!> Definition for block 2 simulation: swirl atomizer
module block2_class
   use precision,         only: WP
   use geometry,          only: Dout
   use config_class,      only: config
   use hypre_str_class,   only: hypre_str
   ! use fourier3d_class,   only: fourier3d
   use incomp_class,      only: incomp
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   public :: block2
   
   !> block 2 object
   type :: block2
      class(config), pointer :: cfg             !< Pointer to config
      type(incomp)           :: fs            !< Single phase incompressible flow solver
      type(hypre_str)        :: ps              !< Structured hypre pressure solver
      ! type(fourier3d)            :: ps            !< Fourier pressure solver
	   type(hypre_str)        :: vs              !< Structured hypre implicit solver
      type(timetracker)      :: time            !< Time tracker                                     
      type(ensight)          :: ens_out         !< Ensight output
      type(event)            :: ens_evt         !< Ensight event
      type(monitor)          :: mfile,cflfile   !< Monitor files
      !> Private work arrays
      real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2

   !> Problem definition
   real(WP) :: visc,bforce,rhoUaxial_tgt,rhoUaxial_avg,rhoUtheta_tgt,rhoUtheta_avg

   
contains
   

   !> Function that localizes the top (x+) of the domain
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   

   !> Function that localizes the bottom (x-) of the domain boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.le.pg%imin) isIn=.true.
   end function xm_locator
   
   
   !> Initialization of problem solver
   subroutine init(b)
      use param, only: param_read
      implicit none
      class(block2), intent(inout) :: b
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(b%resU (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resV (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resW (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Ui   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Vi   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Wi   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         b%time=timetracker(amRoot=b%cfg%amRoot,name='jet')
         call param_read('2 Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use incomp_class, only: dirichlet,clipped_neumann
         use hypre_str_class,  only: pcg_pfmg
         integer :: i,j,k
         ! Create flow solver
         b%fs=incomp(cfg=b%cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Liquid density',b%fs%rho)
         call param_read('Liquid dynamic viscosity',visc); b%fs%visc=visc
         ! Configure pressure solver
         ! b%ps=fourier3d(cfg=b%cfg,name='Pressure',nst=7)
         b%ps=hypre_str(cfg=b%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         b%ps%maxlevel=14
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Configure implicit velocity solver
         b%vs=hypre_str(cfg=b%cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',b%vs%maxit)
         call param_read('Implicit tolerance',b%vs%rcvg)
         ! Inflow on the left of domain
         call b%fs%add_bcond(name='inflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Clipped Neumann outflow on the right of domain
         call b%fs%add_bcond(name='bc_xp' ,type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
			! Setup the solver
			call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%vs)
      end block create_and_initialize_flow_solver
      

      ! Initialize our velocity field 
      initialize_velocity: block
         ! Zero initial field in the domain
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP
         ! Apply all other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)
         ! Compute MFR through all boundary conditions
         call b%fs%get_mfr()
         ! Adjust MFR for global mass balance
         call b%fs%correct_mfr()
         ! Compute cell-centered velocity
         call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
         ! Compute divergence
         call b%fs%get_div()
      end block initialize_velocity
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(cfg=b%cfg,name='jet')
         ! Create event for Ensight output
         b%ens_evt=event(time=b%time,name='Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('Pressure',b%fs%P)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation2')
         call b%mfile%add_column(b%time%n,'Timestep number')
         call b%mfile%add_column(b%time%t,'Time')
         call b%mfile%add_column(b%time%dt,'Timestep size')
         call b%mfile%add_column(b%time%cfl,'Maximum CFL')
         call b%mfile%add_column(b%fs%Umax,'Umax')
         call b%mfile%add_column(b%fs%Vmax,'Vmax')
         call b%mfile%add_column(b%fs%Wmax,'Wmax')
         call b%mfile%add_column(b%fs%Pmax,'Pmax')
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl2')
         call b%cflfile%add_column(b%time%n,'Timestep number')
         call b%cflfile%add_column(b%time%t,'Time')
         call b%cflfile%add_column(b%fs%CFLc_x,'Convective xCFL')
         call b%cflfile%add_column(b%fs%CFLc_y,'Convective yCFL')
         call b%cflfile%add_column(b%fs%CFLc_z,'Convective zCFL')
         call b%cflfile%add_column(b%fs%CFLv_x,'Viscous xCFL')
         call b%cflfile%add_column(b%fs%CFLv_y,'Viscous yCFL')
         call b%cflfile%add_column(b%fs%CFLv_z,'Viscous zCFL')
         call b%cflfile%write()
      end block create_monitor
      
      
   end subroutine init
   
   
   !> Take a time step with block 2
   subroutine step(b,Uinflow,Vinflow,Winflow)
      implicit none
      class(block2), intent(inout) :: b
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Uinflow     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Vinflow     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(b%cfg%imino_:,b%cfg%jmino_:,b%cfg%kmino_:), intent(inout) :: Winflow     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
         
      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Apply time-varying Dirichlet conditions
      reapply_dirichlet: block
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k
         ! Apply velocity from turbulent annular pipe
         call b%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            b%fs%U(i,j,k)=Uinflow(i,j,k)
            b%fs%V(i,j,k)=Vinflow(i,j,k)
            b%fs%W(i,j,k)=Winflow(i,j,k)
         end do
      end block reapply_dirichlet
         
      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W
         
      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)
            
         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)
            
         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)
            
         ! Assemble explicit residual
         b%resU=-2.0_WP*(b%fs%rho*b%fs%U-b%fs%rho*b%fs%Uold)+b%time%dt*b%resU
         b%resV=-2.0_WP*(b%fs%rho*b%fs%V-b%fs%rho*b%fs%Vold)+b%time%dt*b%resV
         b%resW=-2.0_WP*(b%fs%rho*b%fs%W-b%fs%rho*b%fs%Wold)+b%time%dt*b%resW
            
         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)
            
         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW
            
         ! Apply other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)
            
         ! Solve Poisson equation
         call b%fs%correct_mfr()
         call b%fs%get_div()
         b%fs%psolv%rhs=-b%fs%cfg%vol*b%fs%div*b%fs%rho/b%time%dt
         b%fs%psolv%sol=0.0_WP
         call b%fs%psolv%solve()
         call b%fs%shift_p(b%fs%psolv%sol)
            
         ! Correct velocity
         call b%fs%get_pgrad(b%fs%psolv%sol,b%resU,b%resV,b%resW)
         b%fs%P=b%fs%P+b%fs%psolv%sol
         b%fs%U=b%fs%U-b%time%dt*b%resU/b%fs%rho
         b%fs%V=b%fs%V-b%time%dt*b%resV/b%fs%rho
         b%fs%W=b%fs%W-b%time%dt*b%resW/b%fs%rho
            
         ! Increment sub-iteration counter
         b%time%it=b%time%it+1
            
      end do
         
      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

      ! Output to ensight
      if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
         
      ! Perform and output monitoring
      call b%fs%get_max()
      call b%mfile%write()
      call b%cflfile%write()
         

   end subroutine step
      
   
   !> Finalize b2 simulation
   subroutine final(b)
      implicit none
      class(block2), intent(inout) :: b
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi)
      
   end subroutine final
   
   
end module block2_class
