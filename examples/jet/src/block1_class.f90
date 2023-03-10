!> Definition for block 1 simulation: turbulent annular pipe for inflow
module block1_class
   use precision,         only: WP
   use geometry,          only: Dout
   use ibconfig_class,    only: ibconfig
   use hypre_str_class,   only: hypre_str
   ! use fourier3d_class,   only: fourier3d
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   public :: block1
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type :: block1
      class(ibconfig),   pointer :: cfg           !< Pointer to config
      type(incomp)               :: fs            !< Single phase incompressible flow solver
      type(hypre_str)            :: ps
      ! type(fourier3d)            :: ps            !< Fourier pressure solver
      type(hypre_str)            :: vs            !< Structured hypre implicit solver
      type(sgsmodel)             :: sgs           !< SGS model
      type(timetracker)          :: time          !< Time tracker
      type(ensight)              :: ens_out       !< Ensight output
      type(event)                :: ens_evt       !< Ensight event 
      type(monitor)              :: mfile,cflfile !< Monitor files
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable     :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable     :: Ui,Vi,Wi
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block1
   
   !> Problem definition
   real(WP) :: visc,mfr_target,mfr,bforce
   
   
contains
   
   
   !> Compute massflow rate
   function get_bodyforce_mfr(cfg,rho,U,srcU) result(mfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      use ibconfig_class,    only: ibconfig
      class(ibconfig), intent(in) :: cfg
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: U !< Cell-centered velocity field
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(in) :: rho
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol,mfr
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  vol=cfg%dxm(i)*cfg%dy(j)*cfg%dz(k)*cfg%VF(i,j,k)
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*(rho*U(i,j,k)+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_
                  vol=cfg%dxm(i)*cfg%dy(j)*cfg%dz(k)**cfg%VF(i,j,k)
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*rho*U(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr); mfr=mfr/Uvol
   end function get_bodyforce_mfr
   
   
   !> Initialization of problem solver
   subroutine init(b)
      use param, only: param_read
      implicit none
      class(block1), intent(inout) :: b


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(b%resU         (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resV         (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%resW         (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Ui           (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Vi           (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%Wi           (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%gradU(1:3,1:3,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_)) 
      end block allocate_work_arrays


      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         b%time=timetracker(amRoot=b%cfg%amRoot,name='turbulent_pipe')
         call param_read('1 Max timestep size',b%time%dtmax)
         call param_read('Max time',b%time%tmax)
         call param_read('Max cfl number',b%time%cflmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker

      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
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
         ! Setup the solver
         call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Ubulk,amp
         ! Initial fields
         call param_read('Bulk velocity',Ubulk)
         b%fs%U=Ubulk; b%fs%V=0.0_WP; b%fs%W=0.0_WP; b%fs%P=0.0_WP
         ! For faster transition
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
            do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
               do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                  b%fs%U(i,j,k)=b%fs%U(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*b%fs%cfg%zm(k)/b%fs%cfg%zL)*cos(8.0_WP*twoPi*b%fs%cfg%ym(j)/b%fs%cfg%yL)
                  b%fs%V(i,j,k)=b%fs%V(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*b%fs%cfg%xm(i)/b%fs%cfg%xL)
                  b%fs%W(i,j,k)=b%fs%W(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*b%fs%cfg%xm(i)/b%fs%cfg%xL)
               end do
            end do
         end do
         call b%fs%cfg%sync(b%fs%U)
         call b%fs%cfg%sync(b%fs%V)
         call b%fs%cfg%sync(b%fs%W)
         ! Compute cell-centered velocity
         call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
         ! Compute divergence
         call b%fs%get_div()
         ! Get target MFR and zero bodyforce
         mfr=get_bodyforce_mfr(b%cfg,b%fs%rho,b%Ui)
         mfr_target=mfr
         bforce=0.0_WP
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         b%sgs=sgsmodel(cfg=b%fs%cfg,umask=b%fs%umask,vmask=b%fs%vmask,wmask=b%fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(cfg=b%cfg,name='pipe')
         ! Create event for Ensight output
         b%ens_evt=event(time=b%time,name='Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('levelset',b%cfg%Gib)
         call b%ens_out%add_scalar('pressure',b%fs%P)
         call b%ens_out%add_scalar('visc_sgs',b%sgs%visc)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation')
         call b%mfile%add_column(b%time%n,'Timestep number')
         call b%mfile%add_column(b%time%t,'Time')
         call b%mfile%add_column(b%time%dt,'Timestep size')
         call b%mfile%add_column(b%time%cfl,'Maximum CFL')
         call b%mfile%add_column(mfr,'MFR')
         call b%mfile%add_column(bforce,'Body force')
         call b%mfile%add_column(b%fs%Umax,'Umax')
         call b%mfile%add_column(b%fs%Vmax,'Vmax')
         call b%mfile%add_column(b%fs%Wmax,'Wmax')
         call b%mfile%add_column(b%fs%Pmax,'Pmax')
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl')
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
   
   
   !> Take a time step with block 1
   subroutine step(b)
      implicit none
      class(block1), intent(inout) :: b
         
      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()
         
      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W
         
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         b%resU=b%fs%rho
         call b%fs%get_gradu(b%gradU)
         call b%sgs%get_visc(type=vreman,dt=b%time%dtold,rho=b%resU,gradu=b%gradU)
         b%fs%visc=visc+b%sgs%visc
      end block sgs_modeling
         
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
            
         ! Add body forcing
         bodyforcing: block
         real(WP) :: mfr
         mfr=get_bodyforce_mfr(b%cfg,b%fs%rho,b%Ui,b%resU)
         bforce=(mfr_target-mfr)/b%time%dtmid
         b%resU=b%resU+b%time%dt*bforce
         end block bodyforcing
            
         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)
            
         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW
            
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     b%fs%U(i,j,k)=sum(b%fs%itpr_x(:,i,j,k)*b%cfg%VF(i-1:i,j,k))*b%fs%U(i,j,k)
                     b%fs%V(i,j,k)=sum(b%fs%itpr_y(:,i,j,k)*b%cfg%VF(i,j-1:j,k))*b%fs%V(i,j,k)
                     b%fs%W(i,j,k)=sum(b%fs%itpr_z(:,i,j,k)*b%cfg%VF(i,j,k-1:k))*b%fs%W(i,j,k)
                  end do
               end do
            end do
            call b%fs%cfg%sync(b%fs%U)
            call b%fs%cfg%sync(b%fs%V)
            call b%fs%cfg%sync(b%fs%W)
         end block ibforcing
           
         ! Apply other boundary conditions on the resulting fields
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
   
   
   !> Finalize b1 simulation
   subroutine final(b)
      implicit none
      class(block1), intent(inout) :: b
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%gradU)
      
   end subroutine final
   
end module block1_class
