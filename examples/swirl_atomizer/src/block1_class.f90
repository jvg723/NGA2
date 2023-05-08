!> Definition for block 1 simulation: turbulent annular pipe for inflow
module block1_class
   use precision,         only: WP
   use geometry,          only: Dout
   use ibconfig_class,    only: ibconfig
   use ddadi_class,       only: ddadi
   use fftxyz_class,      only: fftxyz
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use fene_class,        only: fene
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
      type(fftxyz)               :: ps            !< Fourier pressure solver
      type(ddadi)                :: vs,ss         !< Velocity/Scalar solver
      type(fene)                 :: nn            !< FENE model for non-newt behavior
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
   real(WP) :: visc,bforce,rhoUaxial_tgt,rhoUaxial_avg,rhoUtheta_tgt,rhoUtheta_avg
   real(WP) :: swirl_number,swirl_number_tgt
   
   
contains
   
   
   !> Function that computes swirl number
   function get_swirl_number(cfg,U,V,W,R,rho) result(SN)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      use ibconfig_class,    only: ibconfig
      class(ibconfig), intent(in) :: cfg
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: U,V,W !< Cell-centered velocity field
      real(WP), intent(in) :: R     !< Typically the outer radius
      real(WP), intent(in) :: rho   !< Fluid density
      real(WP) :: SN
      integer :: i,j,k,ierr
      real(WP) :: theta,radius
      real(WP) :: myaxialflux,mythetaflux,axialflux,thetaflux
      myaxialflux=0.0_WP; mythetaflux=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               radius=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
               theta=atan2(cfg%ym(j),cfg%zm(k))
               myaxialflux=myaxialflux+cfg%vol(i,j,k)*cfg%VF(i,j,k)*rho*U(i,j,k)**2
               mythetaflux=mythetaflux+cfg%vol(i,j,k)*cfg%VF(i,j,k)*rho*U(i,j,k)*(V(i,j,k)*cos(theta)-W(i,j,k)*sin(theta))
            end do
         end do
      end do
      call MPI_ALLREDUCE(myaxialflux,axialflux,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(mythetaflux,thetaflux,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      SN=thetaflux/(R*max(axialflux,epsilon(1.0_WP)))
   end function get_swirl_number
   
   
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
         b%time=timetracker(amRoot=b%cfg%amRoot,name='annular_pipe')
         call param_read('1 Max timestep size',b%time%dtmax)
         call param_read('Max time',b%time%tmax)
         call param_read('Max cfl number',b%time%cflmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker

      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg,gmres_pfmg
         ! Create flow solver
         b%fs=incomp(cfg=b%cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Liquid density',b%fs%rho)
         call param_read('Liquid dynamic viscosity',visc)
         ! Configure pressure solver
         b%ps=fftxyz(cfg=b%cfg,name='Pressure',nst=7)
         ! Configure velocity solver
			b%vs=ddadi(cfg=b%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%vs)
      end block create_flow_solver

      ! Create a FENE model 
      create_fene: block 
         use multiscalar_class, only: bquick
         use fene_class,        only: fenecr
         integer :: i,j,k
         ! Create FENE model solver
         b%nn=fene(cfg=b%cfg,model=fenecr,scheme=bquick,name='FENE')
         ! Assign unity density for simplicity
         b%nn%rho=1.0_WP
         ! Maximum extensibility of polymer chain
         call param_read('Maximum polymer extensibility',b%nn%Lmax)
         ! Relaxation time for polymer
         call param_read('Polymer relaxation time',b%nn%trelax)
         ! Polymer viscosity at zero strain rate
         call param_read('Polymer viscosity',b%nn%visc);b%nn%visc_p=b%nn%visc; b%fs%visc=visc+b%nn%visc_p
         ! Powerlaw coefficient in Carreau model
         call param_read('Carreau powerlaw',b%nn%ncoeff)
         ! Configure implicit scalar solver
         b%ss=ddadi(cfg=b%cfg,name='scalar',nst=13)
         ! Setup the solver
         call b%nn%setup(implicit_solver=b%ss)
         ! Initialize conformation tensor to identity
         b%nn%SC(:,:,:,1)=1.0_WP !< Cxx
         b%nn%SC(:,:,:,4)=1.0_WP !< Cyy
         b%nn%SC(:,:,:,6)=1.0_WP !< Czz
      end block create_fene
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Uaxial,Utheta,amp,theta,radius
         ! Zero out velocity
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP; b%fs%P=0.0_WP
         ! Read velocity field parameters
         call param_read('Bulk axial velocity',Uaxial)
         call param_read('Bulk theta velocity',Utheta)
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         ! Initialize velocity
         do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
            do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
               do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                  ! U velocity
                  b%fs%U(i,j,k)=+Uaxial+Uaxial*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Uaxial*cos(8.0_WP*twoPi*b%fs%cfg%zm(k)/b%fs%cfg%zL)*cos(8.0_WP*twoPi*b%fs%cfg%ym(j)/b%fs%cfg%yL)
                  ! V velocity
                  radius=sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)
                  theta=atan2(b%cfg%y(j),b%cfg%zm(k))
                  b%fs%V(i,j,k)=+Utheta*radius*cos(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*b%fs%cfg%xm(i)/b%fs%cfg%xL)
                  ! W velocity
                  radius=sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)
                  theta=atan2(b%cfg%ym(j),b%cfg%z(k))
                  b%fs%W(i,j,k)=-Utheta*radius*sin(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*b%fs%cfg%xm(i)/b%fs%cfg%xL)
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
         ! Get target rhoUaxial
         call b%cfg%integrate(A=b%fs%rho*b%Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/b%cfg%fluid_vol
         rhoUaxial_tgt=rhoUaxial_avg
         ! Get target rhoUtheta
         b%resU=0.0_WP
         do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
            do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
               do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                  radius=sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)
                  theta=atan2(b%cfg%ym(j),b%cfg%zm(k))
                  if (radius.gt.0.0_WP) b%resU(i,j,k)=(b%Vi(i,j,k)*cos(theta)-b%Wi(i,j,k)*sin(theta))/radius
               end do
            end do
         end do
         call b%cfg%integrate(A=b%resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/b%cfg%fluid_vol
         rhoUtheta_tgt=rhoUtheta_avg
         ! Compute swirl number and coeff
         swirl_number=get_swirl_number(cfg=b%cfg,U=b%Ui,V=b%Vi,W=b%Wi,R=0.5_WP*Dout,rho=b%fs%rho)
         swirl_number_tgt=swirl_number
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
         call b%ens_out%add_scalar('visc_l',b%nn%visc_p)
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
         call b%mfile%add_column(rhoUaxial_avg,'Average rhoUaxial')
         call b%mfile%add_column(rhoUtheta_avg,'Average rhoUtheta')
         call b%mfile%add_column(swirl_number,'Swirl number')
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

      ! fs%visc_l is the solvent viscosity, nn%visc is the zero strainrate polymer viscosity
      shear_thinning: block
         integer :: i,j,k
         real(WP), dimension(:,:,:,:), allocatable :: SR
         ! Allocate SR array
         allocate(SR(1:6,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         ! Calculate strain rate
         call b%fs%get_strainrate(SR)
         ! Update polymer viscosity using Carreau model
         call b%nn%update_visc_p(SR)
         ! Update total liquid viscosity
         b%fs%visc=visc+b%nn%visc_p
         ! Deallocate SR array
         deallocate(SR)
      end block shear_thinning
         
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         b%resU=b%fs%rho
         call b%fs%get_gradu(b%gradU)
         call b%sgs%get_visc(type=vreman,dt=b%time%dtold,rho=b%resU,gradu=b%gradU)
         b%fs%visc=visc+b%sgs%visc
      end block sgs_modeling
         
      ! Calculate body forcing
      calc_bodyforcing: block
         integer :: i,j,k
         real(WP) :: theta,radius
         call b%cfg%integrate(A=b%fs%rho*b%Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/b%cfg%fluid_vol
         b%resU=0.0_WP
         do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
            do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
               do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                  radius=sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2)
                  theta=atan2(b%cfg%ym(j),b%cfg%zm(k))
                  if (radius.gt.0.0_WP) b%resU(i,j,k)=(b%Vi(i,j,k)*cos(theta)-b%Wi(i,j,k)*sin(theta))/radius
               end do
            end do
         end do
         call b%cfg%integrate(A=b%resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/b%cfg%fluid_vol
      end block calc_bodyforcing
         
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
            
         ! Add body forcing (do we need to time it by the volume?)
         add_bodyforcing: block
            integer :: i,j,k
            real(WP) :: theta,radius
            b%resU=b%resU+(rhoUaxial_tgt-rhoUaxial_avg)
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     radius=sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)
                     theta=atan2(b%cfg%y(j),b%cfg%zm(k))
                     b%resV(i,j,k)=b%resV(i,j,k)+(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*cos(theta)
                     radius=sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)
                     theta=atan2(b%cfg%ym(j),b%cfg%z(k))
                     b%resW(i,j,k)=b%resW(i,j,k)-(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*sin(theta)
                  end do
               end do
            end do
            call b%fs%cfg%sync(b%resU)
            call b%fs%cfg%sync(b%resV)
            call b%fs%cfg%sync(b%resW)
         end block add_bodyforcing
            
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
         
      ! Compute swirl number
      swirl_number=get_swirl_number(cfg=b%cfg,U=b%Ui,V=b%Vi,W=b%Wi,R=0.5_WP*Dout,rho=b%fs%rho)
         
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
