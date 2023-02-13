!> Definition for block 1 simulation: turbulent annular pipe for inflow
module block1_class
   use precision,         only: WP
   use geometry,          only: D1,D2,get_VF
   use config_class,      only: config
   use hypre_str_class,   only: hypre_str
   !use pfft3d_class,      only: pfft3d
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   public :: block1
   
   !> block 1 object
   type :: block1
      class(config),    pointer :: cfg           !< Pointer to config
      type(incomp),      public :: fs            !< Single phase incompressible flow solver
      type(hypre_str),   public :: ps            !< Structured hypre pressure solver
      !type(pfft3d),      public :: ps
      type(hypre_str),   public :: vs            !< Structured hypre implicit solver
      type(sgsmodel),    public :: sgs           !< SGS model
      type(timetracker), public :: time          !< Time tracker
      type(ensight)             :: ens_out       !< Ensight output
      type(event)               :: ens_evt       !< Ensight event 
      type(monitor)             :: mfile,cflfile !< Monitor files
      !> Private work arrays
      real(WP), dimension(:,:,:), allocatable     :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable     :: Ui,Vi,Wi
      real(WP), dimension(:,:,:), allocatable     :: G
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block1

   !> Problem definition
   real(WP) :: R1,R2                                  !< Inner and outer radii of annulus
   real(WP) :: Ubulk,SW                               !< Liquid axial velocity and swirl ratio
   real(WP) :: Umfr,Vmfr,Wmfr                         !< Mass flow rates for forcing
   real(WP) :: Umfr_target,Vmfr_target,Wmfr_target    !< Mass flow rate targets
   real(WP) :: Ubforce,Vbforce,Wbforce                !< Body forcing
   real(WP) :: visc                                   !< Fluid Viscosity   
   
   
contains
   
   
   !> Compute massflow rate in x/U direction
   function get_bodyforce_Umfr(srcU) result(Umfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol,Umfr
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  !myRhoU=myRhoU+vol*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))*fs%rho+vol*srcU(i,j,k)
                  myRhoU=myRhoU+vol*(fs%rho*fs%U(i,j,k)+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  !myRhoU=myRhoU+vol*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))*fs%rho
                  myRhoU=myRhoU+vol*fs%rho*fs%U(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,Umfr,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Umfr=Umfr/Uvol
   end function get_bodyforce_Umfr
   
   !> Compute massflow rate in y/V direction
   function get_bodyforce_Vmfr(srcV) result(Vmfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: srcV !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoV,myVvol,Vvol,Vmfr
      myRhoV=0.0_WP; myVvol=0.0_WP
      if (present(srcV)) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'V')
                  myVvol=myVvol+vol
                  !myRhoV=myRhoV+vol*(2.0_WP*fs%V(i,j,k)-fs%Vold(i,j,k))*fs%rho+vol*srcV(i,j,k)
                  myRhoV=myRhoV+vol*(fs%rho*fs%V(i,j,k)+srcV(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'V')
                  myVvol=myVvol+vol
                  !myRhoV=myRhoV+vol*(2.0_WP*fs%V(i,j,k)-fs%Vold(i,j,k))*fs%rho
                  myRhoV=myRhoV+vol*fs%rho*fs%V(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myVvol,Vvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoV,Vmfr,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Vmfr=Vmfr/Vvol
   end function get_bodyforce_Vmfr

   !> Compute massflow rate in z/W direction
   function get_bodyforce_Wmfr(srcW) result(Wmfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: srcW !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoW,myWvol,Wvol,Wmfr
      myRhoW=0.0_WP; myWvol=0.0_WP
      if (present(srcW)) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'W')
                  myWvol=myWvol+vol
                  !myRhoW=myRhoW+vol*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))*fs%rho+vol*srcW(i,j,k)
                  myRhoW=myRhoW+vol*(fs%rho*fs%W(i,j,k)+srcW(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*get_VF(i,j,k,'W')
                  myWvol=myWvol+vol
                  !myRhoW=myRhoW+vol*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))*fs%rho
                  myRhoW=myRhoW+vol*fs%rho*fs%W(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myWvol,Wvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoW,Wmfr,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Wmfr=Wmfr/Wvol
   end function get_bodyforce_Wmfr
   
   
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
         allocate(b%G            (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
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
         use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         b%fs=incomp(cfg=b%cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',b%fs%rho)
         call param_read('Dynamic viscosity',visc); b%fs%visc=visc
         ! Configure pressure solver
         !ps=pfft3d(cfg=cfg,name='Pressure',nst=7)
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
         integer :: i,j,k
         real(WP) :: r,theta
         ! Initial fields
         call param_read('Bulk velocity',Ubulk)
         call param_read('Swirl ratio',SW)
         fs%U=Ubulk; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! Apply axial swirl 
         do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
            do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
               do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                  ! Set U velocity
                  r=sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2) !< Radius in cylindrical coordinates
                  b%fs%U(i,j,k)=Ubulk
                  ! Set V velocity
                  r=sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)  !< Radius in cylindrical coordinates
                  theta=atan2(b%cfg%y(j),b%cfg%zm(k))   !< Angle  in cylindrical coordinates
                  b%fs%V(i,j,k)=+4.0_WP*SW*Ubulk*(-r**2+(D2+D1)*r-D2*D1)/((D2-D1)**2)*cos(theta)
                  ! Set W velocity
                  r=sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)  !< Radius in cylindrical coordinates
                  theta=atan2(cfg%ym(j),cfg%z(k))   !< Angle  in cylindrical coordinates
                  b%fs%W(i,j,k)=-4.0_WP*SW*Ubulk*(-r**2+(D2+D1)*r-D2*D1)/((D2-D1)**2)*sin(theta)
               end do
            end do
         end do
         call b%fs%cfg%sync(b%fs%U)
         call b%fs%cfg%sync(b%fs%V)
         call b%fs%cfg%sync(b%fs%W)
         ! Compute cell-centered velocity
         call b%fs%interp_velb%(Ui,b%Vi,b%Wi)
         ! Compute divergence
         call b%fs%get_div()
         ! Get target MFR and zero bodyforce
         Umfr=get_bodyforce_Umfr()
         Vmfr=get_bodyforce_Vmfr()
         Wmfr=get_bodyforce_Wmfr()
         Umfr_target=Umfr
         Vmfr_target=Vmfr
         Wmfr_target=Wmfr
         Ubforce=0.0_WP
         Vbforce=0.0_WP
         Wbforce=0.0_WP
      end block initialize_velocity
      
      
      ! Initialize IBM fields
      initialize_ibm: block
         integer :: i,j,k
         real(WP) :: r
         do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
            do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
               do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
                  r=sqrt(b%fs%cfg%ym(j)**2+b%fs%cfg%zm(k)**2)
                  if (r.ge.0.5_WP*D1) then 
                     b%G(i,j,k)=0.5_WP*D2-r
                  else
                     b%G(i,j,k)=r-0.5_WP*D1
                  end if
               end do
            end do
         end do
      end block initialize_ibm
      
      
      ! Create an LES model
      create_sgs: block
         b%sgs=sgsmodel(cfg=b%fs%cfg,umask=b%fs%umask,vmask=b%fs%vmask,wmask=b%fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         b%ens_out=ensight(cfg=b%cfg,name='annular_pipe')
         ! Create event for Ensight output
         b%ens_evt=event(time=b%time,name='Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call ens_out%add_scalar('levelset',b%G)
         call ens_out%add_scalar('pressure',b%fs%P)
         call ens_out%add_scalar('visc_sgs',b%sgs%visc)
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         ! Create simulation monitor
         b%mfile=monitor(b%fs%cfg%amRoot,'simulation1')
         call b%mfile%add_column(b%time%n,'Timestep number')
         call b%mfile%add_column(b%time%t,'Time')
         call b%mfile%add_column(b%time%dt,'Timestep size')
         call b%mfile%add_column(b%time%cfl,'Maximum CFL')
         call b%mfile%add_column(Umfr,'U MFR')
         call b%mfile%add_column(Vmfr,'V MFR')
         call b%mfile%add_column(Wmfr,'W MFR')
         call b%mfile%add_column(Ubforce,'U Body force')
         call b%mfile%add_column(Vbforce,'V Body force')
         call b%mfile%add_column(Wbforce,'W Body force')
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
            real(WP) :: Umfr,Vmfr,Wmfr
            ! Forcing of U velocity
            Umfr=get_bodyforce_Umfr(b%resU)
            Ubforce=(Umfr_target-Umfr)/b%time%dtmid
            b%resU=b%resU+b%time%dt*Ubforce
            ! Forcing of V velocity
            Vmfr=get_bodyforce_Vmfr(b%resV)
            Vbforce=(Vmfr_target-Vmfr)/b%time%dtmid
            b%resV=b%resV+b%time%dt*Vbforce
            ! Forcing of W velocity
            Wmfr=get_bodyforce_Wmfr(b%resW)
            Wbforce=(Wmfr_target-Wmfr)/b%time%dtmid
            b%resW=b%resW+b%time%dt*Wbforce
         end block bodyforcing
            
         ! Apply IB forcing to enforce BC at the pipe walls
         !ibforcing: block
         !   integer :: i,j,k
         !   do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
         !      do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
         !         do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
         !            b%resU(i,j,k)=b%resU(i,j,k)-(1.0_WP-get_VF(i,j,k,'U'))*b%fs%rho*b%fs%U(i,j,k)
         !            b%resV(i,j,k)=b%resV(i,j,k)-(1.0_WP-get_VF(i,j,k,'V'))*b%fs%rho*b%fs%V(i,j,k)
         !            b%resW(i,j,k)=b%resW(i,j,k)-(1.0_WP-get_VF(i,j,k,'W'))*b%fs%rho*b%fs%W(i,j,k)
         !         end do
         !      end do
         !   end do
         !end block ibforcing

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
                     b%fs%U(i,j,k)=get_VF(i,j,k,'U')*b%fs%U(i,j,k)
                     b%fs%V(i,j,k)=get_VF(i,j,k,'V')*b%fs%V(i,j,k)
                     b%fs%W(i,j,k)=get_VF(i,j,k,'W')*b%fs%W(i,j,k)
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
      if b%(ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
         
      ! Perform and output monitoring
      call b%fs%get_max()
      call b%mfile%write()
      call b%cflfile%write()
      
   end subroutine step
   
   
   !> Finalize b1 simulation
   subroutine final
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
