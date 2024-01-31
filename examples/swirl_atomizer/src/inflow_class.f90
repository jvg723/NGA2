!> Definition for block 1 simulation: turbulent annular pipe for inflow
module inflow_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use ddadi_class,       only: ddadi
   use fftxyz_class,      only: fftxyz
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use datafile_class,    only: datafile
   use monitor_class,     only: monitor
   implicit none
   private

   public :: inflow
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type :: inflow
      
      !> Provide a datafile and an event tracker for saving restarts
      type(datafile) :: df
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(ibconfig) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs    !< Incompressible flow solver
      type(fftxyz)      :: ps    !< Fourier pressure solver
      type(ddadi)       :: vs    !< Velocity/Scalar solver
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(ensight) :: ens_out  !< Ensight output for flow variables
      type(event)   :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
  
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities


   contains
      procedure, private :: geometry_init          !< Initialize geometry for nozzle
      procedure, private :: simulation_init        !< Initialize simulation for nozzle
      procedure :: init                            !< Initialize block
      procedure :: step                            !< Advance block
      procedure :: final                           !< Finalize block
   end type inflow
   
   !> Fluid definition
   real(WP) :: visc

   !> Forcing terms
   real(WP) :: bforce,rhoUaxial_tgt,rhoUaxial_avg,rhoUtheta_tgt,rhoUtheta_avg
   real(WP) :: swirl_number,swirl_number_tgt

   !> Geometric data for case
   real(WP) :: Dout,Din
   
   
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
   
   
   !> Initialization of inflow simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(inflow), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='input_inflow')
      
      ! Initialize the geometry
      call this%geometry_init()
      
      ! Initialize the simulation
      call this%simulation_init()
      
   end subroutine init

   !> Initialization of problem geometry
   subroutine geometry_init(this)
      use sgrid_class, only: sgrid
      use parallel,    only: comm,group,nproc,rank
      use mpi_f08,     only: MPI_Group,MPI_Group_range_incl
      implicit none
      class(inflow) :: this
      type(sgrid)   :: grid
      

      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call this%input%read('nx',nx); allocate(x(nx+1))
         call this%input%read('ny',ny); allocate(y(ny+1))
         call this%input%read('nz',nz); allocate(z(nz+1))

         ! Read in geometry definition
         call this%input%read('Pipe length',Lx)
         call this%input%read('Outer Pipe diameter',Dout)
         call this%input%read('Inner Pipe diameter',Din)
     
         ! Domain lengths
         dx=Lx/real(nx,WP)
         no=6
         if (ny.gt.1) then
            Ly=Dout+real(2*no,WP)*Dout/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=Dout+real(2*no,WP)*Dout/real(ny-2*no,WP)
         else
            Lz=dx
         end if            
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do 

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='annular_pipe')
      
      end block create_grid
         

      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         use ibconfig_class, only: bigot,sharp
         integer :: i,j,k
         integer, dimension(3) :: partition
         real(WP) :: radius

         ! Read in partition
          call this%input%read('Partition',partition)

         ! Create partitioned grid
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)

         ! Create IB walls for this config
         !> Create IB field
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  radius=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
                  this%cfg%Gib(i,j,k)=max(0.5_WP*Din-radius,radius-0.5_WP*Dout)
               end do
            end do
         end do           
         !> Get normal vector
         call this%cfg%calculate_normal()           
         !> Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)

      end block create_cfg
      
   

   end subroutine geometry_init

   !> Initialize simulation
   subroutine simulation_init(this)
      implicit none
      class(inflow), intent(inout) :: this

      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot,name='annular_pipe')
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW         (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi           (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)) 
      end block allocate_work_arrays

      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call this%input%read('Density',this%fs%rho)
         call this%input%read('Dynamic viscosity',visc); this%fs%visc=visc
         ! Configure pressure solver
         this%ps=fftxyz(cfg=this%cfg,name='Pressure',nst=7)
         ! Configure velocity solver
			this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Uaxial,Utheta,amp,theta,radius
         ! Zero out velocity
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP; this%fs%P=0.0_WP
         ! Read velocity field parameters
         call this%input%read('Bulk axial velocity',Uaxial)
         call this%input%read('Bulk theta velocity',Utheta)
         call this%input%read('Fluctuation amp',amp,default=0.0_WP)
         ! Initialize velocity
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  ! U velocity
                  this%fs%U(i,j,k)=+Uaxial+Uaxial*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Uaxial*cos(8.0_WP*twoPi*this%fs%cfg%zm(k)/this%fs%cfg%zL)*cos(8.0_WP*twoPi*this%fs%cfg%ym(j)/this%fs%cfg%yL)
                  ! V velocity
                  radius=sqrt(this%cfg%y(j)**2+this%cfg%zm(k)**2)
                  theta=atan2(this%cfg%y(j),this%cfg%zm(k))
                  this%fs%V(i,j,k)=+Utheta*radius*cos(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
                  ! W velocity
                  radius=sqrt(this%cfg%ym(j)**2+this%cfg%z(k)**2)
                  theta=atan2(this%cfg%ym(j),this%cfg%z(k))
                  this%fs%W(i,j,k)=-Utheta*radius*sin(theta)+Utheta*radius*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Utheta*radius*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
               end do
            end do
         end do
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
         ! Get target rhoUaxial
         call this%cfg%integrate(A=this%fs%rho*this%Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/this%cfg%fluid_vol
         rhoUaxial_tgt=rhoUaxial_avg
         ! Get target rhoUtheta
         this%resU=0.0_WP
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  radius=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
                  theta=atan2(this%cfg%ym(j),this%cfg%zm(k))
                  if (radius.gt.0.0_WP) this%resU(i,j,k)=(this%Vi(i,j,k)*cos(theta)-this%Wi(i,j,k)*sin(theta))/radius
               end do
            end do
         end do
         call this%cfg%integrate(A=this%resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/this%cfg%fluid_vol
         rhoUtheta_tgt=rhoUtheta_avg
         ! Compute swirl number and coeff
         swirl_number=get_swirl_number(cfg=this%cfg,U=this%Ui,V=this%Vi,W=this%Wi,R=0.5_WP*Dout,rho=this%fs%rho)
         swirl_number_tgt=swirl_number
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='pipe')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('levelset',this%cfg%Gib)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('visc_sgs',this%sgs%visc)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_inflow')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(rhoUaxial_avg,'Average rhoUaxial')
         call this%mfile%add_column(rhoUtheta_avg,'Average rhoUtheta')
         call this%mfile%add_column(swirl_number,'Swirl number')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_inflow')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Take a time step with block 1
   subroutine step(this)
      implicit none
      class(inflow), intent(inout) :: this
         
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
         
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
         
      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         this%resU=this%fs%rho
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         this%fs%visc=visc+this%sgs%visc
      end block sgs_modeling
         
      ! Calculate body forcing
      calc_bodyforcing: block
         integer :: i,j,k
         real(WP) :: theta,radius
         call this%cfg%integrate(A=this%fs%rho*this%Ui,integral=rhoUaxial_avg); rhoUaxial_avg=rhoUaxial_avg/this%cfg%fluid_vol
         this%resU=0.0_WP
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  radius=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
                  theta=atan2(this%cfg%ym(j),this%cfg%zm(k))
                  if (radius.gt.0.0_WP) this%resU(i,j,k)=(this%Vi(i,j,k)*cos(theta)-this%Wi(i,j,k)*sin(theta))/radius
               end do
            end do
         end do
         call this%cfg%integrate(A=this%resU,integral=rhoUtheta_avg); rhoUtheta_avg=rhoUtheta_avg/this%cfg%fluid_vol
      end block calc_bodyforcing
         
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
            
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
            
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
            
         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW
            
         ! Add body forcing (do we need to time it by the volume?)
         add_bodyforcing: block
            integer :: i,j,k
            real(WP) :: theta,radius
            this%resU=this%resU+(rhoUaxial_tgt-rhoUaxial_avg)
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     radius=sqrt(this%cfg%y(j)**2+this%cfg%zm(k)**2)
                     theta=atan2(this%cfg%y(j),this%cfg%zm(k))
                     this%resV(i,j,k)=this%resV(i,j,k)+(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*cos(theta)
                     radius=sqrt(this%cfg%ym(j)**2+this%cfg%z(k)**2)
                     theta=atan2(this%cfg%ym(j),this%cfg%z(k))
                     this%resW(i,j,k)=this%resW(i,j,k)-(swirl_number_tgt-swirl_number)*rhoUaxial_tgt*radius*sin(theta)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%resU)
            call this%fs%cfg%sync(this%resV)
            call this%fs%cfg%sync(this%resW)
         end block add_bodyforcing
            
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
            
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
            
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*this%fs%U(i,j,k)
                     this%fs%V(i,j,k)=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*this%fs%V(i,j,k)
                     this%fs%W(i,j,k)=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*this%fs%W(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing
           
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
            
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
            
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
            
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
            
      end do
         
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
         
      ! Output to ensight
      if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
         
      ! Compute swirl number
      swirl_number=get_swirl_number(cfg=this%cfg,U=this%Ui,V=this%Vi,W=this%Wi,R=0.5_WP*Dout,rho=this%fs%rho)
         
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()

      
   end subroutine step
   
   
   !> Finalize b1 simulation
   subroutine final(this)
      implicit none
      class(inflow), intent(inout) :: this
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
end module inflow_class
