!> Definition for a spray class
module spray_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use iterator_class,    only: iterator
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   implicit none
   private
   
   public :: spray
   
   !> spray object
   type :: spray
      
      !> Provide a pardata and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB based on polygon
      type(config) :: cfg
      
      !> Flow solver
      type(incomp)      :: fs        !< Single flow solver
      type(hypre_str)   :: ps        !< HYPRE linear solver for pressure
      type(ddadi)       :: vs        !< DDADI linear solver for velocity
      type(sgsmodel)    :: sgs       !< SGS model for eddy viscosity
      type(timetracker) :: time      !< Time info
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile     !< General simulation monitoring
      type(monitor) :: cflfile   !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:,:), allocatable :: SR                !< Strain rate tensor
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Timing info
      type(monitor) :: timefile !< Timing monitoring
      type(timer)   :: tstep    !< Timer for step
      type(timer)   :: tsgs     !< Timer for SGS
      type(timer)   :: tvel     !< Timer for velocity
      type(timer)   :: tpres    !< Timer for pressure
      type(timer)   :: tlpadv   !< Timer for advancing particles
      
      !> Drop transfer modeling
      logical :: use_drop_transfer !< Do we use droplet transfer
      type(lpt)      :: lp         !< Lagrangian particle tracking
      type(monitor)  :: pfile      !< Particle monitoring
      type(partmesh) :: pmesh      !< Particle mesh for lpt
      real(WP) :: dmax             !< Maximum diameter for transfer
      real(WP) :: dmin             !< Minimum diameter below which transfer is automatic
      real(WP) :: ddel             !< Minimum diameter below which structure is directly deleted
      real(WP) :: emax             !< Maximum eccentricity for transfer
      real(WP) :: vof_tf_drop      !< Integral of VOF transfered from droplet conversion
      real(WP) :: vof_deleted      !< Integral of VOF deleted
      integer  :: np_drop

      logical ::  use_lig_transfer  !< Do we use ligament transfer
      real(WP) :: lmin
      real(WP) :: lmake
      real(WP) :: lper
      real(WP) :: lstratio
      real(WP) :: dw
      real(WP) :: ldmin
      real(WP) :: size_ratio
      real(WP) :: vof_tf_lig
      integer  :: np_lig
      
   contains
      procedure :: init                            !< Initialize spray simulation
      procedure :: step                            !< Advance spray simulation by one time step
      procedure :: final                           !< Finalize spray simulation
   end type spray
   
   
contains
   
   !> Initialization of spray simulation
   subroutine init(this)
      implicit none
      class(spray), intent(inout) :: this
      
      
      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='spray.input')
      end block read_input
      
      
      ! Initialize config object
      create_config: block
         use parallel,    only: group
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni,x,y,z
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x_uni(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y_uni(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z_uni(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! Add stretching
         call this%input%read('Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call this%input%read('Stretch ratio in yz',sratio_yz)
         call this%input%read('Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call this%input%read('Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x *(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='spray')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create config
         this%cfg=config(grp=group,decomp=partition,grid=grid)
         ! Create walls
         this%cfg%VF=1.0_WP
      end block create_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR       (1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Create a two-phase flow solver with bconds
      create_flow_solver: block
         use incomp_class,      only: clipped_neumann,dirichlet,slip
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Single-Phase NS')
         ! Set the flow properties
         call this%input%read('Gas dynamic viscosity', this%fs%visc)
         call this%input%read('Gas density', this%fs%rho)
         ! Set acceleration of gravity
         call this%input%read('Gravity',this%fs%gravity)
         ! Inlets and coflow on the left
         call this%fs%add_bcond(name='inlets',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=left_boundary)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: rad
         ! Zero velocity except if restarting
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs

      create_lpt: block
         this%lp=lpt(cfg=this%cfg,name='spray')
         call this%input%read('Liquid density', this%lp%rho)
         this%lp%gravity=this%fs%gravity
         this%lp%filter_width=3.5_WP*this%cfg%min_meshsize
         call this%lp%resize(0)
      end block create_lpt
      
      ! ! Prepare Lagrangian drop model
      ! ! id=1, transfer_drops
      ! ! id=2, transfer_drops as ligament
      ! ! id=3, transfer_ligament
      ! prepare_transfer: block
      !    ! Is transfer used?
      !    call this%input%read('Transfer drops',this%use_drop_transfer,default=.true.)
      !    call this%input%read('Transfer ligaments',this%use_lig_transfer,default=.true.)
      !    ! Create CCLs
      !    call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      !    call this%ccl_lig%initialize(pg=this%cfg%pgrid,name='ccl_lig')
      !    ! Setup lpt solver
      !    if (this%use_drop_transfer.or.this%use_lig_transfer) then
      !       ! Create lpt solver
      !       this%lp=lpt(cfg=this%cfg,name='spray')
      !       this%lp%rho=this%fs%rho_l
      !       this%lp%gravity=this%fs%gravity
      !       this%lp%filter_width=3.5_WP*this%cfg%min_meshsize
      !       call this%lp%resize(0)
      !    end if
      !    ! Initialize drop transfer routine
      !    if (this%use_drop_transfer) then
      !       ! Create CCL
      !       ! call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      !       ! Set parameters for transfer
      !       this%ddel=0.2_WP*this%cfg%min_meshsize
      !       this%dmin=1.5_WP*this%cfg%min_meshsize
      !       this%dmax=1.0e-3_WP
      !       this%emax=0.8_WP
      !       ! Zero out transfered volume
      !       this%vof_tf_drop=0.0_WP
      !       this%np_drop=0
      !    end if
      !    if (this%use_lig_transfer) then
      !       ! Create CCL LIG
      !       ! call this%ccl_lig%initialize(pg=this%cfg%pgrid,name='ccl_lig')
      !       this%ldmin=1.0e-2_WP
      !       this%dw =0.697_WP
      !       this%size_ratio=0.707_WP 
      !       this%lmin=1.0_WP
      !       ! this%lmake=1.5_WP
      !       this%lmake=3.0_WP
      !       this%lper=0.9_WP
      !       this%lstratio=1.5_WP
      !       ! Zero out monitoring variables
      !       this%vof_tf_lig=0.0_WP
      !       this%np_lig=0
      !    end if
      ! end block prepare_transfer

      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use incomp_class,            only: bcond
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP) :: rad
         integer :: i,j,k,n
         type(bcond), pointer :: mybc
         logical :: partfile_exists
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! We are restarting, read the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart_spray/data_'//trim(timestamp))
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
            ! Finally, handle particle I/O
            if (this%use_drop_transfer.or.this%use_lig_transfer) then
               ! Check if particle file exists
               inquire(file='restart_spray/part_'//trim(timestamp),exist=partfile_exists)
               ! If so, read it
               if (partfile_exists) call this%lp%read(filename='restart_spray/part_'//trim(timestamp))
            end if
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart_spray')) call makedir('restart_spray')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart


      
      ! Create partmesh object for particle output
      if (this%use_drop_transfer.or.this%use_lig_transfer) then
         create_pmesh: block
            integer :: i
            this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
            this%pmesh%varname(1)='radius'
            this%pmesh%varname(2)='id'
            this%pmesh%vecname(1)='velocity'
            call this%lp%update_partmesh(this%pmesh)
            do i=1,this%lp%np_
               this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
               this%pmesh%var(2,i)=this%lp%p(i)%id
               this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
            end do
         end block create_pmesh
      end if

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='spray')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_scalar('divergence',this%fs%div)
         if (this%use_drop_transfer.or.this%use_lig_transfer) call this%ens_out%add_particle('part',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight

      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_spray')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')  
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_spray')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
         ! Create particle monitor
         if (this%use_drop_transfer.or.this%use_lig_transfer.or.this%use_lig_transfer) then
            call this%lp%get_max()
            this%pfile=monitor(amroot=this%lp%cfg%amRoot,name='particles')
            call this%pfile%add_column(this%time%n,'Timestep number')
            call this%pfile%add_column(this%time%t,'Time')
            call this%pfile%add_column(this%lp%np,'Particle number')
            call this%pfile%add_column(this%lp%vp_tot,'Particle volume')
            call this%pfile%add_column(this%lp%np_new,'Npart new')
            call this%pfile%add_column(this%np_drop,'Npart new drop')
            call this%pfile%add_column(this%np_lig, 'Npart new lig')
            call this%pfile%add_column(this%lp%vp_new,'Vpart new')
            call this%pfile%add_column(this%lp%np_out,'Npart removed')
            call this%pfile%add_column(this%lp%vp_out,'Vpart removed')
            call this%pfile%add_column(this%lp%Umin,'Particle Umin')
            call this%pfile%add_column(this%lp%Umax,'Particle Umax')
            call this%pfile%add_column(this%lp%Vmin,'Particle Vmin')
            call this%pfile%add_column(this%lp%Vmax,'Particle Vmax')
            call this%pfile%add_column(this%lp%Wmin,'Particle Wmin')
            call this%pfile%add_column(this%lp%Wmax,'Particle Wmax')
            call this%pfile%add_column(this%lp%dmin,'Particle dmin')
            call this%pfile%add_column(this%lp%dmax,'Particle dmax')
            call this%pfile%write()
         end if
      end block create_monitor

      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep  =timer(comm=this%cfg%comm,name='Timestep')
         this%tvel   =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres  =timer(comm=this%cfg%comm,name='Pressure')
         this%tsgs   =timer(comm=this%cfg%comm,name='SGSmodel')
         this%tlpadv =timer(comm=this%cfg%comm,name='AdvancePart')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
         call this%timefile%add_column(this%tsgs%time  ,trim(this%tsgs%name))
         call this%timefile%add_column(this%tlpadv%time,trim(this%tlpadv%name))
      end block create_timing
      
      
   contains
      
      !> Function that localizes the right domain boundary
      function right_boundary(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function right_boundary
      
      
      !> Function that localizes the left boundary
      function left_boundary(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function left_boundary
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the top (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the bottom (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      implicit none
      class(spray), intent(inout) :: this
      
      ! Reset all timers and start timestep timer
      call this%tstep%reset()
      call this%tsgs%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%tlpadv%reset()
      call this%tstep%start()
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Advance lagrangian droplets
      if (this%use_drop_transfer.or.this%use_lig_transfer) then
         this%resU=this%fs%rho_g
         this%resV=this%fs%visc_g
         call this%tlpadv%start()
         ! if (this%vf%cfg%amRoot) print *, "pre advance lpt"
         call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)
         ! if (this%vf%cfg%amRoot) print *, "post advance lpt"
         call this%tlpadv%stop()
      end if
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W

      ! Reset here gas viscosity
      this%fs%visc=visc_g
      
      ! ! Turbulence modeling
      ! call this%tsgs%start() ! Start SGS timer
      ! sgs_modeling: block
      !    use sgsmodel_class, only: vreman,dynamic_smag
      !    integer :: i,j,k
      !    ! Get velocity gradient tensor and strain rate tensor
      !    call this%fs%get_gradu(this%gradU)
      !    this%SR(6,:,:,:)=(this%gradU(1,1,:,:,:)+this%gradU(2,2,:,:,:)+this%gradU(3,3,:,:,:))/3.0_WP ! div
      !    this%SR(1,:,:,:)=this%gradU(1,1,:,:,:)-this%SR(6,:,:,:)                                     ! du/dx-div/3
      !    this%SR(2,:,:,:)=this%gradU(2,2,:,:,:)-this%SR(6,:,:,:)                                     ! dv/dy-div/3
      !    this%SR(3,:,:,:)=this%gradU(3,3,:,:,:)-this%SR(6,:,:,:)                                     ! dw/dz-div/3
      !    this%SR(4,:,:,:)=0.5_WP*(this%gradU(1,2,:,:,:)+this%gradU(2,1,:,:,:))                       ! (du/dy+dv/dx)/2
      !    this%SR(5,:,:,:)=0.5_WP*(this%gradU(2,3,:,:,:)+this%gradU(3,2,:,:,:))                       ! (dv/dz+dw/dy)/2
      !    this%SR(6,:,:,:)=0.5_WP*(this%gradU(3,1,:,:,:)+this%gradU(1,3,:,:,:))                       ! (dw/dx+du/dz)/2
      !    ! Get turbulent viscosity
      !    this%resU=this%vf%VF*this%fs%rho_l+(1.0_WP-this%vf%VF)*this%fs%rho_g
      !    !call this%sgs%get_visc(type=dynamic_smag,dt=this%time%dtold,rho=this%resU,Ui=this%Ui,Vi=this%Vi,Wi=this%Wi,SR=this%SR)
      !    call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
      !    ! Add sgs visc to our two-phase viscosities
      !    do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
      !       this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
      !       this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
      !       this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
      !       this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
      !    end do; end do; end do
      ! end block sgs_modeling
      ! call this%tsgs%stop() ! Stop SGS timer
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Start velocity timer
         call this%tvel%start()
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Add momentum source terms
         call this%fs%addsrc_gravity(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop velocity timer and start pressure timer
         call this%tvel%stop()
         call this%tpres%start()
         
         ! Solve Poisson equation
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Stop pressure timer
         call this%tpres%stop()
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()

      ! ! attempt transfter
      ! attempt_transfer : block
      !    ! Zero out monitoring variables
      !    this%lp%np_new=0
      !    this%lp%vp_new=0.0_WP
      !    ! Transfer VOF into droplets
      !    call this%tltrans%start() ! Start burst timer
      !    if (this%use_lig_transfer) call this%transfer_ligs()
      !    call this%tltrans%stop() ! Stop burst timer
      !    call this%tdtrans%start() ! Start transfer timer
      !    if (this%use_drop_transfer) call this%transfer_drops()
      !    call this%tdtrans%stop() ! Stop transfer timer
      ! end block attempt_transfer
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update particle mesh object
         if (this%use_drop_transfer.or.this%use_lig_transfer) then
            update_pmesh: block
               integer :: i
               call this%lp%update_partmesh(this%pmesh)
               do i=1,this%lp%np_
                  this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
                  this%pmesh%var(2,i)=this%lp%p(i)%id
                  this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
               end do
            end block update_pmesh 
         end if
         ! Write ensight files
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Stop timestep timer
      call this%tstep%stop()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      if (this%use_drop_transfer.or.this%use_lig_transfer) then
         call this%lp%get_max()
         call this%pfile%write()
      end if
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            integer :: i,j,k
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%write(fdata='restart_spray/data_'//trim(adjustl(timestamp)))
            ! Finally, handle particle I/O
            if (this%use_drop_transfer.or.this%use_lig_transfer) call this%lp%write(filename='restart_spray/part_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize spray simulation
   subroutine final(this)
      implicit none
      class(spray), intent(inout) :: this
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      deallocate(this%gradU,this%SR)
   end subroutine final
   
   
end module spray_class

