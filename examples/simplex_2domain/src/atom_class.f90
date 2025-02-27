!> Definition for an atomization class
module atom_class
   use precision,         only: WP
   use inputfile_class,   only: inputfile
   use config_class,      only: config
   use iterator_class,    only: iterator
   use surfmesh_class,    only: surfmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   public :: atom
   
   !> Atom object
   type :: atom

      !> Provide a pardata and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      type(ddadi)       :: vs    !< DDADI solver for velocity
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer !< Edge of domain where we actively remove VOF
      real(WP)       :: vof_removed       !< Integral of VOF removed
      integer        :: nlayer=4          !> Hardcode size of buffer layer for VOF removal


   contains
      procedure, private :: geometry_init          !< Initialize geometry for nozzle
      procedure, private :: simulation_init        !< Initialize simulation for nozzle
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
   end type atom

   
   
contains
   
   
   !> Initialization of atom simulation
   subroutine init(this)
      use parallel, only: amRoot
      implicit none
      class(atom), intent(inout) :: this
      
      ! Read the input
      this%input=inputfile(amRoot=amRoot,filename='atomization.input')
      
      ! Initialize the geometry
      call this%geometry_init()
      
      ! Initialize the simulation
      call this%simulation_init()
      
   end subroutine init


   !> Initialize geometry
   subroutine geometry_init(this)
      use sgrid_class, only: sgrid
      implicit none
      class(atom) :: this
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni
         real(WP), dimension(:), allocatable :: x,y,z
         
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
            x(i)=x(i-1)+sratio_x*(x(i-1)-x(i-2))
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
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='atom')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call this%input%read('Partition',partition)
         
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=partition,grid=grid)

         ! Create masks for this config
         this%cfg%VF=1.0_WP
         
      end block create_cfg
      

   end subroutine geometry_init
   
   
   !> Initialize simulation
   subroutine simulation_init(this)
      implicit none
      class(atom), intent(inout) :: this
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: r2pnet,plicnet,remap
         integer :: i,j,k
         real(WP) :: xloc,rad
         ! Create a VOF solver with LVIRA
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         ! Initialize to flat interface in liquid needle
         ! xloc=0.0_WP !< Interface initially at x=0
         ! do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
         !    do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
         !       do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
         !          rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
         !          if (this%vf%cfg%xm(i).lt.xloc.and.rad.le.0.5_WP*dl) then
         !             this%vf%VF(i,j,k)=1.0_WP
         !          else
         !             this%vf%VF(i,j,k)=0.0_WP
         !          end if
         !          this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
         !          this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
         !       end do
         !    end do
         ! end do
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator

      
      ! Create an incompressible flow solver with bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use tpns_class,      only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set the flow properties
         call this%input%read('Liquid dynamic viscosity',this%fs%visc_l)
         call this%input%read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call this%input%read('Liquid density',this%fs%rho_l)
         call this%input%read('Gas density'   ,this%fs%rho_g)
         call this%input%read('Surface tension coefficient',this%fs%sigma)
         ! Set acceleration of gravity
         call this%input%read('Gravity',this%fs%gravity)
         ! ! Inlets and coflow on the left
         ! call this%fs%add_bcond(name='inlets',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=left_boundary)
         ! ! Outflow on the right
         ! call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=20
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure velocity solver
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
      end block create_flow_solver
      

      ! Initialize our velocity field
      initialize_velocity: block
         use mpi_f08,    only: MPI_ALLREDUCE,MPI_SUM
         use parallel,   only: MPI_REAL_WP
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer  :: n,i,j,k,ierr
         ! real(WP) :: Ugas,myAgas,Agas,Uliq,myAliq,Aliq
         ! real(WP) :: Qgas,Qliq,myU
         ! real(WP), parameter :: SLPM2SI=1.66667E-5_WP
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! ! Read in gas flow rate and convert to SI
         ! call this%input%read('Gas flow rate (SLPM)',Qgas)
         ! Qgas=Qgas*SLPM2SI
         ! ! Calculate gas flow area - no overlap here!
         ! myAgas=0.0_WP
         ! call this%fs%get_bcond('gas_inlet',mybc)
         ! do n=1,mybc%itr%n_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    myAgas=myAgas+this%cfg%dy(j)*this%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
         ! end do
         ! call MPI_ALLREDUCE(myAgas,Agas,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! ! Calculate bulk gas velocity
         ! Ugas=Qgas/Agas
         ! ! Apply Dirichlet at gas inlet
         ! call this%fs%get_bcond('gas_inlet',mybc)
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    this%fs%U(i,j,k)=+sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*Ugas
         ! end do
         ! ! Read in liquid flow rate and convert to SI
         ! call this%input%read('Liquid flow rate (SLPM)',Qliq)
         ! Qliq=Qliq*SLPM2SI
         ! ! Calculate liquid flow area - no overlap here!
         ! myAliq=0.0_WP
         ! call this%fs%get_bcond('liq_inlet',mybc)
         ! do n=1,mybc%itr%n_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    myAliq=myAliq+this%cfg%dy(j)*this%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
         ! end do
         ! call MPI_ALLREDUCE(myAliq,Aliq,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! ! Calculate bulk axial velocity
         ! Uliq=Qliq/Aliq
         ! ! Apply Dirichlet at liquid injector port
         ! call this%fs%get_bcond('liq_inlet',mybc)
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    myU=2.0_WP*Uliq*(1.0_WP-min((this%fs%cfg%ym(j)**2+this%fs%cfg%zm(k)**2)/rl**2,1.0_WP))
         !    this%fs%U(i,j,k)=+sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*myU
         ! end do
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
      end block initialize_velocity


      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,np,nplane
         this%smesh=surfmesh(nvar=3,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='edge_sensor'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh_nowall(this%smesh)
         ! Calculate thickness even for plic
         if (.not.this%vf%two_planes) then
            allocate(this%vf%thickness(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%vf%thickness=0.0_WP
         end if
         call this%vf%get_thickness()
         ! Populate surface variables
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                        this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='atom')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_surface('plic',this%smesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_atom')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_atom')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
      end block create_monitor

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
   

      ! !> Function that localizes liquid stream at -x
      ! function liq_inlet(pg,i,j,k) result(isIn)
      !    use pgrid_class, only: pgrid
      !    class(pgrid), intent(in) :: pg
      !    integer, intent(in) :: i,j,k
      !    logical :: isIn
      !    real(WP) :: rad
      !    isIn=.false.
      !    rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      !    if (rad.lt.0.5_WP*dl.and.i.eq.pg%imin) isIn=.true.
      ! end function liq_inlet
      
      
      ! !> Function that localizes gas stream at -x
      ! function gas_inlet(pg,i,j,k) result(isIn)
      !    use pgrid_class, only: pgrid
      !    class(pgrid), intent(in) :: pg
      !    integer, intent(in) :: i,j,k
      !    logical :: isIn
      !    real(WP) :: rad
      !    isIn=.false.
      !    rad=sqrt(pg%ym(j)**2+pg%zm(k)**2)
      !    if (rad.ge.0.5_WP*dl.and.rad.lt.0.5_WP*dg.and.i.eq.pg%imin) isIn=.true.
      ! end function gas_inlet


      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer.or.&
         &   j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator


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
      
      
   end subroutine simulation_init
   

   !> Take one time step
   subroutine step(this)
      implicit none
      class(atom), intent(inout) :: this
      
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
         
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf)

      ! Turbulence modeling
      sgs_modeling: block
         use sgsmodel_class, only: vreman
         integer :: i,j,k
         this%resU=this%fs%rho_g
         call this%fs%get_gradu(this%gradU)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_
            do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_
               do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
                  this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
                  this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
                  this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
                  this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
               end do
            end do
         end do
      end block sgs_modeling
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
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
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         if (this%vf%two_planes) then
            call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         else
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         end if
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
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         real(WP) :: my_vof_removed
         my_vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            my_vof_removed=my_vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(my_vof_removed,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surface mesh
         update_smesh: block
            use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
            integer :: i,j,k,np,nplane
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh_nowall(this%smesh)
            ! Calculate thickness even for plic
            if (.not.this%vf%two_planes) call this%vf%get_thickness()
            ! Populate surface variables
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Write ensight files
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(atom), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%gradU)
      
   end subroutine final
   
   
end module atom_class