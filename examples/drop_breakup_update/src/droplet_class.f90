!> Definition for a droplet atomization class
module droplet_class
   use precision,            only: WP
   use config_class,         only: config
   use iterator_class,       only: iterator
   use ensight_class,        only: ensight
   use surfmesh_class,       only: surfmesh
   use hypre_str_class,      only: hypre_str
   use ddadi_class,          only: ddadi
   use vfs_class,            only: vfs
   use tpns_class,           only: tpns
   use tpviscoelastic_class, only: tpviscoelastic
   use timetracker_class,    only: timetracker
   use lpt_class,            only: lpt
   use transfermodel_class,  only: transfermodels
   use event_class,          only: event
   use monitor_class,        only: monitor
   implicit none
   private
   
   public :: droplet
   
   !> droplet object
   type :: droplet
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)            :: vf    !< Volume fraction solver
      type(tpns)           :: fs    !< Two-phase flow solver
      type(hypre_str)      :: ps    !< Structured Hypre linear solver for pressure
      type(ddadi)          :: vs    !< DDADI solver for velocity 
      type(tpviscoelastic) :: ve
      type(lpt)            :: lp
      type(timetracker)    :: time  !< Time info
      type(transfermodels) :: tm
      type(event)          :: ppevt
     

      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor)  :: mfile       !< General simulation monitoring
      type(monitor)  :: cflfile     !< CFL monitoring
      type(monitor)  :: scfile      !< Scalar monitoring
      type(monitor)  :: filmfile    !< Film monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
      real(WP), dimension(:,:,:,:,:), allocatable :: Atmp
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      
      
   contains
      procedure :: init                    !< Initialize nozzle simulation
      procedure :: step                    !< Advance nozzle simulation by one time step
      procedure :: final                   !< Finalize nozzle simulation
   end type droplet

   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=5

   integer :: id_largest_film
   real(WP) :: initial_volume

   !> Check for stabilization 
   logical :: stabilization 


contains
   
   
   !> Function that defines a level set function for a droplet
   function levelset_droplet(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
   end function levelset_droplet

   
   !> Initialization of droplet simulation
   subroutine init(this)
      implicit none
      class(droplet), intent(inout) :: this
      
      
      ! Create the droplet mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xlig
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); call param_read('X droplet',xlig)
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xlig
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='droplet')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_config
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU      (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV      (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW      (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui        (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi        (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi        (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resSC     (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         allocate(this%SCtmp     (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Atmp(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,elvira,r2p,flux_storage
         use mms_geom,  only: cube_refine_vol
         use param,     only: param_read
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2p,transport_method=flux_storage,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         ! Initialize to a droplet
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_droplet,0.0_WP,amr_ref_lvl)
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
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
      end block create_iterator

      
      ! Create a multiphase flow solver with bconds
      create_flow_solver: block
         use mathtools,       only: Pi
         use param,           only: param_read
         use tpns_class,      only: dirichlet,clipped_neumann,bcond
         use hypre_str_class, only: pcg_pfmg2
         type(bcond), pointer :: mybc
         integer :: n,i,j,k      
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call param_read('Density ratio',this%fs%rho_l)
         call param_read('Reynolds number',this%fs%visc_g); this%fs%visc_g=1.0_WP/this%fs%visc_g
         call param_read('Viscosity ratio',this%fs%visc_l); this%fs%visc_l=this%fs%visc_g*this%fs%visc_l
         call param_read('Weber number',this%fs%sigma); this%fs%sigma=1.0_WP/this%fs%sigma
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=1.0_WP
         end do
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
      end block create_flow_solver

      ! Create a viscoleastic model with log conformation stablization method
      create_viscoelastic: block
         use param,                only: param_read
         use tpviscoelastic_class, only: eptt,oldroydb
         integer :: i,j,k
         ! Create viscoelastic model solver
         call this%ve%init(cfg=this%cfg,phase=0,model=oldroydb,name='viscoelastic')
         ! Relaxation time for polymer
         call param_read('Polymer relaxation time',this%ve%trelax)
         ! Polymer viscosity
         call param_read('Polymer viscosity ratio',this%ve%visc_p);
         this%ve%visc_p=this%fs%visc_l*((1.00_WP-this%ve%visc_p)/this%ve%visc_p)
         ! Setup without an implicit solver
         call this%ve%setup()
         ! Check first if we use stabilization
         call param_read('Stabilization',stabilization,default=.false.)
         ! Initialize C scalar fields
         if (stabilization) then 
            !> Allocate storage fo eigenvalues and vectors
            allocate(this%ve%eigenval    (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%ve%eigenval=0.0_WP
            allocate(this%ve%eigenvec(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%ve%eigenvec=0.0_WP
            !> Allocate storage for reconstructured C and Cold
            allocate(this%ve%SCrec   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6)); this%ve%SCrec=0.0_WP
            allocate(this%ve%SCrecold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6)); this%ve%SCrecold=0.0_WP
            do k=this%cfg%kmino_,this%cfg%kmaxo_
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     if (this%vf%VF(i,j,k).gt.0.0_WP) then
                        this%ve%SCrec(i,j,k,1)=1.0_WP  !< Cxx
                        this%ve%SCrec(i,j,k,4)=1.0_WP  !< Cyy
                        this%ve%SCrec(i,j,k,6)=1.0_WP  !< Czz
                     end if
                  end do
               end do
            end do
            ! Get eigenvalues and eigenvectors
            call this%ve%get_eigensystem(this%vf%VF)
         else
            do k=this%cfg%kmino_,this%cfg%kmaxo_
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     if (this%vf%VF(i,j,k).gt.0.0_WP) then
                        this%ve%SC(i,j,k,1)=1.0_WP  !< Cxx
                        this%ve%SC(i,j,k,4)=1.0_WP  !< Cyy
                        this%ve%SC(i,j,k,6)=1.0_WP  !< Czz
                     end if
                  end do
               end do
            end do
         end if
         ! Apply boundary conditions
         call this%ve%apply_bcond(this%time%t,this%time%dt)
      end block create_viscoelastic

      ! Create a Lagrangian spray tracker
      create_lpt: block
         use param, only: param_read
         ! Create the solver
         this%lp=lpt(cfg=this%cfg,name='spray')
         ! Get particle density from the flow solver
         this%lp%rho=this%fs%rho_l
         ! Turn off drag
         this%lp%drag_model='none'         
         ! Initialize with zero particles
         call this%lp%resize(0)
         ! Get initial particle volume fraction
         call this%lp%update_VF()               
         ! Get particle statistics
         call this%lp%get_max()
      end block create_lpt

      ! Create a transfer model object
      create_filmmodel: block
         use param, only: param_read
         use mathtools, only: Pi
         initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3
         call this%tm%initialize(cfg=this%cfg,vf=this%vf,fs=this%fs,lp=this%lp)
      end block create_filmmodel
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=12,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='curv'
         this%smesh%varname(3)='edge_sensor'
         this%smesh%varname(4)='thin_sensor'
         this%smesh%varname(5)='thickness'
         this%smesh%varname(6)='trC'
         this%smesh%varname(7)='Cxx'
         this%smesh%varname(8)='Cxy'
         this%smesh%varname(9)='Cxz'
         this%smesh%varname(10)='Cyy'
         this%smesh%varname(11)='Cyz'
         this%smesh%varname(12)='Czz'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Also populate nplane variable
         this%smesh%var(1,:)=1.0_WP
         ! Initalize variables to 0
         this%smesh%var(2,:)=0.0_WP
         this%smesh%var(3,:)=0.0_WP
         this%smesh%var(4,:)=0.0_WP
         this%smesh%var(5,:)=0.0_WP
         this%smesh%var(6,:)=0.0_WP
         this%smesh%var(7,:)=0.0_WP
         this%smesh%var(8,:)=0.0_WP
         this%smesh%var(9,:)=0.0_WP
         this%smesh%var(10,:)=0.0_WP
         this%smesh%var(11,:)=0.0_WP
         this%smesh%var(12,:)=0.0_WP
         np=0
         if (stabilization) then
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                           this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                           this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                           this%smesh%var(6,np)=this%ve%SCrec(i,j,k,1)+this%ve%SCrec(i,j,k,4)+this%ve%SCrec(i,j,k,6)
                           this%smesh%var(7,np)=this%ve%SCrec(i,j,k,1)
                           this%smesh%var(8,np)=this%ve%SCrec(i,j,k,2)
                           this%smesh%var(9,np)=this%ve%SCrec(i,j,k,3)
                           this%smesh%var(10,np)=this%ve%SCrec(i,j,k,4)
                           this%smesh%var(11,np)=this%ve%SCrec(i,j,k,5)
                           this%smesh%var(12,np)=this%ve%SCrec(i,j,k,6)
                        end if
                     end do
                  end do
               end do
            end do
         else
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                           this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                           this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                           this%smesh%var(6,np)=this%ve%SC(i,j,k,1)+this%ve%SC(i,j,k,4)+this%ve%SC(i,j,k,6)
                           this%smesh%var(7,np)=this%ve%SC(i,j,k,1)
                           this%smesh%var(8,np)=this%ve%SC(i,j,k,2)
                           this%smesh%var(9,np)=this%ve%SC(i,j,k,3)
                           this%smesh%var(10,np)=this%ve%SC(i,j,k,4)
                           this%smesh%var(11,np)=this%ve%SC(i,j,k,5)
                           this%smesh%var(12,np)=this%ve%SC(i,j,k,6)
                        end if
                     end do
                  end do
               end do
            end do
         end if
      end block create_smesh

      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         integer :: nsc
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='droplet')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_scalar('thin_sensor',this%vf%thin_sensor)
         call this%ens_out%add_scalar('edge_sensor',this%vf%edge_sensor)
         call this%ens_out%add_vector('edge_normal',this%resU,this%resV,this%resW)
         call this%ens_out%add_surface('plic',this%smesh)
         do nsc=1,this%ve%nscalar
            call this%ens_out%add_scalar(trim(this%ve%SCname(nsc)),this%ve%SCrec(:,:,:,nsc))
            !call this%ens_out%add_scalar(trim(this%ve%SCname(nsc)),this%ve%SC(:,:,:,nsc))
         end do
         call this%ens_out%add_scalar('eigval1',this%ve%eigenval(1,:,:,:))
         call this%ens_out%add_scalar('eigval2',this%ve%eigenval(2,:,:,:))
         call this%ens_out%add_scalar('eigval3',this%ve%eigenval(3,:,:,:))
         call this%ens_out%add_scalar('eigvec11',this%ve%eigenvec(1,1,:,:,:))
         call this%ens_out%add_scalar('eigvec12',this%ve%eigenvec(1,2,:,:,:))
         call this%ens_out%add_scalar('eigvec13',this%ve%eigenvec(1,3,:,:,:))
         call this%ens_out%add_scalar('eigvec21',this%ve%eigenvec(2,1,:,:,:))
         call this%ens_out%add_scalar('eigvec22',this%ve%eigenvec(2,2,:,:,:))
         call this%ens_out%add_scalar('eigvec23',this%ve%eigenvec(2,3,:,:,:))
         call this%ens_out%add_scalar('eigvec31',this%ve%eigenvec(3,1,:,:,:))
         call this%ens_out%add_scalar('eigvec32',this%ve%eigenvec(3,2,:,:,:))
         call this%ens_out%add_scalar('eigvec33',this%ve%eigenvec(3,3,:,:,:))
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         integer :: nsc
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         if (stabilization) then
            call this%ve%get_max_reconstructed(this%vf%VF)
            call this%ve%get_max(this%vf%VF)
         else
            call this%ve%get_max(this%vf%VF)
         end if
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
         call this%mfile%add_column(this%vf%flotsam_error,'Flotsam error')
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
         ! Create scalar monitor
         this%scfile=monitor(this%ve%cfg%amRoot,'scalar')
         call this%scfile%add_column(this%time%n,'Timestep number')
         call this%scfile%add_column(this%time%t,'Time')
         if (stabilization) then
            do nsc=1,this%ve%nscalar
               call this%scfile%add_column(this%ve%SCrecmin(nsc),trim(this%ve%SCname(nsc))//'_RCmin')
               call this%scfile%add_column(this%ve%SCrecmax(nsc),trim(this%ve%SCname(nsc))//'_RCmax')
            end do
            do nsc=1,this%ve%nscalar
               call this%scfile%add_column(this%ve%SCmin(nsc),trim(this%ve%SCname(nsc))//'_lnmin')
               call this%scfile%add_column(this%ve%SCmax(nsc),trim(this%ve%SCname(nsc))//'_lnmax')
            end do
         else
            do nsc=1,this%ve%nscalar
               call this%scfile%add_column(this%ve%SCmin(nsc),trim(this%ve%SCname(nsc))//'_min')
               call this%scfile%add_column(this%ve%SCmax(nsc),trim(this%ve%SCname(nsc))//'_max')
            end do
         end if
         call this%scfile%write()
         ! Create film thickness monitor
         this%filmfile=monitor(amroot=this%fs%cfg%amRoot,name='film')
         call this%filmfile%add_column(this%time%n,'Timestep number')
         call this%filmfile%add_column(this%time%t,'Time')
         call this%filmfile%add_column(this%tm%min_thickness,'Min thickness')
         call this%filmfile%add_column(this%tm%film_volume,'Largest Film volume')
         call this%filmfile%add_column(this%tm%film_ratio,'Vb/V0')
         call this%filmfile%add_column(this%tm%converted_volume,'VOF-LPT Converted volume')
         call this%filmfile%add_column(this%tm%xL,'Droplet x-extent')
         call this%filmfile%add_column(this%tm%yL,'Droplet y-extent')
         call this%filmfile%add_column(this%tm%zL,'Droplet z-extent')
         call this%filmfile%add_column(this%tm%xbary,'Droplet x-position')
         call this%filmfile%add_column(this%tm%ybary,'Droplet y-position')
         call this%filmfile%add_column(this%tm%zbary,'Droplet z-position')
         call this%filmfile%add_column(this%tm%ubary,'Droplet x-velocity')
         call this%filmfile%add_column(this%tm%vbary,'Droplet y-velocity')
         call this%filmfile%add_column(this%tm%wbary,'Droplet z-velocity')
         call this%filmfile%write()
      end block create_monitor

   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc,harmonic_visc
      implicit none
      class(droplet), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! Calculate grad(U)
      call this%fs%get_gradU(this%gradU)

      ! Remember old reconstructed conformation tensor
      if (stabilization) this%ve%SCrecold=this%ve%SCrec

      ! Transport our liquid conformation tensor using log conformation
      advance_scalar: block
         integer :: i,j,k,nsc
         ! ! Add source terms for constitutive model
         ! if (stabilization) then 
         !    ! Streching 
         !    call this%ve%get_CgradU_log(this%gradU,this%SCtmp,this%vf%VFold); this%resSC=this%SCtmp
         !    ! Relxation
         !    ! call this%ve%get_relax_log(this%SCtmp,this%vf%VFold);             this%resSC=this%resSC+this%SCtmp
         ! else
         !    call this%ve%get_CgradU(this%gradU,this%SCtmp);    this%resSC=this%SCtmp
         !    call this%ve%get_relax(this%SCtmp,this%time%dt);   this%resSC=this%resSC+this%SCtmp
         ! end if
         ! this%ve%SC=this%ve%SC+this%time%dt*this%resSC
         ! call this%ve%apply_bcond(this%time%t,this%time%dt)
         ! this%ve%SCold=this%ve%SC
         ! Explicit calculation of dSC/dt from scalar equation
         call this%ve%get_dSCdt(dSCdt=this%resSC,U=this%fs%U,V=this%fs%V,W=this%fs%W,VFold=this%vf%VFold,VF=this%vf%VF,detailed_face_flux=this%vf%detailed_face_flux,dt=this%time%dt)
         ! Update our scalars
         do nsc=1,this%ve%nscalar
            where (this%ve%mask.eq.0.and.this%vf%VF.ne.0.0_WP) this%ve%SC(:,:,:,nsc)=(this%vf%VFold*this%ve%SCold(:,:,:,nsc)+this%time%dt*this%resSC(:,:,:,nsc))/this%vf%VF
            where (this%vf%VF.eq.0.0_WP) this%ve%SC(:,:,:,nsc)=0.0_WP
         end do
         ! Add source terms for constitutive model
         if (stabilization) then 
            ! Streching 
            call this%ve%get_CgradU_log(this%gradU,this%SCtmp,this%vf%VF); this%resSC=this%SCtmp
            ! Relxation
            ! call this%ve%get_relax_log(this%SCtmp,this%vf%VF);             this%resSC=this%resSC+this%SCtmp
         else
            call this%ve%get_CgradU(this%gradU,this%SCtmp);    this%resSC=this%SCtmp
            call this%ve%get_relax(this%SCtmp,this%time%dt);   this%resSC=this%resSC+this%SCtmp
         end if
         this%ve%SC=this%ve%SC+this%time%dt*this%resSC
         call this%ve%apply_bcond(this%time%t,this%time%dt)
         ! Apply boundary conditions
         ! call this%ve%apply_bcond(this%time%t,this%time%dt)
      end block advance_scalar

      if (stabilization) then 
         ! Get eigenvalues and eigenvectors
         call this%ve%get_eigensystem(this%vf%VF)
         ! Reconstruct conformation tensor
         call this%ve%reconstruct_conformation(this%vf%VF)
         ! Add in relaxtion source from semi-anlaytical integration
         call this%ve%get_relax_analytical(this%time%dt,this%vf%VF)
         ! Reconstruct lnC for next time step
         !> get eigenvalues and eigenvectors based on reconstructed C
         call this%ve%get_eigensystem_SCrec(this%vf%VF)
         !> Reconstruct lnC from eigenvalues and eigenvectors
         call this%ve%reconstruct_log_conformation(this%vf%VF)
         ! Take exp(eigenvalues) to use in next time-step
         this%ve%eigenval=exp(this%ve%eigenval)
      end if

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
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)

      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! ============ VELOCITY SOLVER ======================

         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)

         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)

         ! ! Add polymer stress term
         ! polymer_stress: block
         !    use tpviscoelastic_class, only: oldroydb,eptt
         !    integer :: i,j,k,nsc,n
         !    real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
         !    real(WP), dimension(:,:,:,:), allocatable :: stress
         !    real(WP) :: coeff,trace
         !    ! Allocate work arrays
         !    allocate(stress(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         !    allocate(Txy   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    allocate(Tyz   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    allocate(Tzx   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    ! Calculate polymer stress for a given model
         !    stress=0.0_WP
         !    if (stabilization) then !< Build stress tensor from reconstructed C
         !       select case (this%ve%model)
         !       case (oldroydb)
         !          coeff=this%ve%visc_p/this%ve%trelax
         !          do k=this%cfg%kmino_,this%cfg%kmaxo_
         !             do j=this%cfg%jmino_,this%cfg%jmaxo_
         !                do i=this%cfg%imino_,this%cfg%imaxo_
         !                   stress(i,j,k,1)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,1)-1.0_WP) !> xx tensor component
         !                   stress(i,j,k,2)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,2)-0.0_WP) !> xy tensor component
         !                   stress(i,j,k,3)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,3)-0.0_WP) !> xz tensor component
         !                   stress(i,j,k,4)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,4)-1.0_WP) !> yy tensor component
         !                   stress(i,j,k,5)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,5)-0.0_WP) !> yz tensor component
         !                   stress(i,j,k,6)=this%vf%VF(i,j,k)*coeff*(this%ve%SCrec(i,j,k,6)-1.0_WP) !> zz tensor component
         !                end do
         !             end do
         !          end do
         !       end select 
         !    else
         !       select case (this%ve%model)
         !       case (oldroydb)
         !          ! Calculate the polymer stress
         !          call this%ve%get_relax(stress,this%time%dt)
         !          ! Build liquid stress tensor
         !          do nsc=1,6
         !             stress(:,:,:,nsc)=-this%ve%visc_p*this%vf%VF*stress(:,:,:,nsc)
         !          end do
         !       end select
         !    end if
         !    ! Interpolate tensor components to cell edges
         !    do k=this%cfg%kmin_,this%cfg%kmax_+1
         !       do j=this%cfg%jmin_,this%cfg%jmax_+1
         !          do i=this%cfg%imin_,this%cfg%imax_+1
         !             Txy(i,j,k)=sum(this%fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
         !             Tyz(i,j,k)=sum(this%fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
         !             Tzx(i,j,k)=sum(this%fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
         !          end do
         !       end do
         !    end do
         !    ! Add divergence of stress to residual
         !    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         !       do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
         !          do i=this%fs%cfg%imin_,this%fs%cfg%imax_
         !             if (this%fs%umask(i,j,k).eq.0) this%resU(i,j,k)=this%resU(i,j,k)+sum(this%fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
         !             &                                               +sum(this%fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))                     &
         !             &                                               +sum(this%fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
         !             if (this%fs%vmask(i,j,k).eq.0) this%resV(i,j,k)=this%resV(i,j,k)+sum(this%fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
         !             &                                               +sum(this%fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))                &
         !             &                                               +sum(this%fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
         !             if (this%fs%wmask(i,j,k).eq.0) this%resW(i,j,k)=this%resW(i,j,k)+sum(this%fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
         !             &                                               +sum(this%fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))                     &                  
         !             &                                               +sum(this%fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
         !          end do
         !       end do
         !    end do
         !    ! Clean up
         !    deallocate(stress,Txy,Tyz,Tzx)
         ! end block polymer_stress
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals - if running implicitly 
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         ! Apply these residuals - if running explicitly
         ! this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho_U
         ! this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho_V
         ! this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho_W
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         !call this%fs%update_laplacian(pinpoint=[this%fs%cfg%imin,this%fs%cfg%jmin,this%fs%cfg%kmin])
         call this%fs%correct_mfr()
         call this%fs%get_div()
         !call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         call this%fs%add_surface_tension_jump_thin(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         !if (this%cfg%amRoot) this%fs%psolv%rhs(this%cfg%imin,this%cfg%jmin,this%cfg%kmin)=0.0_WP
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         integer :: n
         do n=1,this%vof_removal_layer%no_
            this%vf%VF(this%vof_removal_layer%map(1,n),this%vof_removal_layer%map(2,n),this%vof_removal_layer%map(3,n))=0.0_WP
         end do
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surfmesh object
         update_smesh: block
            use irl_fortran_interface
            integer :: i,j,k,nplane,np
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=1.0_WP
            ! Initalize variables to 0
            this%smesh%var(2,:)=0.0_WP
            this%smesh%var(3,:)=0.0_WP
            this%smesh%var(4,:)=0.0_WP
            this%smesh%var(5,:)=0.0_WP
            this%smesh%var(6,:)=0.0_WP
            this%smesh%var(7,:)=0.0_WP
            this%smesh%var(8,:)=0.0_WP
            this%smesh%var(9,:)=0.0_WP
            this%smesh%var(10,:)=0.0_WP
            this%smesh%var(11,:)=0.0_WP
            this%smesh%var(12,:)=0.0_WP
            np=0
            if (stabilization) then
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                              this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                              this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                              this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                              this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                              this%smesh%var(6,np)=this%ve%SCrec(i,j,k,1)+this%ve%SCrec(i,j,k,4)+this%ve%SCrec(i,j,k,6)
                              this%smesh%var(7,np)=this%ve%SCrec(i,j,k,1)
                              this%smesh%var(8,np)=this%ve%SCrec(i,j,k,2)
                              this%smesh%var(9,np)=this%ve%SCrec(i,j,k,3)
                              this%smesh%var(10,np)=this%ve%SCrec(i,j,k,4)
                              this%smesh%var(11,np)=this%ve%SCrec(i,j,k,5)
                              this%smesh%var(12,np)=this%ve%SCrec(i,j,k,6)
                           end if
                        end do
                     end do
                  end do
               end do
            else
               do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
                  do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                     do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                              this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                              this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                              this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                              this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                              this%smesh%var(6,np)=this%ve%SC(i,j,k,1)+this%ve%SC(i,j,k,4)+this%ve%SC(i,j,k,6)
                              this%smesh%var(7,np)=this%ve%SC(i,j,k,1)
                              this%smesh%var(8,np)=this%ve%SC(i,j,k,2)
                              this%smesh%var(9,np)=this%ve%SC(i,j,k,3)
                              this%smesh%var(10,np)=this%ve%SC(i,j,k,4)
                              this%smesh%var(11,np)=this%ve%SC(i,j,k,5)
                              this%smesh%var(12,np)=this%ve%SC(i,j,k,6)
                           end if
                        end do
                     end do
                  end do
               end do
            end if
         end block update_smesh
         ! Transfer edge normal data
         this%resU=this%vf%edge_normal(1,:,:,:)
         this%resV=this%vf%edge_normal(2,:,:,:)
         this%resW=this%vf%edge_normal(3,:,:,:)
         ! Perform ensight output
         call this%ens_out%write_data(this%time%t)
      end if

      ! Calculate film stats
      call this%tm%transfer_vf_to_drops()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      if (stabilization) then
         call this%ve%get_max_reconstructed(this%vf%VF)
         call this%ve%get_max(this%vf%VF)
      else
         call this%ve%get_max(this%vf%VF)
      end if
      call this%mfile%write()
      call this%cflfile%write()
      call this%scfile%write()
      call this%filmfile%write()
      
   end subroutine step
   

   !> Finalize droplet simulation
   subroutine final(this)
      implicit none
      class(droplet), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      deallocate(this%resSC,this%gradU,this%SCtmp)
      
   end subroutine final
   
   
   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator
   
   
   !> Function that localizes region of VOF removal
   function vof_removal_layer_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.ge.pg%imax-nlayer) isIn=.true.
   end function vof_removal_layer_locator   
   
end module droplet_class