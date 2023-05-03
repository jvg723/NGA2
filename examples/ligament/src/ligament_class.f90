!> Definition for a ligament atomization class
module ligament_class
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use fene_class,        only: fene
   use ccl_class,         only: ccl
   use lpt_class,         only: lpt
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: ligament
   
   !> Ligament object
   type :: ligament
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf       !< Volume fraction solver
      type(tpns)        :: fs       !< Two-phase flow solver
      type(hypre_str)   :: ps       !< Structured Hypre linear solver for pressure
      type(ddadi)       :: vs,ss    !< DDADI solver for velocity and scalar
      type(fene)        :: nn       !< FENE model for polymer stress
      type(ccl)         :: cc       !< Connected component labeling class
      type(lpt)         :: lp       !< Particle tracking
      type(timetracker) :: time     !< Time info
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(partmesh) :: pmesh
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor)  :: mfile       !< General simulation monitoring
      type(monitor)  :: cflfile     !< CFL monitoring
      type(monitor)  :: scfile      !< Scalar monitoring
      type(monitor)  :: sprayfile   !< Spray monitoring file
      
      !> Work arrays
      real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      real(WP), dimension(:,:,:),     allocatable :: SRmag
      real(WP), dimension(:,:,:,:),   allocatable :: SR
      real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU

      !> Normalized viscosity array
      real(WP), dimension(:,:,:),   allocatable :: visc_norm
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      
      
   contains
      procedure :: init                    !< Initialize nozzle simulation
      procedure :: step                    !< Advance nozzle simulation by one time step
      procedure :: final                   !< Finalize nozzle simulation
      procedure :: transfer_vf_to_drops    !< Convert VOF to particles
   end type ligament
   
   
   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=5

   !> Shear thinning parameters
	real(WP) :: power_constant,alpha,visc_0,visc_inf,visc_l,visc_s,visc_p0

   !> Transfer model parameters
   real(WP) :: filmthickness_over_dx  =5.0e-1_WP
   real(WP) :: min_filmthickness      =1.0e-3_WP
   real(WP) :: diam_over_filmthickness=1.0e+1_WP
   real(WP) :: max_eccentricity       =5.0e-1_WP
   real(WP) :: d_threshold            =1.0e-1_WP
   

contains
   
   
   !> Function that defines a level set function for a droplet
   function levelset_droplet(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)
   end function levelset_droplet


   !> Function that defines a level set function for a ligament
   function levelset_ligament(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP-sqrt(xyz(1)**2+xyz(2)**2)
   end function levelset_ligament
   
   
   !> Initialization of ligament simulation
   subroutine init(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      
      ! Create the ligament mesh
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
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); call param_read('X ligament',xlig)
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
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=y,xper=.false.,yper=.true.,zper=.true.,name='Ligament')
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
         allocate(this%resU     (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV     (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW     (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui       (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi       (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi       (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SRmag    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR       (6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
			allocate(this%visc_norm(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resSC    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         allocate(this%SCtmp    (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,elvira,r2p,art,lvira
         use mms_geom,  only: cube_refine_vol
         use param,     only: param_read
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         this%vf=vfs(cfg=this%cfg,reconstruction_method=r2p,name='VOF')
         ! Initialize to a ligament
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_ligament,0.0_WP,amr_ref_lvl)
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
         use hypre_str_class, only: pcg_pfmg
         type(bcond), pointer :: mybc
         integer :: n,i,j,k      
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call param_read('Density ratio',this%fs%rho_l)
         call param_read('Reynolds number',this%fs%visc_g); this%fs%visc_g=1.0_WP/this%fs%visc_g
         call param_read('Viscosity ratio',visc_l); this%fs%visc_l=this%fs%visc_g*visc_l; visc_s=this%fs%visc_g*visc_l
         call param_read('Weber number',this%fs%sigma); this%fs%sigma=1.0_WP/this%fs%sigma
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         this%ps%maxlevel=20
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

      ! Create a FENE model 
      create_fene: block 
         use param, only: param_read
         use multiscalar_class, only: bquick
         use fene_class,        only: fenecr
         integer :: i,j,k
         ! Create FENE model solver
         this%nn=fene(cfg=this%cfg,model=fenecr,scheme=bquick,name='FENE')
         ! Assign unity density for simplicity
         this%nn%rho=1.0_WP
         ! Maximum extensibility of polymer chain
         call param_read('Maximum polymer extensibility',this%nn%Lmax)
         ! Polymer viscosity
         call param_read('Polymer viscosity',visc_p0)
         ! Configure implicit scalar solver
         this%ss=ddadi(cfg=this%cfg,name='scalar',nst=13)
         ! Setup the solver
         call this%nn%setup(implicit_solver=this%ss)
         ! Initialize conformation tensor to identity
         this%nn%SC(:,:,:,1)=1.0_WP !< Cxx
         this%nn%SC(:,:,:,4)=1.0_WP !< Cyy
         this%nn%SC(:,:,:,6)=1.0_WP !< Czz
      end block create_fene

      ! ! Set inital viscosity of liquid phase 
		! liquid_viscosity: block
      !    use param, only: param_read
      !    ! Carreau model parameters	
      !    call param_read('Power law constant',power_constant)
      !    call param_read('Shear rate parameter',alpha)
      !    ! !> Carreau model
      !    ! ! No infintie shear rate viscosity
      !    ! visc_inf=0.0_WP
      !    ! ! Zero shear rate viscosity
      !    ! visc_0=visc_l
      !    ! this%visc_norm=visc_0/visc_0
      !    ! Inital hybrid model viscosity
      !    visc_0=visc_s+visc_p0        !< Zero shear rate viscosity
      !    this%nn%visc=visc_p0         !< Polymer viscosity
      !    this%fs%visc_l=visc_0        !< Initial fluid viscosity
      !    this%visc_norm=visc_0/visc_0
      ! end block liquid_viscosity 

       ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%max_interface_planes=2
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%film_classify(Lbary=this%vf%Lbary,Gbary=this%vf%Gbary)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl

      ! Create a Lagrangian spray tracker
      create_lpt: block
         ! Create the solver
         this%lp=lpt(cfg=this%cfg,name='spray')
         ! Get particle density from the flow solver
         this%lp%rho=this%fs%rho_l
      end block create_lpt
   

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=8,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='type'
         this%smesh%varname(3)='x_velocity'
         this%smesh%varname(4)='y_velocity'
         this%smesh%varname(5)='z_velocity'
         this%smesh%varname(6)='visc_l'
         this%smesh%varname(7)='SRmag'
         this%smesh%varname(8)='norm_visc'
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
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%type(i,j,k)
                        this%smesh%var(3,np)=this%Ui(i,j,k)
                        this%smesh%var(4,np)=this%Vi(i,j,k)
                        this%smesh%var(5,np)=this%Wi(i,j,k)
                        this%smesh%var(6,np)=this%fs%visc_l(i,j,k)
                        this%smesh%var(7,np)=this%SRmag(i,j,k)
                        this%smesh%var(8,np)=this%visc_norm(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         ! Include an extra variable for droplet diameter
         this%pmesh=partmesh(nvar=1,nvec=0,name='lpt')
         this%pmesh%varname(1)='diameter'
         ! Transfer particles to pmesh
         call this%lp%update_partmesh(this%pmesh)
         ! Also populate diameter variable
         do i=1,this%lp%np_
            this%pmesh%var(1,i)=this%lp%p(i)%d
         end do
      end block create_pmesh


      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         integer :: nsc
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='ligament')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_surface('plic',this%smesh)
         call this%ens_out%add_scalar('SRmag', this%SRmag)
         call this%ens_out%add_scalar('visc_l',this%fs%visc_l)
         call this%ens_out%add_scalar('visc_vorm',this%visc_norm)
         do nsc=1,this%nn%nscalar
            call this%ens_out%add_scalar(trim(this%nn%SCname(nsc)),this%nn%SC(:,:,:,nsc))
         end do
         call this%ens_out%add_particle('spray',this%pmesh)
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
         call this%lp%get_max()
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
         this%scfile=monitor(this%nn%cfg%amRoot,'scalar')
         call this%scfile%add_column(this%time%n,'Timestep number')
         call this%scfile%add_column(this%time%t,'Time')
         do nsc=1,this%nn%nscalar
            call this%scfile%add_column(this%nn%SCmin(nsc),trim(this%nn%SCname(nsc))//'_min')
            call this%scfile%add_column(this%nn%SCmax(nsc),trim(this%nn%SCname(nsc))//'_max')
         end do
         call this%scfile%write()
         ! Create a spray monitor
         this%sprayfile=monitor(amroot=this%lp%cfg%amRoot,name='spray')
         call this%sprayfile%add_column(this%time%n,'Timestep number')
         call this%sprayfile%add_column(this%time%t,'Time')
         call this%sprayfile%add_column(this%time%dt,'Timestep size')
         call this%sprayfile%add_column(this%lp%np,'Droplet number')
         call this%sprayfile%add_column(this%lp%Umin, 'Umin')
         call this%sprayfile%add_column(this%lp%Umax, 'Umax')
         call this%sprayfile%add_column(this%lp%Umean,'Umean')
         call this%sprayfile%add_column(this%lp%Vmin, 'Vmin')
         call this%sprayfile%add_column(this%lp%Vmax, 'Vmax')
         call this%sprayfile%add_column(this%lp%Vmean,'Vmean')
         call this%sprayfile%add_column(this%lp%Wmin, 'Wmin')
         call this%sprayfile%add_column(this%lp%Wmax, 'Wmax')
         call this%sprayfile%add_column(this%lp%Wmean,'Wmean')
         call this%sprayfile%add_column(this%lp%dmin, 'dmin')
         call this%sprayfile%add_column(this%lp%dmax, 'dmax')
         call this%sprayfile%add_column(this%lp%dmean,'dmean')
         call this%sprayfile%write()
      end block create_monitor

      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! ! Calculate SR
		! call this%fs%get_strainrate(SR=this%SR)

      ! ! Model shear thinning viscosity
		! update_viscosity: block
      !    integer :: i,j,k
      !    do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_
      !      do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_
      !        do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
      !          ! Strain rate magnitude
      !          this%SRmag(i,j,k)=sqrt(2.00_WP*this%SR(1,i,j,k)**2+this%SR(2,i,j,k)**2+this%SR(3,i,j,k)**2+2.0_WP*(this%SR(4,i,j,k)**2+this%SR(5,i,j,k)**2+this%SR(6,i,j,k)**2))
      !          ! Carreau Model
      !          this%fs%visc_l(i,j,k)=visc_inf+(visc_0-visc_inf)*(1.00_WP+(alpha*this%SRmag(i,j,k))**2)**((power_constant-1.00_WP)/2.00_WP)
      !          this%visc_norm(i,j,k)=this%fs%visc_l(i,j,k)/visc_0
      !        end do
      !      end do
      !    end do
      !    call this%fs%cfg%sync(this%fs%visc_l)
      ! end block update_viscosity

      ! Advance our spray
      this%resU=this%fs%rho_g; this%resV=this%fs%visc_g
      call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W

      ! Remember old scalars
      this%nn%SCold=this%nn%SC
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
         
      ! VOF solver step
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)

      ! Calculate grad(U)
      call this%fs%get_gradU(this%gradU)
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)

         ! ! ============= SCALAR SOLVER =======================
            
         ! ! Reset interpolation metrics to QUICK scheme
         ! call this%nn%metric_reset()
            
         ! ! Build mid-time scalar
         ! this%nn%SC=0.5_WP*(this%nn%SC+this%nn%SCold)
         
         ! ! Explicit calculation of drhoSC/dt from scalar equation
         ! call this%nn%get_drhoSCdt(this%resSC,this%fs%Uold,this%fs%Vold,this%fs%Wold)
         
         ! ! Perform bquick procedure
         ! bquick: block
         !    integer :: i,j,k
         !    logical, dimension(:,:,:), allocatable :: flag
         !    ! Allocate work array
         !    allocate(flag(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    ! Assemble explicit residual
         !    this%resSC=-2.0_WP*(this%nn%SC-this%nn%SCold)+this%time%dt*this%resSC
         !    ! Apply it to get explicit scalar prediction
         !    this%SCtmp=2.0_WP*this%nn%SC-this%nn%SCold+this%resSC
         !    ! Check cells that require bquick
         !    do k=this%nn%cfg%kmino_,this%nn%cfg%kmaxo_
         !       do j=this%nn%cfg%jmino_,this%nn%cfg%jmaxo_
         !          do i=this%nn%cfg%imino_,this%nn%cfg%imaxo_
         !             if (this%SCtmp(i,j,k,1).le.0.0_WP.or.this%SCtmp(i,j,k,4).le.0.0_WP.or.this%SCtmp(i,j,k,6).le.0.0_WP.or.&
         !             &   this%SCtmp(i,j,k,1)+this%SCtmp(i,j,k,4)+this%SCtmp(i,j,k,6).ge.this%nn%Lmax**2) then
         !                flag(i,j,k)=.true.
         !             else
         !                flag(i,j,k)=.false.
         !             end if
         !          end do
         !       end do
         !    end do
         !    ! Adjust metrics
         !    call this%nn%metric_adjust(this%SCtmp,flag)
         !    ! Clean up
         !    deallocate(flag)
         !    ! Recompute drhoSC/dt
         !    call this%nn%get_drhoSCdt(this%resSC,this%fs%Uold,this%fs%Vold,this%fs%Wold)
         ! end block bquick
         
         ! ! Add fene sources
         ! call this%nn%addsrc_CgradU(this%gradU,this%resSC)
         ! call this%nn%addsrc_relax(this%resSC)
         
         ! ! Assemble explicit residual
         ! this%resSC=-2.0_WP*(this%nn%SC-this%nn%SCold)+this%time%dt*this%resSC
         
         ! ! Form implicit residual
         ! call this%nn%solve_implicit(this%time%dt,this%resSC,this%fs%Uold,this%fs%Vold,this%fs%Wold)
         
         ! ! Update scalars
         ! this%nn%SC=2.0_WP*this%nn%SC-this%nn%SCold+this%resSC
         
         ! ! Force the gas scalar to identity
         ! gas_scalar_forcing: block
         !    integer :: i,j,k
         !    do k=this%nn%cfg%kmino_,this%nn%cfg%kmaxo_
         !       do j=this%nn%cfg%jmino_,this%nn%cfg%jmaxo_
         !          do i=this%nn%cfg%imino_,this%nn%cfg%imaxo_
         !             if (this%nn%mask(i,j,k).eq.0) then
         !                this%nn%SC(i,j,k,1)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,1)+(1.0_WP-this%vf%VF(i,j,k))*1.0_WP
         !                this%nn%SC(i,j,k,2)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,2)
         !                this%nn%SC(i,j,k,3)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,3)
         !                this%nn%SC(i,j,k,4)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,4)+(1.0_WP-this%vf%VF(i,j,k))*1.0_WP
         !                this%nn%SC(i,j,k,5)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,5)
         !                this%nn%SC(i,j,k,6)=this%vf%VF(i,j,k)*this%nn%SC(i,j,k,6)+(1.0_WP-this%vf%VF(i,j,k))*1.0_WP
         !             end if
         !          end do
         !       end do
         !    end do
         ! end block gas_scalar_forcing

         ! ! Apply all other boundary conditions on the resulting field
         ! call this%nn%apply_bcond(this%time%t,this%time%dt)
         ! ! ===================================================
         
         ! ============ VELOCITY SOLVER ======================

         ! Calculate SR
         call this%fs%get_strainrate(SR=this%SR)

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
         !    integer :: i,j,k,n
         !    real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
         !    real(WP), dimension(:,:,:,:), allocatable :: stress
         !    ! Allocate work arrays
         !    allocate(stress(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
         !    allocate(Txy   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    allocate(Tyz   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    allocate(Tzx   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         !    ! Calculate the polymer relaxation
         !    stress=0.0_WP; call this%nn%addsrc_relax(stress)
         !    ! Build liquid stress tensor - rate depdenent viscosity (Carreau model from Ohta et al. 2019)
         !    do n=1,6
         !       this%SRmag=0.0_WP; this%nn%visc=visc_p0; this%fs%visc_l=visc_0
         !       do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
         !          do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
         !             do i=this%fs%cfg%imin_,this%fs%cfg%imax_
         !                ! Cell Strain rate magnitude
         !                this%SRmag(i,j,k)=sqrt(this%SR(1,i,j,k)**2+this%SR(2,i,j,k)**2+this%SR(3,i,j,k)**2+2.0_WP*(this%SR(4,i,j,k)**2+this%SR(5,i,j,k)**2+this%SR(6,i,j,k)**2))
         !                ! Rate depdentdent polymer viscostiy - polymer stress added to liquid viscosity
         !                this%nn%visc=(visc_0-visc_s)*(1.00_WP+(alpha*this%SRmag(i,j,k))**2)**((power_constant-1.00_WP)/2.00_WP)
         !                this%fs%visc_l(i,j,k)=visc_s+this%nn%visc
         !                ! Viscoelastic stress
         !                stress(i,j,k,n)=-this%nn%visc*this%vf%VF(i,j,k)*stress(i,j,k,n)
         !             end do
         !          end do
         !       end do
         !       ! Sync up stress
         !       call this%fs%cfg%sync(stress(:,:,:,n))
         !    end do
         !    ! Sync up total liquid viscosity
         !    call this%fs%cfg%sync(this%fs%visc_l)
         !    ! ! Build liquid stress tensor
         !    ! do n=1,6
         !    !    stress(:,:,:,n)=-nn%visc*vf%VF*stress(:,:,:,n)
         !    ! end do
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
         !             &                                                               +sum(this%fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
         !             &                                                               +sum(this%fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
         !             if (this%fs%vmask(i,j,k).eq.0) this%resV(i,j,k)=this%resV(i,j,k)+sum(this%fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
         !             &                                                               +sum(this%fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
         !             &                                                               +sum(this%fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
         !             if (this%fs%wmask(i,j,k).eq.0) this%resW(i,j,k)=this%resW(i,j,k)+sum(this%fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
         !             &                                                               +sum(this%fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &                  
         !             &                                                               +sum(this%fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
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
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
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

      ! Perform volume-fraction-to-droplet transfer
      call this%transfer_vf_to_drops()
      
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
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%type(i,j,k)
                           this%smesh%var(2,np)=this%vf%type(i,j,k)
                           this%smesh%var(3,np)=this%Ui(i,j,k)
                           this%smesh%var(4,np)=this%Vi(i,j,k)
                           this%smesh%var(5,np)=this%Wi(i,j,k)
                           this%smesh%var(6,np)=this%fs%visc_l(i,j,k)
                           this%smesh%var(7,np)=this%SRmag(i,j,k)
                           this%smesh%var(8,np)=this%visc_norm(i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Update partmesh object
         update_pmesh: block
            integer :: i
            ! Transfer particles to pmesh
            call this%lp%update_partmesh(this%pmesh)
            ! Also populate diameter variable
            do i=1,this%lp%np_
               this%pmesh%var(1,i)=this%lp%p(i)%d
            end do
         end block update_pmesh
         ! Perform ensight output
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%sprayfile%write() 
      
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi,this%SR,this%SRmag,this%resSC,this%gradU,this%SCtmp)
      
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

   !> Transfer vf to drops
   subroutine transfer_vf_to_drops(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Perform a first pass with simplest CCL
      call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      
      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam
         
         ! Loops over film segments contained locally
         do m=1,this%cc%n_meta_struct
            
            ! Test if sphericity is compatible with transfer
            lmin=this%cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=this%cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=this%cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
            if (eccentricity.gt.max_eccentricity) cycle
            
            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle
            
            ! Create drop from available liquid volume - only one root does that
            if (this%cc%cfg%amRoot) then
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(0,8)                                                                                 !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
               this%lp%p(np)%Acol=0.0_WP                                                                                   !< Give zero collision force (axial)
               this%lp%p(np)%Tcol=0.0_WP                                                                                   !< Give zero collision force (tangential)
               this%lp%p(np)%d   =diam                                                                                     !< Assign diameter to account for full volume
               this%lp%p(np)%pos =[this%cc%meta_structures_list(m)%x,this%cc%meta_structures_list(m)%y,this%cc%meta_structures_list(m)%z] !< Place the drop at the liquid barycenter
               this%lp%p(np)%vel =[this%cc%meta_structures_list(m)%u,this%cc%meta_structures_list(m)%v,this%cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])                !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                                                                        !< Activate it
               ! Increment particle counter
               this%lp%np_=np
            end if
            
            ! Find local structs with matching id
            do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
               if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(m)%id) cycle
               ! Remove liquid in meta-structure cells
               do l=1,this%cc%struct_list(this%cc%struct_map_(n))%nnode ! Loops over cells within local
                  i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
                  j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
                  k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
                  ! Remove liquid in that cell
                  this%vf%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      ! Perform more detailed CCL in a second pass
      this%cc%max_interface_planes=2
      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%cc%get_min_thickness()
      call this%cc%sort_by_thickness()
      
      ! Loop through identified films and remove those that are thin enough
      remove_film: block
         use mathtools, only: pi
         integer :: m,n,i,j,k,np,ip,np_old
         real(WP) :: Vt,Vl,Hl,Vd
         
         ! Loops over film segments contained locally
         do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
            
            ! Skip non-liquid films
            if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
            
            ! Skip films that are still thick enough
            if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.min_filmthickness) cycle
            
            ! We are still here: transfer the film to drops
            Vt=0.0_WP      ! Transferred volume
            Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            np_old=this%lp%np_  ! Remember old number of particles
            do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
               i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
               j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
               k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
               ! Increment liquid volume to remove
               Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               ! Estimate drop size based on local film thickness in current cell
               Hl=max(this%cc%film_thickness(i,j,k),min_filmthickness)
               Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
               ! Create drops from available liquid volume
               do while (Vl-Vd.gt.0.0_WP)
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  this%lp%p(np)%Acol=0.0_WP                                     !< Give zero collision force (axial)
                  this%lp%p(np)%Tcol=0.0_WP                                     !< Give zero collision force (tangential)
                  this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                     !< Place the drop at the liquid barycenter
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
                  this%lp%p(np)%flag=0                                          !< Activate it
                  ! Increment particle counter
                  this%lp%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
               end do
               ! Remove liquid in that cell
               this%vf%VF(i,j,k)=0.0_WP
            end do
            
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
               ! Add one last drop for remaining liquid volume
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(0,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%Acol=0.0_WP                                     !< Give zero collision force (axial)
               this%lp%p(np)%Tcol=0.0_WP                                     !< Give zero collision force (tangential)
               this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                     !< Place the drop at the liquid barycenter
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
            else ! Some particles were created, make them all larger
               do ip=np_old+1,this%lp%np_
                  this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
               end do
            end if
         end do
         
      end block remove_film
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      ! Resync the spray
      call this%lp%sync()
      
   end subroutine transfer_vf_to_drops
   
   
end module ligament_class