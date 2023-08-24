!> Definition for a ligament atomization class
module ligament_class
   use string,            only: str_medium
   use precision,         only: WP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use ccl_class,         only: ccl
   use lpt_class,         only: lpt
   use transfermodel_class,only: transfermodels
   implicit none
   private
   
   public :: ligament
   
   !> Ligament object
   type :: ligament
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      !type(ddadi)       :: vs    !< DDADI solver for velocity
      type(timetracker) :: time  !< Time info
      type(ccl)         :: cc
      type(lpt)         :: lp
      type(transfermodels) :: tm

      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(partmesh) :: pmesh
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      
      !> Transfer model parameters
      real(WP) :: filmthickness_over_dx  =5.0e-1_WP
      real(WP) :: min_filmthickness      =1.0e-4_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP
      real(WP) :: max_eccentricity       =5.0e-1_WP
      real(WP) :: d_threshold            =1.0e-1_WP
      logical  :: debug

   contains
      procedure :: init                            !< Initialize nozzle simulation
      procedure :: step                            !< Advance nozzle simulation by one time step
      procedure :: final                           !< Finalize nozzle simulation
      procedure :: bag_droplet_gamma
      procedure :: transfer_vf_to_drops
   end type ligament
   
   
   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=5
   

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
      
      ! ! this%debug flag
      ! activate_debug: block
      !    use param, only: param_read
      !    call param_read('Print debug',this%debug,default=.false.)   
      ! end block activate_debug
      
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
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: VFlo,VFhi,elvira,r2p
         use mms_geom,  only: cube_refine_vol
         use param,     only: param_read
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         !this%vf=vfs(cfg=this%cfg,reconstruction_method=r2p,name='VOF')
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=r2p,name='VOF')
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
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
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
      

      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         integer :: i,j,k
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%max_interface_planes=2
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl


      ! Create a Lagrangian spray tracker
      create_lpt: block
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
         ! Handle restarts
         ! if (restarted) call this%lp%read(filename=trim(lpt_file))
      end block create_lpt


      ! Create a transfer model object
      create_transfermodel: block
         call this%tm%initialize(cfg=this%cfg,vf=this%vf,fs=this%fs,lp=this%lp)
      end block create_transfermodel


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         this%smesh=surfmesh(nvar=5,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='curv'
         this%smesh%varname(3)='edge_sensor'
         this%smesh%varname(4)='thin_sensor'
         this%smesh%varname(5)='thickness'
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
                        this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                        this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                        this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                        this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
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
         call this%ens_out%add_scalar('thin_sensor',this%vf%thin_sensor)
         call this%ens_out%add_scalar('edge_sensor',this%vf%edge_sensor)
         call this%ens_out%add_vector('edge_normal',this%resU,this%resV,this%resW)
         call this%ens_out%add_surface('plic',this%smesh)
         call this%ens_out%add_particle('spray',this%pmesh)
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
         call this%mfile%add_column(this%lp%np,'Particle number')
         call this%mfile%add_column(this%vf%VFmax,'VOF maximum')
         call this%mfile%add_column(this%vf%VFmin,'VOF minimum')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vf%flotsam_error,'Flotsam error')
         call this%mfile%add_column(this%vf%thinstruct_error,'Film error')
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
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc,harmonic_visc
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! Advance our spray
      this%resU=this%fs%rho_g; this%resV=this%fs%visc_g
      call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)      
   
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
         !call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU/this%fs%rho_U
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV/this%fs%rho_V
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW/this%fs%rho_W
         
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
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Perform volume-fraction-to-droplet transfer
      call this%tm%transfer_vf_to_drops()

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
                           this%smesh%var(2,np)=this%vf%curv2p(nplane,i,j,k)
                           this%smesh%var(3,np)=this%vf%edge_sensor(i,j,k)
                           this%smesh%var(4,np)=this%vf%thin_sensor(i,j,k)
                           this%smesh%var(5,np)=this%vf%thickness  (i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Transfer edge normal data
         this%resU=this%vf%edge_normal(1,:,:,:)
         this%resV=this%vf%edge_normal(2,:,:,:)
         this%resW=this%vf%edge_normal(3,:,:,:)
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
      
      
   end subroutine step
   

   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ligament), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
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
      use vfs_class, only: r2p
      use messager, only: die
      implicit none
      class(ligament), intent(inout) :: this
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,id_largest_film

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
            eccentricity=sqrt(1.0_WP-lmin**2/(lmax**2+tiny(1.0_WP)))
            if (eccentricity.gt.this%max_eccentricity) cycle

            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if ((diam.eq.0.0_WP).or.(diam.gt.this%d_threshold)) cycle

            
            ! Create drop from available liquid volume - only one root does that
            if (this%cc%cfg%amRoot) then
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(1,8)                                                                                 !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
               this%lp%p(np)%Acol =0.0_WP                                                                                   !< Give zero collision force
               this%lp%p(np)%Tcol =0.0_WP                                                                                   !< Give zero collision force
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
      if (this%debug) print *,'ts',this%time%n,'rank',this%cfg%rank,'done remove struct'

      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      if (this%vf%reconstruction_method.eq.r2p) then
         ! Perform more detailed CCL in a second pass
         this%cc%max_interface_planes=2
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%vf%sense_interface()
         call this%cc%get_min_thickness()
         call this%cc%sort_by_thickness()

         if (this%debug) print *,'ts',this%time%n,'rank',this%cfg%rank,'remove film'

         ! Loop through identified films and remove those that are thin enough
         remove_film: block
            use mathtools, only: pi,normalize,cross_product
            use random,    only: random_uniform
            use myrandom,  only: random_gamma
            use parallel,  only: MPI_REAL_WP
            use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_LOGICAL,MPI_LOR
            use irl_fortran_interface
            integer :: m,n,l,i,j,k,np,ip,np_old,ii,jj,kk
            real(WP) :: Vt,Vl,Hl,Vd,diam,d0
            real(WP), dimension(3) :: pref,nref,tref,sref
            real(WP) :: theta
            logical :: has_burst
            real(WP), dimension(2) :: s_bary,c_bary
            integer :: nodes_transferred
            ! Conversion model
            real(WP) :: alpha, beta, curv_sum, ncurv
            real(WP) :: max_film_droplet_diam=4.0e-4_WP

            ! Loops over film segments contained locally
            d0=1.0_WP ! Characteristic diameter

            do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
            
               ! Skip non-liquid films
               if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
            
               ! Skip films that are still thick enough
               if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.this%min_filmthickness) cycle

               ! We are still here: transfer the film to drops
               print *,'ts',this%time%n,'rank',this%cfg%rank,'film m',m,'start hole at min thickness'
               i=this%cc%film_list(this%cc%film_map_(m))%node(1,1)
               j=this%cc%film_list(this%cc%film_map_(m))%node(2,1)
               k=this%cc%film_list(this%cc%film_map_(m))%node(3,1)
               ! Compute cell-averaged curvature for cell containing thinnest film segment
               curv_sum=0.0_WP; ncurv=0.0_WP
               do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                     curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                     ncurv=ncurv+1.0_WP
                  end if
               end do
               ! Determine gamma distribution parameters                        
               call this%bag_droplet_gamma(this%cc%film_thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
               Vt=0.0_WP      ! Transferred volume
               Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
               np_old=this%lp%np_  ! Remember old number of particles
               Vd=pi/6.0_WP*(random_gamma(alpha,.true.)*beta*d0)**3
               do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
                  i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
                  j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
                  k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
                  ! Get local coordinate system
                  nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
                  select case (maxloc(abs(nref),1))
                  case (1)
                     tref=normalize([+nref(2),-nref(1),0.0_WP])
                  case (2)
                     tref=normalize([0.0_WP,+nref(3),-nref(2)])
                  case (3)
                     tref=normalize([-nref(3),0.0_WP,+nref(1)])
                  end select
                  sref=cross_product(nref,tref)
                  ! Increment liquid volume to remove
                  Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)

                  ! Compute cell-averaged curvature
                  curv_sum=0.0_WP; ncurv=0.0_WP
                  do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                        curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                        ncurv=ncurv+1.0_WP
                     end if
                  end do
                  ! Determine gamma distribution parameters
                  call this%bag_droplet_gamma(this%cc%film_thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
                  
                  ! Create drops from available liquid volume
                  do while (Vl-Vd.gt.0.0_WP)
                     ! Make room for new drop
                     np=this%lp%np_+1; call this%lp%resize(np)
                     ! Add the drop
                     this%lp%p(np)%id  =int(1,8)                                   !< Give id (maybe based on break-up model?)
                     this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                     this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
                     this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
                     this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
                     ! theta        =random_uniform(0.0_WP,2.0_WP*pi)           !< Angular position of randomly placed droplet relative to liquid barycenter
                     ! this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(0.0_WP,0.5_WP*this%vf%cfg%meshsize(i,j,k))*(cos(theta)*tref+sin(theta)*sref)
                     this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref              
                     this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                     this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
                     this%lp%p(np)%flag=0                                          !< Activate it
                     ! Increment particle counter
                     this%lp%np_=np
                     ! Update tracked volumes
                     Vl=Vl-Vd
                     Vt=Vt+Vd
                     ! Output diameter, velocity, and position
                     ! Generate new droplet volume
                     Vd=pi/6.0_WP*(random_gamma(alpha,.true.)*beta*d0)**3
                  end do

                  ! Remove liquid in that cell
                  this%vf%VF(i,j,k)=0.0_WP

               end do
               
               ! Based on how many particles were created, decide what to do with left-over volume
               if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                  ! Add one last drop for remaining liquid volume
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(1,8)                                   !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
                  this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
                  this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
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

            end do ! local films
            
            ! Sync VF and clean up IRL and band
            call this%vf%cfg%sync(this%vf%VF)
            call this%vf%clean_irl_and_band()                  

            ! Clean up CCL
            call this%cc%deallocate_lists()

         end block remove_film

      end if ! r2p
     
      ! Resync the spray
      call this%lp%sync()

   end subroutine transfer_vf_to_drops
   

   !> Generate a Gamma distribution for bag droplet formation
   !> where the number PDF is in the form
   !> p_n(x=d;alpha,beta)=x**(alpha-1)*exp(-x/beta)/beta**alpha/gamma(alpha)
   !> Adapted from Jackiw and Ashgriz 2022, JFM
   subroutine bag_droplet_gamma(this,h,R,alpha,beta)
      use myrandom, only: random_gamma
      implicit none
      class(ligament), intent(inout) :: this
      real(WP), intent(in) :: h,R
      real(WP), intent(out) :: alpha,beta
      real(WP) :: d0,Utc,ac,b,lambda_rr,dr,ds,Oh
      real(WP) :: mean, stdev
      
      ! assert h,R != 0
      ! Drop diameter
      d0=1.0_WP
      ! Retraction speed
      Utc=sqrt(2.0_WP*this%fs%sigma/this%fs%rho_l/h)
      ! Centripetal acceleration
      ac=Utc**2/R
      ! Rim diameter
      b=sqrt(this%fs%sigma/this%fs%rho_l/ac)
      ! Receding rim wavelength
      lambda_rr=4.5_WP*b
      ! RP droplet diameter
      dr=1.89_WP*b
      ! Rim Ohnesorge number
      Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b**3*this%fs%sigma)
      ! Satellite droplet diameter
      ds=dr/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
      ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
      mean=0.25_WP*(h+b+dr+ds)/d0
      stdev=sqrt(0.25_WP*sum(([h,b,dr,ds]/d0-mean)**2))
      ! Gamma distribution parameters
      alpha=(mean/stdev)**2
      beta=stdev**2/mean

   end subroutine bag_droplet_gamma


end module ligament_class
