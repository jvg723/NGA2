!> Definition for block 2 simulation: swirl atomizer
module block2_class
   use precision,         only: WP
   use geometry,          only: Dout,Din
   use config_class,      only: config
   use iterator_class,    only: iterator
   use hypre_str_class,   only: hypre_str
   use hypre_uns_class,   only: hypre_uns
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use fene_class,        only: fene
   use ccl_class,         only: ccl
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
      class(config), pointer :: cfg                    !< Pointer to config
      ! type(hypre_str)        :: ps                     !< Structured hypre pressure solver
	   type(hypre_uns)        :: ps                     !< Unstructured hypre pressure solver
      type(ddadi)            :: vs,ss                  !< DDAI implicit and scalar solver
      type(tpns)             :: fs                     !< Two phase incompressible flow solver
      type(vfs)              :: vf                     !< VF solver
      type(fene)             :: nn                     !< FENE model for non-newt behavior
      type(ccl)              :: cc                     !< Connected component labeling class
      type(timetracker)      :: time                   !< Time tracker
      type(surfmesh)         :: smesh                  !< Surfmesh                                       
      type(ensight)          :: ens_out                !< Ensight output
      type(event)            :: ens_evt                !< Ensight event
      type(monitor)          :: mfile,cflfile,scfile   !< Monitor files
      type(event)            :: ppevt                  !< Event for post-processing
      !> Private work arrays
      real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:),     allocatable :: SRmag
      real(WP), dimension(:,:,:,:),   allocatable :: SR
      real(WP), dimension(:,:,:),     allocatable :: visc_norm
      real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2

   !> Eccentricity threshold
   real(WP) :: max_eccentricity=2.0_WP

   !> Shear thinning parameters
	real(WP) :: power_constant,alpha,visc_0,visc_inf,visc_l,beta

   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=4

   !> Structure post processing 
   real(WP) :: x_over_xL   !< Parameter to set location of collection along x-axis
   real(WP) :: volume,x_cg,y_cg,z_cg,U_vel,V_vel,W_vel,length_1,length_2,length_3,eccent
   real(WP), dimension(:,:), allocatable :: mother_dropsize,satellite_dropsize,wave_number !< Arrays for K&M 2020 breakup model

contains
   
   
   !> Function that defines a level set function for an annular liquid region
	function levelset_annulus(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G,r
      r=sqrt(xyz(2)**2+xyz(3)**2)
	   G=min(r-0.5_WP*Din,0.5_WP*Dout-r)
	end function levelset_annulus
   

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

   !> Function that localizes region of VOF removal
   function vof_removal_layer_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.ge.pg%imax-nlayer.or.&
      &   j.le.pg%jmin+nlayer.or.&
      &   j.ge.pg%jmax-nlayer.or.&
      &   k.le.pg%kmin+nlayer.or.&
      &   k.ge.pg%kmax-nlayer) isIn=.true.
   end function vof_removal_layer_locator
   
   
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
         allocate(b%SRmag(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         allocate(b%SR (6,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         ! ! Array for normalized viscosity
			allocate(b%visc_norm(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
         ! Scalar solver arrays
         allocate(b%resSC (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_,1:6))
         allocate(b%SCtmp (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_,1:6))
         allocate(b%gradU (1:3,1:3,b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         b%time=timetracker(amRoot=b%cfg%amRoot,name='swirl_atomizer')
         call param_read('2 Max timestep size',b%time%dtmax)
         call param_read('Max cfl number',b%time%cflmax)
         call param_read('Max time',b%time%tmax)
         b%time%dt=b%time%dtmax
         b%time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,r2p,VFhi,VFlo
         integer :: i,j,k,n,si,sj,sk
			real(WP), dimension(3,8) :: cube_vertex
			real(WP), dimension(3) :: v_cent,a_cent
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         b%vf=vfs(cfg=b%cfg,reconstruction_method=lvira,name='VOF')
         ! Set full domain to gas
         do k=b%vf%cfg%kmino_,b%vf%cfg%kmaxo_
            do j=b%vf%cfg%jmino_,b%vf%cfg%jmaxo_
               do i=b%vf%cfg%imino_,b%vf%cfg%imaxo_
                  b%vf%VF(i,j,k)=0.0_WP
                  b%vf%Lbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
                  b%vf%Gbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Initialize an annular interface at the inlet only
         if (b%vf%cfg%iproc.eq.1) then
            do k=b%vf%cfg%kmino_,b%vf%cfg%kmaxo_
               do j=b%vf%cfg%jmino_,b%vf%cfg%jmaxo_
                  do i=b%vf%cfg%imino,b%vf%cfg%imin-1
                     ! Set cube vertices
				         n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1; cube_vertex(:,n)=[b%vf%cfg%x(i+si),b%vf%cfg%y(j+sj),b%vf%cfg%z(k+sk)]
                           end do
                        end do
                     end do
                     ! Call adaptive refinement code to get volume and barycenters recursively
				         vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_annulus,0.0_WP,amr_ref_lvl)
                     b%vf%VF(i,j,k)=vol/b%vf%cfg%vol(i,j,k)
                     if (b%vf%VF(i,j,k).ge.VFlo.and.b%vf%VF(i,j,k).le.VFhi) then
                        b%vf%Lbary(:,i,j,k)=v_cent
                        b%vf%Gbary(:,i,j,k)=([b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]-b%vf%VF(i,j,k)*b%vf%Lbary(:,i,j,k))/(1.0_WP-b%vf%VF(i,j,k))
                     else
                        b%vf%Lbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
                        b%vf%Gbary(:,i,j,k)=[b%vf%cfg%xm(i),b%vf%cfg%ym(j),b%vf%cfg%zm(k)]
                     end if
                  end do
               end do
            end do
         end if
         ! Update the band
         call b%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call b%vf%build_interface()
         ! Create discontinuous polygon mesh from IRL interface
         call b%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call b%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call b%vf%subcell_vol()
         ! Calculate curvature
         call b%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call b%vf%reset_volume_moments()
      end block create_and_initialize_vof

      ! Create an iterator for removing VOF at edges
      create_iterator: block
         b%vof_removal_layer=iterator(b%cfg,'VOF removal',vof_removal_layer_locator)
         b%vof_removed=0.0_WP
      end block create_iterator

      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,slip
         ! use hypre_str_class, only: pcg_pfmg,gmres_pfmg
         use hypre_uns_class,  only: gmres_amg,pcg_amg
         integer :: i,j,k
         ! Create flow solver
         b%fs=tpns(cfg=b%cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',visc_l) !; b%fs%visc_l=visc_l
         call param_read('Gas dynamic viscosity'   ,b%fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',b%fs%rho_l)
         call param_read('Gas density'   ,b%fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',b%fs%sigma)
         ! Inflow on the left of domain
         call b%fs%add_bcond(name='inflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Clipped Neumann outflow on the right of domain
         call b%fs%add_bcond(name='bc_xp' ,type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=xp_locator)
         ! Slip on the sides
         call b%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call b%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call b%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call b%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         b%ps=hypre_uns(cfg=b%cfg,name='Pressure',method=pcg_amg,nst=7)
			! b%ps=hypre_str(cfg=b%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ! b%ps%maxlevel=10
         ! b%ps=hypre_uns(cfg=b%cfg,name='Pressure',method=gmres_amg,nst=7)
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Configure velocity solver
			b%vs=ddadi(cfg=b%cfg,name='Velocity',nst=7)
			! Setup the solver
			call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%vs)
      end block create_and_initialize_flow_solver

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
         call param_read('Polymer relaxation time',b%nn%lambda)
         ! ! Polymer viscosity
         ! call param_read('Polymer viscosity',nn%visc)
         ! Configure implicit scalar solver
         b%ss=ddadi(cfg=b%cfg,name='scalar',nst=13)
         ! Setup the solver
         call b%nn%setup(implicit_solver=b%ss)
         ! Initialize conformation tensor to identity
         b%nn%SC(:,:,:,1)=1.0_WP !< Cxx
         b%nn%SC(:,:,:,4)=1.0_WP !< Cyy
         b%nn%SC(:,:,:,6)=1.0_WP !< Czz
      end block create_fene

      ! Set inital viscosity of liquid phase 
		liquid_viscosity: block
         ! Carreau model parameters	
         call param_read('Power law constant',power_constant)
         call param_read('Shear rate parameter',alpha)
         call param_read('Viscosity Ratio',beta)
         ! call param_read('Zero shear rate viscosity',visc_0)
         ! !> Carreau model
         ! ! No infintie shear rate viscosity
         ! visc_inf=0.0_WP
         ! ! Zero shear rate viscosity
         ! visc_0=visc_l; b%fs%visc_l=visc_0
         ! b%visc_norm=b%fs%visc_l/visc_0
         !> Hybrid model
         ! Zero/initial shear rate polymer viscosity
         b%nn%visc=visc_l/beta-visc_l
         ! Zero shear rate liquid viscosity
         visc_0=visc_l+b%nn%visc
         b%fs%visc_l=visc_0
      end block liquid_viscosity 
      

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

      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         b%cc=ccl(cfg=b%cfg,name='CCL')
         b%cc%max_interface_planes=2
         b%cc%VFlo=VFlo
         b%cc%dot_threshold=-0.5_WP
         ! Perform CCL step
         call b%cc%build_lists(VF=b%vf%VF,poly=b%vf%interface_polygon,U=b%fs%U,V=b%fs%V,W=b%fs%W)
         call b%cc%film_classify(Lbary=b%vf%Lbary,Gbary=b%vf%Gbary)
         call b%cc%deallocate_lists()
      end block create_and_initialize_ccl
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         b%smesh=surfmesh(nvar=7,name='plic')
         b%smesh%varname(1)='nplane'
         b%smesh%varname(2)='x_velocity'
         b%smesh%varname(3)='y_velocity'
         b%smesh%varname(4)='z_velocity'
         b%smesh%varname(5)='visc_l'
         b%smesh%varname(6)='SRmag'
         b%smesh%varname(7)='norm_viscl'
         ! Transfer polygons to smesh
         call b%vf%update_surfmesh(b%smesh)
         ! Also populate nplane variable
         b%smesh%var(1,:)=1.0_WP
         ! Initalize variables to 0
         b%smesh%var(2,:)=0.0_WP
         b%smesh%var(3,:)=0.0_WP
         b%smesh%var(4,:)=0.0_WP
         b%smesh%var(5,:)=0.0_WP
         b%smesh%var(6,:)=0.0_WP
         b%smesh%var(7,:)=0.0_WP
         np=0
         do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
            do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
               do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1
                        b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                        b%smesh%var(2,np)=b%Ui(i,j,k)
                        b%smesh%var(3,np)=b%Vi(i,j,k)
                        b%smesh%var(4,np)=b%Wi(i,j,k)
                        b%smesh%var(5,np)=b%fs%visc_l(i,j,k)
                        b%smesh%var(6,np)=b%SRmag(i,j,k)
                        b%smesh%var(7,np)=b%visc_norm(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         ! Create Ensight output from cfg
         b%ens_out=ensight(cfg=b%cfg,name='swirl_atomizer')
         ! Create event for Ensight output
         b%ens_evt=event(time=b%time,name='Ensight output')
         call param_read('Ensight output period',b%ens_evt%tper)
         ! Add variables to output
         call b%ens_out%add_vector('velocity',b%Ui,b%Vi,b%Wi)
         call b%ens_out%add_scalar('VOF',b%vf%VF)
         call b%ens_out%add_scalar('Pressure',b%fs%P)
         call b%ens_out%add_scalar('curvature',b%vf%curv)
         call b%ens_out%add_surface('vofplic',b%smesh)
         call b%ens_out%add_scalar('SRmag',b%SRmag)
         call b%ens_out%add_scalar('visc_l',b%fs%visc_l)
         ! call b%ens_out%add_scalar('normvisc_l',b%visc_norm)
         do nsc=1,b%nn%nscalar
            call b%ens_out%add_scalar(trim(b%nn%SCname(nsc)),b%nn%SC(:,:,:,nsc))
         end do
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc   
         ! Prepare some info about fields
         call b%fs%get_cfl(b%time%dt,b%time%cfl)
         call b%fs%get_max()
         call b%vf%get_max()
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
         call b%mfile%add_column(b%vf%VFmax,'VOF maximum')
         call b%mfile%add_column(b%vf%VFmin,'VOF minimum')
         call b%mfile%add_column(b%vf%VFint,'VOF integral')
         call b%mfile%add_column(b%vf%SDint,'SD integral')
         call b%mfile%add_column(b%fs%divmax,'Maximum divergence')
         call b%mfile%add_column(b%fs%psolv%it,'Pressure iteration')
         call b%mfile%add_column(b%fs%psolv%rerr,'Pressure error')
         call b%mfile%write()
         ! Create CFL monitor
         b%cflfile=monitor(b%fs%cfg%amRoot,'cfl2')
         call b%cflfile%add_column(b%time%n,'Timestep number')
         call b%cflfile%add_column(b%time%t,'Time')
         call b%cflfile%add_column(b%fs%CFLst,'STension CFL')
         call b%cflfile%add_column(b%fs%CFLc_x,'Convective xCFL')
         call b%cflfile%add_column(b%fs%CFLc_y,'Convective yCFL')
         call b%cflfile%add_column(b%fs%CFLc_z,'Convective zCFL')
         call b%cflfile%add_column(b%fs%CFLv_x,'Viscous xCFL')
         call b%cflfile%add_column(b%fs%CFLv_y,'Viscous yCFL')
         call b%cflfile%add_column(b%fs%CFLv_z,'Viscous zCFL')
         call b%cflfile%write()
         ! Create scalar monitor
         b%scfile=monitor(b%nn%cfg%amRoot,'scalar')
         call b%scfile%add_column(b%time%n,'Timestep number')
         call b%scfile%add_column(b%time%t,'Time')
         do nsc=1,b%nn%nscalar
            call b%scfile%add_column(b%nn%SCmin(nsc),trim(b%nn%SCname(nsc))//'_min')
            call b%scfile%add_column(b%nn%SCmax(nsc),trim(b%nn%SCname(nsc))//'_max')
         end do
         call b%scfile%write()
      end block create_monitor

      ! Create a specialized post-processing file
      structure_postproc: block
         ! Create directory to write stats to
         if (b%cfg%amRoot) call execute_command_line('mkdir -p stats')
         ! Plane to set binning domain
         call param_read('Analysis plane',x_over_xL)
         ! ! Perform the output
         ! call structure_identification(b)
      end block structure_postproc
      
   end subroutine init
   
   
   !> Take a time step with block 2
   subroutine step(b,Uinflow,Vinflow,Winflow)
      use tpns_class, only: arithmetic_visc
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
         use tpns_class, only: bcond
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

      ! ! Calculate SR
		! call b%fs%get_strainrate(SR=b%SR)

      ! ! Model shear thinning viscosity
		! update_viscosity: block
      !    integer :: i,j,k
      !    do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
      !      do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
      !        do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
      !          ! Strain rate magnitude
      !          b%SRmag(i,j,k)=sqrt(2.00_WP*b%SR(1,i,j,k)**2+b%SR(2,i,j,k)**2+b%SR(3,i,j,k)**2+2.0_WP*(b%SR(4,i,j,k)**2+b%SR(5,i,j,k)**2+b%SR(6,i,j,k)**2))
      !          ! Carreau Model
      !          b%fs%visc_l(i,j,k)=visc_inf+(visc_0-visc_inf)*(1.00_WP+(alpha*b%SRmag(i,j,k))**2)**((power_constant-1.00_WP)/2.00_WP)
      !          b%visc_norm(i,j,k)=b%fs%visc_l(i,j,k)/visc_0
      !        end do
      !      end do
      !    end do
      !    call b%fs%cfg%sync(b%fs%visc_l)
      ! end block update_viscosity
         
      ! Remember old VOF
      b%vf%VFold=b%vf%VF
         
      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W

      ! Remember old scalars
      b%nn%SCold=b%nn%SC
         
      ! Prepare old staggered density (at n)
      call b%fs%get_olddensity(vf=b%vf)
         
      ! VOF solver step
      call b%vf%advance(dt=b%time%dt,U=b%fs%U,V=b%fs%V,W=b%fs%W)
         
      ! Prepare new staggered viscosity (at n+1)
      call b%fs%get_viscosity(vf=b%vf,strat=arithmetic_visc)

      ! Calculate SR
      call b%fs%get_strainrate(SR=b%SR)
         
      ! Calculate grad(U)
      call b%fs%get_gradU(b%gradU)
         
      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)

         ! ============= SCALAR SOLVER =======================
            
         ! Reset interpolation metrics to QUICK scheme
         call b%nn%metric_reset()
            
         ! Build mid-time scalar
         b%nn%SC=0.5_WP*(b%nn%SC+b%nn%SCold)
         
         ! Explicit calculation of drhoSC/dt from scalar equation
         call b%nn%get_drhoSCdt(b%resSC,b%fs%Uold,b%fs%Vold,b%fs%Wold)
         
         ! Perform bquick procedure
         bquick: block
            integer :: i,j,k
            logical, dimension(:,:,:), allocatable :: flag
            ! Allocate work array
            allocate(flag(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
            ! Assemble explicit residual
            b%resSC=-2.0_WP*(b%nn%SC-b%nn%SCold)+b%time%dt*b%resSC
            ! Apply it to get explicit scalar prediction
            b%SCtmp=2.0_WP*b%nn%SC-b%nn%SCold+b%resSC
            ! Check cells that require bquick
            do k=b%nn%cfg%kmino_,b%nn%cfg%kmaxo_
               do j=b%nn%cfg%jmino_,b%nn%cfg%jmaxo_
                  do i=b%nn%cfg%imino_,b%nn%cfg%imaxo_
                     if (b%SCtmp(i,j,k,1).le.0.0_WP.or.b%SCtmp(i,j,k,4).le.0.0_WP.or.b%SCtmp(i,j,k,6).le.0.0_WP.or.&
                     &   b%SCtmp(i,j,k,1)+b%SCtmp(i,j,k,4)+b%SCtmp(i,j,k,6).ge.b%nn%Lmax**2) then
                        flag(i,j,k)=.true.
                     else
                        flag(i,j,k)=.false.
                     end if
                  end do
               end do
            end do
            ! Adjust metrics
            call b%nn%metric_adjust(b%SCtmp,flag)
            ! Clean up
            deallocate(flag)
            ! Recompute drhoSC/dt
            call b%nn%get_drhoSCdt(b%resSC,b%fs%Uold,b%fs%Vold,b%fs%Wold)
         end block bquick
         
         ! Add fene sources
         call b%nn%addsrc_CgradU(b%gradU,b%resSC)
         call b%nn%addsrc_relax(b%resSC)
         
         ! Assemble explicit residual
         b%resSC=-2.0_WP*(b%nn%SC-b%nn%SCold)+b%time%dt*b%resSC
         
         ! Form implicit residual
         call b%nn%solve_implicit(b%time%dt,b%resSC,b%fs%Uold,b%fs%Vold,b%fs%Wold)
         
         ! Update scalars
         b%nn%SC=2.0_WP*b%nn%SC-b%nn%SCold+b%resSC
         
         ! Force the gas scalar to identity
         gas_scalar_forcing: block
            integer :: i,j,k
            do k=b%nn%cfg%kmino_,b%nn%cfg%kmaxo_
               do j=b%nn%cfg%jmino_,b%nn%cfg%jmaxo_
                  do i=b%nn%cfg%imino_,b%nn%cfg%imaxo_
                     if (b%nn%mask(i,j,k).eq.0) then
                        b%nn%SC(i,j,k,1)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,1)+(1.0_WP-b%vf%VF(i,j,k))*1.0_WP
                        b%nn%SC(i,j,k,2)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,2)
                        b%nn%SC(i,j,k,3)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,3)
                        b%nn%SC(i,j,k,4)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,4)+(1.0_WP-b%vf%VF(i,j,k))*1.0_WP
                        b%nn%SC(i,j,k,5)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,5)
                        b%nn%SC(i,j,k,6)=b%vf%VF(i,j,k)*b%nn%SC(i,j,k,6)+(1.0_WP-b%vf%VF(i,j,k))*1.0_WP
                     end if
                  end do
               end do
            end do
         end block gas_scalar_forcing

         ! Apply all other boundary conditions on the resulting field
         call b%nn%apply_bcond(b%time%t,b%time%dt)
         ! ===================================================
            
         ! ============ VELOCITY SOLVER ======================

         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)
            
         ! Preliminary mass and momentum transport step at the interface
         call b%fs%prepare_advection_upwind(dt=b%time%dt)
            
         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)
            
         ! Add momentum source terms
         call b%fs%addsrc_gravity(b%resU,b%resV,b%resW)

         ! Add polymer stress term
         polymer_stress: block
            integer :: i,j,k,n
            real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
            real(WP), dimension(:,:,:,:), allocatable :: stress
            ! Allocate work arrays
            allocate(stress(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_,1:6))
            allocate(Txy   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
            allocate(Tyz   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
            allocate(Tzx   (b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_))
            ! Calculate the polymer relaxation
            stress=0.0_WP; call b%nn%addsrc_relax(stress)
            ! Build liquid stress tensor - rate depdenent viscosity (Carreau model from Ohta et al. 2019)
            do n=1,6
               b%SRmag=0.0_WP; b%nn%visc=visc_0-visc_l;b%fs%visc_l=visc_0
               do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
                  do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
                     do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
                        ! Cell Strain rate magnitude
                        b%SRmag(i,j,k)=sqrt(b%SR(1,i,j,k)**2+b%SR(2,i,j,k)**2+b%SR(3,i,j,k)**2+2.0_WP*(b%SR(4,i,j,k)**2+b%SR(5,i,j,k)**2+b%SR(6,i,j,k)**2))
                        ! Rate depdentdent polymer viscostiy
                        b%nn%visc=(visc_0-visc_l)*(1.00_WP+(alpha*b%SRmag(i,j,k))**2)**((power_constant-1.00_WP)/2.00_WP)
                        b%fs%visc_l(i,j,k)=visc_l+b%nn%visc
                        b%visc_norm(i,j,k)=b%fs%visc_l(i,j,k)/visc_0
                        ! Viscoelastic stress
                        stress(i,j,k,n)=-b%nn%visc*b%vf%VF(i,j,k)*stress(i,j,k,n)
                     end do
                  end do
               end do
               
            end do
            ! ! Build liquid stress tensor - constant viscosity
            ! do n=1,6
            !    stress(:,:,:,n)=-nn%visc*vf%VF*stress(:,:,:,n)
            ! end do
            ! Interpolate tensor components to cell edges
            do k=b%cfg%kmin_,b%cfg%kmax_+1
               do j=b%cfg%jmin_,b%cfg%jmax_+1
                  do i=b%cfg%imin_,b%cfg%imax_+1
                     Txy(i,j,k)=sum(b%fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
                     Tyz(i,j,k)=sum(b%fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
                     Tzx(i,j,k)=sum(b%fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
                  end do
               end do
            end do
            ! Add divergence of stress to residual
            do k=b%fs%cfg%kmin_,b%fs%cfg%kmax_
               do j=b%fs%cfg%jmin_,b%fs%cfg%jmax_
                  do i=b%fs%cfg%imin_,b%fs%cfg%imax_
                     if (b%fs%umask(i,j,k).eq.0) b%resU(i,j,k)=b%resU(i,j,k)+sum(b%fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
                     &                                                      +sum(b%fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
                     &                                                      +sum(b%fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
                     if (b%fs%vmask(i,j,k).eq.0) b%resV(i,j,k)=b%resV(i,j,k)+sum(b%fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
                     &                                                      +sum(b%fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
                     &                                                      +sum(b%fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
                     if (b%fs%wmask(i,j,k).eq.0) b%resW(i,j,k)=b%resW(i,j,k)+sum(b%fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
                     &                                                      +sum(b%fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &                  
                     &                                                      +sum(b%fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
                  end do
               end do
            end do
            ! Clean up
            deallocate(stress,Txy,Tyz,Tzx)
         end block polymer_stress
            
         ! Assemble explicit residual
         b%resU=-2.0_WP*b%fs%rho_U*b%fs%U+(b%fs%rho_Uold+b%fs%rho_U)*b%fs%Uold+b%time%dt*b%resU
         b%resV=-2.0_WP*b%fs%rho_V*b%fs%V+(b%fs%rho_Vold+b%fs%rho_V)*b%fs%Vold+b%time%dt*b%resV
         b%resW=-2.0_WP*b%fs%rho_W*b%fs%W+(b%fs%rho_Wold+b%fs%rho_W)*b%fs%Wold+b%time%dt*b%resW
            
         ! Form implicit residuals
         call b%fs%solve_implicit(b%time%dt,b%resU,b%resV,b%resW)
            
         ! Apply these residuals
         b%fs%U=2.0_WP*b%fs%U-b%fs%Uold+b%resU
         b%fs%V=2.0_WP*b%fs%V-b%fs%Vold+b%resV
         b%fs%W=2.0_WP*b%fs%W-b%fs%Wold+b%resW
            
         ! Apply other boundary conditions
         call b%fs%apply_bcond(b%time%t,b%time%dt)
            
         ! Solve Poisson equation
         call b%fs%update_laplacian()
         call b%fs%correct_mfr()
         call b%fs%get_div()
         call b%fs%add_surface_tension_jump(dt=b%time%dt,div=b%fs%div,vf=b%vf)
         b%fs%psolv%rhs=-b%fs%cfg%vol*b%fs%div/b%time%dt
         b%fs%psolv%sol=0.0_WP
         call b%fs%psolv%solve()
         call b%fs%shift_p(b%fs%psolv%sol)
            
         ! Correct velocity
         call b%fs%get_pgrad(b%fs%psolv%sol,b%resU,b%resV,b%resW)
         b%fs%P=b%fs%P+b%fs%psolv%sol
         b%fs%U=b%fs%U-b%time%dt*b%resU/b%fs%rho_U
         b%fs%V=b%fs%V-b%time%dt*b%resV/b%fs%rho_V
         b%fs%W=b%fs%W-b%time%dt*b%resW/b%fs%rho_W
         ! ===================================================
            
         ! Increment sub-iteration counter
         b%time%it=b%time%it+1
            
      end do
         
      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         real(WP) :: my_vof_removed
         my_vof_removed=0.0_WP
         do n=1,b%vof_removal_layer%no_
            i=b%vof_removal_layer%map(1,n)
            j=b%vof_removal_layer%map(2,n)
            k=b%vof_removal_layer%map(3,n)
            my_vof_removed=my_vof_removed+b%cfg%vol(i,j,k)*b%vf%VF(i,j,k)
            b%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(my_vof_removed,b%vof_removed,1,MPI_REAL_WP,MPI_SUM,b%cfg%comm,ierr)
      end block remove_vof

      ! Output to ensight
      if (b%ens_evt%occurs()) then 
         ! Update surfmesh object
         update_smesh: block
            use irl_fortran_interface
            integer :: nplane,np,i,j,k
            ! Transfer polygons to smesh
            call b%vf%update_surfmesh(b%smesh)
            ! Also populate nplane variable
            b%smesh%var(1,:)=1.0_WP
            ! Initalize variables to 0
            b%smesh%var(2,:)=0.0_WP
            b%smesh%var(3,:)=0.0_WP
            b%smesh%var(4,:)=0.0_WP
            b%smesh%var(5,:)=0.0_WP
            b%smesh%var(6,:)=0.0_WP
            b%smesh%var(7,:)=0.0_WP
            np=0
            do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
               do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                  do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1
                           b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                           b%smesh%var(2,np)=b%Ui(i,j,k)
                           b%smesh%var(3,np)=b%Vi(i,j,k)
                           b%smesh%var(4,np)=b%Wi(i,j,k)
                           b%smesh%var(5,np)=b%fs%visc_l(i,j,k)
                           b%smesh%var(6,np)=b%SRmag(i,j,k)
                           b%smesh%var(7,np)=b%visc_norm(i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Perform ensight output 
         call b%ens_out%write_data(b%time%t)
      end if
         
      ! Perform and output monitoring
      call b%nn%get_max()
      call b%fs%get_max()
      call b%vf%get_max()
      call b%mfile%write()
      call b%cflfile%write()
      call b%scfile%write()

      ! Specialized post-processing
      call structure_identification(b)

      
   end subroutine step

   !> Collect information about liquid structures
   subroutine structure_identification(b)
      implicit none
      class(block2), intent(inout) :: b
      
      ! Perform detailed CCL pass to identify structures
      b%cc%max_interface_planes=2
      call b%cc%build_lists(VF=b%vf%VF,poly=b%vf%interface_polygon,U=b%fs%U,V=b%fs%V,W=b%fs%W)
      call b%cc%get_min_thickness()
      call b%cc%sort_by_thickness()
      
      ! Loop through identified detached structs and sort by ligaments and particles
      print_structures: block
         integer :: m,n,i,j,k

         ! Loops over film segments contained locally
         do m=1,b%cc%n_meta_struct

            ! Check if structure x_cg .ge. x plane location
            if (b%cc%meta_structures_list(m)%x.ge.x_over_xL*b%cfg%xL) then

               ! Call writting subroutine to output structure stats
               call stats_writer(b,b%time%t,b%cc%meta_structures_list(m)%vol,b%cc%meta_structures_list(m)%x,b%cc%meta_structures_list(m)%y,b%cc%meta_structures_list(m)%z,b%cc%meta_structures_list(m)%u,b%cc%meta_structures_list(m)%v,b%cc%meta_structures_list(m)%w,b%cc%meta_structures_list(m)%lengths(1),b%cc%meta_structures_list(m)%lengths(2),b%cc%meta_structures_list(m)%lengths(3)) 
            
            end if

         end do
         
      end block print_structures
      
      ! Clean up CCL
      call b%cc%deallocate_lists()

      
   end subroutine structure_identification
   
   !> Specialized subroutine to print our stats of liquid structures
   subroutine stats_writer(b,time,volume,x_cg,y_cg,z_cg,U_vel,V_vel,W_vel,length_1,length_2,length_3)
      use string,      only: str_medium
      implicit none
      class(block2), intent(inout) :: b
      real(WP) :: time,volume,x_cg,y_cg,z_cg,U_vel,V_vel,W_vel,length_1,length_2,length_3
      integer  :: iunit,ierr
      logical  :: exist
      ! Print stats from root processor
      if (b%cfg%amRoot) then
         inquire(file='./stats/structure_stats',exist=exist)
         if (exist) then
            open(newunit=iunit,file='./stats/structure_stats',form='formatted',status='old',position='append',access='stream',iostat=ierr)
         else
            open(newunit=iunit,file='./stats/structure_stats',form='formatted',status='new',access='stream',iostat=ierr)
            write(iunit,'(a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a12)') 'time','volume','x_cg','y_cg','z_cg','U','V','W','length_1','length_2','length_3'
         end if
         ! Write Data
         write(iunit,'(es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5)') time,volume,x_cg,y_cg,z_cg,U_vel,V_vel,W_vel,length_1,length_2,length_3
         close(iunit)
      end if
   end subroutine stats_writer
   
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
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR,b%SRmag)
      deallocate(b%resSC,b%SCtmp,b%gradU,b%SR,b%SRmag)
      
   end subroutine final
   
   
end module block2_class
