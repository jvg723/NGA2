!> Definition for block 2 simulation: swirl atomizer
module block2_class
   use precision,          only: WP
   use geometry,           only: Dout,Din
   use config_class,       only: config
   use iterator_class,     only: iterator
   use hypre_str_class,    only: hypre_str
   use ddadi_class,        only: ddadi
   use tpns_class,         only: tpns
   use vfs_class,          only: vfs
   use lpt_class,          only: lpt
   use stokeslet_class,    only: stokeslet
   use timetracker_class,  only: timetracker
   use ensight_class,      only: ensight
   use surfmesh_class,     only: surfmesh
   use partmesh_class,     only: partmesh
   use event_class,        only: event
   use cclabel_class,      only: cclabel
   use monitor_class,      only: monitor
   use breakup_class,      only: breakup
   implicit none
   private

   public :: block2
   
   !> block 2 object
   type :: block2
      class(config), pointer :: cfg                    !< Pointer to config
      type(hypre_str)        :: ps                     !< Structured hypre pressure solver
      type(ddadi)            :: vs                     !< DDAI implicit solver
      type(tpns)             :: fs                     !< Two phase incompressible flow solver
      type(vfs)              :: vf                     !< VF solve
      type(timetracker)      :: time                   !< Time tracker
      type(surfmesh)         :: smesh                  !< Surfmesh   
      type(partmesh)         :: pmesh,pmesh_lpt        !< Particle mesh for stokeslets           
      type(ensight)          :: ens_out                !< Ensight output
      type(event)            :: ens_evt                !< Ensight event
      type(cclabel)          :: ccl                    !< Label object for thin regions
      !> Break-up modeling
      type(lpt)         :: lp    !< Lagrangian particle solver
      type(stokeslet)   :: ss
      type(breakup)     :: bu    !< SGS break-up model
      !> Monitoring files
      type(monitor) :: mfile     !< General monitor files
      type(monitor) :: cflfile   !< CFL monitor files
      !> Private work arrays
      real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
      real(WP), dimension(:,:,:),     allocatable :: Uslip,Vslip,Wslip
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2

   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=4

   real(WP) :: VFlolim


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
         allocate(b%Uslip(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_)); b%Uslip=0.0_WP
         allocate(b%Vslip(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_)); b%Vslip=0.0_WP
         allocate(b%Wslip(b%cfg%imino_:b%cfg%imaxo_,b%cfg%jmino_:b%cfg%jmaxo_,b%cfg%kmino_:b%cfg%kmaxo_)); b%Wslip=0.0_WP
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
         use vfs_class, only: VFlo,VFhi,elvira,r2p,lvira,flux,remap,r2plig
         integer :: i,j,k,n,si,sj,sk
			real(WP), dimension(3,8) :: cube_vertex
			real(WP), dimension(3) :: v_cent,a_cent
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call b%vf%initialize(cfg=b%cfg,reconstruction_method=r2plig,transport_method=flux,name='VOF')
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
         ! Min Film thickess
         call param_read('Minimum local thickness',b%vf%thin_thld_min)
         ! call b%input%read('nx',nx)
         ! min_thickness = b%vf%thin_thld_min*b%cfg%xL/(1.0_WP*nx)
         ! b%vf%puncture = .false.
         VFlolim = VFlo 
         ! b%vf%vol_conv = 0.0_WP
         ! b%vf%vol_loss = 0.0_WP
         ! b%vf%vol_left = 0.0_WP
         ! b%vf%punctured = 0.0_WP
      end block create_and_initialize_vof


      ! Create an iterator for removing VOF at edges
      create_iterator: block
         b%vof_removal_layer=iterator(b%cfg,'VOF removal',vof_removal_layer_locator)
      end block create_iterator
      

      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,slip
         use hypre_str_class, only: pcg_pfmg2
         integer :: i,j,k
         ! Create flow solver
         b%fs=tpns(cfg=b%cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',b%fs%visc_l)
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
         b%ps=hypre_str(cfg=b%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         b%ps%maxlevel=16
         call param_read('Pressure iteration',b%ps%maxit)
         call param_read('Pressure tolerance',b%ps%rcvg)
         ! Configure velocity solver
			b%vs=ddadi(cfg=b%cfg,name='Velocity',nst=7)
			! Setup the solver
			call b%fs%setup(pressure_solver=b%ps,implicit_solver=b%vs)
      end block create_and_initialize_flow_solver
      

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

      ! Create lpt solver
      create_lpt_solver: block
         ! Create the solver
         b%lp=lpt(cfg=b%cfg,name='spray')
         ! Get particle density from the flow solver
         b%lp%rho=b%fs%rho_l
         ! Turn off drag
         b%lp%drag_model='none'         
         ! Initialize with zero particles
         call b%lp%resize(0)
         ! Get initial particle volume fraction
         call b%lp%update_VF()               
         ! Get particle statistics
         call b%lp%get_max()
      end block create_lpt_solver

      create_breakup_model: block
         call b%bu%initialize(b%vf,b%fs,b%lp)
         b%vf%thin_thld_min=b%bu%min_filmthickness/b%vf%cfg%min_meshsize
      end block create_breakup_model

      ! Create surfmesh object for interface polygon output 
      create_smesh: block
         use irl_fortran_interface
         use vfs_class, only: lvira,r2p
         integer :: i,j,k,nplane,np
         select case (b%vf%reconstruction_method)
         case (r2p)
            ! Include an extra variable for number of planes
            b%smesh=surfmesh(nvar=13,name='plic')
            b%smesh%varname(1)='nplane'
            b%smesh%varname(2)='curv'
            b%smesh%varname(3)='edge_sensor'
            b%smesh%varname(4)='thin_sensor'
            b%smesh%varname(5)='thickness'
            b%smesh%varname(6)='id_ccl_film'
            b%smesh%varname(7)='edge_sensor'
            b%smesh%varname(8)='Uslip'
            b%smesh%varname(9)='Vslip'
            b%smesh%varname(10)='Wslip'
            b%smesh%varname(11)='film_type'
            b%smesh%varname(12)='id_ccl_edge'
            b%smesh%varname(13)='lig_ind'
            ! Transfer polygons to smesh
            call b%vf%update_surfmesh(b%smesh)
            ! Also populate nplane variable
            b%smesh%var(1,:)=1.0_WP
            np=0
            do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
               do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                  do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                           b%smesh%var(2,np)=b%vf%curv2p(nplane,i,j,k)
                           b%smesh%var(3,np)=b%vf%edge_sensor(i,j,k)
                           b%smesh%var(4,np)=b%vf%thin_sensor(i,j,k)
                           b%smesh%var(5,np)=b%vf%thickness  (i,j,k)
                           b%smesh%var(6,np)=real(b%bu%ccl_film%id(i,j,k),WP)
                           b%smesh%var(7,np)=b%vf%edge_sensor(i,j,k)
                           b%smesh%var(8,np)=b%Uslip(i,j,k)
                           b%smesh%var(9,np)=b%Vslip(i,j,k)
                           b%smesh%var(10,np)=b%Wslip(i,j,k)
                           b%smesh%var(11,np)=real(b%vf%film_type(i,j,k),WP)
                           b%smesh%var(12,np)=real(b%bu%ccl_edge%id(i,j,k),WP)
                           b%smesh%var(13,np)=real(b%vf%lig_ind(i,j,k),WP)
                        end if
                     end do
                  end do
               end do
            end do
         case (lvira)
            ! Include an extra variable for number of planes
            b%smesh=surfmesh(nvar=1,name='plic')
            b%smesh%varname(1)='nplane'
            ! Transfer polygons to smesh
            call b%vf%update_surfmesh(b%smesh)
            ! Also populate nplane variable
            b%smesh%var(1,:)=1.0_WP
            np=0
            do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
               do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                  do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                        end if
                     end do
                  end do
               end do
            end do
         end select
      end block create_smesh

      create_pmesh: block
         integer :: i
         call b%ss%initialize(cfg=b%cfg,name='stokeslet')
         b%pmesh=partmesh(nvar=2,nvec=1,name='stokeslet')
         ! b%pmesh%varname(1)='id'
         b%pmesh%varname(1)='fcoeff'
         b%pmesh%varname(2)='VF'
         b%pmesh%vecname(1)='nedge'
         call b%ss%update_partmesh(b%pmesh)
         do i=1,b%ss%np_
            ! b%pmesh%var(1,i)=ss%p(i)%id
            b%pmesh%var(1,i)=b%ss%p(i)%fcoeff
            b%pmesh%var(1,i)=b%ss%p(i)%VF
            b%pmesh%vec(:,1,i)=b%ss%p(i)%nedge
         end do
      end block create_pmesh

      ! Create partmesh object for Lagrangian particle output
      create_pmesh_lpt: block
         integer :: i
         ! Include an extra variable for droplet diameter
         b%pmesh_lpt=partmesh(nvar=2,nvec=1,name='lpt')
         b%pmesh_lpt%varname(1)='diameter'
         b%pmesh_lpt%varname(2)='id'
         b%pmesh_lpt%vecname(1)='velocity'
         ! Transfer particles to pmesh
         call b%lp%update_partmesh(b%pmesh_lpt)
         ! Also populate diameter variable
         do i=1,b%lp%np_
            b%pmesh_lpt%var(1,i)=b%lp%p(i)%d
            b%pmesh_lpt%var(2,i)=b%lp%p(i)%id
            b%pmesh_lpt%vec(:,1,i)=b%lp%p(i)%vel
         end do
      end block create_pmesh_lpt

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
         call b%ens_out%add_vector('SLIPvelocity',b%Uslip,b%Vslip,b%Wslip)
         call b%ens_out%add_scalar('VOF',b%vf%VF)
         call b%ens_out%add_scalar('Pressure',b%fs%P)
         call b%ens_out%add_scalar('thin_sensor',b%vf%thin_sensor)
         call b%ens_out%add_scalar('curvature',b%vf%curv)
         call b%ens_out%add_surface('vofplic',b%smesh)
         call b%ens_out%add_particle('stokeslets',b%pmesh)
         call b%ens_out%add_particle('spray',b%pmesh_lpt)
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
      end block create_monitor

      
   end subroutine init
   
   
   !> Take a time step with block 2
   subroutine step(b,Uinflow,Vinflow,Winflow)
      use tpns_class, only: arithmetic_visc
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_WTIME
      use parallel, only: MPI_REAL_WP
      implicit none
      class(block2), intent(inout) :: b
      real(WP) :: starttime_step,endtime_step,my_time_step
      integer :: ierr
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

      ! Advance our spray
      b%resU=b%fs%rho_g; b%resV=b%fs%visc_g
      call b%lp%advance(dt=b%time%dt,U=b%fs%U,V=b%fs%V,W=b%fs%W,rho=b%resU,visc=b%resV)   
         
      ! Remember old VOF
      b%vf%VFold=b%vf%VF
         
      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W
         
      ! Prepare old staggered density (at n)
      call b%fs%get_olddensity(vf=b%vf)

         
      ! VOF solver step
      if (b%vf%cfg%amRoot) print *,"ligamentclasshere1"
      call b%fs%add_slipvel(b%vf,b%ss,b%Uslip,b%Vslip,b%Wslip,b%bu%min_filmthickness)
      ! if (b%vf%cfg%amRoot) print *,"herebefore"
      if (b%vf%cfg%amRoot) print *,"ligamentclasshere2"
      b%Uslip =b%fs%U + b%Uslip
      b%Vslip =b%fs%V + b%Vslip
      b%Wslip =b%fs%W + b%Wslip
      if (b%ss%np.gt.0) then
         if (b%vf%cfg%amRoot) print *,"ligamentclasshere3"
         call b%vf%advance(dt=b%time%dt,U=b%Uslip,V=b%Vslip,W=b%Wslip,type_flux=2)
         if (b%vf%cfg%amRoot) print *,"ligamentclasshere4"
      else 
         if (b%vf%cfg%amRoot) print *,"ligamentclasshere5"
         call b%vf%advance(dt=b%time%dt,U=b%Uslip,V=b%Vslip,W=b%Wslip,type_flux=1)
         if (b%vf%cfg%amRoot) print *,"ligamentclasshere6"
      end if
      ! if (b%vf%cfg%amRoot) print *,"here12"
      ! ! if (b%vf%cfg%amRoot) print *,"after add slip"
      
         
      ! Prepare new staggered viscosity (at n+1)
      call b%fs%get_viscosity(vf=b%vf,strat=arithmetic_visc)
         
      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)
            
         ! ============ VELOCITY SOLVER ======================

         ! Build mid-time velocity
         b%fs%U=0.5_WP*(b%fs%U+b%fs%Uold)
         b%fs%V=0.5_WP*(b%fs%V+b%fs%Vold)
         b%fs%W=0.5_WP*(b%fs%W+b%fs%Wold)
            
         ! Preliminary mass and momentum transport step at the interface
         call b%fs%prepare_advection_upwind(dt=b%time%dt)
            
         ! Explicit calculation of drho*u/dt from NS
         call b%fs%get_dmomdt(b%resU,b%resV,b%resW)
            
         ! ! Add momentum source terms
         ! call b%fs%addsrc_gravity(b%resU,b%resV,b%resW)
            
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
         solve_poisson: block
            use vfs_class, only: lvira,r2p
            call b%fs%update_laplacian()
            call b%fs%correct_mfr()
            call b%fs%get_div()
            select case (b%vf%reconstruction_method)
            case (r2p)
               ! call b%fs%add_surface_tension_jump_thin(dt=b%time%dt,div=b%fs%div,vf=b%vf)
               call b%fs%add_surface_tension_jump_twoVF(dt=b%time%dt,div=b%fs%div,vf=b%vf)
            case (lvira)
               call b%fs%add_surface_tension_jump(dt=b%time%dt,div=b%fs%div,vf=b%vf)
            end select
            b%fs%psolv%rhs=-b%fs%cfg%vol*b%fs%div/b%time%dt
            b%fs%psolv%sol=0.0_WP
            call b%fs%psolv%solve()
            call b%fs%shift_p(b%fs%psolv%sol)
         end block solve_poisson
            
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

      ! label film and ligament based on local moment of inertia information
      label_film_ligament: block
         if (b%vf%cfg%amRoot) print *,"ligamentclasshere7"
         ! call b%get_neighbortype()
         ! if (b%vf%cfg%amRoot) print *,"here4"
         call b%bu%attempt_breakup_retraction()
         ! if (b%vf%cfg%amRoot) print *,"here5"
      end block label_film_ligament

      
      ! Remove VOF at edge of domain
      remove_vof: block
         integer :: n
         do n=1,b%vof_removal_layer%no_
            b%vf%VF(b%vof_removal_layer%map(1,n),b%vof_removal_layer%map(2,n),b%vof_removal_layer%map(3,n))=0.0_WP
         end do
      end block remove_vof

      
      ! Output to ensight
      if (b%ens_evt%occurs()) then 
         ! Update surfmesh object 
         update_smesh: block
            use irl_fortran_interface
            use vfs_class, only: lvira,r2p
            integer :: nplane,np,i,j,k
            select case (b%vf%reconstruction_method)
            case (r2p)
               ! Transfer polygons to smesh
               call b%vf%update_surfmesh(b%smesh)
               ! Also populate nplane variable
               b%smesh%var(1,:)=1.0_WP
               np=0
               do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
                  do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                     do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                              b%smesh%var(2,np)=b%vf%curv2p(nplane,i,j,k)
                              b%smesh%var(3,np)=b%vf%edge_sensor(i,j,k)
                              b%smesh%var(4,np)=b%vf%thin_sensor(i,j,k)
                              b%smesh%var(5,np)=b%vf%thickness  (i,j,k)
                              b%smesh%var(6,np)=real(b%bu%ccl_film%id(i,j,k),WP)
                              b%smesh%var(7,np)=b%vf%edge_sensor(i,j,k)
                              b%smesh%var(8,np)=b%Uslip(i,j,k)
                              b%smesh%var(9,np)=b%Vslip(i,j,k)
                              b%smesh%var(10,np)=b%Wslip(i,j,k)
                              b%smesh%var(11,np)=real(b%vf%film_type(i,j,k),WP)
                              b%smesh%var(12,np)=real(b%bu%ccl_edge%id(i,j,k),WP)
                              b%smesh%var(13,np)=real(b%vf%lig_ind(i,j,k),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            case (lvira)
               ! Transfer polygons to smesh
               call b%vf%update_surfmesh(b%smesh)
               ! Also populate nplane variable
               b%smesh%var(1,:)=1.0_WP
               np=0
               do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
                  do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                     do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
                           end if
                        end do
                     end do
                  end do
               end do
            end select
         end block update_smesh

         update_pmesh: block
            integer :: i
            call b%ss%update_partmesh(b%pmesh)
            do i=1,b%ss%np_
               b%pmesh%var(1,i)=b%ss%p(i)%fcoeff
               b%pmesh%var(2,i)=b%ss%p(i)%VF
               b%pmesh%vec(:,1,i)=b%ss%p(i)%nedge
            end do
         end block update_pmesh  

         update_pmesh_lpt: block
            integer :: i
            call b%lp%update_partmesh(b%pmesh_lpt)
            do i=1,b%lp%np_
               b%pmesh_lpt%var(1,i)=b%lp%p(i)%d
               b%pmesh_lpt%var(2,i)=b%lp%p(i)%id
               b%pmesh_lpt%vec(:,1,i)=b%lp%p(i)%vel
            end do
         end block update_pmesh_lpt  
   
         ! Perform ensight output 
         call b%ens_out%write_data(b%time%t)
      end if
         
      ! Perform and output monitoring
      call b%fs%get_max()
      call b%vf%get_max()
      call b%mfile%write()
      call b%cflfile%write()

   end subroutine step
   
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
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%Uslip,b%Vslip,b%Wslip)
      
   end subroutine final

   
end module block2_class
