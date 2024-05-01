!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use iterator_class,    only: iterator
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use cclabel_class,     only: cclabel
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps                     !< Structured hypre pressure solver
   type(ddadi),       public :: vs                     !< DDAI implicit solver
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(cclabel),     public :: ccl                    !< Label object
   type(timetracker), public :: time
   type(surfmesh),    public :: smesh                                              
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: Uslip,Vslip,Wslip

   !> Iterator for VOF removal
   type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
   real(WP) :: vof_removed              !< Integral of VOF removed
   
   !> Problem definition
   real(WP) :: R1,R2   !< Inner and outer radii of annulus
   real(WP) :: Ul,SW   !< Liquid axial velocity and swirl ratio

   !> Hardcode size of buffer layer for VOF removal
   integer, parameter :: nlayer=4

   !> Min film thickness 
   real(WP), parameter :: min_filmthickness=4e-2_WP
   
contains
   
   
   !> Function that defines a level set function for an annular liquid region
	function levelset_annulus(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G,r
      r=sqrt(xyz(2)**2+xyz(3)**2)
	   G=min(r-R1,R2-r)
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

   !> Function that identifies cells that need a label to min thickness region 
   logical function make_label(i,j,k)
      implicit none
      integer, intent(in) :: i,j,k
      if (vf%thickness(i,j,k).le.min_filmthickness.and.vf%thickness(i,j,k).gt.0.0_WP) then
         make_label=.true.
      else
         make_label=.false.
      end if
   end function make_label

   !> Function that identifies if cell pairs have same label
   logical function same_label(i1,j1,k1,i2,j2,k2)
      implicit none
      integer, intent(in) :: i1,j1,k1,i2,j2,k2
      same_label=.true.
   end function same_label
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Uslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Uslip=0.0_WP
         allocate(Vslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vslip=0.0_WP
         allocate(Wslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wslip=0.0_WP
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='pressure_swirl')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
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
         call vf%initialize(cfg=cfg,reconstruction_method=r2p,name='VOF')
         ! Set full domain to gas
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  vf%VF(i,j,k)=0.0_WP
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Initialize an annular interface at the inlet only
         call param_read('Inner radius',R1)
         call param_read('Outer radius',R2)
         if (vf%cfg%iproc.eq.1) then
            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino,vf%cfg%imin-1
                     ! Set cube vertices
				         n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                           end do
                        end do
                     end do
                     ! Call adaptive refinement code to get volume and barycenters recursively
				         vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_annulus,0.0_WP,amr_ref_lvl)
                     vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                     if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                        vf%Lbary(:,i,j,k)=v_cent
                        vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                     else
                        vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                        vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     end if
                  end do
               end do
            end do
         end if
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof

      ! Create an iterator for removing VOF at edges
      create_iterator: block
         vof_removal_layer=iterator(cfg,'VOF removal',vof_removal_layer_locator)
         vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,slip
         use hypre_str_class, only: pcg_pfmg2
         integer :: i,j,k
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity'   ,fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density'   ,fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Inflow on the left of domain
         call fs%add_bcond(name='inflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Clipped Neumann outflow on the right of domain
         call fs%add_bcond(name='bc_xp' ,type=clipped_neumann,face='x',dir=+1,canCorrect=.false. ,locator=xp_locator)
         ! Slip on the sides
         call fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure velocity solver
			vs=ddadi(cfg=cfg,name='Velocity',nst=7)
			! Setup the solver
			call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_and_initialize_flow_solver
      

      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         real(WP) :: r,theta
         integer  :: n,i,j,k
         ! Zero initial field in the domain
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Read in inflow parameters
         call param_read('Liquid velocity',Ul)
         call param_read('Swirl ratio',SW)
         ! Apply axial and swirl component Dirichlet at inlet
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            ! Set U velocity
            r=sqrt(cfg%ym(j)**2+cfg%zm(k)**2) !< Radius in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) fs%U(i,j,k)=Ul
            ! Set V velocity
            r=sqrt(cfg%y(j)**2+cfg%zm(k)**2)  !< Radius in cylindrical coordinates
            theta=atan2(cfg%y(j),cfg%zm(k))   !< Angle  in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) fs%V(i,j,k)=+4.0_WP*SW*Ul*(-r**2+(R2+R1)*r-R2*R1)/((R2-R1)**2)*cos(theta)
            ! Set W velocity
            r=sqrt(cfg%ym(j)**2+cfg%z(k)**2)  !< Radius in cylindrical coordinates
            theta=atan2(cfg%ym(j),cfg%z(k))   !< Angle  in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) fs%W(i,j,k)=-4.0_WP*SW*Ul*(-r**2+(R2+R1)*r-R2*R1)/((R2-R1)**2)*sin(theta)
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity

      ! Create cclabel object
      film_label: block
         call ccl%initialize(pg=cfg%pgrid,name='thin_region_label')
      end block film_label
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=9,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='edge_sensor'
         smesh%varname(4)='thin_sensor'
         smesh%varname(5)='thickness'
         smesh%varname(6)='id_thickenss'
         smesh%varname(7)='Uslip'
         smesh%varname(8)='Vslip'
         smesh%varname(9)='Wslip'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                        smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                        smesh%var(3,np)=vf%edge_sensor(i,j,k)
                        smesh%var(4,np)=vf%thin_sensor(i,j,k)
                        smesh%var(5,np)=vf%thickness  (i,j,k)
                        smesh%var(6,np)=real(ccl%id(i,j,k),WP)
                        smesh%var(7,np)=Uslip(i,j,k)
                        smesh%var(8,np)=Vslip(i,j,k)
                        smesh%var(9,np)=Wslip(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='swirl_atomizer')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('Pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh)
         call ens_out%add_vector('edge_normal',resU,resV,resW)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(vf%SDint,'SD integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: arithmetic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! Add in slip velocity at hole edges
         slip_velocity: block
            integer :: i,j,k
            ! Store current velocity field 
            Uslip=fs%U
            Vslip=fs%V
            Wslip=fs%W
            ! Add in retraction velocities
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     if (vf%edge_sensor(i,j,k).ge.0.10_WP) then
                     ! if (vf%edge_sensor(i,j,k).ge.0.10_WP) then Local threshold
                        Uslip(i  ,j,k)=Uslip(i  ,j,k)+0.5_WP*vf%edge_normal(1,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                        Uslip(i+1,j,k)=Uslip(i+1,j,k)+0.5_WP*vf%edge_normal(1,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                        Vslip(i,j  ,k)=Vslip(i,j  ,k)+0.5_WP*vf%edge_normal(2,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                        Vslip(i,j+1,k)=Vslip(i,j+1,k)+0.5_WP*vf%edge_normal(2,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                        Wslip(i,j,k  )=Wslip(i,j,k  )+0.5_WP*vf%edge_normal(3,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                        Wslip(i,j,k+1)=Wslip(i,j,k+1)+0.5_WP*vf%edge_normal(3,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*vf%thickness(i,j,k)))
                     end if
                  end do 
               end do 
            end do
         end block slip_velocity

         
         ! VOF solver step
         ! call vf%advance(dt=time%dt,U=Uslip,V=Vslip,W=Wslip)
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! ! Add momentum source terms
            ! call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Label thin film regions
         label_thin: block
            ! Label regions and track time
            call ccl%build(make_label,same_label)
         end block label_thin

         ! Puncture a hole in the film based upon the thin region label
         puncture_film: block
            integer :: i,j,k,nn,n
            do n=1,ccl%nstruct
               do nn=1,ccl%struct(n)%n_
                  i=ccl%struct(n)%map(1,nn); j=ccl%struct(n)%map(2,nn); k=ccl%struct(n)%map(3,nn)
                  vf%VF(i,j,k)=0.0_WP
               end do
            end do
         end block puncture_film

         ! Remove VOF at edge of domain
         remove_vof: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
            use parallel, only: MPI_REAL_WP
            integer :: n,i,j,k,ierr
            real(WP) :: my_vof_removed
            my_vof_removed=0.0_WP
            do n=1,vof_removal_layer%no_
               i=vof_removal_layer%map(1,n)
               j=vof_removal_layer%map(2,n)
               k=vof_removal_layer%map(3,n)
               my_vof_removed=my_vof_removed+cfg%vol(i,j,k)*vf%VF(i,j,k)
               vf%VF(i,j,k)=0.0_WP
            end do
            call MPI_ALLREDUCE(my_vof_removed,vof_removed,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
         end block remove_vof

         ! Output to ensight
         if (ens_evt%occurs()) then 
            ! Update surfmesh object
            update_smesh: block
               use irl_fortran_interface
               integer :: nplane,np,i,j,k
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Also populate nplane variable
               smesh%var(1,:)=1.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                              smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                              smesh%var(3,np)=vf%edge_sensor(i,j,k)
                              smesh%var(4,np)=vf%thin_sensor(i,j,k)
                              smesh%var(5,np)=vf%thickness  (i,j,k)
                              smesh%var(6,np)=real(ccl%id(i,j,k),WP)
                              smesh%var(7,np)=Uslip(i,j,k)
                              smesh%var(8,np)=Vslip(i,j,k)
                              smesh%var(9,np)=Wslip(i,j,k)
                           end if
                        end do
                     end do
                  end do
               end do
            end block update_smesh
            ! Transfer edge normal data
            resU=vf%edge_normal(1,:,:,:)
            resV=vf%edge_normal(2,:,:,:)
            resW=vf%edge_normal(3,:,:,:)
            ! Perform ensight output 
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         
      end do
      
      
   end subroutine simulation_run
      
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      deallocate(Uslip,Vslip,Wslip)
      
   end subroutine simulation_final
   
   
end module simulation
