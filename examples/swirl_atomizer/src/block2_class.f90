!> Definition for block 2 simulation: swirl atomizer
module block2_class
   use precision,         only: WP
   use geometry,          only: D1,D2
   use config_class,      only: config
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
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
   class(config),    pointer :: cfg             !< Pointer to config
   type(tpns),        public :: fs              !< Two phase incompressible flow solver
   type(vfs),         public :: vf              !< VF solver
   type(timetracker), public :: time            !< Time tracker
   type(surfmesh),    public :: smesh           !< Surfmesh                                       
   type(ensight)             :: ens_out         !< Ensight output
   type(event)               :: ens_evt         !< Ensight event
   type(monitor)             :: mfile,cflfile   !< Monitor files
   !> Private work arrays
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: SRmag
   real(WP), dimension(:,:,:,:), allocatable :: SR
   contains
      procedure :: init                   !< Initialize block
      procedure :: step                   !< Advance block
      procedure :: final                  !< Finalize block
   end type block2
   
   !> Problem definition
   real(WP) :: R1,R2   !< Inner and outer radii of annulus
   real(WP) :: Ul,SW   !< Liquid axial velocity and swirl ratio
   
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
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         b%time=timetracker(amRoot=b%cfg%amRoot,name='swirl_atomizer')
         call param_read('1 Max timestep size',b%time%dtmax)
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
         call param_read('Inner radius',R1)
         call param_read('Outer radius',R2)
         if (b%vf%cfg%iproc.eq.1) then
            do k=b%vf%cfg%kmino_b%,vf%cfg%kmaxo_
               do j=b%vf%cfg%jmino_b%,vf%cfg%jmaxo_
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
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: dirichlet,clipped_neumann
         use ils_class,  only: pcg_pfmg
         integer :: i,j,k
         real(WP) :: visc_l,visc_g
         ! Create flow solver
         b%fs=tpns(cfg=b%cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',visc_l); b%fs%visc_l=visc_l
         call param_read('Gas dynamic viscosity'   ,visc_g); b%fs%visc_g=visc_g
         ! Assign constant density to each phase
         call param_read('Liquid density',b%fs%rho_l)
         call param_read('Gas density'   ,b%fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',b%fs%sigma)
         ! Inflow on the left of domain
         call b%fs%add_bcond(name='inflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Clipped Neumann outflow on the right of domain
         call b%fs%add_bcond(name='bc_xp' ,type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         ! Configure pressure solver
         call param_read('Pressure iteration',b%fs%psolv%maxit)
         call param_read('Pressure tolerance',b%fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',b%fs%implicit%maxit)
         call param_read('Implicit tolerance',b%fs%implicit%rcvg)
         ! Setup the solver
         call b%fs%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
      end block create_and_initialize_flow_solver
      

      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         real(WP) :: r,theta
         integer  :: n,i,j,k
         ! Zero initial field in the domain
         b%fs%U=0.0_WP; b%fs%V=0.0_WP; b%fs%W=0.0_WP
         ! Read in inflow parameters
         call param_read('Liquid velocity',Ul)
         call param_read('Swirl ratio',SW)
         ! Apply axial and swirl component Dirichlet at inlet
         call b%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            ! Set U velocity
            r=sqrt(b%cfg%ym(j)**2+b%cfg%zm(k)**2) !< Radius in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) b%fs%U(i,j,k)=Ul
            ! Set V velocity
            r=sqrt(b%cfg%y(j)**2+b%cfg%zm(k)**2)  !< Radius in cylindrical coordinates
            theta=atan2(b%cfg%y(j),b%cfg%zm(k))   !< Angle  in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) b%fs%V(i,j,k)=+4.0_WP*SW*Ul*(-r**2+(R2+R1)*r-R2*R1)/((R2-R1)**2)*cos(theta)
            ! Set W velocity
            r=sqrt(b%cfg%ym(j)**2+b%cfg%z(k)**2)  !< Radius in cylindrical coordinates
            theta=atan2(b%cfg%ym(j),b%cfg%z(k))   !< Angle  in cylindrical coordinates
            if (r.ge.R1.and.r.le.R2) b%fs%W(i,j,k)=-4.0_WP*SW*Ul*(-r**2+(R2+R1)*r-R2*R1)/((R2-R1)**2)*sin(theta)
         end do
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
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         b%smesh=surfmesh(nvar=1,name='plic')
         b%smesh%varname(1)='nplane'
         ! Transfer polygons to smesh
         call b%vf%update_surfmesh(smesh)
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
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
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
         ! Output to ensight
         if (b%ens_evt%occurs()) call b%ens_out%write_data(b%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
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
   subroutine step(b)
      implicit none
      class(block2), intent(inout) :: b
         
      ! Increment time
      call b%fs%get_cfl(b%time%dt,b%time%cfl)
      call b%time%adjust_dt()
      call b%time%increment()

      ! Calculate SR
      ! call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)

      ! ! Model shear thinning fluid
      ! nonewt: block
      !    integer :: i,j,k
      !    real(WP), parameter :: C=1.137e-3_WP
      !    real(WP), parameter :: n=0.3_WP
      !    ! Update viscosity
      !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
      !       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
      !          do i=fs%cfg%imino_,fs%cfg%imaxo_
      !             ! Power law model 
      !             SRmag(i,j,k)=sqrt(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
      !             SRmag(i,j,k)=max(SRmag(i,j,k),1000.0_WP**(1.0_WP/(n-1.0_WP)))
      !             fs%visc_l(i,j,k)=C*SRmag(i,j,k)**(n-1.0_WP)
      !             ! Carreau Model
      !             ! SRmag(i,j,k)=sqrt(2.00_WP*SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
      !             ! fs%visc_l(i,j,k)=visc_inf+(visc_0-visc_inf)*(1.00_WP+(lambda*SRmag(i,j,k))**2.00_WP)**((n-1.00_WP)/2.00_WP)
      !          end do
      !       end do
      !    end do
      !    ! call fs%cfg%sync(fs%visc_l)
      ! end block nonewt
         
      ! Remember old VOF
      b%vf%VFold=b%vf%VF
         
      ! Remember old velocity
      b%fs%Uold=b%fs%U
      b%fs%Vold=b%fs%V
      b%fs%Wold=b%fs%W
         
      ! Prepare old staggered density (at n)
      call b%fs%get_olddensity(vf=b%vf)
         
      ! VOF solver step
      call b%vf%advance(dt=b%time%dt,U=b%fs%U,V=b%fs%V,W=b%fs%W)
         
      ! Prepare new staggered viscosity (at n+1)
      call b%fs%get_viscosity(vf=b%vf)
         
      ! Perform sub-iterations
      do while (b%time%it.le.b%time%itmax)
            
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
            
         ! Increment sub-iteration counter
         b%time%it=b%time%it+1
            
      end do
         
      ! Recompute interpolated velocity and divergence
      call b%fs%interp_vel(b%Ui,b%Vi,b%Wi)
      call b%fs%get_div()

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
            np=0
            do k=b%vf%cfg%kmin_,b%vf%cfg%kmax_
               do j=b%vf%cfg%jmin_,b%vf%cfg%jmax_
                  do i=b%vf%cfg%imin_,b%vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVerticesb%(b%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; b%smesh%var(1,np)=real(getNumberOfPlanes(b%vf%liquid_gas_interface(i,j,k)),WP)
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
      call b%fs%get_max()
      call b%vf%get_max()
      call b%mfile%write()
      call b%cflfile%write()
         
      ! After we're done clip all VOF at the exit area and along the sides - hopefully nothing's left
      clip_vof: block
         integer :: i,j,k
         do k=b%fs%cfg%kmino_,b%fs%cfg%kmaxo_
            do j=b%fs%cfg%jmino_,b%fs%cfg%jmaxo_
               do i=b%fs%cfg%imino_,b%fs%cfg%imaxo_
                  if (i.ge.b%vf%cfg%imax-5) b%vf%VF(i,j,k)=0.0_WP
                  if (j.ge.b%vf%cfg%jmax-5) b%vf%VF(i,j,k)=0.0_WP
                  if (j.le.b%vf%cfg%jmin+5) b%vf%VF(i,j,k)=0.0_WP
                  if (k.ge.b%vf%cfg%kmax-5) b%vf%VF(i,j,k)=0.0_WP
                  if (k.le.b%vf%cfg%kmin+5) b%vf%VF(i,j,k)=0.0_WP
               end do
            end do
         end do
      end block clip_vof
      
      
   end subroutine step
      
   
   !> Finalize b2 simulation
   subroutine final
      implicit none
      class(block2), intent(inout) :: b
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(b%resU,b%resV,b%resW,b%Ui,b%Vi,b%Wi,b%SR,b%SRmag)
      
   end subroutine final
   
   
end module block2_class
