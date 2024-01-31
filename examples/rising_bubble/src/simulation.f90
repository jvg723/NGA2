!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,            only: WP
   use geometry,             only: cfg,moving_domain
   use hypre_str_class,      only: hypre_str
   use ddadi_class,          only: ddadi
   use tpns_class,           only: tpns
   use vfs_class,            only: vfs
   use tpviscoelastic_class, only: tpviscoelastic
   use timetracker_class,    only: timetracker
   use ensight_class,        only: ensight
   use surfmesh_class,       only: surfmesh
   use event_class,          only: event
   use monitor_class,        only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),      public :: ps
   type(ddadi),          public :: vs
   type(tpns),           public :: fs
   type(vfs),            public :: vf
   type(tpviscoelastic), public :: ve
   type(timetracker),    public :: time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,bubblefile,scfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp,SR,Eigenvalues,Conf
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU,Eigenvectors
   
   !> Problem definition
   real(WP), dimension(3) :: center,gravity
   real(WP) :: volume,radius,Ycent,Vrise
   real(WP) :: Vin,Vin_old,Vrise_ref,Ycent_ref,G,ti
   
contains


   !> Function that defines a level set function for a rising bubble problem
   function levelset_rising_bubble(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the bubble
      G=-radius+sqrt(sum((xyz-center)**2))
   end function levelset_rising_bubble
   
   
   !> Function that localizes the y+ side of the domain
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator
   
   
   !> Function that localizes the y- side of the domain
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


      !> Function that localizes y- boundary for scalar
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc


   !> Function that localizes y+ boundary for scalar
   function yp_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator_sc 
   
   
   !> Routine that computes rise velocity
   subroutine rise_vel()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i,j,k,ierr
      real(WP) :: myYcent,myVrise,myvol,bubble_vol
      myVrise=0.0_WP
      myvol=0.0_WP
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               myYcent=myYcent+vf%cfg%ym(j)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myVrise=myVrise+Vi(i,j,k)*(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
               myvol=myvol+(1.0_WP-vf%VF(i,j,k))*cfg%vol(i,j,k)
            end do
         end do
      end do
      call MPI_ALLREDUCE(myYcent,Ycent     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myVrise,Vrise     ,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(myvol  ,bubble_vol,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Ycent=Ycent/bubble_vol
      Vrise=Vrise/bubble_vol
   end subroutine
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(SCtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(Conf(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
         allocate(SR   (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Eigenvalues(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3))
         allocate(Eigenvectors(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3,1:3))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: elvira,VFhi,VFlo
         use mathtools, only: Pi
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=elvira,name='VOF',store_detailed_flux=.true.)
         ! Initialize a bubble
         call param_read('Bubble position',center,default=[0.0_WP,0.0_WP,0.0_WP])
         call param_read('Bubble volume',radius)
         radius=(radius*3.0_WP/(4.0_WP*Pi))**(1.0_WP/3.0_WP)
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_rising_bubble,0.0_WP,amr_ref_lvl)
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
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
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
      

      ! Prepare PID controller if domain is moving
      if (moving_domain) then
         prepare_controller: block
            ! Store target data
            Ycent_ref=center(2)
            Vrise_ref=0.0_WP
            ! Controller parameters
            G=0.5_WP
            ti=time%dtmax
         end block prepare_controller
      end if
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2,pcg_pfmg,gmres_pfmg
         use tpns_class,      only: clipped_neumann,dirichlet
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Assign acceleration of gravity
         call param_read('Gravity',gravity); fs%gravity=gravity
         ! Dirichlet inflow at the top and clipped Neumann outflow at the bottom
         if (moving_domain) then
            call fs%add_bcond(name='inflow' ,type=dirichlet      ,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
            call fs%add_bcond(name='outflow',type=clipped_neumann,face='y',dir=-1,canCorrect=.true. ,locator=ym_locator)
         end if
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=gmres_pfmg,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Zero inflow velocity
         Vin=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create a viscoleastic model with log conformation stablization method
      create_viscoelastic: block
         use tpviscoelastic_class, only: eptt
         use tpscalar_class,       only: bcond,neumann
         type(bcond), pointer :: mybc
         integer :: i,j,k
         ! Create viscoelastic model solver
         call ve%init(cfg=cfg,phase=0,model=eptt,name='viscoelastic')
         ! Relaxation time for polymer
         call param_read('Polymer relaxation time',ve%trelax)
         ! Polymer viscosity
         call param_read('Polymer viscosity',ve%visc_p)
         ! Extensional viscosity parameter (ePTT)
         call param_read('Extensional viscosity parameter',ve%elongvisc)
         ! Affine parameter (ePTT)
         call param_read('Extensional viscosity parameter',ve%affinecoeff)
         ! Apply boundary conditions
         if (moving_domain) then
            call ve%add_bcond(name='yp_sc',type=neumann,locator=yp_locator_sc,dir='yp')
            call ve%add_bcond(name='ym_sc',type=neumann,locator=ym_locator_sc,dir='ym')
         end if
         ! Setup without an implicit solver
         call ve%setup()
         ! Initialize C scalar fields (C=I -> 0 intital value in ln(C) feild) 
         Conf=0.0_WP
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  if (vf%VF(i,j,k).gt.0.0_WP) then
                     Conf(i,j,k,1)=1.0_WP  !< Cxx
                     Conf(i,j,k,4)=1.0_WP  !< Cyy
                     Conf(i,j,k,6)=1.0_WP  !< Czz
                     ! ve%SC(i,j,k,1)=1.0_WP  !< Cxx
                     ! ve%SC(i,j,k,4)=1.0_WP  !< Cyy
                     ! ve%SC(i,j,k,6)=1.0_WP  !< Czz
                  end if
               end do
            end do
         end do
         ! Get eigenvalues and eigenvectors
         call ve%get_eigensystem(Eigenvalues,Eigenvectors)
         ! Apply boundary conditions
         call ve%apply_bcond(time%t,time%dt)
      end block create_viscoelastic
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh
      

      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='RisingBubble')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         do nsc=1,ve%nscalar
            call ens_out%add_scalar(trim(ve%SCname(nsc)),ve%SC(:,:,:,nsc))
         end do
         call ens_out%add_scalar('recCxx',Conf(:,:,:,1))
         call ens_out%add_scalar('recCxy',Conf(:,:,:,2))
         call ens_out%add_scalar('recCxz',Conf(:,:,:,3))
         call ens_out%add_scalar('recCyy',Conf(:,:,:,4))
         call ens_out%add_scalar('recCyz',Conf(:,:,:,5))
         call ens_out%add_scalar('recCzz',Conf(:,:,:,6))
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         integer :: nsc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
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
         ! Create bubble monitor
         bubblefile=monitor(fs%cfg%amRoot,'bubble')
         call bubblefile%add_column(time%n,'Timestep number')
         call bubblefile%add_column(time%t,'Time')
         call bubblefile%add_column(Ycent,'Y centroid')
         call bubblefile%add_column(Vrise,'Rise velocity')
         call bubblefile%add_column(Vin,'Inflow velocity')
         call bubblefile%write()
         ! Create scalar monitor
         scfile=monitor(ve%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do nsc=1,ve%nscalar
            call scfile%add_column(ve%SCmin(nsc),trim(ve%SCname(nsc))//'_min')
            call scfile%add_column(ve%SCmax(nsc),trim(ve%SCname(nsc))//'_max')
         end do
         call scfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      use tpns_class, only: harmonic_visc
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Control inflow condition
         if (moving_domain) then
            control_inflow: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer  :: n,i,j,k
               ! Get new inflow velocity
               Vin_old=Vin
               Vin=G*((Vrise_ref-Vrise)+(Ycent_ref-Ycent)/ti)
               ! Apply inflow at top
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%V(i,j,k)=Vin
               end do
            end block control_inflow
         end if

         ! Calculate grad(U)
         call fs%get_gradU(gradU)

         ! Calcualte SR tensor
         call fs%get_strainrate(SR)

         ! Transport our liquid conformation tensor using log conformation
         advance_scalar: block
            integer :: i,j,k,nsc
            ! Add source terms for constitutive model
            call ve%get_CgradU(Eigenvalues,Eigenvectors,gradU,SCtmp,SR); resSC=SCtmp
            call ve%get_relax(Eigenvalues,Eigenvectors,SCtmp);           resSC=resSC+SCtmp
            ! call ve%get_CgradU(gradU,SCtmp,SR); resSC=SCtmp
            ! call ve%get_relax(SCtmp,time%dt);   resSC=resSC+SCtmp
            ve%SC=ve%SC+time%dt*resSC
            call ve%apply_bcond(time%t,time%dt)
            ve%SCold=ve%SC
            ! Explicit calculation of dSC/dt from scalar equation
            call ve%get_dSCdt(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,VFold=vf%VFold,VF=vf%VF,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)
            ! Update our scalars
            do nsc=1,ve%nscalar
               where (ve%mask.eq.0.and.vf%VF.ne.0.0_WP) ve%SC(:,:,:,nsc)=(vf%VFold*ve%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/vf%VF
               where (vf%VF.eq.0.0_WP) ve%SC(:,:,:,nsc)=0.0_WP
            end do
            ! Apply boundary conditions
            call ve%apply_bcond(time%t,time%dt)
         end block advance_scalar
         ! Update eigenvalues and eigenvectors
         call ve%get_eigensystem(Eigenvalues,Eigenvectors)
         ! Build conformation tensor from its diagonalized state
         reconstruct_conformation: block
            integer :: i,j,k,nsc
            Conf=0.0_WP
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     if (vf%VF(i,j,k).gt.0.0_WP) then
                        ! Skip non-solved cells
                        if (ve%mask(i,j,k).ne.0) cycle
                        ! Reconstruct conformation tensor (C=R*exp(ln(Lambda))*R^T={{Cxx,Cxy,Cxz},{Cxy,Cyy,Cyz},{Cxz,Cyz,Czz}})
                        !>xx tensor component
                        Conf(i,j,k,1)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,1,1)**2                     +Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,1,2)**2+Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,1,3)**2
                        !>xy tensor component
                        Conf(i,j,k,2)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)+Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)+Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)
                        !>xz tensor component
                        Conf(i,j,k,3)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)+Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)+Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)
                        !>yy tensor component
                        Conf(i,j,k,4)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,2,1)**2                     +Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,2,2)**2+Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,2,3)**2
                        !>yz tensor component
                        Conf(i,j,k,5)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)+Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,2,2)*Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)
                        !>zz tensor component
                        Conf(i,j,k,6)=Eigenvalues(i,j,k,1)*Eigenvectors(i,j,k,3,1)**2                     +Eigenvalues(i,j,k,2)*Eigenvectors(i,j,k,3,2)**2+Eigenvalues(i,j,k,3)*Eigenvectors(i,j,k,3,3)**2
                     end if
                  end do
               end do
            end do
         end block reconstruct_conformation
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
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
            
            ! Add momentum source terms - adjust gravity if accelerating frame of reference
            if (moving_domain) fs%gravity(2)=gravity(2)+(Vin-Vin_old)/time%dt
            call fs%addsrc_gravity(resU,resV,resW)

            ! Add polymer stress term
            polymer_stress: block
               use tpviscoelastic_class, only: fenecr,oldroydb,lptt,eptt,fenep
               integer :: i,j,k,nsc
               real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
               real(WP), dimension(:,:,:,:), allocatable :: stress
               real(WP) :: coeff
               ! Allocate work arrays
               allocate(stress(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
               allocate(Txy   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tyz   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               allocate(Tzx   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
               ! Calculate polymer stress for a given model
               select case (ve%model)
               case (oldroydb,fenecr,fenep)
               ! Calculate the polymer relaxation
               ! call ve%get_relax(stress,time%dt)
               call ve%get_relax(Eigenvalues,Eigenvectors,stress)
               ! Build liquid stress tensor
               do nsc=1,6
                  ! stress(:,:,:,nsc)=-nn%visc(:,:,:)*vf%VF*stress(:,:,:,nsc)
                  stress(:,:,:,nsc)=-ve%visc_p*vf%VF*stress(:,:,:,nsc)
               end do
               case (eptt,lptt)
                  ! coeff=ve%visc_p/(ve%trelax*(1-ve%affinecoeff))
                  do k=cfg%kmino_,cfg%kmaxo_
                     do j=cfg%jmino_,cfg%jmaxo_
                        do i=cfg%imino_,cfg%imaxo_
                           ! stress(i,j,k,1)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,1)-1.0_WP) !> xx tensor component
                           ! stress(i,j,k,2)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,2)-0.0_WP) !> xy tensor component
                           ! stress(i,j,k,3)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,3)-0.0_WP) !> xz tensor component
                           ! stress(i,j,k,4)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,4)-1.0_WP) !> yy tensor component
                           ! stress(i,j,k,5)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,5)-0.0_WP) !> yz tensor component
                           ! stress(i,j,k,6)=vf%VF(i,j,k)*coeff*(ve%SC(i,j,k,6)-1.0_WP) !> zz tensor component
                           stress(i,j,k,1)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,1,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                           stress(i,j,k,2)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                           stress(i,j,k,3)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                           stress(i,j,k,4)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,2,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                           stress(i,j,k,5)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                           stress(i,j,k,6)=vf%VF(i,j,k)*(-ve%visc_p/(ve%trelax*(1.00_WP-ve%affinecoeff)))*(Eigenvectors(i,j,k,3,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
                        end do
                     end do
                  end do
               end select
               ! Interpolate tensor components to cell edges
               do k=cfg%kmin_,cfg%kmax_+1
                  do j=cfg%jmin_,cfg%jmax_+1
                     do i=cfg%imin_,cfg%imax_+1
                        Txy(i,j,k)=sum(fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
                        Tyz(i,j,k)=sum(fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
                        Tzx(i,j,k)=sum(fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
                     end do
                  end do
               end do
               ! Add divergence of stress to residual
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
                        &                                                +sum(fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
                        &                                                +sum(fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
                        if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
                        &                                                +sum(fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
                        &                                                +sum(fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
                        if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
                        &                                                +sum(fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &                  
                        &                                                +sum(fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
                     end do
                  end do
               end do
               ! Clean up
               deallocate(stress,Txy,Tyz,Tzx)
            end block polymer_stress
            
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
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call rise_vel()
         call mfile%write()
         call cflfile%write()
         call bubblefile%write()
         call scfile%write()
         
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
      deallocate(resSC,SCtmp,SR,gradU,Eigenvalues,Eigenvectors,Conf)
      
   end subroutine simulation_final
   
   
end module simulation