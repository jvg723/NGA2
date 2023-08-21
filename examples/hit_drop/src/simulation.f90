!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,nx,Lx
   use fft3d_class,       only: fft3d
   use tpns_class,        only: tpns
   use tracer_class,      only: tracer
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use irl_fortran_interface
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   type(fft3d),       public :: ps
   type(tpns),        public :: fs
   type(timetracker), public :: time
   type(vfs),         public :: vf
   type(tracer),      public :: pt
 
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   type(surfmesh) :: smesh
  
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lptfile,hitfile,cvgfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: visc,rho,meanU,meanV,meanW
   real(WP) :: Urms0,L0,TKE0,EPS0,ReT_max,ReL_max,eta_min
   real(WP) :: ReL_tgt,ReT_tgt,eta_tgt,taueta_tgt,We_tgt,inject_time
   real(WP) :: TKE,URMS
   real(WP) :: tauinf,G,Gdtau,Gdtaui,dx
   real(WP), dimension(3) :: center
   real(WP) :: radius
   integer :: npart
   logical :: normal_only

   !> For monitoring
   real(WP) :: EPS
   real(WP) :: Re_L,Re_lambda
   real(WP) :: eta,ell
   real(WP) :: dx_eta,ell_Lx,Re_ratio,eps_ratio,tke_ratio,nondtime
   

contains
   
   
   !> Function that defines a level set function for a cylinder
   function levelset_sphere(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=radius-sqrt(sum((xyz-center)**2))
   end function levelset_sphere

   !> Compute turbulence stats
   subroutine compute_stats()
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
      use parallel, only: MPI_REAL_WP
      real(WP) :: myTKE,myEPS
      integer :: i,j,k,ierr
      
      ! Compute mean velocities
      call fs%cfg%integrate(A=Ui,integral=meanU); meanU=meanU/fs%cfg%vol_total
      call fs%cfg%integrate(A=Vi,integral=meanV); meanV=meanV/fs%cfg%vol_total
      call fs%cfg%integrate(A=Wi,integral=meanW); meanW=meanW/fs%cfg%vol_total
      
      ! Compute strainrate and grad(u)
      call fs%get_strainrate(SR=SR)
      call fs%get_gradu(gradu=gradU)
      
      ! Compute current TKE and dissipation rate
      myTKE=0.0_WP
      myEPS=0.0_WP
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               myTKE=myTKE+0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
               myEPS=myEPS+2.0_WP*fs%cfg%vol(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
            end do
         end do
      end do

      call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE=TKE/fs%cfg%vol_total
      call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPS=EPS*visc/fs%cfg%vol_total
      
      ! Compute standard parameters for HIT
      Urms=sqrt(2.0_WP/3.0_WP*TKE)
      Re_L=TKE**2.0_WP/(visc*EPS)
      Re_lambda=sqrt(20.0_WP*Re_L/3.0_WP)
      eta=(visc**3.0_WP/EPS)**0.25_WP
      ell=(2.0_WP*TKE/3.0_WP)**1.5_WP/EPS
 
      ! Some more useful info
      nondtime =time%t/tauinf
      dx_eta   =dx/eta
      eps_ratio=EPS/EPS0
      tke_ratio=TKE/TKE0
      ell_Lx   =ell/Lx
      Re_ratio =Re_lambda/ReL_tgt
      
   end subroutine compute_stats
   
   
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
         allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      ! Prepare initial velocity field
      initialize_hit: block
         use param,     only: param_exists
         use mathtools, only: Pi
         use messager,  only: die
         real(WP) :: max_forcing_estimate,L
         ! Mesh spacing
         dx=Lx/real(nx,WP)
         ! We assume that the integral lengthscale is Lx/5
         L0=Lx/5.0_WP
         ! The density is set to 1
         rho=1.0_WP
         ! The RMS of velocity fluctation is set to 1
         Urms0=1.0_WP
         ! Dissipation rate
         EPS0=Urms0**3.0_WP/L0
         ! Turbulent kinetic energy
         TKE0=1.5_WP*Urms0**2.0_WP
         ! Eddy turnover time
         tauinf=L0/Urms0
         ! Initialize target Re_lambda
         call param_read('Reynolds number',ReL_tgt,default=0.0_WP)
         ! Mininum eta and max Re_lambda based on mesh-spacing
         eta_min=1.5_WP*dx/Pi ! From Pope, p.347
         L=TKE0**1.5_WP/EPS0 ! From Pope, p.200
         ReT_max=(L/eta_min)**(4.0_WP/3.0_WP) ! From Pope, p.200
         ReL_max=sqrt(20.0_WP/3.0_WP*ReT_max) ! From Pope, p.200
         ! We can define now targets based on the requested Re or the maximum achievable
         ReL_tgt=max(ReL_tgt,ReL_max)
         ReT_tgt=3.0_WP/20.0_WP*ReL_tgt**2.0_WP ! From Pope, p.200
         eta_tgt=L/ReT_tgt**(3.0_WP/4.0_WP) ! From Pope, p.200
         ! Having a target Re_lambda, we can now set the viscosity
         visc=rho*TKE0**2.0_WP/(ReT_tgt*EPS0) ! From Pope, p.200
         ! Calculate other target quantities
         taueta_tgt=sqrt(visc/EPS0)
         ! Read in forcing parameter (we need dt<tauinf/forcing)
         max_forcing_estimate=3.0_WP*tauinf*Urms0/dx
         call param_read('Forcing constant',G,default=max_forcing_estimate)
         Gdtau =G/tauinf
         Gdtaui=1.0_WP/Gdtau
         ! Droplet injection time (in terms of eddy turnover)
         call param_read('Droplet injection time',inject_time,default=1.0_WP); inject_time=inject_time*tauinf
      end block initialize_hit

      ! Initialize our VOF solver
      create_vof: block
         use vfs_class,only: r2p,lvira,VFhi,VFlo
         use mpi_f08,  only: MPI_WTIME
         use string,   only: str_medium,lowercase
         integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
         character(len=str_medium) :: read_curvature_method
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=5
         real(WP) :: start, finish
         ! Create a VOF solver with r2p reconstruction
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize two droplets
         call param_read('Droplet diameter',radius); radius=0.5_WP*radius
         call param_read('Droplet position',center)      
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
      end block create_vof

      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         fs%visc_l = visc; fs%visc_g = fs%visc_l
         ! Assign constant density to each phase
         fs%rho_l = rho; fs%rho_g = fs%rho_l
         ! Read in surface tension coefficient
         call param_read('Weber number',We_tgt)
         fs%sigma=2.0_WP*rho*EPS0**(2.0_WP/3.0_WP)*(2.0_WP*radius)**(5.0_WP/3.0_WP)/We_tgt ! Based on Risso and Fabre (1998)
         ! Prepare and configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='curv'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         call add_surfgrid_variable(smesh,1,vf%curv)      
      end block create_smesh

      ! Prepare initial velocity field
      initialize_velocity: block
         use random,    only: random_normal
         use param,     only: param_exists
         use mathtools, only: Pi
         use messager,  only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         integer :: i,j,k
         ! Print out some expected turbulence properties
         if (fs%cfg%amRoot) then
            write(output_unit,'("Target parameters:")')
            write(output_unit,'("[HIT setup]  => visc              =",es12.5)') visc
            write(output_unit,'("[HIT setup]  => Urms              =",es12.5)') Urms0
            write(output_unit,'("[HIT setup]  => Re_lambda         =",es12.5)') ReL_tgt
            write(output_unit,'("[HIT setup]  => Re_turb           =",es12.5)') ReT_tgt
            write(output_unit,'("[HIT setup]  => Kolmogorov Lscale =",es12.5)') eta_tgt
            write(output_unit,'("[HIT setup]  => Kolmogorov Tscale =",es12.5)') taueta_tgt
            write(output_unit,'("[HIT setup]  => Epsilon           =",es12.5)') EPS0
            write(output_unit,'("[HIT setup]  => Eddyturnover time =",es12.5)') tauinf
            write(output_unit,'("[Drop setup] => We                =",es12.5)') We_tgt
            write(output_unit,'("[Drop setup] => sigma             =",es12.5)') fs%sigma
         end if
         ! Gaussian initial field
         ! call random_init(.true., .true.)
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
         fs%U=fs%U-meanU
         fs%V=fs%V-meanV
         fs%W=fs%W-meanW
         ! Project to ensure divergence-free
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div*rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU/rho
         fs%V=fs%V-time%dt*resV/rho
         fs%W=fs%W-time%dt*resW/rho
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Compute turbulence stats
         call compute_stats()
      end block initialize_velocity
      
      ! Initialize our tracer particle tracker
      initialize_pt: block
         ! Get number of particles
         call param_read('Number of tracers',npart)
         ! Remove tangential component
         call param_read('Remove tangential motion',normal_only,default=.true.)
         ! Create solver
         pt=tracer(cfg=cfg,name='PT')
         ! Initialize with zero particles
         call pt%resize(0)
      end block initialize_pt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=0,nvec=2,name='tracers')
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='acceleration'
         call pt%update_partmesh(pmesh)
         do i=1,pt%np_
            pmesh%vec(:,1,i)=pt%p(i)%vel
            pmesh%vec(:,2,i)=pt%p(i)%acc
         end do
      end block create_pmesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='HIT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('signeddist',vf%G)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('vofplic',smesh) 
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
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
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create hit monitor
         hitfile=monitor(fs%cfg%amRoot,'hit')
         call hitfile%add_column(time%n,'Timestep number')
         call hitfile%add_column(time%t,'Time')
         call hitfile%add_column(Re_L,'Re_L')
         call hitfile%add_column(Re_lambda,'Re_lambda')
         call hitfile%add_column(eta,'eta')
         call hitfile%add_column(TKE,'TKE')
         call hitfile%add_column(URMS,'Urms')
         call hitfile%add_column(EPS,'EPS')
         call hitfile%add_column(ell,'L')
         call hitfile%write()
         ! Create hit convergence monitor
         cvgfile=monitor(fs%cfg%amRoot,'convergence')
         call cvgfile%add_column(time%n,'Timestep number')
         call cvgfile%add_column(time%t,'Time')
         call cvgfile%add_column(nondtime,'Time/t_int')
         call cvgfile%add_column(Re_ratio,'Re_ratio')
         call cvgfile%add_column(eps_ratio,'EPS_ratio')
         call cvgfile%add_column(tke_ratio,'TKE_ratio')
         call cvgfile%add_column(dx_eta,'dx/eta')
         call cvgfile%add_column(ell_Lx,'ell/Lx')
         call cvgfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      implicit none
      logical :: droplet_inserted=.false.

      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Insert droplet once
         if (.not.droplet_inserted.and.time%t.gt.inject_time) then
            ! Initialize VOF field of droplet
            call insert_drop(vf=vf)
            ! Seed interface with tracers
            insert_tracers: block
               use mpi_f08, only: MPI_INTEGER8,MPI_MAX
               integer :: i,np0_,ierr
               integer(kind=8) :: maxid_,maxid  !< Keep track of maximum particle id
               real :: hk,thetak,phik,realn
               ! Initial number of particles
               np0_=pt%np_
               ! Determine id to assign to particle
               maxid_=0
               do i=1,pt%np_
                  maxid_=max(maxid_,pt%p(i)%id)
               end do
               call MPI_ALLREDUCE(maxid_,maxid,1,MPI_INTEGER8,MPI_MAX,cfg%comm,ierr)
      
               ! Add new particles
               if (cfg%amRoot) then
                  ! Create space for new particle
                  call pt%resize(np0_+npart)
                  ! Initialize parameters
                  phik=0.0_WP;realn=REAL(npart,WP)
                  do i=1,npart
                     ! Helicoidal seeding based on Saff and Kuijlaars (1997)
                     hk=-1.0_WP+2.0_WP*(REAL(i,WP)-1.0_WP)/(realn-1.0_WP); thetak=acos(hk)
                     if (i.eq.1.or.i.eq.npart) then
                        phik=0.0_WP
                     else
                        phik=phik+3.6_WP/sqrt(realn)/sqrt(1.0_WP-hk**2.0_WP)
                     end if
                     ! Set particle ID
                     pt%p(np0_+i)%id=maxid+int(i,8)
                     ! Seed particle on sphere
                     pt%p(np0_+i)%pos(1)=center(1)+radius*sin(thetak)*cos(phik)
                     pt%p(np0_+i)%pos(2)=center(2)+radius*sin(thetak)*sin(phik)
                     pt%p(np0_+i)%pos(3)=center(3)+radius*cos(thetak)
                     ! Localize the particle
                     pt%p(np0_+i)%ind=cfg%get_ijk_global(pt%p(np0_+i)%pos,[cfg%imin,cfg%jmin,cfg%kmin])
                     ! Make it an "official" particle
                     pt%p(np0_+i)%flag=0
                  end do
               end if
               ! Communicate particles
               call pt%sync()
            end block insert_tracers
            droplet_inserted=.true.
         end if

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

         ! Advance and project tracer particles on the interface
         if (.not.normal_only) then
            call pt%advance(dt=time%dtmid,U=fs%U,V=fs%V,W=fs%W)
         else
            advance_tracers_normally: block
               use mathtools, only: Pi,normalize
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
               use parallel, only: MPI_REAL_WP
               real(WP), dimension(3) :: old_vel,old_pos,mymassvel,massvel
               real(WP), dimension(3) :: normal
               real(WP) :: mymass,mass,norm
               integer :: i,j,k,ierr

               ! Interpolate velocity at cell center
               call fs%interp_vel(Ui,Vi,Wi)

               ! Get gradient of distance function to the interface
               call fs%get_pgrad(vf%G,resU,resV,resW)

               ! Compute center of mass velocity
               mymassvel=[0.0_WP,0.0_WP,0.0_WP]; mymass=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        mymass=mymass+vf%VF(i,j,k)
                        mymassvel=mymassvel+vf%VF(i,j,k)*[Ui(i,j,k),Vi(i,j,k),Wi(i,j,k)]
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(mymass,mass,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr);
               call MPI_ALLREDUCE(mymassvel,massvel,3,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); massvel=massvel/mass
               
               ! Advance the equations
               do i=1,pt%np_
                  ! Avoid particles with id=0
                  if (pt%p(i)%id.eq.0) cycle
                  ! Remember the particle
                  old_vel=pt%p(i)%vel
                  old_pos=pt%p(i)%pos
                  ! Interpolate fluid velocity to particle location
                  normal=normalize([cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resU,bc='n'),&
                  &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resV,bc='n'),&
                  &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resW,bc='n')])
                  pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
                  pt%p(i)%vel=dot_product(pt%p(i)%vel-massvel,normal)*normal+massvel
                  ! Advance with Euler prediction
                  pt%p(i)%pos=old_pos+0.5_WP*time%dtmid*pt%p(i)%vel
                  ! Correct with midpoint rule
                  normal=normalize([cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resU,bc='n'),&
                  &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resV,bc='n'),&
                  &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resW,bc='n')])
                  pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
                  pt%p(i)%vel=dot_product(pt%p(i)%vel-massvel,normal)*normal+massvel
                  pt%p(i)%pos=old_pos+time%dtmid*pt%p(i)%vel
                  ! Relocalize the particle
                  pt%p(i)%ind=cfg%get_ijk_global(pt%p(i)%pos,pt%p(i)%ind)
                  ! Correct the position to take into account periodicity
                  if (cfg%xper) pt%p(i)%pos(1)=cfg%x(cfg%imin)+modulo(pt%p(i)%pos(1)-cfg%x(cfg%imin),cfg%xL)
                  if (cfg%yper) pt%p(i)%pos(2)=cfg%y(cfg%jmin)+modulo(pt%p(i)%pos(2)-cfg%y(cfg%jmin),cfg%yL)
                  if (cfg%zper) pt%p(i)%pos(3)=cfg%z(cfg%kmin)+modulo(pt%p(i)%pos(3)-cfg%z(cfg%kmin),cfg%zL)
                  ! Relocalize the particle
                  pt%p(i)%ind=cfg%get_ijk_global(pt%p(i)%pos,pt%p(i)%ind)
                  ! Interpolate fluid quantities to particle location
                  pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
                  ! Update acceleration
                  pt%p(i)%acc=(pt%p(i)%vel-old_vel)/time%dtmid
               end do
               ! Communicate particles
               call pt%sync()
               ! Log/screen output
               logging: block
                  use, intrinsic :: iso_fortran_env, only: output_unit
                  use param,    only: verbose
                  use messager, only: log
                  use string,   only: str_long
                  character(len=str_long) :: message
                  if (cfg%amRoot) then
                     write(message,'("Tracer solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(pt%name),trim(cfg%name),pt%np
                     if (verbose.gt.1) write(output_unit,'(a)') trim(message)
                     if (verbose.gt.0) call log(message)
                  end if
               end block logging
            end block advance_tracers_normally
         end if
         call project_tracers(vf=vf,pt=pt)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Add linear forcing term based on Bassenne et al. (2016)
            linear_forcing: block
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
               use parallel, only: MPI_REAL_WP
               real(WP) :: myTKE,A,myEPSp,EPSp
               integer :: i,j,k,ierr
               ! Calculate mean velocity
               call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
               ! Calculate TKE and pseudo-EPS
               call fs%interp_vel(Ui,Vi,Wi)
               call fs%get_gradu(gradu=gradU)
               myTKE=0.0_WP; myEPSp=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        myTKE =myTKE +0.5_WP*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
                        myEPSp=myEPSp+fs%cfg%vol(i,j,k)*visc*(gradU(1,1,i,j,k)**2+gradU(1,2,i,j,k)**2+gradU(1,3,i,j,k)**2+&
                        &                                               gradU(2,1,i,j,k)**2+gradU(2,2,i,j,k)**2+gradU(2,3,i,j,k)**2+&
                        &                                               gradU(3,1,i,j,k)**2+gradU(3,2,i,j,k)**2+gradU(3,3,i,j,k)**2)
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myTKE ,TKE ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE =TKE /fs%cfg%vol_total
               call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPSp=EPSp/fs%cfg%vol_total/rho
               A=(EPSp-Gdtau*(TKE-TKE0))/(2.0_WP*TKE)*rho
               resU=resU+time%dt*(fs%U-meanU)*A
               resV=resV+time%dt*(fs%V-meanV)*A
               resW=resW+time%dt*(fs%W-meanW)*A
            end block linear_forcing
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            ! call fs%update_laplacian()
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
            ! Update surfmesh object
            update_smesh: block
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Populate curv variable
               call add_surfgrid_variable(smesh,1,vf%curv)      
            end block update_smesh           
            update_pmesh: block
              integer :: i
              call pt%update_partmesh(pmesh)
              do i=1,pt%np_
                  pmesh%vec(:,1,i)=pt%p(i)%vel
                  pmesh%vec(:,2,i)=pt%p(i)%acc
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call compute_stats()
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call hitfile%write()
         call cvgfile%write()
         
      end do
      
   end subroutine simulation_run

   !> Make a surface scalar variable
   subroutine add_surfgrid_variable(smesh,var_index,A)
      use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
      implicit none
      class(surfmesh), intent(inout) :: smesh
      integer, intent(in) :: var_index
      real(WP), dimension(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_), intent(in) :: A 
      integer :: i,j,k,shape,np,nplane
      
      ! Fill out arrays
      if (smesh%nPoly.gt.0) then
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     shape=getNumberOfVertices(vf%interface_polygon(nplane,i,j,k))
                     if (shape.gt.0) then
                        ! Increment polygon counter
                        np=np+1
                        ! Set nplane variable
                        smesh%var(var_index,np)=A(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      else
         smesh%var(var_index,1)=1
      end if      
      
   end subroutine add_surfgrid_variable
   
   ! Initialize our VOF field
   subroutine insert_drop(vf)
      use mms_geom, only: cube_refine_vol
      use vfs_class,only: r2p,lvira,VFhi,VFlo
      use mpi_f08,  only: MPI_WTIME
      use string,   only: str_medium,lowercase
      implicit none
      class(vfs), intent(inout) :: vf
      integer :: i,j,k,n,si,sj,sk,curvature_method,stencil_size,hf_backup_method
      character(len=str_medium) :: read_curvature_method
      real(WP), dimension(3,8) :: cube_vertex
      real(WP), dimension(3) :: v_cent,a_cent
      real(WP) :: vol,area
      integer, parameter :: amr_ref_lvl=5
      real(WP) :: start, finish
      ! Initialize two droplets
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
               call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
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
   end subroutine insert_drop

   ! Project tracers on the interface
   subroutine project_tracers(vf,pt)
      use mathtools, only: normalize
      implicit none
      class(vfs), intent(inout) :: vf
      class(tracer), intent(inout) :: pt
      integer :: i, j, maxit
      real(WP) :: distx, disty, distz, dist, maxdist
      real(WP), dimension(3) :: oldvel, oldpos, normal
      ! Maximum projection iterations
      maxit=10
      ! Maximum distance from the interface wanted
      maxdist=1.0E-9_WP*dx
      ! Get gradient of distance function to the interface
      call fs%get_pgrad(vf%G,resU,resV,resW)
      ! Loop over local tracers
      do j=1,maxit
         do i=1,pt%np_
            ! Avoid particles with id=0
            if (pt%p(i)%id.eq.0) cycle
            oldpos=pt%p(i)%pos
            ! Interpolate distance function
            dist=cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=vf%G,bc='n')
            if (abs(dist).gt.maxdist) then
               ! Interpolate interface normal
               normal=normalize([cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resU,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resV,bc='n'),&
               &                 cfg%get_scalar(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),S=resW,bc='n')])
               ! Move tracer along distance function gradient
               pt%p(i)%pos=pt%p(i)%pos-dist*normal
               ! Correct the position to take into account periodicity
               pt%p(i)%pos(1)=cfg%x(cfg%imin)+modulo(pt%p(i)%pos(1)-cfg%x(cfg%imin),cfg%xL)
               pt%p(i)%pos(2)=cfg%y(cfg%jmin)+modulo(pt%p(i)%pos(2)-cfg%y(cfg%jmin),cfg%yL)
               pt%p(i)%pos(3)=cfg%z(cfg%kmin)+modulo(pt%p(i)%pos(3)-cfg%z(cfg%kmin),cfg%zL)
               ! Relocalize the tracer
               pt%p(i)%ind=cfg%get_ijk_global(pt%p(i)%pos,pt%p(i)%ind)
            end if
         end do
         ! Communicate tracers across procs
         call pt%sync()
      end do
      do i=1,pt%np_
         ! Interpolate fluid quantities to particle location
         pt%p(i)%vel=cfg%get_velocity(pos=pt%p(i)%pos,i0=pt%p(i)%ind(1),j0=pt%p(i)%ind(2),k0=pt%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
      end do
   end subroutine project_tracers

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradU)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
