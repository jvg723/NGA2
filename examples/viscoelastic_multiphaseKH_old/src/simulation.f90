!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,            only: WP
   use geometry,             only: cfg
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
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
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
   type(monitor) :: mfile,cflfile,scfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   
   !> Problem definition
   real(WP) :: Reg,Weg,r_visc,r_rho,r_vel,delta_l
   integer :: nwaveX,nwaveZ
   real(WP), dimension(:), allocatable :: wnumbX,wshiftX,wampX,wnumbZ,wshiftZ,wampZ

   !> Check for stabilization 
   logical :: stabilization 
   
contains
   
   
   !> Function that defines a level set function for a initial wavy interface
   function levelset_wavy(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      integer :: nX,nZ
      G=-xyz(2)
      do nX=1,nwaveX
         do nZ=1,nwaveZ
            G=G+wampX(nX)*cos(wnumbX(nX)*(xyz(1)-wshiftX(nX)))*wampZ(nZ)*cos(wnumbZ(nZ)*(xyz(3)-wshiftZ(nZ)))
         end do
      end do
   end function levelset_wavy
   
   
   !> Specialized subroutine that outputs wave amplitude information
   !subroutine postproc_data()
   !   use irl_fortran_interface
   !   use mathtools, only: Pi
   !   use string,    only: str_medium
   !   use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
   !   use parallel,  only: MPI_REAL_WP
   !   implicit none
   !   integer :: ierr,i,j,k,my_size
   !   real(WP) :: my_height
   !   real(WP), dimension(:), allocatable :: temp
   !   ! Calculate new amplitude
   !   grate=amp
   !   my_height=0.0_WP
   !   do k=vf%cfg%kmin_,vf%cfg%kmax_
   !      do i=vf%cfg%imin_,vf%cfg%imax_
   !         ! Find closest vertical column to center
   !         if (vf%cfg%x(i).le.0.5_WP*vf%cfg%xL.and.vf%cfg%x(i+1).gt.0.5_WP*vf%cfg%xL.and.vf%cfg%z(k).le.0.0_WP.and.vf%cfg%z(k+1).gt.0.0_WP) then
   !            ! Integrate height
   !            do j=vf%cfg%jmin_,vf%cfg%jmax_
   !               my_height=my_height+vf%VF(i,j,k)*vf%cfg%dy(j)
   !            end do
   !         end if
   !      end do
   !   end do
   !   call MPI_ALLREDUCE(my_height,amp,1,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr)
   !   amp=amp-0.5_WP*vf%cfg%yL
   !   ! Estimate growth rate
   !   if (time%t.gt.0.0_WP) then
   !      grate=(amp-grate)/time%dt
   !   else
   !      grate=0.0_WP
   !   end if
   !   ! Store time and amplitude series
   !   if (.not.allocated(all_time)) then
   !      my_size=0
   !   else
   !      my_size=size(all_time,dim=1)
   !   end if
   !   allocate(temp(my_size+1)); temp(1:my_size)=all_time; temp(my_size+1)=time%t; call MOVE_ALLOC(temp,all_time)
   !   allocate(temp(my_size+1)); temp(1:my_size)=all_amp ; temp(my_size+1)=amp   ; call MOVE_ALLOC(temp,all_amp )
   !end subroutine postproc_data
   
   
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
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         use vfs_class, only: plicnet,VFhi,VFlo,flux_storage
         use mathtools, only: twoPi
         use random,    only: random_uniform
         use parallel,  only: MPI_REAL_WP
         use mpi_f08
         integer :: i,j,k,n,si,sj,sk,ierr
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=plicnet,transport_method=flux_storage,name='VOF')
         !vf%cons_correct=.false.
         ! Prepare initialize interface parameters
         nwaveX=6
         allocate(wnumbX(nwaveX),wshiftX(nwaveX),wampX(nwaveX))
         wampX=0.3_WP/real(nwaveX,WP)
         wnumbX=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%xL
         if (cfg%amRoot) then
            do n=1,nwaveX
               wshiftX(n)=random_uniform(lo=-0.5_WP*cfg%xL,hi=+0.5_WP*cfg%xL)
            end do
         end if
         call MPI_BCAST(wshiftX,nwaveX,MPI_REAL_WP,0,cfg%comm,ierr)
         nwaveZ=6
         allocate(wnumbZ(nwaveZ),wshiftZ(nwaveZ),wampZ(nwaveZ))
         wampZ=0.3_WP/real(nwaveZ,WP)
         wnumbZ=[3.0_WP,4.0_WP,5.0_WP,6.0_WP,7.0_WP,8.0_WP]*twoPi/cfg%zL
         if (cfg%amRoot) then
            do n=1,nwaveZ
               wshiftZ(n)=random_uniform(lo=-0.5_WP*cfg%zL,hi=+0.5_WP*cfg%zL)
            end do
         end if
         call MPI_BCAST(wshiftZ,nwaveZ,MPI_REAL_WP,0,cfg%comm,ierr)
         ! Create the wavy interface
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_wavy,0.0_WP,amr_ref_lvl)
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
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         integer :: i,j,k
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Read in flow conditions
         call param_read('Gas Reynolds number',Reg); fs%visc_g=1.0_WP/(Reg+epsilon(Reg))
         call param_read('Viscosity ratio',r_visc); fs%visc_l=r_visc*fs%visc_g
         call param_read('Density ratio',r_rho); fs%rho_g=1.0_WP; fs%rho_l=r_rho*fs%rho_g
         call param_read('Gas Weber number',Weg); fs%sigma=1.0_WP/(Weg+epsilon(Weg))
         call param_read('Velocity ratio',r_vel); delta_l=r_visc*r_vel
         ! Configure pressure solver
			ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Set initial velocity field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%cfg%ym(j).le.0.0_WP) then
                     ! Use the liquid profile
                     fs%U(i,j,k)=r_vel*erf(fs%cfg%ym(j)/delta_l)
                  else
                     ! Use the gas profile
                     fs%U(i,j,k)=erf(fs%cfg%ym(j))
                  end if
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver

      ! Create a viscoleastic model with log conformation stablization method
      create_viscoelastic: block
         use tpviscoelastic_class, only: oldroydb
         use tpscalar_class,       only: bcond,neumann
         type(bcond), pointer :: mybc
         integer :: i,j,k
         ! Create viscoelastic model solver
         call ve%init(cfg=cfg,phase=0,model=oldroydb,name='viscoelastic')
         ! Relaxation time for polymer
         call param_read('Weissenberg Number',ve%trelax)
         ! Polymer viscosity
         call param_read('Polymer viscosity ratio',ve%visc_p); ve%visc_p=ve%visc_p*fs%visc_l
         ! Setup without an implicit solver
         call ve%setup()
         ! Check first if we use stabilization
         call param_read('Stabilization',stabilization,default=.false.)
         ! Initialize C scalar fields
         if (stabilization) then 
            !> Allocate storage fo eigenvalues and vectors
            allocate(ve%eigenval    (1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); ve%eigenval=0.0_WP
            allocate(ve%eigenvec(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); ve%eigenvec=0.0_WP
            !> Allocate storage for reconstructured C
            allocate(ve%SCrec   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6)); ve%SCrec=0.0_WP
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     if (vf%VF(i,j,k).gt.0.0_WP) then
                        ve%SCrec(i,j,k,1)=1.0_WP  !< Cxx
                        ve%SCrec(i,j,k,4)=1.0_WP  !< Cyy
                        ve%SCrec(i,j,k,6)=1.0_WP  !< Czz
                     end if
                  end do
               end do
            end do
            ! Get eigenvalues and eigenvectors
            call ve%get_eigensystem(vf%VF)
         end if
         ! Apply boundary conditions
         call ve%apply_bcond(time%t,time%dt)
      end block create_viscoelastic
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=7,name='plic')
         smesh%varname(1)='trC'
         smesh%varname(2)='Cxx'
         smesh%varname(3)='Cxy'
         smesh%varname(4)='Cxz'
         smesh%varname(5)='Cyy'
         smesh%varname(6)='Cyz'
         smesh%varname(7)='Czz'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Initalize variables to 0
         smesh%var(1,:)=0.0_WP
         smesh%var(2,:)=0.0_WP
         smesh%var(3,:)=0.0_WP
         smesh%var(4,:)=0.0_WP
         smesh%var(5,:)=0.0_WP
         smesh%var(6,:)=0.0_WP
         smesh%var(7,:)=0.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; 
                        smesh%var(1,np)=ve%SCrec(i,j,k,1)+ve%SCrec(i,j,k,4)+ve%SCrec(i,j,k,6)
                        smesh%var(2,np)=ve%SCrec(i,j,k,1)
                        smesh%var(3,np)=ve%SCrec(i,j,k,2)
                        smesh%var(4,np)=ve%SCrec(i,j,k,3)
                        smesh%var(5,np)=ve%SCrec(i,j,k,4)
                        smesh%var(6,np)=ve%SCrec(i,j,k,5)
                        smesh%var(7,np)=ve%SCrec(i,j,k,6)
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
         ens_out=ensight(cfg=cfg,name='VE_MPKH')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         if (stabilization) then
            do nsc=1,ve%nscalar
               call ens_out%add_scalar(trim(ve%SCname(nsc)),ve%SCrec(:,:,:,nsc))
            end do
         end if
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
         if (stabilization) then
            call ve%get_max_reconstructed(vf%VF)
         end if
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
         ! Create scalar monitor
         scfile=monitor(ve%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         if (stabilization) then
            do nsc=1,ve%nscalar
               call scfile%add_column(ve%SCrecmin(nsc),trim(ve%SCname(nsc))//'_min')
               call scfile%add_column(ve%SCrecmax(nsc),trim(ve%SCname(nsc))//'_max')
            end do
         end if
         call scfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())!.and.amp.lt.0.1_WP*vf%cfg%yL.and.time%t.lt.20.0_WP*tau)
         
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
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)

         ! Calculate grad(U)
         call fs%get_gradU(gradU)

         ! Transport our liquid conformation tensor using log conformation
         advance_scalar: block
            integer :: ierr
            integer :: i,j,k,nsc
            ! Add streching source term for constitutive model
            if (stabilization) then 
               call ve%get_CgradU_log(gradU,SCtmp,vf%VFold); resSC=SCtmp
            end if
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

         ! Add in relaxation forcing and reconstruct C
         if (stabilization) then 
            ! Get eigenvalues and eigenvectors
            call ve%get_eigensystem(vf%VF)
            ! Reconstruct conformation tensor
            call ve%reconstruct_conformation(vf%VF)
            ! Add in relaxtion source from semi-anlaytical integration
            call ve%get_relax_analytical(time%dt,vf%VF)
            ! Reconstruct lnC for next time step
            !> get eigenvalues and eigenvectors based on reconstructed C
            call ve%get_eigensystem_SCrec(vf%VF)
            !> Reconstruct lnC from eigenvalues and eigenvectors
            call ve%reconstruct_log_conformation(vf%VF)
            ! Take exp(eigenvalues) to use in next time-step
            ve%eigenval=exp(ve%eigenval)
         end if
         
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

            ! Add polymer stress term
            polymer_stress: block
               use tpviscoelastic_class, only: oldroydb
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
               stress=0.0_WP
               if (stabilization) then !< Build stress tensor from reconstructed C
                  select case (ve%model)
                  case (oldroydb)
                     coeff=ve%visc_p/ve%trelax
                     do k=cfg%kmino_,cfg%kmaxo_
                        do j=cfg%jmino_,cfg%jmaxo_
                           do i=cfg%imino_,cfg%imaxo_
                              if (ve%mask(i,j,k).ne.0) cycle
                              if (vf%VF(i,j,k).eq.0.0_WP) cycle
                              stress(i,j,k,1)=coeff*(ve%SCrec(i,j,k,1)-1.0_WP) !> xx tensor component
                              stress(i,j,k,2)=coeff*(ve%SCrec(i,j,k,2)-0.0_WP) !> xy tensor component
                              stress(i,j,k,3)=coeff*(ve%SCrec(i,j,k,3)-0.0_WP) !> xz tensor component
                              stress(i,j,k,4)=coeff*(ve%SCrec(i,j,k,4)-1.0_WP) !> yy tensor component
                              stress(i,j,k,5)=coeff*(ve%SCrec(i,j,k,5)-0.0_WP) !> yz tensor component
                              stress(i,j,k,6)=coeff*(ve%SCrec(i,j,k,6)-1.0_WP) !> zz tensor component
                           end do
                        end do
                     end do
                  end select 
               end if
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
            ! Update surfmesh object
            create_smesh: block
               use irl_fortran_interface
               integer :: i,j,k,nplane,np
               ! Transfer polygons to smesh
               call vf%update_surfmesh(smesh)
               ! Initalize variables to 0
               smesh%var(1,:)=0.0_WP
               smesh%var(2,:)=0.0_WP
               smesh%var(3,:)=0.0_WP
               smesh%var(4,:)=0.0_WP
               smesh%var(5,:)=0.0_WP
               smesh%var(6,:)=0.0_WP
               smesh%var(7,:)=0.0_WP
               np=0
               do k=vf%cfg%kmin_,vf%cfg%kmax_
                  do j=vf%cfg%jmin_,vf%cfg%jmax_
                     do i=vf%cfg%imin_,vf%cfg%imax_
                        do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                           if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                              np=np+1; 
                              smesh%var(1,np)=ve%SCrec(i,j,k,1)+ve%SCrec(i,j,k,4)+ve%SCrec(i,j,k,6)
                              smesh%var(2,np)=ve%SCrec(i,j,k,1)
                              smesh%var(3,np)=ve%SCrec(i,j,k,2)
                              smesh%var(4,np)=ve%SCrec(i,j,k,3)
                              smesh%var(5,np)=ve%SCrec(i,j,k,4)
                              smesh%var(6,np)=ve%SCrec(i,j,k,5)
                              smesh%var(7,np)=ve%SCrec(i,j,k,6)
                           end if
                        end do
                     end do
                  end do
               end do
            end block create_smesh
            ! Perform ensight output
            call ens_out%write_data(time%t)
         end if 

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         if (stabilization) then
            call ve%get_max_reconstructed(vf%VF)
         end if
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         
      end do
      
      ! Post-process growth rate using ODRPACK
      ! odr_fit: block
      !    use, intrinsic :: iso_fortran_env, only: output_unit
      !    use mathtools, only: twoPi
      !    use messager,  only: log
      !    use string,    only: str_long
      !    character(len=str_long) :: message
      !    integer :: i
      !    ! ODRPACK variables - explicit model based on exponential of time
      !    integer                       :: N                      !> Number of observations (number of polygons)
      !    integer , parameter           :: M=1                    !> Number of elements per explanatory variables (1 time)
      !    integer , parameter           :: NP=2                   !> Number of parameters in our model (2 for a normalized exponential in time with time shift)
      !    integer , parameter           :: NQ=1                   !> Number of response per observation (only 1, the normalized amplitude)
      !    real(WP), dimension(NP)       :: BETA=0.0_WP            !> Array of model parameter values (the growth rate and time shift)
      !    real(WP), dimension(:,:)  , allocatable :: YY           !> Value of response variable (of size LDYYxNQ)
      !    integer                       :: LDYY                   !> Leading dimension of YY (equals N since an explicit model is used)
      !    real(WP), dimension(:,:)  , allocatable :: XX           !> Value of explanatory variable (of size LDXXxM)
      !    integer                       :: LDXX                   !> Leading dimension of XX (equals N)
      !    real(WP), dimension(:,:,:), allocatable :: WE           !> Weighting of response data (of size LDWExLD2WExNQ)
      !    integer                       :: LDWE                   !> Leading dimension of WE (equals N since an explicit model is used)
      !    integer                       :: LD2WE                  !> Second dimension of WE (equals NQ)
      !    real(WP), dimension(:,:,:), allocatable :: WD           !> Weighting of explanatory data (of size LDWDxLD2WDxM)
      !    integer                       :: LDWD                   !> Leading dimension of WD (equals N)
      !    integer                       :: LD2WD                  !> Second dimension of WD (equals 1)
      !    integer , dimension(NP)       :: IFIXB=-1               !> Whether any model parameters has to be kept constant
      !    integer , parameter           :: LDIFX=1                !> Leading dimension of IFIXX (equals 1)
      !    integer , dimension(LDIFX,M)  :: IFIXX=-1               !> Whether any explanatory variable data is to be treated as "fixed"
      !    integer                       :: JOB=00030              !> 5-digit parameter flag that controls execution (this invokes analytical Jacobian with explicit model)
      !    integer                       :: NDIGIT=1               !> Number of reliable digits in our model - let ODRPACK figure it out on its own
      !    real(WP)                      :: TAUFAC=0.0_WP          !> To control size of first step (ignored here)
      !    real(WP)                      :: SSTOL=-1.0_WP          !> Relative cvg of sum of squares: this sets it to 1e-8             ********* Need to change to sth else
      !    real(WP)                      :: PARTOL=-1.0_WP         !> Relative cvg for model parameters: this sets it to 1e-11         ********* Need to change to sth else
      !    integer                       :: MAXIT=-1               !> Maximum number of iterations                                     ********* Need to change to sth else
      !    integer                       :: IPRINT=0               !> 4-digit parameter flag for controlling printing (default is -1)
      !    integer                       :: LUNERR=10              !> Logical unit for error reporting (6 by default)
      !    integer                       :: LUNRPT=10              !> Logical unit for reporting
      !    real(WP), dimension(NP)       :: STPB=0.0_WP            !> Relative step sizes for Jacobian for model parameters (here, default)
      !    integer , parameter           :: LDSTPD=1               !> Leading dimension of STPD, either 1 or N (here, 1)
      !    real(WP), dimension(LDSTPD,1) :: STPD=0.0_WP            !> Relative step sizes for Jacobian for input errors (here, default)
      !    real(WP), dimension(NP)       :: SCLB=1.0_WP            !> Scaling for the model parameters (here, not default but set to 1.0 to avoid rescaling 0 coefficients)
      !    real(WP), dimension(:,:)  , allocatable :: SCLD         !> Scaling for the input errors (here, not default but set to 1.0 to avoid rescaling 0 coefficients)
      !    integer                       :: LDSCLD                 !> Leading dimension of SCLD, either 1 or N (here, N)
      !    integer                       :: LWORK                  !> Size of WORK array
      !    real(WP), dimension(:)    , allocatable :: WORK         !> WORK array
      !    integer , parameter           :: LiWORK=20+NP+NQ*(NP+M) !> Size of IWORK array
      !    integer , dimension(LiWORK)   :: iWORK                  !> iWORK array
      !    integer                       :: INFO                   !> Why the calculations stopped
      !    ! Copy over data and sizes
      !    N=size(all_time,dim=1)
      !    LDYY=N; allocate(YY(LDYY,NQ)); YY(:,1)=all_amp/amp0
      !    LDXX=N; allocate(XX(LDXX,M )); XX(:,1)=all_time
      !    LDWE=N; LD2WE=NQ; allocate(WE(LDWE,LD2WE,NQ)); WE=1.0_WP
      !    LDWD=N; LD2WD=1 ; allocate(WD(LDWD,LD2WD,M )); WD=1.0_WP
      !    LDSCLD=N; allocate(SCLD(LDSCLD,M)); SCLD=1.0_WP
      !    LWORK=18+11*NP+NP**2+M+M**2+4*N*NQ+6*N*M+2*N*NQ*NP+2*N*NQ*M+NQ**2+5*NQ+NQ*(NP+M)+(LDWE*LD2WE)*NQ; allocate(WORK(LWORK))
      !    ! Call ODRPACK to find time shift
      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)
      !    ! Adjust weights to eliminate the early non-exponential part
      !    do i=1,size(all_time,dim=1)
      !       if (all_time(i).le.2.0_WP*BETA(2)) then
      !          WE(i,1,1)=0.0_WP
      !          WD(i,1,1)=0.0_WP
      !       end if
      !    end do
      !    ! Call ODRPACK again to find growth rate
      !    call DODRC(exponential_model,N,M,NP,NQ,BETA,YY,LDYY,XX,LDXX,WE,LDWE,LD2WE,WD,LDWD,LD2WD,IFIXB,IFIXX,LDIFX,JOB,NDIGIT,TAUFAC,&
      !    &          SSTOL,PARTOL,MAXIT,IPRINT,LUNERR,LUNRPT,STPB,STPD,LDSTPD,SCLB,SCLD,LDSCLD,WORK,LWORK,iWORK,LiWORK,INFO)
      !    ! Get back growth rate
      !    if (fs%cfg%amRoot) then
      !       write(output_unit,'(es12.5,x,es12.5,x,es12.5,x,es12.5)') lc,tau,twoPi/fs%cfg%xL*lc,BETA(1)*tau
      !       write(message    ,'("Reference time scale   = ",es12.5)') tau               ; call log(message)
      !       write(message    ,'("Cut-off length scale   = ",es12.5)') lc                ; call log(message)
      !       write(message    ,'("Normalized growth rate = ",es12.5)') BETA(1)*tau       ; call log(message)
      !       write(message    ,'("Normalized wave number = ",es12.5)') twoPi/fs%cfg%xL*lc; call log(message)
      !    end if
      ! end block odr_fit
      
      
   end subroutine simulation_run
   
   
   !> Definition of our exponential function of time model
   subroutine exponential_model(N,M,NP,NQ,LDN,LDM,LDNP,BETA,XPLUSD,IFIXB,IFIXX,LDFIX,IDEVAL,F,FJACB,FJACD,ISTOP)
      implicit none
      ! Input parameters
      integer , intent(in) :: IDEVAL,LDFIX,LDM,LDN,LDNP,M,N,NP,NQ
      integer , dimension(NP)     , intent(in) :: IFIXB
      integer , dimension(LDFIX,M), intent(in) :: IFIXX
      real(WP), dimension(NP)     , intent(in) :: BETA
      real(WP), dimension(LDN,M)  , intent(in) :: XPLUSD
      ! Output parameters
      real(WP), dimension(LDN,NQ) :: F
      real(WP), dimension(LDN,LDNP,NQ) :: FJACB
      real(WP), dimension(LDN,LDM ,NQ) :: FJACD
      integer :: ISTOP,i
      ! Check stopping condition - all values are acceptable
      ISTOP=0
      ! Compute model value
      if (mod(IDEVAL,10).ge.1) then
         do i=1,N
            F(i,1)=exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
      ! Compute model derivatives with respect to BETA
      if (mod(IDEVAL/10,10).GE.1) then
         do i=1,N
            FJACB(i,1,1)=(XPLUSD(i,1)-BETA(2))*exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
            FJACB(i,2,1)=            -BETA(1) *exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
      ! Compute model derivatives with respect to input
      if (mod(IDEVAL/100,10).GE.1) then
         do i=1,N
            FJACD(i,1,1)=BETA(1)*exp(BETA(1)*(XPLUSD(i,1)-BETA(2)))
         end do
      end if
   end subroutine exponential_model
   
   
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
      deallocate(resSC,SCtmp,gradU)
      
   end subroutine simulation_final
   
   
end module simulation
