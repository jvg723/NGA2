!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use film_class,        only: film
   use config_class,      only: config
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   
   implicit none
   private
   
   !> Film simulation
   type(film) :: dns
   
   !> R2P post-processing
   type(config)   :: cfg
   type(vfs)      :: vf
   type(tpns)        :: fs    !< Two-phase flow solver
   type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
   type(ddadi)       :: vs    !< DDADI solver for velocity
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: pproc_evt
   type(timetracker) :: time

   type(monitor) :: mfile    !< General simulation monitoring
   type(monitor) :: cflfile  !< CFL monitoring

   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Uslip,Vslip,Wslip
   ! real(WP), dimension(:,:,:), allocatable :: Um,Vm,Wm
   ! real(WP), dimension(:,:,:), allocatable :: Ul,Vl,Wl
   ! real(WP), dimension(:,:,:), allocatable :: Ug,Vg,Wg
   real(WP) :: rho_ratio,visc_ratio,Oh,dh

   public :: simulation_init,simulation_run,simulation_final

contains
   
   !> Initialization of our simulation
   subroutine simulation_init
      
      ! Initialize film simulation
      call dns%init()

      ! Create a coarse grid config
      coarse_grid: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,coarse_factor
         real(WP) :: Lx,Ly,Lz
         ! Coarsening factor
         call param_read('Coarsening factor',coarse_factor)
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); nx=nx/coarse_factor; allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny)
         Ly=3.0_WP*Ly; ny=3*ny !< Force the coarse mesh to be 3D
         ny=ny/coarse_factor; allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); nz=nz/coarse_factor; allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General coarse serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.true.,name='coarseFilm')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block coarse_grid
      
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ui=0.0_WP
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vi=0.0_WP
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wi=0.0_WP
         allocate(Uslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Uslip=0.0_WP
         allocate(Vslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vslip=0.0_WP
         allocate(Wslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wslip=0.0_WP
      end block allocate_work_arrays

      initialize_timetracker: block
         use param, only: param_read
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Create a coarse VF solver
      coarse_vof: block
         use vfs_class, only: r2p,remap,flux
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=r2p,transport_method=flux,name='coarseVOF')
      end block coarse_vof

      coarse_fs: block
         use param, only: param_read
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant density and viscosity to each phase
         fs%rho_l =1.0_WP; call param_read('Density ratio'  , rho_ratio); fs%rho_g =fs%rho_l / rho_ratio
         fs%visc_l=1.0_WP; call param_read('Viscosity ratio',visc_ratio); fs%visc_g=fs%visc_l/visc_ratio
         ! ! Read in Ohnesorge number and assign surface tension coefficient
         call param_read('Oh',Oh); fs%sigma=Oh**(-2)
         ! ! Prepare and configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         !ps%maxlevel=16
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! ! ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block coarse_fs

      ! Create a coarse surfmesh object for interface polygon output
      coarse_smesh: block
         smesh=surfmesh(nvar=5,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='edge_sensor'
         smesh%varname(4)='thin_sensor'
         smesh%varname(5)='thickness'
      end block coarse_smesh
      
      ! Also create ensight output
      coarse_ensight: block
         ! Allocate coarse velocity storage
         ! allocate(Um(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Vm(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Wm(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Ul(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Vl(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Wl(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Ug(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Vg(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(Wg(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='coarse')
         ! Add variables to output
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('SLIPvelocity',Uslip,Vslip,Wslip)
         call ens_out%add_scalar('curvature',vf%curv)
         ! call ens_out%add_vector('ST',fs%Pjx,fs%Pjy,fs%Pjz)
         ! call ens_out%add_vector('Um',Um,Vm,Wm)
         ! call ens_out%add_vector('Ul',Ul,Vl,Wl)
         ! call ens_out%add_vector('Ug',Ug,Vg,Wg)
         call ens_out%add_surface('r2p',smesh)
      end block coarse_ensight
      
      ! Finally, create event for postprocessing
      pproc_event: block
         use param, only: param_read
         pproc_evt=event(time=time,name='Postproc')
         call param_read('Postproc period',pproc_evt%tper)
      end block pproc_event

      ! Perform postprocessing
      ! if (pproc_evt%occurs()) call pproc()
      call pproc()
      

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
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k,iunit
      real(WP) :: Hfilm 
      Hfilm = 1.0_WP
      iunit =1
      open(unit=iunit,file='VFint.csv',status='unknown')
      write(iunit,'(a12,3x,a12,3x,a12)') 'Time','VFint'

      do while (time%n.lt.300)
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
         call fs%add_slipvel(vf,Uslip,Vslip,Wslip,Hfilm)
         ! Uslip = Uslip +  fs%U
         ! Vslip = Vslip +  fs%V
         ! Wslip = Wslip +  fs%W
         ! VOF solver step
         ! call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         call vf%advance(dt=time%dt,U=Uslip,V=Vslip,W=Wslip)

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)
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
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
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
            call fs%add_surface_tension_jump_thin(dt=time%dt,div=fs%div,vf=vf)
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
         ! Update surfmesh object
         update_smesh: block
            use irl_fortran_interface
            integer :: nplane,np,i,j,k
            ! Transfer polygons to smesh
            call vf%update_surfmesh(smesh)
            ! Also populate nplane variable
            smesh%var(1,:)=0.0_WP
            np=0
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; 
                           smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                           smesh%var(2,np)=vf%curv2p(nplane,i,j,k)
                           smesh%var(3,np)=vf%edge_sensor(i,j,k)
                           smesh%var(4,np)=vf%thin_sensor(i,j,k)
                           smesh%var(5,np)=vf%thickness  (i,j,k)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Output to ensight
         ! if (ens_evt%occurs()) 
         call ens_out%write_data(time%t)
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()
         write(iunit,'(1x,f24.16,1x,f24.16)') time%t,vf%VFint
      end do

      ! Film drives overall time integration
      ! do while (.not.dns%time%done())

      !    ! Advance film simulation
      !    call dns%step()
         
      !    ! Perform postprocessing
      !    if (pproc_evt%occurs()) call pproc()
         
      ! end do
      
   end subroutine simulation_run


   !> Postprocess phasic moments from DNS
   !> Feed to coarse mesh for R2P reconstruction
   subroutine pproc()
      implicit none
      
      ! Integrate volume moments from fine to coarse mesh
      compute_volume_moments: block
         use vfs_class, only: VFlo,VFhi
         integer :: fi,fj,fk
         integer :: ci,cj,ck
         ! Reset coarse volume moments
         vf%VF=0.0_WP; vf%Lbary=0.0_WP; vf%Gbary=0.0_WP
         ! Loop over fine mesh, find corresponding coarse mesh cell, integrate moments
         do fk=dns%cfg%kmin_,dns%cfg%kmax_
            ck=cfg%kmin_; do while (dns%cfg%zm(fk).ge.cfg%z(ck+1)); ck=ck+1; end do
            do fj=dns%cfg%jmin_,dns%cfg%jmax_
               cj=cfg%jmin_; do while (dns%cfg%ym(fj).ge.cfg%y(cj+1)); cj=cj+1; end do
               do fi=dns%cfg%imin_,dns%cfg%imax_
                  ci=cfg%imin_; do while (dns%cfg%xm(fi).ge.cfg%x(ci+1)); ci=ci+1; end do
                  vf%VF(ci,cj,ck)=vf%VF(ci,cj,ck)+dns%vf%VF(fi,fj,fk)*dns%cfg%vol(fi,fj,fk)
                  vf%Lbary(:,ci,cj,ck)=vf%Lbary(:,ci,cj,ck)+(       dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%vf%Lbary(:,fi,fj,fk)
                  vf%Gbary(:,ci,cj,ck)=vf%Gbary(:,ci,cj,ck)+(1.0_WP-dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%vf%Gbary(:,fi,fj,fk)
               end do
            end do
         end do
         do ck=cfg%kmin_,cfg%kmax_
            do cj=cfg%jmin_,cfg%jmax_
               do ci=cfg%imin_,cfg%imax_
                  vf%VF(ci,cj,ck)=vf%VF(ci,cj,ck)/cfg%vol(ci,cj,ck)
                  if (vf%VF(ci,cj,ck).ge.VFlo) vf%Lbary(:,ci,cj,ck)=vf%Lbary(:,ci,cj,ck)/((       vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
                  if (vf%VF(ci,cj,ck).le.VFhi) vf%Gbary(:,ci,cj,ck)=vf%Gbary(:,ci,cj,ck)/((1.0_WP-vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
               end do
            end do
         end do
         call cfg%sync(vf%VF)
         call vf%sync_and_clean_barycenters()
      end block compute_volume_moments
      
      ! Integrate surface moments from fine to coarse mesh
      compute_surface_moments: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE,MPI_INTEGER
         use parallel, only: MPI_REAL_WP
         use irl_fortran_interface
         use vfs_class, only : VFlo,VFhi
         integer :: fi,fj,fk,n,t,ierr
         integer :: ci,cj,ck
         type(DivPoly_type) :: divided_polygon
         type(Tri_type) :: triangle
         real(IRL_double), dimension(1:4) :: plane_data
         real(IRL_double), dimension(:,:,:), allocatable :: tri,tmp
         integer, dimension(:), allocatable :: ntri_proc
         integer :: ntri,list_size,localizer_id,ii,count
         real(WP), dimension(3) :: bary
         integer, dimension(3) :: ind
         type(TagAccListVM_VMAN_type) :: accumulated_moments_from_tri
         type(ListVM_VMAN_type) :: moments_list_from_tri
         type(VMAN_type) :: volume_moments_and_normal
         ! Allocate IRL data
         call new(divided_polygon)
         call new(triangle)
         call new(accumulated_moments_from_tri)
         call new(moments_list_from_tri)
         call new(volume_moments_and_normal)
         ! Count the number of surface triangles in DNS and communicate
         allocate(ntri_proc(1:dns%cfg%nproc)); ntri_proc=0
         do fk=dns%cfg%kmin_,dns%cfg%kmax_
            do fj=dns%cfg%jmin_,dns%cfg%jmax_
               do fi=dns%cfg%imin_,dns%cfg%imax_
                  ! Skip if no interface
                  if (dns%vf%VF(fi,fj,fk).lt.VFlo.or.dns%vf%VF(fi,fj,fk).gt.VFhi) cycle
                  ! Construct triangulation of each interface plane
                  do n=1,getNumberOfPlanes(dns%vf%liquid_gas_interface(fi,fj,fk))
                     ! Skip planes outside of the cell
                     if (getNumberOfVertices(dns%vf%interface_polygon(n,fi,fj,fk)).eq.0) cycle
                     ! Get DividedPolygon from the plane
                     call constructFromPolygon(divided_polygon,dns%vf%interface_polygon(n,fi,fj,fk))
                     ! Increment triangle counter
                     ntri_proc(dns%cfg%rank+1)=ntri_proc(dns%cfg%rank+1)+getNumberOfSimplicesInDecomposition(divided_polygon)
                  end do
               end do
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,ntri_proc,dns%cfg%nproc,MPI_INTEGER,MPI_SUM,dns%cfg%comm,ierr)
         ntri=sum(ntri_proc)
         ! Allocate global storage of our triangles
         allocate(tri(1:3,1:3,ntri)); tri=0.0_WP
         ! Create global triangulation of the DNS interface and communicate
         count=sum(ntri_proc(1:dns%cfg%rank))
         do fk=dns%cfg%kmin_,dns%cfg%kmax_
            do fj=dns%cfg%jmin_,dns%cfg%jmax_
               do fi=dns%cfg%imin_,dns%cfg%imax_
                  ! Skip if no interface
                  if (dns%vf%VF(fi,fj,fk).lt.VFlo.or.dns%vf%VF(fi,fj,fk).gt.VFhi) cycle
                  ! Construct triangulation of each interface plane
                  do n=1,getNumberOfPlanes(dns%vf%liquid_gas_interface(fi,fj,fk))
                     ! Skip planes outside of the cell
                     if (getNumberOfVertices(dns%vf%interface_polygon(n,fi,fj,fk)).eq.0) cycle
                     ! Get DividedPolygon from the plane
                     call constructFromPolygon(divided_polygon,dns%vf%interface_polygon(n,fi,fj,fk))
                     ! Check if point ordering correct, flip if not
                     plane_data=getPlane(dns%vf%liquid_gas_interface(fi,fj,fk),n-1)
                     if (abs(1.0_WP-dot_product(calculateNormal(divided_polygon),plane_data(1:3))).gt.1.0_WP) call reversePtOrdering(divided_polygon)
                     ! Loop over triangles from DividedPolygon
                     do t=1,getNumberOfSimplicesInDecomposition(divided_polygon)
                        ! Get the triangle vertices
                        call getSimplexFromDecomposition(divided_polygon,t-1,triangle)
                        ! Store at correct global location
                        count=count+1
                        tri(1:3,1:3,count)=getVertices(triangle)
                     end do
                  end do
               end do
            end do
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,tri,9*ntri,MPI_REAL_WP,MPI_SUM,dns%cfg%comm,ierr)
         ! Finally, deal with periodicity - brute force copy
         allocate(tmp(1:3,1:3,9*ntri)); tmp=0.0_WP
         do n=1,ntri
            ! Original triangle and its 8 copies
            tmp(:,:,0*ntri+n)=tri(:,:,n) ! Original
            tmp(:,:,1*ntri+n)=tri(:,:,n) ! x-
            tmp(:,:,2*ntri+n)=tri(:,:,n) ! x+
            tmp(:,:,3*ntri+n)=tri(:,:,n) ! z-
            tmp(:,:,4*ntri+n)=tri(:,:,n) ! z+
            tmp(:,:,5*ntri+n)=tri(:,:,n) ! x-/z-
            tmp(:,:,6*ntri+n)=tri(:,:,n) ! x+/z-
            tmp(:,:,7*ntri+n)=tri(:,:,n) ! x-/z+
            tmp(:,:,8*ntri+n)=tri(:,:,n) ! x+/z+
            ! Apply periodicity in x- (1, 5, 7)
            tmp(1,:,1*ntri+n)=tmp(1,:,1*ntri+n)-dns%cfg%xL
            tmp(1,:,5*ntri+n)=tmp(1,:,5*ntri+n)-dns%cfg%xL
            tmp(1,:,7*ntri+n)=tmp(1,:,7*ntri+n)-dns%cfg%xL
            ! Apply periodicity in x+ (2, 6, 8)
            tmp(1,:,2*ntri+n)=tmp(1,:,2*ntri+n)+dns%cfg%xL
            tmp(1,:,6*ntri+n)=tmp(1,:,6*ntri+n)+dns%cfg%xL
            tmp(1,:,8*ntri+n)=tmp(1,:,8*ntri+n)+dns%cfg%xL
            ! Apply periodicity in z- (3, 5, 6)
            tmp(3,:,3*ntri+n)=tmp(3,:,3*ntri+n)-dns%cfg%zL
            tmp(3,:,5*ntri+n)=tmp(3,:,5*ntri+n)-dns%cfg%zL
            tmp(3,:,6*ntri+n)=tmp(3,:,6*ntri+n)-dns%cfg%zL
            ! Apply periodicity in z+ (4, 7, 8)
            tmp(3,:,4*ntri+n)=tmp(3,:,4*ntri+n)+dns%cfg%zL
            tmp(3,:,7*ntri+n)=tmp(3,:,7*ntri+n)+dns%cfg%zL
            tmp(3,:,8*ntri+n)=tmp(3,:,8*ntri+n)+dns%cfg%zL
         end do
         ntri=9*ntri; call move_alloc(tmp,tri)
         ! Finally, populate the coarse surface data arrays
         do ck=cfg%kmino_,cfg%kmaxo_
            do cj=cfg%jmino_,cfg%jmaxo_
               do ci=cfg%imino_,cfg%imaxo_
                  ! Reset storage and surface area
                  call clear(vf%triangle_moments_storage(ci,cj,ck))
                  vf%SD(ci,cj,ck)=0.0_WP
                  ! Traverse all DNS triangles and store the ones local to our cell
                  do n=1,ntri
                     ! Get barycenter of triangle
                     bary=(tri(:,1,n)+tri(:,2,n)+tri(:,3,n))/3.0_WP
                     ! Test if it's in our cell
                     if (bary(1).gt.cfg%x(ci).and.bary(1).lt.cfg%x(ci+1).and.&
                     &   bary(2).gt.cfg%y(cj).and.bary(2).lt.cfg%y(cj+1).and.&
                     &   bary(3).gt.cfg%z(ck).and.bary(3).lt.cfg%z(ck+1)) then
                        ! Add the triangle moments
                        call construct(triangle,tri(:,:,n))
                        call calculateAndSetPlaneOfExistence(triangle)
                        ! Cut it by the mesh
                        call getMoments(triangle,vf%localizer_link(ci,cj,ck),accumulated_moments_from_tri)
                        ! Loop through each cell and append to triangle_moments_storage in each cell
                        list_size=getSize(accumulated_moments_from_tri)
                        do ii=1,list_size
                           localizer_id=getTagForIndex(accumulated_moments_from_tri,ii-1)
                           ind=cfg%get_ijk_from_lexico(localizer_id)
                           call getListAtIndex(accumulated_moments_from_tri,ii-1,moments_list_from_tri)
                           call append(vf%triangle_moments_storage(ind(1),ind(2),ind(3)),moments_list_from_tri)
                        end do
                     end if
                  end do
                  ! Finally compute surface area
                  do ii=0,getSize(vf%triangle_moments_storage(ci,cj,ck))-1
                     call getMoments(vf%triangle_moments_storage(ci,cj,ck),ii,volume_moments_and_normal)
                     vf%SD(ci,cj,ck)=vf%SD(ci,cj,ck)+getVolume(volume_moments_and_normal)
                  end do
                  vf%SD(ci,cj,ck)=vf%SD(ci,cj,ck)/cfg%vol(ci,cj,ck)
               end do
            end do
         end do
      end block compute_surface_moments
      
      ! Perform R2P reconstruction on coarse mesh
      !----> Needs volume moments: VF, Lbary, and Gbary
      !----> Needs surface moments: triangle_moments_storage and SD
      update_vf: block
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Perform interface sensing
         call vf%sense_interface()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         !call vf%reset_volume_moments()
      end block update_vf

      ! Update smesh data
      update_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
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
                     end if
                  end do
               end do
            end do
         end do
      end block update_smesh
      
      ! ! Integrate velocities from fine to coarse mesh
      ! compute_velocities: block
      !    use vfs_class, only: VFlo,VFhi
      !    integer :: fi,fj,fk
      !    integer :: ci,cj,ck
      !    ! Reset velocities
      !    Um=0.0_WP; Vm=0.0_WP; Wm=0.0_WP
      !    Ul=0.0_WP; Vl=0.0_WP; Wl=0.0_WP
      !    Ug=0.0_WP; Vg=0.0_WP; Wg=0.0_WP
      !    ! Loop over fine mesh, find corresponding coarse mesh cell, integrate velocities
      !    do fk=dns%cfg%kmin_,dns%cfg%kmax_
      !       ck=cfg%kmin_; do while (dns%cfg%zm(fk).ge.cfg%z(ck+1)); ck=ck+1; end do
      !       do fj=dns%cfg%jmin_,dns%cfg%jmax_
      !          cj=cfg%jmin_; do while (dns%cfg%ym(fj).ge.cfg%y(cj+1)); cj=cj+1; end do
      !          do fi=dns%cfg%imin_,dns%cfg%imax_
      !             ci=cfg%imin_; do while (dns%cfg%xm(fi).ge.cfg%x(ci+1)); ci=ci+1; end do
      !             ! Mixture velocity
      !             Um(ci,cj,ck)=Um(ci,cj,ck)+                             dns%cfg%vol(fi,fj,fk)*dns%Ui(fi,fj,fk)
      !             Vm(ci,cj,ck)=Vm(ci,cj,ck)+                             dns%cfg%vol(fi,fj,fk)*dns%Vi(fi,fj,fk)
      !             Wm(ci,cj,ck)=Wm(ci,cj,ck)+                             dns%cfg%vol(fi,fj,fk)*dns%Wi(fi,fj,fk)
      !             ! Liquid velocity
      !             Ul(ci,cj,ck)=Ul(ci,cj,ck)+(       dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Ui(fi,fj,fk)
      !             Vl(ci,cj,ck)=Vl(ci,cj,ck)+(       dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Vi(fi,fj,fk)
      !             Wl(ci,cj,ck)=Wl(ci,cj,ck)+(       dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Wi(fi,fj,fk)
      !             ! Gas velocity
      !             Ug(ci,cj,ck)=Ug(ci,cj,ck)+(1.0_WP-dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Ui(fi,fj,fk)
      !             Vg(ci,cj,ck)=Vg(ci,cj,ck)+(1.0_WP-dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Vi(fi,fj,fk)
      !             Wg(ci,cj,ck)=Wg(ci,cj,ck)+(1.0_WP-dns%vf%VF(fi,fj,fk))*dns%cfg%vol(fi,fj,fk)*dns%Wi(fi,fj,fk)
      !          end do
      !       end do
      !    end do
      !    do ck=cfg%kmin_,cfg%kmax_
      !       do cj=cfg%jmin_,cfg%jmax_
      !          do ci=cfg%imin_,cfg%imax_
      !             ! Mixture velocity
      !             Um(ci,cj,ck)=Um(ci,cj,ck)/cfg%vol(ci,cj,ck)
      !             Vm(ci,cj,ck)=Vm(ci,cj,ck)/cfg%vol(ci,cj,ck)
      !             Wm(ci,cj,ck)=Wm(ci,cj,ck)/cfg%vol(ci,cj,ck)
      !             ! Liquid velocity
      !             if (vf%VF(ci,cj,ck).ge.VFlo) then
      !                Ul(ci,cj,ck)=Ul(ci,cj,ck)/((       vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !                Vl(ci,cj,ck)=Vl(ci,cj,ck)/((       vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !                Wl(ci,cj,ck)=Wl(ci,cj,ck)/((       vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !             end if
      !             ! Gas velocity
      !             if (vf%VF(ci,cj,ck).le.VFhi) then
      !                Ug(ci,cj,ck)=Ug(ci,cj,ck)/((1.0_WP-vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !                Vg(ci,cj,ck)=Vg(ci,cj,ck)/((1.0_WP-vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !                Wg(ci,cj,ck)=Wg(ci,cj,ck)/((1.0_WP-vf%VF(ci,cj,ck))*cfg%vol(ci,cj,ck))
      !             end if
      !          end do
      !       end do
      !    end do
      !    call cfg%sync(Um); call cfg%sync(Vm); call cfg%sync(Wm)
      !    call cfg%sync(Ul); call cfg%sync(Vl); call cfg%sync(Wl)
      !    call cfg%sync(Ug); call cfg%sync(Vg); call cfg%sync(Wg)
      ! end block compute_velocities

      ! Perform ensight output
      call ens_out%write_data(time%t)
      
   end subroutine pproc
   

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize film simulation
      call dns%final()
      
   end subroutine simulation_final
   
   
end module simulation
