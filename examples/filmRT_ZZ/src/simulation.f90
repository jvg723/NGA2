!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use tpns_class,        only: tpns
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(tpns),        public :: fs
   type(hypre_str)   :: ps
   type(ddadi)       :: vs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   type(surfmesh) :: smesh
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Uslip,Vslip,Wslip
   real(WP), dimension(:,:,:,:), allocatable :: bary_temp
   
   !> Problem definition
   real(WP) :: amp0,Hfilm,dh
   integer :: nh
   real(WP), dimension(:,:), allocatable :: ph

contains
   
   
   !> Function that defines a level set function for a perturbed film problem
   function levelset_film(xyz,t) result(G)
      use mathtools, only: twoPi
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G1,G2,G,dist
      real(WP), dimension(3) :: c,pos,myh
      integer :: n
      ! G1=+xyz(3)+0.5_WP*Hfilm!-amp0*cos(twoPi*xyz(1))
      ! G2=-xyz(3)+0.5_WP*Hfilm!+amp0*cos(twoPi*xyz(1))
      ! G=min(G1,G2)
      G=0.5_WP*Hfilm-abs(xyz(3))
      ! Add the holes
      ! do n=1,nh
      !    ! Store 2D positions
      !    pos=[xyz(1) ,0.0_WP,xyz(3) ]
      !    myh=[ph(1,n),0.0_WP,ph(2,n)]
      !    ! Account for periodicity
      !    if (myh(1)-pos(1).gt.+0.5_WP*cfg%xL) myh(1)=myh(1)-cfg%xL
      !    if (myh(1)-pos(1).lt.-0.5_WP*cfg%xL) myh(1)=myh(1)+cfg%xL
      !    if (myh(3)-pos(3).gt.+0.5_WP*cfg%zL) myh(3)=myh(3)-cfg%zL
      !    if (myh(3)-pos(3).lt.-0.5_WP*cfg%zL) myh(3)=myh(3)+cfg%zL
      !    ! Check distance to hole center
      !    dist=norm2(myh-pos)
      !    ! Check if we are right on top of hole, if so move it a bit
      !    if (dist.eq.0.0_WP) then
      !       myh=myh+[tiny(1.0_WP),0.0_WP,tiny(1.0_WP)]
      !       dist=norm2(myh-pos)
      !    end if
      !    ! Check if hole is close enough
      !    if (dist.lt.0.5_WP*dh) then
      !       ! Get closest point on torus
      !       c=myh+0.5_WP*dh*(pos-myh)/dist
      !       G=0.5_WP-norm2(xyz-c)
      !    end if
      ! end do

   end function levelset_film
   
   
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
         allocate(Uslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Uslip=0.0_WP
         allocate(Vslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vslip=0.0_WP
         allocate(Wslip(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wslip=0.0_WP
         allocate(bary_temp(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); bary_temp=0.0_WP
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
       ! Prepare random holes for film perforation
      initialize_holes: block
         use random,   only: random_uniform
         use parallel, only: MPI_REAL_WP
         use mpi_f08,  only: MPI_BCAST
         integer  :: n,nn,ierr
         real(WP) :: xh,zh,xxh,zzh
         logical  :: is_overlap
         real(WP), parameter :: safety_margin=1.5_WP
         ! Read in the hole size
         call param_read('Size of holes',dh)
         ! Read in the number of holes
         call param_read('Number of holes',nh)
         ! Allocate hole position array
         allocate(ph(2,nh)); ph=0.0_WP
         ! If only one hole, leave it in the middle
         if (nh.gt.1) then
            ! Root assigns random non-overlapping hole positions
            if (cfg%amRoot) then
               n=0
               do while (n.lt.nh)
                  ! Draw a random position on the film
                  xh=random_uniform(lo=cfg%x(cfg%imin),hi=cfg%x(cfg%imax+1))
                  zh=random_uniform(lo=cfg%z(cfg%kmin),hi=cfg%z(cfg%kmax+1))
                  ! Compare to all previous holes
                  is_overlap=.false.
                  do nn=1,n
                     ! Get the position of the other hole
                     xxh=ph(1,nn); zzh=ph(2,nn)
                     ! Account for periodicity
                     if (xxh-xh.gt.+0.5_WP*cfg%xL) xxh=xxh-cfg%xL
                     if (xxh-xh.lt.-0.5_WP*cfg%xL) xxh=xxh+cfg%xL
                     if (zzh-zh.gt.+0.5_WP*cfg%zL) zzh=zzh-cfg%zL
                     if (zzh-zh.lt.-0.5_WP*cfg%zL) zzh=zzh+cfg%zL
                     ! Check for overlap
                     if (norm2([xxh-xh,zzh-zh]).lt.safety_margin*dh) is_overlap=.true.
                  end do
                  ! If no overlap was found, add the hole to the list
                  if (.not.is_overlap) then
                     n=n+1
                     ph(1,n)=xh
                     ph(2,n)=zh
                  end if
               end do
            end if
            ! Broadcoast the hole positions
            call MPI_BCAST(ph,2*nh,MPI_REAL_WP,0,cfg%comm,ierr)
         end if
      end block initialize_holes

      
      ! Initialize our VOF solver and field （This creates a single plane interface in each cell based on level set)
      create_and_initialize_vof: block
         use mathtools, only: twoPi
         use mms_geom, only: cube_refine_vol
         use vfs_class, only: r2p,lvira,VFhi,VFlo,flux
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         integer :: nx
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=r2p,transport_method=flux,name='VOF')
         ! vf=vfs(cfg=cfg,reconstruction_method=r2p,name='VOF')
         ! Initialize to a thin film
         call param_read('Film thickness',Hfilm)
         call param_read('nx',nx)
         amp0=1.0e-3_WP
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
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_film,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)

                  if((i.eq.nx/2.and.j.eq.nx/2).or.(i.eq.nx/2.and.j.eq.nx/2+1).or.(i.eq.nx/2+1.and.j.eq.nx/2).or.(i.eq.nx/2+1.and.j.eq.nx/2+1)) then
                     vf%VF(i,j,k)= 0.0_WP
                  end if
                  
                  ! if((i.eq.nx/2-1.and.j.eq.nx/2).or.(i.eq.nx/2-1.and.j.eq.nx/2+1).or.(i.eq.nx/2+2.and.j.eq.nx/2).or.(i.eq.nx/2+2.and.j.eq.nx/2+1)) then
                  !    vf%VF(i,j,k)= vf%VF(i,j,k)/2.0_WP 
                  ! end if

                  ! if((i.eq.nx/2.and.j.eq.nx/2-1).or.(i.eq.nx/2+1.and.j.eq.nx/2-1).or.(i.eq.nx/2.and.j.eq.nx/2+2).or.(i.eq.nx/2+1.and.j.eq.nx/2+2)) then
                  !    vf%VF(i,j,k)= vf%VF(i,j,k)/2.0_WP 
                  ! end if

                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if

                  ! if (i.eq.nx/2.and.j.eq.nx/2-1.and.(k.eq.1.or.k.eq.2)) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i+1),vf%cfg%y(j+1)]-[0.5274_WP,0.2732_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2+1.and.j.eq.nx/2-1) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i),vf%cfg%y(j+1)]+[0.5274_WP,-0.2732_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2-1.and.j.eq.nx/2) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i+1),vf%cfg%y(j+1)]-[0.2732_WP,0.5274_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2+2.and.j.eq.nx/2) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i),vf%cfg%y(j+1)]+[0.2732_WP,-0.5274_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2-1.and.j.eq.nx/2+1) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i+1),vf%cfg%y(j)]+[-0.2732_WP,0.5274_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2+2.and.j.eq.nx/2+1) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i),vf%cfg%y(j)]+[0.2732_WP,0.5274_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2.and.j.eq.nx/2+2) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i+1),vf%cfg%y(j)]+[-0.5274_WP,0.2732_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! else if (i.eq.nx/2+1.and.j.eq.nx/2+2) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  !    vf%Gbary(1:2,i,j,k)= [vf%cfg%x(i),vf%cfg%y(j)]+[0.5274_WP,0.2732_WP]*vf%cfg%dx(i)
                  !    vf%Lbary(1:2,i,j,k)= [vf%cfg%xm(i),vf%cfg%ym(j)]-((1-vf%VF(i,j,k))*vf%Gbary(1:2,i,j,k))/vf%VF(i,j,k)
                  ! end if

                  ! if(  (i.eq.nx/2.and.j.eq.nx/2-1).or.(i.eq.nx/2+1.and.j.eq.nx/2-1).or.(i.eq.nx/2-1.and.j.eq.nx/2).or.(i.eq.nx/2+2.and.j.eq.nx/2)&
                  ! &.or.(i.eq.nx/2-1.and.j.eq.nx/2+1).or.(i.eq.nx/2+2.and.j.eq.nx/2+1).or.(i.eq.nx/2.and.j.eq.nx/2+2).or.(i.eq.nx/2+1.and.j.eq.nx/2+2)) then
                  !    vf%VF(i,j,k)= (1-(twoPi-4.0_WP)/8.0)*Hfilm*0.5_WP/vf%cfg%dx(i)
                  ! end if


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
      

      ! ! Initialize our VOF solver and field （This provides a set of prescribed 2 planes interface for my inital hole)
      ! create_and_initialize_vof: block
      !    use mms_geom, only: cube_refine_vol
      !    use vfs_class, only: r2p,lvira,VFhi,VFlo,flux
      !    use irl_fortran_interface
      !    integer :: i,j,k,n,si,sj,sk
      !    real(WP), dimension(3,8) :: cube_vertex
      !    real(WP), dimension(3) :: v_cent,a_cent,norm
      !    real(WP) :: vol,area,dist
      !    integer, parameter :: amr_ref_lvl=4
      !    integer :: nx
      !    ! Create a VOF solver
      !    call vf%initialize(cfg=cfg,reconstruction_method=r2p,transport_method=flux,name='VOF')
      !    ! vf=vfs(cfg=cfg,reconstruction_method=r2p,name='VOF')
      !    ! Initialize to a thin film
      !    call param_read('Film thickness',Hfilm)
      !    call param_read('nx',nx)

      !    do k=vf%cfg%kmin_,vf%cfg%kmax_
      !       do j=vf%cfg%jmin_,vf%cfg%jmax_
      !          do i=vf%cfg%imin_,vf%cfg%imax_
      !             ! norm=normalize([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                  
      !             if(k.eq.2) then
      !                call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),2)
      !                norm = [0.0_WP,0.0_WP,-1.0_WP]
      !                dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm*1.0_WP)
      !                call setPlane(vf%liquid_gas_interface(i,j,k),0,norm,dist)  
      !                norm = [0.0_WP,0.0_WP,1.0_WP]
      !                dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm*1.0_WP)
      !                call setPlane(vf%liquid_gas_interface(i,j,k),1,norm,dist)  
      !             end if

                  
      !             ! if(k.eq.1) then
      !             !    call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
      !             !    norm = [0.0_WP,0.0_WP,-1.0_WP]
      !             !    dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm)
      !             !    call setPlane(vf%liquid_gas_interface(i,j,k),0,norm,dist)  
      !             !    ! norm = [0.0_WP,0.0_WP,1.0_WP]
      !             !    ! dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm)
      !             !    ! call setPlane(vf%liquid_gas_interface(i,j,k),1,norm,dist)  
      !             ! end if


      !             ! if(k.eq.3) then
      !             !    call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
      !             !    ! norm = [0.0_WP,0.0_WP,-1.0_WP]
      !             !    ! dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm)
      !             !    ! call setPlane(vf%liquid_gas_interface(i,j,k),0,norm,dist)  
      !             !    norm = [0.0_WP,0.0_WP,1.0_WP]
      !             !    dist=dot_product(norm,[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]+norm*Hfilm)
      !             !    call setPlane(vf%liquid_gas_interface(i,j,k),0,norm,dist)  
      !             ! end if

      !          end do
      !       end do
      !    end do
      !    call vf%sync_interface()
      !    ! Reset moments from r2p
      !    call vf%reset_volume_moments()
      !    ! vf%VF = 1-vf%VF
      !    ! bary_temp = vf%Lbary
      !    ! vf%Lbary  =vf%Gbary
      !    ! vf%Gbary =bary_temp
      !    ! Update the band
      !    call vf%update_band()
      !    ! Create discontinuous polygon mesh from IRL interface
      !    call vf%polygonalize_interface()
      !    ! Perform interface sensing
      !    if (vf%two_planes) call vf%sense_interface()
      !    ! Calculate distance from polygons
      !    call vf%distance_from_polygon()
      !    ! Calculate subcell phasic volumes
      !    call vf%subcell_vol()
      !    ! Calculate curvature
      !    call vf%get_curvature()
      ! end block create_and_initialize_vof

      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         ! use ils_class, only: pcg_pfmg
         use hypre_str_class, only: pcg_pfmg2
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
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! call param_read('Pressure iteration',fs%psolv%maxit)
         ! call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! call param_read('Implicit iteration',fs%implicit%maxit)
         ! call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         ! call fs%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=5,name='plic')
         smesh%varname(1)='nplane'
         smesh%varname(2)='curv'
         smesh%varname(3)='edge_sensor'
         smesh%varname(4)='thin_sensor'
         smesh%varname(5)='thickness'
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
      end block create_smesh

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='filmRT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('SLIPvelocity',Uslip,Vslip,Wslip)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_vector('ST',fs%Pjx,fs%Pjy,fs%Pjz)
         call ens_out%add_surface('vofplic',smesh)
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
      implicit none
      integer :: i,j,k
      
      ! Perform time integration
      do while (time%n.lt.100)
         
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
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         Uslip=fs%U
         Vslip=fs%V
         Wslip=fs%W

         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  if (vf%edge_sensor(i,j,k).gt.0.03_WP) then
                     ! print *, i,j,k,"The normal is: ", vf%edge_normal(:,i,j,k)
                     ! Uslip(i,j,k) = Uslip(i,j,k)+vf%edge_normal(1,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     ! Vslip(i,j,k) = Vslip(i,j,k)+vf%edge_normal(2,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     ! Wslip(i,j,k) = Wslip(i,j,k)+vf%edge_normal(3,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     ! print *, i,j,k,"U,V,W SLIP VELOCITY ", Uslip(i,j,k),Vslip(i,j,k),Wslip(i,j,k)
                     Uslip(i  ,j,k) = Uslip(i  ,j,k) + 0.5_WP*vf%edge_normal(1,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     Uslip(i+1,j,k) = Uslip(i+1,j,k) + 0.5_WP*vf%edge_normal(1,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))

                     Vslip(i,j  ,k) = Vslip(i,j  ,k) + 0.5_WP*vf%edge_normal(2,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     Vslip(i,j+1,k) = Vslip(i,j+1,k) + 0.5_WP*vf%edge_normal(2,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))

                     Wslip(i,j,k  ) = Wslip(i,j,k  ) + 0.5_WP*vf%edge_normal(3,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                     Wslip(i,j,k+1) = Wslip(i,j,k+1) + 0.5_WP*vf%edge_normal(3,i,j,k)*sqrt(2.0_WP*fs%sigma/(fs%rho_l*Hfilm))
                  end if


               end do 
            end do 
         end do
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
      
   end subroutine simulation_final
   
   
end module simulation
