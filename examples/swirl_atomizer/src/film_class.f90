!> Film model class:
!> Provides support for film retraction modeling

module film_class
   use precision,      only: WP
   use config_class,   only: config
   use ccl_class,      only: ccl
   use vfs_class,      only: vfs
   use tpns_class,     only: tpns
   use irl_fortran_interface
   implicit none
   private

   ! Expose type/methods
   public :: film

   !> Transfer model object
   type :: film

      ! This is our config
      class(config), pointer :: cfg

      type(ccl)           :: cc
      type(vfs) , pointer :: vf
      type(tpns), pointer :: fs

      !> Transfer model parameters
      real(WP) :: filmthickness_over_dx  =5.0e-1_WP
      real(WP) :: min_filmthickness      =1.0e-4_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP

      ! Hole edge indicator
      real(WP), dimension(:,:,:), allocatable   :: edge    !< Hole edge indicator
      real(WP), dimension(:,:,:), allocatable   :: edgeold !< Hole edge indicator

      ! Film drainage parameter - turned off by default
      real(WP) :: film_drain=0.0_WP                       !< Threshold VF parameter for film drainage procedure (0.0=off)

      ! !> Film measurements
      ! real(WP) :: min_thickness,film_volume,film_ratio

      !> Film retraction
      real(WP) :: edge_threshold !=0.1_WP
      real(WP), dimension(:,:,:), allocatable :: filmU,filmV,filmW ! Face-centered film velocities
      real(WP), dimension(:,:,:), allocatable :: filmUm,filmVm,filmWm ! Velocities at film barycenters      

   contains
      procedure :: initialize
      procedure :: compute_film_velocity
      procedure :: transport_edge_indicator
      ! procedure :: breakup_film_retraction
      ! procedure, private :: bag_droplet_gamma
      ! procedure :: get_min_thickness
      procedure :: detect_edge_regions
   end type film

contains


   !> Initialize
   subroutine initialize(this,cfg,vf,fs)
      use param, only: param_read
      implicit none
      class(film), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      type(vfs), target, intent(in) :: vf
      type(tpns), target, intent(in) :: fs

      this%cfg=>cfg
      this%vf=>vf
      this%fs=>fs

      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%max_interface_planes=2
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl

      ! ! Initialize film statistics
      ! this%min_thickness=0.0_WP
      ! this%film_volume=0.0_WP
      ! this%film_ratio=0.0_WP
      call param_read('Edge threshold',this%edge_threshold,default=0.1_WP)

      ! Allocate film edge indicator
      allocate(this%edge   (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%edge   =0.0_WP
      allocate(this%edgeold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); this%edgeold=0.0_WP

      ! Allocate film velocity
      allocate(this%filmU (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%filmV (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%filmW (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))         
      allocate(this%filmUm(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%filmVm(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(this%filmWm(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))  
      this%filmU=0.0_WP; this%filmV=0.0_WP; this%filmW=0.0_WP

   end subroutine initialize


   !> Taylor-Culick velocity
   subroutine compute_film_velocity(this)
      use vfs_class, only: r2p
      use irl_fortran_interface
      implicit none
      class(film), intent(inout) :: this
      integer :: m,n
      integer :: i,j,k
      real(WP), dimension(3) :: pos ! position where to measure height
      real(WP), dimension(3) :: mynorm ! alternative edge normal
      ! real(WP), dimension(3) :: ph3 ! disk/hole center position in 3D space
      real(WP) :: height, Utc
      real(WP) :: dir
      integer :: count
      integer :: index,ii,jj,kk,ndata
      real(WP) :: ww,wwUsum,wwVsum,wwWsum
      real(WP), dimension(3) :: prefU,prefV,prefW,ploc,buf


      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%cc%film_classify(this%vf%Lbary,this%vf%Gbary)
      if (this%vf%two_planes) then
         call this%vf%sense_interface()      
         call this%detect_edge_regions()
      end if
      !! Prescribe Taylor-Culick velocity based on film thickness near rim
      ! ph3=0.0_WP; ! ph3(1)=ph(1,1); ph3(3)=ph(2,1)
      dir=1.0_WP !-2.0_WP*real(interface_shape,WP) ! positive if hole, negative if disk
      this%filmU=0.0_WP; this%filmV=0.0_WP; this%filmW=0.0_WP
      this%filmUm=0.0_WP; this%filmVm=0.0_WP; this%filmWm=0.0_WP                   
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
         count=0
            do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)                        
            ! if (this%vf%edge_sensor(i,j,k).ge.this%edge_threshold.or.this%edge(i,j,k).ge.2.0_WP) then ! enforce i-2:i continuity too?
               if (this%edge(i,j,k).ge.1.0_WP) then
                  ! Get thickness away from edge
                  if (this%vf%VF(i,j,k).le.0.0_WP) cycle
                  mynorm=this%vf%edge_normal(:,i,j,k)
                  height=max(0.5_WP*this%vf%thickness(i,j,k),0.5_WP*this%min_filmthickness)
                  ! height=0.5e-6_WP
                  Utc=sqrt(this%fs%sigma/this%fs%rho_l/height) !*merge(1.0_WP,0.0_WP,this%edge(i,j,k).ge.1.0_WP)

                  this%filmUm(i,j,k)=-mynorm(1)*Utc
                  this%filmVm(i,j,k)=-mynorm(2)*Utc 
                  this%filmWm(i,j,k)=-mynorm(3)*Utc                              
                  count=count+1
            end if                         
            end do
      end do
      ! Deallocate CCL lists
      call this%cc%deallocate_lists()
      ! Sync velocities
      call this%vf%cfg%sync(this%filmUm)
      call this%vf%cfg%sync(this%filmVm)
      call this%vf%cfg%sync(this%filmWm)

      ! Propagate velocities outward
      do index=1,sum(this%vf%band_count(0:1))
      ! do index=sum(this%vf%band_count(0))+1,sum(this%vf%band_count(1))
         i=this%vf%band_map(1,index)
         j=this%vf%band_map(2,index)
         k=this%vf%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%vf%mask(i,j,k).ne.0) cycle                     
            
         ! ! Skip edge cells with velocities
         ! if (abs(this%filmUm(i,j,k))+abs(this%filmVm(i,j,k))+abs(this%filmWm(i,j,k))+abs(this%filmUm(i-1,j,k))+abs(this%filmVm(i,j-1,k))+abs(this%filmWm(i,j,k-1)).le.0.0_WP) cycle

         ! Get reference locations
         prefU=[this%vf%cfg%x (i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
         prefV=[this%vf%cfg%xm(i),this%vf%cfg%y (j),this%vf%cfg%zm(k)]
         prefW=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%z (k)]
         ndata=0; wwUsum=0.0_WP; wwVsum=0.0_WP; wwWsum=0.0_WP
         do kk=k-1,k+1
            do jj=j-1,j+1
               do ii=i-1,i+1

                  ! Skip the cell if it's a true wall
                  if (this%vf%mask(ii,jj,kk).eq.1) cycle

                  ! Skip the cell if there is zero velocity
                  if (abs(this%filmUm(ii,jj,kk))+abs(this%filmVm(ii,jj,kk))+abs(this%filmWm(ii,jj,kk)).le.0.0_WP) cycle

                  ! Calculate weighted velocity contribution
                  ploc=this%vf%Lbary(:,ii,jj,kk)
                  buf=ploc-prefU
                  ww=wgauss(sqrt(dot_product(buf,buf)),2.5_WP)
                  wwUsum=wwUsum+ww
                  this%filmU(i,j,k)=this%filmU(i,j,k)+this%filmUm(ii,jj,kk)*ww
                  buf=ploc-prefV
                  ww=wgauss(sqrt(dot_product(buf,buf)),2.5_WP)
                  wwVsum=wwVsum+ww
                  this%filmV(i,j,k)=this%filmV(i,j,k)+this%filmVm(ii,jj,kk)*ww
                  buf=ploc-prefW
                  ww=wgauss(sqrt(dot_product(buf,buf)),2.5_WP)
                  wwWsum=wwWsum+ww
                  this%filmW(i,j,k)=this%filmW(i,j,k)+this%filmWm(ii,jj,kk)*ww
                  ndata=ndata+1
               end do
            end do
         end do                
         this%filmU(i,j,k)=this%filmU(i,j,k)/(wwUsum+tiny(1.0_WP))
         this%filmV(i,j,k)=this%filmV(i,j,k)/(wwVsum+tiny(1.0_WP))
         this%filmW(i,j,k)=this%filmW(i,j,k)/(wwWsum+tiny(1.0_WP))
      end do
      ! Sync velocities
      call this%vf%cfg%sync(this%filmU)
      call this%vf%cfg%sync(this%filmV)
      call this%vf%cfg%sync(this%filmW)      

   contains

      ! Quasi-Gaussian weighting function - h=2.5 looks okay
      real(WP) function wgauss(d,h)
         implicit none
         real(WP), intent(in) :: d,h
         if (d.ge.h) then
            wgauss=0.0_WP
         else
            wgauss=(1.0_WP+4.0_WP*d/h)*(1.0_WP-d/h)**4
         end if
      end function wgauss

   end subroutine compute_film_velocity   


   !> Advect edge indicator
   !> Calculate the new VF based on U/V/W and dt
   subroutine transport_edge_indicator(this,dt,U,V,W)
      use mathtools, only: normalize
      implicit none
      class(film), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index,ii,jj,kk
      integer, parameter :: advect_band=1
      real(WP), dimension(3) :: ctr_now,dir
      integer, dimension(3) :: edge_ind

      ! Detect edges
      call this%detect_edge_regions()

      ! Remember old edge indicator
      this%edgeold=this%edge

      ! Transport edge indicator
      do index=1,sum(this%vf%band_count(0:advect_band))
         i=this%vf%band_map(1,index)
         j=this%vf%band_map(2,index)
         k=this%vf%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%vf%mask(i,j,k).ne.0) cycle               

         ! Semi-Lagrangian
         ! Project liquid barycenter backwards in time - get direction first
         ctr_now=this%vf%project(this%vf%Lbary(:,i,j,k),i,j,k,-dt,U,V,W)
         ! Get estimated edge cell index
         edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])         
         ! Make sure neighbor cell is reached by projecting a meshsize back if needed
         if (edge_ind(1).eq.i.and.edge_ind(2).eq.j.and.edge_ind(3).eq.k) then
            dir=normalize(ctr_now-this%vf%Lbary(:,i,j,k))
            ctr_now=this%vf%Lbary(:,i,j,k)+dir*this%cfg%meshsize(i,j,k)
            edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         end if      
         if (this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)).ge.1.0_WP) then
            this%edge(i,j,k)=ceiling(this%vf%VF(i,j,k))
         end if    
         
         ! Clear empty cells
         this%edge(i,j,k)=min(ceiling(this%vf%VF(i,j,k)),int(this%edge(i,j,k)))
      end do

      ! Check neighbors for edge cells
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%edge(i,j,k).eq.1) then
                  this%edge(i,j,k)=0
                  do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
                     this%edge(i,j,k)=max(this%edge(i,j,k),merge(1.0_WP,0.0_WP,this%vf%edge_sensor(ii,jj,kk).ge.0.5_WP*this%edge_threshold))
                  end do; end do; end do
               end if
            end do
         end do
      end do

      ! Synchronize edge field
      call this%cfg%sync(this%edge)

   end subroutine transport_edge_indicator


   !> Detect edge-like regions of the interface
   !> Alternative to vfs_class version
   subroutine detect_edge_regions(this)
      ! use mathtools, only: cross_product
      use vfs_class, only: VFlo,VFhi
      implicit none
      class(film), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,ni
      real(WP) :: w,fvol,myvol,volume_sensor,surface_sensor
      real(WP), dimension(3) :: fvbary,fbary,mybary,bary
      real(WP) :: wwidth=2.5_WP
      real(WP) :: wwidth2=4.5_WP
      real(WP), dimension(:,:,:)  , allocatable :: s_tmp
      real(WP), dimension(:,:,:,:), allocatable :: v_tmp
      real(WP) :: surface_area
      ! Default value is 0
      this%vf%edge_sensor=0.0_WP
      this%vf%edge_normal=0.0_WP
      ! Traverse domain and compute sensors
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               ! Skip wall/bcond cells
               if (this%vf%mask(i,j,k).ne.0) cycle
               ! Skip full cells
               if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
               ! Compute filtered volume barycenter sensor
               fvol=0.0_WP; fvbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        bary=this%vf%Lbary(:,ii,jj,kk)
                        w=wgauss(norm2(bary-this%vf%Lbary(:,i,j,k))/this%vf%cfg%meshsize(i,j,k),wwidth)
                        myvol=this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*w
                        fvol =fvol +myvol
                        fvbary=fvbary+myvol*this%vf%Lbary(:,ii,jj,kk)
                     end do
                  end do
               end do
               fvbary=fvbary/fvol
               volume_sensor=norm2(fvbary-this%vf%Lbary(:,i,j,k))/this%vf%cfg%meshsize(i,j,k)
               ! Compute filtered surface barycenter sensor
               fvol=0.0_WP; mybary=0.0_WP
               do ni=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(ni,i,j,k)).ne.0) then
                     myvol =abs(calculateVolume(this%vf%interface_polygon(ni,i,j,k)))
                     fvol  =fvol  +myvol
                     mybary=mybary+myvol*calculateCentroid(this%vf%interface_polygon(ni,i,j,k))
                  end if
               end do
               mybary=mybary/fvol
               fvol=0.0_WP; fbary=0.0_WP
               do kk=k-2,k+2
                  do jj=j-2,j+2
                     do ii=i-2,i+2
                        do ni=1,getNumberOfPlanes(this%vf%liquid_gas_interface(ii,jj,kk))
                           if (getNumberOfVertices(this%vf%interface_polygon(ni,ii,jj,kk)).ne.0) then
                              bary=calculateCentroid(this%vf%interface_polygon(ni,ii,jj,kk))
                              w=wgauss(norm2(bary-mybary)/this%vf%cfg%meshsize(i,j,k),wwidth)
                              myvol=abs(calculateVolume(this%vf%interface_polygon(ni,ii,jj,kk)))*w
                              fvol =fvol +myvol
                              fbary=fbary+myvol*calculateCentroid(this%vf%interface_polygon(ni,ii,jj,kk))
                           end if
                        end do
                     end do
                  end do
               end do
               fbary=fbary/fvol
               surface_sensor=norm2(fbary-mybary)/this%vf%cfg%meshsize(i,j,k)
               ! Aggregate into an edge sensor
               this%vf%edge_sensor(i,j,k)=volume_sensor*surface_sensor
               ! Clip based on thickness
               if (this%vf%thickness(i,j,k).gt.this%vf%thin_thld_max*this%vf%cfg%meshsize(i,j,k)) this%vf%edge_sensor(i,j,k)=0.0_WP
               ! Finally, store edge orientation
               if (this%vf%edge_sensor(i,j,k).gt.0.0_WP) then
                  ! this%vf%edge_normal(:,i,j,k)=-(fbary-mybary)/(norm2(fbary-mybary)+epsilon(1.0_WP))
                  this%vf%edge_normal(:,i,j,k)=-(fvbary-this%vf%Lbary(:,i,j,k))/(norm2(fvbary-this%vf%Lbary(:,i,j,k))+epsilon(1.0_WP))
               end if
            end do
         end do
      end do
      ! Communicate
      call this%vf%cfg%sync(this%vf%edge_sensor)
      call this%vf%cfg%sync(this%vf%edge_normal)
      ! ! Apply an extra step of surface smoothing to our edge info
      ! allocate(s_tmp(    this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); s_tmp=0.0_WP
      ! allocate(v_tmp(1:3,this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); v_tmp=0.0_WP
      ! do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
      !    do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
      !       do i=this%vf%cfg%imin_,this%vf%cfg%imax_
      !          ! Skip wall/bcond/full cells
      !          if (this%vf%mask(i,j,k).ne.0) cycle
      !          if (this%vf%VF(i,j,k).lt.VFlo.or.this%vf%VF(i,j,k).gt.VFhi) cycle
      !          ! Surface-averaged normal magnitude
      !          surface_area=0.0_WP
      !          do kk=k-1,k+1; do jj=j-1,j+1; do ii=i-1,i+1
      !             bary=this%vf%Lbary(:,ii,jj,kk)
      !             w=wgauss(norm2(bary-this%vf%Lbary(:,i,j,k))/this%vf%cfg%meshsize(i,j,k),wwidth2) !*norm2(cross_product(this%vf%edge_normal(:,i,j,k),this%vf%Lbary(:,ii,jj,kk)-this%vf%Lbary(:,i,j,k)))
      !             surface_area  =  surface_area+this%vf%SD(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*w
      !             s_tmp  (i,j,k)=s_tmp  (i,j,k)+this%vf%SD(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%edge_sensor  (ii,jj,kk)*w
      !             v_tmp(:,i,j,k)=v_tmp(:,i,j,k)+this%vf%SD(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%edge_normal(:,ii,jj,kk)*w
      !          end do; end do; end do
      !          if (surface_area.gt.0.0_WP) then
      !             s_tmp  (i,j,k)=s_tmp  (i,j,k)/surface_area
      !             v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/surface_area
      !             v_tmp(:,i,j,k)=v_tmp(:,i,j,k)/(norm2(v_tmp(:,i,j,k))+epsilon(1.0_WP))
      !          end if
      !       end do
      !    end do
      ! end do
      ! call this%vf%cfg%sync(s_tmp); this%vf%edge_sensor=s_tmp; deallocate(s_tmp)
      ! call this%vf%cfg%sync(v_tmp); this%vf%edge_normal=v_tmp; deallocate(v_tmp)
      
      contains

      ! Quasi-Gaussian weighting function - h=2.5 looks okay
      real(WP) function wgauss(r,h)
         implicit none
         real(WP), intent(in) :: r,h
         if (r.ge.h) then
            wgauss=0.0_WP
         else
            wgauss=(1.0_WP+4.0_WP*r/h)*(1.0_WP-r/h)**4
         end if
         ! wgauss=1.0_WP
      end function wgauss

   end subroutine detect_edge_regions

end module film_class
