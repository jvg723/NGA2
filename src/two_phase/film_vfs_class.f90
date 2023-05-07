!> Volume fraction solver class:
!> Provides support for various BC, semi-Lagrangian geometric advancement,
!> curvature calculation, interface reconstruction.
module film_vfs_class
   use precision,      only: WP
   use vfs_class,      only: vfs
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   use irl_fortran_interface
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: film_vfs,bcond
   
   ! Also expose the min and max VF values
   public :: VFhi,VFlo
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available interface reconstructions schemes for VF
   integer, parameter, public :: lvira=1             !< LVIRA scheme
   integer, parameter, public :: elvira=2            !< ELVIRA scheme
   !integer, parameter, public :: mof=3               !< MOF scheme
   integer, parameter, public :: r2p=4               !< R2P scheme
   integer, parameter, public :: swartz=5            !< Swartz scheme
   integer, parameter, public :: art=6               !< ART scheme
   integer, parameter, public :: youngs=7            !< Youngs' scheme
   integer, parameter, public :: lvlset=8            !< Levelset-based scheme
   
   ! IRL cutting moment calculation method
   integer, parameter, public :: recursive_simplex=0 !< Recursive simplex cutting
   integer, parameter, public :: half_edge=1         !< Half-edge cutting
   integer, parameter, public :: nonrecurs_simplex=2 !< Non-recursive simplex cutting
   
   ! Default parameters for volume fraction solver
   integer,  parameter :: nband=3                                 !< Number of cells around the interfacial cells on which localized work is performed
   integer,  parameter :: advect_band=2                           !< How far we do the transport
   integer,  parameter :: distance_band=2                         !< How far we build the distance
   integer,  parameter :: max_interface_planes=2                  !< Maximum number of interfaces allowed (2 for R2P)
   real(WP), parameter :: VFlo=1.0e-12_WP                         !< Minimum VF value considered
   real(WP), parameter :: VFhi=1.0_WP-VFlo                        !< Maximum VF value considered
   real(WP), parameter :: volume_epsilon_factor =1.0e-15_WP       !< Minimum volume  to consider for computational geometry (normalized by min_meshsize**3)
   real(WP), parameter :: surface_epsilon_factor=1.0e-15_WP       !< Minimum surface to consider for computational geometry (normalized by min_meshsize**2)
   real(WP), parameter :: iterative_distfind_tol=1.0e-12_WP       !< Tolerance for iterative plane distance finding
   
   !> Boundary conditions for the volume fraction solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Volume fraction solver object definition
   type, extends(vfs) :: film_vfs
      
      ! Curvature
      real(WP), dimension(:,:,:,:), allocatable :: curv2  !< Interface mean curvatures for two-plane cells

      ! Hole edge indicator
      real(WP), dimension(:,:,:), allocatable   :: edge    !< Hole edge indicator
      real(WP), dimension(:,:,:), allocatable   :: edgeold !< Hole edge indicator

      ! Film drainage parameter - turned off by default
      real(WP) :: film_drain=0.0_WP                       !< Threshold VF parameter for film drainage procedure (0.0=off)
      
   contains

      procedure :: advance                                !< Advance VF to next step
      procedure :: advance_film                           !< Advance VF to next step with a non-solenoidal velocity
      procedure, private :: crude_phase_test              !< Helper function that rapidly assess if a mixed cell might be present
      procedure :: get_curvature                          !< Compute curvature from IRL surface polygons

   end type film_vfs
   
   !> Declare volume fraction solver constructor
   interface film_vfs
      procedure constructor
   end interface film_vfs
   
contains
   
   
   !> Default constructor for volume fraction solver
   function constructor(cfg,reconstruction_method,name) result(self)
      use messager, only: die
      implicit none
      type(film_vfs) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: reconstruction_method
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Check that we have at least 3 overlap cells - we can push that to 2 with limited work!
      if (cfg%no.lt.3) call die('[film_vfs constructor] The config requires at least 3 overlap cells')
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%VF     (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF     =0.0_WP
      allocate(self%VFold  (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFold  =0.0_WP
      allocate(self%Lbary  (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Lbary  =0.0_WP
      allocate(self%Gbary  (3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Gbary  =0.0_WP
      allocate(self%SD     (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SD     =0.0_WP
      allocate(self%G      (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%G      =0.0_WP
      allocate(self%curv   (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%curv   =0.0_WP
      allocate(self%curv2  (2,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%curv2  =0.0_WP
      allocate(self%edge   (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%edge   =0.0_WP
      allocate(self%edgeold(  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%edgeold=0.0_WP
      allocate(self%type   (  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%type   =0.0_WP

      ! Set clipping distance
      self%Gclip=real(distance_band+1,WP)*self%cfg%min_meshsize
      
      ! Subcell phasic volumes
      allocate(self%Lvol(0:1,0:1,0:1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Lvol=0.0_WP
      allocate(self%Gvol(0:1,0:1,0:1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Gvol=0.0_WP
      
      ! Prepare the band arrays
      allocate(self%band(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%band=0
      if (allocated(self%band_map)) deallocate(self%band_map)
      
      ! Set reconstruction method
      self%reconstruction_method=reconstruction_method
      
      ! Initialize IRL
      call self%initialize_irl()
      
      ! Prepare mask for VF
      allocate(self%mask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(self%cfg%imax+1:,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(:,self%cfg%jmax+1:,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(:,:,self%cfg%kmax+1:)=2
      end if
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask)
      
      ! Prepare mask for vertices
      allocate(self%vmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%vmask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%vmask(               :self%cfg%imin,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%vmask(self%cfg%imax+1:             ,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%vmask(:,               :self%cfg%jmin,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%vmask(:,self%cfg%jmax+1:             ,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%vmask(:,:,               :self%cfg%kmin)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%vmask(:,:,self%cfg%kmax+1:             )=2
      end if
      do k=self%cfg%kmino_+1,self%cfg%kmaxo_
         do j=self%cfg%jmino_+1,self%cfg%jmaxo_
            do i=self%cfg%imino_+1,self%cfg%imaxo_
               if (minval(self%cfg%VF(i-1:i,j-1:j,k-1:k)).eq.0.0_WP) self%vmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%vmask)
      if (.not.self%cfg%xper.and.self%cfg%iproc.eq.1) self%vmask(self%cfg%imino,:,:)=self%vmask(self%cfg%imino+1,:,:)
      if (.not.self%cfg%yper.and.self%cfg%jproc.eq.1) self%vmask(:,self%cfg%jmino,:)=self%vmask(:,self%cfg%jmino+1,:)
      if (.not.self%cfg%zper.and.self%cfg%kproc.eq.1) self%vmask(:,:,self%cfg%kmino)=self%vmask(:,:,self%cfg%kmino+1)
      
   end function constructor
    

   !> Calculate the new VF based on U/V/W and dt
   subroutine advance(this,dt,U,V,W)
      use mathtools, only: normalize
      implicit none
      class(film_vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,index,ii,jj,kk
      real(IRL_double), dimension(3,9) :: face
      type(CapDod_type) :: flux_polyhedron
      real(WP) :: Lvolold,Gvolold
      real(WP) :: Lvolinc,Gvolinc
      real(WP) :: Lvolnew,Gvolnew
      real(WP) :: vol_now,crude_VF
      real(WP), dimension(3) :: ctr_now,dir
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      integer, dimension(3) :: edge_ind

      ! Allocate
      call new(flux_polyhedron)
      
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (minval(abs(this%band(i-1:i,j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(1,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(1,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Y flux
               if (minval(abs(this%band(i,j-1:j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(2,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(2,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Z flux
               if (minval(abs(this%band(i,j,k-1:k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j+1,k  ).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(3,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(3,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
            end do
         end do
      end do
      
      ! Compute transported moments
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle
         
         ! Old liquid and gas volumes
         Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
         Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
         
         ! Compute incoming liquid and gas volumes
         Lvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),0)+getVolumePtr(this%face_flux(1,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),0)+getVolumePtr(this%face_flux(2,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),0)+getVolumePtr(this%face_flux(3,i,j,k),0)
         Gvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),1)+getVolumePtr(this%face_flux(1,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),1)+getVolumePtr(this%face_flux(2,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),1)+getVolumePtr(this%face_flux(3,i,j,k),1)
         
         ! Compute new liquid and gas volumes
         Lvolnew=Lvolold+Lvolinc
         Gvolnew=Gvolold+Gvolinc
         
         ! Compute new liquid volume fraction
         this%VF(i,j,k)=Lvolnew/(Lvolnew+Gvolnew)
         
         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Compute old phase barycenters
            this%Lbary(:,i,j,k)=(this%Lbary(:,i,j,k)*Lvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),0)+getCentroidPtr(this%face_flux(1,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),0)+getCentroidPtr(this%face_flux(2,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),0)+getCentroidPtr(this%face_flux(3,i,j,k),0))/Lvolnew
            this%Gbary(:,i,j,k)=(this%Gbary(:,i,j,k)*Gvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),1)+getCentroidPtr(this%face_flux(1,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),1)+getCentroidPtr(this%face_flux(2,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),1)+getCentroidPtr(this%face_flux(3,i,j,k),1))/Gvolnew
            ! Project forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
      end do
      
      ! Synchronize VF field
      call this%cfg%sync(this%VF)
      
      ! Remove flotsams if needed
      call this%remove_flotsams()
      
      ! Synchronize and clean-up barycenter fields
      call this%sync_and_clean_barycenters()
      
      ! Advect interface polygons
      call this%advect_interface(dt,U,V,W)
      
      ! Update the band
      call this%update_band()
      
      ! Perform interface reconstruction from transported moments
      call this%build_interface()
      
      ! Remove thin sheets if needed
      call this%remove_sheets()
      
      ! Create discontinuous polygon mesh from IRL interface
      call this%polygonalize_interface()
      
      ! Calculate distance from polygons
      call this%distance_from_polygon()
      
      ! Calculate subcell phasic volumes
      call this%subcell_vol()
      
      ! Calculate curvature
      call this%get_curvature()
      
      ! Reset moments to guarantee compatibility with interface reconstruction
      call this%reset_volume_moments()

      ! Remember old edge indicator
      this%edgeold=this%edge

      ! Transport edge indicator
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle               

         ! do kk=k-1,k+1
         !    do jj=j-1,j+1
         !       do ii=i-1,i+1
         ! ! do kk=k-2,k+2
         ! !    do jj=j-2,j+2
         ! !       do ii=i-2,i+2             
         !          if (this%VF(ii,jj,kk).eq.0.0_WP.and.this%VFold(ii,jj,kk).gt.VFlo.and.this%edgeold(ii,jj,kk).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%VF(i,j,k))
         !             ! this%edge(i,j,k)=1
         !          end if                  
         !       end do
         !    end do
         ! end do
         ! Semi-Lagrangian
         ! Project liquid barycenter backwards in time - get direction first
         ctr_now=this%project(this%Lbary(:,i,j,k),i,j,k,-dt,U,V,W)
         ! Get estimated edge cell index
         edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])         
         ! Make sure neighbor cell is reached by projecting a meshsize back if needed
         if (edge_ind(1).eq.i.and.edge_ind(2).eq.j.and.edge_ind(3).eq.k) then
            dir=normalize(ctr_now-this%Lbary(:,i,j,k))
            ctr_now=this%Lbary(:,i,j,k)+dir*this%cfg%meshsize(i,j,k)
            edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         end if      
         if (this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)).ge.1.0_WP) then
            this%edge(i,j,k)=ceiling(this%VF(i,j,k))
            ! this%edge(i,j,k)=1
         end if    
         ! this%edge(i,j,k)=max(this%edgeold(i,j,k),this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)))
         ! if (this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)).ge.1.0_WP) print *,'adv1 ijk',i,j,k,'edge_ind',edge_ind,'lbary',this%Lbary(:,i,j,k),'ctr_now',ctr_now,'edge ijk',this%edge(i,j,k)
         ! ! if (i.ne.edge_ind(1)+1) print *,'ijk',i,j,k,'edge_ind',edge_ind,'Lbary',this%Lbary(:,i,j,k),'proj',ctr_now,'edge ijk',this%edge(i,j,k),'edgeold edge_ind',this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)),'edge i-1jk',this%edge(i-1,j,k),'edgeold i-1jk',this%edgeold(i-1,j,k),'edgeold ijk',this%edgeold(i,j,k)
         ! ! if (this%edge(i,j,k).ne.this%edgeold(i,j,k)) print *,'ijk',i,j,k,'edge_ind',edge_ind,'edgeold ijk',this%edgeold(i,j,k),'edgeold edge_ind',this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)),'edge ijk',this%edge(i,j,k),'edge edge_ind',this%edge(edge_ind(1),edge_ind(2),edge_ind(3))
         ! ! this%edge(i,j,k)=ceiling(this%VF(i,j,k))*maxval(this%edgeold(edge_ind(1)-1:edge_ind(1)+1,edge_ind(2)-1:edge_ind(2)+1,edge_ind(3)-1:edge_ind(3)+1))
         ! this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
         ! ! If any neighbor drains, then use interpolated edgeold in projected direction
         ! do kk=k-1,k+1
         !    do jj=j-1,j+1
         !       do ii=i-1,i+1
         !          if (this%VF(ii,jj,kk).eq.0.0_WP.and.this%VFold(ii,jj,kk).gt.VFlo.and.this%edgeold(ii,jj,kk).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
         !          end if
         !       end do
         !    end do
         ! end do         
         ! Clear empty cells
         this%edge(i,j,k)=min(ceiling(this%VF(i,j,k)),int(this%edge(i,j,k)))
      end do

      ! Synchronize edge field
      call this%cfg%sync(this%edge)
      
   end subroutine advance


   !> Calculate the new VF based on U/V/W and dt, where velocity is not divergence-free
   subroutine advance_film(this,dt,U,V,W,reconstruct_flag)
      use mathtools, only: normalize
      implicit none
      class(film_vfs), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      logical, intent(in), optional :: reconstruct_flag
      integer :: i,j,k,index,ii,jj,kk
      real(IRL_double), dimension(3,9) :: face
      type(CapDod_type) :: flux_polyhedron
      real(WP) :: Lvolold,Gvolold
      real(WP) :: Lvolinc,Gvolinc
      real(WP) :: Lvolnew,Gvolnew
      real(WP) :: vol_now,crude_VF
      real(WP), dimension(3) :: ctr_now,dir
      real(WP), dimension(3,2) :: bounding_pts
      integer, dimension(3,2) :: bb_indices
      real(WP) :: neg_flux_sum
      logical :: reconstruct
      integer, dimension(3) :: edge_ind, forw_ind
      
      ! Allocate
      call new(flux_polyhedron)
      ! print *,'rank',this%cfg%rank,'Start flux'
      ! Loop over the domain and compute fluxes using semi-Lagrangian algorithm
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               
               ! X flux
               if (minval(abs(this%band(i-1:i,j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*U(i,j,k)*this%cfg%dy(j)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(1,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(1,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Y flux
               if (minval(abs(this%band(i,j-1:j,k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k+1).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k+1)]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k+1).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*V(i,j,k)*this%cfg%dx(i)*this%cfg%dz(k))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(2,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(2,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
               ! Z flux
               if (minval(abs(this%band(i,j,k-1:k))).le.advect_band) then
                  ! Construct and project face
                  face(:,1)=[this%cfg%x(i  ),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,5)=this%project(face(:,1),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j+1,k  ).eq.1) face(:,5)=face(:,1)
                  face(:,2)=[this%cfg%x(i  ),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,6)=this%project(face(:,2),i,j,k,-dt,U,V,W); if (this%vmask(i  ,j  ,k  ).eq.1) face(:,6)=face(:,2)
                  face(:,3)=[this%cfg%x(i+1),this%cfg%y(j  ),this%cfg%z(k  )]; face(:,7)=this%project(face(:,3),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j  ,k  ).eq.1) face(:,7)=face(:,3)
                  face(:,4)=[this%cfg%x(i+1),this%cfg%y(j+1),this%cfg%z(k  )]; face(:,8)=this%project(face(:,4),i,j,k,-dt,U,V,W); if (this%vmask(i+1,j+1,k  ).eq.1) face(:,8)=face(:,4)
                  face(:,9)=0.25_WP*[sum(face(1,1:4)),sum(face(2,1:4)),sum(face(3,1:4))]
                  face(:,9)=this%project(face(:,9),i,j,k,-dt,U,V,W)
                  ! Form flux polyhedron
                  call construct(flux_polyhedron,face)
                  ! Add solenoidal correction
                  call adjustCapToMatchVolume(flux_polyhedron,dt*W(i,j,k)*this%cfg%dx(i)*this%cfg%dy(j))
                  ! Get bounds for flux polyhedron
                  call getBoundingPts(flux_polyhedron,bounding_pts(:,1),bounding_pts(:,2))
                  bb_indices(:,1)=this%cfg%get_ijk_local(bounding_pts(:,1),[i,j,k])
                  bb_indices(:,2)=this%cfg%get_ijk_local(bounding_pts(:,2),[i,j,k])
                  ! Crudely check phase information for flux polyhedron
                  crude_VF=this%crude_phase_test(bb_indices)
                  if (crude_VF.lt.0.0_WP) then
                     ! Need full geometric flux
                     call getMoments(flux_polyhedron,this%localized_separator_link(i,j,k),this%face_flux(3,i,j,k))
                  else
                     ! Simpler flux calculation
                     vol_now=calculateVolume(flux_polyhedron); ctr_now=calculateCentroid(flux_polyhedron)
                     call construct(this%face_flux(3,i,j,k),[crude_VF*vol_now,crude_VF*vol_now*ctr_now,(1.0_WP-crude_VF)*vol_now,(1.0_WP-crude_VF)*vol_now*ctr_now])
                  end if
               end if
               
            end do
         end do
      end do
      ! print *,'rank',this%cfg%rank,'Compute transported moments'
      ! Compute transported moments
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle
         
         ! Old liquid and gas volumes
         Lvolold=        this%VFold(i,j,k) *this%cfg%vol(i,j,k)
         Gvolold=(1.0_WP-this%VFold(i,j,k))*this%cfg%vol(i,j,k)
         
         ! Compute incoming liquid and gas volumes
         Lvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),0)+getVolumePtr(this%face_flux(1,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),0)+getVolumePtr(this%face_flux(2,i,j,k),0) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),0)+getVolumePtr(this%face_flux(3,i,j,k),0)
         Gvolinc=-getVolumePtr(this%face_flux(1,i+1,j,k),1)+getVolumePtr(this%face_flux(1,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(2,i,j+1,k),1)+getVolumePtr(this%face_flux(2,i,j,k),1) &
         &       -getVolumePtr(this%face_flux(3,i,j,k+1),1)+getVolumePtr(this%face_flux(3,i,j,k),1)
         
         ! Compute new liquid and gas volumes
         Lvolnew=Lvolold+Lvolinc
         Gvolnew=Gvolold+Gvolinc
         
         ! Compute new liquid volume fraction
         this%VF(i,j,k)=Lvolnew/this%cfg%vol(i,j,k) !(Lvolnew+Gvolnew)
         
         ! Distribute trailing volume
         if (Lvolinc.lt.0.0_WP.and.this%VF(i,j,k).lt.this%film_drain) then
            neg_flux_sum=max(getVolumePtr(this%face_flux(1,i+1,j,k),0),0.0_WP)-min(getVolumePtr(this%face_flux(1,i,j,k),0),0.0_WP) &
            &           +max(getVolumePtr(this%face_flux(2,i,j+1,k),0),0.0_WP)-min(getVolumePtr(this%face_flux(2,i,j,k),0),0.0_WP) &
            &           +max(getVolumePtr(this%face_flux(3,i,j,k+1),0),0.0_WP)-min(getVolumePtr(this%face_flux(3,i,j,k),0),0.0_WP)
            this%VF(i-1,j,k)=this%VF(i-1,j,k)+this%VF(i,j,k)*(-min(getVolumePtr(this%face_flux(1,i,j,k),0),0.0_WP)/neg_flux_sum)            
            this%VF(i+1,j,k)=this%VF(i+1,j,k)+this%VF(i,j,k)*max(getVolumePtr(this%face_flux(1,i+1,j,k),0),0.0_WP)/neg_flux_sum            
            this%VF(i,j-1,k)=this%VF(i,j-1,k)+this%VF(i,j,k)*(-min(getVolumePtr(this%face_flux(2,i,j,k),0),0.0_WP)/neg_flux_sum)            
            this%VF(i,j+1,k)=this%VF(i,j+1,k)+this%VF(i,j,k)*max(getVolumePtr(this%face_flux(2,i,j+1,k),0),0.0_WP)/neg_flux_sum            
            this%VF(i,j,k-1)=this%VF(i,j,k-1)+this%VF(i,j,k)*(-min(getVolumePtr(this%face_flux(3,i,j,k),0),0.0_WP)/neg_flux_sum)            
            this%VF(i,j,k+1)=this%VF(i,j,k+1)+this%VF(i,j,k)*max(getVolumePtr(this%face_flux(3,i,j,k+1),0),0.0_WP)/neg_flux_sum            
            this%VF(i,j,k)=0.0_WP
         end if

         ! Only work on higher order moments if VF is in [VFlo,VFhi]
         if (this%VF(i,j,k).lt.VFlo) then
            this%VF(i,j,k)=0.0_WP
         else if (this%VF(i,j,k).gt.VFhi) then
            this%VF(i,j,k)=1.0_WP
         else
            ! Compute old phase barycenters
            this%Lbary(:,i,j,k)=(this%Lbary(:,i,j,k)*Lvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),0)+getCentroidPtr(this%face_flux(1,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),0)+getCentroidPtr(this%face_flux(2,i,j,k),0) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),0)+getCentroidPtr(this%face_flux(3,i,j,k),0))/Lvolnew
            this%Gbary(:,i,j,k)=(this%Gbary(:,i,j,k)*Gvolold-getCentroidPtr(this%face_flux(1,i+1,j,k),1)+getCentroidPtr(this%face_flux(1,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(2,i,j+1,k),1)+getCentroidPtr(this%face_flux(2,i,j,k),1) &
            &                                               -getCentroidPtr(this%face_flux(3,i,j,k+1),1)+getCentroidPtr(this%face_flux(3,i,j,k),1))/Gvolnew
            ! Project forward in time
            this%Lbary(:,i,j,k)=this%project(this%Lbary(:,i,j,k),i,j,k,dt,U,V,W)
            this%Gbary(:,i,j,k)=this%project(this%Gbary(:,i,j,k),i,j,k,dt,U,V,W)
         end if
      end do
      
      ! Synchronize VF field
      call this%cfg%sync(this%VF)
      
      ! Remove flotsams if needed
      call this%remove_flotsams()
      
      ! Synchronize and clean-up barycenter fields
      call this%sync_and_clean_barycenters()
      ! if (this%cfg%amRoot) print *,'rank',this%cfg%rank,'Advect polygons'

      ! Advect interface polygons
      call this%advect_interface(dt,U,V,W)
      
      ! Update the band
      call this%update_band()
      ! if (this%cfg%amRoot) print *,'rank',this%cfg%rank,'Start reconstruction'

      ! reconstruct=.true.
      ! if (present(reconstruct_flag)) reconstruct=reconstruct_flag
      ! if (reconstruct) then  
         ! Perform interface reconstruction from transported moments
         ! if (present(reconstruct_flag)) then
         !    call this%build_r2pdebug()
         ! else
            call this%build_interface()
         ! end if
         ! if (this%cfg%amRoot) print *,'rank',this%cfg%rank,'End reconstruction'

         ! Remove thin sheets if needed
         call this%remove_sheets()
         
         ! Create discontinuous polygon mesh from IRL interface
         call this%polygonalize_interface()
         
         ! Calculate distance from polygons
         call this%distance_from_polygon()
         
         ! Calculate subcell phasic volumes
         call this%subcell_vol()
         
         ! Calculate curvature
         call this%get_curvature()
         
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%reset_volume_moments()
      ! end if
      ! if (this%cfg%amRoot) print *,'rank',this%cfg%rank,'Start edge transport'

      ! Remember old edge indicator
      this%edgeold=this%edge

      ! Transport edge indicator
      do index=1,sum(this%band_count(0:advect_band))
         i=this%band_map(1,index)
         j=this%band_map(2,index)
         k=this%band_map(3,index)
         
         ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
         if (this%mask(i,j,k).ne.0) cycle               

         ! do kk=k-1,k+1
         !    do jj=j-1,j+1
         !       do ii=i-1,i+1
         ! ! do kk=k-2,k+2
         ! !    do jj=j-2,j+2
         ! !       do ii=i-2,i+2             
         !          if (this%VF(ii,jj,kk).eq.0.0_WP.and.this%VFold(ii,jj,kk).gt.VFlo.and.this%edgeold(ii,jj,kk).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%VF(i,j,k))
         !             ! this%edge(i,j,k)=1
         !          end if                  
         !       end do
         !    end do
         ! end do
         ! Semi-Lagrangian
         ! Project liquid barycenter backward in time
         ctr_now=this%project(this%Lbary(:,i,j,k),i,j,k,-dt,U,V,W)
         ! Get estimated edge cell index
         edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         ! Make sure neighbor cell is reached by projecting a meshsize back if needed
         if (edge_ind(1).eq.i.and.edge_ind(2).eq.j.and.edge_ind(3).eq.k) then
            dir=normalize(ctr_now-this%Lbary(:,i,j,k))
            ctr_now=this%Lbary(:,i,j,k)+dir*this%cfg%meshsize(i,j,k)
            edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         end if
          
         ! this%edge(i,j,k)=max(this%edgeold(i,j,k),this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)))
         ! if (this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)).ge.1.0_WP) print *,'adv2 ijk',i,j,k,'edge_ind',edge_ind,'lbary',this%Lbary(:,i,j,k),'ctr_now',ctr_now,'edge ijk',this%edge(i,j,k)
         ! if (edge_ind(2).ne.j.or.edge_ind(3).ne.k) print *,'adv2 moved ijk',i,j,k,'edge_ind',edge_ind,'lbary',this%Lbary(:,i,j,k),'ctr_now',ctr_now,'edge ijk',this%edge(i,j,k)
         ! if (this%edge(i,j,k).ne.this%edgeold(i,j,k)) print *,'ijk',i,j,k,'edge_ind',edge_ind,'edgeold ijk',this%edgeold(i,j,k),'edgeold edge_ind',this%edgeold(edge_ind(1),edge_ind(2),edge_ind(3)),'edge ijk',this%edge(i,j,k),'edge edge_ind',this%edge(edge_ind(1),edge_ind(2),edge_ind(3))
         ! this%edge(i,j,k)=ceiling(this%VF(i,j,k))*maxval(this%edgeold(edge_ind(1)-1:edge_ind(1)+1,edge_ind(2)-1:edge_ind(2)+1,edge_ind(3)-1:edge_ind(3)+1))
         ! If any neighbor drains, then use interpolated edgeold in projected direction
         do kk=k-1,k+1
            do jj=j-1,j+1
               do ii=i-1,i+1
                  if (this%VF(ii,jj,kk).eq.0.0_WP.and.this%VFold(ii,jj,kk).gt.VFlo.and.this%edgeold(ii,jj,kk).ge.1.0_WP) then
                     this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
                  end if
               end do
            end do
         end do

         ! ! Project liquid barycenter 2x backward in time
         ! ctr_now=this%project(this%Lbary(:,i,j,k),i,j,k,-2.0_WP*dt,U,V,W)
         ! ! Get estimated edge cell index
         ! edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         ! ! Make sure neighbor cell is reached by projecting a meshsize back if needed
         ! if (edge_ind(1).eq.i.and.edge_ind(2).eq.j.and.edge_ind(3).eq.k) then
         !    dir=normalize(ctr_now-this%Lbary(:,i,j,k))
         !    ctr_now=this%Lbary(:,i,j,k)+dir*this%cfg%meshsize(i,j,k)
         !    edge_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         ! end if
         ! ! Project liquid barycenter 2x forward in time
         ! ctr_now=this%project(this%Lbary(:,i,j,k),i,j,k,2.0_WP*dt,U,V,W)
         ! ! Get estimated edge cell index
         ! forw_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         ! ! Make sure neighbor cell is reached by projecting a meshsize forward if needed
         ! if (forw_ind(1).eq.i.and.forw_ind(2).eq.j.and.forw_ind(3).eq.k) then
         !    dir=normalize(ctr_now-this%Lbary(:,i,j,k))
         !    ctr_now=this%Lbary(:,i,j,k)+dir*this%cfg%meshsize(i,j,k)
         !    forw_ind=this%cfg%get_ijk_local(ctr_now,[i,j,k])
         ! end if
         ! ! If forward and backward neighbor cells are not empty, then self is not an edge cell
         ! if (this%VF(edge_ind(1),edge_ind(2),edge_ind(3)).gt.0.0_WP.and.this%VF(forw_ind(1),forw_ind(2),forw_ind(3)).gt.0.0_WP) then
         !    ! this%edge(i,j,k)=ceiling(this%VF(i,j,k))
         !    this%edge(i,j,k)=0.0_WP
         ! end if

         ! ! Directional version
         ! select case (maxloc(abs(dir),1))
         ! case (1)
         !    do kk=k-1,k+1
         !       do jj=j-1,j+1 
         !          if (this%VF(i,jj,kk).eq.0.0_WP.and.this%VFold(i,jj,kk).gt.VFlo.and.this%edgeold(i,jj,kk).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
         !          end if
         !       end do
         !    end do
         ! case (2)
         !    do kk=k-1,k+1
         !       do ii=i-1,i+1
         !          if (this%VF(ii,j,kk).eq.0.0_WP.and.this%VFold(ii,j,kk).gt.VFlo.and.this%edgeold(ii,j,kk).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
         !          end if
         !       end do
         !    end do
         ! case (3)
         !    do jj=j-1,j+1
         !       do ii=i-1,i+1
         !          if (this%VF(ii,jj,k).eq.0.0_WP.and.this%VFold(ii,jj,k).gt.VFlo.and.this%edgeold(ii,jj,k).ge.1.0_WP) then
         !             this%edge(i,j,k)=ceiling(this%cfg%get_scalar(ctr_now,i,j,k,this%edgeold,'n'))
         !          end if
         !       end do
         !    end do
         ! end select        
         ! Clear empty cells
         this%edge(i,j,k)=min(ceiling(this%VF(i,j,k)),int(this%edge(i,j,k)))
      end do

      ! Synchronize edge field
      call this%cfg%sync(this%edge)


   end subroutine advance_film   


   !> Private function to rapidly assess if a mixed cell is possible
   pure function crude_phase_test(this,b_ind) result(crude_phase)
      implicit none
      class(film_vfs), intent(in) :: this
      integer, dimension(3,2), intent(in) :: b_ind
      real(WP) :: crude_phase
      integer :: i,j,k
      ! Check if originating cell is mixed, continue if not
      crude_phase=this%VF(b_ind(1,2),b_ind(2,2),b_ind(3,2))
      if (crude_phase.ge.VFlo.and.crude_phase.le.VFhi) then
         ! Already have a mixed cell, we need the full geometry
         crude_phase=-1.0_WP; return
      end if
      ! Check cells in bounding box
      do k=b_ind(3,1),b_ind(3,2)
         do j=b_ind(2,1),b_ind(2,2)
            do i=b_ind(1,1),b_ind(1,2)
               if (this%VF(i,j,k).ne.crude_phase) then
                  ! We could have changed phase, we need the full geometry
                  crude_phase=-1.0_WP; return
               end if
            end do
         end do
      end do
      ! Ensure proper values
      if (crude_phase.gt.VFhi) then
         crude_phase=1.0_WP
      else if (crude_phase.lt.VFlo) then
         crude_phase=0.0_WP
      end if
   end function crude_phase_test


   !> Compute curvature from a least squares fit of the IRL surface
   subroutine get_curvature(this)
      implicit none
      class(film_vfs), intent(inout) :: this
      integer :: i,j,k,n
      real(WP), dimension(max_interface_planes) :: mycurv,mysurf
      real(WP), dimension(max_interface_planes,3) :: mynorm
      real(WP), dimension(3) :: csn,sn
      ! Reset curvature
      this%curv=0.0_WP
      this%curv2=0.0_WP
      ! Traverse interior domain and compute curvature in cells with polygons
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Zero out curvature and surface storage
               mycurv=0.0_WP; mysurf=0.0_WP; mynorm=0.0_WP
               ! Get a curvature for each plane
               do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
                  ! Skip empty polygon
                  if (getNumberOfVertices(this%interface_polygon(n,i,j,k)).eq.0) cycle
                  ! Perform LSQ PLIC barycenter fitting to get curvature
                  ! call this%paraboloid_fit(i,j,k,n,mycurv(n))
                  ! Perform PLIC surface fitting to get curvature
                  call this%paraboloid_integral_fit(i,j,k,n,mycurv(n))
                  ! Also store surface and normal
                  mysurf(n)  =abs(calculateVolume(this%interface_polygon(n,i,j,k)))
                  mynorm(n,:)=    calculateNormal(this%interface_polygon(n,i,j,k))
                  ! Separate curvatures
                  this%curv2(n,i,j,k)=mycurv(n)                  
               end do
               ! Oriented-surface-average curvature
               !csn=0.0_WP; sn=0.0_WP
               !do n=1,getNumberOfPlanes(this%liquid_gas_interface(i,j,k))
               !   csn=csn+mysurf(n)*mynorm(n,:)*mycurv(n)
               !   sn = sn+mysurf(n)*mynorm(n,:)
               !end do
               !if (dot_product(sn,sn).gt.10.0_WP*tiny(1.0_WP)) this%curv(i,j,k)=dot_product(csn,sn)/dot_product(sn,sn)
               ! Surface-averaged curvature
               if (sum(mysurf).gt.0.0_WP) this%curv(i,j,k)=sum(mysurf*mycurv)/sum(mysurf)
               ! Curvature of largest surface
               !if (mysurf(maxloc(mysurf,1)).gt.0.0_WP) this%curv(i,j,k)=mycurv(maxloc(mysurf,1))
               ! Largest curvature
               !this%curv(i,j,k)=mycurv(maxloc(abs(mycurv),1))
               ! Smallest curvature
               !if (getNumberOfPlanes(this%liquid_gas_interface(i,j,k)).eq.2) then
               !   this%curv(i,j,k)=mycurv(minloc(abs(mycurv),1))
               !else
               !   this%curv(i,j,k)=mycurv(1)
               !end if
               ! Clip curvature - may not be needed if we select polygons carefully
               this%curv(i,j,k)=max(min(this%curv(i,j,k),this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k)),-this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k))
               this%curv2(:,i,j,k)=max(min(this%curv2(:,i,j,k),this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k)),-this%maxcurv_times_mesh/this%cfg%meshsize(i,j,k))
            end do
         end do
      end do
      ! Synchronize boundaries
      call this%cfg%sync(this%curv)
      call this%cfg%sync(this%curv2)
   end subroutine get_curvature

end module film_vfs_class
