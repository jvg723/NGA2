!> Two-phase incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density in each phase.
!> Interface is represented using VOF
module film_tpns_class
   use precision,      only: WP
   use tpns_class,     only: tpns
   use string,         only: str_medium
   use config_class,   only: config
   use linsol_class,   only: linsol
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: film_tpns,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: wall=1              !< Dirichlet at zero condition
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   integer, parameter, public :: convective=4        !< Convective outflow condition
   integer, parameter, public :: clipped_neumann=5   !< Clipped Neumann condition (outflow only)
   
   ! List of available contact line models for this solver
   integer, parameter, public :: static_contact=1    !< Static contact line model
   
   ! Parameter for switching schemes around the interface
   real(WP), parameter :: rhoeps_coeff=1.0e-3_WP     !< Parameter for deciding when to switch to upwinded transport
   
   !> Boundary conditions for the two-phase solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      type(iterator) :: itr                               !< This is the iterator for the bcond - this identifies the (i,j,k)
      character(len=1) :: face                            !< Bcond face (x/y/z)
      integer :: dir                                      !< Bcond direction (+1,-1,0 for interior)
      real(WP) :: rdir                                    !< Bcond direction (real variable)
      logical :: canCorrect                               !< Can this bcond be corrected for global conservation?
   end type bcond
   
   !> Two-phase incompressible solver object definition
   type, extends(tpns) :: film_tpns
           
   contains
      procedure :: add_surface_tension_jump_film          !< Add surface tension jump
      procedure :: interp_pj                              !< Calculate interpolated pressure jump
   end type film_tpns
   
   
   !> Declare two-phase incompressible solver constructor
   interface film_tpns
      procedure constructor
   end interface film_tpns
   
contains
   
   
   !> Default constructor for two-phase incompressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(film_tpns) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate flow variables
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      allocate(self%Pjx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjx=0.0_WP
      allocate(self%Pjy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjy=0.0_WP
      allocate(self%Pjz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Pjz=0.0_WP
      allocate(self%dPjx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjx=0.0_WP
      allocate(self%dPjy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjy=0.0_WP
      allocate(self%dPjz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%dPjz=0.0_WP
      allocate(self%rho_U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_U=0.0_WP
      allocate(self%rho_V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_V=0.0_WP
      allocate(self%rho_W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_W=0.0_WP
      allocate(self%visc   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc   =0.0_WP
      allocate(self%visc_xy(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_xy=0.0_WP
      allocate(self%visc_yz(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_yz=0.0_WP
      allocate(self%visc_zx(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc_zx=0.0_WP
      
      ! Allocate flow divergence
      allocate(self%div(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%div=0.0_WP
      
      ! Allocate old flow variables
      allocate(self%Uold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Uold=0.0_WP
      allocate(self%Vold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Vold=0.0_WP
      allocate(self%Wold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Wold=0.0_WP
      allocate(self%rho_Uold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_Uold=0.0_WP
      allocate(self%rho_Vold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_Vold=0.0_WP
      allocate(self%rho_Wold(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%rho_Wold=0.0_WP
      
      ! Momentum fluxes need to be preallocated
      allocate(self%FUX(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FUX=0.0_WP
      allocate(self%FUY(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FUY=0.0_WP
      allocate(self%FUZ(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FUZ=0.0_WP
      allocate(self%FVX(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FVX=0.0_WP
      allocate(self%FVY(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FVY=0.0_WP
      allocate(self%FVZ(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FVZ=0.0_WP
      allocate(self%FWX(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FWX=0.0_WP
      allocate(self%FWY(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FWY=0.0_WP
      allocate(self%FWZ(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%FWZ=0.0_WP
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare P-cell masks
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
      
      ! Prepare face mask for U
      allocate(self%umask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%umask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%umask(self%cfg%imin  ,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%umask(self%cfg%imax+1,:,:)=2
      end if
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_+1,self%cfg%imaxo_
               if (minval(self%cfg%VF(i-1:i,j,k)).eq.0.0_WP) self%umask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%umask)
      if (.not.self%cfg%xper.and.self%cfg%iproc.eq.1) self%umask(self%cfg%imino,:,:)=self%umask(self%cfg%imino+1,:,:)
      
      ! Prepare face mask for V
      allocate(self%vmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%vmask=0
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%vmask(:,self%cfg%jmin  ,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%vmask(:,self%cfg%jmax+1,:)=2
      end if
      do k=self%cfg%kmino_  ,self%cfg%kmaxo_
         do j=self%cfg%jmino_+1,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j-1:j,k)).eq.0.0_WP) self%vmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%vmask)
      if (.not.self%cfg%yper.and.self%cfg%jproc.eq.1) self%vmask(:,self%cfg%jmino,:)=self%vmask(:,self%cfg%jmino+1,:)
      
      ! Prepare face mask for W
      allocate(self%wmask(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%wmask=0
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%wmask(:,:,self%cfg%kmin  )=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%wmask(:,:,self%cfg%kmax+1)=2
      end if
      do k=self%cfg%kmino_+1,self%cfg%kmaxo_
         do j=self%cfg%jmino_  ,self%cfg%jmaxo_
            do i=self%cfg%imino_  ,self%cfg%imaxo_
               if (minval(self%cfg%VF(i,j,k-1:k)).eq.0.0_WP) self%wmask(i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%wmask)
      if (.not.self%cfg%zper.and.self%cfg%kproc.eq.1) self%wmask(:,:,self%cfg%kmino)=self%wmask(:,:,self%cfg%kmino+1)
      
   end function constructor

   !> Add surface tension jump term using CSF
   !> GFM-style film treatment, precompute alpha
   subroutine add_surface_tension_jump_film(this,dt,div,vf,contact_model)
      use messager,  only: die
      use film_vfs_class, only: film_vfs, VFlo
      use ccl_class, only: ccl
      use mathtools, only: normalize
      use irl_fortran_interface
      implicit none
      class(film_tpns), intent(inout) :: this
      real(WP), intent(inout) :: dt     !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      class(film_vfs), intent(inout) :: vf
      integer, intent(in), optional :: contact_model
      integer :: i,j,k
      real(WP) :: mycurv,mysurf
      ! Alpha array creation
      type(ccl) :: cc
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_) :: alpha
      ! Curvature determination
      integer :: nplanem, nplanen      !< nplane my cell, nplane neighbor cell
      integer :: n,ind
      real(WP), dimension(3) :: n0,n1,n2
      real(WP), parameter :: edge_threshold=2.0_WP
      real(WP), parameter :: h_over_dx_thresh=0.5_WP
      logical, parameter :: use_surface_average=.true.
      ! Reconstruction
      type(LVIRANeigh_RectCub_type) :: neighborhood
      integer :: ii,jj,kk,icenter
      type(RectCub_type), dimension(0:26) :: neighborhood_cells
      real(IRL_double)  , dimension(0:26) :: liquid_volume_fraction
      real(IRL_double), dimension(3) :: initial_norm
      real(IRL_double) :: initial_dist
      ! type(RectCub_type) :: cell
      ! real(WP), dimension(1:4) :: plane

      ! Create the CCL object
      cc=ccl(cfg=this%cfg,name='CCL')
      cc%max_interface_planes=2
      cc%VFlo=VFlo
      cc%dot_threshold=-0.5_WP
      ! cc%dot_threshold=0.0_WP
      cc%thickness_cutoff=1.0_WP ! want to capture all films that require special ST treatment
      ! Perform CCL step
      call cc%build_lists(VF=vf%VF,poly=vf%interface_polygon,U=this%U,V=this%V,W=this%W)
      call cc%film_classify(Lbary=vf%Lbary,Gbary=vf%Gbary)
      call cc%deallocate_lists()
      
      ! ! Initialize IRL objects
      ! ! call new(cell)
      ! ! Give ourselves an LVIRA neighborhood of 27 cells
      ! call new(neighborhood)
      ! do i=0,26
      !    call new(neighborhood_cells(i))
      ! end do

      alpha=0.0_WP
      ! Precompute alpha
      do k=this%cfg%kmin_-1,this%cfg%kmax_+1
         do j=this%cfg%jmin_-1,this%cfg%jmax_+1
            do i=this%cfg%imin_-1,this%cfg%imax_+1
               ! if (cc%film_phase(i,j,k).eq.1) then
               !    alpha(i,j,k)=1.0_WP
               ! elseif (cc%film_phase(i,j,k).eq.2) then
               !    alpha(i,j,k)=0.0_WP
               ! else
               !    alpha(i,j,k)=vf%VF(i,j,k)
               ! endif
               if ((cc%film_phase(i,j,k).eq.1).and.(cc%film_type(i,j,k).eq.2)) then
                  alpha(i,j,k)=1.0_WP
               elseif ((cc%film_phase(i,j,k).eq.2).and.(cc%film_type(i,j,k).eq.2)) then
                  alpha(i,j,k)=0.0_WP
               else
                  alpha(i,j,k)=vf%VF(i,j,k)
               endif
               !! Correct edge cell alpha
               if (cc%film_type(i,j,k).eq.1) then
                  if (cc%film_edge(i,j,k).gt.3.0_WP) then
                     ! if (cc%film_thickness(i,j,k).lt.h_over_dx_thresh*this%cfg%meshsize(i,j,k)) then
                        vf%curv2(:,i,j,k)=0.5_WP/cc%film_thickness(i,j,k)
                        vf%curv2(:,i,j,k)=min(max(vf%curv2(:,i,j,k),0.0_WP),1.0_WP/(this%CFLst*vf%cfg%meshsize(i,j,k)))
                     !    ! vf%curv2(:,i,j,k)=1.0e4_WP
                     ! else
                     !    ! vf%curv2(:,i,j,k)=maxval(vf%curv2(:,i,j,k))
                     !    call vf%paraboloid_integral_fit(i,j,k,vf%curv2(1,i,j,k))
                     !    vf%curv2(:,i,j,k)=max(vf%curv2(1,i,j,k),0.0_WP)
                     ! end if
                  !    if (cc%film_thickness(i,j,k).gt.0.0_WP) vf%curv2(:,i,j,k)=0.1_WP/cc%film_thickness(i,j,k)
                     vf%curv(i,j,k)=vf%curv2(1,i,j,k)
                  elseif (cc%film_edge(i,j,k).gt.edge_threshold) then
                  ! if (cc%film_edge(i,j,k).gt.edge_threshold) then
                     !    ! 
                  !    ! ! alpha(i,j,k)=0.5_WP
                     ! if (cc%film_thickness(i,j,k).lt.h_over_dx_thresh*this%cfg%meshsize(i,j,k)) then
                     !    vf%curv2(:,i,j,k)=0.1_WP/cc%film_thickness(i,j,k)
                     !    vf%curv2(:,i,j,k)=min(max(vf%curv2(:,i,j,k),0.0_WP),1.0_WP/(this%CFLst*vf%cfg%meshsize(i,j,k)))
                     ! !    ! vf%curv2(:,i,j,k)=1.0e4_WP
                     ! else
                        ! alpha(i,j,k)=vf%VF(i,j,k)
                     !    vf%curv2(:,i,j,k)=maxval(vf%curv2(:,i,j,k))
                     ! end if
                     vf%curv2(:,i,j,k)=maxval(vf%curv2(:,i,j,k))
                     vf%curv2(:,i,j,k)=max(vf%curv2(:,i,j,k),0.0_WP)
                     ! vf%curv(i,j,k)=maxval(vf%curv2(:,i,j,k))
                     vf%curv(i,j,k)=vf%curv2(1,i,j,k)
                  elseif (cc%film_edge(i,j,k).gt.1.2_WP) then
                     vf%curv2(:,i,j,k)=max(vf%curv2(:,i,j,k),0.0_WP)
                     ! vf%curv(i,j,k)=max(vf%curv(i,j,k),0.0_WP)
                     vf%curv(i,j,k)=maxval(vf%curv2(:,i,j,k))
                  endif            
               endif      
            end do
         end do
      end do
      call this%cfg%sync(alpha)

      ! Store old jump
      this%DPjx=this%Pjx
      this%DPjy=this%Pjy
      this%DPjz=this%Pjz

      ! Zero pressure jump
      this%Pjx=0.0_WP
      this%Pjy=0.0_WP
      this%Pjz=0.0_WP
      ! Calculate pressure jump
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               nplanem=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j,k)).gt.0) then
                     nplanem=nplanem+1
                  end if
               end do            
               ! X face
               nplanen=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i-1,j,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i-1,j,k)).gt.0) then
                     nplanen=nplanen+1
                  end if
               end do
               mycurv=0.0_WP
               select case (nplanem+3*nplanen)
               case (0) ! non-R2P
                  mycurv=0.0_WP
               case (1) ! nplanem=1,nplanen=0
                  mycurv=vf%curv(i,j,k)
               case (3) ! 0+1
                  mycurv=vf%curv(i-1,j,k)
               case (4) ! 1+1
                  mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
                  mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
                  ! mycurv=vf%curv2(1,i,j,k)
               case (2) ! nplanem=2,nplanen=0
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  ind=maxloc([dot_product(n1,[-1.0_WP,0.0_WP,0.0_WP]),dot_product(n2,[-1.0_WP,0.0_WP,0.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i,j,k)
               case (5) ! 2+1
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
                     mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                     n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                     ind=maxloc([dot_product(n1,[-1.0_WP,0.0_WP,0.0_WP]),dot_product(n2,[-1.0_WP,0.0_WP,0.0_WP])],dim=1)
                     ! pair_curv(2)=vf%curv2(ind,i,j,k)
                     ! pair_curv(1)=vf%curv2(1,i-1,j,k)
                     mycurv=vf%curv2(ind,i,j,k) ! could surface average
                  end if
               case (6) ! 0+2
                  n1=calculateNormal(vf%interface_polygon(1,i-1,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i-1,j,k))
                  ind=maxloc([dot_product(n1,[+1.0_WP,0.0_WP,0.0_WP]),dot_product(n2,[+1.0_WP,0.0_WP,0.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i-1,j,k)
               case (7) ! 1+2
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
                     mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i-1,j,k))
                     n2=calculateNormal(vf%interface_polygon(2,i-1,j,k))
                     ind=maxloc([dot_product(n1,[+1.0_WP,0.0_WP,0.0_WP]),dot_product(n2,[+1.0_WP,0.0_WP,0.0_WP])],dim=1)
                     ! pair_curv(1)=vf%curv2(ind,i-1,j,k)
                     ! pair_curv(2)=vf%curv2(1,i,j,k)
                     mycurv=vf%curv2(ind,i-1,j,k)
                  end if
               case (8) ! 2+2
                  mysurf=sum(vf%SD(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))
                  mycurv=sum(vf%SD(i-1:i,j,k)*vf%curv(i-1:i,j,k)*this%cfg%vol(i-1:i,j,k))/mysurf
               end select
               this%Pjx(i,j,k)=this%sigma*mycurv*sum(this%divu_x(:,i,j,k)*alpha(i-1:i,j,k))
               ! Y face
               nplanen=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j-1,k))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j-1,k)).gt.0) then
                     nplanen=nplanen+1
                  end if
               end do
               mycurv=0.0_WP
               select case (nplanem+3*nplanen)
               case (0) ! non-R2P
                  mycurv=0.0_WP
               case (1) ! nplanem=1,nplanen=0
                  mycurv=vf%curv(i,j,k)
               case (3) ! 0+1
                  mycurv=vf%curv(i,j-1,k)
               case (4) ! 1+1         
                  mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
                  mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
                  ! mycurv=vf%curv2(1,i,j,k)
               case (2) ! nplanem=2,nplanen=0
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  ind=maxloc([dot_product(n1,[0.0_WP,-1.0_WP,0.0_WP]),dot_product(n2,[0.0_WP,-1.0_WP,0.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i,j,k)
               case (5) ! 2+1
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
                     mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                     n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                     ind=maxloc([dot_product(n1,[0.0_WP,-1.0_WP,0.0_WP]),dot_product(n2,[0.0_WP,-1.0_WP,0.0_WP])],dim=1)
                     mycurv=vf%curv2(ind,i,j,k)
                  end if
               case (6) ! 0+2
                  n1=calculateNormal(vf%interface_polygon(1,i,j-1,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j-1,k))
                  ind=maxloc([dot_product(n1,[0.0_WP,+1.0_WP,0.0_WP]),dot_product(n2,[0.0_WP,+1.0_WP,0.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i,j-1,k)
               case (7) ! 1+2
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
                     mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i,j-1,k))
                     n2=calculateNormal(vf%interface_polygon(2,i,j-1,k))
                     ind=maxloc([dot_product(n1,[0.0_WP,+1.0_WP,0.0_WP]),dot_product(n2,[0.0_WP,+1.0_WP,0.0_WP])],dim=1)
                     mycurv=vf%curv2(ind,i,j-1,k)
                  end if         
               case (8) ! 2+2
                  mysurf=sum(vf%SD(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))
                  mycurv=sum(vf%SD(i,j-1:j,k)*vf%curv(i,j-1:j,k)*this%cfg%vol(i,j-1:j,k))/mysurf
               end select
               this%Pjy(i,j,k)=this%sigma*mycurv*sum(this%divv_y(:,i,j,k)*alpha(i,j-1:j,k))
               ! Z face
               nplanen=0
               do n=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k-1))
                  if (getNumberOfVertices(vf%interface_polygon(n,i,j,k-1)).gt.0) then
                     nplanen=nplanen+1
                  end if
               end do
               mycurv=0.0_WP
               select case (nplanem+3*nplanen)
               case (0) ! non-R2P
                  mycurv=0.0_WP
               case (1) ! nplanem=1,nplanen=0
                  mycurv=vf%curv(i,j,k)
               case (3) ! 0+1
                  mycurv=vf%curv(i,j,k-1)
               case (4) ! 1+1         
                  mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
                  mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
                  ! mycurv=vf%curv2(1,i,j,k)
               case (2) ! nplanem=2,nplanen=0
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                  ind=maxloc([dot_product(n1,[0.0_WP,0.0_WP,-1.0_WP]),dot_product(n2,[0.0_WP,0.0_WP,-1.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i,j,k)
               case (5) ! 2+1
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
                     mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                     n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                     ind=maxloc([dot_product(n1,[0.0_WP,0.0_WP,-1.0_WP]),dot_product(n2,[0.0_WP,0.0_WP,-1.0_WP])],dim=1)
                     mycurv=vf%curv2(ind,i,j,k)
                  end if
               case (6) ! 0+2
                  n1=calculateNormal(vf%interface_polygon(1,i,j,k-1))
                  n2=calculateNormal(vf%interface_polygon(2,i,j,k-1))
                  ind=maxloc([dot_product(n1,[0.0_WP,0.0_WP,+1.0_WP]),dot_product(n2,[0.0_WP,0.0_WP,+1.0_WP])],dim=1)
                  mycurv=vf%curv2(ind,i,j,k-1)
               case (7) ! 1+2
                  if (use_surface_average) then
                     mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
                     mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
                  else
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k-1))
                     n2=calculateNormal(vf%interface_polygon(2,i,j,k-1))
                     ind=maxloc([dot_product(n1,[0.0_WP,0.0_WP,+1.0_WP]),dot_product(n2,[0.0_WP,0.0_WP,+1.0_WP])],dim=1)
                     mycurv=vf%curv2(ind,i,j,k-1)
                  end if
               case (8) ! 2+2
                  mysurf=sum(vf%SD(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))
                  mycurv=sum(vf%SD(i,j,k-1:k)*vf%curv(i,j,k-1:k)*this%cfg%vol(i,j,k-1:k))/mysurf
               end select
               this%Pjz(i,j,k)=this%sigma*mycurv*sum(this%divw_z(:,i,j,k)*alpha(i,j,k-1:k))
            end do
         end do
      end do
      
      ! Add wall contact force to pressure jump
      if (present(contact_model)) then
         select case (contact_model)
         case (static_contact)
            call this%add_static_contact(vf=vf)
         case default
            call die('[tpns: add_surface_tension_jump_film] Unknown contact model!')
         end select
      end if
      
      ! Compute jump of DP
      this%DPjx=this%Pjx-this%DPjx
      this%DPjy=this%Pjy-this%DPjy
      this%DPjz=this%Pjz-this%DPjz
      
      ! Add div(Pjump) to RP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               div(i,j,k)=div(i,j,k)+dt*(sum(this%divp_x(:,i,j,k)*this%DPjx(i:i+1,j,k)/this%rho_U(i:i+1,j,k))&
               &                        +sum(this%divp_y(:,i,j,k)*this%DPjy(i,j:j+1,k)/this%rho_V(i,j:j+1,k))&
               &                        +sum(this%divp_z(:,i,j,k)*this%DPjz(i,j,k:k+1)/this%rho_W(i,j,k:k+1)))
            end do
         end do
      end do 

   end subroutine add_surface_tension_jump_film

   
   !> Calculate the interpolated pressure jump, including overlap and ghosts
   subroutine interp_pj(this,Pjxi,Pjyi,Pjzi)
      implicit none
      class(film_tpns), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pjxi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pjyi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: Pjzi !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      ! Calculate as far as possible each component
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_-1
               Pjxi(i,j,k)=sum(this%itpu_x(:,i,j,k)*this%Pjx(i:i+1,j,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_-1
            do i=this%cfg%imino_,this%cfg%imaxo_
               Pjyi(i,j,k)=sum(this%itpv_y(:,i,j,k)*this%Pjy(i,j:j+1,k))
            end do
         end do
      end do
      do k=this%cfg%kmino_,this%cfg%kmaxo_-1
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               Pjzi(i,j,k)=sum(this%itpw_z(:,i,j,k)*this%Pjz(i,j,k:k+1))
            end do
         end do
      end do
      ! Add last layer in each direction
      if (.not.this%cfg%xper.and.this%cfg%iproc.eq.this%cfg%npx) Pjxi(this%cfg%imaxo,:,:)=this%Pjx(this%cfg%imaxo,:,:)
      if (.not.this%cfg%yper.and.this%cfg%jproc.eq.this%cfg%npy) Pjyi(:,this%cfg%jmaxo,:)=this%Pjy(:,this%cfg%jmaxo,:)
      if (.not.this%cfg%zper.and.this%cfg%kproc.eq.this%cfg%npz) Pjzi(:,:,this%cfg%kmaxo)=this%Pjz(:,:,this%cfg%kmaxo)
      ! Sync it
      call this%cfg%sync(Pjxi)
      call this%cfg%sync(Pjyi)
      call this%cfg%sync(Pjzi)
   end subroutine interp_pj

end module film_tpns_class
