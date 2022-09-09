!> FENE model solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes constant diffusivity and density.
module fene_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: fene,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available advection schemes for fene transport
   integer, parameter, public :: quick =1             !< Quick scheme
   integer, parameter, public :: bquick=2             !< Bquick scheme
   integer, parameter, public :: upwind=3             !< Upwind scheme
   
   
   !> Boundary conditions for the incompressible solver
   type :: bcond
      type(bcond), pointer :: next                        !< Linked list of bconds
      character(len=str_medium) :: name='UNNAMED_BCOND'   !< Bcond name (default=UNNAMED_BCOND)
      integer :: type                                     !< Bcond type
      integer :: dir                                      !< Bcond direction (1 to 6)
      type(iterator) :: itr                               !< This is the iterator for the bcond
   end type bcond
   
   !> Bcond shift value
   integer, dimension(3,6), parameter :: shift=reshape([+1,0,0,-1,0,0,0,+1,0,0,-1,0,0,0,+1,0,0,-1],shape(shift))
   
   !> Constant density fene solver object definition
   type :: fene
      
      ! This is our config
      class(config), pointer :: cfg                         !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_FENE'      !< Solver name (default=UNNAMED_FENE)
      
      ! Constant property fluid
      real(WP) :: rho                                       !< This is our constant fluid density
      
      ! Boundary condition list
      integer :: nbc                                        !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                      !< List of bcond for our solver
      
      ! fene conformation tensor variable
      real(WP), dimension(:,:,:,:), allocatable :: CT       !< CT array
      
      ! Old fene conformation tensor variable
      real(WP), dimension(:,:,:,:), allocatable :: CTold    !< CTold array

      ! Conformation tensor used in 1st order upwind scheme
      real(WP), dimension(:,:,:,:), allocatable :: CTupwind !< CT array used with upwind

      ! Source term arrays
      real(WP), dimension(:,:,:,:), allocatable :: CTgradU  !< Product of CT and velocity gradient
      real(WP), dimension(:,:,:,:), allocatable :: CTgradUT !< Transpose of CT dot gradU
      real(WP), dimension(:,:,:,:), allocatable :: T        !< Viscoelastic stress tensor

      ! Viscoleastic tensor divergence
      real(WP), dimension(:,:,:,:), allocatable :: divT     !< Divergence array

      ! Number of elements in conformation tensor
      integer :: CTelem=9                                   !< Symmetric tensor with 6 unique elements
      
      ! Implicit fene solver
      type(ils) :: implicit                                 !< Iterative linear solver object for an implicit prediction of the fene residual
      integer, dimension(:,:,:), allocatable :: stmap       !< Inverse map from stencil shift to index location
      
      ! Metrics
      integer :: scheme                                                    !< Advection scheme for fene model
      integer :: nst,qnst                                                  !< Scheme order (and elemental stencil size)
      integer :: stp1,stp2,qstp1,qstp2                                     !< Plus interpolation stencil extent for fene model advection
      integer :: stm1,stm2,qstm1,qstm2                                     !< Minus interpolation stencil extent for fene model advection
      real(WP), dimension(:,:,:,:), allocatable :: itp_xp,itp_yp,itp_zp    !< Plus interpolation coefficents
      real(WP), dimension(:,:,:,:), allocatable :: itp_xm,itp_ym,itp_zm    !< Minus interpolation coefficents
      real(WP), dimension(:,:,:,:), allocatable :: qitp_xp,qitp_yp,qitp_zp !< Stored Plus interpolation coefficents for QUICK
      real(WP), dimension(:,:,:,:), allocatable :: qitp_xm,qitp_ym,qitp_zm !< Stored Minus interpolation coefficents for QUICK
      real(WP), dimension(:,:,:,:), allocatable :: div_x ,div_y ,div_z     !< Divergence for CT
      real(WP), dimension(:,:,:,:), allocatable :: grd_x ,grd_y ,grd_z     !< Velocity gradient  
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for modifying CT metrics
      
      
   contains
      procedure :: print=>fene_print                      !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the fene solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: upwind_metrics                         !< Calculate upwind metrics for a given i,j,k
      procedure :: get_CTgradU                            !< Calculate product and transpose of CT_dot_gradU
      procedure :: get_stressTensor                       !< Calculate the viscoelastic stress tensor
      procedure :: get_drhoCTdt                           !< Calculate drhoCT/dt
      procedure :: get_divT                               !< Calculate viscoelastic tensor divergence
      procedure :: solve_implicit                         !< Solve for the fene residuals implicitly
   end type fene
   
   
   !> Declare fene solver constructor
   interface fene
      procedure constructor
   end interface fene
   
contains
   
   
   !> Default constructor for fene solver
   function constructor(cfg,scheme,name) result(self)
      use messager, only: die
      implicit none
      type(fene) :: self
      class(config), target, intent(in) :: cfg
      integer,  intent(in) :: scheme
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Nullify bcond list
      self%nbc=0
      self%first_bc=>NULL()
      
      ! Allocate variables
      allocate(self%CT      (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%CT      =0.0_WP
      allocate(self%CTold   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%CTold   =0.0_WP
      allocate(self%CTgradU (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%CTgradU =0.0_WP
      allocate(self%CTgradUT(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%CTgradUT=0.0_WP
      allocate(self%T       (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%T       =0.0_WP
      allocate(self%divT    (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3          )); self%divT    =0.0_WP

      ! Prepare advection scheme
      self%scheme=scheme
      select case (self%scheme)
      case (quick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[fene constructor] fene transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case (bquick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[fene constructor] fene transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
         ! Allocate array for computing CT with upwind scheme
         allocate(self%CTupwind(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%CTelem)); self%CTgradUT=0.0_WP
      case default
         call die('[fene constructor] Unknown fene transport scheme selected')
      end select
      
      ! Create implicit fene solver object
      self%implicit=ils(cfg=self%cfg,name='fene',nst=1+6*abs(self%stp1))
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for CT
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
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(fene), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference fene interpolation coefficients
      allocate(this%itp_xp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_xm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_yp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_ym(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_zp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%itp_zm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create fene interpolation coefficients to cell faces
      select case (this%scheme)
      case (quick)
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itp_xp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itp_xm(:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itp_yp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itp_ym(:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itp_zp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itp_zm(:,i,j,k))
               end do
            end do
         end do
      case (bquick)
         allocate(this%qitp_xp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
         allocate(this%qitp_xm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
         allocate(this%qitp_yp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
         allocate(this%qitp_ym(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
         allocate(this%qitp_zp(this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
         allocate(this%qitp_zm(this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itp_xp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itp_xm(:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itp_yp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itp_ym(:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itp_zp(:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itp_zm(:,i,j,k))
               end do
            end do
         end do
         ! Store QUICK stencil for reset
         this%qnst=this%nst
         this%stp1=this%qstp1
         this%stp2=this%qstp2
         this%stm1=this%qstm1
         this%stm2=this%qstm2
         ! Store QUICK parameters for reset
         this%qitp_xp=this%itp_xp
         this%qitp_xm=this%itp_xm
         this%qitp_yp=this%itp_yp
         this%qitp_ym=this%itp_ym
         this%qitp_zp=this%itp_zp
         this%qitp_zm=this%itp_zm
      end select
      
      ! Allocate finite volume divergence operators
      allocate(this%div_x(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_y(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%div_z(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%div_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%div_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%div_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do

      ! Allocate finite difference velocity gradient operators
      allocate(this%grd_x(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%grd_y(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%grd_z(0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create gradient coefficients to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%grd_x(:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FD gradient in x from [x ,ym,zm]
               this%grd_y(:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FD gradient in y from [xm,y ,zm]
               this%grd_z(:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FD gradient in z from [xm,ym,z ]
            end do
         end do
      end do

   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(fene), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask)
      
      ! Adjust fene interpolation to reflect Dirichlet boundaries
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face
               if (this%mask(i-1,j,k).eq.2) then
                  this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm(-1,i,j,k)=1.0_WP
                  this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i  ,j,k).eq.2) then
                  this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm( 0,i,j,k)=1.0_WP
                  this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp( 0,i,j,k)=1.0_WP
               end if
               ! Y face
               if (this%mask(i,j-1,k).eq.2) then
                  this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym(-1,i,j,k)=1.0_WP
                  this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j  ,k).eq.2) then
                  this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym( 0,i,j,k)=1.0_WP
                  this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp( 0,i,j,k)=1.0_WP
               end if
               ! Z face
               if (this%mask(i,j,k-1).eq.2) then
                  this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm(-1,i,j,k)=1.0_WP
                  this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp(-1,i,j,k)=1.0_WP
               end if
               if (this%mask(i,j,k  ).eq.2) then
                  this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm( 0,i,j,k)=1.0_WP
                  this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp( 0,i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to CT divergence
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).gt.0) then
                  this%div_x(:,i,j,k)=0.0_WP
                  this%div_y(:,i,j,k)=0.0_WP
                  this%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(i,j,k).eq.1.or.this%mask(i-1,j,k).eq.1) this%grd_x(:,i,j,k)=0.0_WP     !< FD gradient in x of CT
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j-1,k).eq.1) this%grd_y(:,i,j,k)=0.0_WP     !< FD gradient in y of CT
               if (this%mask(i,j,k).eq.1.or.this%mask(i,j,k-1).eq.1) this%grd_z(:,i,j,k)=0.0_WP     !< FD gradient in z of CT
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%div_x=0.0_WP
         this%grd_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%div_y=0.0_WP
         this%grd_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%div_z=0.0_WP
         this%grd_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the fene solver now that bconds have been defined
   subroutine setup(this,implicit_ils)
      implicit none
      class(fene), intent(inout) :: this
      integer, intent(in) :: implicit_ils
      integer :: count,st
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
      ! Set dynamic stencil map for the fene solver
      count=1; this%implicit%stc(count,:)=[0,0,0]
      do st=1,abs(this%stp1)
         count=count+1; this%implicit%stc(count,:)=[+st,0,0]
         count=count+1; this%implicit%stc(count,:)=[-st,0,0]
         count=count+1; this%implicit%stc(count,:)=[0,+st,0]
         count=count+1; this%implicit%stc(count,:)=[0,-st,0]
         count=count+1; this%implicit%stc(count,:)=[0,0,+st]
         count=count+1; this%implicit%stc(count,:)=[0,0,-st]
      end do
      
      ! Set the diagonal to 1 to make sure all cells participate in solver
      this%implicit%opr(1,:,:,:)=1.0_WP
      
      ! Initialize the implicit fene solver
      call this%implicit%init(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(fene), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer,  intent(in) :: type
      external :: locator
      interface
         logical function locator(pargrid,ind1,ind2,ind3)
            use pgrid_class, only: pgrid
            class(pgrid), intent(in) :: pargrid
            integer, intent(in) :: ind1,ind2,ind3
         end function locator
      end interface
      character(len=2), optional :: dir
      type(bcond), pointer :: new_bc
      integer :: i,j,k,n
      
      ! Prepare new bcond
      allocate(new_bc)
      new_bc%name=trim(adjustl(name))
      new_bc%type=type
      if (present(dir)) then
         select case (lowercase(dir))
         case ('+x','x+','xp','px'); new_bc%dir=1
         case ('-x','x-','xm','mx'); new_bc%dir=2
         case ('+y','y+','yp','py'); new_bc%dir=3
         case ('-y','y-','ym','my'); new_bc%dir=4
         case ('+z','z+','zp','pz'); new_bc%dir=5
         case ('-z','z-','zm','mz'); new_bc%dir=6
         case default; call die('[fene add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[fene apply_bcond] Neumann requires a direction')
         new_bc%dir=0
      end if
      new_bc%itr=iterator(this%cfg,new_bc%name,locator)
      
      ! Insert it up front
      new_bc%next=>this%first_bc
      this%first_bc=>new_bc
      
      ! Increment bcond counter
      this%nbc=this%nbc+1
      
      ! Now adjust the metrics accordingly
      select case (new_bc%type)
      case (dirichlet)
         do n=1,new_bc%itr%n_
            i=new_bc%itr%map(1,n); j=new_bc%itr%map(2,n); k=new_bc%itr%map(3,n)
            this%mask(i,j,k)=2
         end do
      case (neumann)
         ! No modification - this assumes Neumann is only applied at walls or domain boundaries
      case default
         call die('[fene apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(fene), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[fene get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(fene), intent(inout) :: this
      real(WP), intent(in) :: t,dt
      integer :: i,j,k,n
      type(bcond), pointer :: my_bc
      
      ! Traverse bcond list
      my_bc=>this%first_bc
      do while (associated(my_bc))
         
         ! Only processes inside the bcond work here
         if (my_bc%itr%amIn) then
            
            ! Select appropriate action based on the bcond type
            select case (my_bc%type)
               
            case (dirichlet)           ! Apply Dirichlet conditions
               
               ! This is done by the user directly
               ! Unclear whether we want to do this within the solver...
               
            case (neumann)             ! Apply Neumann condition
               
               ! Implement based on bcond direction
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  this%CT(:,i,j,k)=this%CT(:,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[fene apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%CT)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond

   !> Calculate interpolation metrics for 1st order upwind scheme
   subroutine upwind_metrics(this,i,j,k)
      use mathtools, only: fv_itp_build
      implicit none
      class(fene), intent(inout) :: this
      integer, intent(in) :: i,j,k
      
      ! Set interpolation stencil sizes for 1st order upwind
      this%nst=1
      this%stp1=-(this%nst+1)/2; this%stp2=this%nst+this%stp1-1
      this%stm1=-(this%nst-1)/2; this%stm2=this%nst+this%stm1-1

      ! Interpolation to x-face
      call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itp_xp(:,i,j,k))
      call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itp_xm(:,i,j,k))
      ! Interpolation to y-face
      call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itp_yp(:,i,j,k))
      call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itp_ym(:,i,j,k))
      ! Interpolation to z-face
      call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itp_zp(:,i,j,k))
      call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itp_zm(:,i,j,k))

      ! Sync up masks
      call this%cfg%sync(this%mask)

      ! X face
      if (this%mask(i-1,j,k).eq.2) then
         this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm(-1,i,j,k)=1.0_WP
         this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp(-1,i,j,k)=1.0_WP
      end if
      if (this%mask(i  ,j,k).eq.2) then
         this%itp_xm(:,i,j,k)=0.0_WP; this%itp_xm( 0,i,j,k)=1.0_WP
         this%itp_xp(:,i,j,k)=0.0_WP; this%itp_xp( 0,i,j,k)=1.0_WP
      end if
      ! Y face
      if (this%mask(i,j-1,k).eq.2) then
         this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym(-1,i,j,k)=1.0_WP
         this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp(-1,i,j,k)=1.0_WP
      end if
      if (this%mask(i,j  ,k).eq.2) then
         this%itp_ym(:,i,j,k)=0.0_WP; this%itp_ym( 0,i,j,k)=1.0_WP
         this%itp_yp(:,i,j,k)=0.0_WP; this%itp_yp( 0,i,j,k)=1.0_WP
      end if
      ! Z face
      if (this%mask(i,j,k-1).eq.2) then
         this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm(-1,i,j,k)=1.0_WP
         this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp(-1,i,j,k)=1.0_WP
      end if
      if (this%mask(i,j,k  ).eq.2) then
         this%itp_zm(:,i,j,k)=0.0_WP; this%itp_zm( 0,i,j,k)=1.0_WP
         this%itp_zp(:,i,j,k)=0.0_WP; this%itp_zp( 0,i,j,k)=1.0_WP
      end if

   end subroutine upwind_metrics

   !> Calculate prduct of CT dot gradU and its transpose
   subroutine get_CTgradU(this,scheme,gradU,ii,jj,kk)
      implicit none
      class(fene), intent(inout) :: this
      integer,  intent(in) :: scheme
      integer,  intent(in) :: ii,jj,kk
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: gradU
      integer :: i,j,k
      select case (scheme)
      case(quick)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_          
                  ! Product of C dot gradU
                  this%CTgradU(i,j,k,1)=this%CT(i,j,k,1)*gradU(i,j,k,1)+this%CT(i,j,k,2)*gradU(i,j,k,4)+this%CT(i,j,k,3)*gradU(i,j,k,7); this%CTgradU(i,j,k,2)=this%CT(i,j,k,1)*gradU(i,j,k,2)+this%CT(i,j,k,2)*gradU(i,j,k,5)+this%CT(i,j,k,3)*gradU(i,j,k,8); this%CTgradU(i,j,k,3)=this%CT(i,j,k,1)*gradU(i,j,k,3)+this%CT(i,j,k,2)*gradU(i,j,k,6)+this%CT(i,j,k,3)*gradU(i,j,k,9)
                  this%CTgradU(i,j,k,4)=this%CT(i,j,k,4)*gradU(i,j,k,1)+this%CT(i,j,k,5)*gradU(i,j,k,4)+this%CT(i,j,k,6)*gradU(i,j,k,7); this%CTgradU(i,j,k,5)=this%CT(i,j,k,4)*gradU(i,j,k,2)+this%CT(i,j,k,5)*gradU(i,j,k,5)+this%CT(i,j,k,6)*gradU(i,j,k,8); this%CTgradU(i,j,k,6)=this%CT(i,j,k,4)*gradU(i,j,k,3)+this%CT(i,j,k,5)*gradU(i,j,k,6)+this%CT(i,j,k,6)*gradU(i,j,k,9)
                  this%CTgradU(i,j,k,7)=this%CT(i,j,k,7)*gradU(i,j,k,1)+this%CT(i,j,k,8)*gradU(i,j,k,4)+this%CT(i,j,k,9)*gradU(i,j,k,7); this%CTgradU(i,j,k,8)=this%CT(i,j,k,7)*gradU(i,j,k,2)+this%CT(i,j,k,8)*gradU(i,j,k,5)+this%CT(i,j,k,9)*gradU(i,j,k,8); this%CTgradU(i,j,k,9)=this%CT(i,j,k,7)*gradU(i,j,k,3)+this%CT(i,j,k,8)*gradU(i,j,k,6)+this%CT(i,j,k,9)*gradU(i,j,k,9)
                  ! Transpose it
                  this%CTgradUT(i,j,k,1)=this%CTgradU(i,j,k,1); this%CTgradUT(i,j,k,2)=this%CTgradU(i,j,k,4); this%CTgradUT(i,j,k,3)=this%CTgradU(i,j,k,7)
                  this%CTgradUT(i,j,k,4)=this%CTgradU(i,j,k,2); this%CTgradUT(i,j,k,5)=this%CTgradU(i,j,k,5); this%CTgradUT(i,j,k,6)=this%CTgradU(i,j,k,8)
                  this%CTgradUT(i,j,k,7)=this%CTgradU(i,j,k,3); this%CTgradUT(i,j,k,8)=this%CTgradU(i,j,k,6); this%CTgradUT(i,j,k,9)=this%CTgradU(i,j,k,9)
               end do
            end do
         end do
      case(upwind)
         ! Product of C dot gradU
         this%CTgradU(ii,jj,kk,1)=this%CT(ii,jj,kk,1)*gradU(ii,jj,kk,1)+this%CT(ii,jj,kk,2)*gradU(ii,jj,kk,4)+this%CT(ii,jj,kk,3)*gradU(ii,jj,kk,7); this%CTgradU(ii,jj,kk,2)=this%CT(ii,jj,kk,1)*gradU(ii,jj,kk,2)+this%CT(ii,jj,kk,2)*gradU(ii,jj,kk,5)+this%CT(ii,jj,kk,3)*gradU(ii,jj,kk,8); this%CTgradU(ii,jj,kk,3)=this%CT(ii,jj,kk,1)*gradU(ii,jj,kk,3)+this%CT(ii,jj,kk,2)*gradU(ii,jj,kk,6)+this%CT(ii,jj,kk,3)*gradU(ii,jj,kk,9)
         this%CTgradU(ii,jj,kk,4)=this%CT(ii,jj,kk,4)*gradU(ii,jj,kk,1)+this%CT(ii,jj,kk,5)*gradU(ii,jj,kk,4)+this%CT(ii,jj,kk,6)*gradU(ii,jj,kk,7); this%CTgradU(ii,jj,kk,5)=this%CT(ii,jj,kk,4)*gradU(ii,jj,kk,2)+this%CT(ii,jj,kk,5)*gradU(ii,jj,kk,5)+this%CT(ii,jj,kk,6)*gradU(ii,jj,kk,8); this%CTgradU(ii,jj,kk,6)=this%CT(ii,jj,kk,4)*gradU(ii,jj,kk,3)+this%CT(ii,jj,kk,5)*gradU(ii,jj,kk,6)+this%CT(ii,jj,kk,6)*gradU(ii,jj,kk,9)
         this%CTgradU(ii,jj,kk,7)=this%CT(ii,jj,kk,7)*gradU(ii,jj,kk,1)+this%CT(ii,jj,kk,8)*gradU(ii,jj,kk,4)+this%CT(ii,jj,kk,9)*gradU(ii,jj,kk,7); this%CTgradU(ii,jj,kk,8)=this%CT(ii,jj,kk,7)*gradU(ii,jj,kk,2)+this%CT(ii,jj,kk,8)*gradU(ii,jj,kk,5)+this%CT(ii,jj,kk,9)*gradU(ii,jj,kk,8); this%CTgradU(ii,jj,kk,9)=this%CT(ii,jj,kk,7)*gradU(ii,jj,kk,3)+this%CT(ii,jj,kk,8)*gradU(ii,jj,kk,6)+this%CT(ii,jj,kk,9)*gradU(ii,jj,kk,9)
         ! Transpose it
         this%CTgradUT(ii,jj,kk,1)=this%CTgradU(ii,jj,kk,1); this%CTgradUT(ii,jj,kk,2)=this%CTgradU(ii,jj,kk,4); this%CTgradUT(ii,jj,kk,3)=this%CTgradU(ii,jj,kk,7)
         this%CTgradUT(ii,jj,kk,4)=this%CTgradU(ii,jj,kk,2); this%CTgradUT(ii,jj,kk,5)=this%CTgradU(ii,jj,kk,5); this%CTgradUT(ii,jj,kk,6)=this%CTgradU(ii,jj,kk,8)
         this%CTgradUT(ii,jj,kk,7)=this%CTgradU(ii,jj,kk,3); this%CTgradUT(ii,jj,kk,8)=this%CTgradU(ii,jj,kk,6); this%CTgradUT(ii,jj,kk,9)=this%CTgradU(ii,jj,kk,9)
      end select
      ! Sync it
      call this%cfg%sync(this%CTgradU)
      call this%cfg%sync(this%CTgradUT)
   end subroutine get_CTgradU

   !> Calculate the viscoelastic stress tensor
   subroutine get_stressTensor(this,scheme,We,Lmax,ii,jj,kk)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), intent(in) :: Lmax
      real(WP), intent(in) :: We
      integer,  intent(in) :: scheme
      integer,  intent(in) :: ii,jj,kk
      real(WP) :: a,psi
      integer :: i,j,k
      ! a parameter
      a=1.00_WP-3.00_WP/(Lmax**2)
      select case (scheme)
      case(quick)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! psi parameter
               psi=1.00_WP-(this%CT(i,j,k,1)+this%CT(i,j,k,5)+this%CT(i,j,k,9))/(Lmax**2)   
               ! Build visoelastic stress tensor
               this%T(i,j,k,1)=this%CT(i,j,k,1)/psi-1.00_WP/a; this%T(i,j,k,2)=this%CT(i,j,k,2)/psi-0.00_WP/a; this%T(i,j,k,3)=this%CT(i,j,k,3)/psi-0.00_WP/a;          
               this%T(i,j,k,4)=this%CT(i,j,k,4)/psi-0.00_WP/a; this%T(i,j,k,5)=this%CT(i,j,k,5)/psi-1.00_WP/a; this%T(i,j,k,6)=this%CT(i,j,k,6)/psi-0.00_WP/a;
               this%T(i,j,k,7)=this%CT(i,j,k,7)/psi-0.00_WP/a; this%T(i,j,k,8)=this%CT(i,j,k,8)/psi-0.00_WP/a; this%T(i,j,k,9)=this%CT(i,j,k,9)/psi-1.00_WP/a;
            end do
         end do
      end do
      case(upwind)
         ! psi parameter
         psi=1.00_WP-(this%CT(ii,jj,kk,1)+this%CT(ii,jj,kk,5)+this%CT(ii,jj,kk,9))/(Lmax**2)   
         ! Build visoelastic stress tensor
         this%T(ii,jj,kk,1)=this%CT(ii,jj,kk,1)/psi-1.00_WP/a; this%T(ii,jj,kk,2)=this%CT(ii,jj,kk,2)/psi-0.00_WP/a; this%T(ii,jj,kk,3)=this%CT(ii,jj,kk,3)/psi-0.00_WP/a;          
         this%T(ii,jj,kk,4)=this%CT(ii,jj,kk,4)/psi-0.00_WP/a; this%T(ii,jj,kk,5)=this%CT(ii,jj,kk,5)/psi-1.00_WP/a; this%T(ii,jj,kk,6)=this%CT(ii,jj,kk,6)/psi-0.00_WP/a;
         this%T(ii,jj,kk,7)=this%CT(ii,jj,kk,7)/psi-0.00_WP/a; this%T(ii,jj,kk,8)=this%CT(ii,jj,kk,8)/psi-0.00_WP/a; this%T(ii,jj,kk,9)=this%CT(ii,jj,kk,9)/psi-1.00_WP/a;
      end select
      ! Sync it
      call this%cfg%sync(this%T)
   end subroutine get_stressTensor
   
   !> Calculate the explicit rhoCT time derivative based on rhoU/rhoV/rhoW
   subroutine get_drhoCTdt(this,scheme,drhoCTdt,rhoU,rhoV,rhoW,ii,jj,kk)
      implicit none
      class(fene), intent(inout) :: this
      integer,  intent(in) :: scheme
      integer,  intent(in) :: ii,jj,kk
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(out) :: drhoCTdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)  :: rhoU     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)  :: rhoV     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)  :: rhoW     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,n
      real(WP), dimension(:,:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,this%CTelem))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,this%CTelem))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,this%CTelem))
      select case (scheme)
      case(quick)
      print *,'In QUICK, pre flux'
         ! Flux of rhoCT
         do n=1,this%CTelem
            do k=this%cfg%kmin_,this%cfg%kmax_+1
               do j=this%cfg%jmin_,this%cfg%jmax_+1
                  do i=this%cfg%imin_,this%cfg%imax_+1
                     ! Fluxes on x-face
                     FX(i,j,k,n)=-0.5_WP*(rhoU(i,j,k)+abs(rhoU(i,j,k)))*sum(this%itp_xp(:,i,j,k)*this%CT(i+this%stp1:i+this%stp2,j,k,n)) &
                     &           -0.5_WP*(rhoU(i,j,k)-abs(rhoU(i,j,k)))*sum(this%itp_xm(:,i,j,k)*this%CT(i+this%stm1:i+this%stm2,j,k,n))
                     ! Fluxes on y-face
                     FY(i,j,k,n)=-0.5_WP*(rhoV(i,j,k)+abs(rhoV(i,j,k)))*sum(this%itp_yp(:,i,j,k)*this%CT(i,j+this%stp1:j+this%stp2,k,n)) &
                     &           -0.5_WP*(rhoV(i,j,k)-abs(rhoV(i,j,k)))*sum(this%itp_ym(:,i,j,k)*this%CT(i,j+this%stm1:j+this%stm2,k,n)) 
                     ! Fluxes on z-face
                     FZ(i,j,k,n)=-0.5_WP*(rhoW(i,j,k)+abs(rhoW(i,j,k)))*sum(this%itp_zp(:,i,j,k)*this%CT(i,j,k+this%stp1:k+this%stp2,n)) &
                     &           -0.5_WP*(rhoW(i,j,k)-abs(rhoW(i,j,k)))*sum(this%itp_zm(:,i,j,k)*this%CT(i,j,k+this%stm1:k+this%stm2,n)) 
                  end do
               end do
            end do
         end do
         print *, 'pre time derivative'    
         ! Time derivative of rhoCT
         do n=1,this%CTelem
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%cfg%imin_,this%cfg%imax_
                     drhoCTdt(i,j,k,n)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k,n))+&
                     &                 sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k,n))+&
                     &                 sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1,n))+&
                     &                 this%CTgradU(i,j,k,n)+this%CTgradUT(i,j,k,n)-this%T(i,j,k,n)
                  end do
               end do
            end do
         end do
      case(upwind)
         ! Flux of rhoCT
         do n=1,this%CTelem
            ! Fluxes on x-face
            FX(ii,jj,kk,n)=-0.5_WP*(rhoU(ii,jj,kk)+abs(rhoU(ii,jj,kk)))*sum(this%itp_xp(:,ii,jj,kk)*this%CT(ii+this%stp1:ii+this%stp2,jj,kk,n)) &
            &              -0.5_WP*(rhoU(ii,jj,kk)-abs(rhoU(ii,jj,kk)))*sum(this%itp_xm(:,ii,jj,kk)*this%CT(ii+this%stm1:ii+this%stm2,jj,kk,n))
            ! Fluxes on y-face
            FY(ii,jj,kk,n)=-0.5_WP*(rhoV(ii,jj,kk)+abs(rhoV(ii,jj,kk)))*sum(this%itp_yp(:,ii,jj,kk)*this%CT(ii,jj+this%stp1:jj+this%stp2,kk,n)) &
            &              -0.5_WP*(rhoV(ii,jj,kk)-abs(rhoV(ii,jj,kk)))*sum(this%itp_ym(:,ii,jj,kk)*this%CT(ii,jj+this%stm1:jj+this%stm2,kk,n)) 
            ! Fluxes on z-face
            FZ(ii,jj,kk,n)=-0.5_WP*(rhoW(ii,jj,kk)+abs(rhoW(ii,jj,kk)))*sum(this%itp_zp(:,ii,jj,kk)*this%CT(ii,jj,kk+this%stp1:kk+this%stp2,n)) &
            &              -0.5_WP*(rhoW(ii,jj,kk)-abs(rhoW(ii,jj,kk)))*sum(this%itp_zm(:,ii,jj,kk)*this%CT(ii,jj,kk+this%stm1:kk+this%stm2,n)) 
         end do
         ! Time derivative of rhoCT
         do n=1,this%CTelem
            drhoCTdt(ii,jj,kk,n)=sum(this%div_x(:,ii,jj,kk)*FX(ii:ii+1,jj,kk,n))+&
            &                    sum(this%div_y(:,ii,jj,kk)*FY(ii,jj:jj+1,kk,n))+&
            &                    sum(this%div_z(:,ii,jj,kk)*FZ(ii,jj,kk:kk+1,n))+&
            &                    this%CTgradU(ii,jj,kk,n)+this%CTgradUT(ii,jj,kk,n)-this%T(ii,jj,kk,n)
         end do 
      end select
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync residual
      call this%cfg%sync(drhoCTdt)
   end subroutine get_drhoCTdt

   !> Calculate the viscoelastic tensor divergence
   subroutine get_divT(this)
      implicit none
      class(fene), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%divT(i,j,k,1)=sum(this%div_x(:,i,j,k)*this%T(i:i+1,j,k,1))+&
               &                  sum(this%div_y(:,i,j,k)*this%T(i,j:j+1,k,2))+&
               &                  sum(this%div_z(:,i,j,k)*this%T(i,j,k:k+1,3))
               this%divT(i,j,k,2)=sum(this%div_x(:,i,j,k)*this%T(i:i+1,j,k,4))+&
               &                  sum(this%div_y(:,i,j,k)*this%T(i,j:j+1,k,5))+&
               &                  sum(this%div_z(:,i,j,k)*this%T(i,j,k:k+1,6))
               this%divT(i,j,k,3)=sum(this%div_x(:,i,j,k)*this%T(i:i+1,j,k,7))+&
               &                  sum(this%div_y(:,i,j,k)*this%T(i,j:j+1,k,8))+&
               &                  sum(this%div_z(:,i,j,k)*this%T(i,j,k:k+1,9))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(this%divT)
   end subroutine get_divT
   
   
   !> Solve for implicit fene residual
   subroutine solve_implicit(this,scheme,dt,resCT,rhoU,rhoV,rhoW,ii,jj,kk)
      implicit none
      class(fene), intent(inout) :: this
      integer,  intent(in) :: scheme
      integer,  intent(in) :: ii,jj,kk
      real(WP), intent(in) :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resCT !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)    :: rhoU  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)    :: rhoV  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:   ), intent(in)    :: rhoW  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,sti,std,n
      
      select case (scheme)
      case(quick)
         ! Prepare convective operator
         do n=1,this%CTelem
            this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  do i=this%cfg%imin_,this%cfg%imax_
                     ! Loop over divergence stencil
                     do std=0,1
                        ! Loop over plus interpolation stencil
                        do sti=this%stp1,this%stp2
                           this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)+abs(rhoU(i+std,j,k)))*this%itp_xp(sti,i+std,j,k)
                           this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)+abs(rhoV(i,j+std,k)))*this%itp_yp(sti,i,j+std,k)
                           this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)+abs(rhoW(i,j,k+std)))*this%itp_zp(sti,i,j,k+std)
                        end do
                        ! Loop over minus interpolation stencil
                        do sti=this%stm1,this%stm2
                           this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%div_x(std,i,j,k)*0.5_WP*(rhoU(i+std,j,k)-abs(rhoU(i+std,j,k)))*this%itp_xm(sti,i+std,j,k)
                           this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%div_y(std,i,j,k)*0.5_WP*(rhoV(i,j+std,k)-abs(rhoV(i,j+std,k)))*this%itp_ym(sti,i,j+std,k)
                           this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%div_z(std,i,j,k)*0.5_WP*(rhoW(i,j,k+std)-abs(rhoW(i,j,k+std)))*this%itp_zm(sti,i,j,k+std)
                        end do
                     end do
                  end do
               end do
            end do
            ! Solve the linear system
            call this%implicit%setup()
            this%implicit%rhs=resCT(:,:,:,n)
            this%implicit%sol=0.0_WP
            call this%implicit%solve()
            resCT(:,:,:,n)=this%implicit%sol
         end do
      case(upwind)
         ! Prepare convective operator
         do n=1,this%CTelem
            this%implicit%opr(1,ii,jj,kk)=this%rho; this%implicit%opr(2:,ii,jj,kk)=0.0_WP
            ! Loop over divergence stencil
            do std=0,1
               ! Loop over plus interpolation stencil
               do sti=this%stp1,this%stp2
                  this%implicit%opr(this%implicit%stmap(sti+std,0,0),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),ii,jj,kk)+0.5_WP*dt*this%div_x(std,ii,jj,kk)*0.5_WP*(rhoU(ii+std,jj,kk)+abs(rhoU(ii+std,jj,kk)))*this%itp_xp(sti,ii+std,jj,kk)
                  this%implicit%opr(this%implicit%stmap(0,sti+std,0),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),ii,jj,kk)+0.5_WP*dt*this%div_y(std,ii,jj,kk)*0.5_WP*(rhoV(ii,jj+std,kk)+abs(rhoV(ii,jj+std,kk)))*this%itp_yp(sti,ii,jj+std,kk)
                  this%implicit%opr(this%implicit%stmap(0,0,sti+std),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),ii,jj,kk)+0.5_WP*dt*this%div_z(std,ii,jj,kk)*0.5_WP*(rhoW(ii,jj,kk+std)+abs(rhoW(ii,jj,kk+std)))*this%itp_zp(sti,ii,jj,kk+std)
               end do
               ! Loop over minus interpolation stencil
               do sti=this%stm1,this%stm2
                  this%implicit%opr(this%implicit%stmap(sti+std,0,0),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),ii,jj,kk)+0.5_WP*dt*this%div_x(std,ii,jj,kk)*0.5_WP*(rhoU(ii+std,jj,kk)-abs(rhoU(ii+std,jj,kk)))*this%itp_xm(sti,ii+std,jj,kk)
                  this%implicit%opr(this%implicit%stmap(0,sti+std,0),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),ii,jj,kk)+0.5_WP*dt*this%div_y(std,ii,jj,kk)*0.5_WP*(rhoV(ii,jj+std,kk)-abs(rhoV(ii,jj+std,kk)))*this%itp_ym(sti,ii,jj+std,kk)
                  this%implicit%opr(this%implicit%stmap(0,0,sti+std),ii,jj,kk)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),ii,jj,kk)+0.5_WP*dt*this%div_z(std,ii,jj,kk)*0.5_WP*(rhoW(ii,jj,kk+std)-abs(rhoW(ii,jj,kk+std)))*this%itp_zm(sti,ii,jj,kk+std)
               end do
            end do
            ! Solve the linear system
            call this%implicit%setup()
            this%implicit%rhs=resCT(ii,jj,kk,n)
            this%implicit%sol=0.0_WP
            call this%implicit%solve()
            resCT(ii,jj,kk,n)=this%implicit%sol(ii,jj,kk)
         end do
      end select
      
      ! Sync up residual
      call this%cfg%sync(resCT)
      
   end subroutine solve_implicit
   
   
   !> Print out info for fene solver
   subroutine fene_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(fene), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Constant density FENE model solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine fene_print
   
   
end module fene_class
