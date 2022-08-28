!> Tensor solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes constant diffusivity and density.
module tensor_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ils_class,      only: ils
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: tensor,bcond
   
   ! List of known available bcond for this solver
   integer, parameter, public :: dirichlet=2         !< Dirichlet condition
   integer, parameter, public :: neumann=3           !< Zero normal gradient
   
   ! List of available advection schemes for tensor transport
   integer, parameter, public :: quick =1             !< Quick scheme
   integer, parameter, public :: upwind=2             !< Upwind scheme
   
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
   
   !> Constant density tensor solver object definition
   type :: tensor
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_tensor'  !< Solver name (default=UNNAMED_tensor)
      
      ! Constant property fluid, but diffusivity is still a field due to LES modeling
      real(WP) :: rho                                       !< This is our constant fluid density
      real(WP), dimension(:,:,:,:), allocatable :: diff     !< These is our constant+SGS dynamic diffusivity for the tensor
      
      ! Boundary condition list
      integer :: nbc                                      !< Number of bcond for our solver
      type(bcond), pointer :: first_bc                    !< List of bcond for our solver
      
      ! tensor variable
      real(WP), dimension(:,:,:,:), allocatable :: TN       !< TN array
      
      ! Old tensor variable
      real(WP), dimension(:,:,:,:), allocatable :: TNold    !< TNold array
      
      ! Implicit tensor solver
      type(ils) :: implicit                               !< Iterative linear solver object for an implicit prediction of the tensor residual
      integer, dimension(:,:,:), allocatable :: stmap     !< Inverse map from stencil shift to index location
      
      ! Metrics
      integer :: scheme                                   !< Advection scheme for tensor
      integer :: nst                                      !< Scheme order (and elemental stencil size)
      integer :: stp1,stp2                                !< Plus interpolation stencil extent for tensor advection
      integer :: stm1,stm2                                !< Minus interpolation stencil extent for tensor advection
      real(WP), dimension(:,:,:,:,:), allocatable :: itptn_xp,itptn_yp,itptn_zp   !< Plus interpolation for TN
      real(WP), dimension(:,:,:,:,:), allocatable :: itptn_xm,itptn_ym,itptn_zm   !< Minus interpolation for TN
      real(WP), dimension(:,:,:,:,:), allocatable :: divtn_x ,divtn_y ,divtn_z    !< Divergence for TN
      real(WP), dimension(:,:,:,:,:), allocatable :: grdtn_x ,grdtn_y ,grdtn_z    !< tensor gradient for TN
      real(WP), dimension(:,:,:,:,:), allocatable :: itp_x   ,itp_y   ,itp_z      !< Second order interpolation for TN diffusivity
      
      ! Masking info for metric modification
      integer, dimension(:,:,:,:), allocatable :: mask    !< Integer array used for modifying TN metrics
      
      ! Monitoring quantities
      real(WP) :: TNmax,TNmin,TNint                       !< Maximum and minimum, integral tensor
      
   contains
      procedure :: print=>tensor_print                    !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the tensor solver
      procedure :: add_bcond                              !< Add a boundary condition
      procedure :: get_bcond                              !< Get a boundary condition
      procedure :: apply_bcond                            !< Apply all boundary conditions
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_drhoTNdt                           !< Calculate drhoTN/dt
      procedure :: get_max                                !< Calculate maximum field values
      ! procedure :: get_int                                !< Calculate integral field values
      procedure :: solve_implicit                         !< Solve for the tensor residuals implicitly
   end type tensor
   
   
   !> Declare tensor solver constructor
   interface tensor
      procedure constructor
   end interface tensor
   
contains
   
   
   !> Default constructor for tensor solver
   function constructor(cfg,scheme,name) result(self)
      use messager, only: die
      implicit none
      type(tensor) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: scheme
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
      allocate(self%TN   (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_),6); self%TN   =0.0_WP
      allocate(self%TNold(6,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%TNold=0.0_WP
      allocate(self%diff (6,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%diff =0.0_WP
      
      ! Prepare advection scheme
      self%scheme=scheme
      select case (self%scheme)
      case (quick)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[tensor constructor] tensor transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=3
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case (upwind)
         ! Check current overlap
         if (self%cfg%no.lt.2) call die('[tensor constructor] tensor transport scheme requires larger overlap')
         ! Set interpolation stencil sizes
         self%nst=1
         self%stp1=-(self%nst+1)/2; self%stp2=self%nst+self%stp1-1
         self%stm1=-(self%nst-1)/2; self%stm2=self%nst+self%stm1-1
      case default
         call die('[tensor constructor] Unknown tensor transport scheme selected')
      end select
      
      ! Create implicit tensor solver object
      self%implicit=ils(cfg=self%cfg,name='tensor',nst=1+6*abs(self%stp1))
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for TN
      allocate(self%mask(6,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%mask=0
      if (.not.self%cfg%xper) then
         if (self%cfg%iproc.eq.           1) self%mask(1,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.           1) self%mask(2,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.           1) self%mask(3,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.           1) self%mask(4,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.           1) self%mask(5,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.           1) self%mask(6,:self%cfg%imin-1,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(1,self%cfg%imax+1:,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(2,self%cfg%imax+1:,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(3,self%cfg%imax+1:,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(4,self%cfg%imax+1:,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(5,self%cfg%imax+1:,:,:)=2
         if (self%cfg%iproc.eq.self%cfg%npx) self%mask(6,self%cfg%imax+1:,:,:)=2
      end if
      if (.not.self%cfg%yper) then
         if (self%cfg%jproc.eq.           1) self%mask(1,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.           1) self%mask(2,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.           1) self%mask(3,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.           1) self%mask(4,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.           1) self%mask(5,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.           1) self%mask(6,:,:self%cfg%jmin-1,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(1,:,self%cfg%jmax+1:,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(2,:,self%cfg%jmax+1:,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(3,:,self%cfg%jmax+1:,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(4,:,self%cfg%jmax+1:,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(5,:,self%cfg%jmax+1:,:)=2
         if (self%cfg%jproc.eq.self%cfg%npy) self%mask(6,:,self%cfg%jmax+1:,:)=2
      end if
      if (.not.self%cfg%zper) then
         if (self%cfg%kproc.eq.           1) self%mask(1,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.           1) self%mask(2,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.           1) self%mask(3,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.           1) self%mask(4,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.           1) self%mask(5,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.           1) self%mask(6,:,:,:self%cfg%kmin-1)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(1,:,:,self%cfg%kmax+1:)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(2,:,:,self%cfg%kmax+1:)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(3,:,:,self%cfg%kmax+1:)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(4,:,:,self%cfg%kmax+1:)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(5,:,:,self%cfg%kmax+1:)=2
         if (self%cfg%kproc.eq.self%cfg%npz) self%mask(6,:,:,self%cfg%kmax+1:)=2
      end if
      do k=self%cfg%kmino_,self%cfg%kmaxo_
         do j=self%cfg%jmino_,self%cfg%jmaxo_
            do i=self%cfg%imino_,self%cfg%imaxo_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(1,i,j,k)=1
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(2,i,j,k)=1
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(3,i,j,k)=1
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(4,i,j,k)=1
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(5,i,j,k)=1
               if (self%cfg%VF(i,j,k).eq.0.0_WP) self%mask(6,i,j,k)=1
            end do
         end do
      end do
      call self%cfg%sync(self%mask(1,:,:,:))
      call self%cfg%sync(self%mask(2,:,:,:))
      call self%cfg%sync(self%mask(3,:,:,:))
      call self%cfg%sync(self%mask(4,:,:,:))
      call self%cfg%sync(self%mask(5,:,:,:))
      call self%cfg%sync(self%mask(6,:,:,:))
      
   end function constructor
      
   
   !> Metric initialization with no awareness of walls nor bcond
   subroutine init_metrics(this)
      use mathtools, only: fv_itp_build
      implicit none
      class(tensor), intent(inout) :: this
      integer :: i,j,k
      
      ! Allocate finite difference diffusivity interpolation coefficients
      allocate(this%itp_x(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itp_y(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itp_z(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create diffusivity interpolation coefficients to cell face
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%itp_x(1,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_x(2,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_x(3,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_x(4,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_x(5,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_x(6,:,i,j,k)=this%cfg%dxmi(i)*[this%cfg%xm(i)-this%cfg%x(i),this%cfg%x(i)-this%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
               this%itp_y(1,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_y(2,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_y(3,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_y(4,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_y(5,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_y(6,:,i,j,k)=this%cfg%dymi(j)*[this%cfg%ym(j)-this%cfg%y(j),this%cfg%y(j)-this%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
               this%itp_z(1,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
               this%itp_z(2,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
               this%itp_z(3,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
               this%itp_z(4,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
               this%itp_z(5,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
               this%itp_z(6,:,i,j,k)=this%cfg%dzmi(k)*[this%cfg%zm(k)-this%cfg%z(k),this%cfg%z(k)-this%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
      ! Allocate finite difference tensor interpolation coefficients
      allocate(this%itptn_xp(6,this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itptn_xm(6,this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%itptn_yp(6,this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itptn_ym(6,this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%itptn_zp(6,this%stp1:this%stp2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      allocate(this%itptn_zm(6,this%stm1:this%stm2,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create tensor interpolation coefficients to cell faces
      select case (this%scheme)
      case (quick)
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(6,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(6,:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(6,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(6,:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(6,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(1,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(2,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(3,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(4,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(5,:,i,j,k))
                  call fv_itp_build(n=3,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(6,:,i,j,k))
               end do
            end do
         end do
      case (upwind)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Interpolation to x-face
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stp1:i+this%stp2+1),xp=this%cfg%x(i),coeff=this%itptn_xp(6,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%x(i+this%stm1:i+this%stm2+1),xp=this%cfg%x(i),coeff=this%itptn_xm(6,:,i,j,k))
                  ! Interpolation to y-face
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stp1:j+this%stp2+1),xp=this%cfg%y(j),coeff=this%itptn_yp(6,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%y(j+this%stm1:j+this%stm2+1),xp=this%cfg%y(j),coeff=this%itptn_ym(6,:,i,j,k))
                  ! Interpolation to z-face
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stp1:k+this%stp2+1),xp=this%cfg%z(k),coeff=this%itptn_zp(6,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(1,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(2,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(3,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(4,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(5,:,i,j,k))
                  call fv_itp_build(n=1,x=this%cfg%z(k+this%stm1:k+this%stm2+1),xp=this%cfg%z(k),coeff=this%itptn_zm(6,:,i,j,k))
               end do
            end do
         end do
      end select
      
      ! Allocate finite volume divergence operators
      allocate(this%divtn_x(6,0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%divtn_y(6,0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      allocate(this%divtn_z(6,0:+1,this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%divtn_x(1,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_x(2,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_x(3,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_x(4,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_x(5,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_x(6,:,i,j,k)=this%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FV divergence from [x ,ym,zm]
               this%divtn_y(1,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_y(2,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_y(3,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_y(4,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_y(5,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_y(6,:,i,j,k)=this%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,y ,zm]
               this%divtn_z(1,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               this%divtn_z(2,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               this%divtn_z(3,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               this%divtn_z(4,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               this%divtn_z(5,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
               this%divtn_z(6,:,i,j,k)=this%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FV divergence from [xm,ym,z ]
            end do
         end do
      end do
      
      ! Allocate finite difference velocity gradient operators
      allocate(this%grdtn_x(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< X-face-centered
      allocate(this%grdtn_y(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Y-face-centered
      allocate(this%grdtn_z(6,-1:0,this%cfg%imin_:this%cfg%imax_+1,this%cfg%jmin_:this%cfg%jmax_+1,this%cfg%kmin_:this%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               this%grdtn_x(1,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_x(2,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_x(3,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_x(4,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_x(5,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_x(6,:,i,j,k)=this%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in x from [xm,ym,zm] to [x,ym,zm]
               this%grdtn_y(1,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_y(2,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_y(3,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_y(4,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_y(5,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_y(6,:,i,j,k)=this%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in y from [xm,ym,zm] to [xm,y,zm]
               this%grdtn_z(1,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
               this%grdtn_z(2,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
               this%grdtn_z(3,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
               this%grdtn_z(4,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
               this%grdtn_z(5,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
               this%grdtn_z(6,:,i,j,k)=this%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of TN in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do
      
   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(tensor), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask(1,:,:,:))
      call this%cfg%sync(this%mask(2,:,:,:))
      call this%cfg%sync(this%mask(3,:,:,:))
      call this%cfg%sync(this%mask(4,:,:,:))
      call this%cfg%sync(this%mask(5,:,:,:))
      call this%cfg%sync(this%mask(6,:,:,:))
      
      ! Adjust interpolation coefficients to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Linear interpolation in x
               if (this%mask(1,i,j,k).eq.0.and.this%mask(1,i-1,j,k).gt.0) this%itp_x(1,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(2,i,j,k).eq.0.and.this%mask(2,i-1,j,k).gt.0) this%itp_x(2,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(3,i,j,k).eq.0.and.this%mask(3,i-1,j,k).gt.0) this%itp_x(3,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(4,i,j,k).eq.0.and.this%mask(4,i-1,j,k).gt.0) this%itp_x(4,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(5,i,j,k).eq.0.and.this%mask(5,i-1,j,k).gt.0) this%itp_x(5,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(6,i,j,k).eq.0.and.this%mask(6,i-1,j,k).gt.0) this%itp_x(6,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(1,i,j,k).gt.0.and.this%mask(1,i-1,j,k).eq.0) this%itp_x(1,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(2,i,j,k).gt.0.and.this%mask(2,i-1,j,k).eq.0) this%itp_x(2,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(3,i,j,k).gt.0.and.this%mask(3,i-1,j,k).eq.0) this%itp_x(3,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(4,i,j,k).gt.0.and.this%mask(4,i-1,j,k).eq.0) this%itp_x(4,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(5,i,j,k).gt.0.and.this%mask(5,i-1,j,k).eq.0) this%itp_x(5,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(6,i,j,k).gt.0.and.this%mask(6,i-1,j,k).eq.0) this%itp_x(6,:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in y
               if (this%mask(1,i,j,k).eq.0.and.this%mask(1,i,j-1,k).gt.0) this%itp_y(1,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(2,i,j,k).eq.0.and.this%mask(2,i,j-1,k).gt.0) this%itp_y(2,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(3,i,j,k).eq.0.and.this%mask(3,i,j-1,k).gt.0) this%itp_y(3,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(4,i,j,k).eq.0.and.this%mask(4,i,j-1,k).gt.0) this%itp_y(4,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(5,i,j,k).eq.0.and.this%mask(5,i,j-1,k).gt.0) this%itp_y(5,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(6,i,j,k).eq.0.and.this%mask(6,i,j-1,k).gt.0) this%itp_y(6,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(1,i,j,k).gt.0.and.this%mask(1,i,j-1,k).eq.0) this%itp_y(1,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(2,i,j,k).gt.0.and.this%mask(2,i,j-1,k).eq.0) this%itp_y(2,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(3,i,j,k).gt.0.and.this%mask(3,i,j-1,k).eq.0) this%itp_y(3,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(4,i,j,k).gt.0.and.this%mask(4,i,j-1,k).eq.0) this%itp_y(4,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(5,i,j,k).gt.0.and.this%mask(5,i,j-1,k).eq.0) this%itp_y(5,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(6,i,j,k).gt.0.and.this%mask(6,i,j-1,k).eq.0) this%itp_y(6,:,i,j,k)=[1.0_WP,0.0_WP]
               ! Linear interpolation in z
               if (this%mask(1,i,j,k).eq.0.and.this%mask(1,i,j,k-1).gt.0) this%itp_z(1,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(2,i,j,k).eq.0.and.this%mask(2,i,j,k-1).gt.0) this%itp_z(2,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(3,i,j,k).eq.0.and.this%mask(3,i,j,k-1).gt.0) this%itp_z(3,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(4,i,j,k).eq.0.and.this%mask(4,i,j,k-1).gt.0) this%itp_z(4,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(5,i,j,k).eq.0.and.this%mask(5,i,j,k-1).gt.0) this%itp_z(5,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(6,i,j,k).eq.0.and.this%mask(6,i,j,k-1).gt.0) this%itp_z(6,:,i,j,k)=[0.0_WP,1.0_WP]
               if (this%mask(1,i,j,k).gt.0.and.this%mask(1,i,j,k-1).eq.0) this%itp_z(1,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(2,i,j,k).gt.0.and.this%mask(2,i,j,k-1).eq.0) this%itp_z(2,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(3,i,j,k).gt.0.and.this%mask(3,i,j,k-1).eq.0) this%itp_z(3,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(4,i,j,k).gt.0.and.this%mask(4,i,j,k-1).eq.0) this%itp_z(4,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(5,i,j,k).gt.0.and.this%mask(5,i,j,k-1).eq.0) this%itp_z(5,:,i,j,k)=[1.0_WP,0.0_WP]
               if (this%mask(6,i,j,k).gt.0.and.this%mask(6,i,j,k-1).eq.0) this%itp_z(6,:,i,j,k)=[1.0_WP,0.0_WP]
            end do
         end do
      end do
      
      ! Adjust tensor interpolation to reflect Dirichlet boundaries
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! X face - 1
               if (this%mask(1,i-1,j,k).eq.2) then
                  this%itptn_xm(1,:,i,j,k)=0.0_WP; this%itptn_xm(1,-1,i,j,k)=1.0_WP
                  this%itptn_xp(1,:,i,j,k)=0.0_WP; this%itptn_xp(1,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(1,i  ,j,k).eq.2) then
                  this%itptn_xm(1,:,i,j,k)=0.0_WP; this%itptn_xm(1, 0,i,j,k)=1.0_WP
                  this%itptn_xp(1,:,i,j,k)=0.0_WP; this%itptn_xp(1, 0,i,j,k)=1.0_WP
               end if
               ! X face - 2
               if (this%mask(2,i-1,j,k).eq.2) then
                  this%itptn_xm(2,:,i,j,k)=0.0_WP; this%itptn_xm(2,-1,i,j,k)=1.0_WP
                  this%itptn_xp(2,:,i,j,k)=0.0_WP; this%itptn_xp(2,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(2,i  ,j,k).eq.2) then
                  this%itptn_xm(2,:,i,j,k)=0.0_WP; this%itptn_xm(2, 0,i,j,k)=1.0_WP
                  this%itptn_xp(2,:,i,j,k)=0.0_WP; this%itptn_xp(2, 0,i,j,k)=1.0_WP
               end if
               ! X face - 3
               if (this%mask(3,i-1,j,k).eq.2) then
                  this%itptn_xm(3,:,i,j,k)=0.0_WP; this%itptn_xm(3,-1,i,j,k)=1.0_WP
                  this%itptn_xp(3,:,i,j,k)=0.0_WP; this%itptn_xp(3,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(3,i  ,j,k).eq.2) then
                  this%itptn_xm(3,:,i,j,k)=0.0_WP; this%itptn_xm(3, 0,i,j,k)=1.0_WP
                  this%itptn_xp(3,:,i,j,k)=0.0_WP; this%itptn_xp(3, 0,i,j,k)=1.0_WP
               end if
               ! X face - 4
               if (this%mask(4,i-1,j,k).eq.2) then
                  this%itptn_xm(4,:,i,j,k)=0.0_WP; this%itptn_xm(4,-1,i,j,k)=1.0_WP
                  this%itptn_xp(4,:,i,j,k)=0.0_WP; this%itptn_xp(4,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(4,i  ,j,k).eq.2) then
                  this%itptn_xm(4,:,i,j,k)=0.0_WP; this%itptn_xm(4, 0,i,j,k)=1.0_WP
                  this%itptn_xp(4,:,i,j,k)=0.0_WP; this%itptn_xp(4, 0,i,j,k)=1.0_WP
               end if
               ! X face - 5
               if (this%mask(5,i-1,j,k).eq.2) then
                  this%itptn_xm(5,:,i,j,k)=0.0_WP; this%itptn_xm(5,-1,i,j,k)=1.0_WP
                  this%itptn_xp(5,:,i,j,k)=0.0_WP; this%itptn_xp(5,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(5,i  ,j,k).eq.2) then
                  this%itptn_xm(5,:,i,j,k)=0.0_WP; this%itptn_xm(5, 0,i,j,k)=1.0_WP
                  this%itptn_xp(5,:,i,j,k)=0.0_WP; this%itptn_xp(5, 0,i,j,k)=1.0_WP
               end if
               ! X face - 6
               if (this%mask(6,i-1,j,k).eq.2) then
                  this%itptn_xm(6,:,i,j,k)=0.0_WP; this%itptn_xm(6,-1,i,j,k)=1.0_WP
                  this%itptn_xp(6,:,i,j,k)=0.0_WP; this%itptn_xp(6,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(6,i  ,j,k).eq.2) then
                  this%itptn_xm(6,:,i,j,k)=0.0_WP; this%itptn_xm(6, 0,i,j,k)=1.0_WP
                  this%itptn_xp(6,:,i,j,k)=0.0_WP; this%itptn_xp(6, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 1
               if (this%mask(1,i,j-1,k).eq.2) then
                  this%itptn_ym(1,:,i,j,k)=0.0_WP; this%itptn_ym(1,-1,i,j,k)=1.0_WP
                  this%itptn_yp(1,:,i,j,k)=0.0_WP; this%itptn_yp(1,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(1,i,j  ,k).eq.2) then
                  this%itptn_ym(1,:,i,j,k)=0.0_WP; this%itptn_ym(1, 0,i,j,k)=1.0_WP
                  this%itptn_yp(1,:,i,j,k)=0.0_WP; this%itptn_yp(1, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 2
               if (this%mask(2,i,j-1,k).eq.2) then
                  this%itptn_ym(2,:,i,j,k)=0.0_WP; this%itptn_ym(2,-1,i,j,k)=1.0_WP
                  this%itptn_yp(2,:,i,j,k)=0.0_WP; this%itptn_yp(2,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(2,i,j  ,k).eq.2) then
                  this%itptn_ym(2,:,i,j,k)=0.0_WP; this%itptn_ym(2, 0,i,j,k)=1.0_WP
                  this%itptn_yp(2,:,i,j,k)=0.0_WP; this%itptn_yp(2, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 3
               if (this%mask(3,i,j-1,k).eq.2) then
                  this%itptn_ym(3,:,i,j,k)=0.0_WP; this%itptn_ym(3,-1,i,j,k)=1.0_WP
                  this%itptn_yp(3,:,i,j,k)=0.0_WP; this%itptn_yp(3,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(3,i,j  ,k).eq.2) then
                  this%itptn_ym(3,:,i,j,k)=0.0_WP; this%itptn_ym(3, 0,i,j,k)=1.0_WP
                  this%itptn_yp(3,:,i,j,k)=0.0_WP; this%itptn_yp(3, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 4
               if (this%mask(4,i,j-1,k).eq.2) then
                  this%itptn_ym(4,:,i,j,k)=0.0_WP; this%itptn_ym(4,-1,i,j,k)=1.0_WP
                  this%itptn_yp(4,:,i,j,k)=0.0_WP; this%itptn_yp(4,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(4,i,j  ,k).eq.2) then
                  this%itptn_ym(4,:,i,j,k)=0.0_WP; this%itptn_ym(4, 0,i,j,k)=1.0_WP
                  this%itptn_yp(4,:,i,j,k)=0.0_WP; this%itptn_yp(4, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 5
               if (this%mask(5,i,j-1,k).eq.2) then
                  this%itptn_ym(5,:,i,j,k)=0.0_WP; this%itptn_ym(5,-1,i,j,k)=1.0_WP
                  this%itptn_yp(5,:,i,j,k)=0.0_WP; this%itptn_yp(5,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(5,i,j  ,k).eq.2) then
                  this%itptn_ym(5,:,i,j,k)=0.0_WP; this%itptn_ym(5, 0,i,j,k)=1.0_WP
                  this%itptn_yp(5,:,i,j,k)=0.0_WP; this%itptn_yp(5, 0,i,j,k)=1.0_WP
               end if
               ! Y face - 6
               if (this%mask(6,i,j-1,k).eq.2) then
                  this%itptn_ym(6,:,i,j,k)=0.0_WP; this%itptn_ym(6,-1,i,j,k)=1.0_WP
                  this%itptn_yp(6,:,i,j,k)=0.0_WP; this%itptn_yp(6,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(6,i,j  ,k).eq.2) then
                  this%itptn_ym(6,:,i,j,k)=0.0_WP; this%itptn_ym(6, 0,i,j,k)=1.0_WP
                  this%itptn_yp(6,:,i,j,k)=0.0_WP; this%itptn_yp(6, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 1
               if (this%mask(1,i,j,k-1).eq.2) then
                  this%itptn_zm(1,:,i,j,k)=0.0_WP; this%itptn_zm(1,-1,i,j,k)=1.0_WP
                  this%itptn_zp(1,:,i,j,k)=0.0_WP; this%itptn_zp(1,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(1,i,j,k  ).eq.2) then
                  this%itptn_zm(1,:,i,j,k)=0.0_WP; this%itptn_zm(1, 0,i,j,k)=1.0_WP
                  this%itptn_zp(1,:,i,j,k)=0.0_WP; this%itptn_zp(1, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 2
               if (this%mask(2,i,j,k-1).eq.2) then
                  this%itptn_zm(2,:,i,j,k)=0.0_WP; this%itptn_zm(2,-1,i,j,k)=1.0_WP
                  this%itptn_zp(2,:,i,j,k)=0.0_WP; this%itptn_zp(2,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(2,i,j,k  ).eq.2) then
                  this%itptn_zm(2,:,i,j,k)=0.0_WP; this%itptn_zm(2, 0,i,j,k)=1.0_WP
                  this%itptn_zp(2,:,i,j,k)=0.0_WP; this%itptn_zp(2, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 3
               if (this%mask(3,i,j,k-1).eq.2) then
                  this%itptn_zm(3,:,i,j,k)=0.0_WP; this%itptn_zm(3,-1,i,j,k)=1.0_WP
                  this%itptn_zp(3,:,i,j,k)=0.0_WP; this%itptn_zp(3,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(3,i,j,k  ).eq.2) then
                  this%itptn_zm(3,:,i,j,k)=0.0_WP; this%itptn_zm(3, 0,i,j,k)=1.0_WP
                  this%itptn_zp(3,:,i,j,k)=0.0_WP; this%itptn_zp(3, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 4
               if (this%mask(4,i,j,k-1).eq.2) then
                  this%itptn_zm(4,:,i,j,k)=0.0_WP; this%itptn_zm(4,-1,i,j,k)=1.0_WP
                  this%itptn_zp(4,:,i,j,k)=0.0_WP; this%itptn_zp(4,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(4,i,j,k  ).eq.2) then
                  this%itptn_zm(4,:,i,j,k)=0.0_WP; this%itptn_zm(4, 0,i,j,k)=1.0_WP
                  this%itptn_zp(4,:,i,j,k)=0.0_WP; this%itptn_zp(4, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 5
               if (this%mask(5,i,j,k-1).eq.2) then
                  this%itptn_zm(5,:,i,j,k)=0.0_WP; this%itptn_zm(5,-1,i,j,k)=1.0_WP
                  this%itptn_zp(5,:,i,j,k)=0.0_WP; this%itptn_zp(5,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(5,i,j,k  ).eq.2) then
                  this%itptn_zm(5,:,i,j,k)=0.0_WP; this%itptn_zm(5, 0,i,j,k)=1.0_WP
                  this%itptn_zp(5,:,i,j,k)=0.0_WP; this%itptn_zp(5, 0,i,j,k)=1.0_WP
               end if
               ! Z face - 6
               if (this%mask(6,i,j,k-1).eq.2) then
                  this%itptn_zm(6,:,i,j,k)=0.0_WP; this%itptn_zm(6,-1,i,j,k)=1.0_WP
                  this%itptn_zp(6,:,i,j,k)=0.0_WP; this%itptn_zp(6,-1,i,j,k)=1.0_WP
               end if
               if (this%mask(6,i,j,k  ).eq.2) then
                  this%itptn_zm(6,:,i,j,k)=0.0_WP; this%itptn_zm(6, 0,i,j,k)=1.0_WP
                  this%itptn_zp(6,:,i,j,k)=0.0_WP; this%itptn_zp(6, 0,i,j,k)=1.0_WP
               end if
            end do
         end do
      end do
      
      ! Loop over the domain and apply masked conditions to TN divergence
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(1,i,j,k).gt.0) then
                  this%divtn_x(1,:,i,j,k)=0.0_WP
                  this%divtn_y(1,:,i,j,k)=0.0_WP
                  this%divtn_z(1,:,i,j,k)=0.0_WP
               end if
               if (this%mask(2,i,j,k).gt.0) then
                  this%divtn_x(2,:,i,j,k)=0.0_WP
                  this%divtn_y(2,:,i,j,k)=0.0_WP
                  this%divtn_z(2,:,i,j,k)=0.0_WP
               end if
               if (this%mask(3,i,j,k).gt.0) then
                  this%divtn_x(3,:,i,j,k)=0.0_WP
                  this%divtn_y(3,:,i,j,k)=0.0_WP
                  this%divtn_z(3,:,i,j,k)=0.0_WP
               end if
               if (this%mask(4,i,j,k).gt.0) then
                  this%divtn_x(4,:,i,j,k)=0.0_WP
                  this%divtn_y(4,:,i,j,k)=0.0_WP
                  this%divtn_z(4,:,i,j,k)=0.0_WP
               end if
               if (this%mask(5,i,j,k).gt.0) then
                  this%divtn_x(5,:,i,j,k)=0.0_WP
                  this%divtn_y(5,:,i,j,k)=0.0_WP
                  this%divtn_z(5,:,i,j,k)=0.0_WP
               end if
               if (this%mask(6,i,j,k).gt.0) then
                  this%divtn_x(6,:,i,j,k)=0.0_WP
                  this%divtn_y(6,:,i,j,k)=0.0_WP
                  this%divtn_z(6,:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Adjust gradient coefficients to cell faces for walls (assume Neumann at wall)
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               if (this%mask(1,i,j,k).eq.1.or.this%mask(1,i-1,j,k).eq.1) this%grdtn_x(1,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(2,i,j,k).eq.1.or.this%mask(2,i-1,j,k).eq.1) this%grdtn_x(2,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(3,i,j,k).eq.1.or.this%mask(3,i-1,j,k).eq.1) this%grdtn_x(3,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(4,i,j,k).eq.1.or.this%mask(4,i-1,j,k).eq.1) this%grdtn_x(4,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(5,i,j,k).eq.1.or.this%mask(5,i-1,j,k).eq.1) this%grdtn_x(5,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(6,i,j,k).eq.1.or.this%mask(6,i-1,j,k).eq.1) this%grdtn_x(6,:,i,j,k)=0.0_WP     !< FD gradient in x of TN
               if (this%mask(1,i,j,k).eq.1.or.this%mask(1,i,j-1,k).eq.1) this%grdtn_y(1,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(2,i,j,k).eq.1.or.this%mask(2,i,j-1,k).eq.1) this%grdtn_y(2,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(3,i,j,k).eq.1.or.this%mask(3,i,j-1,k).eq.1) this%grdtn_y(3,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(4,i,j,k).eq.1.or.this%mask(4,i,j-1,k).eq.1) this%grdtn_y(4,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(5,i,j,k).eq.1.or.this%mask(5,i,j-1,k).eq.1) this%grdtn_y(5,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(6,i,j,k).eq.1.or.this%mask(6,i,j-1,k).eq.1) this%grdtn_y(6,:,i,j,k)=0.0_WP     !< FD gradient in y of TN
               if (this%mask(1,i,j,k).eq.1.or.this%mask(1,i,j,k-1).eq.1) this%grdtn_z(1,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
               if (this%mask(2,i,j,k).eq.1.or.this%mask(2,i,j,k-1).eq.1) this%grdtn_z(2,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
               if (this%mask(3,i,j,k).eq.1.or.this%mask(3,i,j,k-1).eq.1) this%grdtn_z(3,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
               if (this%mask(4,i,j,k).eq.1.or.this%mask(4,i,j,k-1).eq.1) this%grdtn_z(4,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
               if (this%mask(5,i,j,k).eq.1.or.this%mask(5,i,j,k-1).eq.1) this%grdtn_z(5,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
               if (this%mask(6,i,j,k).eq.1.or.this%mask(6,i,j,k-1).eq.1) this%grdtn_z(6,:,i,j,k)=0.0_WP     !< FD gradient in z of TN
            end do
         end do
      end do
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%divtn_x=0.0_WP
         this%grdtn_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%divtn_y=0.0_WP
         this%grdtn_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%divtn_z=0.0_WP
         this%grdtn_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the tensor solver now that bconds have been defined
   subroutine setup(this,implicit_ils)
      implicit none
      class(tensor), intent(inout) :: this
      integer, intent(in) :: implicit_ils
      integer :: count,st
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
      ! Set dynamic stencil map for the tensor solver
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
      
      ! Initialize the implicit tensor solver
      call this%implicit%init(implicit_ils)
      
   end subroutine setup
   
   
   !> Add a boundary condition
   subroutine add_bcond(this,name,type,locator,dir)
      use string,   only: lowercase
      use messager, only: die
      implicit none
      class(tensor), intent(inout) :: this
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
         case default; call die('[tensor add_bcond] Unknown bcond direction')
         end select
      else
         if (new_bc%type.eq.neumann) call die('[tensor apply_bcond] Neumann requires a direction')
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
            this%mask(1,i,j,k)=2
            this%mask(2,i,j,k)=2
            this%mask(3,i,j,k)=2
            this%mask(4,i,j,k)=2
            this%mask(5,i,j,k)=2
            this%mask(6,i,j,k)=2
         end do
      case (neumann)
         ! No modification - this assumes Neumann is only applied at walls or domain boundaries
      case default
         call die('[tensor apply_bcond] Unknown bcond type')
      end select
   
   end subroutine add_bcond
   
   
   !> Get a boundary condition
   subroutine get_bcond(this,name,my_bc)
      use messager, only: die
      implicit none
      class(tensor), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(bcond), pointer, intent(out) :: my_bc
      my_bc=>this%first_bc
      search: do while (associated(my_bc))
         if (trim(my_bc%name).eq.trim(name)) exit search
         my_bc=>my_bc%next
      end do search
      if (.not.associated(my_bc)) call die('[tensor get_bcond] Boundary condition was not found')
   end subroutine get_bcond
   
   
   !> Enforce boundary condition
   subroutine apply_bcond(this,t,dt)
      use messager, only: die
      use mpi_f08,  only: MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tensor), intent(inout) :: this
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
                  this%TN(1,i,j,k)=this%TN(1,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%TN(2,i,j,k)=this%TN(2,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%TN(3,i,j,k)=this%TN(3,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%TN(4,i,j,k)=this%TN(4,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%TN(5,i,j,k)=this%TN(5,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
                  this%TN(6,i,j,k)=this%TN(6,i-shift(1,my_bc%dir),j-shift(2,my_bc%dir),k-shift(3,my_bc%dir))
               end do
               
            case default
               call die('[tensor apply_bcond] Unknown bcond type')
            end select
            
         end if
         
         ! Sync full fields after each bcond - this should be optimized
         call this%cfg%sync(this%TN)
         
         ! Move on to the next bcond
         my_bc=>my_bc%next
         
      end do
      
   end subroutine apply_bcond
   
   
   !> Calculate the explicit rhoTN time derivative based on rhoU/rhoV/rhoW
   subroutine get_drhoTNdt(this,drhoTNdt,rhoU,rhoV,rhoW)
      implicit none
      class(tensor), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(out) :: drhoTNdt !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoU     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoV     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)  :: rhoW     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k
      real(WP), dimension(:,:,:,:), allocatable :: FX,FY,FZ
      ! Allocate flux arrays
      allocate(FX(6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Flux of rhoTN
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               ! Fluxes on x-face - 1
               FX(1,i,j,k)=-0.5_WP*(rhoU(1,i,j,k)+abs(rhoU(1,i,j,k)))*sum(this%itptn_xp(1,:,i,j,k)*this%TN(1,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(1,i,j,k)-abs(rhoU(1,i,j,k)))*sum(this%itptn_xm(1,:,i,j,k)*this%TN(1,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(1,:,i,j,k)*this%diff(1,i-1:i,j,k))*sum(this%grdtn_x(1,:,i,j,k)*this%TN(1,i-1:i,j,k))
               ! Fluxes on x-face - 2
               FX(2,i,j,k)=-0.5_WP*(rhoU(2,i,j,k)+abs(rhoU(2,i,j,k)))*sum(this%itptn_xp(2,:,i,j,k)*this%TN(2,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(2,i,j,k)-abs(rhoU(2,i,j,k)))*sum(this%itptn_xm(2,:,i,j,k)*this%TN(2,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(2,:,i,j,k)*this%diff(2,i-1:i,j,k))*sum(this%grdtn_x(2,:,i,j,k)*this%TN(2,i-1:i,j,k))
               ! Fluxes on x-face - 3
               FX(3,i,j,k)=-0.5_WP*(rhoU(3,i,j,k)+abs(rhoU(3,i,j,k)))*sum(this%itptn_xp(3,:,i,j,k)*this%TN(3,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(3,i,j,k)-abs(rhoU(3,i,j,k)))*sum(this%itptn_xm(3,:,i,j,k)*this%TN(3,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(3,:,i,j,k)*this%diff(3,i-1:i,j,k))*sum(this%grdtn_x(3,:,i,j,k)*this%TN(3,i-1:i,j,k))
               ! Fluxes on x-face - 4
               FX(4,i,j,k)=-0.5_WP*(rhoU(4,i,j,k)+abs(rhoU(4,i,j,k)))*sum(this%itptn_xp(4,:,i,j,k)*this%TN(4,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(4,i,j,k)-abs(rhoU(4,i,j,k)))*sum(this%itptn_xm(4,:,i,j,k)*this%TN(4,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(4,:,i,j,k)*this%diff(4,i-1:i,j,k))*sum(this%grdtn_x(4,:,i,j,k)*this%TN(4,i-1:i,j,k))
               ! Fluxes on x-face - 5
               FX(5,i,j,k)=-0.5_WP*(rhoU(5,i,j,k)+abs(rhoU(5,i,j,k)))*sum(this%itptn_xp(5,:,i,j,k)*this%TN(5,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(5,i,j,k)-abs(rhoU(5,i,j,k)))*sum(this%itptn_xm(5,:,i,j,k)*this%TN(5,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(5,:,i,j,k)*this%diff(5,i-1:i,j,k))*sum(this%grdtn_x(5,:,i,j,k)*this%TN(5,i-1:i,j,k))
               ! Fluxes on x-face - 6
               FX(6,i,j,k)=-0.5_WP*(rhoU(6,i,j,k)+abs(rhoU(6,i,j,k)))*sum(this%itptn_xp(6,:,i,j,k)*this%TN(6,i+this%stp1:i+this%stp2,j,k)) &
               &           -0.5_WP*(rhoU(6,i,j,k)-abs(rhoU(6,i,j,k)))*sum(this%itptn_xm(6,:,i,j,k)*this%TN(6,i+this%stm1:i+this%stm2,j,k)) &
               &           +sum(this%itp_x(6,:,i,j,k)*this%diff(6,i-1:i,j,k))*sum(this%grdtn_x(6,:,i,j,k)*this%TN(6,i-1:i,j,k))
               ! Fluxes on y-face - 1
               FY(1,i,j,k)=-0.5_WP*(rhoV(1,i,j,k)+abs(rhoV(1,i,j,k)))*sum(this%itptn_yp(1,:,i,j,k)*this%TN(1,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(1,i,j,k)-abs(rhoV(1,i,j,k)))*sum(this%itptn_ym(1,:,i,j,k)*this%TN(1,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(1,:,i,j,k)*this%diff(1,i,j-1:j,k))*sum(this%grdtn_y(1,:,i,j,k)*this%TN(1,i,j-1:j,k))
               ! Fluxes on y-face - 2
               FY(2,i,j,k)=-0.5_WP*(rhoV(2,i,j,k)+abs(rhoV(2,i,j,k)))*sum(this%itptn_yp(2,:,i,j,k)*this%TN(2,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(2,i,j,k)-abs(rhoV(2,i,j,k)))*sum(this%itptn_ym(2,:,i,j,k)*this%TN(2,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(2,:,i,j,k)*this%diff(2,i,j-1:j,k))*sum(this%grdtn_y(2,:,i,j,k)*this%TN(2,i,j-1:j,k))
               ! Fluxes on y-face - 3
               FY(3,i,j,k)=-0.5_WP*(rhoV(3,i,j,k)+abs(rhoV(3,i,j,k)))*sum(this%itptn_yp(3,:,i,j,k)*this%TN(3,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(3,i,j,k)-abs(rhoV(3,i,j,k)))*sum(this%itptn_ym(3,:,i,j,k)*this%TN(3,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(3,:,i,j,k)*this%diff(3,i,j-1:j,k))*sum(this%grdtn_y(3,:,i,j,k)*this%TN(3,i,j-1:j,k))
               ! Fluxes on y-face - 4
               FY(4,i,j,k)=-0.5_WP*(rhoV(4,i,j,k)+abs(rhoV(4,i,j,k)))*sum(this%itptn_yp(4,:,i,j,k)*this%TN(4,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(4,i,j,k)-abs(rhoV(4,i,j,k)))*sum(this%itptn_ym(4,:,i,j,k)*this%TN(4,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(4,:,i,j,k)*this%diff(4,i,j-1:j,k))*sum(this%grdtn_y(4,:,i,j,k)*this%TN(4,i,j-1:j,k))
               ! Fluxes on y-face - 5
               FY(5,i,j,k)=-0.5_WP*(rhoV(5,i,j,k)+abs(rhoV(5,i,j,k)))*sum(this%itptn_yp(5,:,i,j,k)*this%TN(5,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(5,i,j,k)-abs(rhoV(5,i,j,k)))*sum(this%itptn_ym(5,:,i,j,k)*this%TN(5,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(5,:,i,j,k)*this%diff(5,i,j-1:j,k))*sum(this%grdtn_y(5,:,i,j,k)*this%TN(5,i,j-1:j,k))
               ! Fluxes on y-face - 6
               FY(6,i,j,k)=-0.5_WP*(rhoV(6,i,j,k)+abs(rhoV(6,i,j,k)))*sum(this%itptn_yp(6,:,i,j,k)*this%TN(6,i,j+this%stp1:j+this%stp2,k)) &
               &           -0.5_WP*(rhoV(6,i,j,k)-abs(rhoV(6,i,j,k)))*sum(this%itptn_ym(6,:,i,j,k)*this%TN(6,i,j+this%stm1:j+this%stm2,k)) &
               &           +sum(this%itp_y(6,:,i,j,k)*this%diff(6,i,j-1:j,k))*sum(this%grdtn_y(6,:,i,j,k)*this%TN(6,i,j-1:j,k))
               ! Fluxes on z-face - 1
               FZ(1,i,j,k)=-0.5_WP*(rhoW(1,i,j,k)+abs(rhoW(1,i,j,k)))*sum(this%itptn_zp(1,:,i,j,k)*this%TN(1,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(1,i,j,k)-abs(rhoW(1,i,j,k)))*sum(this%itptn_zm(1,:,i,j,k)*this%TN(1,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(1,:,i,j,k)*this%diff(1,i,j,k-1:k))*sum(this%grdtn_z(1,:,i,j,k)*this%TN(1,i,j,k-1:k))
               ! Fluxes on z-face - 2
               FZ(2,i,j,k)=-0.5_WP*(rhoW(2,i,j,k)+abs(rhoW(2,i,j,k)))*sum(this%itptn_zp(2,:,i,j,k)*this%TN(2,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(2,i,j,k)-abs(rhoW(2,i,j,k)))*sum(this%itptn_zm(2,:,i,j,k)*this%TN(2,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(2,:,i,j,k)*this%diff(2,i,j,k-1:k))*sum(this%grdtn_z(2,:,i,j,k)*this%TN(2,i,j,k-1:k))
               ! Fluxes on z-face - 3
               FZ(3,i,j,k)=-0.5_WP*(rhoW(3,i,j,k)+abs(rhoW(3,i,j,k)))*sum(this%itptn_zp(3,:,i,j,k)*this%TN(3,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(3,i,j,k)-abs(rhoW(3,i,j,k)))*sum(this%itptn_zm(3,:,i,j,k)*this%TN(3,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(3,:,i,j,k)*this%diff(3,i,j,k-1:k))*sum(this%grdtn_z(3,:,i,j,k)*this%TN(3,i,j,k-1:k))
               ! Fluxes on z-face - 4
               FZ(4,i,j,k)=-0.5_WP*(rhoW(4,i,j,k)+abs(rhoW(4,i,j,k)))*sum(this%itptn_zp(4,:,i,j,k)*this%TN(4,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(4,i,j,k)-abs(rhoW(4,i,j,k)))*sum(this%itptn_zm(4,:,i,j,k)*this%TN(4,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(4,:,i,j,k)*this%diff(4,i,j,k-1:k))*sum(this%grdtn_z(4,:,i,j,k)*this%TN(4,i,j,k-1:k))
               ! Fluxes on z-face - 5
               FZ(5,i,j,k)=-0.5_WP*(rhoW(5,i,j,k)+abs(rhoW(5,i,j,k)))*sum(this%itptn_zp(5,:,i,j,k)*this%TN(5,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(5,i,j,k)-abs(rhoW(5,i,j,k)))*sum(this%itptn_zm(5,:,i,j,k)*this%TN(5,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(5,:,i,j,k)*this%diff(5,i,j,k-1:k))*sum(this%grdtn_z(5,:,i,j,k)*this%TN(5,i,j,k-1:k))
               ! Fluxes on z-face - 6
               FZ(6,i,j,k)=-0.5_WP*(rhoW(6,i,j,k)+abs(rhoW(6,i,j,k)))*sum(this%itptn_zp(6,:,i,j,k)*this%TN(6,i,j,k+this%stp1:k+this%stp2)) &
               &           -0.5_WP*(rhoW(6,i,j,k)-abs(rhoW(6,i,j,k)))*sum(this%itptn_zm(6,:,i,j,k)*this%TN(6,i,j,k+this%stm1:k+this%stm2)) &
               &           +sum(this%itp_z(6,:,i,j,k)*this%diff(6,i,j,k-1:k))*sum(this%grdtn_z(6,:,i,j,k)*this%TN(6,i,j,k-1:k))
            end do
         end do
      end do
      ! Time derivative of rhoTN
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               drhoTNdt(1,i,j,k)=sum(this%divtn_x(1,:,i,j,k)*FX(1,i:i+1,j,k))+&
               &                 sum(this%divtn_y(1,:,i,j,k)*FY(1,i,j:j+1,k))+&
               &                 sum(this%divtn_z(1,:,i,j,k)*FZ(1,i,j,k:k+1))
               drhoTNdt(2,i,j,k)=sum(this%divtn_x(2,:,i,j,k)*FX(2,i:i+1,j,k))+&
               &                 sum(this%divtn_y(2,:,i,j,k)*FY(2,i,j:j+1,k))+&
               &                 sum(this%divtn_z(2,:,i,j,k)*FZ(2,i,j,k:k+1))
               drhoTNdt(3,i,j,k)=sum(this%divtn_x(3,:,i,j,k)*FX(3,i:i+1,j,k))+&
               &                 sum(this%divtn_y(3,:,i,j,k)*FY(3,i,j:j+1,k))+&
               &                 sum(this%divtn_z(3,:,i,j,k)*FZ(3,i,j,k:k+1))
               drhoTNdt(4,i,j,k)=sum(this%divtn_x(4,:,i,j,k)*FX(4,i:i+1,j,k))+&
               &                 sum(this%divtn_y(4,:,i,j,k)*FY(4,i,j:j+1,k))+&
               &                 sum(this%divtn_z(4,:,i,j,k)*FZ(4,i,j,k:k+1))
               drhoTNdt(5,i,j,k)=sum(this%divtn_x(5,:,i,j,k)*FX(5,i:i+1,j,k))+&
               &                 sum(this%divtn_y(5,:,i,j,k)*FY(5,i,j:j+1,k))+&
               &                 sum(this%divtn_z(5,:,i,j,k)*FZ(5,i,j,k:k+1))
               drhoTNdt(6,i,j,k)=sum(this%divtn_x(6,:,i,j,k)*FX(6,i:i+1,j,k))+&
               &                 sum(this%divtn_y(6,:,i,j,k)*FY(6,i,j:j+1,k))+&
               &                 sum(this%divtn_z(6,:,i,j,k)*FZ(6,i,j,k:k+1))
            end do
         end do
      end do
      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)
      ! Sync residual
      call this%cfg%sync(drhoTNdt)
   end subroutine get_drhoTNdt
   
   
   !> Calculate the min and max of our TN field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tensor), intent(inout) :: this
      integer :: ierr
      real(WP) :: my_TNmax,my_TNmin
      my_TNmax=maxval(this%TN); call MPI_ALLREDUCE(my_TNmax,this%TNmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_TNmin=minval(this%TN); call MPI_ALLREDUCE(my_TNmin,this%TNmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max
   
   
   ! !> Calculate the integral of our TN field
   ! subroutine get_int(this)
   !    implicit none
   !    class(tensor), intent(inout) :: this
   !    call this%cfg%integrate(this%TN(1,:,:,:),integral=this%TNint(1,:,:,:))
   ! end subroutine get_int
   
   
   !> Solve for implicit tensor residual
   subroutine solve_implicit(this,dt,resTN,rhoU,rhoV,rhoW)
      implicit none
      class(tensor), intent(inout) :: this
      real(WP), intent(in) :: dt
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: resTN !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoU  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoV  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: rhoW  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,n,sti,std
      
      ! Prescribe density
      this%implicit%opr(1,:,:,:)=this%rho; this%implicit%opr(2:,:,:,:)=0.0_WP

      !> Tensor component 1
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(1,std,i,j,k)*0.5_WP*(rhoU(1,i+std,j,k)+abs(rhoU(1,i+std,j,k)))*this%itptn_xp(1,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(1,std,i,j,k)*0.5_WP*(rhoV(1,i,j+std,k)+abs(rhoV(1,i,j+std,k)))*this%itptn_yp(1,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(1,std,i,j,k)*0.5_WP*(rhoW(1,i,j,k+std)+abs(rhoW(1,i,j,k+std)))*this%itptn_zp(1,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(1,std,i,j,k)*0.5_WP*(rhoU(1,i+std,j,k)-abs(rhoU(1,i+std,j,k)))*this%itptn_xm(1,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(1,std,i,j,k)*0.5_WP*(rhoV(1,i,j+std,k)-abs(rhoV(1,i,j+std,k)))*this%itptn_ym(1,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(1,std,i,j,k)*0.5_WP*(rhoW(1,i,j,k+std)-abs(rhoW(1,i,j,k+std)))*this%itptn_zm(1,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(1,+1,i,j,k)*sum(this%itp_x(1,:,i+1,j,k)*this%diff(1,i  :i+1,j,k))*this%grdtn_x(1,-1,i+1,j,k)+&
               &                                                                this%divtn_x(1, 0,i,j,k)*sum(this%itp_x(1,:,i  ,j,k)*this%diff(1,i-1:i  ,j,k))*this%grdtn_x(1, 0,i  ,j,k)+&
               &                                                                this%divtn_y(1,+1,i,j,k)*sum(this%itp_y(1,:,i,j+1,k)*this%diff(1,i,j  :j+1,k))*this%grdtn_y(1,-1,i,j+1,k)+&
               &                                                                this%divtn_y(1, 0,i,j,k)*sum(this%itp_y(1,:,i,j  ,k)*this%diff(1,i,j-1:j  ,k))*this%grdtn_y(1, 0,i,j  ,k)+&
               &                                                                this%divtn_z(1,+1,i,j,k)*sum(this%itp_z(1,:,i,j,k+1)*this%diff(1,i,j,k  :k+1))*this%grdtn_z(1,-1,i,j,k+1)+&
               &                                                                this%divtn_z(1, 0,i,j,k)*sum(this%itp_z(1,:,i,j,k  )*this%diff(1,i,j,k-1:k  ))*this%grdtn_z(1, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(1,+1,i,j,k)*sum(this%itp_x(1,:,i+1,j,k)*this%diff(1,i  :i+1,j,k))*this%grdtn_x(1, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(1, 0,i,j,k)*sum(this%itp_x(1,:,i  ,j,k)*this%diff(1,i-1:i  ,j,k))*this%grdtn_x(1,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(1,+1,i,j,k)*sum(this%itp_y(1,:,i,j+1,k)*this%diff(1,i,j  :j+1,k))*this%grdtn_y(1, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(1, 0,i,j,k)*sum(this%itp_y(1,:,i,j  ,k)*this%diff(1,i,j-1:j  ,k))*this%grdtn_y(1,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(1,+1,i,j,k)*sum(this%itp_z(1,:,i,j,k+1)*this%diff(1,i,j,k  :k+1))*this%grdtn_z(1, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(1, 0,i,j,k)*sum(this%itp_z(1,:,i,j,k  )*this%diff(1,i,j,k-1:k  ))*this%grdtn_z(1,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(1,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(1,:,:,:)=this%implicit%sol

      !> Tensor component 2
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(2,std,i,j,k)*0.5_WP*(rhoU(2,i+std,j,k)+abs(rhoU(2,i+std,j,k)))*this%itptn_xp(2,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(2,std,i,j,k)*0.5_WP*(rhoV(2,i,j+std,k)+abs(rhoV(2,i,j+std,k)))*this%itptn_yp(2,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(2,std,i,j,k)*0.5_WP*(rhoW(2,i,j,k+std)+abs(rhoW(2,i,j,k+std)))*this%itptn_zp(2,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(2,std,i,j,k)*0.5_WP*(rhoU(2,i+std,j,k)-abs(rhoU(2,i+std,j,k)))*this%itptn_xm(2,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(2,std,i,j,k)*0.5_WP*(rhoV(2,i,j+std,k)-abs(rhoV(2,i,j+std,k)))*this%itptn_ym(2,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(2,std,i,j,k)*0.5_WP*(rhoW(2,i,j,k+std)-abs(rhoW(2,i,j,k+std)))*this%itptn_zm(2,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(2,+1,i,j,k)*sum(this%itp_x(2,:,i+1,j,k)*this%diff(2,i  :i+1,j,k))*this%grdtn_x(2,-1,i+1,j,k)+&
               &                                                                this%divtn_x(2, 0,i,j,k)*sum(this%itp_x(2,:,i  ,j,k)*this%diff(2,i-1:i  ,j,k))*this%grdtn_x(2, 0,i  ,j,k)+&
               &                                                                this%divtn_y(2,+1,i,j,k)*sum(this%itp_y(2,:,i,j+1,k)*this%diff(2,i,j  :j+1,k))*this%grdtn_y(2,-1,i,j+1,k)+&
               &                                                                this%divtn_y(2, 0,i,j,k)*sum(this%itp_y(2,:,i,j  ,k)*this%diff(2,i,j-1:j  ,k))*this%grdtn_y(2, 0,i,j  ,k)+&
               &                                                                this%divtn_z(2,+1,i,j,k)*sum(this%itp_z(2,:,i,j,k+1)*this%diff(2,i,j,k  :k+1))*this%grdtn_z(2,-1,i,j,k+1)+&
               &                                                                this%divtn_z(2, 0,i,j,k)*sum(this%itp_z(2,:,i,j,k  )*this%diff(2,i,j,k-1:k  ))*this%grdtn_z(2, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(2,+1,i,j,k)*sum(this%itp_x(2,:,i+1,j,k)*this%diff(2,i  :i+1,j,k))*this%grdtn_x(2, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(2, 0,i,j,k)*sum(this%itp_x(2,:,i  ,j,k)*this%diff(2,i-1:i  ,j,k))*this%grdtn_x(2,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(2,+1,i,j,k)*sum(this%itp_y(2,:,i,j+1,k)*this%diff(2,i,j  :j+1,k))*this%grdtn_y(2, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(2, 0,i,j,k)*sum(this%itp_y(2,:,i,j  ,k)*this%diff(2,i,j-1:j  ,k))*this%grdtn_y(2,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(2,+1,i,j,k)*sum(this%itp_z(2,:,i,j,k+1)*this%diff(2,i,j,k  :k+1))*this%grdtn_z(2, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(2, 0,i,j,k)*sum(this%itp_z(2,:,i,j,k  )*this%diff(2,i,j,k-1:k  ))*this%grdtn_z(2,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(2,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(2,:,:,:)=this%implicit%sol

      !> Tensor component 3
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(3,std,i,j,k)*0.5_WP*(rhoU(3,i+std,j,k)+abs(rhoU(3,i+std,j,k)))*this%itptn_xp(3,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(3,std,i,j,k)*0.5_WP*(rhoV(3,i,j+std,k)+abs(rhoV(3,i,j+std,k)))*this%itptn_yp(3,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(3,std,i,j,k)*0.5_WP*(rhoW(3,i,j,k+std)+abs(rhoW(3,i,j,k+std)))*this%itptn_zp(3,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(3,std,i,j,k)*0.5_WP*(rhoU(3,i+std,j,k)-abs(rhoU(3,i+std,j,k)))*this%itptn_xm(3,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(3,std,i,j,k)*0.5_WP*(rhoV(3,i,j+std,k)-abs(rhoV(3,i,j+std,k)))*this%itptn_ym(3,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(3,std,i,j,k)*0.5_WP*(rhoW(3,i,j,k+std)-abs(rhoW(3,i,j,k+std)))*this%itptn_zm(3,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(3,+1,i,j,k)*sum(this%itp_x(3,:,i+1,j,k)*this%diff(3,i  :i+1,j,k))*this%grdtn_x(3,-1,i+1,j,k)+&
               &                                                                this%divtn_x(3, 0,i,j,k)*sum(this%itp_x(3,:,i  ,j,k)*this%diff(3,i-1:i  ,j,k))*this%grdtn_x(3, 0,i  ,j,k)+&
               &                                                                this%divtn_y(3,+1,i,j,k)*sum(this%itp_y(3,:,i,j+1,k)*this%diff(3,i,j  :j+1,k))*this%grdtn_y(3,-1,i,j+1,k)+&
               &                                                                this%divtn_y(3, 0,i,j,k)*sum(this%itp_y(3,:,i,j  ,k)*this%diff(3,i,j-1:j  ,k))*this%grdtn_y(3, 0,i,j  ,k)+&
               &                                                                this%divtn_z(3,+1,i,j,k)*sum(this%itp_z(3,:,i,j,k+1)*this%diff(3,i,j,k  :k+1))*this%grdtn_z(3,-1,i,j,k+1)+&
               &                                                                this%divtn_z(3, 0,i,j,k)*sum(this%itp_z(3,:,i,j,k  )*this%diff(3,i,j,k-1:k  ))*this%grdtn_z(3, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(3,+1,i,j,k)*sum(this%itp_x(3,:,i+1,j,k)*this%diff(3,i  :i+1,j,k))*this%grdtn_x(3, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(3, 0,i,j,k)*sum(this%itp_x(3,:,i  ,j,k)*this%diff(3,i-1:i  ,j,k))*this%grdtn_x(3,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(3,+1,i,j,k)*sum(this%itp_y(3,:,i,j+1,k)*this%diff(3,i,j  :j+1,k))*this%grdtn_y(3, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(3, 0,i,j,k)*sum(this%itp_y(3,:,i,j  ,k)*this%diff(3,i,j-1:j  ,k))*this%grdtn_y(3,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(3,+1,i,j,k)*sum(this%itp_z(3,:,i,j,k+1)*this%diff(3,i,j,k  :k+1))*this%grdtn_z(3, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(3, 0,i,j,k)*sum(this%itp_z(3,:,i,j,k  )*this%diff(3,i,j,k-1:k  ))*this%grdtn_z(3,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(3,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(3,:,:,:)=this%implicit%sol

      !> Tensor component 4
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(4,std,i,j,k)*0.5_WP*(rhoU(4,i+std,j,k)+abs(rhoU(4,i+std,j,k)))*this%itptn_xp(4,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(4,std,i,j,k)*0.5_WP*(rhoV(4,i,j+std,k)+abs(rhoV(4,i,j+std,k)))*this%itptn_yp(4,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(4,std,i,j,k)*0.5_WP*(rhoW(4,i,j,k+std)+abs(rhoW(4,i,j,k+std)))*this%itptn_zp(4,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(4,std,i,j,k)*0.5_WP*(rhoU(4,i+std,j,k)-abs(rhoU(4,i+std,j,k)))*this%itptn_xm(4,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(4,std,i,j,k)*0.5_WP*(rhoV(4,i,j+std,k)-abs(rhoV(4,i,j+std,k)))*this%itptn_ym(4,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(4,std,i,j,k)*0.5_WP*(rhoW(4,i,j,k+std)-abs(rhoW(4,i,j,k+std)))*this%itptn_zm(4,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(4,+1,i,j,k)*sum(this%itp_x(4,:,i+1,j,k)*this%diff(4,i  :i+1,j,k))*this%grdtn_x(4,-1,i+1,j,k)+&
               &                                                                this%divtn_x(4, 0,i,j,k)*sum(this%itp_x(4,:,i  ,j,k)*this%diff(4,i-1:i  ,j,k))*this%grdtn_x(4, 0,i  ,j,k)+&
               &                                                                this%divtn_y(4,+1,i,j,k)*sum(this%itp_y(4,:,i,j+1,k)*this%diff(4,i,j  :j+1,k))*this%grdtn_y(4,-1,i,j+1,k)+&
               &                                                                this%divtn_y(4, 0,i,j,k)*sum(this%itp_y(4,:,i,j  ,k)*this%diff(4,i,j-1:j  ,k))*this%grdtn_y(4, 0,i,j  ,k)+&
               &                                                                this%divtn_z(4,+1,i,j,k)*sum(this%itp_z(4,:,i,j,k+1)*this%diff(4,i,j,k  :k+1))*this%grdtn_z(4,-1,i,j,k+1)+&
               &                                                                this%divtn_z(4, 0,i,j,k)*sum(this%itp_z(4,:,i,j,k  )*this%diff(4,i,j,k-1:k  ))*this%grdtn_z(4, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(4,+1,i,j,k)*sum(this%itp_x(4,:,i+1,j,k)*this%diff(4,i  :i+1,j,k))*this%grdtn_x(4, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(4, 0,i,j,k)*sum(this%itp_x(4,:,i  ,j,k)*this%diff(4,i-1:i  ,j,k))*this%grdtn_x(4,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(4,+1,i,j,k)*sum(this%itp_y(4,:,i,j+1,k)*this%diff(4,i,j  :j+1,k))*this%grdtn_y(4, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(4, 0,i,j,k)*sum(this%itp_y(4,:,i,j  ,k)*this%diff(4,i,j-1:j  ,k))*this%grdtn_y(4,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(4,+1,i,j,k)*sum(this%itp_z(4,:,i,j,k+1)*this%diff(4,i,j,k  :k+1))*this%grdtn_z(4, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(4, 0,i,j,k)*sum(this%itp_z(4,:,i,j,k  )*this%diff(4,i,j,k-1:k  ))*this%grdtn_z(4,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(4,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(4,:,:,:)=this%implicit%sol

      !> Tensor component 5
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(5,std,i,j,k)*0.5_WP*(rhoU(5,i+std,j,k)+abs(rhoU(5,i+std,j,k)))*this%itptn_xp(5,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(5,std,i,j,k)*0.5_WP*(rhoV(5,i,j+std,k)+abs(rhoV(5,i,j+std,k)))*this%itptn_yp(5,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(5,std,i,j,k)*0.5_WP*(rhoW(5,i,j,k+std)+abs(rhoW(5,i,j,k+std)))*this%itptn_zp(5,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(5,std,i,j,k)*0.5_WP*(rhoU(5,i+std,j,k)-abs(rhoU(5,i+std,j,k)))*this%itptn_xm(5,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(5,std,i,j,k)*0.5_WP*(rhoV(5,i,j+std,k)-abs(rhoV(5,i,j+std,k)))*this%itptn_ym(5,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(5,std,i,j,k)*0.5_WP*(rhoW(5,i,j,k+std)-abs(rhoW(5,i,j,k+std)))*this%itptn_zm(5,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(5,+1,i,j,k)*sum(this%itp_x(5,:,i+1,j,k)*this%diff(5,i  :i+1,j,k))*this%grdtn_x(5,-1,i+1,j,k)+&
               &                                                                this%divtn_x(5, 0,i,j,k)*sum(this%itp_x(5,:,i  ,j,k)*this%diff(5,i-1:i  ,j,k))*this%grdtn_x(5, 0,i  ,j,k)+&
               &                                                                this%divtn_y(5,+1,i,j,k)*sum(this%itp_y(5,:,i,j+1,k)*this%diff(5,i,j  :j+1,k))*this%grdtn_y(5,-1,i,j+1,k)+&
               &                                                                this%divtn_y(5, 0,i,j,k)*sum(this%itp_y(5,:,i,j  ,k)*this%diff(5,i,j-1:j  ,k))*this%grdtn_y(5, 0,i,j  ,k)+&
               &                                                                this%divtn_z(5,+1,i,j,k)*sum(this%itp_z(5,:,i,j,k+1)*this%diff(5,i,j,k  :k+1))*this%grdtn_z(5,-1,i,j,k+1)+&
               &                                                                this%divtn_z(5, 0,i,j,k)*sum(this%itp_z(5,:,i,j,k  )*this%diff(5,i,j,k-1:k  ))*this%grdtn_z(5, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(5,+1,i,j,k)*sum(this%itp_x(5,:,i+1,j,k)*this%diff(5,i  :i+1,j,k))*this%grdtn_x(5, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(5, 0,i,j,k)*sum(this%itp_x(5,:,i  ,j,k)*this%diff(5,i-1:i  ,j,k))*this%grdtn_x(5,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(5,+1,i,j,k)*sum(this%itp_y(5,:,i,j+1,k)*this%diff(5,i,j  :j+1,k))*this%grdtn_y(5, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(5, 0,i,j,k)*sum(this%itp_y(5,:,i,j  ,k)*this%diff(5,i,j-1:j  ,k))*this%grdtn_y(5,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(5,+1,i,j,k)*sum(this%itp_z(5,:,i,j,k+1)*this%diff(5,i,j,k  :k+1))*this%grdtn_z(5, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(5, 0,i,j,k)*sum(this%itp_z(5,:,i,j,k  )*this%diff(5,i,j,k-1:k  ))*this%grdtn_z(5,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(5,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(5,:,:,:)=this%implicit%sol

      !> Tensor component 6
      ! Prepare convective operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Loop over divergence stencil
               do std=0,1
                  ! Loop over plus interpolation stencil
                  do sti=this%stp1,this%stp2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(6,std,i,j,k)*0.5_WP*(rhoU(6,i+std,j,k)+abs(rhoU(6,i+std,j,k)))*this%itptn_xp(6,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(6,std,i,j,k)*0.5_WP*(rhoV(6,i,j+std,k)+abs(rhoV(6,i,j+std,k)))*this%itptn_yp(6,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(6,std,i,j,k)*0.5_WP*(rhoW(6,i,j,k+std)+abs(rhoW(6,i,j,k+std)))*this%itptn_zp(6,sti,i,j,k+std)
                  end do
                  ! Loop over minus interpolation stencil
                  do sti=this%stm1,this%stm2
                     this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)=this%implicit%opr(this%implicit%stmap(sti+std,0,0),i,j,k)+0.5_WP*dt*this%divtn_x(6,std,i,j,k)*0.5_WP*(rhoU(6,i+std,j,k)-abs(rhoU(6,i+std,j,k)))*this%itptn_xm(6,sti,i+std,j,k)
                     this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)=this%implicit%opr(this%implicit%stmap(0,sti+std,0),i,j,k)+0.5_WP*dt*this%divtn_y(6,std,i,j,k)*0.5_WP*(rhoV(6,i,j+std,k)-abs(rhoV(6,i,j+std,k)))*this%itptn_ym(6,sti,i,j+std,k)
                     this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)=this%implicit%opr(this%implicit%stmap(0,0,sti+std),i,j,k)+0.5_WP*dt*this%divtn_z(6,std,i,j,k)*0.5_WP*(rhoW(6,i,j,k+std)-abs(rhoW(6,i,j,k+std)))*this%itptn_zm(6,sti,i,j,k+std)
                  end do
               end do
            end do
         end do
      end do
         
      ! Prepare diffusive operator
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%opr(1,i,j,k)=this%implicit%opr(1,i,j,k)-0.5_WP*dt*(this%divtn_x(6,+1,i,j,k)*sum(this%itp_x(6,:,i+1,j,k)*this%diff(6,i  :i+1,j,k))*this%grdtn_x(6,-1,i+1,j,k)+&
               &                                                                this%divtn_x(6, 0,i,j,k)*sum(this%itp_x(6,:,i  ,j,k)*this%diff(6,i-1:i  ,j,k))*this%grdtn_x(6, 0,i  ,j,k)+&
               &                                                                this%divtn_y(6,+1,i,j,k)*sum(this%itp_y(6,:,i,j+1,k)*this%diff(6,i,j  :j+1,k))*this%grdtn_y(6,-1,i,j+1,k)+&
               &                                                                this%divtn_y(6, 0,i,j,k)*sum(this%itp_y(6,:,i,j  ,k)*this%diff(6,i,j-1:j  ,k))*this%grdtn_y(6, 0,i,j  ,k)+&
               &                                                                this%divtn_z(6,+1,i,j,k)*sum(this%itp_z(6,:,i,j,k+1)*this%diff(6,i,j,k  :k+1))*this%grdtn_z(6,-1,i,j,k+1)+&
               &                                                                this%divtn_z(6, 0,i,j,k)*sum(this%itp_z(6,:,i,j,k  )*this%diff(6,i,j,k-1:k  ))*this%grdtn_z(6, 0,i,j,k  ))
               this%implicit%opr(2,i,j,k)=this%implicit%opr(2,i,j,k)-0.5_WP*dt*(this%divtn_x(6,+1,i,j,k)*sum(this%itp_x(6,:,i+1,j,k)*this%diff(6,i  :i+1,j,k))*this%grdtn_x(6, 0,i+1,j,k))
               this%implicit%opr(3,i,j,k)=this%implicit%opr(3,i,j,k)-0.5_WP*dt*(this%divtn_x(6, 0,i,j,k)*sum(this%itp_x(6,:,i  ,j,k)*this%diff(6,i-1:i  ,j,k))*this%grdtn_x(6,-1,i  ,j,k))
               this%implicit%opr(4,i,j,k)=this%implicit%opr(4,i,j,k)-0.5_WP*dt*(this%divtn_y(6,+1,i,j,k)*sum(this%itp_y(6,:,i,j+1,k)*this%diff(6,i,j  :j+1,k))*this%grdtn_y(6, 0,i,j+1,k))
               this%implicit%opr(5,i,j,k)=this%implicit%opr(5,i,j,k)-0.5_WP*dt*(this%divtn_y(6, 0,i,j,k)*sum(this%itp_y(6,:,i,j  ,k)*this%diff(6,i,j-1:j  ,k))*this%grdtn_y(6,-1,i,j  ,k))
               this%implicit%opr(6,i,j,k)=this%implicit%opr(6,i,j,k)-0.5_WP*dt*(this%divtn_z(6,+1,i,j,k)*sum(this%itp_z(6,:,i,j,k+1)*this%diff(6,i,j,k  :k+1))*this%grdtn_z(6, 0,i,j,k+1))
               this%implicit%opr(7,i,j,k)=this%implicit%opr(7,i,j,k)-0.5_WP*dt*(this%divtn_z(6, 0,i,j,k)*sum(this%itp_z(6,:,i,j,k  )*this%diff(6,i,j,k-1:k  ))*this%grdtn_z(6,-1,i,j,k  ))
            end do
         end do
      end do
         
      ! Solve the linear system
      call this%implicit%setup()
      this%implicit%rhs=resTN(6,:,:,:)
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      resTN(6,:,:,:)=this%implicit%sol
      
      ! Sync up residual
      call this%cfg%sync(resTN)
      
   end subroutine solve_implicit
   
   
   !> Print out info for tensor solver
   subroutine tensor_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(tensor), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Constant density tensor solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine tensor_print
   
   
end module tensor_class
