!> FENE model solver class:
!> Provides support for various BC, RHS calculation, implicit solver
!> Assumes constant diffusivity and density.
module fene_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use iterator_class, only: iterator
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: fene
   
   !> Constant density fene solver object definition
   type :: fene
      
      ! This is our config
      class(config), pointer :: cfg                         !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_FENE'      !< Solver name (default=UNNAMED_FENE)

      ! Source term arrays
      real(WP), dimension(:,:,:,:), allocatable :: CgradU  !< Sum of CdotU and (CdotU)T
      real(WP), dimension(:,:,:,:), allocatable :: T       !< Viscoelastic stress tensor

      ! Viscoleastic tensor divergence
      real(WP), dimension(:,:,:,:), allocatable :: divT     !< Divergence array

      ! Number of elements in conformation tensor
      integer :: Celem=6                                   !< Symmetric tensor with 6 unique elements
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z     !< Divergence for viscoelatic stress tensor 
      
      ! Masking info for metric modification
      integer, dimension(:,:,:), allocatable :: mask      !< Integer array used for modifying C metrics
      
      
   contains
      procedure :: print=>fene_print                      !< Output solver to the screen
      procedure :: setup                                  !< Finish configuring the fene solver
      procedure :: init_metrics                           !< Initialize metrics
      procedure :: adjust_metrics                         !< Adjust metrics
      procedure :: get_CgradU                             !< Calculate product and transpose of C_dot_gradU
      procedure :: get_stressTensor                       !< Calculate the viscoelastic stress tensor
      procedure :: get_divT                               !< Calculate viscoelastic tensor divergence
   end type fene
   
   
   !> Declare fene solver constructor
   interface fene
      procedure constructor
   end interface fene
   
contains
   
   
   !> Default constructor for fene solver
   function constructor(cfg,name) result(self)
      implicit none
      type(fene) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Allocate variables
      allocate(self%CgradU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%Celem)); self%CgradU=0.0_WP
      allocate(self%T     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%Celem)); self%T     =0.0_WP
      allocate(self%divT  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3         )); self%divT  =0.0_WP
      
      ! Prepare default metrics
      call self%init_metrics()
      
      ! Prepare mask for C
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
      implicit none
      class(fene), intent(inout) :: this
      integer :: i,j,k
      
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

   end subroutine init_metrics
   
   
   !> Metric adjustment accounting for bconds and walls - zero out div at bcond and walls
   subroutine adjust_metrics(this)
      implicit none
      class(fene), intent(inout) :: this
      integer :: i,j,k
      
      ! Sync up masks
      call this%cfg%sync(this%mask)
      
      ! Loop over the domain and apply masked conditions to T divergence
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
      
      ! Adjust metrics to account for lower dimensionality
      if (this%cfg%nx.eq.1) then
         this%div_x=0.0_WP
      end if
      if (this%cfg%ny.eq.1) then
         this%div_y=0.0_WP
      end if
      if (this%cfg%nz.eq.1) then
         this%div_z=0.0_WP
      end if
      
   end subroutine adjust_metrics
   
   
   !> Finish setting up the fene solver now that bconds have been defined
   subroutine setup(this)
      implicit none
      class(fene), intent(inout) :: this
      
      ! Adjust metrics based on mask array
      call this%adjust_metrics()
      
   end subroutine setup

   !> Calculate components of tensor (c*graduT)+(C*gradu)^T
   subroutine get_CgradU(this,C,gradu)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:   ), intent(in) :: C
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradu
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! 11 tensor component
               this%CgradU(i,j,k,1)=2.00_WP*(C(i,j,k,1)*gradu(1,1,i,j,k)+C(i,j,k,2)*gradu(1,2,i,j,k)+C(i,j,k,3)*gradu(1,3,i,j,k))
               ! 21/12 tensor component
               this%CgradU(i,j,k,2)=(C(i,j,k,2)*gradu(1,1,i,j,k)+C(i,j,k,4)*gradu(1,2,i,j,k)+C(i,j,k,5)*gradu(1,3,i,j,k))+&
               &                    (C(i,j,k,1)*gradu(2,1,i,j,k)+C(i,j,k,2)*gradu(2,2,i,j,k)+C(i,j,k,3)*gradu(2,3,i,j,k))
               ! 31/13 tensor component
               this%CgradU(i,j,k,3)=(C(i,j,k,3)*gradu(1,1,i,j,k)+C(i,j,k,5)*gradu(1,2,i,j,k)+C(i,j,k,6)*gradu(1,3,i,j,k))+&
               &                    (C(i,j,k,1)*gradu(3,1,i,j,k)+C(i,j,k,2)*gradu(3,2,i,j,k)+C(i,j,k,3)*gradu(3,3,i,j,k))
               ! 22 tensor component
               this%CgradU(i,j,k,4)=2.00_WP*(C(i,j,k,2)*gradu(2,1,i,j,k)+C(i,j,k,4)*gradu(2,2,i,j,k)+C(i,j,k,5)*gradu(2,3,i,j,k))
               ! 32/23 tensor component
               this%CgradU(i,j,k,5)=(C(i,j,k,2)*gradu(3,1,i,j,k)+C(i,j,k,4)*gradu(3,2,i,j,k)+C(i,j,k,5)*gradu(3,3,i,j,k))+&
               &                    (C(i,j,k,3)*gradu(2,1,i,j,k)+C(i,j,k,5)*gradu(2,2,i,j,k)+C(i,j,k,6)*gradu(2,3,i,j,k))
               ! 33 tensor component
               this%CgradU(i,j,k,6)=2.00_WP*(C(i,j,k,3)*gradu(3,1,i,j,k)+C(i,j,k,5)*gradu(3,2,i,j,k)+C(i,j,k,6)*gradu(3,3,i,j,k))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%pgrid_rsync_right_array(this%CgradU)
   end subroutine get_CgradU

   !> Calculate the viscoelastic stress tensor
   subroutine get_stressTensor(this,C,We,Lmax)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(in) :: C
      real(WP), intent(in) :: Lmax,We
      real(WP) :: a,psi
      integer :: i,j,k
      ! a parameter
      a=1.00_WP-3.00_WP/(Lmax**2)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! psi parameter
               psi=1.00_WP-(C(i,j,k,1)+C(i,j,k,2)+C(i,j,k,6))/(Lmax**2)
               ! 11 tensor component
               this%T(i,j,k,1)=(1.00_WP/We)*C(i,j,k,1)/psi-1.00_WP/a
               ! 21/12 tensor component
               this%T(i,j,k,2)=(1.00_WP/We)*C(i,j,k,2)/psi-0.00_WP/a
               ! 31/13 tensor component
               this%T(i,j,k,3)=(1.00_WP/We)*C(i,j,k,3)/psi-0.00_WP/a
               ! 22 tensor component
               this%T(i,j,k,4)=(1.00_WP/We)*C(i,j,k,4)/psi-1.00_WP/a
               ! 32/23 tensor component
               this%T(i,j,k,5)=(1.00_WP/We)*C(i,j,k,5)/psi-0.00_WP/a
               ! 33 tensor component
               this%T(i,j,k,6)=(1.00_WP/We)*C(i,j,k,6)/psi-1.00_WP/a 
            end do
         end do
      end do
      ! Sync it
      call this%cfg%pgrid_rsync_right_array(this%T)
   end subroutine get_stressTensor
   

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
               this%divT(i,j,k,2)=sum(this%div_x(:,i,j,k)*this%T(i:i+1,j,k,2))+&
               &                  sum(this%div_y(:,i,j,k)*this%T(i,j:j+1,k,4))+&
               &                  sum(this%div_z(:,i,j,k)*this%T(i,j,k:k+1,5))
               this%divT(i,j,k,3)=sum(this%div_x(:,i,j,k)*this%T(i:i+1,j,k,3))+&
               &                  sum(this%div_y(:,i,j,k)*this%T(i,j:j+1,k,5))+&
               &                  sum(this%div_z(:,i,j,k)*this%T(i,j,k:k+1,6))
            end do
         end do
      end do
      ! Sync it
      call this%cfg%pgrid_rsync_right_array(this%divT)
   end subroutine get_divT
   
   
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
