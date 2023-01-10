!> FENE model class
!> Extends multiscalar class for the calculation of source terms
!> Assumes constant diffusivity and density.
module fene_class
   use multiscalar_class, only: multiscalar
   use config_class,      only: config
   use precision,         only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: fene
   
   !> Constant density fene solver object definition
   type, extends(multiscalar) :: fene

      ! Source term arrays
      real(WP), dimension(:,:,:,:), allocatable :: CgradU  !< Sum of distortion terms (CdotU and (CdotU)T)
      real(WP), dimension(:,:,:,:), allocatable :: T       !< Stress tensor
      real(WP), dimension(:,:,:,:), allocatable :: divT    !< Stress tensor divergence

   contains
      procedure :: get_CgradU                             !< Calculate product and transpose of C_dot_gradU
      procedure :: get_stressTensor                       !< Calculate the stress tensor
      procedure :: get_divT                               !< Calculate stress tensor divergence
   end type fene
   
   !> Declare fene model constructor
   interface fene
      procedure construct_fene_from_args
   end interface fene
   
contains
   
   
   !> FENE model constructor from multiscalar
   function construct_fene_from_args(cfg,scheme,name) result(self)
      implicit none
      type(fene) :: self                        !< FENE model
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: scheme
      character(len=*), optional :: name

      ! Create a six-scalar solver
      self%multiscalar=multiscalar(cfg=cfg,scheme=scheme,nscalar=6,name=name)

      ! Allocate variables
      allocate(self%CgradU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%nscalar)); self%CgradU=0.0_WP
      allocate(self%T     (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,self%nscalar)); self%T     =0.0_WP
      allocate(self%divT  (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,3           )); self%divT  =0.0_WP
   end function construct_fene_from_args

   !> Calculate components of tensor (c*graduT)+(C*gradu)^T
   subroutine get_CgradU(this,gradu)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradu
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! xx tensor component
               this%CgradU(i,j,k,1)=2.00_WP*(this%SC(i,j,k,1)*gradu(1,1,i,j,k)+this%SC(i,j,k,2)*gradu(2,1,i,j,k)+this%SC(i,j,k,3)*gradu(3,1,i,j,k))
               ! yx/xy tensor component
               this%CgradU(i,j,k,2)=(this%SC(i,j,k,2)*gradu(1,1,i,j,k)+this%SC(i,j,k,4)*gradu(2,1,i,j,k)+this%SC(i,j,k,5)*gradu(3,1,i,j,k))+&
               &                    (this%SC(i,j,k,1)*gradu(1,2,i,j,k)+this%SC(i,j,k,2)*gradu(2,2,i,j,k)+this%SC(i,j,k,3)*gradu(3,2,i,j,k))
               ! zx/xz tensor component
               this%CgradU(i,j,k,3)=(this%SC(i,j,k,3)*gradu(1,1,i,j,k)+this%SC(i,j,k,5)*gradu(2,1,i,j,k)+this%SC(i,j,k,6)*gradu(3,1,i,j,k))+&
               &                    (this%SC(i,j,k,1)*gradu(1,3,i,j,k)+this%SC(i,j,k,2)*gradu(2,3,i,j,k)+this%SC(i,j,k,3)*gradu(3,3,i,j,k))
               ! yy tensor component
               this%CgradU(i,j,k,4)=2.00_WP*(this%SC(i,j,k,2)*gradu(1,2,i,j,k)+this%SC(i,j,k,4)*gradu(2,2,i,j,k)+this%SC(i,j,k,5)*gradu(3,2,i,j,k))
               ! zy/yz tensor component
               this%CgradU(i,j,k,5)=(this%SC(i,j,k,2)*gradu(1,3,i,j,k)+this%SC(i,j,k,4)*gradu(2,3,i,j,k)+this%SC(i,j,k,5)*gradu(3,3,i,j,k))+&
               &                    (this%SC(i,j,k,3)*gradu(1,2,i,j,k)+this%SC(i,j,k,5)*gradu(2,2,i,j,k)+this%SC(i,j,k,6)*gradu(3,2,i,j,k))
               ! zz tensor component
               this%CgradU(i,j,k,6)=2.00_WP*(this%SC(i,j,k,3)*gradu(1,3,i,j,k)+this%SC(i,j,k,5)*gradu(2,3,i,j,k)+this%SC(i,j,k,6)*gradu(3,3,i,j,k))
            end do
         end do
      end do
   end subroutine get_CgradU

   !> Calculate the viscoelastic stress tensor
   subroutine get_stressTensor(this,lambda,Lmax,visc_p)
      implicit none
      class(fene), intent(inout) :: this
      real(WP), intent(in)       :: Lmax,lambda,visc_p
      real(WP), dimension(:,:,:), allocatable :: trC,f_r        
      integer  :: i,j,k
      
      ! Build stress tensor
      trC=this%SC(:,:,:,1)+this%SC(:,:,:,4)+this%SC(:,:,:,6)   !< trace of C
      f_r=(Lmax**2.00_WP-3.00_WP)/(Lmax**2.00_WP-trC)          !< Peterlin Function
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! xx tensor component
               this%T(i,j,k,1)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,1)-1.00_WP)
               ! yx/xy tensor component
               this%T(i,j,k,2)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,2)-0.00_WP)
               ! zx/xz tensor component
               this%T(i,j,k,3)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,3)-0.00_WP)
               ! yy tensor component
               this%T(i,j,k,4)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,4)-1.00_WP)
               ! zy/yz tensor component
               this%T(i,j,k,5)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,5)-0.00_WP)
               ! zz tensor component
               this%T(i,j,k,6)=(visc_p/lambda)*(f_r(i,j,k)*this%SC(i,j,k,6)-1.00_WP)
            end do
         end do
      end do

      ! Deallocate trace(C) array
	   deallocate(trC,f_r)

   end subroutine get_stressTensor
   

   !> Calculate the viscoelastic tensor divergence
   subroutine get_divT(this,fs)
      use incomp_class, only: incomp
      implicit none
      class(fene), intent(inout) :: this
      class(incomp), intent(in)  :: fs
      integer :: i,j,k
      real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx

      ! Allocate tensor components
	   allocate(Txy(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(Tyz(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(Tzx(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      
      ! Interpolate tensor components to cell faces
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               Txy(i,j,k)=sum(fs%itp_xy(:,:,i,j,k)*this%T(i-1:i,j-1:j,k,2))
               Tyz(i,j,k)=sum(fs%itp_yz(:,:,i,j,k)*this%T(i,j-1:j,k-1:k,5))
               Tzx(i,j,k)=sum(fs%itp_xz(:,:,i,j,k)*this%T(i-1:i,j,k-1:k,3))
            end do
         end do
      end do

      ! Compute divergence
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%divT(i,j,k,1)=sum(fs%divu_x(:,i,j,k)*this%T(i-1:i,j,k,1))+&
               &                  sum(fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     +&
               &                  sum(fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
               this%divT(i,j,k,2)=sum(fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     +&
               &                  sum(fs%divv_y(:,i,j,k)*this%T(i,j-1:j,k,4))+&
               &                  sum(fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
               this%divT(i,j,k,3)=sum(fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     +&
               &                  sum(fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     +&                  
               &                  sum(fs%divw_z(:,i,j,k)*this%T(i,j,k-1:k,6))        
            end do
         end do
      end do

      ! Deallocate arrays
      deallocate(Txy,Tyz,Tzx)

   end subroutine get_divT
   
end module fene_class
