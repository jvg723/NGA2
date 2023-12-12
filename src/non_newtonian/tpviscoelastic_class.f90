!> Two-phase viscoelastic model class
!> Extends tpscalar class for calculation of source terms
module tpviscoelastic_class
   use tpscalar_class, only: tpscalar
   use config_class,   only: config
   use precision,      only: WP
   implicit none
   private
   

   ! Expose type/constructor/methods
   public :: tpviscoelastic
   

   ! List of available viscoelastic models
   integer, parameter, public :: fenep =0            !< FENE-P model
   integer, parameter, public :: fenecr=1            !< FENE-CR model
   integer, parameter, public :: oldroydb=2          !< Oldroyd-B model
   integer, parameter, public :: lptt=3              !< Linear Phan-Thien-Tanner model
   integer, parameter, public :: eptt=4              !< Exponential Phan-Thien-Tanner model
   

   !> Constant density viscoelastic solver object definition
   type, extends(tpscalar) :: tpviscoelastic
      
      ! Model parameters
      integer  :: model                                      !< Closure model
      real(WP) :: trelax                                     !< Polymer relaxation timescale
      real(WP) :: Lmax                                       !< Polymer maximum extensibility in FENE model
      real(WP) :: affinecoeff                                !< Parameter for affine motion in PTT model
      real(WP) :: elongvisc                                  !< Extensional parameter for elognational viscosity in PTT model

      ! Arrays for log-conformation stabilization
      logical :: use_stabilization=.false.                   !< Flag to use log-conformation approach (default=.false.)
      real(WP), dimension(:,:,:,:), allocatable :: extension !< Extension matrix decomposed from gradU
      real(WP), dimension(:,:,:,:), allocatable :: rotation  !< Rotation matrix decomposed from gradU

   contains
      procedure :: init                                    !< Initialization of tpviscoelastic class (different name is used because of extension...)
      procedure :: get_CgradU                              !< Calculate streching and distrortion term
      procedure :: get_relax                               !< Calculate relaxation term
      procedure :: get_affine                              !< Source term in PTT equation for non-affine motion
      procedure :: decompose_gradu                         !< Deompose velocity gradient for log-conformation stabilization technique
   end type tpviscoelastic
   
   
contains
   
   
   !> Viscoelastic model initialization
   subroutine init(this,cfg,phase,model,name,use_stabilization)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: phase
      integer, intent(in) :: model
      character(len=*), optional :: name
      logical, optional :: use_stabilization
      ! Create a six-scalar solver for conformation tensor in the liquid
      call this%tpscalar%initialize(cfg=cfg,nscalar=6,name=name)
      this%phase=phase; this%SCname=['Cxx','Cxy','Cxz','Cyy','Cyz','Czz']
      ! Assign closure model for viscoelastic fluid
      this%model=model
      ! Set optional detailed flux info
      if (present(use_stabilization)) then
         allocate(this%extension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6))
      end if
   end subroutine init
   
   
   !> Get CgradU source terms to add to multiscalar residual
   subroutine get_CgradU(this,gradU,resSC)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=2.0_WP*(this%SC(i,j,k,1)*gradU(1,1,i,j,k)+this%SC(i,j,k,2)*gradU(2,1,i,j,k)+this%SC(i,j,k,3)*gradU(3,1,i,j,k))
               ! xy tensor component
               resSC(i,j,k,2)=(this%SC(i,j,k,2)*gradU(1,1,i,j,k)+this%SC(i,j,k,4)*gradU(2,1,i,j,k)+this%SC(i,j,k,5)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,2,i,j,k)+this%SC(i,j,k,2)*gradU(2,2,i,j,k)+this%SC(i,j,k,3)*gradU(3,2,i,j,k))
               ! xz tensor component
               resSC(i,j,k,3)=(this%SC(i,j,k,3)*gradU(1,1,i,j,k)+this%SC(i,j,k,5)*gradU(2,1,i,j,k)+this%SC(i,j,k,6)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,3,i,j,k)+this%SC(i,j,k,2)*gradU(2,3,i,j,k)+this%SC(i,j,k,3)*gradU(3,3,i,j,k))
               ! yy tensor component
               resSC(i,j,k,4)=2.0_WP*(this%SC(i,j,k,2)*gradU(1,2,i,j,k)+this%SC(i,j,k,4)*gradU(2,2,i,j,k)+this%SC(i,j,k,5)*gradU(3,2,i,j,k))
               ! yz tensor component
               resSC(i,j,k,5)=(this%SC(i,j,k,2)*gradU(1,3,i,j,k)+this%SC(i,j,k,4)*gradU(2,3,i,j,k)+this%SC(i,j,k,5)*gradU(3,3,i,j,k))+&
               &              (this%SC(i,j,k,3)*gradU(1,2,i,j,k)+this%SC(i,j,k,5)*gradU(2,2,i,j,k)+this%SC(i,j,k,6)*gradU(3,2,i,j,k))
               ! zz tensor component
               resSC(i,j,k,6)=2.0_WP*(this%SC(i,j,k,3)*gradU(1,3,i,j,k)+this%SC(i,j,k,5)*gradU(2,3,i,j,k)+this%SC(i,j,k,6)*gradU(3,3,i,j,k))
            end do
         end do
      end do
   end subroutine get_CgradU
   

   !> Get S*C terms for PTT equation
   subroutine get_affine(this,SR,resSC)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=-this%affinecoeff*2.0_WP*(SR(1,i,j,k)*this%SC(i,j,k,1)+SR(4,i,j,k)*this%SC(i,j,k,2)+SR(6,i,j,k)*this%SC(i,j,k,3))
               ! xy tensor component
               resSC(i,j,k,2)=-this%affinecoeff*((SR(4,i,j,k)*this%SC(i,j,k,1)+SR(2,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,3))+&
               &                                 (SR(1,i,j,k)*this%SC(i,j,k,2)+SR(4,i,j,k)*this%SC(i,j,k,4)+SR(6,i,j,k)*this%SC(i,j,k,5)))
               ! xz tensor component
               resSC(i,j,k,3)=-this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,1)+SR(5,i,j,k)*this%SC(i,j,k,2)+SR(3,i,j,k)*this%SC(i,j,k,3))+&
               &                                 (SR(1,i,j,k)*this%SC(i,j,k,3)+SR(4,i,j,k)*this%SC(i,j,k,5)+SR(6,i,j,k)*this%SC(i,j,k,6)))
               ! yy tensor component
               resSC(i,j,k,4)=-this%affinecoeff*2.0_WP*(SR(4,i,j,k)*this%SC(i,j,k,2)+SR(2,i,j,k)*this%SC(i,j,k,4)+SR(5,i,j,k)*this%SC(i,j,k,5))
               ! yz tensor component
               resSC(i,j,k,5)=-this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,4)+SR(3,i,j,k)*this%SC(i,j,k,5))+&
               &                                 (SR(4,i,j,k)*this%SC(i,j,k,3)+SR(2,i,j,k)*this%SC(i,j,k,5)+SR(5,i,j,k)*this%SC(i,j,k,6)))
               ! zz tensor component
               resSC(i,j,k,6)=-this%affinecoeff*2.0_WP*(SR(6,i,j,k)*this%SC(i,j,k,3)+SR(5,i,j,k)*this%SC(i,j,k,5)+SR(3,i,j,k)*this%SC(i,j,k,6))
            end do
         end do
      end do
   end subroutine get_affine
   

   !> Add viscoelastic relaxation source
   subroutine get_relax(this,resSC,dt)
      use messager, only: die
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      real(WP), intent(in) :: dt
      integer :: i,j,k
      real(WP) :: coeff
      real(WP), parameter :: safety_margin=10.0_WP
      resSC=0.0_WP
      select case (this%model)
      case (fenep) ! Add relaxation source for FENE-P (1/lambda)(f(r)*C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6)))
                  resSC(i,j,k,1)=-(coeff*this%SC(i,j,k,1)-1.0_WP)/this%trelax   !< xx tensor component
                  resSC(i,j,k,2)=-(coeff*this%SC(i,j,k,2)-0.0_WP)/this%trelax   !< xy tensor component
                  resSC(i,j,k,3)=-(coeff*this%SC(i,j,k,3)-0.0_WP)/this%trelax   !< xz tensor component
                  resSC(i,j,k,4)=-(coeff*this%SC(i,j,k,4)-1.0_WP)/this%trelax   !< yy tensor component
                  resSC(i,j,k,5)=-(coeff*this%SC(i,j,k,5)-0.0_WP)/this%trelax   !< yz tensor component
                  resSC(i,j,k,6)=-(coeff*this%SC(i,j,k,6)-1.0_WP)/this%trelax   !< zz tensor component
               end do
            end do
         end do
      case (fenecr) ! Add relaxation source for FENE-CR (f(r)/lambda*(C-I))
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.0_WP-(this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))/this%Lmax**2
                  coeff=max(epsilon(1.0_WP),min(coeff,1.0_WP))*this%trelax      !< Build a safe adjusted relaxation time scale
                  coeff=max(coeff,safety_margin*dt)                             !< Further clip based on current time step size for stability
                  coeff=1.0_WP/coeff                                            !< Inverse coeff
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP)               !< xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP)               !< xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP)               !< xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP)               !< yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP)               !< yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP)               !< zz tensor component
               end do
            end do
         end do
      case (oldroydb) ! Add relaxation source for Oldroyd-B (1/t_relax)(C-I)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.0_WP/this%trelax                                      !< Inverse of relaxation time
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP)               !< xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP)               !< xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP)               !< xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP)               !< yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP)               !< yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP)               !< zz tensor component
               end do
            end do
         end do
      case (lptt) ! Add relaxation source term for lPTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=1.00_WP+(this%elongvisc/(1.0_WP-this%affinecoeff))*((this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))-3.0_WP)
                  coeff=coeff/this%trelax                                       !< Divide by relaxation time scale
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP)               !< xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP)               !< xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP)               !< xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP)               !< yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP)               !< yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP)               !< zz tensor component
               end do
            end do
         end do
      case (eptt) ! Add relaxation source term for ePTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  coeff=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*((this%SC(i,j,k,1)+this%SC(i,j,k,4)+this%SC(i,j,k,6))-3.0_WP))
                  coeff=coeff/this%trelax                                       !< Divide by relaxation time scale
                  resSC(i,j,k,1)=-coeff*(this%SC(i,j,k,1)-1.0_WP)               !< xx tensor component
                  resSC(i,j,k,2)=-coeff*(this%SC(i,j,k,2)-0.0_WP)               !< xy tensor component
                  resSC(i,j,k,3)=-coeff*(this%SC(i,j,k,3)-0.0_WP)               !< xz tensor component
                  resSC(i,j,k,4)=-coeff*(this%SC(i,j,k,4)-1.0_WP)               !< yy tensor component
                  resSC(i,j,k,5)=-coeff*(this%SC(i,j,k,5)-0.0_WP)               !< yz tensor component
                  resSC(i,j,k,6)=-coeff*(this%SC(i,j,k,6)-1.0_WP)               !< zz tensor component
               end do
            end do
         end do
      case default
         call die('[tpviscoelastic get_relax] Unknown viscoelastic model selected')
      end select
   end subroutine get_relax

   !> Decompose velocity gradient for log-conformation stabilization
   !> Assumes scalar being transported is log(C)
   !> Passes out B and Omega matrices
   subroutine decompose_gradu(this,gradu)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      integer :: i,j,k
      ! Eigenvalues/eigenvectors
      real(WP), dimension(3,3) :: A
      real(WP), dimension(3,3) :: M
      real(WP), dimension(3) :: d
      real(WP)               :: omega12,omega21,omega13,omega31,omega23,omega32
      integer , parameter :: order = 3             !< Conformation tensor is 3x3
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info,n

      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               
               ! Eigenvalues/eigenvectors of C tensor
               !>Assemble C matrix (from C=exp(log(C)))
               A(1,1)=exp(this%SC(i,j,k,1)); A(1,2)=exp(this%SC(i,j,k,2)); A(1,3)=exp(this%SC(i,j,k,3))
               A(2,1)=exp(this%SC(i,j,k,2)); A(2,2)=exp(this%SC(i,j,k,4)); A(2,3)=exp(this%SC(i,j,k,5))
               A(3,1)=exp(this%SC(i,j,k,3)); A(3,2)=exp(this%SC(i,j,k,5)); A(3,3)=exp(this%SC(i,j,k,6))
               !>On exit, A contains eigenvectors, and d contains eigenvalues in ascending order
               call dsyev('V','U',order,A,order,d,work,lwork,info)

               ! Change of base for the velocity gradient (A*gradU*At=M)
               !> M={{m11,m12,m13},{m21,m22,m23},{m31,m32,m33}}
               M(1,1)=A(1,1)*(A(1,1)*gradU(1,1,i,j,k)+A(1,2)*gradU(2,1,i,j,k)+A(1,3)*gradU(3,1,i,j,k))+A(1,2)*(A(1,1)*gradU(1,2,i,j,k)+A(1,2)*gradU(2,2,i,j,k)+A(1,3)*gradU(3,2,i,j,k))+A(1,3)*(A(1,1)*gradU(1,3,i,j,k)+A(1,2)*gradU(2,3,i,j,k)+A(1,3)*gradU(3,3,i,j,k))
               M(1,2)=A(2,1)*(A(1,1)*gradU(1,1,i,j,k)+A(1,2)*gradU(2,1,i,j,k)+A(1,3)*gradU(3,1,i,j,k))+A(2,2)*(A(1,1)*gradU(1,2,i,j,k)+A(1,2)*gradU(2,2,i,j,k)+A(1,3)*gradU(3,2,i,j,k))+A(2,3)*(A(1,1)*gradU(1,3,i,j,k)+A(1,2)*gradU(2,3,i,j,k)+A(1,3)*gradU(3,3,i,j,k))
               M(1,3)=A(3,1)*(A(1,1)*gradU(1,1,i,j,k)+A(1,2)*gradU(2,1,i,j,k)+A(1,3)*gradU(3,1,i,j,k))+A(3,2)*(A(1,1)*gradU(1,2,i,j,k)+A(1,2)*gradU(2,2,i,j,k)+A(1,3)*gradU(3,2,i,j,k))+A(3,3)*(A(1,1)*gradU(1,3,i,j,k)+A(1,2)*gradU(2,3,i,j,k)+A(1,3)*gradU(3,3,i,j,k))
               M(2,1)=A(1,1)*(A(2,1)*gradU(1,1,i,j,k)+A(2,2)*gradU(2,1,i,j,k)+A(2,3)*gradU(3,1,i,j,k))+A(1,2)*(A(2,1)*gradU(1,2,i,j,k)+A(2,2)*gradU(2,2,i,j,k)+A(2,3)*gradU(3,2,i,j,k))+A(1,3)*(A(2,1)*gradU(1,3,i,j,k)+A(2,2)*gradU(2,3,i,j,k)+A(2,3)*gradU(3,3,i,j,k))
               M(2,2)=A(2,1)*(A(2,1)*gradU(1,1,i,j,k)+A(2,2)*gradU(2,1,i,j,k)+A(2,3)*gradU(3,1,i,j,k))+A(2,2)*(A(2,1)*gradU(1,2,i,j,k)+A(2,2)*gradU(2,2,i,j,k)+A(2,3)*gradU(3,2,i,j,k))+A(2,3)*(A(2,1)*gradU(1,3,i,j,k)+A(2,2)*gradU(2,3,i,j,k)+A(2,3)*gradU(3,3,i,j,k))
               M(2,3)=A(3,1)*(A(2,1)*gradU(1,1,i,j,k)+A(2,2)*gradU(2,1,i,j,k)+A(2,3)*gradU(3,1,i,j,k))+A(3,2)*(A(2,1)*gradU(1,2,i,j,k)+A(2,2)*gradU(2,2,i,j,k)+A(2,3)*gradU(3,2,i,j,k))+A(3,3)*(A(2,1)*gradU(1,3,i,j,k)+A(2,2)*gradU(2,3,i,j,k)+A(2,3)*gradU(3,3,i,j,k))
               M(3,1)=A(1,1)*(A(3,1)*gradU(1,1,i,j,k)+A(3,2)*gradU(2,1,i,j,k)+A(3,3)*gradU(3,1,i,j,k))+A(1,2)*(A(3,1)*gradU(1,2,i,j,k)+A(3,2)*gradU(2,2,i,j,k)+A(3,3)*gradU(3,2,i,j,k))+A(1,3)*(A(3,1)*gradU(1,3,i,j,k)+A(3,2)*gradU(2,3,i,j,k)+A(3,3)*gradU(3,3,i,j,k))
               M(3,2)=A(2,1)*(A(3,1)*gradU(1,1,i,j,k)+A(3,2)*gradU(2,1,i,j,k)+A(3,3)*gradU(3,1,i,j,k))+A(2,2)*(A(3,1)*gradU(1,2,i,j,k)+A(3,2)*gradU(2,2,i,j,k)+A(3,3)*gradU(3,2,i,j,k))+A(2,3)*(A(3,1)*gradU(1,3,i,j,k)+A(3,2)*gradU(2,3,i,j,k)+A(3,3)*gradU(3,3,i,j,k))
               M(3,3)=A(3,1)*(A(3,1)*gradU(1,1,i,j,k)+A(3,2)*gradU(2,1,i,j,k)+A(3,3)*gradU(3,1,i,j,k))+A(3,2)*(A(3,1)*gradU(1,2,i,j,k)+A(3,2)*gradU(2,2,i,j,k)+A(3,3)*gradU(3,2,i,j,k))+A(3,3)*(A(3,1)*gradU(1,3,i,j,k)+A(3,2)*gradU(2,3,i,j,k)+A(3,3)*gradU(3,3,i,j,k))


               ! Compute extension matrix (symmetric)
               !>E11
               this%extension(i,j,k,1)=A(1,1)**2*M(1,1)+A(1,2)**2*M(2,2)+A(1,3)**2*M(3,3)
               !>E21
               this%extension(i,j,k,2)=A(1,1)*A(2,1)*M(1,1)+A(1,2)*A(2,2)*M(2,2)+A(1,3)*A(2,3)*M(3,3)
               !>E31
               this%extension(i,j,k,3)=A(1,1)*A(3,1)*M(1,1)+A(1,2)*A(3,2)*M(2,2)+A(1,3)*A(3,3)*M(3,3)
               !>E22
               this%extension(i,j,k,4)=A(2,1)**2*M(1,1)+A(2,2)**2*M(2,2)+A(2,3)**2*M(3,3)
               !>E32
               this%extension(i,j,k,5)=A(2,1)*A(3,1)*M(1,1)+A(2,2)*A(3,2)*M(2,2)+A(2,3)*A(3,3)*M(3,3)
               !>E33
               this%extension(i,j,k,6)=A(3,1)**2*M(1,1)+A(3,2)**2*M(2,2)+A(3,3)**2*M(3,3)

               ! Compute rotation matrix
               !> Matrix terms
               omega12=(d(2)*M(1,2)+d(1)*M(2,1))/(d(1)+d(2)); omega21=(d(1)*M(2,1)+d(2)*M(1,2))/(d(2)-d(1))
               omega23=(d(3)*M(2,3)+d(2)*M(3,2))/(d(2)-d(3)); omega32=(d(2)*M(3,2)+d(3)*M(2,3))/(d(3)-d(2))
               omega13=(d(3)*M(1,3)+d(1)*M(3,1))/(d(1)-d(3)); omega31=(d(1)*M(3,1)+d(3)*M(1,3))/(d(3)-d(1))

               !>Omega11
               this%rotation(i,j,k,1)=0.0_WP
               !>Omega21
               this%rotation(i,j,k,2)=0.0_WP
               !>Omega31
               this%rotation(i,j,k,3)=0.0_WP
               !>Omega22
               this%rotation(i,j,k,4)=0.0_WP
               !>Omega32
               this%rotation(i,j,k,5)=0.0_WP
               !>Omega33
               this%rotation(i,j,k,6)=0.0_WP
            end do 
         end do
      end do

   end subroutine decompose_gradu
   
end module tpviscoelastic_class