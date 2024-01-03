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
      procedure :: stabilization                           !< log-conformation stabilization for upper convected derivative
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

   !> Log-conformation stabilization method for upper convective derivative
   !> Assumes scalar being transported is log(C)
   subroutine stabilization(this,gradu,resSC)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      ! Temp scalar values for matrix multiplication
      real(WP) :: Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz                                        !< Matrix containing eigenvectors of C
      real(WP) :: mxx,mxy,mxz,myx,myy,myz,mzx,mzy,mzz                                        !< Components of M tensor
      real(WP) :: Lambdax,Lambday,Lambdaz                                                    !< Eigenvalues of C
      real(WP) :: omega_xy,omega_xz,omega_yz                                                 !< Components for anti-symetric matric
      real(WP) :: Bxx,Bxy,Bxz,Byy,Byz,Bzz                                                    !< Extensional component of conformation tensor
      real(WP) :: Omegaxx,Omegaxy,Omegaxz,Omegayx,Omegayy,Omegayz,Omegazx,Omegazy,Omegazz    !< Rotational component of conformation tensor
      ! Used for calculation of conformation tensor eigenvalues and eigenvectors
      real(WP), dimension(3,3) :: A
      real(WP), dimension(3) :: d
      integer , parameter :: order = 3             !< Conformation tensor is 3x3
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info,n

      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle

               ! Eigenvalues/eigenvectors of conformation tensor
               !>Assemble conformation tensor (from C=exp(log(C)))
               A(1,1)=exp(this%SC(i,j,k,1)); A(1,2)=exp(this%SC(i,j,k,2)); A(1,3)=exp(this%SC(i,j,k,3))
               A(2,1)=exp(this%SC(i,j,k,2)); A(2,2)=exp(this%SC(i,j,k,4)); A(2,3)=exp(this%SC(i,j,k,5))
               A(3,1)=exp(this%SC(i,j,k,3)); A(3,2)=exp(this%SC(i,j,k,5)); A(3,3)=exp(this%SC(i,j,k,6))
               !>On exit, A contains eigenvectors, and d contains eigenvalues in ascending order
               call dsyev('V','U',order,A,order,d,work,lwork,info)
               !>Form eigenvector tensor (R={{Rxx,Rxy,Rxz},{Ryx,Ryy,Ryz},{Rzx,Rzy,Rzz}})
               Rxx=A(1,1); Rxy=A(1,2); Rxz=A(1,3)
               Ryx=A(2,1); Ryy=A(2,2); Ryz=A(2,3)
               Rzx=A(3,1); Rzy=A(3,2); Rzz=A(3,3)
               !>Eigenvalues for conformation tensor (Lambdax,Lambday,Lambdaz)
               Lambdax=d(1); Lambday=d(2); Lambdaz=d(3)

               ! Form M tensor (M=R^T*gradU^T*R={{mxx,mxy,mxz},{myx,myy,myz},{mzx,mzy,mzz}})
               mxx=Rxx*(gradU(1,1,i,j,k)*Rxx+gradU(1,2,i,j,k)*Ryx+gradU(1,3,i,j,k)*Rzx)+Ryx*(gradU(2,1,i,j,k)*Rxx+gradU(2,2,i,j,k)*Ryx+gradU(2,3,i,j,k)*Rzx)+Rzx*(gradU(3,1,i,j,k)*Rxx+gradU(3,2,i,j,k)*Ryx+gradU(3,3,i,j,k)*Rzx)
               mxy=Rxy*(gradU(1,1,i,j,k)*Rxx+gradU(1,2,i,j,k)*Ryx+gradU(1,3,i,j,k)*Rzx)+Ryy*(gradU(2,1,i,j,k)*Rxx+gradU(2,2,i,j,k)*Ryx+gradU(2,3,i,j,k)*Rzx)+Rzy*(gradU(3,1,i,j,k)*Rxx+gradU(3,2,i,j,k)*Ryx+gradU(3,3,i,j,k)*Rzx)
               mxz=Rxz*(gradU(1,1,i,j,k)*Rxx+gradU(1,2,i,j,k)*Ryx+gradU(1,3,i,j,k)*Rzx)+Ryz*(gradU(2,1,i,j,k)*Rxx+gradU(2,2,i,j,k)*Ryx+gradU(2,3,i,j,k)*Rzx)+Rzz*(gradU(3,1,i,j,k)*Rxx+gradU(3,2,i,j,k)*Ryx+gradU(3,3,i,j,k)*Rzx)
               myx=Rxx*(gradU(1,1,i,j,k)*Rxy+gradU(1,2,i,j,k)*Ryy+gradU(1,3,i,j,k)*Rzy)+Ryx*(gradU(2,1,i,j,k)*Rxy+gradU(2,2,i,j,k)*Ryy+gradU(2,3,i,j,k)*Rzy)+Rzx*(gradU(3,1,i,j,k)*Rxy+gradU(3,2,i,j,k)*Ryy+gradU(3,3,i,j,k)*Rzy)
               myy=Rxy*(gradU(1,1,i,j,k)*Rxy+gradU(1,2,i,j,k)*Ryy+gradU(1,3,i,j,k)*Rzy)+Ryy*(gradU(2,1,i,j,k)*Rxy+gradU(2,2,i,j,k)*Ryy+gradU(2,3,i,j,k)*Rzy)+Rzy*(gradU(3,1,i,j,k)*Rxy+gradU(3,2,i,j,k)*Ryy+gradU(3,3,i,j,k)*Rzy)
               myz=Rxz*(gradU(1,1,i,j,k)*Rxy+gradU(1,2,i,j,k)*Ryy+gradU(1,3,i,j,k)*Rzy)+Ryz*(gradU(2,1,i,j,k)*Rxy+gradU(2,2,i,j,k)*Ryy+gradU(2,3,i,j,k)*Rzy)+Rzz*(gradU(3,1,i,j,k)*Rxy+gradU(3,2,i,j,k)*Ryy+gradU(3,3,i,j,k)*Rzy)
               mzx=Rxx*(gradU(1,1,i,j,k)*Rxz+gradU(1,2,i,j,k)*Ryz+gradU(1,3,i,j,k)*Rzz)+Ryx*(gradU(2,1,i,j,k)*Rxz+gradU(2,2,i,j,k)*Ryz+gradU(2,3,i,j,k)*Rzz)+Rzx*(gradU(3,1,i,j,k)*Rxz+gradU(3,2,i,j,k)*Ryz+gradU(3,3,i,j,k)*Rzz)
               mzy=Rxy*(gradU(1,1,i,j,k)*Rxz+gradU(1,2,i,j,k)*Ryz+gradU(1,3,i,j,k)*Rzz)+Ryy*(gradU(2,1,i,j,k)*Rxz+gradU(2,2,i,j,k)*Ryz+gradU(2,3,i,j,k)*Rzz)+Rzy*(gradU(3,1,i,j,k)*Rxz+gradU(3,2,i,j,k)*Ryz+gradU(3,3,i,j,k)*Rzz)
               mzz=Rxz*(gradU(1,1,i,j,k)*Rxz+gradU(1,2,i,j,k)*Ryz+gradU(1,3,i,j,k)*Rzz)+Ryz*(gradU(2,1,i,j,k)*Rxz+gradU(2,2,i,j,k)*Ryz+gradU(2,3,i,j,k)*Rzz)+Rzz*(gradU(3,1,i,j,k)*Rxz+gradU(3,2,i,j,k)*Ryz+gradU(3,3,i,j,k)*Rzz)
            
               ! Form symmetric extension component of confomration tensor (B=R*{{mxx,0,0},{0,myy,0},{0,0,mzz}}*R^T={{Bxx,Bxy,Bxz},{Bxy,Byy,Byz},{Bxz,Byz,Bzz}})
               !>xx tensor component
               Bxx=mxx*Rxx**2+myy*Rxy**2+mzz*Rxz**2
               !>xy tensor component
               Bxy=mxx*Rxx*Ryx+myy*Rxy*Ryy+mzz*Rxz*Ryz
               !>xz tensor component
               Bxz=mxx*Rxx*Rzx+myy*Rxy*Rzy+mzz*Rxz*Rzz
               !>yy tensor component
               Byy=mxx*Rxy**2+myy*Ryy**2+mzz*Ryz**2
               !>yz tensor component
               Byz=mxx*Ryx*Rzx+myy*Ryy*Rzy+mzz*Ryz*Rzz
               !>zz tensor component
               Bzz=mxx*Rzx**2+myy*Rzy**2+mzz*Rzz**2

               ! Form rotation component of conformation tensor (Omega=R*{{0,omega_xy,omega_xz},{-omega_xy,0,omega_yz},{-omega_xz,-omega_yz,0}}*R^T={{Omegaxx,Omegaxy,Omegaxz},{Omegayx,Omegayy,Omegayz},{Omegazx,Omegazy,Omegazz}})
               !>Antisymmetric components
               omega_xy=(Lambday*mxy+Lambdax*myx)/(Lambday-Lambdax); omega_xz=(Lambdaz*mxz+Lambdax*mzx)/(Lambdaz-Lambdax); omega_yz=(Lambdaz*myz+Lambday*mzy)/(Lambdaz-Lambday)
               !>Tensor components
               Omegaxx=Rxz*(Rxx*omega_xz+Rxy*omega_yz)-Rxx*(Rxy*omega_xy+Rxz*omega_xz)+Rxy*(Rxx*omega_xy-Rxz*omega_yz)
               Omegaxy=Ryz*(Rxx*omega_xz+Rxy*omega_yz)-Ryx*(Rxy*omega_xy+Rxz*omega_xz)+Ryy*(Rxx*omega_xy-Rxz*omega_yz)
               Omegaxz=Rzz*(Rxx*omega_xz+Rxy*omega_yz)-Rzx*(Rxy*omega_xy+Rxz*omega_xz)+Rzy*(Rxx*omega_xy-Rxz*omega_yz)
               Omegayx=Rxz*(Ryx*omega_xz+Ryy*omega_yz)-Rxx*(Ryy*omega_xy+Ryz*omega_xz)+Rxy*(Ryx*omega_xy-Ryz*omega_yz)
               Omegayy=Ryz*(Ryx*omega_xz+Ryy*omega_yz)-Ryx*(Ryy*omega_xy+Ryz*omega_xz)+Ryy*(Ryx*omega_xy-Ryz*omega_yz)
               Omegayz=Rzz*(Ryx*omega_xz+Ryy*omega_yz)-Rzx*(Ryy*omega_xy+Ryz*omega_xz)+Rzy*(Ryx*omega_xy-Ryz*omega_yz)
               Omegazx=Rxz*(Rzx*omega_xz+Rzy*omega_yz)-Rxx*(Rzy*omega_xy+Rzz*omega_xz)+Rxy*(Rzx*omega_xy-Rzz*omega_yz)
               Omegazy=Ryz*(Rzx*omega_xz+Rzy*omega_yz)-Ryx*(Rzy*omega_xy+Rzz*omega_xz)+Ryy*(Rzx*omega_xy-Rzz*omega_yz)
               Omegazz=Rzz*(Rzx*omega_xz+Rzy*omega_yz)-Rzx*(Rzy*omega_xy+Rzz*omega_xz)+Rzy*(Rzx*omega_xy-Rzz*omega_yz)

               ! Add extension and rotation components to resSC (Omega*log(C)-log(C)*Omega+2B)
               !>xx tensor component
               resSC(i,j,k,1)=Omegaxy*this%SC(i,j,k,2)-Omegayx*this%SC(i,j,k,2)+Omegaxz*this%SC(i,j,k,3)-Omegazx*this%SC(i,j,k,3)+2.00_WP*Bxx
               !>xy tensor component
               resSC(i,j,k,2)=Omegaxx*this%SC(i,j,k,2)-Omegaxy*this%SC(i,j,k,1)-Omegayy*this%SC(i,j,k,2)-Omegazy*this%SC(i,j,k,3)+Omegaxy*this%SC(i,j,k,4)+Omegaxz*this%SC(i,j,k,5)+2.00_WP*Bxy
               !>xz tensor component
               resSC(i,j,k,3)=Omegaxx*this%SC(i,j,k,3)-Omegaxz*this%SC(i,j,k,1)-Omegayz*this%SC(i,j,k,2)-Omegazz*this%SC(i,j,k,3)+Omegaxy*this%SC(i,j,k,5)+Omegaxz*this%SC(i,j,k,6)+2.00_WP*Bxz
               !>yy tensor component
               resSC(i,j,k,4)=Omegaxy*this%SC(i,j,k,2)-Omegaxy*this%SC(i,j,k,2)+Omegayz*this%SC(i,j,k,5)-Omegazy*this%SC(i,j,k,5)+2.00_WP*Byy
               !>yz tensor component
               resSC(i,j,k,5)=Omegayx*this%SC(i,j,k,3)-Omegaxz*this%SC(i,j,k,2)-Omegayz*this%SC(i,j,k,4)+Omegayy*this%SC(i,j,k,5)-Omegazz*this%SC(i,j,k,5)+Omegayz*this%SC(i,j,k,6)+2.00_WP*Byz
               !>zz tensor component
               resSC(i,j,k,6)=Omegazx*this%SC(i,j,k,3)-Omegaxz*this%SC(i,j,k,3)-Omegayz*this%SC(i,j,k,5)+Omegazy*this%SC(i,j,k,5)+2.00_WP*Bzz

            end do 
         end do
      end do

   end subroutine stabilization
   
end module tpviscoelastic_class