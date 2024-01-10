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
      real(WP) :: visc_p                                     !< viscosity of polymer

   contains
      procedure :: init                                    !< Initialization of tpviscoelastic class (different name is used because of extension...)
      generic :: get_CgradU=>get_CgradU_standard,get_CgradU_log_conformation
      procedure, private :: get_CgradU_standard            !< Calculate upper-convected time derivative for viscoelastic constitutive model
      procedure, private :: get_CgradU_log_conformation    !< Calculate upper-convected time derivative for viscoelastic constitutive model with log conformation stabilization
      generic :: get_relax=>get_relax_standard,get_relax_log_conformation
      procedure, private :: get_relax_standard             !< Calculate relaxation term for viscoelastic constitutive model
      procedure, private :: get_relax_log_conformation     !< Calculate relaxation term for viscoelastic constitutive model with log conformation stabilization
      procedure :: get_eigensystem                         !< Procedure to calculate eigenvalues and eigenvectors for Conformation tensor
   end type tpviscoelastic
   
   
contains
   
   
   !> Viscoelastic model initialization
   subroutine init(this,cfg,phase,model,name)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: phase
      integer, intent(in) :: model
      character(len=*), optional :: name
      ! Create a six-scalar solver for conformation tensor in the liquid
      call this%tpscalar%initialize(cfg=cfg,nscalar=6,name=name)
      this%phase=phase; this%SCname=['Cxx','Cxy','Cxz','Cyy','Cyz','Czz']
      ! Assign closure model for viscoelastic fluid
      this%model=model
   end subroutine init
   
   
   !> Add upper-convected time derivative
   subroutine get_CgradU_standard(this,gradU,resSC,SR)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),      optional :: SR
      integer :: i,j,k
      resSC=0.0_WP
      select case (this%model)
      case (oldroydb,fenep,fenecr) !< For all other models
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
      case (lptt,eptt) !< Upper convected derivative in PTT models has additional term for affine motion
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! xx tensor component
               resSC(i,j,k,1)=2.0_WP*(this%SC(i,j,k,1)*gradU(1,1,i,j,k)+this%SC(i,j,k,2)*gradU(2,1,i,j,k)+this%SC(i,j,k,3)*gradU(3,1,i,j,k))-&
               &              this%affinecoeff*2.0_WP*(SR(1,i,j,k)*this%SC(i,j,k,1)+SR(4,i,j,k)*this%SC(i,j,k,2)+SR(6,i,j,k)*this%SC(i,j,k,3))
               ! xy tensor component
               resSC(i,j,k,2)=(this%SC(i,j,k,2)*gradU(1,1,i,j,k)+this%SC(i,j,k,4)*gradU(2,1,i,j,k)+this%SC(i,j,k,5)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,2,i,j,k)+this%SC(i,j,k,2)*gradU(2,2,i,j,k)+this%SC(i,j,k,3)*gradU(3,2,i,j,k))-&
               &              this%affinecoeff*((SR(4,i,j,k)*this%SC(i,j,k,1)+SR(2,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,3))+&
               &                                (SR(1,i,j,k)*this%SC(i,j,k,2)+SR(4,i,j,k)*this%SC(i,j,k,4)+SR(6,i,j,k)*this%SC(i,j,k,5)))
               ! xz tensor component
               resSC(i,j,k,3)=(this%SC(i,j,k,3)*gradU(1,1,i,j,k)+this%SC(i,j,k,5)*gradU(2,1,i,j,k)+this%SC(i,j,k,6)*gradU(3,1,i,j,k))+&
               &              (this%SC(i,j,k,1)*gradU(1,3,i,j,k)+this%SC(i,j,k,2)*gradU(2,3,i,j,k)+this%SC(i,j,k,3)*gradU(3,3,i,j,k))-&
               &              this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,1)+SR(5,i,j,k)*this%SC(i,j,k,2)+SR(3,i,j,k)*this%SC(i,j,k,3))+&
               &                                (SR(1,i,j,k)*this%SC(i,j,k,3)+SR(4,i,j,k)*this%SC(i,j,k,5)+SR(6,i,j,k)*this%SC(i,j,k,6)))
               ! yy tensor component
               resSC(i,j,k,4)=2.0_WP*(this%SC(i,j,k,2)*gradU(1,2,i,j,k)+this%SC(i,j,k,4)*gradU(2,2,i,j,k)+this%SC(i,j,k,5)*gradU(3,2,i,j,k))-&
               &              this%affinecoeff*2.0_WP*(SR(4,i,j,k)*this%SC(i,j,k,2)+SR(2,i,j,k)*this%SC(i,j,k,4)+SR(5,i,j,k)*this%SC(i,j,k,5))
               ! yz tensor component
               resSC(i,j,k,5)=(this%SC(i,j,k,2)*gradU(1,3,i,j,k)+this%SC(i,j,k,4)*gradU(2,3,i,j,k)+this%SC(i,j,k,5)*gradU(3,3,i,j,k))+&
               &              (this%SC(i,j,k,3)*gradU(1,2,i,j,k)+this%SC(i,j,k,5)*gradU(2,2,i,j,k)+this%SC(i,j,k,6)*gradU(3,2,i,j,k))-&
               &              this%affinecoeff*((SR(6,i,j,k)*this%SC(i,j,k,2)+SR(5,i,j,k)*this%SC(i,j,k,4)+SR(3,i,j,k)*this%SC(i,j,k,5))+&
               &                                (SR(4,i,j,k)*this%SC(i,j,k,3)+SR(2,i,j,k)*this%SC(i,j,k,5)+SR(5,i,j,k)*this%SC(i,j,k,6)))
               ! zz tensor component
               resSC(i,j,k,6)=2.0_WP*(this%SC(i,j,k,3)*gradU(1,3,i,j,k)+this%SC(i,j,k,5)*gradU(2,3,i,j,k)+this%SC(i,j,k,6)*gradU(3,3,i,j,k))-&
               &              this%affinecoeff*2.0_WP*(SR(6,i,j,k)*this%SC(i,j,k,3)+SR(5,i,j,k)*this%SC(i,j,k,5)+SR(3,i,j,k)*this%SC(i,j,k,6))
            end do
         end do
      end do
      end select
   end subroutine get_CgradU_standard


   !> Add upper-convected time derivative using log-conformation stabilization 
   !> Assumes scalar being transported is ln(C)
   subroutine get_CgradU_log_conformation(this,Eigenvalues,Eigenvectors,gradu,resSC,SR)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:),    intent(in)    :: Eigenvalues   
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:,1:), intent(in)    :: Eigenvectors
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: gradU
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(in)    :: SR
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:),    intent(inout) :: resSC
      integer :: i,j,k
      ! Temp scalar values for matrix multiplication
      real(WP), dimension(3,3) :: tmpMat,M,B,Omega   !< Matrices for diagonalization 
      real(WP), dimension(3,3) :: D,L                !< Matrices for diagonalization in PTT models
      real(WP) :: omega_xy,omega_xz,omega_yz         !< Components for anti-symetric matric
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Check if C is proportional to I based upon C's eigenvalues (i.e., Lambda_ii=Lambda_jj)
               if (abs(Eigenvalues(i,j,k,1)-Eigenvalues(i,j,k,2)).le.1.0e-15_WP.or.abs(Eigenvalues(i,j,k,2)-Eigenvalues(i,j,k,3)).le.1.0e-15_WP.or.abs(Eigenvalues(i,j,k,3)-Eigenvalues(i,j,k,1)).le.1.0e-15_WP) then
                  !>Set B equal to the strain rate tensor
                  B(1,1)=SR(1,i,j,k); B(1,2)=SR(4,i,j,k); B(1,3)=SR(6,i,j,k)
                  B(2,1)=SR(4,i,j,k); B(2,2)=SR(2,i,j,k); B(2,3)=SR(5,i,j,k)
                  B(3,1)=SR(6,i,j,k); B(3,2)=SR(5,i,j,k); B(3,3)=SR(3,i,j,k)
                  !>Set Omega = 0
                  Omega=0.0_WP
               else
                  ! Form M tensor (M=R^T*gradU^T*R={{mxx,mxy,mxz},{myx,myy,myz},{mzx,mzy,mzz}})
                  if (this%model.eq.lptt.or.this%model.eq.eptt) then
                     !>Strain rate tensor useed for calculating L
                     D(1,1)=SR(1,i,j,k); D(1,2)=SR(4,i,j,k); D(1,3)=SR(6,i,j,k)
                     D(2,1)=SR(4,i,j,k); D(2,2)=SR(2,i,j,k); D(2,3)=SR(5,i,j,k)
                     D(3,1)=SR(6,i,j,k); D(3,2)=SR(5,i,j,k); D(3,3)=SR(3,i,j,k)
                     ! Form L tensor L=gradU^T-zeta*D
                     L=transpose(gradU(:,:,i,j,k))-this%affinecoeff*D
                     M=matmul(transpose(Eigenvectors(i,j,k,:,:)),matmul(transpose(L),Eigenvectors(i,j,k,:,:)))
                  else 
                     M=matmul(transpose(Eigenvectors(i,j,k,:,:)),matmul(transpose(gradU(:,:,i,j,k)),Eigenvectors(i,j,k,:,:)))
                  end if
                  ! Temp matrix for calculating B
                  tmpMat=reshape((/ M(1,1),0.0_WP,0.0_WP,0.0_WP,M(2,2),0.0_WP,0.0_WP,0.0_WP,M(3,3) /),shape(tmpMat))
                  ! Form symmetric extension component of confomration tensor (B=R*{{mxx,0,0},{0,myy,0},{0,0,mzz}}*R^T={{Bxx,Bxy,Bxz},{Bxy,Byy,Byz},{Bxz,Byz,Bzz}})
                  B=matmul(Eigenvectors(i,j,k,:,:),matmul(tmpMat,transpose(Eigenvectors(i,j,k,:,:))))
                  ! Antisymmetric components
                  omega_xy=(Eigenvalues(i,j,k,2)*M(1,2)+Eigenvalues(i,j,k,1)*M(2,1))/(Eigenvalues(i,j,k,2)-Eigenvalues(i,j,k,1)) 
                  omega_xz=(Eigenvalues(i,j,k,3)*M(1,3)+Eigenvalues(i,j,k,1)*M(3,1))/(Eigenvalues(i,j,k,3)-Eigenvalues(i,j,k,1)) 
                  omega_yz=(Eigenvalues(i,j,k,3)*M(2,3)+Eigenvalues(i,j,k,2)*M(3,2))/(Eigenvalues(i,j,k,3)-Eigenvalues(i,j,k,2))
                  ! Temp matrix for calculating Omega
                  tmpMat=0.0_WP
                  tmpMat=reshape((/ 0.0_WP,-omega_xy,-omega_xz,omega_xy,0.0_WP,-omega_yz,omega_xz,omega_yz,0.0_WP /),shape(tmpMat))
                  ! Form rotation component of conformation tensor (Omega=R*{{0,omega_xy,omega_xz},{-omega_xy,0,omega_yz},{-omega_xz,-omega_yz,0}}*R^T={{Omegaxx,Omegaxy,Omegaxz},{Omegayx,Omegayy,Omegayz},{Omegazx,Omegazy,Omegazz}})
                  Omega=matmul(Eigenvectors(i,j,k,:,:),matmul(tmpMat,transpose(Eigenvectors(i,j,k,:,:))))
               end if
               ! Add extension and rotation components to resSC (Omega*log(C)-log(C)*Omega+2B)
               !>xx tensor component
               resSC(i,j,k,1)=Omega(1,2)*this%SC(i,j,k,2)-Omega(2,1)*this%SC(i,j,k,2)+Omega(1,3)*this%SC(i,j,k,3)-Omega(3,1)*this%SC(i,j,k,3)+2.00_WP*B(1,1)
               !>xy tensor component
               resSC(i,j,k,2)=Omega(1,1)*this%SC(i,j,k,2)-Omega(1,2)*this%SC(i,j,k,1)-Omega(2,2)*this%SC(i,j,k,2)-Omega(3,2)*this%SC(i,j,k,3)+Omega(1,2)*this%SC(i,j,k,4)+Omega(1,3)*this%SC(i,j,k,5)+2.00_WP*B(2,1)
               !>xz tensor component
               resSC(i,j,k,3)=Omega(1,1)*this%SC(i,j,k,3)-Omega(1,3)*this%SC(i,j,k,1)-Omega(2,3)*this%SC(i,j,k,2)-Omega(3,3)*this%SC(i,j,k,3)+Omega(1,2)*this%SC(i,j,k,5)+Omega(1,3)*this%SC(i,j,k,6)+2.00_WP*B(3,1)
               !>yy tensor component
               resSC(i,j,k,4)=Omega(2,1)*this%SC(i,j,k,2)-Omega(1,2)*this%SC(i,j,k,2)+Omega(2,3)*this%SC(i,j,k,5)-Omega(3,2)*this%SC(i,j,k,5)+2.00_WP*B(2,2)
               !>yz tensor component
               resSC(i,j,k,5)=Omega(2,1)*this%SC(i,j,k,3)-Omega(1,3)*this%SC(i,j,k,2)-Omega(2,3)*this%SC(i,j,k,4)+Omega(2,2)*this%SC(i,j,k,5)-Omega(3,3)*this%SC(i,j,k,5)+Omega(2,3)*this%SC(i,j,k,6)+2.00_WP*B(2,3)
               !>zz tensor component
               resSC(i,j,k,6)=Omega(3,1)*this%SC(i,j,k,3)-Omega(1,3)*this%SC(i,j,k,3)-Omega(2,3)*this%SC(i,j,k,5)+Omega(3,2)*this%SC(i,j,k,5)+2.00_WP*B(3,3)
            end do 
         end do
      end do
   end subroutine get_CgradU_log_conformation


   !> Add viscoelastic relaxation source
   subroutine get_relax_standard(this,resSC,dt)
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
   end subroutine get_relax_standard


   !> Add viscoelastic relaxation source using log-conformation stabilization 
   !> Assumes scalar being transported is ln(C)
   subroutine get_relax_log_conformation(this,Eigenvalues,Eigenvectors,resSC)
      use messager, only: die
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:),    intent(in) :: Eigenvalues   
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:,1:), intent(in) :: Eigenvectors
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      real(WP) :: coeff,trace
      resSC=0.0_WP
      select case (this%model)
      case (fenep) ! Add relaxation source for FENE-P 
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=Eigenvectors(i,j,k,1,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,2,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,3,1)**2*Eigenvalues(i,j,k,1)+&
                  &     Eigenvectors(i,j,k,1,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,2,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,3,2)**2*Eigenvalues(i,j,k,2)+&
                  &     Eigenvectors(i,j,k,1,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,2,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,3,3)**2*Eigenvalues(i,j,k,3)
                  !>Relaxation function coefficent
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-trace)
                  !>Add source term to residual
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)**2                     *((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,2,1)**2                     *((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,3,1)**2                     *((coeff/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((coeff/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((coeff/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (fenecr) ! Add relaxation source for FENE-CR 
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=Eigenvectors(i,j,k,1,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,2,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,3,1)**2*Eigenvalues(i,j,k,1)+&
                  &     Eigenvectors(i,j,k,1,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,2,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,3,2)**2*Eigenvalues(i,j,k,2)+&
                  &     Eigenvectors(i,j,k,1,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,2,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,3,3)**2*Eigenvalues(i,j,k,3)
                  !>Relaxation function coefficent
                  coeff=1.00_WP/(1.0_WP-trace/this%Lmax**2)
                  !>Add source term to residual
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,3,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (oldroydb) ! Add relaxation source for Oldroyd-B 
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,2,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*(Eigenvectors(i,j,k,3,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (lptt) ! Add relaxation source for lPTT
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=Eigenvectors(i,j,k,1,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,2,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,3,1)**2*Eigenvalues(i,j,k,1)+&
                  &     Eigenvectors(i,j,k,1,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,2,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,3,2)**2*Eigenvalues(i,j,k,2)+&
                  &     Eigenvectors(i,j,k,1,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,2,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,3,3)**2*Eigenvalues(i,j,k,3)
                  !>Relaxation function coefficent
                  coeff=1.00_WP+(this%elongvisc/(1.0_WP-this%affinecoeff))*(trace-3.0_WP)
                  !>Add source term to residual
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,3,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (eptt) ! Add relaxation source for ePTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=Eigenvectors(i,j,k,1,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,2,1)**2*Eigenvalues(i,j,k,1)+Eigenvectors(i,j,k,3,1)**2*Eigenvalues(i,j,k,1)+&
                  &     Eigenvectors(i,j,k,1,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,2,2)**2*Eigenvalues(i,j,k,2)+Eigenvectors(i,j,k,3,2)**2*Eigenvalues(i,j,k,2)+&
                  &     Eigenvectors(i,j,k,1,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,2,3)**2*Eigenvalues(i,j,k,3)+Eigenvectors(i,j,k,3,3)**2*Eigenvalues(i,j,k,3)
                  !>Relaxation function coefficent
                  coeff=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*(trace-3.0_WP))
                  !>Add source term to residual
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,2,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,2,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,2,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,1,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,1,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,1,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,2,1)*Eigenvectors(i,j,k,3,1)*((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,2,2)*Eigenvectors(i,j,k,3,2)*((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,2,3)*Eigenvectors(i,j,k,3,3)*((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*coeff*(Eigenvectors(i,j,k,3,1)**2                     *((1.00_WP/Eigenvalues(i,j,k,1))-1.00_WP)+Eigenvectors(i,j,k,3,2)**2                     *((1.00_WP/Eigenvalues(i,j,k,2))-1.00_WP)+Eigenvectors(i,j,k,3,3)**2                     *((1.00_WP/Eigenvalues(i,j,k,3))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case default
         call die('[tpviscoelastic get_relax] Unknown viscoelastic model selected')
      end select
   end subroutine get_relax_log_conformation


   !> Calculate the eigenvalues and eigenevtors of the conformation tensor for the whole domain
   !> Assumes scalar being transported is ln(C)
   subroutine get_eigensystem(this,Eigenvalues,Eigenvectors)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:),    intent(out) :: Eigenvalues   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:,1:), intent(out) :: Eigenvectors  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3,1:3)
      integer :: i,j,k
      ! Used for calculation of conformation tensor eigenvalues and eigenvectors
      real(WP), dimension(3,3) :: A
      real(WP), dimension(3) :: d
      integer , parameter :: order = 3             !< Conformation tensor is 3x3
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      integer  :: lwork,info,n
      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      Eigenvalues=0.0_WP
      Eigenvectors=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               !>Assemble ln(C) tensor
               A(1,1)=this%SC(i,j,k,1); A(1,2)=this%SC(i,j,k,2); A(1,3)=this%SC(i,j,k,3)
               A(2,1)=this%SC(i,j,k,2); A(2,2)=this%SC(i,j,k,4); A(2,3)=this%SC(i,j,k,5)
               A(3,1)=this%SC(i,j,k,3); A(3,2)=this%SC(i,j,k,5); A(3,3)=this%SC(i,j,k,6)
               !>On exit, A contains eigenvectors, and d contains ln(eigenvalues) in ascending order
               call dsyev('V','U',order,A,order,d,work,lwork,info)
               !>Eigenvalues (Eigenvalues(1)=Lambda_xx,Eigenvalues(2)=Lambda_yy,Eigenvalues(3)=Lambda_zz)
               Eigenvalues(i,j,k,:)=exp(d(:))
               !>Eigenvectors
               Eigenvectors(i,j,k,:,:)=A(:,:)
            end do 
         end do
      end do
   end subroutine get_eigensystem


end module tpviscoelastic_class