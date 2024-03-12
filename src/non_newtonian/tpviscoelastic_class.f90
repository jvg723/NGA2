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
      integer  :: model                                    !< Closure model
      real(WP) :: trelax                                   !< Polymer relaxation timescale
      real(WP) :: Lmax                                     !< Polymer maximum extensibility in FENE model
      real(WP) :: affinecoeff                              !< Parameter for affine motion in PTT model
      real(WP) :: elongvisc                                !< Extensional parameter for elongational viscosity in PTT model
      real(WP) :: visc_p                                   !< Viscosity of polymer
      ! Eigensystem storage
      real(WP), dimension(:,:,:,:),   allocatable :: eigenval  !< Field of eigenvalues of the conformation tensor
      real(WP), dimension(:,:,:,:,:), allocatable :: eigenvec  !< Field of eigenvectors of the conformation tensor
      ! Storage for reconstructed conformation tensor
      real(WP), dimension(:,:,:,:),   allocatable ::SCrec
      real(WP), dimension(:,:,:,:),   allocatable ::SCrecold
      ! Monitoring quantities
      real(WP), dimension(:), allocatable :: SCrecmax,SCrecmin,SCrecint   !< Maximum and minimum, integral reconstructed scalar feild 
      real(WP) :: eigval1_max,eigval2_max,eigval3_max
      contains
      procedure :: init                                    !< Initialization of tpviscoelastic class (different name is used because of extension...)
      procedure :: get_CgradU                              !< Calculate streching and distortion term
      procedure :: get_relax                               !< Calculate relaxation term
      procedure :: get_relax_analytical                    !< Calculate relaxation term based on a semi-analytical integration 
      procedure :: get_CgradU_log                          !< Calculate streching and distortion term for log-conformation tensor
      procedure :: get_relax_log                           !< Calculate relaxation term for log-conformation tensor
      procedure :: get_eigensystem                         !< Calculate eigenvalues and eigenvectors for conformation tensor
      procedure :: get_eigensystem_SCrec                   !< Calculate eigenvalues and eigenvectors for conformation tensor
      procedure :: reconstruct_conformation                !< Reconstruct conformation tensor for decomposed eigenvalues and eigenvectors
      procedure :: reconstruct_log_conformation            !< Reconstruct log conformation tensor for decomposed eigenvalues and eigenvectors
      procedure :: get_max_reconstructed                   !< Calculate maximum and integral field value for reconstructed C field
      procedure :: get_max_eig
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
   subroutine get_CgradU(this,gradU,resSC)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(inout) :: resSC
      integer :: i,j,k
      real(WP), dimension(6) :: SR
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
      ! Upper convected derivative in PTT models has additional term for affine motion
      select case (this%model)
      case (lptt,eptt)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Skip non-solved cells
                  if (this%mask(i,j,k).ne.0) cycle
                  ! Build strain rate tensor in same order as conformation tensor
                  SR(1)=gradU(1,1,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP
                  SR(2)=0.5_WP*(gradU(1,2,i,j,k)+gradU(2,1,i,j,k))
                  SR(3)=0.5_WP*(gradU(1,3,i,j,k)+gradU(3,1,i,j,k))
                  SR(4)=gradU(2,2,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP
                  SR(5)=0.5_WP*(gradU(2,3,i,j,k)+gradU(3,2,i,j,k))
                  SR(6)=gradU(3,3,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP
                  ! xx tensor component
                  resSC(i,j,k,1)=resSC(i,j,k,1)-this%affinecoeff*2.0_WP*(SR(1)*this%SC(i,j,k,1)+SR(2)*this%SC(i,j,k,2)+SR(3)*this%SC(i,j,k,3))
                  ! xy tensor component
                  resSC(i,j,k,2)=resSC(i,j,k,2)-this%affinecoeff*((SR(2)*this%SC(i,j,k,1)+SR(4)*this%SC(i,j,k,2)+SR(5)*this%SC(i,j,k,3))+&
                  &                                               (SR(1)*this%SC(i,j,k,2)+SR(2)*this%SC(i,j,k,4)+SR(3)*this%SC(i,j,k,5)))
                  ! xz tensor component
                  resSC(i,j,k,3)=resSC(i,j,k,3)-this%affinecoeff*((SR(3)*this%SC(i,j,k,1)+SR(5)*this%SC(i,j,k,2)+SR(6)*this%SC(i,j,k,3))+&
                  &                                               (SR(1)*this%SC(i,j,k,3)+SR(2)*this%SC(i,j,k,5)+SR(3)*this%SC(i,j,k,6)))
                  ! yy tensor component
                  resSC(i,j,k,4)=resSC(i,j,k,4)-this%affinecoeff*2.0_WP*(SR(2)*this%SC(i,j,k,2)+SR(4)*this%SC(i,j,k,4)+SR(5)*this%SC(i,j,k,5))
                  ! yz tensor component
                  resSC(i,j,k,5)=resSC(i,j,k,5)-this%affinecoeff*((SR(3)*this%SC(i,j,k,2)+SR(5)*this%SC(i,j,k,4)+SR(6)*this%SC(i,j,k,5))+&
                  &                                               (SR(2)*this%SC(i,j,k,3)+SR(4)*this%SC(i,j,k,5)+SR(5)*this%SC(i,j,k,6)))
                  ! zz tensor component
                  resSC(i,j,k,6)=resSC(i,j,k,6)-this%affinecoeff*2.0_WP*(SR(3)*this%SC(i,j,k,3)+SR(5)*this%SC(i,j,k,5)+SR(6)*this%SC(i,j,k,6))
               end do
            end do
         end do
      end select
   end subroutine get_CgradU
   
   
   !> Get CgradU source terms to add to multiscalar residual
   !> Assumes scalar being transported is ln(C)
   subroutine get_CgradU_log(this,gradu,resSC,time_step)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in)    :: gradU
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:),    intent(inout) :: resSC
      integer, intent(in) :: time_step
      integer :: i,j,k
      ! Temp scalar values for matrix multiplication
      real(WP), dimension(3,3) :: tmpMat,M,B,Omega   !< Matrices for diagonalization 
      real(WP), dimension(3,3) :: D                  !< Strain rate tensor for diagonalization in PTT models
      real(WP) :: omega_xy,omega_xz,omega_yz         !< Components for anti-symetric matrix
      resSC=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Zero out diagonaliztion matrices
               B=0.0_WP; Omega=0.0_WP; D=0.0_WP; M=0.0_WP; tmpMat=0.0_WP
               omega_xy=0.0_WP; omega_xz=0.0_WP; omega_yz=0.0_WP
               ! Check  C is proportional to I based upon C's eigenvalues (i.e., Lambda_ii=Lambda_jj)
               if (abs(this%eigenval(1,i,j,k)-this%eigenval(2,i,j,k)).le.1.0e-15_WP.or.abs(this%eigenval(2,i,j,k)-this%eigenval(3,i,j,k)).le.1.0e-15_WP) then
                  !>Set B equal to the strain rate tensor
                  B(1,1)=gradU(1,1,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP; B(1,2)=0.5_WP*(gradU(1,2,i,j,k)+gradU(2,1,i,j,k));                                   B(1,3)=0.5_WP*(gradU(1,3,i,j,k)+gradU(3,1,i,j,k))
                  B(2,1)=0.5_WP*(gradU(1,2,i,j,k)+gradU(2,1,i,j,k));                                   B(2,2)=gradU(2,2,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP; B(2,3)=0.5_WP*(gradU(2,3,i,j,k)+gradU(3,2,i,j,k))
                  B(3,1)=0.5_WP*(gradU(1,3,i,j,k)+gradU(3,1,i,j,k));                                   B(3,2)=0.5_WP*(gradU(2,3,i,j,k)+gradU(3,2,i,j,k));                                   B(3,3)=gradU(3,3,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP  
               else
                  select case (this%model)
                  case (lptt,eptt)
                     !>Strain rate tensor used for calculating L
                     D(1,1)=gradU(1,1,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP; D(1,2)=0.5_WP*(gradU(1,2,i,j,k)+gradU(2,1,i,j,k));                                   D(1,3)=0.5_WP*(gradU(1,3,i,j,k)+gradU(3,1,i,j,k))
                     D(2,1)=0.5_WP*(gradU(1,2,i,j,k)+gradU(2,1,i,j,k));                                   D(2,2)=gradU(2,2,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP; D(2,3)=0.5_WP*(gradU(2,3,i,j,k)+gradU(3,2,i,j,k))
                     D(3,1)=0.5_WP*(gradU(1,3,i,j,k)+gradU(3,1,i,j,k));                                   D(3,2)=0.5_WP*(gradU(2,3,i,j,k)+gradU(3,2,i,j,k));                                   D(3,3)=gradU(3,3,i,j,k)-(gradU(1,1,i,j,k)+gradU(2,2,i,j,k)+gradU(3,3,i,j,k))/3.0_WP
                     ! Use L tensor L=gradU^T-zeta*D for PTT models
                     M=matmul(transpose(this%eigenvec(:,:,i,j,k)),matmul((transpose(gradU(:,:,i,j,k))-this%affinecoeff*D),this%eigenvec(:,:,i,j,k)))
                  case default
                     M=matmul(transpose(this%eigenvec(:,:,i,j,k)),matmul(transpose(gradU(:,:,i,j,k)),this%eigenvec(:,:,i,j,k)))
                  end select
                  ! Temp matrix for calculating B
                  tmpMat=reshape((/ M(1,1),0.0_WP,0.0_WP,0.0_WP,M(2,2),0.0_WP,0.0_WP,0.0_WP,M(3,3) /),shape(tmpMat))
                  ! Form symmetric extension component of confomration tensor (B=R*{{mxx,0,0},{0,myy,0},{0,0,mzz}}*R^T={{Bxx,Bxy,Bxz},{Bxy,Byy,Byz},{Bxz,Byz,Bzz}})
                  B=matmul(this%eigenvec(:,:,i,j,k),matmul(tmpMat,transpose(this%eigenvec(:,:,i,j,k))))
                  ! Antisymmetric components
                  omega_xy=(this%eigenval(2,i,j,k)*M(1,2)+this%eigenval(1,i,j,k)*M(2,1))/(this%eigenval(2,i,j,k)-this%eigenval(1,i,j,k)) 
                  omega_xz=(this%eigenval(3,i,j,k)*M(1,3)+this%eigenval(1,i,j,k)*M(3,1))/(this%eigenval(3,i,j,k)-this%eigenval(1,i,j,k))
                  omega_yz=(this%eigenval(3,i,j,k)*M(2,3)+this%eigenval(2,i,j,k)*M(3,2))/(this%eigenval(3,i,j,k)-this%eigenval(2,i,j,k))
                  ! Temp matrix for calculating Omega
                  tmpMat=0.0_WP
                  tmpMat=reshape((/ 0.0_WP,-omega_xy,-omega_xz,omega_xy,0.0_WP,-omega_yz,omega_xz,omega_yz,0.0_WP /),shape(tmpMat))
                  ! Form rotation component of conformation tensor (Omega=R*{{0,omega_xy,omega_xz},{-omega_xy,0,omega_yz},{-omega_xz,-omega_yz,0}}*R^T={{Omegaxx,Omegaxy,Omegaxz},{Omegayx,Omegayy,Omegayz},{Omegazx,Omegazy,Omegazz}})
                  Omega=matmul(this%eigenvec(:,:,i,j,k),matmul(tmpMat,transpose(this%eigenvec(:,:,i,j,k))))
               end if
               ! Add extension and rotation components to resSC (Omega*log(C)-log(C)*Omega+2B)
               !>xx tensor component
               resSC(i,j,k,1)=2.00_WP*B(1,1)+Omega(1,2)*this%SC(i,j,k,2)-Omega(2,1)*this%SC(i,j,k,2)+Omega(1,3)*this%SC(i,j,k,3)-Omega(3,1)*this%SC(i,j,k,3)
               !>xy tensor component
               resSC(i,j,k,2)=2.00_WP*B(2,1)+Omega(2,1)*this%SC(i,j,k,1)-Omega(1,1)*this%SC(i,j,k,2)+Omega(2,2)*this%SC(i,j,k,2)+Omega(2,3)*this%SC(i,j,k,3)-Omega(2,1)*this%SC(i,j,k,4)-Omega(3,1)*this%SC(i,j,k,5)
               !>xz tensor component
               resSC(i,j,k,3)=2.00_WP*B(3,1)+Omega(3,1)*this%SC(i,j,k,1)+Omega(3,2)*this%SC(i,j,k,2)-Omega(1,1)*this%SC(i,j,k,3)+Omega(3,3)*this%SC(i,j,k,3)-Omega(2,1)*this%SC(i,j,k,5)-Omega(3,1)*this%SC(i,j,k,6)
               !>yy tensor component
               resSC(i,j,k,4)=2.00_WP*B(2,2)-Omega(1,2)*this%SC(i,j,k,2)+Omega(2,1)*this%SC(i,j,k,2)+Omega(2,3)*this%SC(i,j,k,5)-Omega(3,2)*this%SC(i,j,k,5)
               !>yz tensor component
               resSC(i,j,k,5)=2.00_WP*B(3,2)+Omega(3,1)*this%SC(i,j,k,2)-Omega(1,2)*this%SC(i,j,k,3)+Omega(3,2)*this%SC(i,j,k,4)-Omega(2,2)*this%SC(i,j,k,5)+Omega(3,3)*this%SC(i,j,k,5)-Omega(3,2)*this%SC(i,j,k,6)
               !>zz tensor component
               resSC(i,j,k,6)=2.00_WP*B(3,3)-Omega(1,3)*this%SC(i,j,k,3)+Omega(3,1)*this%SC(i,j,k,3)-Omega(2,3)*this%SC(i,j,k,5)+Omega(3,2)*this%SC(i,j,k,5)
            end do 
         end do
      end do
   end subroutine get_CgradU_log
   
   
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

   !> Add viscoelastic relaxation source based on semi analtical integration
   subroutine get_relax_analytical(this,dt)
      use messager, only: die
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: i,j,k
      real(WP) :: f,coeff
      select case (this%model)
      case (oldroydb) ! Add relaxation source for Oldroyd-B (1/t_relax)(C-I)
         coeff=-1.00_WP*(dt/this%trelax)  
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  this%SCrec(i,j,k,1)=this%SCrec(i,j,k,1)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< xx tensor component
                  this%SCrec(i,j,k,2)=this%SCrec(i,j,k,2)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xy tensor component
                  this%SCrec(i,j,k,3)=this%SCrec(i,j,k,3)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xz tensor component
                  this%SCrec(i,j,k,4)=this%SCrec(i,j,k,4)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< yy tensor component
                  this%SCrec(i,j,k,5)=this%SCrec(i,j,k,5)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< yz tensor component
                  this%SCrec(i,j,k,6)=this%SCrec(i,j,k,6)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< zz tensor component
               end do
            end do
         end do
      case (eptt) ! Add relaxation source for ePTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  ! f=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*((this%SCold(i,j,k,1)+this%SCold(i,j,k,4)+this%SCold(i,j,k,6))-3.0_WP))
                  ! coeff=-1.00_WP*((f*dt)/this%trelax)    
                  ! this%SC(i,j,k,1)=this%SC(i,j,k,1)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< xx tensor component
                  ! this%SC(i,j,k,2)=this%SC(i,j,k,2)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xy tensor component
                  ! this%SC(i,j,k,3)=this%SC(i,j,k,3)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xz tensor component
                  ! this%SC(i,j,k,4)=this%SC(i,j,k,4)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< yy tensor component
                  ! this%SC(i,j,k,5)=this%SC(i,j,k,5)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< yz tensor component
                  ! this%SC(i,j,k,6)=this%SC(i,j,k,6)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< zz tensor component
                  f=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*((this%SCrecold(i,j,k,1)+this%SCrecold(i,j,k,4)+this%SCrecold(i,j,k,6))-3.0_WP))
                  coeff=-1.00_WP*((f*dt)/this%trelax)    
                  this%SCrec(i,j,k,1)=this%SCrec(i,j,k,1)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< xx tensor component
                  this%SCrec(i,j,k,2)=this%SCrec(i,j,k,2)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xy tensor component
                  this%SCrec(i,j,k,3)=this%SCrec(i,j,k,3)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< xz tensor component
                  this%SCrec(i,j,k,4)=this%SCrec(i,j,k,4)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< yy tensor component
                  this%SCrec(i,j,k,5)=this%SCrec(i,j,k,5)*exp(coeff)+(1.00_WP-exp(coeff))*0.0_WP !< yz tensor component
                  this%SCrec(i,j,k,6)=this%SCrec(i,j,k,6)*exp(coeff)+(1.00_WP-exp(coeff))*1.0_WP !< zz tensor component
               end do
            end do
         end do
      case default
         call die('[tpviscoelastic get_relax] Unknown viscoelastic model selected')
      end select
   end subroutine get_relax_analytical


   !> Add viscoelastic relaxation source using log-conformation stabilization 
   !> Assumes scalar being transported is ln(C)
   subroutine get_relax_log(this,resSC)
      use messager, only: die
      implicit none
      class(tpviscoelastic), intent(inout) :: this
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
                  trace=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,j)**2+this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
                  !>Relaxation function coefficent
                  coeff=(this%Lmax**2-3.00_WP)/(this%Lmax**2-trace)
                  !>Add source term to residual
                  resSC(i,j,k,1)=(1.00_WP/this%trelax)*(this%eigenvec(1,1,i,j,k)**2                      *((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)**2                      *((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)**2                      *((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=(1.00_WP/this%trelax)*(this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)*((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)*((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)*((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=(1.00_WP/this%trelax)*(this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=(1.00_WP/this%trelax)*(this%eigenvec(2,1,i,j,k)**2                      *((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)**2                      *((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)**2                      *((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=(1.00_WP/this%trelax)*(this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=(1.00_WP/this%trelax)*(this%eigenvec(3,1,i,j,k)**2                      *((coeff/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(3,2,i,j,k)**2                      *((coeff/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(3,3,i,j,k)**2                      *((coeff/this%eigenval(3,i,j,k))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (fenecr) ! Add relaxation source for FENE-CR 
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,j)**2+this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
                  !>Relaxation function coefficent
                  coeff=1.00_WP/(1.0_WP-trace/this%Lmax**2)
                  coeff=coeff/this%trelax
                  !>Add source term to residual        
                  resSC(i,j,k,1)=coeff*(this%eigenvec(1,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=coeff*(this%eigenvec(2,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=coeff*(this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=coeff*(this%eigenvec(3,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(3,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(3,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (oldroydb) ! Add relaxation source for Oldroyd-B 
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  coeff=1.00_WP/this%trelax
                  resSC(i,j,k,1)=coeff*(this%eigenvec(1,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=coeff*(this%eigenvec(2,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=coeff*(this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=coeff*(this%eigenvec(3,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(3,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(3,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case (lptt) ! Add relaxation source for lPTT
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle                              !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,j)**2+this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
                  !>Relaxation function coefficent
                  coeff=1.00_WP+(this%elongvisc/(1.0_WP-this%affinecoeff))*(trace-3.0_WP)
                  coeff=coeff/this%trelax
                  !>Add source term to residual
                  resSC(i,j,k,1)=coeff*(this%eigenvec(1,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=coeff*(this%eigenvec(2,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=coeff*(this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=coeff*(this%eigenvec(3,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(3,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(3,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< zz tensor component
               
               end do
            end do
         end do
      case (eptt) ! Add relaxation source for ePTT model
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (this%mask(i,j,k).ne.0) cycle !< Skip non-solved cells
                  !>Trace of reconstructed conformation tensor
                  trace=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2+&
                  &     this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,j)**2+this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2+this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
                  !>Relaxation function coefficent
                  coeff=exp(this%elongvisc/(1.0_WP-this%affinecoeff)*(trace-3.0_WP))
                  coeff=coeff/this%trelax
                  !>Add source term to residual
                  resSC(i,j,k,1)=coeff*(this%eigenvec(1,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xx tensor component
                  resSC(i,j,k,2)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xy tensor component
                  resSC(i,j,k,3)=coeff*(this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< xz tensor component
                  resSC(i,j,k,4)=coeff*(this%eigenvec(2,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yy tensor component
                  resSC(i,j,k,5)=coeff*(this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)*((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)*((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)*((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< yz tensor component
                  resSC(i,j,k,6)=coeff*(this%eigenvec(3,1,i,j,k)**2                      *((1.00_WP/this%eigenval(1,i,j,k))-1.00_WP)+this%eigenvec(3,2,i,j,k)**2                      *((1.00_WP/this%eigenval(2,i,j,k))-1.00_WP)+this%eigenvec(3,3,i,j,k)**2                      *((1.00_WP/this%eigenval(3,i,j,k))-1.00_WP))  !< zz tensor component
               end do
            end do
         end do
      case default
         call die('[tpviscoelastic get_relax] Unknown viscoelastic model selected')
      end select
   end subroutine get_relax_log

   !> Calculate the this%eigenval and eigenvectors of the conformation tensor for the whole domain
   !> Assumes scalar being transported is ln(C)
   subroutine get_eigensystem(this,VF)
      use mathtools, only: eigensolve3
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: VF
      integer :: i,j,k
      real(WP), dimension(3,3) :: A
      real(WP) :: trace
      ! First ensure storage is allocated
      ! if (.not.allocated(this%eigenval)) allocate(this%eigenval    (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! if (.not.allocated(this%eigenvec)) allocate(this%eigenvec(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Empty storage
      this%eigenval=0.0_WP
      this%eigenvec=0.0_WP
      ! Loop over the domain
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (VF(i,j,k).gt.0.0_WP) then
                  A=0.0_WP
                  ! Form local matrix to diagonalize
                  A(1,1)=this%SC(i,j,k,1); A(1,2)=this%SC(i,j,k,2); A(1,3)=this%SC(i,j,k,3)
                  A(2,1)=this%SC(i,j,k,2); A(2,2)=this%SC(i,j,k,4); A(2,3)=this%SC(i,j,k,5)
                  A(3,1)=this%SC(i,j,k,3); A(3,2)=this%SC(i,j,k,5); A(3,3)=this%SC(i,j,k,6)
                  trace=A(1,1)+A(2,2)+A(3,3)
                  if (A(2,1).le.1e-15.and.A(3,1).le.1e-15.and.A(3,2).le.1e-15) then
                  ! >If C is propotional to I, set eigenvalues and eigenvectors to that of the identifty tensor
                  ! Eigenvectors
                  this%eigenvec(1,1,i,j,k)=1.00_WP; this%eigenvec(1,2,i,j,k)=0.00_WP; this%eigenvec(1,3,i,j,k)=0.00_WP
                  this%eigenvec(2,1,i,j,k)=0.00_WP; this%eigenvec(2,2,i,j,k)=1.00_WP; this%eigenvec(2,3,i,j,k)=0.00_WP
                  this%eigenvec(3,1,i,j,k)=0.00_WP; this%eigenvec(3,2,i,j,k)=0.00_WP; this%eigenvec(3,3,i,j,k)=1.00_WP
                  ! Eigenvalues
                  this%eigenval(1,i,j,k)=1.00_WP
                  this%eigenval(2,i,j,k)=1.00_WP
                  this%eigenval(3,i,j,k)=1.00_WP
                  else
                  ! Diagonalize it
                     call eigensolve3(A,this%eigenvec(:,:,i,j,k),this%eigenval(:,i,j,k))
                     ! Take the exponential
                     this%eigenval(:,i,j,k)=exp(this%eigenval(:,i,j,k))
                     ! ! Scale by alpha
                     ! this%eigenvec(:,:,i,j,k)=VF(i,j,k)*this%eigenvec(:,:,i,j,k)
                     ! this%eigenval(:,i,j,k)=VF(i,j,k)*this%eigenval(:,i,j,k)
                  end if
               end if
            end do
         end do
      end do
   end subroutine get_eigensystem

   ! !> Calculate the this%eigenval and eigenvectors of the conformation tensor for the whole domain
   ! subroutine get_eigensystem(this,Atmp)
   !    use mathtools, only: eigensolve3
   !    implicit none
   !    class(tpviscoelastic), intent(inout) :: this
   !    real(WP), dimension(1:,1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: Atmp
   !    integer :: i,j,k
   !    ! real(WP), dimension(3,3) :: A
   !    ! First ensure storage is allocated
   !    ! if (.not.allocated(this%eigenval)) allocate(this%eigenval    (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
   !    ! if (.not.allocated(this%eigenvec)) allocate(this%eigenvec(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
   !    ! Empty storage
   !    this%eigenval=0.0_WP
   !    this%eigenvec=0.0_WP
   !    ! Loop over the domain
   !    do k=this%cfg%kmino_,this%cfg%kmaxo_
   !       do j=this%cfg%jmino_,this%cfg%jmaxo_
   !          do i=this%cfg%imino_,this%cfg%imaxo_
   !             ! Form local matrix to diagonalize
   !             ! A(1,1)=this%SC(i,j,k,1); A(1,2)=this%SC(i,j,k,2); A(1,3)=this%SC(i,j,k,3)
   !             ! A(2,1)=this%SC(i,j,k,2); A(2,2)=this%SC(i,j,k,4); A(2,3)=this%SC(i,j,k,5)
   !             ! A(3,1)=this%SC(i,j,k,3); A(3,2)=this%SC(i,j,k,5); A(3,3)=this%SC(i,j,k,6)
   !             ! Check if A is proprtional to I
   !             ! if (Atmp(2,1,i,j,k).le.1e-15.and.Atmp(3,1,i,j,k).le.1e-15.and.Atmp(3,2,i,j,k).le.1e-15) then
   !                !>If C is propotional to I, set eigenvalues and eigenvectors to that of the identifty tensor
   !                ! Eigenvectors
   !                ! this%eigenvec(1,1,i,j,k)=1.00_WP; this%eigenvec(1,2,i,j,k)=0.00_WP; this%eigenvec(1,3,i,j,k)=0.00_WP
   !                ! this%eigenvec(2,1,i,j,k)=0.00_WP; this%eigenvec(2,2,i,j,k)=1.00_WP; this%eigenvec(2,3,i,j,k)=0.00_WP
   !                ! this%eigenvec(3,1,i,j,k)=0.00_WP; this%eigenvec(3,2,i,j,k)=0.00_WP; this%eigenvec(3,3,i,j,k)=1.00_WP
   !                ! Eigenvalues
   !                ! this%eigenval(1,i,j,k)=1.00_WP
   !                ! this%eigenval(2,i,j,k)=1.00_WP
   !                ! this%eigenval(3,i,j,k)=1.00_WP
   !             ! else
   !                ! Diagonalize it
   !                call eigensolve3(Atmp,this%eigenvec(:,:,i,j,k),this%eigenval(:,i,j,k))
   !                ! Take the exponential
   !                ! this%eigenval(:,i,j,k)=exp(this%eigenval(:,i,j,k))
   !                ! this%eigenval(:,i,j,k)=this%eigenval(:,i,j,k)
   !             ! end if
   !          end do
   !       end do
   !    end do
   ! end subroutine get_eigensystem

   
   !> Calculate the this%eigenval and eigenvectors of the conformation tensor for the whole domain
   !> Assumes scalar being transported is C
   subroutine get_eigensystem_SCrec(this)
      use mathtools, only: eigensolve3
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      integer :: i,j,k
      real(WP), dimension(3,3) :: A
      real(WP) :: trace
      ! First ensure storage is allocated
      ! if (.not.allocated(this%eigenval)) allocate(this%eigenval    (1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! if (.not.allocated(this%eigenvec)) allocate(this%eigenvec(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! Empty storage
      this%eigenval=0.0_WP
      this%eigenvec=0.0_WP
      ! Loop over full domain
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Form local matrix to diagonalize
               A=0.0_WP
               A(1,1)=this%SCrec(i,j,k,1); A(1,2)=this%SCrec(i,j,k,2); A(1,3)=this%SCrec(i,j,k,3)
               A(2,1)=this%SCrec(i,j,k,2); A(2,2)=this%SCrec(i,j,k,4); A(2,3)=this%SCrec(i,j,k,5)
               A(3,1)=this%SCrec(i,j,k,3); A(3,2)=this%SCrec(i,j,k,5); A(3,3)=this%SCrec(i,j,k,6)
               trace=A(1,1)+A(2,2)+A(3,3)
               if (A(2,1).le.1e-15.and.A(3,1).le.1e-15.and.A(3,2).le.1e-15.and.trace.le.3.10_WP.and.trace.ge.2.90_WP) then
               else
               ! Diagonalize it
               call eigensolve3(A,this%eigenvec(:,:,i,j,k),this%eigenval(:,i,j,k))
               end if
            end do
         end do
      end do
   end subroutine get_eigensystem_SCrec

   !> Reconstruct the conformation tensor for its decomposed eigenvalues and eigenvectors
   subroutine reconstruct_conformation(this)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      integer :: i,j,k
      ! First ensure storage is allocated
      ! if (.not.allocated(this%SCrec))    allocate(this%SCrec(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6)) 
      ! if (.not.allocated(this%SCrecold)) allocate(this%SCrecold(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:6)) 
      this%SCrec=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! Reconstruct conformation tensor (C=R*exp(ln(Lambda))*R^T={{Cxx,Cxy,Cxz},{Cxy,Cyy,Cyz},{Cxz,Cyz,Czz}})
               !>xx tensor component
               this%SCrec(i,j,k,1)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2
               !>xy tensor component
               this%SCrec(i,j,k,2)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)
               !>xz tensor component
               this%SCrec(i,j,k,3)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)
               !>yy tensor component
               this%SCrec(i,j,k,4)=this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2
               !>yz tensor component
               this%SCrec(i,j,k,5)=this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)
               !>zz tensor component
               this%SCrec(i,j,k,6)=this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
            end do
         end do
      end do
   end subroutine reconstruct_conformation

   !> Reconstruct the log conformation tensor for its decomposed eigenvalues and eigenvectors
   subroutine reconstruct_log_conformation(this)
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      integer :: i,j,k
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               ! Skip non-solved cells
               if (this%mask(i,j,k).ne.0) cycle
               ! ! Converteigenvalues into log
               ! this%eigenval(:,i,j,k)=log(this%eigenval(:,i,j,k))
               ! Reconstruct conformation tensor (C=R*exp(ln(Lambda))*R^T={{Cxx,Cxy,Cxz},{Cxy,Cyy,Cyz},{Cxz,Cyz,Czz}})
               !>xx tensor component
               this%SC(i,j,k,1)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)**2
               !>xy tensor component
               this%SC(i,j,k,2)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)*this%eigenvec(2,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)*this%eigenvec(2,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)*this%eigenvec(2,3,i,j,k)
               !>xz tensor component
               this%SC(i,j,k,3)=this%eigenval(1,i,j,k)*this%eigenvec(1,1,i,j,k)*this%eigenvec(3,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(1,2,i,j,k)*this%eigenvec(3,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(1,3,i,j,k)*this%eigenvec(3,3,i,j,k)
               !>yy tensor component
               this%SC(i,j,k,4)=this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)**2
               !>yz tensor component
               this%SC(i,j,k,5)=this%eigenval(1,i,j,k)*this%eigenvec(2,1,i,j,k)*this%eigenvec(3,1,i,j,k)+this%eigenval(2,i,j,k)*this%eigenvec(2,2,i,j,k)*this%eigenvec(3,2,i,j,k)+this%eigenval(3,i,j,k)*this%eigenvec(2,3,i,j,k)*this%eigenvec(3,3,i,j,k)
               !>zz tensor component
               this%SC(i,j,k,6)=this%eigenval(1,i,j,k)*this%eigenvec(3,1,i,j,k)**2                      +this%eigenval(2,i,j,k)*this%eigenvec(3,2,i,j,k)**2                      +this%eigenval(3,i,j,k)*this%eigenvec(3,3,i,j,k)**2
            end do
         end do
      end do
   end subroutine reconstruct_log_conformation

   !> Calculate the min, max, and int of our reconstructed C field
   subroutine get_max_reconstructed(this,VF)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: VF        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: ierr,nsc
      real(WP) :: my_SCmax,my_SCmin
      real(WP), dimension(:,:,:), allocatable :: tmp
      allocate(tmp(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      ! First ensure storage is allocated
      if (.not.allocated(this%SCrecint)) allocate(this%SCrecint(1:6)) 
      if (.not.allocated(this%SCrecmin)) allocate(this%SCrecmin(1:6)) 
      if (.not.allocated(this%SCrecmax)) allocate(this%SCrecmax(1:6)) 
      do nsc=1,this%nscalar
         my_SCmax=maxval(this%SCrec(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmax,this%SCrecmax(nsc),1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         my_SCmin=minval(this%SCrec(:,:,:,nsc)); call MPI_ALLREDUCE(my_SCmin,this%SCrecmin(nsc),1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         if      (this%phase(nsc).eq.0) then ! Liquid scalar
            tmp=this%SCrec(:,:,:,nsc)*(       VF(:,:,:))
         else if (this%phase(nsc).eq.1) then ! Gas scalar
            tmp=this%SCrec(:,:,:,nsc)*(1.0_WP-VF(:,:,:))
         end if
         call this%cfg%integrate(A=tmp,integral=this%SCrecint(nsc))
      end do
      deallocate(tmp)
   end subroutine get_max_reconstructed

   !> Calculate the max of our fields
   subroutine get_max_eig(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(tpviscoelastic), intent(inout) :: this
      integer :: i,j,k,ierr
      real(WP) :: my_Eigval1_max,my_Eigval2_max,my_Eigval3_max
      
      ! Set all to zero
      my_Eigval1_max=0.0_WP; my_Eigval2_max=0.0_WP; my_Eigval3_max=0.0_WP;
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Eigval1_max=max(my_Eigval1_max,abs(this%eigenval(1,i,j,k)))
               my_Eigval2_max=max(my_Eigval2_max,abs(this%eigenval(2,i,j,k)))
               my_Eigval3_max=max(my_Eigval3_max,abs(this%eigenval(3,i,j,k)))
            end do
         end do
      end do
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Eigval1_max  ,this%eigenval(1,i,j,k)  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Eigval2_max  ,this%eigenval(2,i,j,k)  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Eigval3_max  ,this%eigenval(3,i,j,k)  ,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_max_eig
   
   

   
   
end module tpviscoelastic_class