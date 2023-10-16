!> Inselastic fluid model class
!> Calculates properties such as shear rate dependent viscosity,...
module inelastic_class
   use config_class,      only: config
   use precision,         only: WP
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: inelastic
   
   ! List of available inelastic models
   integer, parameter, public :: powerlaw=0                !< Power law fluid model
   integer, parameter, public :: carreau =1                !< Carreau shear thinning fluid

   
   !> Inelastic fluid object definition
   type :: inelastic
      ! Model parameters
      integer  :: model                                    !< inelastic model
      real(WP) :: ncoeff                                   !< Powerlaw coefficient
      real(WP) :: tchar                                    !< Characteristic fluid timescale
      real(WP) :: visc_zero                                !< Zero shear rate fluid viscosity
      real(WP) :: flow_index                               !< Flow consistency index of the fluid for a power law fluid
      ! Fluid viscosity
      real(WP), dimension(:,:,:), allocatable :: visc      !< Rate dependent fluid viscosity
      ! SRmag
      real(WP), dimension(:,:,:), allocatable :: SRmag     !< Strain rate magnitude
      ! Monitoring info
      real(WP) :: visc_min,visc_max                        !< Min and max fluid viscosity
   contains
      procedure :: update_visc                             !< Update fluid viscsity given strain rate tensor
      procedure :: get_max                                 !< Monitoring for inelastic class
   end type inelastic
   
   !> Declare inelastic fluid model constructor
   interface inelastic
      procedure constructor
   end interface inelastic
   
   
contains
   
   
   !> Inelastic fluid model constructor from multiscalar
   function constructor(cfg,model,name) result(self)
      implicit none
      type(inelastic) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: model
      integer, intent(in) :: scheme
      character(len=*), optional :: name
      ! Assign closure model for inelastic fluid
      self%model=model
      ! Allocate storage for rate dependent viscosity and SR magnitude
      allocate(self%visc (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%visc =0.0_WP
      allocate(self%SRmag(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%SRmag=0.0_WP
   end function constructor
   
   
   !> Compute visc from SR
   subroutine update_visc(this,SR)
      implicit none
      class(inelastic), intent(inout) :: this
      real(WP), dimension(1:,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: SR
      integer :: i,j,k
      select case (this%model)
      case (powerlaw) ! Power law fluid
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Compute second invariant of strain rate tensor = sqrt(2*SR**2)
                  this%SRmag(i,j,k)=sqrt(2.0_WP*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2)))
                  ! Compute viscosity
                  this%visc(i,j,k)=this%flow_index*this%SRmag(i,j,k)**this%ncoeff
               end do
            end do
         end do
      case (carreau) ! Carreau shear thinning fluid
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Compute second invariant of strain rate tensor = sqrt(2*SR**2)
                  this%SRmag(i,j,k)=sqrt(2.0_WP*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2)))
                  ! Compute viscosity
                  this%visc(i,j,k)=this%visc_zero*(1.0_WP+(this%tchar*this%SRmag(i,j,k))**2)**(0.5_WP*this%ncoeff-0.5_WP)
               end do
            end do
         end do
      case default
         call die('[inelastic update_visc] Unknown rate depdenent viscosity model selected')
      end select
   end subroutine update_visc
   

   !> Calculate the min and max of our viscosity field
   subroutine get_max(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(inelastic), intent(inout) :: this
      integer :: ierr
      real(WP) :: my_visc_pmax,my_visc_pmin
      my_visc_pmax=maxval(this%visc_p); call MPI_ALLREDUCE(my_visc_pmax,this%visc_pmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_visc_pmin=minval(this%visc_p); call MPI_ALLREDUCE(my_visc_pmin,this%visc_pmin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine get_max


end module inelastic_class