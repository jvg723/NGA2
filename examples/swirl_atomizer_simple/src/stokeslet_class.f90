!> Smooth particle hydrodynamics solver
module stokeslet_class
    use precision,      only: WP
    use string,         only: str_medium
    use config_class,   only: config
    use mpi_f08,        only: MPI_Datatype,MPI_INTEGER,MPI_DOUBLE_PRECISION
    implicit none
    private
 
 
    ! Expose type/constructor/methods
    public :: stokeslet
 
 
    ! List of known available smoothing kernels
    integer, parameter, public :: cubic=1
    integer, parameter, public :: quintic=2
    integer, parameter, public :: gaussian=3
    integer, parameter, public :: cubicBS=4
    integer, parameter, public :: wendland=5
 
    !> Memory adaptation parameter
    real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
    real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
 
 
    !> I/O chunk size to read at a time
    integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing
 
 
    !> Particle definition
    type :: part
       !> MPI_DOUBLE_PRECISION data
       real(WP) :: fcoeff                     !< Forcing coefficient
       real(WP) :: VF
       real(WP), dimension(3) :: pos          !< Particle center coordinates
       real(WP), dimension(3) :: nedge        !< Edge retraction orientation
       !> MPI_INTEGER data
       integer :: id                          !< ID the object is associated with
       integer , dimension(3) :: ind          !< Index of cell containing particle center
       integer :: flag                        !< Control parameter (0=normal, 1=done->will be removed)
    end type part
    !> Number of blocks, block length, and block types in a particle
    integer, parameter                         :: part_nblock=2
    integer           , dimension(part_nblock) :: part_lblock=[8,5]
    type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_DOUBLE_PRECISION,MPI_INTEGER]
    !> MPI_PART derived datatype and size
    type(MPI_Datatype) :: MPI_PART
    integer :: MPI_PART_SIZE
 
    !> Smooth particle hydrodynamics object definition
    type :: stokeslet
 
       ! This config is used for parallelization and for neighbor detection
       class(config), pointer :: cfg
 
       ! This is the name of the solver
       character(len=str_medium) :: name='UNNAMED_SPH'
 
       ! Kernel properties
       real(WP) :: a                                       !< Forcing range ratio
       real(WP) :: b                                       !< Modified regularized ratio
 
       ! Global and local particle data
       integer :: np                                       !< Global number of particles
       integer :: np_                                      !< Local number of particles
       integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
       type(part), dimension(:), allocatable :: p          !< Array of particles of type part
 
 
    contains
       procedure :: initialize                        !< Initialize the stokeslet object
       procedure :: update_partmesh                   !< Update a partmesh object using current particles
    !    procedure :: prepare_mpi_part
       procedure :: resize                            !< Resize particle array to given size
    !    procedure :: sync                              !< Synchronize particles across interprocessor boundaries
    !    procedure :: update_periodic                   !< Account for periodic conditions
    end type stokeslet
 
 contains
 

 
 
    !> Default constructor for SPH solver
    ! function initialize(cfg,name) result(self)
    subroutine initialize(this,cfg,name)
       implicit none
       class(stokeslet) :: this
       class(config), target, intent(in) :: cfg
       character(len=*), optional :: name
       integer :: i,j,k,l
       ! Set the name for the solver
       if (present(name)) this%name=trim(adjustl(name))
 
       ! Point to pgrid object
       this%cfg=>cfg
 
       ! Allocate variables
       allocate(this%np_proc(1:this%cfg%nproc)); this%np_proc=0
       this%np_=0; this%np=0
       call this%resize(0)
 
       ! Initialize MPI derived datatype for a particle
    !    call prepare_mpi_part()
 
       ! Log/screen output
       logging: block
          use, intrinsic :: iso_fortran_env, only: output_unit
          use param,    only: verbose
          use messager, only: log
          use string,   only: str_long
          character(len=str_long) :: message
          if (this%cfg%amRoot) then
             write(message,'("Stokeslet object [",a,"] on partitioned grid [",a,"]")') trim(this%name),trim(this%cfg%name)
             if (verbose.gt.1) write(output_unit,'(a)') trim(message)
             if (verbose.gt.0) call log(message)
          end if
       end block logging
 
    end subroutine initialize
 
 
    !> Update particle mesh using our current particles
    subroutine update_partmesh(this,pmesh)
       use partmesh_class, only: partmesh
       implicit none
       class(stokeslet), intent(inout) :: this
       class(partmesh), intent(inout) :: pmesh
       integer :: i
       ! Reset particle mesh storage
       call pmesh%reset()
       ! Nothing else to do if no particle is present
       if (this%np_.eq.0) return
       ! Copy particle info
       call pmesh%set_size(this%np_)
       do i=1,this%np_
          pmesh%pos(:,i)=this%p(i)%pos
       end do
    end subroutine update_partmesh
 
 
    !> Adaptation of particle array size
    subroutine resize(this,n)
       implicit none
       class(stokeslet), intent(inout) :: this
       integer, intent(in) :: n
       type(part), dimension(:), allocatable :: tmp
       integer :: size_now,size_new
       ! Resize particle array to size n
       if (.not.allocated(this%p)) then
          ! Allocate directly to size n
          allocate(this%p(n))
          this%p(1:n)%flag=1
       else
          ! Update from a non-zero size to another non-zero size
          size_now=size(this%p,dim=1)
          if (n.gt.size_now) then
             size_new=max(n,int(real(size_now,WP)*coeff_up))
             allocate(tmp(size_new))
             tmp(1:size_now)=this%p
             tmp(size_now+1:)%flag=1
             call move_alloc(tmp,this%p)
          else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
             allocate(tmp(n))
             tmp(1:n)=this%p(1:n)
             call move_alloc(tmp,this%p)
          end if
       end if
    end subroutine resize
 
 
 
 end module stokeslet_class
 