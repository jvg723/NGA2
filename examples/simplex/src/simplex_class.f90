!> Definition for a simplex class
module simplex_class
   use precision,         only: WP
   use string,            only: str_medium
   use inputfile_class,   only: inputfile
   use ibconfig_class,    only: ibconfig
   use polygon_class,     only: polygon
   use surfmesh_class,    only: surfmesh
   use partmesh_class,    only: partmesh
   use ensight_class,     only: ensight
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use lpt_class,         only: lpt
   use cclabel_class,     only: cclabel
   use iterator_class,    only: iterator
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use pardata_class,     only: pardata
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   implicit none
   private
   
   public :: simplex
   
   !> Simplex object
   type :: simplex
      
      !> Provide a pardata and an event tracker for saving restarts
      type(event)    :: save_evt
      type(pardata)  :: df
      logical :: restarted
      
      !> Input file for the simulation
      type(inputfile) :: input
      
      !> Config with IB based on polygon
      type(polygon)  :: poly
      type(ibconfig) :: cfg
      
      !> Flow solver
      type(vfs)         :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< HYPRE linear solver for pressure
      !type(ddadi)       :: vs    !< DDADI linear solver for velocity
      type(sgsmodel)    :: sgs   !< SGS model for eddy viscosity
      type(timetracker) :: time  !< Time info
      type(cclabel)     :: ccl   !< CCLabel to transfer droplets
      
      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile     !< General simulation monitoring
      type(monitor) :: cflfile   !< CFL monitoring
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU           !< Velocity gradient
      real(WP), dimension(:,:,:,:), allocatable :: SR                !< Strain rate tensor
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib         !< IB slip velocity
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP)       :: vof_removed        !< Integral of VOF removed
      integer        :: nlayer=4           !< Size of buffer layer for VOF removal
      
      !> Timing info
      type(monitor) :: timefile     !< Timing monitoring
      type(timer)   :: tstep        !< Timer for step
      type(timer)   :: tsgs         !< Timer for SGS
      type(timer)   :: tvel         !< Timer for velocity
      type(timer)   :: tpres        !< Timer for pressure
      type(timer)   :: tvof         !< Timer for VOF
      type(timer)   :: ttrans       !< Timer for VOF transfer
      type(timer)   :: ttrans_lig   !< Timer for ligament VOF transfer

      !> Event for flow rate analysis
      type(event) :: flowrate_evt  !< Event trigger for flow rate analysis
      
      !> Drop transfer modeling
      logical :: use_drop_transfer    !< Do we use droplet transfer
      logical :: use_ligament_breakup !< Do we breakup ligaments
      type(lpt)      :: lp            !< Lagrangian particle tracking
      type(monitor)  :: pfile         !< Particle monitoring
      type(partmesh) :: pmesh         !< Particle mesh for lpt
      real(WP) :: dmax                !< Maximum diameter for transfer
      real(WP) :: dmin                !< Minimum diameter below which transfer is automatic
      real(WP) :: ddel                !< Minimum diameter below which structure is directly deleted
      real(WP) :: emax                !< Maximum eccentricity for transfer
      real(WP) :: vof_transfered      !< Integral of VOF transfered
      real(WP) :: vof_deleted         !< Integral of VOF deleted

      !> Ligament transfer
      type(cclabel) :: ccl_ligament
      real(WP), dimension(:,:,:), allocatable :: unfiltered_thickness            !< Tmp film_type for output purposes

      !> Additional CCL for transfering droplets in buffer layer
      type(cclabel) :: ccl_buffer
      
      !> Inlet pipes geometry and flow rates
      real(WP) :: Rinlet=0.002_WP
      real(WP) :: Rexit=0.00143_WP
      real(WP) :: Rpipe=0.000185_WP
      real(WP) :: Rcoflow=0.003_WP
      real(WP), dimension(3) :: p1=[-0.00442_WP,0.0_WP,+0.001245_WP]
      real(WP), dimension(3) :: p2=[-0.00442_WP,0.0_WP,-0.001245_WP]
      real(WP), dimension(3) :: n1=[+0.6_WP,-0.8_WP,0.0_WP]
      real(WP), dimension(3) :: n2=[+0.6_WP,+0.8_WP,0.0_WP]
      real(WP) :: Ucoflow,mfr,Apipe

      !> Break-up model parameters
      real(WP) :: vol_convert=0.0_WP
      real(WP), dimension(:,:,:), allocatable :: film_type            !< Tmp film_type for output purposes
      real(WP), dimension(:,:,:), allocatable :: thickness            !< Tmp film_type for output purposes
      real(WP), dimension(:,:,:), allocatable :: d1,d2,d3            !< Tmp film_type for output purposes
      real(WP), dimension(:,:,:), allocatable :: l1,l2,l3            !< Tmp film_type for output purposes
      real(WP), dimension(:,:,:), allocatable :: d3_d1,d3_d2,d2_d1            !< Tmp film_type for output purposes
      
   contains
      procedure :: init                            !< Initialize simplex simulation
      procedure :: step                            !< Advance simplex simulation by one time step
      procedure :: final                           !< Finalize simplex simulation
      procedure :: transfer_drops                  !< Transfer drops to a Lagrangian representation
      procedure :: transfer_ligaments              !< Transfer ligaments to a Lagrangian representation
      procedure :: transfer_buffer                 !< Transfer structs in the buffer layer
      procedure :: analyze_flowrate                !< Compute and output flow rate through the nozzle
      procedure :: get_thickness_unfiltered
      procedure :: get_cclstats
      procedure :: get_structminthickness
      procedure :: get_localfilmtype
   end type simplex

   ! Temp arrays for ligament transfer
   real(WP), dimension(:,:,:), allocatable :: tmpthickness
   integer, dimension(:,:,:), allocatable :: tmpfilm_type
   real(WP) :: min_ligamentthickness=1.1_WP
   
   
contains
   
   
   !> Compute and output flow rate through the nozzle
   subroutine analyze_flowrate(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      use filesys,  only: makedir,isdir
      use string,   only: str_medium
      implicit none
      class(simplex), intent(inout) :: this
      real(WP), dimension(:), allocatable :: CSA_s,CSA_f,CSA_l,CSA_g !< Solid, fluid, liquid, and gas cross-sectional areas
      real(WP), dimension(:), allocatable :: VFR_s,VFR_f,VFR_l,VFR_g !< Solid, fluid, liquid, and gas volume flow rates
      character(len=str_medium) :: filename,timestamp
      integer :: i,j,k,ierr,iunit
      real(WP) :: area,rad
      ! Allocate and initialize 1D storage for areas and flow rates
      allocate(CSA_s(this%cfg%imin:this%cfg%imax)); CSA_s=0.0_WP
      allocate(CSA_f(this%cfg%imin:this%cfg%imax)); CSA_f=0.0_WP
      allocate(CSA_l(this%cfg%imin:this%cfg%imax)); CSA_l=0.0_WP
      allocate(CSA_g(this%cfg%imin:this%cfg%imax)); CSA_g=0.0_WP
      allocate(VFR_s(this%cfg%imin:this%cfg%imax)); VFR_s=0.0_WP
      allocate(VFR_f(this%cfg%imin:this%cfg%imax)); VFR_f=0.0_WP
      allocate(VFR_l(this%cfg%imin:this%cfg%imax)); VFR_l=0.0_WP
      allocate(VFR_g(this%cfg%imin:this%cfg%imax)); VFR_g=0.0_WP
      ! Integrate various areas and flow rates
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               ! Cell cross-sectional area and radial location
               area=this%cfg%dy(j)*this%cfg%dz(k)
               rad=sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)
               ! Increment solid area and flow rate
               CSA_s(i)=CSA_s(i)+area*(1.0_WP-this%cfg%VF(i,j,k))
               VFR_s(i)=VFR_s(i)+area*(1.0_WP-this%cfg%VF(i,j,k))*this%Ui(i,j,k)
               ! Increment fluid areas and flow rates only inside the nozzle
               if ((this%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet).or.(this%cfg%xm(i).ge.-0.0015_WP.and.this%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit).or.this%cfg%xm(i).gt.0.0_WP) then
                  CSA_f(i)=CSA_f(i)+area*this%cfg%VF(i,j,k)
                  CSA_l(i)=CSA_l(i)+area*this%cfg%VF(i,j,k)*        this%vf%VF(i,j,k)
                  CSA_g(i)=CSA_g(i)+area*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))
                  VFR_f(i)=VFR_f(i)+area*this%cfg%VF(i,j,k)                           *this%Ui(i,j,k)
                  VFR_l(i)=VFR_l(i)+area*this%cfg%VF(i,j,k)*        this%vf%VF(i,j,k) *this%Ui(i,j,k)
                  VFR_g(i)=VFR_g(i)+area*this%cfg%VF(i,j,k)*(1.0_WP-this%vf%VF(i,j,k))*this%Ui(i,j,k)
               end if
            end do
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_s,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_f,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_l,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,CSA_g,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_s,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_f,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_l,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFR_g,this%cfg%nx,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
      ! Only root process outputs to a file
      if (this%cfg%amRoot) then
         if (.not.isdir('flowrate')) call makedir('flowrate')
         filename='flowrate_'; write(timestamp,'(es12.5)') this%time%t
         open(newunit=iunit,file='flowrate/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(999999(a12,x))') 'xm','CSA_s','CSA_f','CSA_l','CSA_g','VFR_s','VFR_f','VFR_l','VFR_g'
         do i=this%cfg%imin,this%cfg%imax
            write(iunit,'(999999(es12.5,x))') this%cfg%xm(i),CSA_s(i),CSA_f(i),CSA_l(i),CSA_g(i),VFR_s(i),VFR_f(i),VFR_l(i),VFR_g(i)
         end do
         close(iunit)
      end if
   end subroutine analyze_flowrate
   
   
   !> Transfer droplet to Lagrangian representation
   subroutine transfer_drops(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi
      class(simplex), intent(inout) :: this
      real(WP), dimension(:)    , allocatable :: dvol
      real(WP), dimension(:,:)  , allocatable :: dpos
      real(WP), dimension(:,:)  , allocatable :: dvel
      real(WP), dimension(:,:,:), allocatable :: dmoi
      real(WP), dimension(:)    , allocatable :: drem
      integer :: n,m,ierr,i,j,k,nmax
      real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
      logical :: transfer
      ! Moment of inertia calculation using lapack
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                            !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info
      
      ! Query optimal work array size
      if (.not.allocated(work)) then
         call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
         lwork=int(lwork_query(1)); allocate(work(lwork))
      end if
      
      ! Start by performing a CCL
      call this%ccl%build(make_label,same_label)
      
      ! Allocate droplet stats arrays
      allocate(dvol(1:this%ccl%nstruct        )); dvol=0.0_WP
      allocate(dpos(1:this%ccl%nstruct,1:3    )); dpos=0.0_WP
      allocate(dvel(1:this%ccl%nstruct,1:3    )); dvel=0.0_WP
      allocate(dmoi(1:this%ccl%nstruct,1:3,1:3)); dmoi=0.0_WP
      allocate(drem(1:this%ccl%nstruct        )); drem=0.0_WP
      
      ! First pass to accumulate volume, position, and velocity
      do n=1,this%ccl%nstruct
         ! Loop over cells in structure
         do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m)
            j=this%ccl%struct(n)%map(2,m)
            k=this%ccl%struct(n)%map(3,m)
            ! Get cell position, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL
            y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL
            z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL
            ! Accumulate volume, position, and velocity
            dvol(n  )=dvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            dpos(n,:)=dpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
            dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
            ! ! Check if struct touches auto-transfer layer
            ! if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
            ! &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
            ! &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
            ! &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
            ! &   k.ge.this%vf%cfg%kmax-this%nlayer) drem(n)=1.0_WP
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dpos,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvel,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,drem,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      
      ! Second pass to accumulate moment of inertia
      do n=1,this%ccl%nstruct
         ! Get drop barycenter
         x0=dpos(n,1)/dvol(n)
         y0=dpos(n,2)/dvol(n)
         z0=dpos(n,3)/dvol(n)
         ! Loop over cells in structure
         do m=1,this%ccl%struct(n)%n_
            ! Get cell indices
            i=this%ccl%struct(n)%map(1,m)
            j=this%ccl%struct(n)%map(2,m)
            k=this%ccl%struct(n)%map(3,m)
            ! Get cell position relative to drop barycenter, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL-x0
            y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL-y0
            z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL-z0
            ! Accumulate moment of inertia
            dmoi(n,1,1)=dmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
            dmoi(n,2,2)=dmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
            dmoi(n,3,3)=dmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
            dmoi(n,1,2)=dmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
            dmoi(n,1,3)=dmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
            dmoi(n,2,3)=dmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      
      ! Third pass to generate normalized drop stats
      do n=1,this%ccl%nstruct
         ! Get drop barycenter, accounting for periodicity
         dpos(n,:)=dpos(n,:)/dvol(n)
         if (this%vf%cfg%xper.and.dpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) dpos(n,1)=dpos(n,1)+this%vf%cfg%xL
         if (this%vf%cfg%yper.and.dpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%vf%cfg%yL
         if (this%vf%cfg%zper.and.dpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%vf%cfg%zL
         ! Get drop velocity
         dvel(n,:)=dvel(n,:)/dvol(n)
      end do
      
      ! Find the liquid core
      nmax=maxloc(dvol,dim=1)
      
      ! Zero out monitoring variables
      this%vof_transfered=0.0_WP
      this%vof_deleted=0.0_WP
      this%lp%np_new=0
      this%lp%vp_new=0.0_WP
      
      ! Transfer drops based on our criteria
      do n=1,this%ccl%nstruct
         
         ! Compute diameter
         diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)
         
         ! Decide whether to transfer based on diameter
         if (diam.gt.this%dmax) then
            ! Too big to transfer
            transfer=.false.
         else if (diam.le.this%ddel) then
            ! Too small to track, delete immediately
            transfer=.false.
            ! Zero out VF in the structure
            do m=1,this%ccl%struct(n)%n_
               this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            ! Increment monitoring variables
            this%vof_deleted=this%vof_deleted+dvol(n)
         else if (diam.gt.this%ddel.and.diam.le.this%dmin) then
            ! Small enough to transfer automatically
            transfer=.true.
         else
            ! In between, check eccentricity from moment of inertia tensor
            A=dmoi(n,:,:)
            call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
            d=max(0.0_WP,d)                             !< Get rid of very small negative values (due to machine accuracy)
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            if (ecc.gt.this%emax) then
               ! Too eccentric to transfer yet
               transfer=.false.
            else
               ! Spherical enough to transfer
               transfer=.true.
            end if
         end if
         
         ! Force transfer if drop touches auto-transfer layer
         ! if (drem(n).gt.0.0_WP) transfer=.true.
         
         ! But prevent transfer if that's the core
         if (n.eq.nmax) transfer=.false.

         ! Perform transfer
         if (transfer) then
            
            ! Root creates a new Lagrangian drop
            if (this%vf%cfg%amRoot) then
               ! Increment particle counter
               this%lp%np_=this%lp%np_+1
               ! Make room for new drop
               call this%lp%resize(this%lp%np_)
               ! Add the drop
               this%lp%p(this%lp%np_)%id  =int(1,8)
               this%lp%p(this%lp%np_)%d   =diam
               this%lp%p(this%lp%np_)%pos =dpos(n,:)
               this%lp%p(this%lp%np_)%vel =dvel(n,:)
               this%lp%p(this%lp%np_)%ind =this%lp%cfg%get_ijk_global(dpos(n,:),[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])
               this%lp%p(this%lp%np_)%flag=0
               this%lp%p(this%lp%np_)%dt  =0.0_WP
               this%lp%p(this%lp%np_)%Acol=0.0_WP
               this%lp%p(this%lp%np_)%Tcol=0.0_WP
            end if
            
            ! Zero out VF in the structure
            do m=1,this%ccl%struct(n)%n_
               this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
            end do
            
            ! Increment monitoring variables
            this%vof_transfered=this%vof_transfered+dvol(n)
            this%lp%np_new=this%lp%np_new+1
            this%lp%vp_new=this%lp%vp_new+dvol(n)

         end if
         
      end do
      
      ! Synchronize VF fields
      call this%vf%sync_interface()
      call this%vf%clean_irl_and_band()
      
      ! Synchronize particles
      call this%lp%sync()
      
      ! Deallocate all but work array
      deallocate(dvol,dpos,dvel,dmoi,drem)
      
   contains
      
      !> Function that identifies cells that need a label
      logical function make_label(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
            make_label=.true.
         else
            make_label=.false.
         end if
      end function make_label
      

      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
      
   end subroutine transfer_drops

   !> Transfer ligaments to Lagrangian representation
   subroutine transfer_ligaments(this)
      use mathtools, only: pi,twoPi
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use messager, only: die
      use irl_fortran_interface
      implicit none
      class(simplex), intent(inout) :: this
      integer  :: n,nn,i,j,k,ii,jj,kk,ierr,np,ip,m,iunit,rank
      ! Stats of the ccl objects
      real(WP), dimension(:)    , allocatable :: x,y,z,u,v,w,vol,maxlength
      real(WP), dimension(:,:,:), allocatable :: axes
      real(WP), dimension(:,:)  , allocatable :: lengths
      real(WP), dimension(:)    , allocatable :: min_thickness,f_ligament
      real(WP), dimension(:)    , allocatable :: drem
      real(WP) :: myint,integral
      logical :: transfer
      
      ! Varaibles determing transfer
      real(WP) :: Vt,Vl,Vd,minor_radius,diam,Vrim,Lrim
      integer  :: nmain,nsat,np_old,np_start
      ! Prescribe inviscid breakup parameters
      real(WP), parameter :: min_diam=1.0e-6_WP
      real(WP), parameter :: dimless_wavenumber=0.697_WP
      real(WP), parameter :: size_ratio=0.707_WP !0.015_WP

      ! Start by performing a CCL
      call this%ccl_ligament%build(make_label_ligament,same_label_ligament)

      ! stats of the ligaments that will be used for modeling breakup
      allocate(x(1:this%ccl_ligament%nstruct),y(1:this%ccl_ligament%nstruct),z(1:this%ccl_ligament%nstruct));x=0.0_WP;y=0.0_WP;z=0.0_WP
      allocate(u(1:this%ccl_ligament%nstruct),v(1:this%ccl_ligament%nstruct),w(1:this%ccl_ligament%nstruct));u=0.0_WP;v=0.0_WP;w=0.0_WP
      allocate(maxlength(1:this%ccl_ligament%nstruct),lengths(1:this%ccl_ligament%nstruct,1:3));maxlength=0.0_WP;lengths=0.0_WP
      allocate(vol(1:this%ccl_ligament%nstruct),axes(1:this%ccl_ligament%nstruct,1:3,1:3));vol=0.0_WP;axes=0.0_WP
      allocate(min_thickness(1:this%ccl_ligament%nstruct),f_ligament(1:this%ccl_ligament%nstruct)); min_thickness = this%vf%cfg%min_meshsize; f_ligament=0.0_WP
      allocate(drem(1:this%ccl_ligament%nstruct)); drem=0.0_WP
      
      myint =0.0_WP; integral =0.0_WP
   
      ! if (this%vf%cfg%amRoot) print *, "breakup1"
      call this%get_cclstats(this%ccl_ligament,x,y,z,u,v,w,vol,lengths,maxlength,axes,f_ligament)
      ! if (this%vf%cfg%amRoot) print *, "breakup2"
      call this%get_structminthickness(this%ccl_ligament,min_thickness,1)

      ! Tag ligaments touching buffer layer
      do n=1,this%ccl_ligament%nstruct
         ! Loop over cells in structure
         do m=1,this%ccl_ligament%struct(n)%n_
            ! Get cell indices
            i=this%ccl_ligament%struct(n)%map(1,m)
            j=this%ccl_ligament%struct(n)%map(2,m)
            k=this%ccl_ligament%struct(n)%map(3,m)
            ! Check if struct touches auto-transfer layer
            if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
            &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
            &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
            &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
            &   k.ge.this%vf%cfg%kmax-this%nlayer) drem(n)=1.0_WP
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,drem,1*this%ccl_ligament%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

      np_start=this%lp%np_
      do n=1,this%ccl_ligament%nstruct
         
         ! Tag if struct is eligible to be broken up
         transfer=.true.
         if (min_thickness(n).gt.min_ligamentthickness*this%vf%cfg%min_meshsize) then 
            transfer=.false.
         else if (vol(n).lt.1.0_WP*this%vf%cfg%min_meshsize**3) then
            transfer=.false.
         else if (drem(n).gt.0.0_WP) then ! If its in the butter automatically break it
            transfer=.true.
         else if (f_ligament(n).lt.0.9_WP) then 
            transfer=.false.
         end if 

         ! Cycle if not being transfered
         if (transfer.eqv..false.) cycle 

         ! Drop size method from Kim & Moin (2011)
         Lrim=maxlength(n) ! Assume a cylinder ligament
         Vrim=vol(n)
         minor_radius=sqrt(Vrim/pi/Lrim)                 
         nmain=floor(dimless_wavenumber*Lrim/twoPi/minor_radius)
         ! Skip if not a droplet is formed
         if (nmain.lt.1) cycle
         nsat=nmain+1
         diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)
         ! Restriction on the smallest droplet diameter via breakup
         diam=max(diam,min_diam)

         if (nmain.gt.1) then
            Vd=pi/6.0_WP*(diam**3+(size_ratio*diam)**3)
            Vt=0.0_WP; Vl=0.0_WP     
            np_old=this%lp%np_  ! Remember old number of particles
            do nn=1,this%ccl_ligament%struct(n)%n_
               i=this%ccl_ligament%struct(n)%map(1,nn); j=this%ccl_ligament%struct(n)%map(2,nn); k=this%ccl_ligament%struct(n)%map(3,nn)
               ! Increment liquid volume to remove
               Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               ! Create drops from available liquid volume
               do while (Vl-Vd.gt.0.0_WP)            
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(2,8)                                                                                          !< Give id 
                  this%lp%p(np)%dt  =0.0_WP                                                                                            !< Let the drop find it own integration time
                  this%lp%p(np)%Acol=0.0_WP                                                                                            !< Give zero collision force
                  this%lp%p(np)%Tcol=0.0_WP                                                                                            !< Give zero collision force
                  this%lp%p(np)%d   =diam                                                                                              !< Assign diameter to account for full volume
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                                                                            !< Place the drop at the liquid barycenter
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)!< Interpolate local cell velocity as drop velocity
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])!< Place the drop in the proper cell for the this%lp%cfg
                  this%lp%p(np)%flag=0                                                                                                 !< Activate it
                  ! Increment particle counter
                  this%lp%np_=np
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(2,8)                                                                                   
                  this%lp%p(np)%dt  =0.0_WP                                                                                     
                  this%lp%p(np)%Acol=0.0_WP                                                                                     
                  this%lp%p(np)%Tcol=0.0_WP                                                                                     
                  this%lp%p(np)%d   =size_ratio*diam                                                                                       
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+2.0_WP*diam*axes(n,:,1)
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    
                  this%lp%p(np)%flag=0                                                                                          
                  ! Increment particle counter
                  this%lp%np_=np
                  ! Update tracked volumes
                  Vl=Vl-Vd
                  Vt=Vt+Vd
               end do
               ! Remove liquid in that cell
               myint = myint + this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               this%vf%VF(i,j,k)=0.0_WP
            end do
            ! Based on how many particles were created, decide what to do with left-over volume
            if (Vl.gt.0.0_WP) then
               if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                  ! Add one last drop for remaining liquid volume
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(2,8)                                 
                  this%lp%p(np)%dt  =0.0_WP                                    
                  this%lp%p(np)%Acol=0.0_WP                                    
                  this%lp%p(np)%Tcol=0.0_WP                                    
                  this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)           
                  this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                    
                  this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) 
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) 
                  this%lp%p(np)%flag=0                                  
                  ! Increment particle counter
                  this%lp%np_=np
               else ! Some particles were created, make them all larger
                  do ip=np_old+1,this%lp%np_
                     this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                  end do
               end if
            end if
   
         else ! nmain=1
            if (this%vf%cfg%amRoot) then
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(3,8)                                                                               
               this%lp%p(np)%dt  =0.0_WP                                                                                  
               this%lp%p(np)%Acol=0.0_WP                                                                                  
               this%lp%p(np)%Tcol=0.0_WP                                                                                  
               this%lp%p(np)%d   =diam                                                                                    
               this%lp%p(np)%pos =[x(n),y(n),z(n)] 
               this%lp%p(np)%vel =[u(n),v(n),w(n)] 
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])     
               this%lp%p(np)%flag=0                                                                                        
               ! Increment particle counter
               this%lp%np_=np
               do nn=1,2
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(3,8)                                                                               
                  this%lp%p(np)%dt  =0.0_WP                                                                                  
                  this%lp%p(np)%Acol=0.0_WP                                                                                  
                  this%lp%p(np)%Tcol=0.0_WP                                                                                  
                  this%lp%p(np)%d   =size_ratio*diam                                                                         
                  this%lp%p(np)%pos =[x(n),y(n),z(n)] &
                                    +sign(1.0_WP,real(nn,WP)-1.5_WP)*2.0_WP*diam*axes(n,:,1)
                  this%lp%p(np)%vel =[u(n),v(n),w(n)] 
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) 
                  this%lp%p(np)%flag=0                                                                                          
                  ! Increment particle counter
                  this%lp%np_=np
               end do
            end if
   
            do nn=1,this%ccl_ligament%struct(n)%n_
               i=this%ccl_ligament%struct(n)%map(1,nn); j=this%ccl_ligament%struct(n)%map(2,nn); k=this%ccl_ligament%struct(n)%map(3,nn)
               myint = myint + this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               this%vf%VF(i,j,k)=0.0_WP
            end do    
         end if

      end do
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
   
      call MPI_ALLREDUCE(myint,integral,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)    
      ! this%vol_convert = this%vol_convert + integral 

      call this%lp%sync()
      deallocate(x,y,z,u,v,w,vol,lengths,maxlength,axes,drem)
      
      
      contains


      !> Function that identifies cells that need a label to min thickness region
      logical function make_label_ligament(i,j,k)
         use vfs_class, only: VFlo
         implicit none
         integer, intent(in) :: i,j,k
         ! ZZ original value = 1.5
         ! Try 3.584
         ! if ((this%vf%VF(i,j,k).gt.VFlo).and.this%unfiltered_thickness(i,j,k).lt.1.5_WP*this%vf%cfg%min_meshsize) then
         if ((this%vf%VF(i,j,k).gt.VFlo).and.this%unfiltered_thickness(i,j,k).lt.3.584_WP*this%vf%cfg%min_meshsize) then
            make_label_ligament=.true.
         else
            make_label_ligament=.false.
         end if
      end function make_label_ligament
      
      
      !> Function that identifies if cell pairs have same label
      logical function same_label_ligament(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label_ligament=.true.
      end function same_label_ligament

   end subroutine transfer_ligaments

   !> Transfer structs in the buffer layer
   subroutine transfer_buffer(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MAX,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: pi
      class(simplex), intent(inout) :: this
      real(WP), dimension(:)    , allocatable :: dvol
      real(WP), dimension(:,:)  , allocatable :: dpos
      real(WP), dimension(:,:)  , allocatable :: dvel
      real(WP), dimension(:,:,:), allocatable :: dmoi
      real(WP), dimension(:)    , allocatable :: drem
      integer :: n,m,ierr,i,j,k,nmax
      real(WP) :: x,y,z,x0,y0,z0,diam,ecc,lmax,lmid,lmin
      logical :: transfer, write_stats
      ! Moment of inertia calculation using lapack
      real(WP), dimension(:), allocatable, save :: work !< Saved!
      integer, save :: lwork                            !< Saved!
      real(WP), dimension(1) :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer :: info
      
      ! Query optimal work array size
      if (.not.allocated(work)) then
         call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
         lwork=int(lwork_query(1)); allocate(work(lwork))
      end if
      
      ! Start by performing a CCL
      call this%ccl_buffer%build(make_label_buffer,same_label)
      
      ! Allocate droplet stats arrays
      allocate(dvol(1:this%ccl_buffer%nstruct        )); dvol=0.0_WP
      allocate(dpos(1:this%ccl_buffer%nstruct,1:3    )); dpos=0.0_WP
      allocate(dvel(1:this%ccl_buffer%nstruct,1:3    )); dvel=0.0_WP
      allocate(dmoi(1:this%ccl_buffer%nstruct,1:3,1:3)); dmoi=0.0_WP
      allocate(drem(1:this%ccl_buffer%nstruct        )); drem=0.0_WP
      
      ! First pass to accumulate volume, position, and velocity
      do n=1,this%ccl_buffer%nstruct
         ! Loop over cells in structure
         do m=1,this%ccl_buffer%struct(n)%n_
            ! Get cell indices
            i=this%ccl_buffer%struct(n)%map(1,m)
            j=this%ccl_buffer%struct(n)%map(2,m)
            k=this%ccl_buffer%struct(n)%map(3,m)
            ! Get cell position, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl_buffer%struct(n)%per(1)*this%vf%cfg%xL
            y=this%vf%cfg%ym(j)-this%ccl_buffer%struct(n)%per(2)*this%vf%cfg%yL
            z=this%vf%cfg%zm(k)-this%ccl_buffer%struct(n)%per(3)*this%vf%cfg%zL
            ! Accumulate volume, position, and velocity
            dvol(n  )=dvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            dpos(n,:)=dpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
            dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
            ! Check if struct touches auto-transfer layer
            if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
            &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
            &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
            &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
            &   k.ge.this%vf%cfg%kmax-this%nlayer) drem(n)=1.0_WP
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl_buffer%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dpos,3*this%ccl_buffer%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvel,3*this%ccl_buffer%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,drem,1*this%ccl_buffer%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

      
      ! Second pass to accumulate moment of inertia
      do n=1,this%ccl_buffer%nstruct
         ! Get drop barycenter
         x0=dpos(n,1)/dvol(n)
         y0=dpos(n,2)/dvol(n)
         z0=dpos(n,3)/dvol(n)
         ! Loop over cells in structure
         do m=1,this%ccl_buffer%struct(n)%n_
            ! Get cell indices
            i=this%ccl_buffer%struct(n)%map(1,m)
            j=this%ccl_buffer%struct(n)%map(2,m)
            k=this%ccl_buffer%struct(n)%map(3,m)
            ! Get cell position relative to drop barycenter, accounting for periodicity
            x=this%vf%cfg%xm(i)-this%ccl_buffer%struct(n)%per(1)*this%vf%cfg%xL-x0
            y=this%vf%cfg%ym(j)-this%ccl_buffer%struct(n)%per(2)*this%vf%cfg%yL-y0
            z=this%vf%cfg%zm(k)-this%ccl_buffer%struct(n)%per(3)*this%vf%cfg%zL-z0
            ! Accumulate moment of inertia
            dmoi(n,1,1)=dmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
            dmoi(n,2,2)=dmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
            dmoi(n,3,3)=dmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
            dmoi(n,1,2)=dmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
            dmoi(n,1,3)=dmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
            dmoi(n,2,3)=dmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*this%ccl_buffer%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      
      ! Third pass to generate normalized drop stats
      do n=1,this%ccl_buffer%nstruct
         ! Get drop barycenter, accounting for periodicity
         dpos(n,:)=dpos(n,:)/dvol(n)
         if (this%vf%cfg%xper.and.dpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) dpos(n,1)=dpos(n,1)+this%vf%cfg%xL
         if (this%vf%cfg%yper.and.dpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%vf%cfg%yL
         if (this%vf%cfg%zper.and.dpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%vf%cfg%zL
         ! Get drop velocity
         dvel(n,:)=dvel(n,:)/dvol(n)
      end do
      
      ! Find the liquid core
      nmax=maxloc(dvol,dim=1)
      
      ! Zero out monitoring variables
      this%vof_transfered=0.0_WP
      this%vof_deleted=0.0_WP
      this%lp%np_new=0
      this%lp%vp_new=0.0_WP
      
      ! Transfer drops based on our criteria
      do n=1,this%ccl_buffer%nstruct

         ! Cycle if struct in NOT in buffer layer
         if (drem(n).lt.1.0_WP) cycle
         
         ! Compute diameter
         diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)
         
         ! Decide whether to transfer based on diameter
         if (diam.gt.this%dmax) then
            ! Too big to transfer
            transfer=.false.
         else if (diam.le.this%ddel) then
            ! Too small to track, delete immediately
            transfer=.false.
            ! Zero out VF in the structure
            do m=1,this%ccl_buffer%struct(n)%n_
               this%vf%VF(this%ccl_buffer%struct(n)%map(1,m),this%ccl_buffer%struct(n)%map(2,m),this%ccl_buffer%struct(n)%map(3,m))=0.0_WP
            end do
            ! Increment monitoring variables
            this%vof_deleted=this%vof_deleted+dvol(n)
         else if (diam.gt.this%ddel.and.diam.le.this%dmin) then
            ! Small enough to transfer automatically
            transfer=.true.
         else
            ! In between, check eccentricity from moment of inertia tensor
            A=dmoi(n,:,:)
            call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
            d=max(0.0_WP,d)                             !< Get rid of very small negative values (due to machine accuracy)
            ! Get characteristic lengths of drop
            lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
            lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
            lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
            if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
            ecc=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
            if (ecc.gt.this%emax) then
               ! Too eccentric to transfer yet
               transfer=.false.
               write_stats=.true.
            else
               ! Spherical enough to transfer
               transfer=.true.
            end if
         end if
         
         ! But prevent transfer if that's the core
         if (n.eq.nmax) transfer=.false.

         ! Perform transfer
         if (transfer) then
            
            ! Root creates a new Lagrangian drop
            if (this%vf%cfg%amRoot) then
               ! Increment particle counter
               this%lp%np_=this%lp%np_+1
               ! Make room for new drop
               call this%lp%resize(this%lp%np_)
               ! Add the drop
               this%lp%p(this%lp%np_)%id  =int(4,8)
               this%lp%p(this%lp%np_)%d   =diam
               this%lp%p(this%lp%np_)%pos =dpos(n,:)
               this%lp%p(this%lp%np_)%vel =dvel(n,:)
               this%lp%p(this%lp%np_)%ind =this%lp%cfg%get_ijk_global(dpos(n,:),[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])
               this%lp%p(this%lp%np_)%flag=0
               this%lp%p(this%lp%np_)%dt  =0.0_WP
               this%lp%p(this%lp%np_)%Acol=0.0_WP
               this%lp%p(this%lp%np_)%Tcol=0.0_WP
            end if
            
            ! Zero out VF in the structure
            do m=1,this%ccl_buffer%struct(n)%n_
               this%vf%VF(this%ccl_buffer%struct(n)%map(1,m),this%ccl_buffer%struct(n)%map(2,m),this%ccl_buffer%struct(n)%map(3,m))=0.0_WP
            end do
            
            ! Increment monitoring variables
            this%vof_transfered=this%vof_transfered+dvol(n)
            this%lp%np_new=this%lp%np_new+1
            this%lp%vp_new=this%lp%vp_new+dvol(n)

         end if

         ! Write stats if struct cannot be transfered and zero out VF in struct
         if (write_stats) then

            ! Zero out VF in the structure
            do m=1,this%ccl_buffer%struct(n)%n_
               this%vf%VF(this%ccl_buffer%struct(n)%map(1,m),this%ccl_buffer%struct(n)%map(2,m),this%ccl_buffer%struct(n)%map(3,m))=0.0_WP
            end do
            
            ! Increment monitoring variables
            this%vof_transfered=this%vof_transfered+dvol(n)
            
         end if
         
      end do
      
      ! Synchronize VF fields
      call this%vf%sync_interface()
      call this%vf%clean_irl_and_band()
      
      ! Synchronize particles
      call this%lp%sync()
      
      ! Deallocate all but work array
      deallocate(dvol,dpos,dvel,dmoi,drem)
      
   contains
      
      !> Function that identifies cells that need a label
      logical function make_label_buffer(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.0_WP) then
            make_label_buffer=.true.
         else
            make_label_buffer=.false.
         end if
      end function make_label_buffer
      

      !> Function that identifies if cell pairs have same label
      logical function same_label(i1,j1,k1,i2,j2,k2)
         implicit none
         integer, intent(in) :: i1,j1,k1,i2,j2,k2
         same_label=.true.
      end function same_label
      
   end subroutine transfer_buffer

   !> Measure local thickness of multiphasic structure
   subroutine get_thickness_unfiltered(this)
      use vfs_class, only: VFlo,VFhi
      implicit none
      class(simplex), intent(inout) :: this
      integer :: i,j,k,ii,jj,kk,nneigh
      real(WP) :: lvol,gvol,area
      ! Reset thickness
      tmpthickness=1.0_WP;nneigh=3!nneigh=1
      ! First compute thickness based on current surface and volume moments (SD and VF)
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               lvol=0.0_WP; area=0.0_WP
               do kk = k-nneigh,k+nneigh
                  do jj = j-nneigh,j+nneigh
                     do ii = i-nneigh,i+nneigh
                        lvol = lvol + this%vf%VF(ii,jj,kk)
                        area = area + this%vf%SD(ii,jj,kk)
                     end do
                  end do
               end do
               if (this%vf%VF(i,j,k).lt.VFlo) then
                  tmpthickness(i,j,k) = 0.0_WP
               else if (area .gt. 0.0_WP) then    
                  tmpthickness(i,j,k) = 2.0_WP*lvol/(area+tiny(1.0_WP))
               else
                  tmpthickness(i,j,k) = 3.0_WP*this%cfg%min_meshsize
               end if
            end do
         end do
      end do
      call this%vf%cfg%sync(tmpthickness)
      this%unfiltered_thickness = tmpthickness
   end subroutine get_thickness_unfiltered
   
   
   !> Initialization of simplex simulation
   subroutine init(this)
      implicit none
      class(simplex), intent(inout) :: this
      
      
      ! Setup an input file
      read_input: block
         use parallel, only: amRoot
         this%input=inputfile(amRoot=amRoot,filename='simplex.input')
      end block read_input
      
      
      ! Initialize ibconfig object
      create_config: block
         use parallel,    only: group
         use sgrid_class, only: cartesian,sgrid
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz,ns_yz,ns_x
         real(WP) :: Lx,Ly,Lz,xshift,sratio_yz,sratio_x
         real(WP), dimension(:), allocatable :: x_uni,y_uni,z_uni,x,y,z
         integer, dimension(3) :: partition
         ! Read in grid definition
         call this%input%read('Lx',Lx); call this%input%read('nx',nx); allocate(x_uni(nx+1)); call this%input%read('X shift',xshift)
         call this%input%read('Ly',Ly); call this%input%read('ny',ny); allocate(y_uni(ny+1))
         call this%input%read('Lz',Lz); call this%input%read('nz',nz); allocate(z_uni(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x_uni(i)=real(i-1,WP)/real(nx,WP)*Lx-xshift
         end do
         do j=1,ny+1
            y_uni(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z_uni(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! Add stretching
         call this%input%read('Stretched cells in yz',ns_yz,default=0)
         if (ns_yz.gt.0) call this%input%read('Stretch ratio in yz',sratio_yz)
         call this%input%read('Stretched cells in x' ,ns_x ,default=0)
         if (ns_x .gt.0) call this%input%read('Stretch ratio in x' ,sratio_x )
         allocate(x(nx+1+1*ns_x )); x(      1:      1+nx)=x_uni
         allocate(y(ny+1+2*ns_yz)); y(ns_yz+1:ns_yz+1+ny)=y_uni
         allocate(z(nz+1+2*ns_yz)); z(ns_yz+1:ns_yz+1+nz)=z_uni
         do i=nx+2,nx+1+ns_x
            x(i)=x(i-1)+sratio_x *(x(i-1)-x(i-2))
         end do
         do j=ns_yz,1,-1
            y(j)=y(j+1)+sratio_yz*(y(j+1)-y(j+2))
         end do
         do j=ns_yz+2+ny,ny+1+2*ns_yz
            y(j)=y(j-1)+sratio_yz*(y(j-1)-y(j-2))
         end do
         do k=ns_yz,1,-1
            z(k)=z(k+1)+sratio_yz*(z(k+1)-z(k+2))
         end do
         do k=ns_yz+2+nz,nz+1+2*ns_yz
            z(k)=z(k-1)+sratio_yz*(z(k-1)-z(k-2))
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='simplex')
         ! Read in partition
         call this%input%read('Partition',partition)
         ! Create ibconfig
         this%cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_config
      
      
      ! Now initialize simplex nozzle geometry
      create_simplex: block
         use ibconfig_class, only: sharp
         integer :: i,j,k
         ! Create polygon
         call this%poly%initialize(nvert=10,name='simplex')
         this%poly%vert(:, 1)=[-0.01000_WP,0.00000_WP]
         this%poly%vert(:, 2)=[-0.00442_WP,0.00000_WP]
         this%poly%vert(:, 3)=[-0.00442_WP,0.00160_WP]
         this%poly%vert(:, 4)=[-0.00385_WP,0.00160_WP]
         this%poly%vert(:, 5)=[-0.00175_WP,0.00039_WP]
         this%poly%vert(:, 6)=[-0.00114_WP,0.00039_WP]
         this%poly%vert(:, 7)=[ 0.00000_WP,0.00143_WP]
         this%poly%vert(:, 8)=[ 0.00000_WP,0.00177_WP]
         this%poly%vert(:, 9)=[-0.00122_WP,0.00279_WP]
         this%poly%vert(:,10)=[-0.01000_WP,0.00279_WP]
         ! Initialize IB distance field
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  ! Calculate distance from object obtained by revolution of polygon for the plenum
                  this%cfg%Gib(i,j,k)=-this%poly%get_distance([this%cfg%xm(i),sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2)])
               end do
            end do
         end do
         ! Get normal vector
         call this%cfg%calculate_normal()
         ! Get VF field
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         ! Carve out inlet pipes
         create_inlet_pipes: block
            use mms_geom, only: cube_refine_vol
            integer :: si,sj,sk,n
            real(WP), dimension(3,8) :: cube_vertex
            real(WP), dimension(3) :: v_cent,a_cent
            real(WP) :: vol,area,contact
            integer, parameter :: amr_ref_lvl=4
            do k=this%cfg%kmino_,this%cfg%kmaxo_
               do j=this%cfg%jmino_,this%cfg%jmaxo_
                  do i=this%cfg%imino_,this%cfg%imaxo_
                     ! Only work to the left of the plenum
                     if (this%cfg%xm(i)-this%p1(1).gt.this%cfg%min_meshsize) cycle
                     ! Set cube vertices
                     n=0
                     do sk=0,1
                        do sj=0,1
                           do si=0,1
                              n=n+1; cube_vertex(:,n)=[this%cfg%x(i+si),this%cfg%y(j+sj),this%cfg%z(k+sk)]
                           end do
                        end do
                     end do
                     ! Call adaptive refinement code to get volume fraction recursively - inlet pipe 1
                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_inlet_pipe_1,0.0_WP,amr_ref_lvl)
                     this%cfg%VF(i,j,k)=max(this%cfg%VF(i,j,k), vol/this%cfg%vol(i,j,k))
                     this%cfg%SD(i,j,k)=max(this%cfg%SD(i,j,k),area/this%cfg%vol(i,j,k))
                     ! Call adaptive refinement code to get volume fraction recursively - inlet pipe 2
                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_inlet_pipe_2,0.0_WP,amr_ref_lvl)
                     this%cfg%VF(i,j,k)=max(this%cfg%VF(i,j,k), vol/this%cfg%vol(i,j,k))
                     this%cfg%SD(i,j,k)=max(this%cfg%SD(i,j,k),area/this%cfg%vol(i,j,k))
                  end do
               end do
            end do
         end block create_inlet_pipes
         ! Apply Neumann on VF and apply stair-stepping at entrance
         if (this%cfg%iproc.eq.1) then
            ! Stair-step entrance
            this%cfg%VF(this%cfg%imin,:,:)=max(real(nint(this%cfg%VF(this%cfg%imin,:,:)),WP),epsilon(1.0_WP))
            ! Copy into overlap layer
            do i=this%cfg%imino,this%cfg%imin-1
               this%cfg%VF(i,:,:)=this%cfg%VF(this%cfg%imin,:,:)
            end do
         end if
         ! Recompute domain volume
         call this%cfg%calc_fluid_vol()
      end block create_simplex
      
      
      ! Initialize flow rate
      set_flowrate: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         use string,   only: str_long
         use messager, only: log
         integer :: j,k,ierr
         character(len=str_long) :: message
         ! Read mass flow rate
         call this%input%read('Mass flow rate',this%mfr)
         ! Read coflow velocity
         call this%input%read('Coflow velocity',this%Ucoflow)
         ! Integrate inlet pipe surface area
         this%Apipe=0.0_WP
         if (this%cfg%iproc.eq.1) then
            do k=this%cfg%kmin_,this%cfg%kmax_
               do j=this%cfg%jmin_,this%cfg%jmax_
                  if (sqrt(this%cfg%ym(j)**2+this%cfg%zm(k)**2).lt.this%Rinlet) then
                     this%Apipe=this%Apipe+this%cfg%VF(this%cfg%imin,j,k)*this%cfg%dy(j)*this%cfg%dz(k)
                  end if
               end do
            end do
         end if
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%Apipe,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         ! Log this info
         if (this%cfg%amRoot) then
            write(message,'("Inlet pipe area is ",es12.5)') this%Apipe; call log(message)
         end if
      end block set_flowrate
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call this%input%read('Max timestep size',this%time%dtmax)
         call this%input%read('Max cfl number',this%time%cflmax)
         call this%input%read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%SR       (1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Uib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wib (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(tmpthickness(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(tmpfilm_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%unfiltered_thickness(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%unfiltered_thickness=0.0_WP
         allocate(this%film_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%film_type=0.0_WP
         allocate(this%d1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d1=0.0_WP
         allocate(this%d2(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d2=0.0_WP
         allocate(this%d3(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d3=0.0_WP
         allocate(this%l1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%l1=0.0_WP
         allocate(this%l2(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%l2=0.0_WP
         allocate(this%l3(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%l3=0.0_WP
         allocate(this%d3_d1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d3_d1=0.0_WP
         allocate(this%d3_d2(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d3_d2=0.0_WP
         allocate(this%d2_d1(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));this%d2_d1=0.0_WP
      end block allocate_work_arrays
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: remap,plicnet,r2pnet
         integer :: i,j,k
         real(WP) :: rad
         ! Create a VOF solver with plicnet
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         this%vf%flotsam_thld=0.0_WP
         this%vf%maxcurv_times_mesh=1.0_WP
         ! Initialize to flat interface at exit
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                  ! Ensure the nozzle is filled with liquid up to the throat with wet walls
                  if (this%vf%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else if (this%vf%cfg%xm(i).ge.-0.0015_WP.and.this%vf%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else
                     this%vf%VF(i,j,k)=0.0_WP
                  end if
                  ! Initialize phasic barycenters
                  this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
         this%vof_removed=0.0_WP
      end block create_iterator
      
      
      ! Create a two-phase flow solver with bconds
      create_flow_solver: block
         use tpns_class,      only: clipped_neumann,dirichlet,slip
         use hypre_str_class, only: pcg_pfmg2
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-Phase NS')
         ! Set the flow properties
         call this%input%read('Liquid dynamic viscosity',this%fs%visc_l)
         call this%input%read('Gas dynamic viscosity'   ,this%fs%visc_g)
         call this%input%read('Liquid density',this%fs%rho_l)
         call this%input%read('Gas density'   ,this%fs%rho_g)
         call this%input%read('Surface tension coefficient',this%fs%sigma)
         ! Set acceleration of gravity
         call this%input%read('Gravity',this%fs%gravity)
         ! Inlets on the left
         call this%fs%add_bcond(name='inlets',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=pipe_inlets)
         call this%fs%add_bcond(name='coflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=coflow_inlet)
         ! Outflow on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.false.,locator=right_boundary)
         ! Slip on the sides
         call this%fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call this%fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call this%fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         call this%fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call this%input%read('Pressure iteration',this%ps%maxit)
         call this%input%read('Pressure tolerance',this%ps%rcvg)
         ! Configure velocity solver
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         ! Zero velocity except if restarting
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply Dirichlet condition at pipe inlets
         call this%fs%get_bcond('inlets',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=this%cfg%VF(i,j,k)*this%mfr/(this%fs%rho_l*this%Apipe)
         end do
         ! Apply Dirichlet condition at coflow inlet
         call this%fs%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=this%Ucoflow
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block initialize_velocity
      
      
      ! Create an LES model
      create_sgs: block
         this%sgs=sgsmodel(cfg=this%fs%cfg,umask=this%fs%umask,vmask=this%fs%vmask,wmask=this%fs%wmask)
      end block create_sgs
      
      
      ! Prepare Lagrangian drop model
      prepare_transfer: block
         ! Is transfer used?
         call this%input%read('Transfer drops',this%use_drop_transfer,default=.true.)
         ! Only initialize transfer model if used
         if (this%use_drop_transfer) then
            ! Create CCL
            call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
            ! Create lpt solver
            this%lp=lpt(cfg=this%cfg,name='spray')
            this%lp%rho=this%fs%rho_l
            this%lp%gravity=this%fs%gravity
            this%lp%filter_width=3.5_WP*this%cfg%min_meshsize
            call this%lp%resize(0)
            ! Set parameters for transfer
            this%ddel=0.2_WP*this%cfg%min_meshsize
            this%dmin=1.5_WP*this%cfg%min_meshsize
            this%dmax=1.0e-3_WP
            this%emax=0.8_WP
            ! Zero out transfered volume
            this%vof_transfered=0.0_WP
         end if
      end block prepare_transfer

      ! Prepare ligament breakup model
      prepare_ligament_breakup: block
         call this%input%read('Breakup ligaments',this%use_ligament_breakup,default=.true.)
         call this%ccl_ligament%initialize(pg=this%cfg%pgrid,name='ccl_ligament')
         call this%ccl_buffer%initialize(pg=this%cfg%pgrid,name='ccl_buffer')
      end block prepare_ligament_breakup
      
      
      ! Handle restart/saves here
      handle_restart: block
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         use tpns_class,            only: bcond
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         real(WP) :: rad
         integer :: i,j,k,n
         type(bcond), pointer :: mybc
         logical :: partfile_exists
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call this%input%read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call this%input%read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call this%input%read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! We are restarting, read the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(timestamp))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Ensure that boundaries are correct
            if (this%vf%cfg%iproc.eq.1)               this%vf%VF(this%vf%cfg%imino:this%vf%cfg%imin-1,:,:)=0.0_WP
            if (this%vf%cfg%iproc.eq.this%vf%cfg%npx) this%vf%VF(this%vf%cfg%imax+1:this%vf%cfg%imaxo,:,:)=0.0_WP
            if (this%vf%cfg%jproc.eq.1)               this%vf%VF(:,this%vf%cfg%jmino:this%vf%cfg%jmin-1,:)=0.0_WP
            if (this%vf%cfg%jproc.eq.this%vf%cfg%npy) this%vf%VF(:,this%vf%cfg%jmax+1:this%vf%cfg%jmaxo,:)=0.0_WP
            if (this%vf%cfg%kproc.eq.1)               this%vf%VF(:,:,this%vf%cfg%kmino:this%vf%cfg%kmin-1)=0.0_WP
            if (this%vf%cfg%kproc.eq.this%vf%cfg%npz) this%vf%VF(:,:,this%vf%cfg%kmax+1:this%vf%cfg%kmaxo)=0.0_WP
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                     if (i.lt.this%vf%cfg%imin.and.rad.le.this%Rinlet) this%vf%VF(i,j,k)=1.0_WP
                  end do
               end do
            end do
            ! Update the band
            call this%vf%update_band()
            ! Set interface planes at the boundaries
            call this%vf%set_full_bcond()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Reapply inflow boundary conditions in case input has changed
            call this%fs%get_bcond('inlets',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               this%fs%U(i,j,k)=this%cfg%VF(i,j,k)*this%mfr/(this%fs%rho_l*this%Apipe)
            end do
            call this%fs%get_bcond('coflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               this%fs%U(i,j,k)=this%Ucoflow
            end do
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
            ! Finally, handle particle I/O
            if (this%use_drop_transfer) then
               ! Check if particle file exists
               inquire(file='restart/part_'//trim(timestamp),exist=partfile_exists)
               ! If so, read it
               if (partfile_exists) call this%lp%read(filename='restart/part_'//trim(timestamp))
            end if
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         use vfs_class, only: r2p,plicnet,r2pnet,lvira
         integer :: i,j,k,np,nplane
         this%smesh=surfmesh(nvar=6,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='id_ccl_ligament'
         this%smesh%varname(4)='unf_filt'
         this%smesh%varname(5)='film_type'
         this%smesh%varname(6)='id_ccl_drop'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh_nowall(this%smesh)
         ! Calculate thickness even for plic
         if (.not.this%vf%two_planes) then
            allocate(this%vf%thickness(this%vf%cfg%imino_:this%vf%cfg%imaxo_,this%vf%cfg%jmino_:this%vf%cfg%jmaxo_,this%vf%cfg%kmino_:this%vf%cfg%kmaxo_)); this%vf%thickness=0.0_WP
         end if
         call this%vf%get_thickness()
         ! Populate surface variables
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                        this%smesh%var(3,np)=real(this%ccl_ligament%id(i,j,k),WP)
                        this%smesh%var(4,np)=this%unfiltered_thickness(i,j,k)
                        this%smesh%var(5,np)=real(this%film_type(i,j,k),WP)
                        this%smesh%var(6,np)=real(this%ccl%id(i,j,k),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Create partmesh object for particle output
      ! id=1, transfer drops
      ! id=2, transfer_ligaments (nmain>1)
      ! id=3, transfer_ligaments (nmain=1)
      ! id=4, Auto transfer droplets in buffer
      if (this%use_drop_transfer) then
         create_pmesh: block
            integer :: i
            this%pmesh=partmesh(nvar=2,nvec=1,name='lpt')
            this%pmesh%varname(1)='radius'
            this%pmesh%vecname(1)='velocity'
            this%pmesh%varname(2)='id'
            call this%lp%update_partmesh(this%pmesh)
            do i=1,this%lp%np_
               this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
               this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
               this%pmesh%var(2,i)=this%lp%p(i)%id
            end do
         end block create_pmesh
      end if
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='simplex')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call this%input%read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('divergence',this%fs%div)
         call this%ens_out%add_surface('plic',this%smesh)
         if (this%use_drop_transfer) call this%ens_out%add_particle('part',this%pmesh)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_simplex')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vof_deleted,'VOF deleted')
         call this%mfile%add_column(this%vof_transfered,'VOF transfered')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_simplex')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
         ! Create particle monitor
         if (this%use_drop_transfer) then
            call this%lp%get_max()
            this%pfile=monitor(amroot=this%lp%cfg%amRoot,name='particles')
            call this%pfile%add_column(this%time%n,'Timestep number')
            call this%pfile%add_column(this%time%t,'Time')
            call this%pfile%add_column(this%lp%np,'Particle number')
            call this%pfile%add_column(this%lp%vp_tot,'Particle volume')
            call this%pfile%add_column(this%lp%np_new,'Npart new')
            call this%pfile%add_column(this%lp%vp_new,'Vpart new')
            call this%pfile%add_column(this%lp%np_out,'Npart removed')
            call this%pfile%add_column(this%lp%vp_out,'Vpart removed')
            call this%pfile%add_column(this%lp%Umin,'Particle Umin')
            call this%pfile%add_column(this%lp%Umax,'Particle Umax')
            call this%pfile%add_column(this%lp%Vmin,'Particle Vmin')
            call this%pfile%add_column(this%lp%Vmax,'Particle Vmax')
            call this%pfile%add_column(this%lp%Wmin,'Particle Wmin')
            call this%pfile%add_column(this%lp%Wmax,'Particle Wmax')
            call this%pfile%add_column(this%lp%dmin,'Particle dmin')
            call this%pfile%add_column(this%lp%dmax,'Particle dmax')
            call this%pfile%write()
         end if
      end block create_monitor
      
      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep     =timer(comm=this%cfg%comm,name='Timestep')
         this%tvof      =timer(comm=this%cfg%comm,name='VOFsolve')
         this%tvel      =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres     =timer(comm=this%cfg%comm,name='Pressure')
         this%tsgs      =timer(comm=this%cfg%comm,name='SGSmodel')
         this%ttrans    =timer(comm=this%cfg%comm,name='Transfer')
         this%ttrans_lig=timer(comm=this%cfg%comm,name='Transfer_lig')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvof%time  ,trim(this%tvof%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
         call this%timefile%add_column(this%tsgs%time  ,trim(this%tsgs%name))
         call this%timefile%add_column(this%ttrans%time,trim(this%ttrans%name))
         call this%timefile%add_column(this%ttrans_lig%time,trim(this%ttrans_lig%name))
      end block create_timing

      
      ! Create an event for flow rate analysis
      flowrate_analysis_prep: block
         this%flowrate_evt=event(time=this%time,name='Flow rate output')
         call this%input%read('Flow rate output period',this%flowrate_evt%tper,default=huge(1.0_WP))
      end block flowrate_analysis_prep
      
   contains
      
      !> Function that localizes the right domain boundary
      function right_boundary(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function right_boundary
      
      
      !> Function that localizes the pipe inlets
      function pipe_inlets(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).lt.this%Rinlet) isIn=.true.
      end function pipe_inlets
      
      
      !> Function that localizes the coflow inlet
      function coflow_inlet(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin.and.sqrt(pg%ym(j)**2+pg%zm(k)**2).gt.this%Rcoflow) isIn=.true.
      end function coflow_inlet
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer.or.&
         &   j.le.pg%jmin+this%nlayer.or.&
         &   j.ge.pg%jmax-this%nlayer.or.&
         &   k.le.pg%kmin+this%nlayer.or.&
         &   k.ge.pg%kmax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that localizes the top (y+) of the domain
      function yp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmax+1) isIn=.true.
      end function yp_locator
      
      
      !> Function that localizes the bottom (y-) of the domain
      function ym_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.eq.pg%jmin) isIn=.true.
      end function ym_locator
      
      
      !> Function that localizes the top (z+) of the domain
      function zp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmax+1) isIn=.true.
      end function zp_locator
      
      
      !> Function that localizes the bottom (z-) of the domain
      function zm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (k.eq.pg%kmin) isIn=.true.
      end function zm_locator
      
      
      !> Function that defines a level set function for inlet pipe 1
      function levelset_inlet_pipe_1(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP), dimension(3) :: v,p
         real(WP) :: G
         v=xyz-this%p1
         p=v-this%n1*dot_product(v,this%n1)
         G=this%Rpipe-sqrt(dot_product(p,p))
      end function levelset_inlet_pipe_1
      
      
      !> Function that defines a level set function for inlet pipe 2
      function levelset_inlet_pipe_2(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP), dimension(3) :: v,p
         real(WP) :: G
         v=xyz-this%p2
         p=v-this%n2*dot_product(v,this%n2)
         G=this%Rpipe-sqrt(dot_product(p,p))
      end function levelset_inlet_pipe_2
      
   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(simplex), intent(inout) :: this
      
      ! Reset all timers and start timestep timer
      call this%tstep%reset()
      call this%tvof%reset()
      call this%tsgs%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%ttrans%reset()
      call this%ttrans_lig%reset()
      call this%tstep%start()
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()
      
      ! Advance lagrangian droplets
      if (this%use_drop_transfer) then
         this%resU=this%fs%rho_g
         this%resV=this%fs%visc_g
         call this%lp%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W,rho=this%resU,visc=this%resV)
      end if
      
      ! Remember old VOF
      this%vf%VFold=this%vf%VF
      
      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old staggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)
      
      ! VOF solver step
      call this%tvof%start() ! Start VOF timer
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%tvof%stop() ! Stop VOF timer
      
      ! Prepare new staggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Turbulence modeling
      call this%tsgs%start() ! Start SGS timer
      sgs_modeling: block
         use sgsmodel_class, only: vreman,dynamic_smag
         use ibconfig_class, only: VFlo,VFhi
         real(WP), parameter :: Cwm=1.0_WP ! Whitmore, Bose, and Moin, as well as Hausmann and van Wachem
         integer :: i,j,k
         ! Get velocity gradient tensor and strain rate tensor
         call this%fs%get_gradu(this%gradU)
         this%SR(6,:,:,:)=(this%gradU(1,1,:,:,:)+this%gradU(2,2,:,:,:)+this%gradU(3,3,:,:,:))/3.0_WP ! div
         this%SR(1,:,:,:)=this%gradU(1,1,:,:,:)-this%SR(6,:,:,:)                                     ! du/dx-div/3
         this%SR(2,:,:,:)=this%gradU(2,2,:,:,:)-this%SR(6,:,:,:)                                     ! dv/dy-div/3
         this%SR(3,:,:,:)=this%gradU(3,3,:,:,:)-this%SR(6,:,:,:)                                     ! dw/dz-div/3
         this%SR(4,:,:,:)=0.5_WP*(this%gradU(1,2,:,:,:)+this%gradU(2,1,:,:,:))                       ! (du/dy+dv/dx)/2
         this%SR(5,:,:,:)=0.5_WP*(this%gradU(2,3,:,:,:)+this%gradU(3,2,:,:,:))                       ! (dv/dz+dw/dy)/2
         this%SR(6,:,:,:)=0.5_WP*(this%gradU(3,1,:,:,:)+this%gradU(1,3,:,:,:))                       ! (dw/dx+du/dz)/2
         ! Get turbulent viscosity
         this%resU=this%vf%VF*this%fs%rho_l+(1.0_WP-this%vf%VF)*this%fs%rho_g
         !call this%sgs%get_visc(type=dynamic_smag,dt=this%time%dtold,rho=this%resU,Ui=this%Ui,Vi=this%Vi,Wi=this%Wi,SR=this%SR)
         call this%sgs%get_visc(type=vreman,dt=this%time%dtold,rho=this%resU,gradu=this%gradU)
         ! Add sgs visc to our two-phase viscosities
         do k=this%fs%cfg%kmino_+1,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_+1,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_+1,this%fs%cfg%imaxo_
            this%fs%visc(i,j,k)   =this%fs%visc(i,j,k)   +this%sgs%visc(i,j,k)
            this%fs%visc_xy(i,j,k)=this%fs%visc_xy(i,j,k)+sum(this%fs%itp_xy(:,:,i,j,k)*this%sgs%visc(i-1:i,j-1:j,k))
            this%fs%visc_yz(i,j,k)=this%fs%visc_yz(i,j,k)+sum(this%fs%itp_yz(:,:,i,j,k)*this%sgs%visc(i,j-1:j,k-1:k))
            this%fs%visc_zx(i,j,k)=this%fs%visc_zx(i,j,k)+sum(this%fs%itp_xz(:,:,i,j,k)*this%sgs%visc(i-1:i,j,k-1:k))
         end do; end do; end do
         ! Compute slip velocity using Cslip*delta*du/dn
         do k=this%fs%cfg%kmino_,this%fs%cfg%kmaxo_; do j=this%fs%cfg%jmino_,this%fs%cfg%jmaxo_; do i=this%fs%cfg%imino_,this%fs%cfg%imaxo_
            this%Uib(i,j,k)=-Cwm*this%fs%cfg%meshsize(i,j,k)*sum(this%gradU(:,1,i,j,k)*this%cfg%Nib(:,i,j,k))
            this%Vib(i,j,k)=-Cwm*this%fs%cfg%meshsize(i,j,k)*sum(this%gradU(:,2,i,j,k)*this%cfg%Nib(:,i,j,k))
            this%Wib(i,j,k)=-Cwm*this%fs%cfg%meshsize(i,j,k)*sum(this%gradU(:,3,i,j,k)*this%cfg%Nib(:,i,j,k))
         end do; end do; end do
      end block sgs_modeling
      call this%tsgs%stop() ! Stop SGS timer
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Start velocity timer
         call this%tvel%start()
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)
         
         ! Add momentum source terms
         call this%fs%addsrc_gravity(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply IB forcing to enforce wall boundary conditions
         ibforcing: block
            integer :: i,j,k
            real(WP) :: sd
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_; do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_; do i=this%fs%cfg%imin_,this%fs%cfg%imax_
               ! Use different model at inlet vs. inside the nozzle
               if (this%cfg%xm(i).le.this%p1(1)+this%cfg%min_meshsize) then
                  ! U cell
                  if (this%fs%umask(i,j,k).eq.0) then
                     sd=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                     this%fs%U(i,j,k)=sd*this%fs%U(i,j,k)
                  end if
                  ! V cell
                  if (this%fs%vmask(i,j,k).eq.0) then
                     sd=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))
                     this%fs%V(i,j,k)=sd*this%fs%V(i,j,k)
                  end if
                  ! W cell
                  if (this%fs%wmask(i,j,k).eq.0) then
                     sd=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))
                     this%fs%W(i,j,k)=sd*this%fs%W(i,j,k)
                  end if
               else
                  ! U cell
                  if (this%fs%umask(i,j,k).eq.0) then
                     sd=min(sum(this%fs%itpr_x(:,i,j,k)*this%cfg%SD(i-1:i,j,k)*this%cfg%meshsize(i-1:i,j,k)),1.0_WP)
                     this%fs%U(i,j,k)=(1.0_WP-sd)*this%fs%U(i,j,k)+sd*sum(this%fs%itpr_x(:,i,j,k)*this%Uib(i-1:i,j,k))
                  end if
                  ! V cell
                  if (this%fs%vmask(i,j,k).eq.0) then
                     sd=min(sum(this%fs%itpr_y(:,i,j,k)*this%cfg%SD(i,j-1:j,k)*this%cfg%meshsize(i,j-1:j,k)),1.0_WP)
                     this%fs%V(i,j,k)=(1.0_WP-sd)*this%fs%V(i,j,k)+sd*sum(this%fs%itpr_y(:,i,j,k)*this%Vib(i,j-1:j,k))
                  end if
                  ! W cell
                  if (this%fs%wmask(i,j,k).eq.0) then
                     sd=min(sum(this%fs%itpr_z(:,i,j,k)*this%cfg%SD(i,j,k-1:k)*this%cfg%meshsize(i,j,k-1:k)),1.0_WP)
                     this%fs%W(i,j,k)=(1.0_WP-sd)*this%fs%W(i,j,k)+sd*sum(this%fs%itpr_z(:,i,j,k)*this%Wib(i,j,k-1:k))
                  end if
               end if
            end do; end do; end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing
         
         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop velocity timer and start pressure timer
         call this%tvel%stop()
         call this%tpres%start()
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         if (this%vf%two_planes) then
            call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         else
            call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         end if
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho_U
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho_V
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho_W
         
         ! Stop pressure timer
         call this%tpres%stop()
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Transfer VOF into droplets
      call this%ttrans%start() ! Start transfer timer
      if (this%use_drop_transfer) call this%transfer_drops()
      call this%ttrans%stop() ! Stop transfer timer

      ! Transfer ligament to droplets
      call this%ttrans_lig%start() ! Start transfer timer
      if (this%use_ligament_breakup) then 
         call this%get_localfilmtype()
         call this%get_thickness_unfiltered()
         call this%transfer_ligaments()
      end if
      call this%ttrans_lig%stop() ! Stop transfer timer
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
      end block remove_vof
      
      ! Reset VOF in IBs
      reset_VOF_IB: block
         use ibconfig_class, only: VFlo
         integer :: i,j,k
         real(WP) :: rad
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Only work in pure wall/IB cells
                  if (this%vf%cfg%VF(i,j,k).gt.VFlo) cycle
                  ! Compute our local radius
                  rad=sqrt(this%vf%cfg%ym(j)**2+this%vf%cfg%zm(k)**2)
                  ! Ensure the nozzle is filled with liquid up to the throat with wet walls
                  if (this%vf%cfg%xm(i).lt.-0.0015_WP.and.rad.le.this%Rinlet) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else if (this%vf%cfg%xm(i).ge.-0.0015_WP.and.this%vf%cfg%xm(i).lt.0.0_WP.and.rad.le.this%Rexit) then
                     this%vf%VF(i,j,k)=1.0_WP
                  else
                     this%vf%VF(i,j,k)=0.0_WP
                  end if
                  ! Initialize phasic barycenters
                  this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
               end do
            end do
         end do
         call this%vf%clean_irl_and_band()
      end block reset_VOF_IB
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surface mesh
         update_smesh: block
            use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
            use vfs_class, only: plicnet,r2p,r2pnet,lvira
            integer :: i,j,k,np,nplane
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh_nowall(this%smesh)
            ! Calculate thickness even for plic
            if (.not.this%vf%two_planes) call this%vf%get_thickness()
            ! Populate surface variables
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=this%vf%thickness(i,j,k)
                           this%smesh%var(3,np)=real(this%ccl_ligament%id(i,j,k),WP)
                           this%smesh%var(4,np)=this%unfiltered_thickness(i,j,k)
                           this%smesh%var(5,np)=real(this%film_type(i,j,k),WP)
                           this%smesh%var(6,np)=real(this%ccl%id(i,j,k),WP)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Update particle mesh object
         if (this%use_drop_transfer) then
            update_pmesh: block
               integer :: i
               call this%lp%update_partmesh(this%pmesh)
               do i=1,this%lp%np_
                  this%pmesh%var(1,i)=0.5_WP*this%lp%p(i)%d
                  this%pmesh%vec(:,1,i)=this%lp%p(i)%vel
                  this%pmesh%var(2,i)=this%lp%p(i)%id
               end do
            end block update_pmesh 
         end if
         ! Write ensight files
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Output flow rate
      if (this%flowrate_evt%occurs()) call this%analyze_flowrate()

      ! Stop timestep timer
      call this%tstep%stop()
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      if (this%use_drop_transfer) then
         call this%lp%get_max()
         call this%pfile%write()
      end if
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Finally, handle particle I/O
            if (this%use_drop_transfer) call this%lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
         end block save_restart
      end if
      
   end subroutine step
   
   
   !> Finalize simplex simulation
   subroutine final(this)
      implicit none
      class(simplex), intent(inout) :: this
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      deallocate(this%gradU,this%Uib,this%Vib,this%Wib,this%SR)
   end subroutine final

   !> Detect edge-like regions of the interface
   subroutine get_localfilmtype(this)
      implicit none
      class(simplex), intent(inout) :: this
      integer , parameter        :: order = 3
      integer :: ii,jj,kk,i,j,k,lwork,info
      real(WP) :: lvol,ratio
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      real(WP), dimension(3) :: lx_vol, Ltmp,d
      real(WP), dimension(3,3) :: Imom
      this%film_type=0.0_WP
      this%d1= 0.0_WP;this%d2=0.0_WP;this%d3=0.0_WP
      this%l1= 0.0_WP;this%l2=0.0_WP;this%l3=0.0_WP
      this%d3_d1= 0.0_WP;this%d3_d2=0.0_WP;this%d2_d1=0.0_WP
      tmpfilm_type = 0; ratio=1.5_WP
      call dsyev('V','U',order,Imom,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      ! Traverse domain and compute sensors
      do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
         do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
            do i=this%vf%cfg%imin_,this%vf%cfg%imax_
               Imom = 0.0_WP; lvol = 0.0_WP; lx_vol=0.0_WP
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        ! Volume
                        lvol = lvol + this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        ! Center of gravity
                        lx_vol = lx_vol + this%vf%Lbary(:,ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                     end do
                  end do
               end do

               lx_vol = lx_vol/lvol
               do kk = k-2,k+2
                  do jj = j-2,j+2
                     do ii = i-2,i+2
                        ! Location of film node
                        Ltmp = this%vf%Lbary(:,ii,jj,kk) - lx_vol
                        Imom(1,1) = Imom(1,1) + (Ltmp(2)**2 + Ltmp(3)**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        Imom(2,2) = Imom(2,2) + (Ltmp(1)**2 + Ltmp(3)**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        Imom(3,3) = Imom(3,3) + (Ltmp(1)**2 + Ltmp(2)**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        
                        Imom(1,2) = Imom(1,2) - Ltmp(1)*Ltmp(2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        Imom(1,3) = Imom(1,3) - Ltmp(1)*Ltmp(3)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
                        Imom(2,3) = Imom(2,3) - Ltmp(2)*Ltmp(3)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)   
                     end do
                  end do
               end do
               call dsyev('V','U',order,Imom,order,d,work,lwork,info)
               d = max(0.0_WP,d)
      
               ! if (lvol .gt. 0.0_WP) then
               if (d(3).gt.(ratio*d(1))) tmpfilm_type(i,j,k) = tmpfilm_type(i,j,k) + 1
               if (d(3).gt.(ratio*d(2))) tmpfilm_type(i,j,k) = tmpfilm_type(i,j,k) + 1

               if (d(3).gt.(ratio*d(1))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1.0_WP
               if (d(3).gt.(ratio*d(2))) this%film_type(i,j,k) = this%film_type(i,j,k) + 1.0_WP

               this%d1(i,j,k)=d(1)
               this%d2(i,j,k)=d(2)
               this%d3(i,j,k)=d(3)
               this%l1(i,j,k)=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/(lvol+tiny(1.0_WP)))
               this%l2(i,j,k)=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/(lvol+tiny(1.0_WP)))
               this%l3(i,j,k)=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/(lvol+tiny(1.0_WP)))
               ! print *,d(3)/(d(1)+tiny(1.0_WP)),d(3)/(d(2)+tiny(1.0_WP)),d(2)/(d(1)+tiny(1.0_WP))
               this%d3_d1(i,j,k)=d(3)/(d(1)+tiny(1.0_WP))
               this%d3_d2(i,j,k)=d(3)/(d(2)+tiny(1.0_WP))
               this%d2_d1(i,j,k)=d(2)/(d(1)+tiny(1.0_WP))
               ! this%film_type(i,j,k) = tmpfilm_type(i,j,k)
               ! end if
            end do 
         end do 
      end do 
      deallocate(work)
      call this%vf%cfg%sync(tmpfilm_type)
      call this%vf%cfg%sync(this%film_type)
   end subroutine get_localfilmtype


   subroutine get_cclstats(this,ccl,x,y,z,u,v,w,vol,lengths,maxlength,axes,f_ligament)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX,MPI_INTEGER
      use parallel,  only: MPI_REAL_WP
      use irl_fortran_interface
      implicit none
      class(simplex), intent(inout) :: this
      class(cclabel), intent(in)::ccl
      real(WP), dimension(1:), intent(inout) :: x,y,z,u,v,w,vol,maxlength,f_ligament
      real(WP), dimension(1:,1:), intent(inout) :: lengths
      real(WP), dimension(1:,1:,1:), intent(inout) :: axes
      integer :: n,nn,i,j,k,ii,jj,kk,ierr,np,ip,m,iunit,rank,per_x,per_y,per_z
      ! Allocate variables to get stats
      real(WP), dimension(:), allocatable :: vol_,x_vol_,y_vol_,z_vol_
      real(WP), dimension(:), allocatable :: u_vol_,v_vol_,w_vol_
      real(WP), dimension(:), allocatable :: x_min_,x_min,x_max_,x_max
      real(WP), dimension(:), allocatable :: y_min_,y_min,y_max_,y_max
      real(WP), dimension(:), allocatable :: z_min_,z_min,z_max_,z_max
      real(WP), dimension(:,:,:), allocatable :: Imom_,Imom
      ! For ligament type
      real(WP), dimension(:), allocatable :: ncell_,ncell,n_ligament_,n_ligament
      real(WP) :: xtmp,ytmp,ztmp
      ! Moment of inertia variable
      real(WP), dimension(:), allocatable :: work
      real(WP), dimension(1)   :: lwork_query
      real(WP), dimension(3) :: d
      real(WP), dimension(3,3) :: A
      integer , parameter :: order = 3
      integer  :: lwork,info
      allocate(vol_(1:ccl%nstruct));vol_=0.0_WP
      allocate(x_min(1:ccl%nstruct),x_min_(1:ccl%nstruct));x_min=0.0_WP;x_min_= 10000.0_WP
      allocate(x_max(1:ccl%nstruct),x_max_(1:ccl%nstruct));x_max=0.0_WP;x_max_=-10000.0_WP
      allocate(y_min(1:ccl%nstruct),y_min_(1:ccl%nstruct));y_min=0.0_WP;y_min_= 10000.0_WP
      allocate(y_max(1:ccl%nstruct),y_max_(1:ccl%nstruct));y_max=0.0_WP;y_max_=-10000.0_WP
      allocate(z_min(1:ccl%nstruct),z_min_(1:ccl%nstruct));z_min=0.0_WP;z_min_= 10000.0_WP
      allocate(z_max(1:ccl%nstruct),z_max_(1:ccl%nstruct));z_max=0.0_WP;z_max_=-10000.0_WP
      allocate(x_vol_(1:ccl%nstruct),y_vol_(1:ccl%nstruct),z_vol_(1:ccl%nstruct));x_vol_=0.0_WP;y_vol_=0.0_WP;z_vol_=0.0_WP
      allocate(u_vol_(1:ccl%nstruct),v_vol_(1:ccl%nstruct),w_vol_(1:ccl%nstruct));u_vol_=0.0_WP;v_vol_=0.0_WP;w_vol_=0.0_WP
      allocate(Imom(1:ccl%nstruct,3,3),Imom_(1:ccl%nstruct,3,3));Imom=0.0_WP;Imom_=0.0_WP
      allocate(ncell(1:ccl%nstruct),ncell_(1:ccl%nstruct),n_ligament_(1:ccl%nstruct),n_ligament(1:ccl%nstruct))
      ncell=0.0_WP;ncell_=0.0_WP;n_ligament=0.0_WP;n_ligament_=0.0_WP
      ! Query optimal work array size
      call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))
      do i=1,ccl%nstruct
         ! Periodicity
         per_x = ccl%struct(i)%per(1); per_y = ccl%struct(i)%per(2); per_z = ccl%struct(i)%per(3)
         ! get number of local cells
         ncell_(i) = ccl%struct(i)%n_*1.0_WP
         do j=1,ccl%struct(i)%n_
            ii=ccl%struct(i)%map(1,j); jj=ccl%struct(i)%map(2,j); kk=ccl%struct(i)%map(3,j)
            ! Location of struct node
            xtmp = this%vf%cfg%xm(ii)-per_x*this%vf%cfg%xL
            ytmp = this%vf%cfg%ym(jj)-per_y*this%vf%cfg%yL
            ztmp = this%vf%cfg%zm(kk)-per_z*this%vf%cfg%zL
            ! Volume
            vol_(i) = vol_(i) + this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            ! Center of gravity
            x_vol_(i) = x_vol_(i) + xtmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            y_vol_(i) = y_vol_(i) + ytmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            z_vol_(i) = z_vol_(i) + ztmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            ! Average gas velocity inside struct
            u_vol_(i) = u_vol_(i) + this%fs%U(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            v_vol_(i) = v_vol_(i) + this%fs%V(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            w_vol_(i) = w_vol_(i) + this%fs%W(ii,jj,kk)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            if(tmpfilm_type(ii,jj,kk).eq.1) n_ligament_(i)=n_ligament_(i)+1.0_WP
         end do
      end do
      ! Sum parallel stats
      call MPI_ALLREDUCE(vol_,vol,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(x_vol_,x,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_vol_,y,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_vol_,z,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(u_vol_,u,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(v_vol_,v,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(w_vol_,w,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(ncell_,ncell,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(n_ligament_,n_ligament,ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      do i=1,ccl%nstruct
         ! Periodicity
         per_x = ccl%struct(i)%per(1); per_y = ccl%struct(i)%per(2); per_z = ccl%struct(i)%per(3)
         do j=1,ccl%struct(i)%n_
            ! Indices of struct node
            ii=ccl%struct(i)%map(1,j); jj=ccl%struct(i)%map(2,j); kk=ccl%struct(i)%map(3,j)
            xtmp = this%vf%cfg%xm(ii)-per_x*this%vf%cfg%xL-x(i)/vol(i)
            ytmp = this%vf%cfg%ym(jj)-per_y*this%vf%cfg%yL-y(i)/vol(i)
            ztmp = this%vf%cfg%zm(kk)-per_z*this%vf%cfg%zL-z(i)/vol(i)
            ! Moment of Inertia
            Imom_(i,1,1) = Imom_(i,1,1) + (ytmp**2 + ztmp**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            Imom_(i,2,2) = Imom_(i,2,2) + (xtmp**2 + ztmp**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            Imom_(i,3,3) = Imom_(i,3,3) + (xtmp**2 + ytmp**2)*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            Imom_(i,1,2) = Imom_(i,1,2) - xtmp*ytmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            Imom_(i,1,3) = Imom_(i,1,3) - xtmp*ztmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            Imom_(i,2,3) = Imom_(i,2,3) - ytmp*ztmp*this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)
            do n=1,2
               if (getNumberOfVertices(this%vf%interface_polygon(n,ii,jj,kk)).gt.0) then
                  d = calculateCentroid(this%vf%interface_polygon(n,ii,jj,kk))
                  x_min_(i) = min(x_min_(i),d(1)); x_max_(i) = max(x_max_(i),d(1))
                  y_min_(i) = min(y_min_(i),d(2)); y_max_(i) = max(y_max_(i),d(2))
                  z_min_(i) = min(z_min_(i),d(3)); z_max_(i) = max(z_max_(i),d(3))
               end if
            end do
            ! ! Min thickness
            ! min_thickness_(i) = min(min_thickness_(i),this%struct_thickness(ii,jj,kk))
         end do
      end do
      ! Sum parallel stat on Imom
      do i=1,3
         do j=1,3
            call MPI_ALLREDUCE(Imom_(:,i,j),Imom(:,i,j),ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
         end do
      end do
      ! Get extents
      call MPI_ALLREDUCE(x_min_,x_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(x_max_,x_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_min_,y_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(y_max_,y_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_min_,z_min,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      call MPI_ALLREDUCE(z_max_,z_max,ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      !!  Get min thickness
      !  call MPI_ALLREDUCE(min_thickness_,min_thickness,this%ccl_ligament%nstruct,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      ! Store data
      do i=1,ccl%nstruct
         ! Center of gravity
         x(i) = x(i)/vol(i); y(i) = y(i)/vol(i); z(i) = z(i)/vol(i)
         ! Periodicity: transport back inside domain if needed
         if (x(i).lt.this%vf%cfg%x(this%vf%cfg%imin)) x(i) = x(i)+this%vf%cfg%xL
         if (y(i).lt.this%vf%cfg%y(this%vf%cfg%jmin)) y(i) = y(i)+this%vf%cfg%yL
         if (z(i).lt.this%vf%cfg%z(this%vf%cfg%kmin)) z(i) = z(i)+this%vf%cfg%zL
         u(i)=u(i)/vol(i); v(i)=v(i)/vol(i); w(i)=w(i)/vol(i)
         ! maxlength(i) = hypot(hypot(x_max(i)-x_min(i),y_max(i)-y_min(i))**2,z_max(i)-z_min(i))
         maxlength(i) = sqrt((x_max(i)-x_min(i))**2+(y_max(i)-y_min(i))**2+(z_max(i)-z_min(i))**2)
         ! Eigenvalues/eigenvectors of moments of inertia tensor
         A = Imom(i,:,:); n = 3
         ! On exit, A contains eigenvectors, and d contains eigenvalues in ascending order
         call dsyev('V','U',n,A,n,d,work,lwork,info)
         ! Get rid of very small negative values (due to machine accuracy)
         d = max(0.0_WP,d)
         ! Store characteristic lengths
         lengths(i,1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol(i))
         lengths(i,2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol(i))
         lengths(i,3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol(i))
         ! Zero out length in 3rd dimension if 2D
         if (this%vf%cfg%nx.eq.1.or.this%vf%cfg%ny.eq.1.or.this%vf%cfg%nz.eq.1) lengths(i,3)=0.0_WP
         ! Store principal axes
         axes(i,:,:) = A
         ! Use max of bounding box and MoI-derived lengths as length
         maxlength(i) = max(maxlength(i),lengths(i,1))
         if (ncell(i).eq.0) then
            f_ligament(i) = 0.0_WP
         else
            f_ligament(i) = 1.0_WP*n_ligament(i)/(1.0_WP*ncell(i))
         end if
      end do
      ! Deallocate arrays
      deallocate(vol_,x_vol_,y_vol_,z_vol_,u_vol_,v_vol_,w_vol_,Imom_,Imom)
      deallocate(x_min_,y_min_,z_min_,x_max_,y_max_,z_max_)
      deallocate(x_min,y_min,z_min,x_max,y_max,z_max)
      deallocate(ncell,ncell_,n_ligament,n_ligament_)
   end subroutine get_cclstats

   subroutine get_structminthickness(this,ccl,min_thickness,type)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MIN,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      implicit none
      class(simplex), intent(inout) :: this 
      class(cclabel), intent(in)::ccl
      real(WP) , dimension(1:), intent(inout) :: min_thickness
      integer, intent(in) :: type
      integer :: n,nn,i,j,k,ierr
      real(WP), dimension(:), allocatable :: min_thickness_
      allocate(min_thickness_(1:ccl%nstruct)); min_thickness_ = 5.0_WP*this%vf%cfg%min_meshsize
      call this%vf%get_thickness()
      do n=1,ccl%nstruct
         do nn=1,ccl%struct(n)%n_
            i=ccl%struct(n)%map(1,nn); j=ccl%struct(n)%map(2,nn); k=ccl%struct(n)%map(3,nn)
            if (type.eq.1) then
               min_thickness_(n) = min(min_thickness_(n),tmpthickness(i,j,k))
            else 
               min_thickness_(n) = min(min_thickness_(n),this%vf%thickness(i,j,k))
            end if
         end do
      end do
      call MPI_ALLREDUCE(min_thickness_,min_thickness,ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
      deallocate(min_thickness_)
   end subroutine get_structminthickness
   
   
end module simplex_class