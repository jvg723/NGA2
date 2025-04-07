!> Various definitions and tools for running an NGA2 simulation
module simulation
   use mpi_f08,         only: MPI_Group
   use precision,       only: WP
   use simplex_class,   only: simplex
   use atom_class,      only: atom
   use coupler_class,   only: coupler
   use inputfile_class, only: inputfile
   implicit none
   private
   
   !> Simplex simulation
   type(simplex) :: spx

   !> Atomization simulation
   type(atom) :: atomization

   !> Input files to read in patitions
   type(inputfile) :: input_spx, input_atom

   !> Couplers from simplex to atomization
   type(coupler) :: xcpl_s2a,ycpl_s2a,zcpl_s2a !> Velocity
   type(coupler) :: vfcpl_s2a

   !> for VOF coupling
   integer :: vof_couple_max=5    !<number of cells over which VOF is coupled between domains
   integer :: vof_couple_min=1    !<number of cells over which VOF is coupled between domains

   !> For MPI groups
   public :: group_spx,isInGrp_spx
   public :: group_atom,isInGrp_atom
   integer, dimension(3) :: partition_spx,partition_atom
   logical :: isInGrp_spx,isInGrp_atom
   type(MPI_Group) :: group_spx,group_atom
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use parallel,    only: comm,group,nproc,rank,amRoot
      use mpi_f08,     only: MPI_Group,MPI_Group_range_incl
      implicit none

      ! mpi_groups: block
      !    integer, dimension(3,1) :: grange
      !    integer :: ierr
      !    ! Read in partition
      !    input_spx =inputfile(amRoot=amRoot,filename='simplex.input')
      !    input_atom=inputfile(amRoot=amRoot,filename='atomization.input')
      !    call input_spx%read('Partition',partition_spx)
      !    call input_atom%read('Partition',partition_atom)
      !    ! Create an MPI group along with logical for the simplex nozzle on the lowest ranks
      !    grange(:,1)=[0,product(partition_spx)-1,1]
      !    call MPI_Group_range_incl(group,1,grange,group_spx,ierr)
      !    isInGrp_spx=.false.; if (rank.le.product(partition_spx)-1) isInGrp_spx=.true.
      !    ! Create an MPI group along with logical for the atomization domain on the highest ranks
      !    grange(:,1)=[nproc-product(partition_atom),nproc-1,1]
      !    call MPI_Group_range_incl(group,1,grange,group_atom,ierr)
      !    isInGrp_atom=.false.; if (rank.ge.nproc-product(partition_atom)) isInGrp_atom=.true.
      ! end block mpi_groups
      
      ! Initialize simplex simulation
      call spx%init()

      ! Initialize atomization simulation
      call atomization%init()

      ! Initialize couplers from injector to atomization
      create_coupler_s2a: block
         use parallel, only: group
         ! Setup couplers
         xcpl_s2a=coupler(src_grp=group,dst_grp=group,name='simplex_to_atom');  call xcpl_s2a%set_src(spx%cfg,'x');  call xcpl_s2a%set_dst(atomization%cfg,'x');  call xcpl_s2a%initialize()
         ycpl_s2a=coupler(src_grp=group,dst_grp=group,name='simplex_to_atom');  call ycpl_s2a%set_src(spx%cfg,'y');  call ycpl_s2a%set_dst(atomization%cfg,'y');  call ycpl_s2a%initialize()
         zcpl_s2a=coupler(src_grp=group,dst_grp=group,name='simplex_to_atom');  call zcpl_s2a%set_src(spx%cfg,'z');  call zcpl_s2a%set_dst(atomization%cfg,'z');  call zcpl_s2a%initialize()
         vfcpl_s2a=coupler(src_grp=group,dst_grp=group,name='simplex_to_atom'); call vfcpl_s2a%set_src(spx%cfg,'c'); call vfcpl_s2a%set_dst(atomization%cfg,'c'); call vfcpl_s2a%initialize()
      end block create_coupler_s2a

      ! Create region to couple VOF between domains
      create_coupler_region: block
         atomization%vfcouple_xmin=atomization%cfg%xm(atomization%cfg%imin+vof_couple_min)
         atomization%vfcouple_xmax=atomization%cfg%xm(atomization%cfg%imin+vof_couple_max)
      end block create_coupler_region

      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none


      ! ! Simplex drives overall time integration
      ! do while (.not.spx%time%done())
         
      !    call spx%step()

      ! end do
      
      ! Simplex drives overall time integration
      do while (.not.spx%time%done())
         
         ! Advance simplex simulation
         call spx%step()

         ! Handle coupling velocity between simplex and atomization
         coupling_velocity_s2a: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cpl12x/y/z couplers
            call xcpl_s2a%push(spx%fs%U);   call xcpl_s2a%transfer();  call xcpl_s2a%pull(atomization%resU)
            call ycpl_s2a%push(spx%fs%V);   call ycpl_s2a%transfer();  call ycpl_s2a%pull(atomization%resV)
            call zcpl_s2a%push(spx%fs%W);   call zcpl_s2a%transfer();  call zcpl_s2a%pull(atomization%resW)
            call atomization%fs%get_bcond('inlets',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               atomization%fs%U(i  ,j,k)=atomization%resU(i  ,j,k)*sum(atomization%fs%itpr_x(:,i  ,j,k)*atomization%cfg%VF(i-1:i,    j,    k))
               atomization%fs%V(i-1,j,k)=atomization%resV(i-1,j,k)*sum(atomization%fs%itpr_y(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j-1:j,    k))
               atomization%fs%W(i-1,j,k)=atomization%resW(i-1,j,k)*sum(atomization%fs%itpr_z(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j    ,k-1:k))
            end do
         end block coupling_velocity_s2a


         ! Handle coupling VOF between simplex and atomization
         coupling_vof_s2a: block
            integer :: i,j,k
            ! Exchange data using cell center coupler
            atomization%tempVF=0.0_WP
            call vfcpl_s2a%push(spx%vf%VF); call vfcpl_s2a%transfer(); call vfcpl_s2a%pull(atomization%tempVF)
            ! Exchange VOF based upon x/y/z position
            do k=atomization%vf%cfg%kmino_,atomization%vf%cfg%kmaxo_
               do j=atomization%vf%cfg%jmino_,atomization%vf%cfg%jmaxo_
                  do i=atomization%vf%cfg%imino_,atomization%vf%cfg%imaxo_
                     if (atomization%vf%cfg%xm(i).ge.atomization%vfcouple_xmin.and.atomization%vf%cfg%xm(i).le.atomization%vfcouple_xmax) then
                        atomization%vf%VF(i,j,k)=atomization%tempVF(i,j,k)
                     end if
                  end do
               end do
            end do
            !> sync arrays
            call atomization%cfg%sync(atomization%vf%VF)
            !> Reconstruct interface to get updated VF and moments
            call atomization%vf%build_interface()
            call atomization%vf%reset_volume_moments()
         end block coupling_vof_s2a
      
         ! Advance atomization simulation until it's caught up
         do while (atomization%time%t.le.spx%time%t)
            call atomization%step()
         end do

      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize simplex simulation
      call spx%final()

      ! Finalize atomization simulation
      call atomization%final()
      
   end subroutine simulation_final

      
   ! !> Transfer droplet to Lagrangian representation
   ! subroutine transfer_drops(this)
   !    use mpi_f08,   
   !    use parallel,  only: MPI_REAL_WP
   !    use mathtools, only: pi,twoPi
   !    use irl_fortran_interface
   !    class(simplex), intent(inout) :: this
   !    real(WP), dimension(:)    , allocatable :: dvol
   !    real(WP), dimension(:,:)  , allocatable :: dpos
   !    real(WP), dimension(:,:)  , allocatable :: dvel
   !    real(WP), dimension(:,:,:), allocatable :: dmoi
   !    real(WP), dimension(:)    , allocatable :: drem
   !    real(WP), dimension(:)    , allocatable :: dlen
   !    real(WP), dimension(:)    , allocatable :: decc
   !    real(WP), dimension(:)    , allocatable :: xmin,xmax,ymin,ymax,zmin,zmax
   !    integer :: n,m,ierr,i,j,k,l,nmax,np_start
   !    real(WP) :: x,y,z,x0,y0,z0,diam,lmax,lmid,lmin
   !    logical :: transfer
   !    integer :: nmain,nsat
   !    real(WP) :: Vt,Vl,Vd,minor_radius,Vrim,Lrim,Lrp
   !    ! Moment of inertia calculation using lapack
   !    real(WP), dimension(:), allocatable, save :: work !< Saved!
   !    integer, save :: lwork                            !< Saved!
   !    real(WP), dimension(1) :: lwork_query
   !    real(WP), dimension(3) :: d
   !    real(WP), dimension(3,3) :: A
   !    integer :: info
      
   !    ! Query optimal work array size
   !    if (.not.allocated(work)) then
   !       call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
   !       lwork=int(lwork_query(1)); allocate(work(lwork))
   !    end if
      
   !    ! Start by performing a CCL
   !    call this%ccl%build(make_label,same_label)
      
   !    ! Allocate droplet stats arrays
   !    allocate(dvol(1:this%ccl%nstruct        )); dvol=0.0_WP
   !    allocate(dpos(1:this%ccl%nstruct,1:3    )); dpos=0.0_WP
   !    allocate(dvel(1:this%ccl%nstruct,1:3    )); dvel=0.0_WP
   !    allocate(dmoi(1:this%ccl%nstruct,1:3,1:3)); dmoi=0.0_WP
   !    allocate(drem(1:this%ccl%nstruct        )); drem=0.0_WP
   !    allocate(dlen(1:this%ccl%nstruct        )); dlen=0.0_WP
   !    allocate(decc(1:this%ccl%nstruct)); decc=0.0_WP
   !    allocate(xmin(1:this%ccl%nstruct),xmax(1:this%ccl%nstruct)); xmin=HUGE(x);xmax=-HUGE(x)
   !    allocate(ymin(1:this%ccl%nstruct),ymax(1:this%ccl%nstruct)); ymin=HUGE(x);ymax=-HUGE(x)
   !    allocate(zmin(1:this%ccl%nstruct),zmax(1:this%ccl%nstruct)); zmin=HUGE(x);zmax=-HUGE(x)
      
   !    ! First pass to accumulate volume, position, and velocity
   !    do n=1,this%ccl%nstruct
   !       ! Loop over cells in structure
   !       do m=1,this%ccl%struct(n)%n_
   !          ! Get cell indices
   !          i=this%ccl%struct(n)%map(1,m)
   !          j=this%ccl%struct(n)%map(2,m)
   !          k=this%ccl%struct(n)%map(3,m)
   !          ! Get cell position, accounting for periodicity
   !          x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL
   !          y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL
   !          z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL
   !          ! Accumulate volume, position, and velocity
   !          dvol(n  )=dvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
   !          dpos(n,:)=dpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
   !          dvel(n,:)=dvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
   !          ! Check if drop touches auto-transfer layer
   !          if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
   !          &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
   !          &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
   !          &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
   !          &   k.ge.this%vf%cfg%kmax-this%nlayer) drem(n)=1.0_WP
   !          ! Get the structures's locally largest and smallest x,y,z locations
   !          do l=1,2
   !             if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
   !                d = calculateCentroid(this%vf%interface_polygon(l,i,j,k))
   !                xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
   !                ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
   !                zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
   !             end if
   !          end do
   !       end do
   !    end do
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,1*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,dpos,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,dvel,3*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,drem,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1*this%ccl%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
      
   !    ! Second pass to accumulate moment of inertia
   !    do n=1,this%ccl%nstruct
   !       ! Get drop barycenter
   !       x0=dpos(n,1)/dvol(n)
   !       y0=dpos(n,2)/dvol(n)
   !       z0=dpos(n,3)/dvol(n)
   !       ! Loop over cells in structure
   !       do m=1,this%ccl%struct(n)%n_
   !          ! Get cell indices
   !          i=this%ccl%struct(n)%map(1,m)
   !          j=this%ccl%struct(n)%map(2,m)
   !          k=this%ccl%struct(n)%map(3,m)
   !          ! Get cell position relative to drop barycenter, accounting for periodicity
   !          x=this%vf%cfg%xm(i)-this%ccl%struct(n)%per(1)*this%vf%cfg%xL-x0
   !          y=this%vf%cfg%ym(j)-this%ccl%struct(n)%per(2)*this%vf%cfg%yL-y0
   !          z=this%vf%cfg%zm(k)-this%ccl%struct(n)%per(3)*this%vf%cfg%zL-z0
   !          ! Accumulate moment of inertia
   !          dmoi(n,1,1)=dmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
   !          dmoi(n,2,2)=dmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
   !          dmoi(n,3,3)=dmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
   !          dmoi(n,1,2)=dmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
   !          dmoi(n,1,3)=dmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
   !          dmoi(n,2,3)=dmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
   !       end do
   !    end do
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,dmoi,9*this%ccl%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      
   !    ! Third pass to generate normalized drop stats
   !    do n=1,this%ccl%nstruct
   !       ! Get drop barycenter, accounting for periodicity
   !       dpos(n,:)=dpos(n,:)/dvol(n)
   !       if (this%vf%cfg%xper.and.dpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) dpos(n,1)=dpos(n,1)+this%vf%cfg%xL
   !       if (this%vf%cfg%yper.and.dpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) dpos(n,2)=dpos(n,2)+this%vf%cfg%yL
   !       if (this%vf%cfg%zper.and.dpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) dpos(n,3)=dpos(n,3)+this%vf%cfg%zL
   !       ! Get drop velocity
   !       dvel(n,:)=dvel(n,:)/dvol(n)
   !       ! Calculate maximum length of the structure
   !       A=dmoi(n,:,:)
   !       call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
   !       d=max(0.0_WP,d)    
   !       ! Replace with corrected eigenvectors for future ligament droplet placement
   !       dmoi(n,:,:)=A
   !       ! Get characteristic lengths of drop
   !       lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/dvol(n))
   !       lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/dvol(n))
   !       lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/dvol(n))
   !       if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
   !       ! Compute eccentricity
   !       decc(n)=sqrt(1.0_WP-lmin**2/(lmax**2+epsilon(1.0_WP)))
   !       ! Use max of bounding box and MoI-derived lengths as length
   !       dlen(n)=max(sqrt((xmax(n)-xmin(n))**2+(ymax(n)-ymin(n))**2+(zmax(n)-zmin(n))**2),lmax)
   !    end do
      
   !    ! Find the liquid core
   !    nmax=maxloc(dvol,dim=1)
      
   !    ! Zero out monitoring variables
   !    this%vof_tf_drop=0.0_WP
   !    this%vof_deleted=0.0_WP
   !    this%np_drop=0
      
   !    ! Transfer drops based on our criteria
   !    do n=1,this%ccl%nstruct
         
   !       ! Compute diameter
   !       diam=(6.0_WP*dvol(n)/pi)**(1.0_WP/3.0_WP)
         
   !       ! Decide if struct is being transfered
   !       if (diam.gt.this%dmax) then
   !          ! Too big to transfer
   !          transfer=.false.
   !       else if (diam.le.this%ddel) then
   !          ! Too small to track, delete immediately
   !          transfer=.false.
   !          ! Zero out VF in the structure
   !          do m=1,this%ccl%struct(n)%n_
   !             this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
   !          end do
   !          ! Increment monitoring variables
   !          this%vof_deleted=this%vof_deleted+dvol(n)
   !       else if (diam.gt.this%ddel.and.diam.le.this%dmin) then
   !          ! Small enough to transfer automatically
   !          transfer=.true.
   !       else
   !          ! In between, check eccentricity from moment of inertia tensor
   !          if (decc(n).gt.this%emax) then
   !             ! Too eccentric to transfer yet
   !             transfer=.false.
   !          else
   !             ! Spherical enough to transfer
   !             transfer=.true.
   !          end if
   !       end if

   !       ! Force transfer if drop touches auto-transfer layer
   !       if (drem(n).gt.0.0_WP) transfer=.true.
         
   !       ! But prevent transfer if that's the core
   !       if (n.eq.nmax) transfer=.false.
         
   !       ! Perform transfer
   !       if (transfer) then

   !          if (dlen(n).le.0.0_WP) cycle
            
   !          if (decc(n).gt.this%emax) then !> convert as a ligament
   !             !>Break-up as ligament (Drop size method from Kim & Moin (2020))
   !             Lrim=dlen(n) 
   !             Vrim=dvol(n)
   !             minor_radius=sqrt(Vrim/pi/Lrim) 
   !             nmain=floor(this%dw*Lrim/(twoPi*minor_radius))
   !             nsat=nmain+1
   !             diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)
   !             ! Only the main processor is in charge of creating droplets
   !             if (this%cfg%amRoot) then
   !                Lrp = twoPi*minor_radius/this%dw
   !                ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " Lrp=",Lrp," Lrim=", Lrim, " Vrim=",Vrim, " minor_radius=", minor_radius, " diam=", diam
   !                ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " pre-conversion"," Lrp=",Lrp
   !                do l=1,nsat+nmain
   !                   ! Increment particle counter
   !                   this%lp%np_=this%lp%np_+1
   !                   ! Make room for new drop
   !                   call this%lp%resize(this%lp%np_)
   !                   ! Add the drop
   !                   this%lp%p(this%lp%np_)%id  =int(2,8)                                                                               
   !                   if (mod(l,2).eq.1) then
   !                      this%lp%p(this%lp%np_)%d=diam*this%size_ratio                                                                                    
   !                   else
   !                      this%lp%p(this%lp%np_)%d=diam                                                                                    
   !                   end if
   !                   this%lp%p(this%lp%np_)%pos=dpos(n,:)+0.5_WP*Lrp*(l-(nmain+1))*dmoi(n,:,1)  
   !                   ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " dpos(n,1)=",dpos(n,1), " dpos(n,2)=",dpos(n,2), " dpos(n,3)=",dpos(n,3)
   !                   ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " dmoi(n,1,1)=",dmoi(n,1,1), " dmoi(n,2,1)=",dmoi(n,2,1), " dmoi(n,3,1)=",dmoi(n,3,1)
   !                   ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " pos(1)=",this%lp%p(this%lp%np_)%pos(1), " pos(2)=",this%lp%p(this%lp%np_)%pos(2), " pos(3)=",this%lp%p(this%lp%np_)%pos(3)   
   !                   this%lp%p(this%lp%np_)%vel=dvel(n,:)
   !                   ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " dvel(n,1)=",dvel(n,1), " dvel(n,2)=",dvel(n,2), " dvel(n,3)=",dvel(n,3)
   !                   this%lp%p(this%lp%np_)%ind=this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])  
   !                   ! if (this%vf%cfg%amRoot) print *, "ligament This paritcle id id=", n, " ind(1)=",this%lp%p(this%lp%np_)%ind(1), " ind(2)=",this%lp%p(this%lp%np_)%ind(2), " ind(3)=",this%lp%p(this%lp%np_)%ind(3)   
   !                   if (ABS(this%lp%p(this%lp%np_)%pos(1)).ge.this%cfg%xL/2.00_WP.or.ABS(this%lp%p(this%lp%np_)%pos(2)).ge.this%cfg%yL/2.00_WP.or.ABS(this%lp%p(this%lp%np_)%pos(3)).ge.this%cfg%zL/2.00_WP) then
   !                      this%lp%p(this%lp%np_)%flag=1 
   !                   else
   !                      this%lp%p(this%lp%np_)%flag=0 
   !                   end if                                                                                
   !                   this%lp%p(this%lp%np_)%dt  =0.0_WP                                                                                  
   !                   this%lp%p(this%lp%np_)%Acol=0.0_WP                                                                                  
   !                   this%lp%p(this%lp%np_)%Tcol=0.0_WP
   !                end do
   !                ! Increment monitoring variables
   !                this%vof_tf_drop=this%vof_tf_drop+dvol(n)
   !                this%np_drop=this%np_drop+nmain+nsat
   !                this%lp%np_new=this%lp%np_new+nmain+nsat
   !                this%lp%vp_new=this%lp%vp_new+dvol(n)
   !             end if
   !          else !> convert to drop
   !             ! Root creates a new Lagrangian drop
   !             if (this%vf%cfg%amRoot) then
   !                np_start=this%lp%np_
   !                ! Increment particle counter
   !                this%lp%np_=this%lp%np_+1
   !                ! Make room for new drop
   !                call this%lp%resize(this%lp%np_)
   !                ! Add the drop
   !                this%lp%p(this%lp%np_)%id  =int(1,8)
   !                this%lp%p(this%lp%np_)%d   =diam
   !                this%lp%p(this%lp%np_)%pos =dpos(n,:)
   !                ! if (this%vf%cfg%amRoot) print *, "drop This paritcle id id=", n, " dpos(n,1)=",dpos(n,1), " dpos(n,2)=",dpos(n,2), " dpos(n,3)=",dpos(n,3)
   !                ! if (this%vf%cfg%amRoot) print *, "drop This paritcle id id=", n, " dmoi(n,1,1)=",dmoi(n,1,1), " dmoi(n,2,1)=",dmoi(n,2,1), " dmoi(n,3,1)=",dmoi(n,3,1)
   !                ! if (this%vf%cfg%amRoot) print *, "drop This paritcle id id=", n, " pos(1)=",this%lp%p(this%lp%np_)%pos(1), " pos(2)=",this%lp%p(this%lp%np_)%pos(2), " pos(3)=",this%lp%p(this%lp%np_)%pos(3)   
   !                this%lp%p(this%lp%np_)%vel =dvel(n,:)
   !                ! if (this%vf%cfg%amRoot) print *, "drop This paritcle id id=", n, " dvel(n,1)=",dvel(n,1), " dvel(n,2)=",dvel(n,2), " dvel(n,3)=",dvel(n,3)
   !                this%lp%p(this%lp%np_)%ind =this%lp%cfg%get_ijk_global(dpos(n,:),[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])
   !                ! if (this%vf%cfg%amRoot) print *, "drop This paritcle id id=", n, " ind(1)=",this%lp%p(this%lp%np_)%ind(1), " ind(2)=",this%lp%p(this%lp%np_)%ind(2), " ind(3)=",this%lp%p(this%lp%np_)%ind(3)   
   !                this%lp%p(this%lp%np_)%flag=0
   !                this%lp%p(this%lp%np_)%dt  =0.0_WP
   !                this%lp%p(this%lp%np_)%Acol=0.0_WP
   !                this%lp%p(this%lp%np_)%Tcol=0.0_WP
   !             end if
   !             ! Increment monitoring variables
   !             this%vof_tf_drop=this%vof_tf_drop+dvol(n)
   !             this%np_drop=this%np_drop+1
   !             this%lp%np_new=this%lp%np_new+1
   !             this%lp%vp_new=this%lp%vp_new+dvol(n)
   !          end if
            
   !          ! Zero out VF in the structure
   !          do m=1,this%ccl%struct(n)%n_
   !             this%vf%VF(this%ccl%struct(n)%map(1,m),this%ccl%struct(n)%map(2,m),this%ccl%struct(n)%map(3,m))=0.0_WP
   !          end do

   !       end if
         
   !    end do
      
   !    ! Synchronize VF fields
   !    call this%vf%sync_interface()
   !    call this%vf%clean_irl_and_band()
      
   !    ! Synchronize particles
   !    call this%lp%sync()
      
   !    ! Deallocate all but work array
   !    deallocate(dvol,dpos,dvel,dmoi,drem)
   !    deallocate(dlen,xmin,ymin,zmin,decc)
      
   ! contains
      
   !    !> Function that identifies cells that need a label
   !    logical function make_label(i,j,k)
   !       implicit none
   !       integer, intent(in) :: i,j,k
   !       if (this%vf%VF(i,j,k).gt.0.0_WP) then
   !          make_label=.true.
   !       else
   !          make_label=.false.
   !       end if
   !    end function make_label
      
   !    !> Function that identifies if cell pairs have same label
   !    logical function same_label(i1,j1,k1,i2,j2,k2)
   !       implicit none
   !       integer, intent(in) :: i1,j1,k1,i2,j2,k2
   !       same_label=.true.
   !    end function same_label
      
   ! end subroutine transfer_drops

   ! subroutine transfer_ligs(this)
   !    use vfs_class, only: VFlo,VFhi
   !    use mathtools, only: pi,twoPi
   !    use mpi_f08
   !    use parallel,  only: MPI_REAL_WP
   !    use messager, only: die
   !    use irl_fortran_interface
   !    implicit none
   !    class(simplex), intent(inout) :: this
   !    real(WP), dimension(:)    , allocatable :: lvol
   !    real(WP), dimension(:)    , allocatable :: lthc
   !    real(WP), dimension(:)    , allocatable :: llen
   !    real(WP), dimension(:)    , allocatable :: lnum
   !    real(WP), dimension(:)    , allocatable :: lper
   !    real(WP), dimension(:,:)  , allocatable :: lpos
   !    real(WP), dimension(:,:)  , allocatable :: lvel
   !    real(WP), dimension(:,:,:), allocatable :: lmoi
   !    real(WP), dimension(:)    , allocatable :: lrem
   !    real(WP), dimension(:)    , allocatable :: lSR
   !    real(WP), dimension(:)    , allocatable :: xmin,xmax,ymin,ymax,zmin,zmax
   !    integer :: n,m,ierr,i,j,k,l,ii,jj,kk,iunit,totalnewp,np_start,np_old,count,ip,rank
   !    real(WP) :: x,y,z,x0,y0,z0,lmax,lmid,lmin
   !    character(len=str_medium) :: filename
   !    integer, dimension(:), allocatable ::  plist,dispels
   !    real(WP), dimension(:,:), allocatable :: pinfo,pinfo_
   !    real(WP) :: Vt,Vl,Vd,minor_radius,diam,Vrim,Lrim
   !    real(WP) :: Oh,Trp,Lrp,Tsr,SR_tmp
   !    real(WP), dimension(1:3) :: tangent
   !    real(WP), dimension(:,:,:,:), allocatable :: SR
   !    integer :: nmain,nsat
   !    integer :: nmax
   !    real(WP), dimension(:,:,:), allocatable :: thickness
   !    integer,  dimension(:,:,:), allocatable :: struct_type
   !    ! Moment of inertia calculation using lapack
   !    real(WP), dimension(:), allocatable, save :: work !< Saved!
   !    integer, save :: lwork                            !< Saved!
   !    real(WP), dimension(1) :: lwork_query
   !    real(WP), dimension(3) :: d
   !    real(WP), dimension(3,3) :: A
   !    integer :: info

   !    ! Query optimal work array size
   !    if (.not.allocated(work)) then
   !    call dsyev('V','U',3,A,3,d,lwork_query,-1,info)
   !    lwork=int(lwork_query(1)); allocate(work(lwork))
   !    end if

   !    ! Get thickness and local struct_type for global information calculation
   !    allocate(thickness  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));thickness=0.0_WP
   !    allocate(struct_type(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));struct_type=0
   !    call get_liginfo()

   !    this%thickness=thickness
   !    this%struct_type=struct_type*1.0_WP
   !    ! Start by performing a CCL based on ligament criteria
   !    call this%ccl_lig%build(make_label,same_label)

   !    if (this%ccl_lig%nstruct.ge.1) then

   !    ! Allocate ligament stats arrays
   !    allocate(lvol(1:this%ccl_lig%nstruct        )); lvol=0.0_WP
   !    allocate(lthc(1:this%ccl_lig%nstruct        )); lthc=HUGE(x)
   !    allocate(llen(1:this%ccl_lig%nstruct        )); llen=0.0_WP
   !    allocate(lnum(1:this%ccl_lig%nstruct        )); lnum=0.0_WP
   !    allocate(lper(1:this%ccl_lig%nstruct        )); lper=0.0_WP
   !    allocate(lpos(1:this%ccl_lig%nstruct,1:3    )); lpos=0.0_WP
   !    allocate(lvel(1:this%ccl_lig%nstruct,1:3    )); lvel=0.0_WP
   !    allocate(lmoi(1:this%ccl_lig%nstruct,1:3,1:3)); lmoi=0.0_WP
   !    allocate(lrem(1:this%ccl_lig%nstruct        )); lrem=0.0_WP
   !    allocate(lSR(1:this%ccl_lig%nstruct        )); lSR=-HUGE(x)
   !    allocate(xmin(1:this%ccl_lig%nstruct),xmax(1:this%ccl_lig%nstruct)); xmin=HUGE(x);xmax=-HUGE(x)
   !    allocate(ymin(1:this%ccl_lig%nstruct),ymax(1:this%ccl_lig%nstruct)); ymin=HUGE(x);ymax=-HUGE(x)
   !    allocate(zmin(1:this%ccl_lig%nstruct),zmax(1:this%ccl_lig%nstruct)); zmin=HUGE(x);zmax=-HUGE(x)
   !    allocate(SR(1:6,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_));SR=0.0_WP
   !    call this%fs%get_strainrate(SR)
   !    ! First pass to accumulate volume, position, min thickness and ligament percentage
   !    do n=1,this%ccl_lig%nstruct
   !       ! Loop over cells in structure
   !       lnum(n)=lnum(n)+1.0_WP*this%ccl_lig%struct(n)%n_
   !       do m=1,this%ccl_lig%struct(n)%n_
   !           ! Get cell indices
   !           i=this%ccl_lig%struct(n)%map(1,m)
   !           j=this%ccl_lig%struct(n)%map(2,m)
   !           k=this%ccl_lig%struct(n)%map(3,m)
   !           ! Get cell position, accounting for periodicity
   !           x=this%vf%cfg%xm(i)-this%ccl_lig%struct(n)%per(1)*this%vf%cfg%xL
   !           y=this%vf%cfg%ym(j)-this%ccl_lig%struct(n)%per(2)*this%vf%cfg%yL
   !           z=this%vf%cfg%zm(k)-this%ccl_lig%struct(n)%per(3)*this%vf%cfg%zL
   !           ! Accumulate volume and position. Get min thickness and ligament percentage
   !           lvol(n  )=lvol(n  )+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
   !           lpos(n,:)=lpos(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[x,y,z]
   !           lvel(n,:)=lvel(n,:)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*[this%Ui(i,j,k),this%Vi(i,j,k),this%Wi(i,j,k)]
   !           lthc(n)=min(lthc(n),thickness(i,j,k))
   !           if (struct_type(i,j,k).eq.1) lper(n)=lper(n)+1.0_WP
   !           ! Check if ligament touches auto-transfer layer
   !           if (i.ge.this%vf%cfg%imax-this%nlayer.or.&
   !           &   j.le.this%vf%cfg%jmin+this%nlayer.or.&
   !           &   j.ge.this%vf%cfg%jmax-this%nlayer.or.&
   !           &   k.le.this%vf%cfg%kmin+this%nlayer.or.&
   !           &   k.ge.this%vf%cfg%kmax-this%nlayer) lrem(n)=1.0_WP
   !           ! Get the structures's locally largest and smallest x,y,z locations
   !           do l=1,2
   !             if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
   !                d = calculateCentroid(this%vf%interface_polygon(l,i,j,k))
   !                xmin(n)=min(xmin(n),d(1)); xmax(n)=max(xmax(n),d(1))
   !                ymin(n)=min(ymin(n),d(2)); ymax(n)=max(ymax(n),d(2))
   !                zmin(n)=min(zmin(n),d(3)); zmax(n)=max(zmax(n),d(3))
   !             end if
   !          end do
   !       end do
   !    end do
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lvol,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lpos,3*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lvel,3*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lrem,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lthc,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lnum,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lper,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,xmin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,ymin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,zmin,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MIN,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,xmax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,ymax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,zmax,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)
        
   !    ! Second pass to accumulate moment of inertia
   !    do n=1,this%ccl_lig%nstruct
   !       ! Get ligament barycenter
   !       x0=lpos(n,1)/lvol(n)
   !       y0=lpos(n,2)/lvol(n)
   !       z0=lpos(n,3)/lvol(n)
   !       ! Loop over cells in structure
   !       do m=1,this%ccl_lig%struct(n)%n_
   !           ! Get cell indices
   !           i=this%ccl_lig%struct(n)%map(1,m)
   !           j=this%ccl_lig%struct(n)%map(2,m)
   !           k=this%ccl_lig%struct(n)%map(3,m)
   !           ! Get cell position relative to drop barycenter, accounting for periodicity
   !           x=this%vf%cfg%xm(i)-this%ccl_lig%struct(n)%per(1)*this%vf%cfg%xL-x0
   !           y=this%vf%cfg%ym(j)-this%ccl_lig%struct(n)%per(2)*this%vf%cfg%yL-y0
   !           z=this%vf%cfg%zm(k)-this%ccl_lig%struct(n)%per(3)*this%vf%cfg%zL-z0
   !           ! Accumulate moment of inertia
   !           lmoi(n,1,1)=lmoi(n,1,1)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y**2+z**2)
   !           lmoi(n,2,2)=lmoi(n,2,2)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(z**2+x**2)
   !           lmoi(n,3,3)=lmoi(n,3,3)+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x**2+y**2)
   !           lmoi(n,1,2)=lmoi(n,1,2)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*y)
   !           lmoi(n,1,3)=lmoi(n,1,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(x*z)
   !           lmoi(n,2,3)=lmoi(n,2,3)-this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)*(y*z)
   !       end do
   !    end do
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lmoi,9*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)

   !    ! Third pass to generalize ligament stats
   !    do n=1,this%ccl_lig%nstruct
   !       ! Get ligament, accounting for periodicity
   !       lpos(n,:)=lpos(n,:)/lvol(n)
   !       if (this%vf%cfg%xper.and.lpos(n,1).lt.this%vf%cfg%x(this%vf%cfg%imin)) lpos(n,1)=lpos(n,1)+this%vf%cfg%xL
   !       if (this%vf%cfg%yper.and.lpos(n,2).lt.this%vf%cfg%y(this%vf%cfg%jmin)) lpos(n,2)=lpos(n,2)+this%vf%cfg%yL
   !       if (this%vf%cfg%zper.and.lpos(n,3).lt.this%vf%cfg%z(this%vf%cfg%kmin)) lpos(n,3)=lpos(n,3)+this%vf%cfg%zL
   !       ! Get drop velocity
   !       lvel(n,:)=lvel(n,:)/lvol(n)
   !       ! Calculate the percentage of ligament structure type
   !       lper(n)=lper(n)/lnum(n)
   !       ! Calculate maximum length of the structure
   !       A=lmoi(n,:,:)
   !       call dsyev('V','U',3,A,3,d,work,lwork,info) !< On exit, A contains eigenvectors and d contains eigenvalues in ascending order
   !       d=max(0.0_WP,d)    
   !       ! Replace with corrected eigenvectors for future ligament droplet placement
   !       lmoi(n,:,:)=A
   !       ! Get characteristic lengths of drop
   !       lmax=sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/lvol(n))
   !       ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n," lmax=",lmax," d(2)=",d(2)," d(1)=",d(1)," d(3)=",d(3)," lvol(n)=",lvol(n)
   !       lmid=sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/lvol(n))
   !       lmin=sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/lvol(n))
   !       if (lmin.eq.0.0_WP) lmin=lmid ! Handle 2D case
   !       ! Use max of bounding box and MoI-derived lengths as length
   !       !hypot(hypot(xmax(n)-xmin(n),ymax(n)-ymin(n))**2,zmax(n)-zmin(n))
   !       llen(n) = max(sqrt((xmax(n)-xmin(n))**2+(ymax(n)-ymin(n))**2+(zmax(n)-zmin(n))**2),lmax)
   !       ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " llen(n)=",llen(n), " lmax=",lmax
   !       ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " xmax(n)=",xmax(n), " xmin(n)=",xmin(n), " ymax(n)=",ymax(n), " ymin(n)=",ymin(n), " zmax(n)=",zmax(n), " zmin(n)=",zmin(n)

   !       ! With the tangent direction of the ligament, we can evaluate the strain rate of each cell of the ligament
   !       tangent = lmoi(n,:,1)
   !       do m=1,this%ccl_lig%struct(n)%n_
   !          ! Get cell indices
   !          i=this%ccl_lig%struct(n)%map(1,m)
   !          j=this%ccl_lig%struct(n)%map(2,m)
   !          k=this%ccl_lig%struct(n)%map(3,m)
   !          SR_tmp =SR(1,i,j,k)*tangent(1)**2        +SR(2,i,j,k)*tangent(2)**2        +SR(3,i,j,k)*tangent(3)**2 + &
   !        & 2.0_WP*(SR(4,i,j,k)*tangent(1)*tangent(2)+SR(5,i,j,k)*tangent(2)*tangent(3)+SR(6,i,j,k)*tangent(1)*tangent(3))
   !          lSR(n) = max(lSR(n),abs(SR_tmp))
   !       end do
   !    end do
   !    ! Find the maximum tangential strain rate of each ligament
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,lSR,1*this%ccl_lig%nstruct,MPI_REAL_WP,MPI_MAX,this%vf%cfg%comm,ierr)

   !    ! Find the liquid core
   !    nmax=maxloc(lvol,dim=1)

   !    ! Zero out monitoring variables
   !    this%vof_tf_lig=0.0_WP
   !    this%np_lig=0
   !    ! Record initial droplets in each processor for future outputing purpose
   !    np_start=this%lp%np_
   !    ! Perform transfer
   !    do n=1,this%ccl_lig%nstruct
   !       ! Cycle if struct is core
   !       if (n.eq.nmax) cycle
   !       ! Assume a cylinder ligament
   !       Lrim=llen(n) !< this is 0 sometimes
   !       Vrim=lvol(n)
   !       minor_radius=sqrt(Vrim/pi/Lrim)  
   !       ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " minor_radius=",minor_radius, " Lrim=",Lrim, " Vrim=",Vrim
   !       ! Drop size method from Kim & Moin (2020)
   !       nmain=floor(this%dw*Lrim/(twoPi*minor_radius))
   !       ! Calculate breakup time scale based on inviscid RP instability analysis
   !       Trp=2.91258_WP*sqrt(this%fs%rho_l*minor_radius**3/this%fs%sigma)
   !       ! Calcuate time scale based on maximum local strainrate
   !       Tsr=1.0_WP/lSR(n)
   !       ! Only breakup if minimum thickness is reached, sufficient volume of the ligament, enough of local ligament-like structures,
   !       ! local time scale asscoiated with strain rate is on par or bigger than the RP time scale, and its length is longer than the inviscid most unstable wavelength
   !       if ((lthc(n).le.this%lmin*this%cfg%min_meshsize).and.(lvol(n).ge.this%cfg%min_meshsize**3).and.(lper(n).ge.this%lper).and.(Trp.le.Tsr).and.(nmain.ge.1)) then
   !       ! else if(lrem(n).gt.0.0_WP) then
   !          ! if (this%vf%cfg%amRoot) print *, "lig is in buffer with nmain=", nmain, " and id=", n
   !       else
   !          cycle
   !       end if

   !       if (llen(n).le.0.0_WP) cycle
      
   !       ! if (this%vf%cfg%amRoot) print *, "This is the min_thickness", lthc(n), ",lig percentage:", lper(n),"max length:",llen(n),&
   !       ! & "how many cells",lnum(n), "vol:",lvol(n),"nmain", nmain, "Trp:", Trp, "Tsr:", Tsr, "Trp/Tsr", Trp/Tsr,"and id:", n
      
   !       nsat=nmain+1
   !       diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+this%size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)
         
   !       ! if (this%vf%cfg%amRoot) print *, "This is the diameter: ", diam,"and id: ", n

   !       ! Only the main processor is in charge of creating droplets
   !       if (this%cfg%amRoot) then
   !          Lrp = twoPi*minor_radius/this%dw
   !          ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " Lrp=",Lrp, " twoPi=",twoPi, " minor_radius=",minor_radius, " this%dw=",this%dw
   !          do l=1,nsat+nmain
   !             ! Increment particle counter
   !             this%lp%np_=this%lp%np_+1
   !             ! Make room for new drop
   !             call this%lp%resize(this%lp%np_)
   !             ! Add the drop
   !             this%lp%p(this%lp%np_)%id  =int(3,8)                                                                               
   !             if (mod(l,2).eq.1) then
   !                this%lp%p(this%lp%np_)%d=diam*this%size_ratio                                                                                    
   !             else
   !                this%lp%p(this%lp%np_)%d=diam                                                                                    
   !             end if
   !             ! if (llen(n).eq.0.0_WP) then
   !             !    this%lp%p(this%lp%np_)%pos=lpos(n,:)
   !             ! else
   !             !    this%lp%p(this%lp%np_)%pos=lpos(n,:)+0.5_WP*Lrp*(l-(nmain+1))*lmoi(n,:,1)
   !             ! end if
   !             this%lp%p(this%lp%np_)%pos=lpos(n,:)+0.5_WP*Lrp*(l-(nmain+1))*lmoi(n,:,1)
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " lpos(n,1)=",lpos(n,1), " lpos(n,2)=",lpos(n,2), " lpos(n,3)=",lpos(n,3)
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " Lrp=",Lrp
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " lmoi(n,1,1)=",lmoi(n,1,1), " lmoi(n,2,1)=",lmoi(n,2,1), " lmoi(n,3,1)=",lmoi(n,3,1)
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " pos(1)=",this%lp%p(this%lp%np_)%pos(1), " pos(2)=",this%lp%p(this%lp%np_)%pos(2), " pos(3)=",this%lp%p(this%lp%np_)%pos(3)   
   !             this%lp%p(this%lp%np_)%vel =lvel(n,:)
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " lvel(n,1)=",lvel(n,1), " lvel(n,2)=",lvel(n,2), " lvel(n,3)=",lvel(n,3)
   !             this%lp%p(this%lp%np_)%ind =this%cfg%get_ijk_global(this%lp%p(this%lp%np_)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])
   !             ! if (this%vf%cfg%amRoot) print *, "This paritcle id id=", n, " ind(1)=",this%lp%p(this%lp%np_)%ind(1), " ind(2)=",this%lp%p(this%lp%np_)%ind(2), " ind(3)=",this%lp%p(this%lp%np_)%ind(3)   
   !             this%lp%p(this%lp%np_)%flag=0                                                                                        
   !             this%lp%p(this%lp%np_)%dt  =0.0_WP                                                                                  
   !             this%lp%p(this%lp%np_)%Acol=0.0_WP                                                                                  
   !             this%lp%p(this%lp%np_)%Tcol=0.0_WP
   !          end do
   !          ! Increment monitoring variables
   !          this%lp%np_new=this%lp%np_new+nmain+nsat
   !          this%np_lig=this%np_lig+nmain+nsat
   !          this%vof_tf_lig=this%vof_tf_lig+lvol(n)
   !          this%lp%vp_new=this%lp%vp_new+lvol(n)
   !       end if

   !       ! empty out the VF
   !       do m=1,this%ccl_lig%struct(n)%n_
   !          i=this%ccl_lig%struct(n)%map(1,m); j=this%ccl_lig%struct(n)%map(2,m); k=this%ccl_lig%struct(n)%map(3,m)
   !          this%vf%VF(i,j,k)=0.0_WP
   !       end do    

   !    end do

   !    ! Synchronize VF fields
   !    call this%vf%cfg%sync(this%vf%VF)
   !    call this%vf%clean_irl_and_band()

   !    ! Synchronize particles
   !    call this%lp%sync()
      
   !    ! Integrate monitoring variables 
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_tf_lig,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(MPI_IN_PLACE,this%np_lig    ,1,MPI_INTEGER,MPI_SUM,this%vf%cfg%comm,ierr)

   !    deallocate(thickness,struct_type)
   !    deallocate(lvol,lthc,llen,lnum,lper,lpos,lvel,lmoi,lrem,lSR,xmin,ymin,zmin,SR)

   !    end if

   ! contains

   !    ! Calculate thickness and struct_type based on moment of inertia
   !    subroutine get_liginfo()
   !       implicit none 
   !       real(WP) :: tmpvol,tmparea
   !       real(WP), dimension(1:3) :: tmpxvol, tmpL
   !       integer :: nneigh_moi, nneigh_thickness
   !       nneigh_moi=2; nneigh_thickness=3
   !       do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
   !          do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
   !             do i=this%vf%cfg%imin_,this%vf%cfg%imax_
   !                ! calculate thickness
   !                tmpvol=0.0_WP; tmparea=0.0_WP
   !                do kk = k-nneigh_thickness,k+nneigh_thickness
   !                   do jj = j-nneigh_thickness,j+nneigh_thickness
   !                      do ii = i-nneigh_thickness,i+nneigh_thickness
   !                         tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
   !                         tmparea = tmparea + this%vf%SD(ii,jj,kk)*this%cfg%vol(i,j,k)
   !                      end do
   !                   end do
   !                end do
   !                ! Calculate thickness
   !                if (this%vf%VF(i,j,k).le.VFlo) then
   !                   thickness(i,j,k) = 0.0_WP
   !                else if (tmparea .gt. 0.0_WP) then    
   !                   thickness(i,j,k) = 2.0_WP*tmpvol/(tmparea+tiny(1.0_WP))
   !                else
   !                   thickness(i,j,k) = 3.5_WP*this%cfg%min_meshsize
   !                end if

   !                ! Calculate moi
   !                tmpvol=0.0_WP; tmpxvol=0.0_WP; A=0.0_WP
   !                ! First pass to accumulate volume, surface area, and position
   !                do kk = k-nneigh_moi,k+nneigh_moi
   !                   do jj = j-nneigh_moi,j+nneigh_moi
   !                      do ii = i-nneigh_moi,i+nneigh_moi
   !                         tmpvol = tmpvol + this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
   !                         tmpxvol = tmpxvol + this%vf%Lbary(:,ii,jj,kk)*this%vf%VF(ii,jj,kk)*this%cfg%vol(i,j,k)
   !                      end do
   !                   end do
   !                end do
   !                ! Second pass to accumulate moment of inertia
   !                tmpxvol = tmpxvol/tmpvol
   !                do kk = k-nneigh_moi,k+nneigh_moi
   !                   do jj = j-nneigh_moi,j+nneigh_moi
   !                      do ii = i-nneigh_moi,i+nneigh_moi
   !                         ! Location of film node
   !                         tmpL = this%vf%Lbary(:,ii,jj,kk) - tmpxvol
   !                         A(1,1)=A(1,1)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(2)**2+tmpL(3)**2)
   !                         A(2,2)=A(2,2)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(3)**2)
   !                         A(3,3)=A(3,3)+this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*(tmpL(1)**2+tmpL(2)**2)
   !                         A(1,2)=A(1,2)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(2)
   !                         A(1,3)=A(1,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(1)*tmpL(3)
   !                         A(2,3)=A(2,3)-this%vf%cfg%vol(ii,jj,kk)*this%vf%VF(ii,jj,kk)*tmpL(2)*tmpL(3)   
   !                      end do
   !                   end do
   !                end do
   !                ! Calculate local struct type
   !                call dsyev('V','U',3,A,3,d,work,lwork,info)
   !                d=max(0.0_WP,d)
   !                if (d(3).gt.(this%lstratio*d(1))) struct_type(i,j,k)=struct_type(i,j,k)+1
   !                if (d(3).gt.(this%lstratio*d(2))) struct_type(i,j,k)=struct_type(i,j,k)+1
   !             end do 
   !          end do 
   !       end do
   !       call this%vf%cfg%sync(thickness)
   !       call this%vf%cfg%sync(struct_type)
   !    end subroutine

   !    !> Function that identifies cells that need a label
   !    logical function make_label(i,j,k)
   !       implicit none
   !       integer, intent(in) :: i,j,k
   !       if ((this%vf%VF(i,j,k).gt.VFlo).and.(thickness(i,j,k).lt.this%lmake*this%cfg%min_meshsize))then
   !          make_label=.true.
   !       else
   !        make_label=.false.
   !       end if
   !    end function make_label

   !    !> Function that identifies if cell pairs have same label
   !    logical function same_label(i1,j1,k1,i2,j2,k2)
   !       implicit none
   !       integer, intent(in) :: i1,j1,k1,i2,j2,k2
   !       same_label=.true.
   !    end function same_label

   ! end subroutine transfer_ligs
   
   
end module simulation