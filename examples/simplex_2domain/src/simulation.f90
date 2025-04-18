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


      ! Simplex drives overall time integration
      do while (.not.spx%time%done())
         
         call spx%step()

      end do
      
      ! ! Simplex drives overall time integration
      ! do while (.not.spx%time%done())
         
      !    ! Advance simplex simulation
      !    call spx%step()

      !    ! Handle coupling velocity between simplex and atomization
      !    coupling_velocity_s2a: block
      !       use tpns_class, only: bcond
      !       integer :: n,i,j,k
      !       type(bcond), pointer :: mybc
      !       ! Exchange data using cpl12x/y/z couplers
      !       call xcpl_s2a%push(spx%fs%U);   call xcpl_s2a%transfer();  call xcpl_s2a%pull(atomization%resU)
      !       call ycpl_s2a%push(spx%fs%V);   call ycpl_s2a%transfer();  call ycpl_s2a%pull(atomization%resV)
      !       call zcpl_s2a%push(spx%fs%W);   call zcpl_s2a%transfer();  call zcpl_s2a%pull(atomization%resW)
      !       call atomization%fs%get_bcond('inlets',mybc)
      !       do n=1,mybc%itr%no_
      !          i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
      !          atomization%fs%U(i  ,j,k)=atomization%resU(i  ,j,k)*sum(atomization%fs%itpr_x(:,i  ,j,k)*atomization%cfg%VF(i-1:i,    j,    k))
      !          atomization%fs%V(i-1,j,k)=atomization%resV(i-1,j,k)*sum(atomization%fs%itpr_y(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j-1:j,    k))
      !          atomization%fs%W(i-1,j,k)=atomization%resW(i-1,j,k)*sum(atomization%fs%itpr_z(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j    ,k-1:k))
      !       end do
      !    end block coupling_velocity_s2a


      !    ! Handle coupling VOF between simplex and atomization
      !    coupling_vof_s2a: block
      !       integer :: i,j,k
      !       ! Exchange data using cell center coupler
      !       atomization%tempVF=0.0_WP
      !       call vfcpl_s2a%push(spx%vf%VF); call vfcpl_s2a%transfer(); call vfcpl_s2a%pull(atomization%tempVF)
      !       ! Exchange VOF based upon x/y/z position
      !       do k=atomization%vf%cfg%kmino_,atomization%vf%cfg%kmaxo_
      !          do j=atomization%vf%cfg%jmino_,atomization%vf%cfg%jmaxo_
      !             do i=atomization%vf%cfg%imino_,atomization%vf%cfg%imaxo_
      !                if (atomization%vf%cfg%xm(i).ge.atomization%vfcouple_xmin.and.atomization%vf%cfg%xm(i).le.atomization%vfcouple_xmax) then
      !                   atomization%vf%VF(i,j,k)=atomization%tempVF(i,j,k)
      !                end if
      !             end do
      !          end do
      !       end do
      !       !> sync arrays
      !       call atomization%cfg%sync(atomization%vf%VF)
      !       !> Reconstruct interface to get updated VF and moments
      !       call atomization%vf%build_interface()
      !       call atomization%vf%reset_volume_moments()
      !    end block coupling_vof_s2a
      
      !    ! Advance atomization simulation until it's caught up
      !    do while (atomization%time%t.le.spx%time%t)
      !       call atomization%step()
      !    end do

      ! end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize simplex simulation
      call spx%final()

      ! Finalize atomization simulation
      call atomization%final()
      
   end subroutine simulation_final
   
   
end module simulation