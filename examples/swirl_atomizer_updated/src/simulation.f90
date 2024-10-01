!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
   use geometry,          only: cfg1,group1,isInGrp1
   use geometry,          only: cfg2,group2,isInGrp2
   use block1_class,      only: block1
   use block2_class,      only: block2
   use coupler_class,     only: coupler
   use timetracker_class, only: timetracker
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final

   !> Block 1 and 2 objects
   type(block1) :: b1
   type(block2) :: b2

   !> Couplers between blocks
   type(coupler) :: cpl12x,cpl12y,cpl12z

   !> Time when annular pipe is turbulent
   real(WP) :: transition_time

   !> Storage for coupled fields
   real(WP), dimension(:,:,:), allocatable :: U1on2,V1on2,W1on2

   
contains
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_read
      implicit none
      
      ! Initialize blocks from mpi groups
      if (isInGrp1) then 
         b1%cfg=>cfg1
         call b1%init()
         ! Time for annular pipe to become turbulent 
         call param_read('Transistion',transition_time)
      end if
      if (isInGrp2) then 
         b2%cfg=>cfg2
         call b2%init()
      end if

      ! Initialize the couplers
      coupler_prep: block
         use parallel, only: group
         ! Both groups prepare the coupler (couple pipe to swirl)
         if (isInGrp1.or.isInGrp2) then
            ! Block 1 to block 2
            cpl12x=coupler(src_grp=group1,dst_grp=group2,name='pipe_to_swirl_x') 
            cpl12y=coupler(src_grp=group1,dst_grp=group2,name='pipe_to_swirl_y') 
            cpl12z=coupler(src_grp=group1,dst_grp=group2,name='pipe_to_swirl_z') 
            ! Set the grids
            if (isInGrp1) then
               call cpl12x%set_src(cfg1,'x')  
               call cpl12y%set_src(cfg1,'y')  
               call cpl12z%set_src(cfg1,'z')
            end if
            if (isInGrp2) then 
               call cpl12x%set_dst(cfg2,'x')
               call cpl12y%set_dst(cfg2,'y')
               call cpl12z%set_dst(cfg2,'z')
            end if
            ! Initialize the metrics
            call cpl12x%initialize()
            call cpl12y%initialize()
            call cpl12z%initialize()
            ! Allocate arrays for velocity coupling
            if (isInGrp2) then
               allocate(U1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U1on2=0.0_WP
               allocate(V1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); V1on2=0.0_WP
               allocate(W1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); W1on2=0.0_WP
            end if
         end if
      end block coupler_prep


   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation 
   subroutine simulation_run
      implicit none

      ! Advance block 1 until annular pipe is turbulent
      if (isInGrp1) then 
         
         do while (b1%time%t.le.transition_time)
            
            call b1%step()
         
         end do
         
         ! Reset block 1 time to 0
         b1%time%t=0.0_WP

      end if
      
      ! Perform time integration - block 2 is the main driver here
      if (isInGrp1.or.isInGrp2) then
         
         do while (.not.b2%time%done()) 
            
            ! Exchange velocity field using cpl12x/y/z couplers
            coupling_step: block
               ! Zero the fields
               U1on2=0.0_WP; V1on2=0.0_WP; W1on2=0.0_WP 
               ! Push the data from B1
               if (isInGrp1) then
                  call cpl12x%push(b1%fs%U) 
                  call cpl12y%push(b1%fs%V) 
                  call cpl12z%push(b1%fs%W) 
               end if
               ! Transfer the data
               call cpl12x%transfer() 
               call cpl12y%transfer() 
               call cpl12z%transfer()
               ! Pull the data to B2
               if (isInGrp2) then
                  call cpl12x%pull(U1on2)
                  call cpl12y%pull(V1on2)
                  call cpl12z%pull(W1on2)
               end if
            end block coupling_step
            
            ! Advance block 2
            if (isInGrp2) call b2%step(U1on2,V1on2,W1on2)
            
            ! Broadcast current B2 time and done condition
            b2_broadcast: block
               use mpi_f08,  only: MPI_BCAST,MPI_LOGICAL
               use parallel, only: comm,nproc,MPI_REAL_WP
               integer :: ierr
               ! Broadcast B2 time to all processors
               call MPI_BCAST(b2%time%t,1,MPI_REAL_WP,nproc-1,comm,ierr)
               ! Broadcast B2 done condition
               call MPI_BCAST(b2%time%done(),1,MPI_LOGICAL,nproc-1,comm,ierr)
            end block b2_broadcast

            ! Advance block 1 until we've caught up
            if (isInGrp1) then
               
               do while (b1%time%t.lt.b2%time%t)

                  ! Advance block 1
                  call b1%step()

               end do

            end if
   
         end do

      end if
      
   end subroutine simulation_run
      
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays for both blocks
      if (isInGrp1) call b1%final()
      if (isInGrp2) call b2%final()
      
   end subroutine simulation_final
   
   
end module simulation
