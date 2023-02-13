!> Various definitions and tools for running an NGA2 simulation
module simulation
   use string,            only: str_medium
   use precision,         only: WP
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
      use geometry, only: cfg1,cfg2
      use param,    only: param_read
      implicit none
      
      ! Initialize both blocks
      b1%cfg=>cfg1; call b1%init()
      b2%cfg=>cfg2; call b2%init()

      ! Initialize the couplers
      coupler_prep: block
         use parallel, only: group
         ! Block 1 to block 2
         cpl12x=coupler(src_grp=group,dst_grp=group,name='pipe_to_swirl_x'); call cpl12x%set_src(cfg1,'x'); call cpl12x%set_dst(cfg2,'x'); call cpl12x%initialize()
         cpl12y=coupler(src_grp=group,dst_grp=group,name='pipe_to_swirl_y'); call cpl12y%set_src(cfg1,'y'); call cpl12y%set_dst(cfg2,'y'); call cpl12y%initialize()
         cpl12z=coupler(src_grp=group,dst_grp=group,name='pipe_to_swirl_z'); call cpl12z%set_src(cfg1,'z'); call cpl12z%set_dst(cfg2,'z'); call cpl12z%initialize()
         allocate(U1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); U1on2=0.0_WP
         allocate(V1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); V1on2=0.0_WP
         allocate(W1on2(cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); W1on2=0.0_WP
      end block coupler_prep

      ! Setup nudging region in block 2
      b2%nudge_trans=20.0_WP*b2%cfg%min_meshsize
      b2%nudge_xmin =b1%cfg%x(3*b1%cfg%nx/4)
      b2%nudge_xmax =b1%cfg%x(3*b1%cfg%nx/4)
      b2%nudge_ymin =b1%cfg%y(b2%cfg%jmin)
      b2%nudge_ymax =b1%cfg%y(b2%cfg%jmax+1)
      b2%nudge_zmin =b1%cfg%z(b2%cfg%kmin)
      b2%nudge_zmax =b1%cfg%z(b2%cfg%kmax+1)

      ! Time for annular pipe to become turbulent 
      call param_read('Transistion',transition_time)
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation 
   subroutine simulation_run
      implicit none
      integer :: i,j,k
      
      ! Perform time integration - block 2 is the main driver here
      do while (.not.b2%time%done())
         
         if(b1%time.le.transistion_time) then 

            ! Advance block 1 until annular pipe is turbulent
            call b1%step()

         else
            
            ! Advance block 2
            call b2%step()

            ! Advance block 1 until we've caught up
            do while (b1%time%t.lt.b2%time%t)

               ! Advance block 1
               call b1%step()

               ! Exchange data using cpl12x/y/z couplers and the most recent velocity
               U1on2=0.0_WP; call cpl12x%push(b1%fs%U); call cpl12x%transfer(); call cpl12x%pull(U1on2)
               V1on2=0.0_WP; call cpl12y%push(b1%fs%V); call cpl12y%transfer(); call cpl12y%pull(V1on2)
               W1on2=0.0_WP; call cpl12z%push(b1%fs%W); call cpl12z%transfer(); call cpl12z%pull(W1on2)

            end do

         end if 
 
      end do
      
      
   end subroutine simulation_run
      
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Deallocate work arrays for both blocks
      call b1%final()
      call b2%final()
      
   end subroutine simulation_final
   
   
end module simulation
