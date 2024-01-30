!> Various definitions and tools for running an NGA2 simulation
module simulation
   use inflow_class,   only: inflow
   implicit none
   private
   
   !> pipe simulation
   type(inflow) :: pipe
   
   
   public :: simulation_init,simulation_run,simulation_final

contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize pipe simulation
      call pipe%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.pipe%time%done())

         call pipe%step()

      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Finalize pipe simulation
      call pipe%final()
      
   end subroutine simulation_final
   

end module simulation
