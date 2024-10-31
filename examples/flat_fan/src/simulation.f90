!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use flat_fan_class, only: flat_fan
   implicit none
   private
   
   !> flat fan simulation
   type(flat_fan) :: ff
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
      ! Initialize flat fan simulation
      call ff%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! flat fan drives overall time integration
      do while (.not.ff%time%done())
         ! Advance flat fan simulation
         call ff%step()
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize flat fan simulation
      call ff%final()
      
   end subroutine simulation_final
   
   
end module simulation