!> Various definitions and tools for running an NGA2 simulation
module simulation
   use mpi_f08,         only: MPI_Group
   use precision,       only: WP
   use simplex_class,   only: simplex
   use coupler_class,   only: coupler
   use inputfile_class, only: inputfile
   implicit none
   private
   
   !> Simplex simulation
   type(simplex) :: spx

   !> Input files to read in patitions
   type(inputfile) :: input_spx

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
      implicit none
      
      ! Initialize simplex simulation
      call spx%init()
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none


      ! Simplex drives overall time integration
      do while (.not.spx%time%done())
         
         call spx%step()

      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize simplex simulation
      call spx%final()
      
   end subroutine simulation_final
   
   
end module simulation