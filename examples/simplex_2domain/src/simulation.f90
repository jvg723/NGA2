!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,     only: WP
   use simplex_class, only: simplex
   use atom_class,    only: atom
   use coupler_class, only: coupler
   implicit none
   private
   
   !> Simplex simulation
   type(simplex) :: spx

   !> Atomization simulation
   type(atom) :: atomization

   !> Couplers from simplex to atomization
   type(coupler) :: xcpl_s2a,ycpl_s2a,zcpl_s2a !> Velocity
   type(coupler) :: vfcpl_s2a

   !> Storage for passing VOF
   real(WP), dimension(:,:,:), allocatable :: tempVF
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      implicit none
      
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
         ! allocate storage for temp VOF
         allocate(tempVF(atomization%cfg%imino_:atomization%cfg%imaxo_,atomization%cfg%jmino_:atomization%cfg%jmaxo_,atomization%cfg%kmino_:atomization%cfg%kmaxo_)); tempVF=0.0_WP 
      end block create_coupler_s2a

      ! Create region to couple VOF between domains
      create_coupler_region: block
         atomization%vfcouple_xmin=-1.00_WP*atomization%xshift
         atomization%vfcouple_xmax=spx%cfg%xm(spx%cfg%imax-spx%nlayer)
         atomization%vfcouple_ymin=spx%cfg%ym(spx%cfg%jmin)
         atomization%vfcouple_ymax=spx%cfg%ym(spx%cfg%jmax)
         atomization%vfcouple_zmin=spx%cfg%zm(spx%cfg%kmin)
         atomization%vfcouple_zmax=spx%cfg%zm(spx%cfg%kmax)
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
      do while (.not.atomization%time%done())
         
         ! Advance simplex simulation until it's caught up
         do while (spx%time%t.le.atomization%time%t)
            call spx%step()
         end do

         ! Handle coupling velocity between simplex and atomization
         coupling_velocity_s2a: block
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cpl12x/y/z couplers
            call xcpl_s2a%push(spx%fs%U);   call xcpl_s2a%transfer();  call xcpl_s2a%pull(atomization%resU)
            call ycpl_s2a%push(spx%fs%V);   call ycpl_s2a%transfer();  call ycpl_s2a%pull(atomization%resV)
            call zcpl_s2a%push(spx%fs%W);   call zcpl_s2a%transfer();  call zcpl_s2a%pull(atomization%resW)
            !>Pass VF field to atomization domain
            call vfcpl_s2a%push(spx%vf%VF); call vfcpl_s2a%transfer(); call vfcpl_s2a%pull(atomization%vf%VF)
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
            use tpns_class, only: bcond
            integer :: n,i,j,k
            type(bcond), pointer :: mybc
            ! Exchange data using cell center coupler
            tempVF=0.0_WP
            call vfcpl_s2a%push(spx%vf%VF); call vfcpl_s2a%transfer(); call vfcpl_s2a%pull(tempVF)
            ! call atomization%fs%get_bcond('inlets',mybc)
            ! do n=1,mybc%itr%no_
            !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            !    atomization%fs%U(i  ,j,k)=atomization%resU(i  ,j,k)*sum(atomization%fs%itpr_x(:,i  ,j,k)*atomization%cfg%VF(i-1:i,    j,    k))
            !    atomization%fs%V(i-1,j,k)=atomization%resV(i-1,j,k)*sum(atomization%fs%itpr_y(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j-1:j,    k))
            !    atomization%fs%W(i-1,j,k)=atomization%resW(i-1,j,k)*sum(atomization%fs%itpr_z(:,i-1,j,k)*atomization%cfg%VF(i-1  ,j    ,k-1:k))
            ! end do
         end block coupling_vof_s2a
      
         ! Advance atomization simulation
         call atomization%step()

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
   
   
end module simulation