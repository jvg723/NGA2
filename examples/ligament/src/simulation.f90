!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,      only: WP
   use hit_class,      only: hit
   use ligament_class, only: ligament
   implicit none
   private
   
   !> HIT simulation
   type(hit) :: turb
   logical :: isInHITGrp
   
   !> Ligament atomization simulation
   type(ligament) :: atom
   
   public :: simulation_init,simulation_run,simulation_final

   !> Transfer model parameters
   real(WP) :: filmthickness_over_dx  =5.0e-1_WP
   real(WP) :: min_filmthickness      =1.0e-3_WP
   real(WP) :: diam_over_filmthickness=1.0e+1_WP
   real(WP) :: max_eccentricity       =5.0e-1_WP
   real(WP) :: d_threshold            =1.0e-1_WP
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize atomization simulation
      call atom%init()
      
      ! Create an MPI group using leftmost processors only
      create_hit_group: block
         use parallel, only: group,comm
         use mpi_f08,  only: MPI_Group_incl
         integer, dimension(:), allocatable :: ranks
         integer, dimension(3) :: coord
         integer :: n,ngrp,ierr,ny,nz
         ngrp=atom%cfg%npy*atom%cfg%npz
         allocate(ranks(ngrp))
         ngrp=0
         do nz=1,atom%cfg%npz
            do ny=1,atom%cfg%npy
               ngrp=ngrp+1
               coord=[0,ny-1,nz-1]
               call MPI_CART_RANK(atom%cfg%comm,coord,ranks(ngrp),ierr)
            end do
         end do
         call MPI_Group_incl(group,ngrp,ranks,hit_group,ierr)
         if (atom%cfg%iproc.eq.1) then
            isInHITGrp=.true.
         else
            isInHITGrp=.false.
         end if
      end block create_hit_group
      
      ! Prepare HIT simulation
      if (isInHITGrp) then
         prepare_hit: block
            real(WP) :: dt
            ! Initialize HIT
            call turb%init(group=hit_group)
            ! Run HIT until t/tau_eddy=20
            dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt !< Estimate maximum stable dt
            do while (turb%time%t.lt.20.0_WP*turb%tau_tgt)
               call turb%step(dt)
            end do
         end block prepare_hit
      end if
      
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
      implicit none
      
      ! Atomization drives overall time integration
      do while (.not.atom%time%done())
         
         ! Advance atomization simulation
         call atom%step()

         ! Perform breakup
         call transfer_vf_to_drops()
         
         ! Advance HIT simulation and transfer velocity info
         if (isInHITGrp) then
            ! Advance HIT
            call turb%step(atom%time%dt)
            ! Transfer turbulent velocity from hit to rta
            apply_boundary_condition: block
               use tpns_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ihit
               real(WP) :: rescaling
               rescaling=turb%ti/turb%Urms_tgt
               call atom%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  ihit=i-atom%fs%cfg%imin+turb%fs%cfg%imax+1
                  atom%fs%U(i  ,j,k)=1.0_WP+turb%fs%U(ihit  ,j,k)*rescaling
                  atom%fs%V(i-1,j,k)=       turb%fs%V(ihit-1,j,k)*rescaling
                  atom%fs%W(i-1,j,k)=       turb%fs%W(ihit-1,j,k)*rescaling
               end do
            end block apply_boundary_condition
         end if
         
      end do
      
   end subroutine simulation_run

   !> Transfer vf to drops
   subroutine transfer_vf_to_drops()
      implicit none
      
      ! Perform a first pass with simplest CCL
      call atom%cc%build_lists(VF=atom%vf%VF,U=atom%fs%U,V=atom%fs%V,W=atom%fs%W)
      
      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam
         
         ! Loops over film segments contained locally
         do m=1,atom%cc%n_meta_struct
            
            ! Test if sphericity is compatible with transfer
            lmin=atom%cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=atom%cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=atom%cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/lmax**2)
            if (eccentricity.gt.max_eccentricity) cycle
            
            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*atom%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if (diam.eq.0.0_WP.or.diam.gt.d_threshold) cycle
            
            ! Find local structs with matching id
            do n=atom%cc%sync_offset+1,atom%cc%sync_offset+atom%cc%n_struct
               if (atom%cc%struct_list(atom%cc%struct_map_(n))%parent.ne.atom%cc%meta_structures_list(m)%id) cycle
               ! Remove liquid in meta-structure cells
               do l=1,atom%cc%struct_list(atom%cc%struct_map_(n))%nnode ! Loops over cells within local
                  i=atom%cc%struct_list(atom%cc%struct_map_(n))%node(1,l)
                  j=atom%cc%struct_list(atom%cc%struct_map_(n))%node(2,l)
                  k=atom%cc%struct_list(atom%cc%struct_map_(n))%node(3,l)
                  ! Remove liquid in that cell
                  atom%vf%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct
      
      ! Sync VF and clean up IRL and band
      call atom%vf%cfg%sync(atom%vf%VF)
      call atom%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call atom%cc%deallocate_lists()
      
      ! Perform more detailed CCL in a second pass
      atom%cc%max_interface_planes=2
      call atom%cc%build_lists(VF=atom%vf%VF,poly=atom%vf%interface_polygon,U=atom%fs%U,V=atom%fs%V,W=atom%fs%W)
      call atom%cc%get_min_thickness()
      call atom%cc%sort_by_thickness()
      
      ! Loop through identified films and remove those that are thin enough
      remove_film: block
         use mathtools, only: pi
         integer :: m,n,i,j,k,np,ip,np_old
         real(WP) :: Vt,Vl,Hl,Vd
         
         ! Loops over film segments contained locally
         do m=atom%cc%film_sync_offset+1,atom%cc%film_sync_offset+atom%cc%n_film
            
            ! Skip non-liquid films
            if (atom%cc%film_list(atom%cc%film_map_(m))%phase.ne.1) cycle
            
            ! Skip films that are still thick enough
            if (atom%cc%film_list(atom%cc%film_map_(m))%min_thickness.gt.min_filmthickness) cycle
            
            ! We are still here: transfer the film to drops
            Vt=0.0_WP      ! Transferred volume
            Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            do n=1,atom%cc%film_list(atom%cc%film_map_(m))%nnode ! Loops over cells within local film segment
               i=atom%cc%film_list(atom%cc%film_map_(m))%node(1,n)
               j=atom%cc%film_list(atom%cc%film_map_(m))%node(2,n)
               k=atom%cc%film_list(atom%cc%film_map_(m))%node(3,n)
               ! Increment liquid volume to remove
               Vl=Vl+atom%vf%VF(i,j,k)*atom%vf%cfg%vol(i,j,k)
               ! Estimate drop size based on local film thickness in current cell
               Hl=max(atom%cc%film_thickness(i,j,k),min_filmthickness)
               Vd=pi/6.0_WP*(diam_over_filmthickness*Hl)**3
               ! Remove liquid in that cell
               atom%vf%VF(i,j,k)=0.0_WP
            end do

         end do
         
      end block remove_film
      
      ! Sync VF and clean up IRL and band
      call atom%vf%cfg%sync(atom%vf%VF)
      call atom%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call atom%cc%deallocate_lists()
      
   end subroutine transfer_vf_to_drops
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize atomization simulation
      call atom%final()
      
      ! Finalize HIT simulation
      if (isInHITGrp) call turb%final()
      
   end subroutine simulation_final
   
   
end module simulation