!> Transfer model class:
!> Provides support for conversion of subgrid films and ligaments to Lagrangian droplets

module transfermodel_class
   use precision,      only: WP
   use config_class,   only: config
   use ccl_class,      only: ccl
   use lpt_class,      only: lpt
   use vfs_class,      only: vfs
   use tpns_class,     only: tpns
   use string,         only: str_medium
   use irl_fortran_interface
   implicit none
   private

   ! Expose type/methods
   public :: transfermodels

   !> Transfer model object
   type :: transfermodels

      ! This is our config
      class(config), pointer :: cfg

      type(ccl)           :: cc
      type(lpt) , pointer :: lp
      type(vfs) , pointer :: vf
      type(tpns), pointer :: fs

      !> Transfer model parameters
      real(WP) :: filmthickness_over_dx  =5.0e-1_WP
      real(WP) :: min_filmthickness      !=1.0e-3_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP
      real(WP) :: max_eccentricity       =8.0e-1_WP
      real(WP) :: d_threshold            =6.0e-1_WP
      real(WP) :: ligament_film_ratio_th =1.0e-2_WP

      !> Film breakup with retraction
      integer  :: max_nodes              =10
      logical  :: instantaneous_burst
      logical  :: has_sampled
      real(WP) :: volume_to_convert
      real(WP) :: current_sampled_volume
      
      !> Film retracting rim Rayleigh-Plateau droplet diameter
      real(WP) :: d_film_rp

      !> Droplet d^3 mean, Sauter mean, and mass median diameters
      real(WP) :: d30,d32,mmd

      !> Drop measurements
      real(WP) :: xL,yL,zL

      !> Drop position
      real(WP) :: xbary,ybary,zbary,ubary,vbary,wbary

      !> Film measurements
      real(WP) :: min_thickness,film_volume,film_ratio,my_converted_volume,converted_volume

      !> Film breakup
      logical :: burst

   contains

      procedure :: initialize
      procedure :: transfer_vf_to_drops
      procedure :: transfer_vf_to_drops_retraction
      procedure :: transfer_detached_struct
      procedure :: breakup_film_instantaneous
      procedure :: puncture_film
      procedure :: breakup_film_retraction
      procedure :: breakup_ligament
      procedure, private :: bag_droplet_gamma
      procedure :: spray_statistics_setup
      procedure :: spray_statistics
      procedure :: get_min_thickness

   end type transfermodels

contains


   !> Initialize
   subroutine initialize(this,cfg,vf,fs,lp)
      use param, only: param_read
      implicit none
      class(transfermodels), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      type(vfs), target, intent(in) :: vf
      type(tpns), target, intent(in) :: fs
      type(lpt), target, intent(in) ::  lp

      this%cfg=>cfg
      this%lp=>lp
      this%vf=>vf
      this%fs=>fs

      ! Create a connected-component labeling object
      create_and_initialize_ccl: block
         use vfs_class, only: VFlo
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl

      ! Initialize film statistics
      this%min_thickness=0.0_WP
      this%film_volume=0.0_WP
      this%film_ratio=0.0_WP
      this%converted_volume=0.0_WP
      this%my_converted_volume=0.0_WP
      this%volume_to_convert=0.0_WP
      this%current_sampled_volume=0.0_WP

      ! Input film breakup parameters
      call param_read('Instantaneous burst',this%instantaneous_burst,default=.true.)
      if (this%instantaneous_burst) then
         this%max_nodes=HUGE(1)
      else
         call param_read('Max nodes to burst',this%max_nodes,default=10)
      end if
      call param_read('Threshold film thickness',this%min_filmthickness,default=1.0e-3_WP)

      ! Initialize film burst indicator
      call param_read('Initialize as burst',this%burst,default=.false.)

      ! Setup spray statistics output
      call this%spray_statistics_setup()
      call this%spray_statistics()

   end subroutine initialize


   !> Transfer vf to drops
   subroutine transfer_vf_to_drops(this)
      use vfs_class, only: r2p
      implicit none
      class(transfermodels), intent(inout) :: this

      ! Transfer detached structs
      call this%transfer_detached_struct()
      
      if (this%vf%reconstruction_method.eq.r2p) then

         if (.not.this%burst) call this%breakup_film_instantaneous()

         if (this%burst) call this%breakup_ligament()

      end if ! r2p
     
      ! Sum converted volume
      sum_transfer: block
         use mpi_f08,  only: MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: ierr
         call MPI_REDUCE(this%my_converted_volume,this%converted_volume,1,MPI_REAL_WP,MPI_SUM,0,this%cfg%comm,ierr) 
      end block sum_transfer

      ! Resync the spray
      call this%lp%sync()

   end subroutine transfer_vf_to_drops
   

   !> Transfer vf to drops
   subroutine transfer_vf_to_drops_retraction(this,fm,dt)
      use film_class,only: film
      use vfs_class, only: r2p
      implicit none
      class(transfermodels), intent(inout) :: this
      class(film), intent(inout) :: fm
      real(WP), intent(in) :: dt

      ! Transfer detached structs
      call this%transfer_detached_struct()
      
      if (this%vf%reconstruction_method.eq.r2p) then

         if (.not.this%burst) call this%puncture_film(fm)

         if (this%burst) then
            call fm%compute_film_velocity()
            call this%breakup_film_retraction(fm,dt)
            call this%breakup_ligament()
         end if

      end if ! r2p
     
      ! Sum converted volume
      sum_transfer: block
         use mpi_f08,  only: MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: ierr
         call MPI_REDUCE(this%my_converted_volume,this%converted_volume,1,MPI_REAL_WP,MPI_SUM,0,this%cfg%comm,ierr) 
      end block sum_transfer

      ! Resync the spray
      call this%lp%sync()

   end subroutine transfer_vf_to_drops_retraction


   !> Transfer detached structures
   subroutine transfer_detached_struct(this)
      use messager,  only: die
      use string,    only: str_medium
      use mathtools, only: Pi
      implicit none
      class(transfermodels), intent(inout) :: this
      real(WP) :: lmin,lmax,eccentricity,diam
      integer :: m,n,l,i,j,k,np
      character(len=str_medium) :: filename
      integer :: iunit,ierr
      logical :: autotransfer
      real(WP), parameter :: initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3

      ! Perform a first pass with simplest CCL
      this%cc%use_struct_cutoff=.false.
      call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)

      ! Loops over global list of structures
      ! if (this%cfg%amRoot) print *,'start detached'
      do m=1,this%cc%n_meta_struct
         if (this%cc%meta_structures_list(m)%vol.gt.0.5_WP*initial_volume) then
            this%xbary=this%cc%meta_structures_list(m)%x
            this%ybary=this%cc%meta_structures_list(m)%y
            this%zbary=this%cc%meta_structures_list(m)%z
            this%ubary=this%cc%meta_structures_list(m)%u
            this%vbary=this%cc%meta_structures_list(m)%v
            this%wbary=this%cc%meta_structures_list(m)%w
         end if
         diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
         autotransfer=.false.
         ! Test if structure is at end of domain
         if (this%cc%meta_structures_list(m)%x.gt.this%cc%cfg%x(this%cc%cfg%imax-10)) autotransfer=.true.
         if (.not.autotransfer) then
            ! Test if sphericity is compatible with transfer
            lmin=this%cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=this%cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=this%cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/(lmax**2+tiny(1.0_WP)))
            ! if (this%cfg%amRoot) print *,'id',this%cc%meta_structures_list(m)%id,'eccentricity',eccentricity,'lengths',this%cc%meta_structures_list(m)%lengths

            ! Test if diameter is compatible with transfer
            ! if (this%cfg%amRoot) print *,'id',this%cc%meta_structures_list(m)%id,'diam',diam
            if (eccentricity.gt.this%max_eccentricity) cycle
            if ((diam.eq.0.0_WP).or.(diam.gt.this%d_threshold)) cycle
         end if
         
         ! Create drop from available liquid volume - only one root does that
         if (this%cc%cfg%amRoot) then
            ! Make room for new drop
            np=this%lp%np_+1; call this%lp%resize(np)
            ! Add the drop
            this%lp%p(np)%id  =int(1,8)                                                                                 !< Give id (maybe based on break-up model?)
            this%lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
            this%lp%p(np)%Acol =0.0_WP                                                                                   !< Give zero collision force
            this%lp%p(np)%Tcol =0.0_WP                                                                                   !< Give zero collision force
            this%lp%p(np)%d   =diam                                                                                     !< Assign diameter to account for full volume
            this%lp%p(np)%pos =[this%cc%meta_structures_list(m)%x,this%cc%meta_structures_list(m)%y,this%cc%meta_structures_list(m)%z] !< Place the drop at the liquid barycenter
            this%lp%p(np)%vel =[this%cc%meta_structures_list(m)%u,this%cc%meta_structures_list(m)%v,this%cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
            this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])                !< Place the drop in the proper cell for the this%lp%cfg
            this%lp%p(np)%flag=0                                                                                        !< Activate it
            ! Increment particle counter
            this%lp%np_=np

            !!! Write to droplet list !!!
            ! Open the file
            filename='spray-all/droplets'
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameter, velocity, and position
            write(iunit,*) this%lp%p(np)%d,this%lp%p(np)%vel(1),this%lp%p(np)%vel(2),this%lp%p(np)%vel(3),norm2([this%lp%p(np)%vel(1),this%lp%p(np)%vel(2),this%lp%p(np)%vel(3)]),this%lp%p(np)%pos(1),this%lp%p(np)%pos(2),this%lp%p(np)%pos(3),'detached',this%lp%p(np)%id  
            ! Close the file
            close(iunit)

         end if

         ! Find local structs with matching id
         do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
            if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(m)%id) cycle
            ! Remove liquid in meta-structure cells
            do l=1,this%cc%struct_list(this%cc%struct_map_(n))%nnode ! Loops over cells within local
               i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
               j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
               k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
               ! Remove liquid in that cell
               this%my_converted_volume=this%my_converted_volume+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
               this%vf%VF(i,j,k)=0.0_WP
            end do
         end do
         
      end do
         
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()

   end subroutine transfer_detached_struct


   !> Transfer vf to drops
   subroutine breakup_film_instantaneous(this)
      use messager,  only: die
      use string,    only: str_medium
      use mathtools, only: Pi,normalize,cross_product
      use random,    only: random_uniform,random_gamma
      use parallel,  only: MPI_REAL_WP
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_LOGICAL,MPI_LOR
      use irl_fortran_interface
      implicit none
      class(transfermodels), intent(inout) :: this
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,id_largest_film
      real(WP), parameter :: initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3
      integer :: m,n,l,i,j,k,np,ip,np_old,np_start
      real(WP) :: Vt,Vl,Vd,diam,d0
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP) :: theta
      logical :: has_burst
      real(WP), dimension(2) :: s_bary,c_bary
      integer :: nodes_transferred
      ! Conversion model
      real(WP) :: alpha, beta, curv_sum, ncurv
      real(WP) :: max_film_droplet_diam

      ! Perform more detailed CCL in a second pass
      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      ! call this%cc%film_classify(this%vf%Lbary,this%vf%Gbary)
      call this%vf%sense_interface()
      ! call this%cc%get_min_thickness()
      call this%get_min_thickness(this%min_thickness,this%film_volume,id_largest_film)
      this%film_ratio=this%film_volume/initial_volume
      call this%cc%sort_by_thickness()

      ! Assign droplet extents (assuming there is only one metastruct)
      droplet_extent: block
         this%xL=0.0_WP;this%yL=0.0_WP;this%zL=0.0_WP
         do m=1,this%cc%n_meta_struct
            if (this%cc%meta_structures_list(m)%vol.le.epsilon(1.0_WP)) cycle
            this%xL=max(this%xL,this%cc%meta_structures_list(m)%xL)
            this%yL=max(this%yL,this%cc%meta_structures_list(m)%yL)
            this%zL=max(this%zL,this%cc%meta_structures_list(m)%zL)
         end do
      end block droplet_extent

      ! Loops over film segments contained locally
      d0=1.0_WP ! Characteristic diameter
      np_start=this%lp%np_  ! Remember old number of particles
      has_burst=.false.

      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
      
         ! Skip non-liquid films
         if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
      
         ! Skip films that are still thick enough
         if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.this%min_filmthickness) cycle

         ! We are still here: transfer the film to drops
         has_burst=.true.
         i=this%cc%film_list(this%cc%film_map_(m))%node(1,1)
         j=this%cc%film_list(this%cc%film_map_(m))%node(2,1)
         k=this%cc%film_list(this%cc%film_map_(m))%node(3,1)
         ! Compute cell-averaged curvature for cell containing thinnest film segment
         curv_sum=0.0_WP; ncurv=0.0_WP
         do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
            if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
               curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
               ncurv=ncurv+1.0_WP
            end if
         end do
         ! Determine gamma distribution parameters                        
         call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
         Vt=0.0_WP      ! Transferred volume
         Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
         np_old=this%lp%np_  ! Remember old number of particles
         max_film_droplet_diam=2.0_WP*this%d_film_rp
         Vd=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            ! Get local coordinate system
            nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
            select case (maxloc(abs(nref),1))
            case (1)
               tref=normalize([+nref(2),-nref(1),0.0_WP])
            case (2)
               tref=normalize([0.0_WP,+nref(3),-nref(2)])
            case (3)
               tref=normalize([-nref(3),0.0_WP,+nref(1)])
            end select
            sref=cross_product(nref,tref)
            ! Increment liquid volume to remove
            Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)

            ! Compute cell-averaged curvature
            curv_sum=0.0_WP; ncurv=0.0_WP
            do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                  curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                  ncurv=ncurv+1.0_WP
               end if
            end do
            ! Determine gamma distribution parameters
            call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
            
            ! Create drops from available liquid volume
            do while (Vl-Vd.gt.0.0_WP)
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(2,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
               ! theta        =random_uniform(0.0_WP,2.0_WP*pi)           !< Angular position of randomly placed droplet relative to liquid barycenter
               ! this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(0.0_WP,0.5_WP*this%vf%cfg%meshsize(i,j,k))*(cos(theta)*tref+sin(theta)*sref)
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref              
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
               ! Update tracked volumes
               Vl=Vl-Vd
               Vt=Vt+Vd
               ! Output diameter, velocity, and position
               ! Generate new droplet volume
               max_film_droplet_diam=2.0_WP*this%d_film_rp
               Vd=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
            end do

            ! Remove liquid in that cell
            this%my_converted_volume=this%my_converted_volume+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP

         end do
         
         ! Based on how many particles were created, decide what to do with left-over volume
         if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
            ! Add one last drop for remaining liquid volume
            np=this%lp%np_+1; call this%lp%resize(np)
            ! Add the drop
            this%lp%p(np)%id  =int(3,8)                                   !< Give id (maybe based on break-up model?)
            this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
            this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
            this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
            this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
            this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                     !< Place the drop at the liquid barycenter
            this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
            this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
            this%lp%p(np)%flag=0                                          !< Activate it
            ! Increment particle counter
            this%lp%np_=np
         else ! Some particles were created, make them all larger
            do ip=np_old+1,this%lp%np_
               this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
            end do
         end if

      end do ! local films
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()                  

      ! Mark as burst if a film has burst
      call MPI_ALLREDUCE(has_burst,this%burst,1,MPI_LOGICAL,MPI_LOR,this%cfg%comm,ierr)

      ! Clean up CCL
      call this%cc%deallocate_lists()

      ! Write the diameters and velocities
      filename='spray-all/droplets'
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            ! open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameters and velocities
            do i=np_start+1,this%lp%np_
               write(iunit,*) this%lp%p(i)%d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),norm2([this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]),this%lp%p(i)%pos(1),this%lp%p(i)%pos(2),this%lp%p(i)%pos(3),'film',this%lp%p(np)%id  
            end do
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine breakup_film_instantaneous


   !> Puncture film
   subroutine puncture_film(this,fm)
      use film_class,only: film
      use messager,  only: die
      use string,    only: str_medium
      use mathtools, only: Pi,normalize,cross_product
      use random,    only: random_uniform,random_gamma
      use parallel,  only: MPI_REAL_WP
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_LOGICAL,MPI_LOR
      use irl_fortran_interface
      implicit none
      class(transfermodels), intent(inout) :: this
      class(film), intent(inout) :: fm
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,id_largest_film
      real(WP), parameter :: initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3
      integer :: m,n,l,i,j,k,np,ip,np_old,np_start,ii,jj,kk
      real(WP) :: Vt,Vl,Vd,diam,d0
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP) :: theta
      logical :: has_burst
      real(WP), dimension(2) :: s_bary,c_bary
      integer :: nodes_transferred
      ! Conversion model
      real(WP) :: alpha, beta, curv_sum, ncurv
      real(WP) :: max_film_droplet_diam

      ! Perform more detailed CCL in a second pass
      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      ! call this%cc%film_classify(this%vf%Lbary,this%vf%Gbary)
      call this%vf%sense_interface()
      ! call this%cc%get_min_thickness()
      call this%get_min_thickness(this%min_thickness,this%film_volume,id_largest_film)
      this%film_ratio=this%film_volume/initial_volume
      call this%cc%sort_by_thickness()

      ! Clear temporary edge indicator
      fm%edgeold=0.0_WP

      ! Loops over film segments contained locally
      d0=1.0_WP ! Characteristic diameter
      np_start=this%lp%np_  ! Remember old number of particles
      has_burst=.false.
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
      
         ! Skip non-liquid films
         if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
      
         ! Skip films that are still thick enough
         if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.this%min_filmthickness) cycle

         ! Skip film if processor does not contain minimum thickness
         ! This step should be refined to account for punch hole radius
         if (this%cc%film_list(this%cc%film_map_(m))%minloc_rank.ne.this%cfg%rank) cycle

         ! We are still here: transfer the film to drops
         has_burst=.true.
         i=this%cc%film_list(this%cc%film_map_(m))%node(1,1)
         j=this%cc%film_list(this%cc%film_map_(m))%node(2,1)
         k=this%cc%film_list(this%cc%film_map_(m))%node(3,1)
         ! Compute cell-averaged curvature for cell containing thinnest film segment
         curv_sum=0.0_WP; ncurv=0.0_WP
         do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
            if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
               curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
               ncurv=ncurv+1.0_WP
            end if
         end do
         ! Determine gamma distribution parameters                        
         call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
         Vt=0.0_WP      ! Transferred volume
         Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
         np_old=this%lp%np_  ! Remember old number of particles
         max_film_droplet_diam=2.0_WP*this%d_film_rp
         Vd=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            if (n.gt.this%max_nodes) exit
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            ! Get local coordinate system
            nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
            select case (maxloc(abs(nref),1))
            case (1)
               tref=normalize([+nref(2),-nref(1),0.0_WP])
            case (2)
               tref=normalize([0.0_WP,+nref(3),-nref(2)])
            case (3)
               tref=normalize([-nref(3),0.0_WP,+nref(1)])
            end select
            sref=cross_product(nref,tref)
            ! Increment liquid volume to remove
            Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)

            ! Compute cell-averaged curvature
            curv_sum=0.0_WP; ncurv=0.0_WP
            do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                  curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                  ncurv=ncurv+1.0_WP
               end if
            end do
            ! Determine gamma distribution parameters
            call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
            
            ! Create drops from available liquid volume
            do while (Vl-Vd.gt.0.0_WP)
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(4,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
               this%lp%p(np)%d   =(6.0_WP*Vd/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref              
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
               ! Update tracked volumes
               Vl=Vl-Vd
               Vt=Vt+Vd
               ! Output diameter, velocity, and position
               ! Generate new droplet volume
               max_film_droplet_diam=2.0_WP*this%d_film_rp
               Vd=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
            end do

            ! Remove liquid in that cell
            this%my_converted_volume=this%my_converted_volume+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
            ! Initialize edge indicator
            fm%edgeold(i,j,k)=1.0_WP ! for now, just to propagate to neighbors after sync

         end do
         
         ! Based on how many particles were created, decide what to do with left-over volume
         if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
            ! Add one last drop for remaining liquid volume
            np=this%lp%np_+1; call this%lp%resize(np)
            ! Add the drop
            this%lp%p(np)%id  =int(5,8)                                   !< Give id (maybe based on break-up model?)
            this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
            this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
            this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
            this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
            this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                     !< Place the drop at the liquid barycenter
            this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
            this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
            this%lp%p(np)%flag=0                                          !< Activate it
            ! Increment particle counter
            this%lp%np_=np
         else ! Some particles were created, make them all larger
            do ip=np_old+1,this%lp%np_
               this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
            end do
         end if

      end do ! local films
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%cfg%sync(fm%edgeold)
      call this%vf%clean_irl_and_band()                  

      ! Propagate hole edge
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            do kk=k-1,k+1
               do jj=j-1,j+1
                  do ii=i-1,i+1                  
                     if (fm%edgeold(ii,jj,kk).ge.1.0_WP) then
                        fm%edge(i,j,k)=2.0_WP*ceiling(this%vf%VF(i,j,k))
                     end if
                  end do
               end do
            end do
         end do 
      end do

      ! Sync edge indicator
      call this%vf%cfg%sync(fm%edge)

      ! Mark as burst if a film has burst
      call MPI_ALLREDUCE(has_burst,this%burst,1,MPI_LOGICAL,MPI_LOR,this%cfg%comm,ierr)

      ! Clean up CCL
      call this%cc%deallocate_lists()

      ! Write the diameters and velocities
      filename='spray-all/droplets'
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            ! open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameters and velocities
            do i=np_start+1,this%lp%np_
               write(iunit,*) this%lp%p(i)%d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),norm2([this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]),this%lp%p(i)%pos(1),this%lp%p(i)%pos(2),this%lp%p(i)%pos(3),'film',this%lp%p(np)%id  
            end do
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine puncture_film


   !> Transfer vf to drops
   subroutine breakup_film_retraction(this,fm,dt)
      use film_class,only: film
      use messager,  only: die
      use string,    only: str_medium
      use mathtools, only: Pi,normalize,cross_product
      use random,    only: random_uniform,random_gamma
      use parallel,  only: MPI_REAL_WP
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_LOGICAL,MPI_LOR
      use irl_fortran_interface
      implicit none
      class(transfermodels), intent(inout) :: this
      class(film), intent(inout) :: fm
      real(WP), intent(in)       :: dt
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,id_largest_film
      real(WP), parameter :: initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3
      integer :: m,n,l,i,j,k,np,ip,np_old,np_start,ii,jj,kk
      real(WP) :: Vt,Vl,Vd,diam,d0
      real(WP), dimension(3) :: pref,nref,tref,sref
      real(WP) :: theta
      logical :: has_burst
      real(WP), dimension(2) :: s_bary,c_bary
      integer :: nodes_transferred
      ! Conversion model
      real(WP) :: alpha, beta, curv_sum, ncurv
      real(WP) :: max_film_droplet_diam
      ! Droplet shedding
      real(WP) :: qout,qin,b,a,length,Qdt
      real(WP) :: total_edge_volume,local_volume_to_convert,Vtsum

      ! Perform more detailed CCL in a second pass
      call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%vf%sense_interface()
      call fm%detect_edge_regions()
      call this%get_min_thickness(this%min_thickness,this%film_volume,id_largest_film)
      this%film_ratio=this%film_volume/initial_volume
      call this%cc%sort_by_thickness()

      ! Assign droplet extents (assuming there is only one metastruct)
      droplet_extent: block
         this%xL=0.0_WP;this%yL=0.0_WP;this%zL=0.0_WP
         do m=1,this%cc%n_meta_struct
            if (this%cc%meta_structures_list(m)%vol.le.epsilon(1.0_WP)) cycle
            this%xL=max(this%xL,this%cc%meta_structures_list(m)%xL)
            this%yL=max(this%yL,this%cc%meta_structures_list(m)%yL)
            this%zL=max(this%zL,this%cc%meta_structures_list(m)%zL)
         end do
      end block droplet_extent

      ! Clear temporary edge indicator
      fm%edgeold=0.0_WP
      
      ! Loops over film segments contained locally
      d0=1.0_WP ! Characteristic diameter
      np_start=this%lp%np_  ! Remember old number of particles
      has_burst=.false.
      Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
      ! Determine how much edge volume is contained locally
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
         ! Skip film if not largest film
         if (this%cc%film_list(this%cc%film_map_(m))%parent.ne.id_largest_film) cycle
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            ! if ((this%vf%edge_sensor(i,j,k).lt.fm%edge_threshold.or.fm%edge(i,j,k).le.0.0_WP).and.(fm%edge(i,j,k).lt.2.0_WP)) cycle
            if (fm%edge(i,j,k).lt.1.0_WP) cycle
            ! Increment liquid volume to remove
            qin=this%vf%thickness(i,j,k)*norm2([fm%filmUm(i,j,k),fm%filmVm(i,j,k),fm%filmWm(i,j,k)])
            qout=qin ! assume constant rim thickness for now
            ! length=this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)/(this%vf%thickness(i,j,k)+epsilon(1.0_WP))/this%vf%cfg%meshsize(i,j,k)
            length=this%vf%cfg%meshsize(i,j,k)
            Qdt=qout*length*dt
            ! Qdt=0.5_WP*this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
            Vl=Vl+min(this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k),Qdt)
            ! Remove liquid in that cell
            this%my_converted_volume=this%my_converted_volume+min(this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k),Qdt)
            this%vf%VF(i,j,k)=max(0.0_WP,this%vf%VF(i,j,k)-Qdt/this%vf%cfg%vol(i,j,k))
            ! Initialize edge indicator
            if (this%vf%VF(i,j,k).le.0.0_WP) fm%edgeold(i,j,k)=1.0_WP ! for now, just to propagate to neighbors after sync
         end do
      end do
      ! Sum edge volume across procs
      call MPI_ALLREDUCE(Vl,total_edge_volume,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      this%volume_to_convert=this%volume_to_convert+total_edge_volume
      ! local_volume_to_convert=Vl/(total_edge_volume+epsilon(1.0_WP))*this%volume_to_convert

      ! Create droplets
      ! if (local_volume_to_convert.gt.0.0_WP) then
      Vt=0.0_WP      ! Transferred volume
      if (Vl.gt.0.0_WP) then
         local_volume_to_convert=Vl/total_edge_volume*this%volume_to_convert

         do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
            ! Skip film if not largest film
            if (this%cc%film_list(this%cc%film_map_(m))%parent.ne.id_largest_film) cycle
            ! We are still here: transfer the film edge to drops
            print *,'rank',this%cfg%rank,'film m',m,'shed droplets'
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,1)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,1)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,1)
            ! Compute cell-averaged curvature for cell containing thinnest film segment
            curv_sum=0.0_WP; ncurv=0.0_WP
            do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                  curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                  ncurv=ncurv+1.0_WP
               end if
            end do
            ! Determine gamma distribution parameters                        
            ! Vt=0.0_WP      ! Transferred volume
            ! Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
            np_old=this%lp%np_  ! Remember old number of particles
            if (.not.this%has_sampled) then
               ! call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
               call this%bag_droplet_gamma(this%min_filmthickness,2.0_WP*ncurv/curv_sum,alpha,beta)
               ! call this%bag_droplet_gamma(max(this%min_filmthickness,this%vf%thickness(i,j,k)),2.0_WP*ncurv/curv_sum,alpha,beta)
               max_film_droplet_diam=2.0_WP*this%d_film_rp
               this%current_sampled_volume=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
               this%has_sampled=.true.
            end if
            do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            ! do n=this%cc%film_list(this%cc%film_map_(m))%nnode,1,-1 ! Loops over cells within local film segment in reverse

               if (local_volume_to_convert.lt.this%current_sampled_volume) then
                  print *,'rank',this%cfg%rank,'film m',m,'node',n,'local_volume_to_convert < Vd'
                  exit
               end if                        
               i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
               j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
               k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
               ! if ((this%vf%edge_sensor(i,j,k).lt.fm%edge_threshold.or.fm%edge(i,j,k).le.0.0_WP).and.(fm%edge(i,j,k).lt.2.0_WP)) cycle
               if (fm%edge(i,j,k).lt.1.0_WP) cycle
               ! Get local coordinate system
               nref=calculateNormal(this%vf%interface_polygon(1,i,j,k))
               select case (maxloc(abs(nref),1))
               case (1)
                  tref=normalize([+nref(2),-nref(1),0.0_WP])
               case (2)
                  tref=normalize([0.0_WP,+nref(3),-nref(2)])
               case (3)
                  tref=normalize([-nref(3),0.0_WP,+nref(1)])
               end select
               sref=cross_product(nref,tref)

               ! Compute cell-averaged curvature
               curv_sum=0.0_WP; ncurv=0.0_WP
               do l=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                  if (getNumberOfVertices(this%vf%interface_polygon(l,i,j,k)).gt.0) then
                     curv_sum=curv_sum+abs(this%vf%curv2p(l,i,j,k))
                     ncurv=ncurv+1.0_WP
                  end if
               end do
               ! Determine gamma distribution parameters
               ! call this%bag_droplet_gamma(this%vf%thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
               call this%bag_droplet_gamma(this%min_filmthickness,2.0_WP*ncurv/curv_sum,alpha,beta)
               ! call this%bag_droplet_gamma(max(this%min_filmthickness,this%vf%thickness(i,j,k)),2.0_WP*ncurv/curv_sum,alpha,beta)
               ! Create drops from available liquid volume
               ! Make room for new drop
               np=this%lp%np_+1; call this%lp%resize(np)
               ! Add the drop
               this%lp%p(np)%id  =int(6,8)                                   !< Give id (maybe based on break-up model?)
               this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
               this%lp%p(np)%Acol=0.0_WP                                     !< Give zero collision force
               this%lp%p(np)%Tcol=0.0_WP                                     !< Give zero collision force
               this%lp%p(np)%d   =(6.0_WP*this%current_sampled_volume/pi)**(1.0_WP/3.0_WP)            !< Assign diameter from model above
               ! theta        =random_uniform(0.0_WP,2.0_WP*pi)           !< Angular position of randomly placed droplet relative to liquid barycenter
               ! this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(0.0_WP,0.5_WP*this%vf%cfg%meshsize(i,j,k))*(cos(theta)*tref+sin(theta)*sref)
               this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*tref+random_uniform(-0.5_WP*this%vf%cfg%meshsize(i,j,k),0.5_WP*this%vf%cfg%meshsize(i,j,k))*sref
               this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
               this%lp%p(np)%vel =this%lp%p(np)%vel+this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=fm%filmU,V=fm%filmV,W=fm%filmW) 
               this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])    !< Place the drop in the proper cell for the this%lp%cfg
               this%lp%p(np)%flag=0                                          !< Activate it
               ! Increment particle counter
               this%lp%np_=np
               ! Update tracked volumes
               ! Vl=Vl-this%current_sampled_volume
               Vt=Vt+this%current_sampled_volume
               local_volume_to_convert=local_volume_to_convert-this%current_sampled_volume
               this%my_converted_volume=this%my_converted_volume+this%current_sampled_volume
               ! Output diameter, velocity, and position
               ! Generate new droplet volume
               max_film_droplet_diam=2.0_WP*this%d_film_rp
               this%current_sampled_volume=pi/6.0_WP*(min(random_gamma(alpha)*beta*d0,max_film_droplet_diam))**3
               ! if (local_volume_to_convert.lt.Vd) then
               !    print *,'ts',time%n,'rank',this%vf%cfg%rank,'film m',m,'local_volume_to_convert < Vd'
               !    exit
               ! end if

               ! ! Remove liquid in that cell
               ! this%vf%VF(i,j,k)=max(0.0_WP,this%vf%VF(i,j,k)-Qdt/this%vf%cfg%vol(i,j,k))
               ! ! Initialize edge indicator
               ! if (this%vf%VF(i,j,k).le.0.0_WP) fm%edgeold(i,j,k)=1.0_WP ! for now, just to propagate to neighbors after sync
            end do ! nnode

            ! ! Based on how many particles were created, decide what to do with left-over volume
            ! ! only when film is completely drained
            ! if (n.eq.this%cc%film_list(this%cc%film_map_(m))%nnode) then
            !    if (Vt.eq.0.0_WP.and.local_volume_to_convert.gt.0.0_WP) then ! No particle was created, we need one...
            !       ! Add one last drop for remaining liquid volume
            !       np=this%lp%np_+1; call this%lp%resize(np)
            !       ! Add the drop
            !       this%lp%p(np)%id  =int(7,8)                                   !< Give id (maybe based on break-up model?)
            !       this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
            !       this%lp%p(np)%Acol =0.0_WP                                     !< Give zero collision force
            !       this%lp%p(np)%Tcol =0.0_WP                                                                                   !< Give zero collision force
            !       this%lp%p(np)%d   =(6.0_WP*local_volume_to_convert/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
            !       this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
            !       this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
            !       this%lp%p(np)%vel =this%lp%p(np)%vel+[fm%filmUm(i,j,k),fm%filmVm(i,j,k),fm%filmWm(i,j,k)]
            !       this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
            !       this%lp%p(np)%flag=0                                          !< Activate it
            !       ! Increment particle counter
            !       this%lp%np_=np
            !    else ! Some particles were created, make them all larger
            !       do ip=np_old+1,this%lp%np_
            !          this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+local_volume_to_convert)/Vt)**(1.0_WP/3.0_WP)
            !       end do
            !    end if
            ! end if

         end do ! local films
      end if ! local_volume_to_convert.gt.0.0_WP

      ! Sync converted volume
      call MPI_ALLREDUCE(Vt,Vtsum,1,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      this%volume_to_convert=this%volume_to_convert-Vtsum

      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%cfg%sync(fm%edgeold)
      call this%vf%clean_irl_and_band()

      ! ! Project barycenters
      ! do m=1,sum(this%vf%band_count(0:1))
      !    i=this%vf%band_map(1,m)
      !    j=this%vf%band_map(2,m)
      !    k=this%vf%band_map(3,m)
         
      !    ! Skip wall/bcond cells - bconds need to be provided elsewhere directly!
      !    if (this%vf%mask(i,j,k).ne.0) cycle
      !    ! Project forward in time
      !    this%vf%Lbary(:,i,j,k)=this%vf%project(this%vf%Lbary(:,i,j,k),i,j,k,dt,fm%filmU,fm%filmV,fm%filmW)
      !    this%vf%Gbary(:,i,j,k)=this%vf%project(this%vf%Gbary(:,i,j,k),i,j,k,dt,fm%filmU,fm%filmV,fm%filmW)
      ! end do
      
      ! ! Synchronize and clean-up barycenter fields
      ! call this%vf%sync_and_clean_barycenters()
      
      ! Propagate hole edge
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            ! if (this%vf%edge_sensor(i,j,k).lt.0.5_WP*fm%edge_threshold) cycle
            do kk=k-1,k+1
               do jj=j-1,j+1
                  do ii=i-1,i+1                  
                     if (fm%edgeold(ii,jj,kk).ge.1.0_WP) then
                        fm%edge(i,j,k)=ceiling(this%vf%VF(i,j,k))
                     end if
                  end do
               end do
            end do
         end do 
      end do

      ! Sync edge indicator
      call this%vf%cfg%sync(fm%edge)

      ! Clean up CCL
      call this%cc%deallocate_lists()

      ! Write the diameters and velocities
      filename='spray-all/droplets'
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            ! open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr)
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameters and velocities
            do i=np_start+1,this%lp%np_
               write(iunit,*) this%lp%p(i)%d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),norm2([this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]),this%lp%p(i)%pos(1),this%lp%p(i)%pos(2),this%lp%p(i)%pos(3),'film',this%lp%p(np)%id  
            end do
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine breakup_film_retraction


   !> Break up ligaments with a nonlinear capillary breakup model
   !> Adapted from Kim and Moin 2020, Combustion Science and Technology
   subroutine breakup_ligament(this)
      use mathtools, only: pi,twoPi,normalize,cross_product
      use random,    only: random_uniform
      use messager,  only: die
      use irl_fortran_interface
      implicit none
      class(transfermodels), intent(inout) :: this
      integer :: m,n,l,i,j,k,np,ip
      real(WP) :: Vt,Vl,Vd
      integer  :: nmain,nsat,np_old,np_start
      real(WP) :: minor_radius,diam,Vrim,Lrim
      ! real(WP) :: lmin,lmax,eccentricity
      ! Prescribe inviscid breakup parameters
      real(WP), parameter :: min_diam=1.0e-2_WP
      integer , parameter :: max_drops=1.0e2
      real(WP), parameter :: dimless_wavenumber=0.697_WP
      real(WP), parameter :: size_ratio=0.707_WP !0.015_WP
      integer :: ierr,iunit,rank,id_largest_film
      character(len=str_medium) :: filename
      real(WP), parameter :: initial_volume=4.0_WP/3.0_WP*Pi*0.5_WP**3
      real(WP), parameter :: min_sphericity=2.0_WP

      ! Perform ligament-identifying CCL      
      this%cc%use_struct_cutoff=.true.
      this%cc%struct_thickness_cutoff=3.0_WP
      call this%cc%build_lists(this%vf%VF,this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%cc%struct_classify(this%vf%Lbary,this%vf%Gbary)
      call this%get_min_thickness(this%min_thickness,this%film_volume,id_largest_film)
      this%film_ratio=this%film_volume/initial_volume
      ! call this%cc%sort_struct_by_thickness()

      ! Remember old number of particles
      np_start=this%lp%np_

      ! if (this%film_ratio.lt.1.0e-5_WP) then
      if (.true.) then
         ! if (this%cfg%amRoot) print *,'start ligament breakup tests'
         do m=1,this%cc%n_meta_struct
            ! Avoids erroneous conversion of ligament tips; small detached structures are converted elsewhere
            ! if (this%cfg%amRoot) print *,'m',m,'vol barrier',this%cc%meta_structures_list(m)%vol.lt.2.0_WP*this%cfg%min_meshsize**3,'vol',this%cc%meta_structures_list(m)%vol
            if (this%cc%meta_structures_list(m)%vol.lt.2.0_WP*this%cfg%min_meshsize**3) cycle

            ! Test if ligament is long enough to break up
            Lrim=this%cc%meta_structures_list(m)%length !s(1)
            Vrim=this%cc%meta_structures_list(m)%vol
            minor_radius=sqrt(Vrim/pi/Lrim)                  
            ! Drop size method from Kim & Moin (2011)
            nmain=floor(dimless_wavenumber*Lrim/twoPi/minor_radius)
            ! if (this%cfg%amRoot) print *,'m',m,'nmain barrier',nmain.lt.1,'nmain',nmain
            if (nmain.lt.1) cycle
            ! if (this%cc%meta_structures_list(m)%type.ne.1) cycle
            ! Sphericity test
            ! diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            ! if (this%cfg%amRoot) print *,'eq. sph. diam',diam
            ! if (this%cc%meta_structures_list(m)%lengths(1)/diam.lt.min_sphericity) cycle
            ! if (this%cc%meta_structures_list(m)%lengths(2)/diam.gt.1.0_WP) cycle
            ! if (this%cc%meta_structures_list(m)%length/diam.lt.min_sphericity) cycle
            ! if (2.0_WP*minor_radius/diam.gt.1.0_WP) cycle

            ! Thickness test
            ! if (this%cfg%amRoot) print *,'m',m,'thickness barrier',this%cc%meta_structures_list(m)%min_thickness.gt.1.0_WP*this%cfg%min_meshsize,'thickness',this%cc%meta_structures_list(m)%min_thickness
            if (this%cc%meta_structures_list(m)%min_thickness.gt.2.0_WP*this%cfg%min_meshsize) cycle
      
            ! Type ratio test
            ! if (this%cfg%amRoot) print *,'m',m,'type ratio barrier',this%cc%meta_structures_list(m)%f_ligament.lt.0.9_WP,'flig',this%cc%meta_structures_list(m)%f_ligament,'nnode',this%cc%meta_structures_list(m)%nnode
            if (this%cc%meta_structures_list(m)%f_ligament.lt.0.9_WP) cycle

            ! Compute main droplet diameter
            nsat=nmain+1
            diam=(6.0_WP*Vrim/pi/(real(nmain,WP)+size_ratio**3*real(nsat,WP)))**(1.0_WP/3.0_WP)
            diam=max(diam,min_diam)

            if (nmain.gt.1) then
               do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
                  if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(m)%id) cycle
                  Vd=pi/6.0_WP*(diam**3+size_ratio*diam**3)
                  Vt=0.0_WP      ! Transferred volume
                  Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation            
                  np_old=this%lp%np_  ! Remember old number of particles
                  do l=this%cc%struct_list(this%cc%struct_map_(n))%nnode,1,-1 ! Loops over cells within local film segment
                     i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
                     j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
                     k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
                     ! Increment liquid volume to remove
                     Vl=Vl+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
                     ! Create drops from available liquid volume
                     do while (Vl-Vd.gt.0.0_WP.and.this%lp%np_-np_start.lt.max_drops)            
                        ! Make room for new drop
                        np=this%lp%np_+1; call this%lp%resize(np)
                        ! Add the drop
                        this%lp%p(np)%id  =int(8,8)                                                                                   !< Give id (maybe based on break-up model?)
                        this%lp%p(np)%dt  =0.0_WP                                                                                     !< Let the drop find it own integration time
                        this%lp%p(np)%Acol=0.0_WP                                                                                     !< Give zero collision force
                        this%lp%p(np)%Tcol=0.0_WP                                                                                     !< Give zero collision force
                        this%lp%p(np)%d   =diam                                                                                       !< Assign diameter to account for full volume
                        this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                                                                     !< Place the drop at the liquid barycenter
                        this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np-1)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                        this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np-1)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])                  !< Place the drop in the proper cell for the this%lp%cfg
                        this%lp%p(np)%flag=0                                                                                          !< Activate it
                        ! Increment particle counter
                        this%lp%np_=np
                        ! Make room for new drop
                        np=this%lp%np_+1; call this%lp%resize(np)
                        ! Add the drop
                        this%lp%p(np)%id  =int(9,8)                                                                                   !< Give id (maybe based on break-up model?)
                        this%lp%p(np)%dt  =0.0_WP                                                                                     !< Let the drop find it own integration time
                        this%lp%p(np)%Acol=0.0_WP                                                                                     !< Give zero collision force
                        this%lp%p(np)%Tcol=0.0_WP                                                                                     !< Give zero collision force
                        this%lp%p(np)%d   =size_ratio*diam                                                                                       !< Assign diameter to account for full volume
                        this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)+2.0_WP*diam*this%cc%meta_structures_list(m)%axes(:,1)          
                        this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W)    !< Interpolate local cell velocity as drop velocity
                        this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])                  !< Place the drop in the proper cell for the this%lp%cfg
                        this%lp%p(np)%flag=0                                                                                          !< Activate it
                        ! Increment particle counter
                        this%lp%np_=np
                        ! Update tracked volumes
                        Vl=Vl-Vd
                        Vt=Vt+Vd
                        if (this%lp%np_-np_start.ge.max_drops) print *,'MAX DROPS EXCEEDED',this%lp%np_,np_start,'Vrim',Vrim,'nmain',nmain,'Lrim',Lrim,'min rad',minor_radius
                     end do
                     
                     ! Remove liquid in that cell
                     this%my_converted_volume=this%my_converted_volume+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
                     this%vf%VF(i,j,k)=0.0_WP
                     
                  end do

                  ! Based on how many particles were created, decide what to do with left-over volume
                  if (Vl.gt.0.0_WP) then
                     if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                        ! Add one last drop for remaining liquid volume
                        np=this%lp%np_+1; call this%lp%resize(np)
                        ! Add the drop
                        this%lp%p(np)%id  =int(10,8)                                   !< Give id (maybe based on break-up model?)
                        this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                        this%lp%p(np)%Acol=0.0_WP                                    !< Give zero collision force
                        this%lp%p(np)%Tcol=0.0_WP                                    !< Give zero collision force
                        this%lp%p(np)%d   =(6.0_WP*Vl/pi)**(1.0_WP/3.0_WP)            !< Assign diameter based on remaining liquid volume
                        this%lp%p(np)%pos =this%vf%Lbary(:,i,j,k)                          !< Place the drop at the liquid barycenter
                        this%lp%p(np)%vel =this%fs%cfg%get_velocity(pos=this%lp%p(np)%pos,i0=i,j0=j,k0=k,U=this%fs%U,V=this%fs%V,W=this%fs%W) !< Interpolate local cell velocity as drop velocity
                        this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin]) !< Place the drop in the proper cell for the this%lp%cfg
                        this%lp%p(np)%flag=0                                          !< Activate it
                        ! Increment particle counter
                        this%lp%np_=np
                     else ! Some particles were created, make them all larger
                        do ip=np_old+1,this%lp%np_
                           this%lp%p(ip)%d=this%lp%p(ip)%d*((Vt+Vl)/Vt)**(1.0_WP/3.0_WP)
                        end do
                     end if
                  end if
               end do ! n_struct
            else ! nmain=1
               ! Create drop from available liquid volume - only one root does that
               if (this%cc%cfg%amRoot) then
                  ! Make room for new drop
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(11,8)                                                                                !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                                                                   !< Let the drop find it own integration time
                  this%lp%p(np)%Acol=0.0_WP                                                                                  !< Give zero collision force
                  this%lp%p(np)%Tcol=0.0_WP                                                                                  !< Give zero collision force
                  this%lp%p(np)%d   =diam                                                                                     !< Assign diameter to account for full volume
                  this%lp%p(np)%pos =[this%cc%meta_structures_list(m)%x,this%cc%meta_structures_list(m)%y,this%cc%meta_structures_list(m)%z] !< Place the drop at the liquid barycenter
                  this%lp%p(np)%vel =[this%cc%meta_structures_list(m)%u,this%cc%meta_structures_list(m)%v,this%cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
                  this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])      !< Place the drop in the proper cell for the this%lp%cfg
                  this%lp%p(np)%flag=0                                                                                        !< Activate it
                  ! Increment particle counter
                  this%lp%np_=np
                  do n=1,2
                     ! Make room for new drop
                     np=this%lp%np_+1; call this%lp%resize(np)
                     ! Add the drop
                     this%lp%p(np)%id  =int(12,8)                                                                                   !< Give id (maybe based on break-up model?)
                     this%lp%p(np)%dt  =0.0_WP                                                                                     !< Let the drop find it own integration time
                     this%lp%p(np)%Acol=0.0_WP                                                                                     !< Give zero collision force
                     this%lp%p(np)%Tcol=0.0_WP                                                                                     !< Give zero collision force
                     this%lp%p(np)%d   =size_ratio*diam                                                                                       !< Assign diameter to account for full volume
                     this%lp%p(np)%pos =[this%cc%meta_structures_list(m)%x,this%cc%meta_structures_list(m)%y,this%cc%meta_structures_list(m)%z] &
                                       +sign(1.0_WP,real(n,WP)-1.5_WP)*2.0_WP*diam*this%cc%meta_structures_list(m)%axes(:,1)          
                     this%lp%p(np)%vel =[this%cc%meta_structures_list(m)%u,this%cc%meta_structures_list(m)%v,this%cc%meta_structures_list(m)%w] !< Assign mean structure velocity as drop velocity
                     this%lp%p(np)%ind =this%lp%cfg%get_ijk_global(this%lp%p(np)%pos,[this%lp%cfg%imin,this%lp%cfg%jmin,this%lp%cfg%kmin])      !< Place the drop in the proper cell for the this%lp%cfg
                     this%lp%p(np)%flag=0                                                                                          !< Activate it
                     ! Increment particle counter
                     this%lp%np_=np
                  end do
               end if

               ! Find local structs with matching id
               do n=this%cc%sync_offset+1,this%cc%sync_offset+this%cc%n_struct
                  if (this%cc%struct_list(this%cc%struct_map_(n))%parent.ne.this%cc%meta_structures_list(m)%id) cycle
                  ! Remove liquid in meta-structure cells
                  do l=1,this%cc%struct_list(this%cc%struct_map_(n))%nnode ! Loops over cells within local
                     i=this%cc%struct_list(this%cc%struct_map_(n))%node(1,l)
                     j=this%cc%struct_list(this%cc%struct_map_(n))%node(2,l)
                     k=this%cc%struct_list(this%cc%struct_map_(n))%node(3,l)
                     ! Remove liquid in that cell
                     this%my_converted_volume=this%my_converted_volume+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
                     this%vf%VF(i,j,k)=0.0_WP
                  end do
               end do
            end if
               
         end do ! n_meta_struct
      end if
      
      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()

      ! Clean up CCL
      call this%cc%deallocate_lists()

      ! Write the diameters and velocities
      filename='spray-all/droplets'
      do rank=0,this%cfg%nproc-1
         if (rank.eq.this%cfg%rank) then
            ! Open the file
            ! open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',position='append',iostat=ierr) ! need to write(iunit) without formatting spec
            open(newunit=iunit,file=trim(filename),form='formatted',status='old',access='stream',position='append',iostat=ierr)
            if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
            ! Output diameters and velocities
            do i=np_start+1,this%lp%np_
               write(iunit,*) this%lp%p(i)%d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),norm2([this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]),this%lp%p(i)%pos(1),this%lp%p(i)%pos(2),this%lp%p(i)%pos(3),'ligament',this%lp%p(i)%id  
            end do
            ! Close the file
            close(iunit)
         end if
         ! Force synchronization
         call MPI_BARRIER(this%cfg%comm,ierr)
      end do

   end subroutine breakup_ligament


   !> Generate a Gamma distribution for bag droplet formation
   !> where the number PDF is in the form
   !> p_n(x=d;alpha,beta)=x**(alpha-1)*exp(-x/beta)/beta**alpha/gamma(alpha)
   !> Adapted from Jackiw and Ashgriz 2022, JFM
   subroutine bag_droplet_gamma(this,h,R,alpha,beta)
      implicit none
      class(transfermodels), intent(inout) :: this
      real(WP), intent(in) :: h,R
      real(WP), intent(out) :: alpha,beta
      real(WP) :: d0,Utc,ac,b,lambda_rr,dr,ds,Oh
      real(WP) :: mean, stdev
      
      ! assert h,R != 0
      ! Drop diameter
      d0=1.0_WP
      ! Retraction speed
      Utc=sqrt(2.0_WP*this%fs%sigma/this%fs%rho_l/h)
      ! Centripetal acceleration
      ac=Utc**2/R
      ! Rim diameter
      b=sqrt(this%fs%sigma/this%fs%rho_l/ac)
      ! Receding rim wavelength
      lambda_rr=4.5_WP*b
      ! RP droplet diameter
      this%d_film_rp=1.89_WP*b
      ! Rim Ohnesorge number
      Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b*this%fs%sigma)
      ! Satellite droplet diameter
      ds=this%d_film_rp/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
      ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
      mean=0.25_WP*(h+b+this%d_film_rp+ds)/d0
      stdev=sqrt(0.25_WP*sum(([h,b,this%d_film_rp,ds]/d0-mean)**2))
      ! Gamma distribution parameters
      alpha=(mean/stdev)**2
      beta=stdev**2/mean

   end subroutine bag_droplet_gamma
    

   !> Setup spray statistics folders
   subroutine spray_statistics_setup(this)
      use c_interface, only: mkdir
      use messager, only: die
      implicit none
      class(transfermodels), intent(inout) :: this
      character(len=str_medium) :: filename
      integer :: iunit,ierr

      ! Create directory
      if (this%lp%cfg%amroot) then
         ! call execute_command_line('mkdir -p spray-all')
         call mkdir('spray-all')
         filename='spray-all/droplets'
         open(newunit=iunit,file=trim(filename),form='formatted',status='replace',access='stream',iostat=ierr)
         if (ierr.ne.0) call die('[transfermodel write spray stats] Could not open file: '//trim(filename))
         ! Write the header
         write(iunit,*) 'Diameter ','U ','V ','W ','Total velocity ','X ','Y ','Z ','origin '
         ! Close the file
         close(iunit)         
      end if
   end subroutine spray_statistics_setup


   !> Output spray statistics (velocity vs. diameter, mean/median diameters)
   subroutine spray_statistics(this)
      use mpi_f08,  only: MPI_REDUCE,MPI_SUM,MPI_BARRIER,MPI_GATHERV
      use messager, only: die
      use parallel, only: MPI_REAL_WP
      implicit none
      class(transfermodels), intent(inout) :: this
      real(WP) :: buf,d0,d2,d3
      integer :: rank,i,ierr
      integer, dimension(:), allocatable :: displacements
      real(WP), dimension(:), allocatable :: vollist_,vollist
      integer, parameter :: col_len=14
      
      ! Create safe np/d0
      d0=real(max(this%lp%np,1),WP)
      ! Initialize others
      d2=0.0_WP;d3=0.0_WP;this%mmd=0.0_WP

      ! Calculate diameter moments
      ! d3 is used as total droplet volume
      allocate(vollist_(1:this%lp%np_))
      d2=0.0_WP
      d3=0.0_WP
      do i=1,this%lp%np_
         d2=d2+this%lp%p(i)%d**2
         d3=d3+this%lp%p(i)%d**3
         vollist_(i)=this%lp%p(i)%d**3
      end do
      ! call MPI_ALLREDUCE(d2,buf,1,MPI_REAL_WP,MPI_SUM,this%lp%cfg%comm,ierr); d2=buf
      ! call MPI_ALLREDUCE(d3,buf,1,MPI_REAL_WP,MPI_SUM,this%lp%cfg%comm,ierr); d3=buf
      call MPI_REDUCE(d2,buf,1,MPI_REAL_WP,MPI_SUM,0,this%lp%cfg%comm,ierr); d2=buf
      call MPI_REDUCE(d3,buf,1,MPI_REAL_WP,MPI_SUM,0,this%lp%cfg%comm,ierr); d3=buf

      ! Gather droplet volumes
      ! if (this%lp%cfg%amroot) then
         allocate(vollist(1:this%lp%np))
         allocate(displacements(this%lp%cfg%nproc)); displacements=0
         do rank=1,this%lp%cfg%nproc-1
            displacements(rank+1)=displacements(rank)+this%lp%np_proc(rank)
         end do
      ! end if
      call MPI_GATHERV(vollist_,this%lp%np_,MPI_REAL_WP,vollist,this%lp%np_proc,displacements,MPI_REAL_WP,0,this%lp%cfg%comm,ierr)
      
      if (this%lp%cfg%amroot) then
         ! Sort volumes
         call quick_sort(vollist)
         ! Calculate d_30, d_32 (Sauter mean), and mass median diameters
         if (d2.le.0.0_WP) then
            this%d32=0.0_WP
            this%d30=0.0_WP
         else
            this%d32=d3/d2
            this%d30=(d3/d0)**(1.0_WP/3.0_WP)
         end if
         buf=0.0_WP
         volloop: do i=1,this%lp%np
            buf=buf+vollist(i)
            if (buf.ge.0.5_WP*d3) then
               this%mmd=vollist(i)**(1.0_WP/3.0_WP)
               exit volloop
            end if
         end do volloop
      end if
   contains
      ! Volume sorting
      recursive subroutine quick_sort(vol)
         implicit none
         real(WP), dimension(:)   :: vol
         integer :: imark
         if (size(vol).gt.1) then
            call quick_sort_partition(vol,imark)
            call quick_sort(vol(     :imark-1))
            call quick_sort(vol(imark:       ))
         end if
      end subroutine quick_sort
      subroutine quick_sort_partition(vol,marker)
         implicit none
         real(WP), dimension(  :) :: vol
         integer , intent(out)    :: marker
         integer :: ii,jj
         real(WP) :: dtmp,x
         x=vol(1)
         ii=0; jj=size(vol)+1
         do
            jj=jj-1
            do
               if (vol(jj).le.x) exit
               jj=jj-1
            end do
            ii=ii+1
            do
               if (vol(ii).ge.x) exit
               ii=ii+1
            end do
            if (ii.lt.jj) then
               dtmp =vol(  ii); vol(  ii)=vol(  jj); vol(  jj)=dtmp
            else if (ii.eq.jj) then
               marker=ii+1
               return
            else
               marker=ii
               return
            endif
         end do
      end subroutine quick_sort_partition
   end subroutine spray_statistics


   !> Find the minimum thickness in the largest volume film
   subroutine get_min_thickness(this,my_min_thickness,my_volume,id_largest_film)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_MINLOC
      use parallel,  only: MPI_REAL_WP,MPI_2REAL_WP
      use mathtools, only: Pi
      implicit none
      class(transfermodels), intent(inout) :: this
      real(WP), intent(out) :: my_min_thickness,my_volume
      integer, intent(out) :: id_largest_film
      integer  :: id,m,n,i,j,k,ierr
      real(WP), dimension(1:this%cc%cfg%nproc*this%cc%n_film_max) :: volume_,volume
      real(WP), dimension(1:2,1:this%cc%cfg%nproc*this%cc%n_film_max) :: min_thickness_,min_thickness_all

      if (this%cc%n_film_max.eq.0) return ! If there are no films globally
      volume_=0.0_WP
      min_thickness_(1,:)=huge(1.0_WP); min_thickness_(2,:)=this%cc%cfg%rank
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film ! Loops over film segments contained locally
         id=this%cc%film_list(this%cc%film_map_(m))%parent
         do n=1,this%cc%film_list(this%cc%film_map_(m))%nnode ! Loops over cells within local film segment
            i=this%cc%film_list(this%cc%film_map_(m))%node(1,n)
            j=this%cc%film_list(this%cc%film_map_(m))%node(2,n)
            k=this%cc%film_list(this%cc%film_map_(m))%node(3,n)
            volume_(id)=volume_(id)+this%vf%VF(i,j,k)*this%vf%cfg%vol(i,j,k)
            min_thickness_(1,id)=min(min_thickness_(1,id),this%vf%thickness(i,j,k))
         end do
      end do
      call MPI_ALLREDUCE(min_thickness_,min_thickness_all,this%cc%cfg%nproc*this%cc%n_film_max,MPI_2REAL_WP,MPI_MINLOC,this%cc%cfg%comm,ierr)
      do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film ! Loops over film segments contained locally
         id=this%cc%film_list(this%cc%film_map_(m))%parent
         this%cc%film_list(this%cc%film_map_(m))%min_thickness = min_thickness_all(1,id)
         this%cc%film_list(this%cc%film_map_(m))%minloc_rank   = int(min_thickness_all(2,id))
      end do
      call MPI_ALLREDUCE(volume_,volume,this%cfg%nproc*this%cc%n_film_max,MPI_REAL_WP,MPI_SUM,this%cc%cfg%comm,ierr)
      my_volume=maxval(volume)
      id_largest_film=maxloc(volume,DIM=1)
      my_min_thickness=min_thickness_all(1,id_largest_film)
      if (my_volume.lt.5.0e-3*4.0_WP/3.0_WP*Pi*0.5_WP**3) then
         my_volume=0.0_WP
         my_min_thickness=0.0_WP
      end if
   end subroutine get_min_thickness


end module transfermodel_class
