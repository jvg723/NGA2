!> Transfer model class:
!> Provides support for conversion of subgrid films and ligaments to Lagrangian droplets

module transfermodel_class
   use precision,      only: WP
   use config_class,   only: config
   use ccl_class,      only: ccl
   use lpt_class,      only: lpt
   use vfs_class,      only: vfs
   use tpns_class,     only: tpns
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
      real(WP) :: min_filmthickness      =1.0e-4_WP
      real(WP) :: diam_over_filmthickness=1.0e+1_WP
      real(WP) :: max_eccentricity       =5.0e-1_WP
      real(WP) :: d_threshold            =1.0e-1_WP

   contains

      procedure :: initialize
      procedure :: transfer_vf_to_drops
      procedure, private :: bag_droplet_gamma


   end type transfermodels

contains


   !> Initialize
   subroutine initialize(this,cfg,vf,fs,lp)
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
         integer :: i,j,k
         ! Create the CCL object
         this%cc=ccl(cfg=this%cfg,name='CCL')
         this%cc%max_interface_planes=2
         this%cc%VFlo=VFlo
         this%cc%dot_threshold=-0.5_WP
         
         ! Perform CCL step
         call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%cc%deallocate_lists()
      end block create_and_initialize_ccl

   end subroutine initialize

   !> Transfer vf to drops
   subroutine transfer_vf_to_drops(this)
      use vfs_class, only: r2p
      use messager, only: die
      use string, only: str_medium
      implicit none
      class(transfermodels), intent(inout) :: this
      character(len=str_medium) :: filename
      integer :: iunit,ierr,rank,id_largest_film

      ! Perform a first pass with simplest CCL
      call this%cc%build_lists(VF=this%vf%VF,U=this%fs%U,V=this%fs%V,W=this%fs%W)

      ! Loop through identified detached structs and remove those that are spherical enough
      remove_struct: block
         use mathtools, only: pi
         integer :: m,n,l,i,j,k,np
         real(WP) :: lmin,lmax,eccentricity,diam

         ! Loops over film segments contained locally
         do m=1,this%cc%n_meta_struct

            ! Test if sphericity is compatible with transfer
            lmin=this%cc%meta_structures_list(m)%lengths(3)
            if (lmin.eq.0.0_WP) lmin=this%cc%meta_structures_list(m)%lengths(2) ! Handle 2D case
            lmax=this%cc%meta_structures_list(m)%lengths(1)
            eccentricity=sqrt(1.0_WP-lmin**2/(lmax**2+tiny(1.0_WP)))
            if (eccentricity.gt.this%max_eccentricity) cycle

            ! Test if diameter is compatible with transfer
            diam=(6.0_WP*this%cc%meta_structures_list(m)%vol/pi)**(1.0_WP/3.0_WP)
            if ((diam.eq.0.0_WP).or.(diam.gt.this%d_threshold)) cycle

            
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
                  this%vf%VF(i,j,k)=0.0_WP
               end do
            end do
            
         end do
         
      end block remove_struct

      ! Sync VF and clean up IRL and band
      call this%vf%cfg%sync(this%vf%VF)
      call this%vf%clean_irl_and_band()
      
      ! Clean up CCL
      call this%cc%deallocate_lists()
      
      if (this%vf%reconstruction_method.eq.r2p) then
         ! Perform more detailed CCL in a second pass
         this%cc%max_interface_planes=2
         call this%cc%build_lists(VF=this%vf%VF,poly=this%vf%interface_polygon,U=this%fs%U,V=this%fs%V,W=this%fs%W)
         call this%vf%sense_interface()
         call this%cc%get_min_thickness()
         call this%cc%sort_by_thickness()

         ! Loop through identified films and remove those that are thin enough
         remove_film: block
            use mathtools, only: pi,normalize,cross_product
            use random,    only: random_uniform
            use myrandom,  only: random_gamma
            use parallel,  only: MPI_REAL_WP
            use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_LOGICAL,MPI_LOR
            use irl_fortran_interface
            integer :: m,n,l,i,j,k,np,ip,np_old,ii,jj,kk
            real(WP) :: Vt,Vl,Hl,Vd,diam,d0
            real(WP), dimension(3) :: pref,nref,tref,sref
            real(WP) :: theta
            logical :: has_burst
            real(WP), dimension(2) :: s_bary,c_bary
            integer :: nodes_transferred
            ! Conversion model
            real(WP) :: alpha, beta, curv_sum, ncurv
            real(WP) :: max_film_droplet_diam=4.0e-4_WP

            ! Loops over film segments contained locally
            d0=1.0_WP ! Characteristic diameter

            do m=this%cc%film_sync_offset+1,this%cc%film_sync_offset+this%cc%n_film
            
               ! Skip non-liquid films
               if (this%cc%film_list(this%cc%film_map_(m))%phase.ne.1) cycle
            
               ! Skip films that are still thick enough
               if (this%cc%film_list(this%cc%film_map_(m))%min_thickness.gt.this%min_filmthickness) cycle

               ! We are still here: transfer the film to drops
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
               call this%bag_droplet_gamma(this%cc%film_thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
               Vt=0.0_WP      ! Transferred volume
               Vl=0.0_WP      ! We will keep track incrementally of the liquid volume to transfer to ensure conservation
               np_old=this%lp%np_  ! Remember old number of particles
               Vd=pi/6.0_WP*(random_gamma(alpha,.true.)*beta*d0)**3
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
                  call this%bag_droplet_gamma(this%cc%film_thickness(i,j,k),2.0_WP*ncurv/curv_sum,alpha,beta)
                  
                  ! Create drops from available liquid volume
                  do while (Vl-Vd.gt.0.0_WP)
                     ! Make room for new drop
                     np=this%lp%np_+1; call this%lp%resize(np)
                     ! Add the drop
                     this%lp%p(np)%id  =int(1,8)                                   !< Give id (maybe based on break-up model?)
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
                     Vd=pi/6.0_WP*(random_gamma(alpha,.true.)*beta*d0)**3
                  end do

                  ! Remove liquid in that cell
                  this%vf%VF(i,j,k)=0.0_WP

               end do
               
               ! Based on how many particles were created, decide what to do with left-over volume
               if (Vt.eq.0.0_WP) then ! No particle was created, we need one...
                  ! Add one last drop for remaining liquid volume
                  np=this%lp%np_+1; call this%lp%resize(np)
                  ! Add the drop
                  this%lp%p(np)%id  =int(1,8)                                   !< Give id (maybe based on break-up model?)
                  this%lp%p(np)%dt  =0.0_WP                                     !< Let the drop find it own integration time
                  this%lp%p(np)%Acol =0.0_WP                                    !< Give zero collision force
                  this%lp%p(np)%Tcol =0.0_WP                                    !< Give zero collision force
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

            end do ! local films
            
            ! Sync VF and clean up IRL and band
            call this%vf%cfg%sync(this%vf%VF)
            call this%vf%clean_irl_and_band()                  

            ! Clean up CCL
            call this%cc%deallocate_lists()

         end block remove_film

      end if ! r2p
     
      ! Resync the spray
      call this%lp%sync()

   end subroutine transfer_vf_to_drops
   

   !> Generate a Gamma distribution for bag droplet formation
   !> where the number PDF is in the form
   !> p_n(x=d;alpha,beta)=x**(alpha-1)*exp(-x/beta)/beta**alpha/gamma(alpha)
   !> Adapted from Jackiw and Ashgriz 2022, JFM
   subroutine bag_droplet_gamma(this,h,R,alpha,beta)
      use myrandom, only: random_gamma
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
      dr=1.89_WP*b
      ! Rim Ohnesorge number
      Oh=this%fs%visc_l/sqrt(this%fs%rho_l*b**3*this%fs%sigma)
      ! Satellite droplet diameter
      ds=dr/sqrt(2.0_WP+3.0_WP*Oh/sqrt(2.0_WP))
      ! Mean and standard deviation of diameter of all modes, normalized by drop diameter
      mean=0.25_WP*(h+b+dr+ds)/d0
      stdev=sqrt(0.25_WP*sum(([h,b,dr,ds]/d0-mean)**2))
      ! Gamma distribution parameters
      alpha=(mean/stdev)**2
      beta=stdev**2/mean

   end subroutine bag_droplet_gamma
    

end module transfermodel_class
