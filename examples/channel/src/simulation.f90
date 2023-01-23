!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp
   use fene_class,        only: fene
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single-phase incompressible flow solver, fene model,sgs model and corresponding time tracker
   type(incomp),      public :: fs
   type(fene),        public :: fm
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs         
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:),   allocatable :: SR,resSC,SC_
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu

   
   !> Fluid viscosity (solvent,polymer,total)
   real(WP) :: visc_s,visc_p,visc_0

   !> Artifical diffusivity for conformation tensor
   real(WP) :: stress_diff

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW
   real(WP) :: H,Q,px

   !> Event for post-processing
   type(event) :: ppevt    
   
   !> FENE-P model parameters
   real(WP) :: Lmax,lambda,Beta

   !> CFL numbers
   real(WP) :: cflc,cflp

contains
   
   
   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_vel()
      use string,      only: str_medium
      use mpi_f08,     only: MPI_ALLREDUCE,MPI_SUM
      use parallel,    only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(Uavg    (fs%cfg%jmin:fs%cfg%jmax)); Uavg =0.0_WP
      allocate(Uavg_   (fs%cfg%jmin:fs%cfg%jmax)); Uavg_=0.0_WP
      allocate(vol_    (fs%cfg%jmin:fs%cfg%jmax)); vol_ =0.0_WP
      allocate(vol     (fs%cfg%jmin:fs%cfg%jmax)); vol  =0.0_WP
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j) =vol_(j) +fs%cfg%vol(i,j,k)
               Uavg_(j)=Uavg_(j)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE( vol_, vol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      do j=fs%cfg%jmin,fs%cfg%jmax
         if (vol(j).gt.0.0_WP) then
            Uavg(j)=Uavg(j)/vol(j)
         else
            Uavg(j)=0.0_WP
         end if
      end do
      ! If root, print it out
      if (fs%cfg%amRoot) then
         filename='./velocity/Uavg_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12)') 'Height','Uavg'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5)') fs%cfg%ym(j),Uavg(j)
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(Uavg,Uavg_,vol,vol_)
   end subroutine postproc_vel

   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_ct()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: Cxxavg,Cxxavg_,Cxyavg,Cxyavg_,Cyyavg,Cyyavg_
      real(WP), dimension(:), allocatable :: Txxavg,Txxavg_,Txyavg,Txyavg_,Tyyavg,Tyyavg_
      real(WP), dimension(:), allocatable :: vol,vol_
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(Cxxavg (fs%cfg%jmin:fs%cfg%jmax)); Cxxavg =0.0_WP
      allocate(Cxxavg_(fs%cfg%jmin:fs%cfg%jmax)); Cxxavg_=0.0_WP
      allocate(Cxyavg (fs%cfg%jmin:fs%cfg%jmax)); Cxyavg =0.0_WP
      allocate(Cxyavg_(fs%cfg%jmin:fs%cfg%jmax)); Cxyavg_=0.0_WP
      allocate(Cyyavg (fs%cfg%jmin:fs%cfg%jmax)); Cyyavg =0.0_WP
      allocate(Cyyavg_(fs%cfg%jmin:fs%cfg%jmax)); Cyyavg_=0.0_WP
      allocate(Txxavg (fs%cfg%jmin:fs%cfg%jmax)); Txxavg =0.0_WP
      allocate(Txxavg_(fs%cfg%jmin:fs%cfg%jmax)); Txxavg_=0.0_WP
      allocate(Txyavg (fs%cfg%jmin:fs%cfg%jmax)); Txyavg =0.0_WP
      allocate(Txyavg_(fs%cfg%jmin:fs%cfg%jmax)); Txyavg_=0.0_WP
      allocate(Tyyavg (fs%cfg%jmin:fs%cfg%jmax)); Tyyavg =0.0_WP
      allocate(Tyyavg_(fs%cfg%jmin:fs%cfg%jmax)); Tyyavg_=0.0_WP
      allocate(vol_   (fs%cfg%jmin:fs%cfg%jmax)); vol_   =0.0_WP
      allocate(vol    (fs%cfg%jmin:fs%cfg%jmax)); vol    =0.0_WP
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j)   =vol_(j)+fs%cfg%vol(i,j,k)
               Cxxavg_(j)=Cxxavg_(j)+fs%cfg%vol(i,j,k)*fm%SC(i,j,k,1)
               Cxyavg_(j)=Cxyavg_(j)+fs%cfg%vol(i,j,k)*fm%SC(i,j,k,2)
               Cyyavg_(j)=Cyyavg_(j)+fs%cfg%vol(i,j,k)*fm%SC(i,j,k,4)
               Txxavg_(j)=Txxavg_(j)+fs%cfg%vol(i,j,k)*fm%T(i,j,k,1)
               Txyavg_(j)=Txyavg_(j)+fs%cfg%vol(i,j,k)*fm%T(i,j,k,2)
               Tyyavg_(j)=Tyyavg_(j)+fs%cfg%vol(i,j,k)*fm%T(i,j,k,4)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(vol_   ,vol   ,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Cxxavg_,Cxxavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Cxyavg_,Cxyavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Cyyavg_,Cyyavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Txxavg_,Txxavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Txyavg_,Txyavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Tyyavg_,Tyyavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      do j=fs%cfg%jmin,fs%cfg%jmax
         if (vol(j).gt.0.0_WP) then
            Cxxavg(j)=Cxxavg(j)/vol(j)
            Cxyavg(j)=Cxyavg(j)/vol(j)
            Cyyavg(j)=Cyyavg(j)/vol(j)
            Txxavg(j)=Txxavg(j)/vol(j)
            Txyavg(j)=Txyavg(j)/vol(j)
            Tyyavg(j)=Tyyavg(j)/vol(j)
         else
            Cxxavg(j)=0.0_WP
            Cxyavg(j)=0.0_WP
            Cyyavg(j)=0.0_WP
            Txxavg(j)=0.0_WP
            Txyavg(j)=0.0_WP
            Tyyavg(j)=0.0_WP
         end if
      end do
      ! If root, print it out
      if (fs%cfg%amRoot) then
         filename='./stress/Cavg_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'Height','Cxx_avg','Cxy_avg','Cyy_avg','Txx_avg','Txy_avg','Tyy_avg'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') fs%cfg%ym(j),Cxxavg(j),Cxyavg(j),Cyyavg(j),Txxavg(j),Txyavg(j),Tyyavg(j)
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(Cxxavg,Cxxavg_,Cxyavg,Cxyavg_,Cyyavg,Cyyavg_,Txxavg,Txxavg_,Txyavg,Txyavg_,Tyyavg,Tyyavg_,vol,vol_)
   end subroutine postproc_ct

   !> Specialized subroutine to plot simulation vs theory curves
   subroutine plotter()
      use string,      only: str_medium
      use mpi_f08,     only: MPI_ALLREDUCE,MPI_SUM
      use parallel,    only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr
      character(len=str_medium) :: timestamp 
      character(len=str_medium) :: vel_file
      character(len=str_medium) :: strs_file
      character(len=str_medium), parameter :: plt_file='~/Builds/NGA2/examples/channel/plot.gp'    
      ! Plot comparison from root processor
      if (fs%cfg%amRoot) then
         ! Time step
         write(timestamp,'(es12.5)') time%t
         ! Velocity
         vel_file=trim('Uavg_')//trim(adjustl(timestamp))
         ! Stress
         strs_file=trim('Cavg_')//trim(adjustl(timestamp))
         ! Store timestep and array naming for reading in gnuplot
         open(newunit=iunit,file='./plots/gp_input',form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,5x,a12,5x,a12,5x,a12,5x,a12,5x,a16,5x,a16)') 'timestep','beta','Uc','H','visc_0','vel_file','strs_file'
         write(iunit,'(es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,es12.5,5x,a16,5x,a16)') time%t,beta,Ubulk,fs%cfg%yL,visc_s+visc_p,vel_file,strs_file
         close(iunit)
         ! Plot the curves using gnuplot
         call execute_command_line('gnuplot ' // plt_file)
      end if
   end subroutine plotter

   !> Calculate F+ values for FENE-P theory velocity curve
   function fpvalues(x,A,C) result(Fplus)
      real(WP), intent(in) :: x,A,C
      real(WP)             :: Fplus
      ! F^+(x)
      Fplus=(C*x+sqrt(A**3.00_WP+(C*x)**2.00_WP))**(1.00_WP/3.00_WP)
   end function fpvalues

   !> Calculate F- values for FENE-P theory velocity curve
   function fmvalues(x,A,C) result(Fminus)
      real(WP), intent(in) :: x,A,C
      real(WP)             :: Fminus
      ! F^-(x)
      Fminus=-abs(C*x-sqrt(A**3.00_WP+(C*x)**2.00_WP))**(1.00_WP/3.00_WP)
   end function fmvalues

   !> Calculate G+ values for FENE-P theory velocity curve
   function gpvalues(x,A,C) result(Gplus)
      real(WP), intent(in) :: x,A,C
      real(WP)             :: Gplus
      ! G^+(x)
      Gplus=3.00_WP*C*x+sqrt(A**3.00_WP+(C*x)**2.00_WP)
   end function gpvalues

   !> Calculate G- values for FENE-P theory velocity curve
   function gmvalues(x,A,C) result(Gminus)
      real(WP), intent(in) :: x,A,C
      real(WP)             :: Gminus
      ! G^-(x)
      Gminus=3.00_WP*C*x-sqrt(A**3.00_WP+(C*x)**2.00_WP)
   end function gmvalues        
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(resU (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (    cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR   (6,  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradu(3,3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6))
         allocate(SC_  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6)) !< Temp SC array for checking bquick bound
      end block allocate_work_arrays
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Channel Flow')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use ils_class, only: gmres_amg,pcg_pfmg,pcg_amg,gmres_pfmg,pfmg,smg,pcg_smg
         use mathtools, only: twoPi
         integer :: i,j,k
         real(WP) :: amp,vel,omega
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Solvent dynamic viscosity',visc_s); fs%visc=visc_s
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=gmres_amg,implicit_ils=smg)      !> 3D case
         ! call fs%setup(pressure_ils=pcg_smg,implicit_ils=smg)      !> 2D Case
         ! Initialize velocity based on specified bulk
         call param_read('Ubulk',Ubulk)
         call param_read('Wbulk',Wbulk)
         where (fs%umask.eq.0) fs%U=Ubulk
         where (fs%wmask.eq.0) fs%W=Wbulk
         meanU=Ubulk
         meanW=Wbulk
         ! To facilitate transition
         call param_read('Perturbation',amp)
         call param_read('Angular velocity',omega)
         vel=sqrt(Ubulk**2+Wbulk**2)
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)
                  if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  ! if (fs%vmask(i,j,k).eq.0) fs%V(i,j,k)=fs%V(i,j,k)+2.00_WP*twoPi*omega*fs%cfg%zm(k)*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  ! if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)-2.00_WP*twoPi*omega*fs%cfg%ym(j)*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      ! Newtonian Solvers: pressure=gmres_amg,implicit=gmres_amg
      ! Non-Newtonian Solvers: pressure=gmres_amg,implicit=pcg_smg

      ! Create a FENE model 
      create_fene: block
         use ils_class,  only: gmres_pilut,bbox,pcg_bbox,gmres_amg
         use multiscalar_class, only: bquick
         integer :: j
         ! Create FENE model solver
         fm=fene(cfg=cfg,scheme=bquick,name='FENE model')
         ! Assign aritifical stress diffusivisty
         call param_read('Stress diffusivisty',stress_diff)
         fm%diff=stress_diff
         ! Assign constant density
         fm%rho=fs%rho
         ! Maximum extensibility of polymer chain
         call param_read('Maximum extension of polymer chain',Lmax)
         ! Relaxation time for polymer
         call param_read('Polymer Relaxation Time',lambda)
         ! Solvent/polymer viscosity ratio
         call param_read('Beta',Beta)
         ! Polymer viscosity
         visc_p=visc_s*((1.00_WP-Beta)/Beta)
         ! Configure implicit scalar solver
         fm%implicit%maxit=fs%implicit%maxit; fm%implicit%rcvg=fs%implicit%rcvg
         ! Setup the solver
         call fm%setup(implicit_ils=gmres_pilut) 
      end block create_fene

      ! Pressure gradient to drive periodic flow
      pressure_grad: block
         ! Total dynamic viscosity (solvent+polymer)
         visc_0=visc_s+visc_p  
         ! Channel height
         H=fs%cfg%yL
         ! Flow rate through the channel
         Q=Ubulk*H 
         ! Contant pressure gradient
         px=-12.00_WP*(visc_0/H**3.00_WP)*Q 
      end block pressure_grad

      ! Initalize FEBNE model tensors
      init_fene: block
         ! real(WP) :: H,Lc,Uc,Q,epsilon,visc,px,A,B,C   !> Terms for theoretical stress tensor calculation in laminar flow
         ! integer :: i,j,k,iunit,ierr
         ! use mpi_f08,  only: MPI_BCAST
         ! use parallel, only: comm,amRoot,MPI_REAL_WP 
         ! real(WP) :: omega,phi,F                  !> Terms for mean shear flow assumption
         ! real(WP), dimension(128) :: dudy=0.0_WP  !> Turbulent velocity gradient for initializing C tensor
         
         ! !> Initalize C and T tensor for laminar channel flow using analytical solution from D.O.A. Cruz et al. (2005)
         ! ! Set conformation tensor to 0
         ! fm%SC=0.0_WP
         ! ! Characteristic scales
         ! H=fs%cfg%yL
         ! Lc=fs%cfg%yL
         ! Uc=Ubulk
         ! ! Flow rate through the channel
         ! Q=Uc*H 
         ! ! Polymer terms
         ! epsilon=1.00_WP-(3.00_WP/Lmax**2.00_WP)
         ! ! Total dynamic viscosity (solvent+polymer)
         ! visc=visc_s+visc_p                                            
         ! ! Pressure gradient in channel
         ! px=-12.00_WP*(visc/H**3.00_WP)*Q 
         ! ! Constant coefficents 
         ! A=(visc_p**2.00_WP/(6.00_WP*epsilon*lambda**2.00_WP))*(1.00_WP+visc_p/visc_s)
         ! C=(visc_p**2.00_WP/(4.00_WP*epsilon*lambda**2.00_WP))*(visc_p/visc_s)*px 
         ! ! For 2D flow set stress yz,xz,yy and zz terms to 0
         ! fm%T(:,:,:,3)=0.00_WP   !Txz
         ! fm%T(:,:,:,4)=0.00_WP   !Tyy
         ! fm%T(:,:,:,5)=0.00_WP   !Tyz
         ! fm%T(:,:,:,6)=0.00_WP   !Tzz
         ! ! Calculate shear and normal stress
         ! do k=fs%cfg%kmin_,fs%cfg%kmax_
         !    do j=fs%cfg%jmin_,fs%cfg%jmax_
         !       do i=fs%cfg%imin_,fs%cfg%imax_
         !          ! Position dependent coefficent being solved for in cubic equation
         !          B=C*fs%cfg%ym(j)
         !          ! Position dependent shear and normmal stress
         !          fm%T(i,j,k,2)=(B+sqrt(A**3.00_WP+B**2.00_WP))**(1.00_WP/3.00_WP)-abs(B-sqrt(A**3.00_WP+B**2.00_WP))**(1.00_WP/3.00_WP) !Txy
         !          fm%T(i,j,k,1)=2.00_WP*(lambda/visc_p)*fm%T(i,j,k,2)**2.00_WP                                                                       !Txx
         !       end do
         !    end do
         ! end do

         ! !> Initalize C and T tensor for turbulent channel flow using procedure of R. Sureshkummar et. al (1997)
         ! ! Read inital dudy from root processor
         ! if (amRoot) then
         !    open(newunit=iunit,file='./initial_conditions/velgrad.txt',form='formatted',status='old',access='stream',iostat=ierr)  
         !       do j=1,128   
         !          read(iunit,'(es12.5)') dudy(j) 
         !       end do 
         !    close(iunit)
         ! end if
         ! ! Then the root broadcasts
         ! call MPI_BCAST(dudy,128,MPI_REAL_WP,0,comm,ierr)
         ! ! Calculate conformation tensor
         ! do j=fs%cfg%jmino_,fs%cfg%jmaxo_
         !    ! Analytical solution parameters
         !    omega=((sqrt(2.00_WP)*Wei)/Lmax)*dudy(j)
         !    phi=log((3.00_WP*sqrt(3.00_WP)*omega/2.00_WP)+sqrt(1.00_WP+(3.00_WP*sqrt(3.00_WP)*omega/2.00_WP)**2.00_WP))
         !    F=(sqrt(3.00_WP)*omega)/(2.00_WP*sinh((phi+epsilon(phi))/3.00_WP))
         !    ! Prescribe analytical solution
         !    fm%SC(:,j,:,1)=(1.00_WP/(F+epsilon(F)))*(1.00_WP+(2.00_WP*Wei**2.00_WP/(F+epsilon(F))**2.00_WP)*(dudy(j)**2.00_WP)) ! Cxx
         !    fm%SC(:,j,:,2)=(Wei/(F+epsilon(F))**2.00_WP)*dudy(j)                                                                ! Cxy
         !    fm%SC(:,j,:,3)=0.00_WP                                                                                              ! Cxz
         !    fm%SC(:,j,:,4)=1.00_WP/(F+epsilon(F))                                                                               ! Cyy
         !    fm%SC(:,j,:,5)=0.00_WP                                                                                              ! Cyz
         !    fm%SC(:,j,:,6)=1.00_WP/(F+epsilon(F))                                                                               ! Czz
         ! end do  
         ! ! Calculate stress tensor
         ! call fm%get_stressTensor(fm%SC,Wei,Lmax) 
         
         ! > Initalize C to 0
         fm%SC=0.0_WP
         call fm%get_stressTensor(lambda,Lmax,visc_p)

      end block init_fene
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('visc_t',sgs%visc)
         call ens_out%add_scalar('Cxx',fm%SC(:,:,:,1))
         call ens_out%add_scalar('Cxy',fm%SC(:,:,:,2))
         call ens_out%add_scalar('Czx',fm%SC(:,:,:,3))
         call ens_out%add_scalar('Cyy',fm%SC(:,:,:,4))
         call ens_out%add_scalar('Czy',fm%SC(:,:,:,5))
         call ens_out%add_scalar('Czz',fm%SC(:,:,:,6))
         call ens_out%add_scalar('Txx',fm%T (:,:,:,1))
         call ens_out%add_scalar('Txy',fm%T (:,:,:,2))
         call ens_out%add_scalar('Tzx',fm%T (:,:,:,3))
         call ens_out%add_scalar('Tyy',fm%T (:,:,:,4))
         call ens_out%add_scalar('Tzy',fm%T (:,:,:,5))
         call ens_out%add_scalar('Tzz',fm%T (:,:,:,6))
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,cflc)
         call fm%get_cfl(time%dt,fs%rho,lambda,visc_p,cflp)
         time%cfl=cflc
         call fs%get_max()
         call fm%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%add_column(fm%CFLp_x,'Polymer Stress xCFL')
         call cflfile%add_column(fm%CFLp_y,'Polymer Stress yCFL')
         call cflfile%add_column(fm%CFLp_z,'Polymer Stress zCFL')
         call cflfile%write()
         ! Create forcing monitor
         forcefile=monitor(fs%cfg%amRoot,'forcing')
         call forcefile%add_column(time%n,'Timestep number')
         call forcefile%add_column(time%t,'Time')
         call forcefile%add_column(meanU,'Bulk U')
         call forcefile%add_column(meanW,'Bulk W')
         call forcefile%write()
      end block create_monitor

      ! Theory solution for 2D FENE-P channel flow from D.O.A. Cruz et al. (2005)
      theory: block
         use string, only: str_medium
         real(WP) :: b2,eps,A,B,C,lam                        !> Terms for theoretical stress tensor calculation in laminar flow
         real(WP) :: Fp_H,Fm_H,Gp_H,Gm_H,Fp_y,Fm_y,Gp_y,Gm_y !> Terms for velocity curves
         integer :: i,j,k,iunit,ierr
         real(WP), dimension(:), allocatable :: Txy,Txx,u
         ! Allocate vertical line storage
         allocate(Txy(fs%cfg%jmin:fs%cfg%jmax)); Txy=0.0_WP
         allocate(Txx(fs%cfg%jmin:fs%cfg%jmax)); Txx=0.0_WP
         allocate(u  (fs%cfg%jmin:fs%cfg%jmax)); u  =0.0_WP
         ! Polymer terms
         b2=Lmax**2.00_WP-3.00_WP
         eps=1.00_WP/((b2+5.00_WP))
         lam=((b2+2.00_WP)/(b2+5.00_WP))*lambda
         ! Constant coefficents 
         A=(visc_p**2.00_WP/(6.00_WP*eps*lam**2.00_WP))*(1.00_WP+visc_p/visc_s)
         C=(visc_p**2.00_WP/(4.00_WP*eps*lam**2.00_WP))*(visc_p/visc_s)*px 
         ! Loop over channel height to calculate velocity and stress
            do j=fs%cfg%jmin,fs%cfg%jmax
               ! Position dependent coefficent being solved for in cubic equation
               B=C*fs%cfg%ym(j)
               ! Shear stress
               Txy(j)=(B+sqrt(A**3.00_WP+B**2.00_WP))**(1.00_WP/3.00_WP)-abs(B-sqrt(A**3.00_WP+B**2.00_WP))**(1.00_WP/3.00_WP) 
               ! Normal stress
               Txx(j)=2.00_WP*(lam/visc_p)*Txy(j)**2.00_WP
               !> Velocity curve parameters
               if (fs%cfg%ym(j).lt.0.00_WP) then
                  ! @ -H/2
                  Fp_H=fpvalues(-H/2.00_WP,A,C)
                  Fm_H=fmvalues(-H/2.00_WP,A,C)
                  Gp_H=gpvalues(-H/2.00_WP,A,C)
                  Gm_H=gmvalues(-H/2.00_WP,A,C)
               else if (fs%cfg%ym(j).ge.0.00_WP) then
                  ! @ H/2
                  Fp_H=fpvalues(H/2.00_WP,A,C)
                  Fm_H=fmvalues(H/2.00_WP,A,C)
                  Gp_H=gpvalues(H/2.00_WP,A,C)
                  Gm_H=gmvalues(H/2.00_WP,A,C)
               end if
               ! @ y(j)
               Fp_y=fpvalues(fs%cfg%ym(j),A,C)
               Fm_y=fmvalues(fs%cfg%ym(j),A,C)
               Gp_y=gpvalues(fs%cfg%ym(j),A,C)
               Gm_y=gmvalues(fs%cfg%ym(j),A,C)
               ! Velocity
               if (Beta.eq.1.00_WP) then
                  u(j)=(-(px*(H/2.00_WP)**2.00_WP)/(2.00_WP*visc_s))*(1.00_WP-(fs%cfg%ym(j)/(H/2.00_WP))**2.00_WP)
               else 
                  u(j)=(-(px*(H/2.00_WP)**2.00_WP)/(2.00_WP*visc_s))*(1.00_WP-(fs%cfg%ym(j)/(H/2.00_WP))**2.00_WP)+(3.00_WP/(8.00_WP*C*visc_s))*(Fp_H*Gm_H-Fp_y*Gm_y+Fm_H*Gp_H-Fm_y*Gp_y)                                                           
               end if
            end do
            ! Create directory to write to
            if (cfg%amRoot) call execute_command_line('mkdir -p plots')
            ! If root, print it out
            if (fs%cfg%amRoot) then
               open(newunit=iunit,file='./plots/theory',form='formatted',status='replace',access='stream',iostat=ierr)
               write(iunit,'(a12,3x,a12,3x,a12,3x,a12)') 'Height','u','Txx','Txy'
               do j=fs%cfg%jmin,fs%cfg%jmax
                  write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') fs%cfg%ym(j),u(j),Txx(j),Txy(j)
               end do
               close(iunit)
            end if 
         ! Deallocate arrays
         deallocate(Txy,Txx,u)
      end block theory

      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Create directory to write to
         if (cfg%amRoot) call execute_command_line('mkdir -p velocity')
         if (cfg%amRoot) call execute_command_line('mkdir -p stress')
         ! Perform the output
         if (ppevt%occurs()) call postproc_vel()
         if (ppevt%occurs()) call postproc_ct()
         if (ppevt%occurs()) call plotter()
      end block create_postproc

   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Calcualte CFL numbers
         call fs%get_cfl(time%dt,cflc)                        !> Connvective CFL
         call fm%get_cfl(time%dt,fs%rho,lambda,visc_p,cflp)   !> Polymer stress CFL
         time%cfl=cflc
         
         ! Increment time
         call time%adjust_dt()
         call time%increment()
         
         ! ! Model non-Newtonian fluid
         ! nonewt: block
         !    integer :: i,j,k
         !    real(WP) :: SRmag
         !    real(WP), parameter :: C=1.0e-2_WP
         !    real(WP), parameter :: n=0.3_WP
         !    ! Calculate SR
         !    call fs%get_strainrate(SR)
         !    ! Update viscosity
         !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         !       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
         !          do i=fs%cfg%imino_,fs%cfg%imaxo_
         !             SRmag=sqrt(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
         !             SRmag=max(SRmag,1000.0_WP**(1.0_WP/(n-1.0_WP)))
         !             fs%visc(i,j,k)=C*SRmag**(n-1.0_WP)
         !          end do
         !       end do
         !    end do
         !    call fs%cfg%sync(fs%visc)
         ! end block nonewt

         ! Remember old scalars
         fm%SCold=fm%SC

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W     

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Build mid-time scalar
            fm%SC=0.5_WP*(fm%SC+fm%SCold)

            ! Form velocity gradient
            call fs%get_gradu(gradu)

            ! ============= SCALAR SOLVER =======================  
            ! Reset interpolation metrics to QUICK scheme
            call fm%metric_reset()

            ! Calculate explicit SC prior to checking bounds
            pre_check: block
               ! Explicit calculation of resSC from scalar equation
               call fm%get_drhoSCdt(resSC,fs%U,fs%V,fs%W)
               ! Assemble explicit residual
               resSC=-2.0_WP*(fm%rho*fm%SC-fm%rho*fm%SCold)+time%dt*resSC
               ! Apply this residual
               SC_=2.0_WP*fm%SC-fm%SCold+resSC
            end block pre_check
            
            ! Check boundedess of explicit SC calculation
            call fm%metric_modification(SC=SC_,SCmin=0.0_WP)

            ! Calculate explicit SC post checking bounds
            post_check: block
               ! Explicit calculation of resSC from scalar equation
               call fm%get_drhoSCdt(resSC,fs%U,fs%V,fs%W)
               ! Assemble explicit residual
               resSC=-2.0_WP*(fm%rho*fm%SC-fm%rho*fm%SCold)+time%dt*resSC
            end block post_check
            
            ! Add FENE source terms
            fene: block
               ! Calculate CgradU terms
               call fm%get_CgradU(gradu)    
               ! Calculate T terms
               call fm%get_stressTensor(lambda,Lmax,visc_p)     
               ! Add source terms to calculated residual
               resSC=resSC+(fm%CgradU-(1.00_WP/visc_p)*fm%T)*time%dt
            end block fene

            ! Form implicit residual
            call fm%solve_implicit(time%dt,resSC,fs%U,fs%V,fs%W)

            ! Apply this residual
            fm%SC=2.0_WP*fm%SC-fm%SCold+resSC

            ! Apply other boundary conditions on the resulting field
            call fm%apply_bcond(time%t,time%dt)
            ! ===================================================

            ! ============= VELOCITY SOLVER ======================  
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Add body forcing
            forcing: block
               ! ! Enforce constant flow rate
               ! use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
               ! use parallel, only: MPI_REAL_WP
               ! integer :: i,j,k,ierr
               ! real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
               ! myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
               ! do k=fs%cfg%kmin_,fs%cfg%kmax_
               !    do j=fs%cfg%jmin_,fs%cfg%jmax_
               !       do i=fs%cfg%imin_,fs%cfg%imax_
               !          if (fs%umask(i,j,k).eq.0) then
               !             myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
               !             myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
               !          end if
               !          if (fs%wmask(i,j,k).eq.0) then
               !             myW   =myW   +fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))
               !             myWvol=myWvol+fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)
               !          end if
               !       end do
               !    end do
               ! end do
               ! call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               ! call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               ! where (fs%umask.eq.0) resU=resU+Ubulk-meanU
               ! call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               ! call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanW=meanW/Wvol
               ! where (fs%wmask.eq.0) resW=resW+Wbulk-meanW
               ! Enforce constant pressure gradient
               where (fs%umask.eq.0) resU=resU+(-px*time%dt)
            end block forcing

            ! Add in polymer stress
            polymer: block
               ! Calculate updated elastic tensor terms
               call fm%get_stressTensor(lambda,Lmax,visc_p)
               ! Get its divergence 
               call fm%get_divT(fs) 
               ! Add divT to momentum equation 
               where (fs%umask.eq.0) resU=resU+fm%divT(:,:,:,1)*time%dt !> x face/U velocity
               where (fs%vmask.eq.0) resV=resV+fm%divT(:,:,:,2)*time%dt !> y face/V velocity
               where (fs%wmask.eq.0) resW=resW+fm%divT(:,:,:,3)*time%dt !> z face/W velocity
            end block polymer

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            ! ====================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call fm%get_max()
         call mfile%write()
         call cflfile%write()
         call forcefile%write()

         ! Specialized post-processing
         if (ppevt%occurs()) call postproc_vel()
         if (ppevt%occurs()) call postproc_ct()
         if (ppevt%occurs()) call plotter()

      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu,resSC,SC_)
      
   end subroutine simulation_final
   
end module simulation
