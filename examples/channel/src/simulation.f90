!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use incomp_class,      only: incomp
   ! use tensor_class,      only: tensor
   use fene_class,        only: fene
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single-phase incompressible flow solver, conformation tensor solvers and corresponding time tracker
   type(incomp),      public :: fs
   type(fene),        public :: fm
   ! type(tensor),      public :: ctquick
   ! type(tensor),      public :: ctupwind    
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile,velfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   ! real(WP), dimension(:,:,:,:), allocatable :: rhoU,rhoV,rhoW
   ! real(WP), dimension(:,:,:,:), allocatable :: resCT,resCTupwind,resCTquick
   ! real(WP), dimension(:,:,:,:), allocatable :: CT,CTstar,CTold
   
   !> Fluid viscosity
   real(WP) :: visc

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW

   !> Event for post-processing
   type(event) :: ppevt

   ! !> Conformation tensor max and min thresholds and magnitude
   ! real(WP) :: CTmag_max,CTmag_min
   
contains
   
   
   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_vel()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(Uavg (fs%cfg%jmin:fs%cfg%jmax)); Uavg =0.0_WP
      allocate(Uavg_(fs%cfg%jmin:fs%cfg%jmax)); Uavg_=0.0_WP
      allocate(vol_ (fs%cfg%jmin:fs%cfg%jmax)); vol_ =0.0_WP
      allocate(vol  (fs%cfg%jmin:fs%cfg%jmax)); vol  =0.0_WP
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j) = vol_(j)+fs%cfg%vol(i,j,k)
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
         filename='Uavg_'
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
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU       (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV       (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW       (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui         (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi         (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi         (  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR         (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(rhoU       (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(rhoV       (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(rhoW       (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(resCT      (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(resCTquick (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(resCTupwind(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(CT         (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(CTold      (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(CTstar     (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use ils_class, only: gmres_amg,pcg_pfmg
         use mathtools, only: twoPi
         integer :: i,j,k
         real(WP) :: amp,vel
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
         ! Initialize velocity based on specified bulk
         call param_read('Ubulk',Ubulk)
         call param_read('Wbulk',Wbulk)
         where (fs%umask.eq.0) fs%U=Ubulk
         where (fs%wmask.eq.0) fs%W=Wbulk
         meanU=Ubulk
         meanW=Wbulk
         ! To facilitate transition
         call param_read('Perturbation',amp)
         vel=sqrt(Ubulk**2+Wbulk**2)
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)
                  if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
               end do
            end do
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      ! Newtonian Solvers: pressure=gmres_amg,implicit=pcg_amg
      ! Newtonian Solvers: pressure=gmres_amg,implicit=pcg_smg

      ! Create a fene model solver
      create_fene_model_solver: block
         use ils_class,  only: pcg_pfmg
         use fene_class, only: bquick
         integer  :: i,j,k
         real(WP) :: Lmax,Wi
         ! Maximum extensibility of polymer chain
         call param_read('Maximum extension of polymer chain',Lmax)
         ! Weissenberg number for polymer
         call param_read('Weissenberg number',Wi)
         ! Create FENE model solver with bquick scheme
         fm=fene(cfg=cfg, scheme=bquick, extension=Lmax,weisnum=Wi, name= 'fene model solver')
         ! Configure implicit tensor solvers for QUICK scheme
         call param_read('FENE model iteration',fm%implicit%maxit)
         call param_read('FENE model tolerance',fm%implicit%rcvg)
         ! Assign constant density for tensor solver
         fm%rho=fs%rho
         ! Setup the solvers
         call fm%setup(implicit_ils=pcg_pfmg)
      end block create_fene_model_solver
      
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
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
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
         call cflfile%write()
         ! Create forcing monitor
         forcefile=monitor(fs%cfg%amRoot,'forcing')
         call forcefile%add_column(time%n,'Timestep number')
         call forcefile%add_column(time%t,'Time')
         call forcefile%add_column(meanU,'Bulk U')
         call forcefile%add_column(meanW,'Bulk W')
         call forcefile%write()
      end block create_monitor
      
      
      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Perform the output
         if (ppevt%occurs()) call postproc_vel()
      end block create_postproc
      
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! ! Model non-Newtonian fluid
         ! nonewt: block
         !    integer :: i,j,k
         !    real(WP) :: SRmag
         !    real(WP), parameter :: C=1.0e-2_WP
         !    real(WP), parameter :: n=0.3_WP
         !    ! Calculate SR
         !    call fs%get_strainrate(Ui,Vi,Wi,SR)
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

         ! ! Remember old tensor
         ! ctquick%TNold=ctquick%TN
         ! ctupwind%TNold=ctupwind%TN
         ! CTold=CT

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! ! Advance the conformation tensor using BQUICK
         ! advance_tensor: block
         !    integer :: i,j,k
         !    real(WP) :: CTmag
         !    ! Remember old conformation tensor
            
         !    ! ! Explicit calculation of drhoTNdt (based on QUICK scheme)
         !    ! call ctquick%get_drhoTNdt(drhoTNdtquick,rhoU,rhoV,rhoW)
         !    ! ! Explicit calculation of drhoTNdt (based on Upwind scheme)
         !    ! call ctupwind%get_drhoTNdt(drhoTNdtupwind,rhoU,rhoV,rhoW)
         !    ! Form implicit residuals for tensor (based on QUICK scheme)
         !    call ctquick%solve_implicit(time%dt,resCTquick,resU,resV,resW)
         !    ! Form implicit residuals for tensor (based on Upwind scheme)
         !    call ctupwind%solve_implicit(time%dt,resCTupwind,resU,resV,resW)
         !    ! Apply the residual to get predictor conformation tensor (based on QUICK scheme)
         !    CTstar=2.0_WP*ctquick%TN-CTold+resCTquick
         !    ! Check boundedness of conformation tensor and adjust CT calculation
         !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         !       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
         !          do i=fs%cfg%imino_,fs%cfg%imaxo_
         !             ! Calculate the magnitude of the conformation tensor
         !             CTmag=sqrt(CTstar(1,i,j,k)**2+CTstar(2,i,j,k)**2+CTstar(3,i,j,k)**2+2.0_WP*(CTstar(4,i,j,k)**2+CTstar(5,i,j,k)**2+CTstar(6,i,j,k)**2))
         !             if (CTmag.gt.CTmag_max.or.CTmag.lt.CTmag_min) then
         !                ! Advance CT in cell with 1st order upwind scheme
         !                CT(:,i,j,k)=2.0_WP*ctupwind%TN(:,i,j,k)-CTold(:,i,j,k)+resCTupwind(:,i,j,k)
         !                resCT(:,i,j,k)=resCTupwind(:,i,j,k)
         !             else 
         !                ! Keep CT value in cell based on QUICK scheme
         !                CT(:,i,j,k)=CTstar(:,i,j,k)
         !                resCT(:,i,j,k)=resCTquick(:,i,j,k)
         !             end if 
         !          end do
         !       end do
         !    end do

         !    call fs%cfg%sync(CT)
            
         !    !Assign current CT to arrays stored in solvers
         !    ctquick%TN =CT
         !    ctupwind%TN=CT
         !    call fs%cfg%sync( ctquick%TN)
         !    call fs%cfg%sync(ctupwind%TN)
            
         !    ! Assign current resCT to parameters used in solvers
         !    resCTquick =resCT
         !    resCTupwind=resCT
         !    call fs%cfg%sync( resCTquick)
         !    call fs%cfg%sync(resCTupwind)
            
         ! end block advance_tensor

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)


            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Add body forcing
            forcing: block
               use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
               use parallel, only: MPI_REAL_WP
               integer :: i,j,k,ierr
               real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
               myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) then
                           myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                           myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                        end if
                        if (fs%wmask(i,j,k).eq.0) then
                           myW   =myW   +fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))
                           myWvol=myWvol+fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)
                        end if
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               where (fs%umask.eq.0) resU=resU+Ubulk-meanU
               call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanW=meanW/Wvol
               where (fs%wmask.eq.0) resW=resW+Wbulk-meanW
            end block forcing   

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
         call mfile%write()
         call cflfile%write()
         call forcefile%write()

         ! Specialized post-processing
         if (ppevt%occurs()) call postproc_vel()

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
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR)
      
   end subroutine simulation_final
   
end module simulation
