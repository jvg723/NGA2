!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,Lx,Ly,Lz,diam
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single two-phase flow solver and volume fraction solver and corresponding time tracker
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(timetracker), public :: time
   type(surfmesh),    public :: smesh                                              
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: SRmag
   real(WP), dimension(:,:,:,:), allocatable :: SR
   
   !> Problem definition
   real(WP) :: Rel,Weg,r_visc,r_rho
   real(WP) :: n,lambda,visc_0,visc_inf,visc_l,visc_g

   !> Post-processing
   type(event) :: ppevt

   !> Liquid axial flow velocity at inlet
   real(WP) :: U0
   
contains
   
   
    !> Function that localizes the center of an annulus located at x=z=y=0
   function annulus(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j)**2.0_WP+pg%zm(k)**2.0_WP.ge.(diam/4.0_WP)**2.0_WP.and.pg%ym(j)**2.0_WP+pg%zm(k)**2.0_WP.lt.(diam/2.0_WP)**2.0_WP) isIn=.true.
   end function annulus

   !> Function that localizes the center of an annulus located at x=z=y=0 for bulk axial flow
   function bulk_axial(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j)**2.0_WP+pg%zm(k)**2.0_WP.lt.(diam/2.0_WP)**2.0_WP) isIn=.true.
   end function bulk_axial

    !> Function that localizes the leftmost domain boundary around the annulus
   function left_boundary_outflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin.and.pg%ym(j)**2.0_WP+pg%zm(k)**2.0_WP.gt.((diam/2.0_WP)**2.0_WP)+(diam/2.0_WP)) isIn=.true.
   end function left_boundary_outflow

   !> Function that localizes the rightmost domain boundary
   function right_boundary_outflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_boundary_outflow
   
   
  !> Specialized subroutine that outputs the vertical liquid distribution
   subroutine postproc_data()
      ! use mathtools, only: Pi
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: myVOF,VOF
      real(WP), dimension(:), allocatable :: myVEL,VEL
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(myVOF(vf%cfg%jmin:vf%cfg%jmax)); myVOF=0.0_WP
      allocate(myVEL(vf%cfg%jmin:vf%cfg%jmax)); myVEL=0.0_WP
      allocate(  VOF(vf%cfg%jmin:vf%cfg%jmax)); VOF=0.0_WP
      allocate(  VEL(vf%cfg%jmin:vf%cfg%jmax)); VEL=0.0_WP
      ! Initialize local data to zero
      myVOF=0.0_WP; myVEL=0.0_WP
      ! Integrate all data over x and z
      do k=vf%cfg%kmin_,vf%cfg%kmax_
         do j=vf%cfg%jmin_,vf%cfg%jmax_
            do i=vf%cfg%imin_,vf%cfg%imax_
               myVOF(j)=myVOF(j)+vf%VF(i,j,k)
               myVEL(j)=myVEL(j)+fs%U(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE(myVOF,VOF,vf%cfg%ny,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr); VOF=VOF/real(vf%cfg%nx*vf%cfg%nz,WP)
      call MPI_ALLREDUCE(myVEL,VEL,vf%cfg%ny,MPI_REAL_WP,MPI_SUM,vf%cfg%comm,ierr); VEL=VEL/real(vf%cfg%nx*vf%cfg%nz,WP)
      ! If root, print it out
      if (vf%cfg%amRoot) then
         ! call execute_command_line('mkdir -p stats')
         ! filename='profile_'
         filename='./stats/profile_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         ! open(newunit=iunit,file='stats/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12)') 'Height','VOF','VEL'
         do j=vf%cfg%jmin,vf%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5)') vf%cfg%ym(j),VOF(j),VEL(j)
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(myVOF,VOF)
      deallocate(myVEL,VEL)
   end subroutine postproc_data
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SRmag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR   (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='pressure_swirl')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,r2p,VFhi,VFlo
         use mathtools, only: twoPi
         use random,    only: random_uniform
         use parallel,  only: MPI_REAL_WP
         use mpi_f08
         integer :: i,j,k
         ! Create a VOF solver
         vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
         ! Initialize to flat interface in the annulus
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  if (i.eq.vf%cfg%imin.and.vf%cfg%ym(j)**2.0_WP+vf%cfg%zm(k)**2.0_WP.ge.(diam/4.0_WP)**2.0_WP.and.vf%cfg%ym(j)**2.0_WP+vf%cfg%zm(k)**2.0_WP.lt.(diam/2.0_WP)**2.0_WP) then
                     vf%VF(i,j,k)=1.0_WP
                  else
                     vf%VF(i,j,k)=0.0_WP
                  end if
                  vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use tpns_class, only: dirichlet,clipped_neumann,neumann
         use ils_class, only: pcg_pfmg,gmres_amg
         use mathtools, only: Pi
         use random,    only: random_uniform
         integer :: i,j,k
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! ! Read in flow conditions
         ! call param_read('Liquid Reynolds number',Rel); fs%visc_l=1.0_WP/(Rel+epsilon(Rel))
         ! call param_read('Viscosity ratio',r_visc); fs%visc_g=fs%visc_l/r_visc
         ! call param_read('Density ratio',r_rho); fs%rho_l=1.0_WP; fs%rho_g=fs%rho_l/r_rho
         ! call param_read('Gas Weber number',Weg); fs%sigma=1.0_WP/(Weg+epsilon(Weg))
         ! ! Read in power law constant for non-newtonian liquid
         ! call param_read('Power law constant',n)
         ! ! Set low and high SR viscosity values
         ! visc_0=r_visc*visc_g
         ! visc_inf=0.0_WP
         ! ! Characteristic time scale
         ! lambda=r_visc
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',visc_l); fs%visc_l=visc_l
         call param_read('Gas dynamic viscosity'   ,visc_g); fs%visc_g=visc_g
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density'   ,fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         ! Inflow on the left
         call fs%add_bcond(name='inflow', type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=annulus)
         call fs%add_bcond(name='bulk', type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=bulk_axial)
         ! Neumann on the left around annulus
         call fs%add_bcond(name='lhs_neumann', type=clipped_neumann,face='x',dir=-1,canCorrect=.false.,locator=left_boundary_outflow)
         ! Outflow on the right
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=right_boundary_outflow)
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_pfmg,implicit_ils=pcg_pfmg)
      end block create_and_initialize_flow_solver

      ! Initialize our velocity field
      initialize_velocity: block
         use tpns_class, only: bcond
         use random,     only: random_uniform
         use mathtools,  only: twoPi
         type(bcond), pointer :: mybc
         real(WP)             :: omega,a,theta,rad,Utheta,Utheta_max
         integer  :: n,i,j,k
         ! Zero initial field in the domain
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Rankine vortex parameters for inflow
         call param_read('Axial velocity',U0)
         a=diam/4.0_WP
         Utheta_max=1.2_WP*U0
         omega=Utheta_max/a
         call fs%get_bcond('inflow',mybc)
         ! Apply swirl component dirichlet at inlet
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            ! Calculate tangential velocity along annulus radius
            rad=0.0_WP;theta=0.0_WP
            rad=sqrt(fs%cfg%zm(k)**2.0_WP+fs%cfg%ym(j)**2.0_WP)
            theta=atan2(fs%cfg%zm(k),fs%cfg%ym(j))
            if (rad.le.a) then
               Utheta=Utheta_max
            elseif (rad.gt.a) then
               Utheta=(omega*a**2.0_WP)/rad
            end if
            ! Put in tangential velocity into cartesian components
            fs%V(i,j,k)= Utheta*cos(theta)
            fs%W(i,j,k)=-Utheta*sin(theta)
         end do
         call fs%get_bcond('bulk',mybc)
         ! Apply bulk component dirichlet at inlet
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            ! Bulk injection velocity
            fs%U(i,j,k)=U0
         end do
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr()
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
      end block initialize_velocity

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface
         integer :: i,j,k,nplane,np
         ! Include an extra variable for number of planes
         smesh=surfmesh(nvar=1,name='plic')
         smesh%varname(1)='nplane'
         ! Transfer polygons to smesh
         call vf%update_surfmesh(smesh)
         ! Also populate nplane variable
         smesh%var(1,:)=1.0_WP
         np=0
         do k=vf%cfg%kmin_,vf%cfg%kmax_
            do j=vf%cfg%jmin_,vf%cfg%jmax_
               do i=vf%cfg%imin_,vf%cfg%imax_
                  do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pressure_swirl')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_scalar('SRmag',SRmag)
         call ens_out%add_scalar('visc_l',fs%visc_l)
         call ens_out%add_surface('vofplic',smesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
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
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(vf%SDint,'SD integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      

      ! Create specialized post-processing
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Create directory to write to
         if (cfg%amRoot) call execute_command_line('mkdir -p stats')
         ! Perform the output
         if (ppevt%occurs()) call postproc_data()
      end block create_postproc


   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      ! use tpns_class, only: static_contact
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         ! Calculate SR
         call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)

         ! ! Model non-Newtonian fluid
         ! nonewt: block
         !    integer :: i,j,k
         !    ! real(WP), parameter :: C=1.137e-3_WP
         !    ! Update viscosity
         !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
         !       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
         !          do i=fs%cfg%imino_,fs%cfg%imaxo_
         !             ! Power law model 
         !             ! SRmag(i,j,k)=sqrt(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
         !             ! SRmag(i,j,k)=max(SRmag(i,j,k),1000.0_WP**(1.0_WP/(n-1.0_WP)))
         !             ! fs%visc_l(i,j,k)=C*SRmag(i,j,k)**(n-1.0_WP)
         !             ! Carreau Model
         !             SRmag(i,j,k)=sqrt(2.00_WP*SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
         !             fs%visc_l(i,j,k)=visc_inf+(visc_0-visc_inf)*(1.00_WP+(lambda*SRmag(i,j,k))**2.00_WP)**((n-1.00_WP)/2.00_WP)
         !          end do
         !       end do
         !    end do
         !    ! call fs%cfg%sync(fs%visc_l)
         ! end block nonewt
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         
         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) then 
            ! Update surfmesh object
            update_smesh: block
            use irl_fortran_interface
            integer :: nplane,np,i,j,k
            ! Transfer polygons to smesh
            call vf%update_surfmesh(smesh)
            ! Also populate nplane variable
            smesh%var(1,:)=1.0_WP
            np=0
            do k=vf%cfg%kmin_,vf%cfg%kmax_
               do j=vf%cfg%jmin_,vf%cfg%jmax_
                  do i=vf%cfg%imin_,vf%cfg%imax_
                     do nplane=1,getNumberOfPlanes(vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; smesh%var(1,np)=real(getNumberOfPlanes(vf%liquid_gas_interface(i,j,k)),WP)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         ! Perform ensight output 
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call mfile%write()
         call cflfile%write()

         ! Specialized post-processing
         if (ppevt%occurs()) call postproc_data()
         
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final
   
   
end module simulation
