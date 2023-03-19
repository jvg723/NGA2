!> Various definitions and tools for running an NGA2 simulation
module simulation
	use precision,         only: WP
	use geometry,          only: cfg
	use hypre_str_class,   only: hypre_str
	use tpns_class,        only: tpns
	use vfs_class,         only: vfs
	! use fene_class,        only: fene
	use timetracker_class, only: timetracker
	use ensight_class,     only: ensight
   	use surfmesh_class,    only: surfmesh
	use event_class,       only: event
	use monitor_class,     only: monitor
	implicit none
	private
	
	!> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
	type(hypre_str),   public :: ps
	type(hypre_str),   public :: vs
	type(tpns),        public :: fs
	type(vfs),         public :: vf
	! type(fene),        public :: fm
	type(timetracker), public :: time
	
	!> Ensight postprocessing
   	type(surfmesh) :: smesh
	type(ensight)  :: ens_out
	type(event)    :: ens_evt
	
	!> Simulation monitor file
	type(monitor) :: mfile,cflfile
	
	public :: simulation_init,simulation_run,simulation_final
	
	!> Private work arrays
	real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
	real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
	real(WP), dimension(:,:,:,:),   allocatable :: resSC,SC_
    real(WP), dimension(:,:,:,:,:), allocatable :: gradu
	
	!> Problem definition
	real(WP), dimension(3) :: center
	real(WP) :: radius
	
	 !> Fluid viscosity (solvent,polymer,total)
	real(WP) :: visc_s,visc_l,visc_p,visc_g

	!> Artifical diffusivity for conformation tensor
	real(WP) :: stress_diff

	!> FENE-P model parameters
	real(WP) :: Lmax,lambda,Beta

contains


   !> Function that defines a level set function for a rising bubble
	function levelset_rising_bubble(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G
		! Create the bubble
	   	G=radius-sqrt(sum((xyz-center)**2))
	end function levelset_rising_bubble

	!> Function that localizes the top (x+) of the domain
	function xp_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (i.eq.pg%imax+1) isIn=.true.
	 end function xp_locator
	 
  
	!> Function that localizes the bottom (x-) of the domain boundary
	function xm_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (i.le.pg%imin) isIn=.true.
	end function xm_locator
  
	!> Function that localizes the top (y+) of the domain
	function yp_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		implicit none
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (j.eq.pg%jmax+1) isIn=.true.
	end function yp_locator
	 
	 
	 !> Function that localizes the bottom (y-) of the domain
	 function ym_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		implicit none
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (j.eq.pg%jmin) isIn=.true.
	 end function ym_locator
	 
	 
	 !> Function that localizes the top (z+) of the domain
	 function zp_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		implicit none
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (k.eq.pg%kmax+1) isIn=.true.
	 end function zp_locator
	 
	 
	 !> Function that localizes the bottom (z-) of the domain
	 function zm_locator(pg,i,j,k) result(isIn)
		use pgrid_class, only: pgrid
		implicit none
		class(pgrid), intent(in) :: pg
		integer, intent(in) :: i,j,k
		logical :: isIn
		isIn=.false.
		if (k.eq.pg%kmin) isIn=.true.
	 end function zm_locator
  
	
	
	!> Initialization of problem solver
	subroutine simulation_init
		use param, only: param_read
		implicit none
		
		
		! Allocate work arrays
	   	allocate_work_arrays: block
		   	allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(gradu(3,3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			! Scalar solver
			allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6))
			allocate(SC_  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6)) !< Temp SC array for checking bquick bound
		end block allocate_work_arrays
		
		
		! Initialize time tracker with 2 subiterations
	   initialize_timetracker: block
		   	time=timetracker(amRoot=cfg%amRoot)
			call param_read('Max timestep size',time%dtmax)
			call param_read('Max cfl number',time%cflmax)
			call param_read('Max time',time%tmax)
			time%dt=time%dtmax
			time%itmax=2
		end block initialize_timetracker
		
		
		! Initialize our VOF solver and field
	   	create_and_initialize_vof: block
			use mms_geom, only: cube_refine_vol
			use vfs_class, only: lvira,VFhi,VFlo
			integer :: i,j,k,n,si,sj,sk
			real(WP), dimension(3,8) :: cube_vertex
			real(WP), dimension(3) :: v_cent,a_cent
			real(WP) :: diam
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
			! Create a VOF solver
			vf=vfs(cfg=cfg,reconstruction_method=lvira,name='VOF')
			! Initialize to a bubble
			call param_read('Diameter',diam)
			center=[0.0_WP,0.0075_WP,0.0_WP]
			radius=diam/2.0_WP
			do k=vf%cfg%kmino_,vf%cfg%kmaxo_
				do j=vf%cfg%jmino_,vf%cfg%jmaxo_
					do i=vf%cfg%imino_,vf%cfg%imaxo_
						! Set cube vertices
				      n=0
						do sk=0,1
							do sj=0,1
								do si=0,1
									n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
								end do
							end do
						end do
						! Call adaptive refinement code to get volume and barycenters recursively
				      	vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
						call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_rising_bubble,0.0_WP,amr_ref_lvl)
						vf%VF(i,j,k)=1.0_WP-(vol/vf%cfg%vol(i,j,k))
						if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
							vf%Lbary(:,i,j,k)=v_cent
							vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
						else
							vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
							vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
						end if
					end do
				end do
			end do
			! Update the band
		   	call vf%update_band()
			! Perform interface reconstruction from VOF field
		   	call vf%build_interface()
			! Set interface planes at the boundaries
         	call vf%set_full_bcond()
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
	   		use tpns_class, only: clipped_neumann,slip,dirichlet
		   	use hypre_str_class, only: pcg_pfmg
			! Create flow solver
			fs=tpns(cfg=cfg,name='Two-phase NS')
			! Assign constant viscosity to each phase
			call param_read('Liquid dynamic viscosity',visc_l); fs%visc_l=visc_l
			call param_read('Gas dynamic viscosity',   visc_g); fs%visc_g=visc_g
			! Assign constant density to each phase
		   	call param_read('Liquid density',fs%rho_l)
			call param_read('Gas density',fs%rho_g)
			! Assign acceleration of gravity
			call param_read('Gravity',fs%gravity)
			! Read in surface tension coefficient
			call param_read('Surface tension coefficient',fs%sigma)
			! Clipped Neumann outflow on the top and bottom of the domain
			call fs%add_bcond(name='bc_yp',type=slip,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
			call fs%add_bcond(name='bc_ym',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
			! ! Clip Neuman on sides
			! call fs%add_bcond(name='bc_xp',type=slip,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
			! call fs%add_bcond(name='bc_xm',type=slip,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
			! call fs%add_bcond(name='bc_zp',type=slip,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
			! call fs%add_bcond(name='bc_zm',type=slip,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
			! Configure pressure solver
			ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
			call param_read('Pressure iteration',ps%maxit)
			call param_read('Pressure tolerance',ps%rcvg)
			! Configure implicit velocity solver
			vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
			call param_read('Implicit iteration',vs%maxit)
			call param_read('Implicit tolerance',vs%rcvg)
			! Setup the solver
			call fs%setup(pressure_solver=ps,implicit_solver=vs)
			! Zero initial field
			fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
			! Calculate cell-centered velocities and divergence
			call fs%interp_vel(Ui,Vi,Wi)
			call fs%get_div()
		end block create_and_initialize_flow_solver

		! ! Create a FENE model 
		! create_fene: block
		! 	use ils_class,  only: gmres
		! 	use multiscalar_class, only: bquick
		! 	! Create FENE model solver
		! 	fm=fene(cfg=cfg,scheme=bquick,name='FENE model')
		! 	! Assign aritifical stress diffusivisty
		! 	call param_read('Stress diffusivisty',stress_diff)
		! 	fm%diff=stress_diff
		! 	! Assign constant density
		! 	fm%rho=1.00_WP
		! 	! Maximum extensibility of polymer chain
		! 	call param_read('Maximum extension of polymer chain',Lmax)
		! 	! Relaxation time for polymer
		! 	call param_read('Polymer Relaxation Time',lambda)
		! 	! Solvent/polymer viscosity ratio
		! 	call param_read('Beta',Beta)
		! 	! Polymer viscosity
		! 	visc_p=visc_s*((1.00_WP-Beta)/Beta)
		! 	! Configure implicit scalar solver
		! 	fm%implicit%maxit=fs%implicit%maxit; fm%implicit%rcvg=fs%implicit%rcvg
		! 	! Setup the solver
		! 	call fm%setup(implicit_ils=gmres) 
		! 	! Intalize conformation tensor to identity matrix
		! 	fm%SC(:,:,:,1)=1.00_WP !Cxx
		! 	fm%SC(:,:,:,2)=0.00_WP !Cyx
		! 	fm%SC(:,:,:,3)=0.00_WP !Czx
		! 	fm%SC(:,:,:,4)=1.00_WP !Cyy
		! 	fm%SC(:,:,:,5)=0.00_WP !Czy
		! 	fm%SC(:,:,:,6)=1.00_WP !Czz
		! 	! Initalize stress tensor
		! 	call fm%get_stressTensor(lambda,Lmax,visc_p)
	 	! end block create_fene
	   
	   
	   ! Add Ensight output
	   create_ensight: block
		   ! Create Ensight output from cfg
		   ens_out=ensight(cfg=cfg,name='RisingBubble')
			! Create event for Ensight output
		   ens_evt=event(time=time,name='Ensight output')
		   call param_read('Ensight output period',ens_evt%tper)
		   ! Add variables to output
		   call ens_out%add_vector('velocity',Ui,Vi,Wi)
		   call ens_out%add_scalar('VOF',vf%VF)
		   call ens_out%add_scalar('pressure',fs%P)
		   call ens_out%add_scalar('curvature',vf%curv)
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
	   
	   
	end subroutine simulation_init
	
	
	!> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
	subroutine simulation_run
      use tpns_class, only: static_contact
		implicit none
		
		! Perform time integration
	   	do while (.not.time%done())
			
			! Increment time
			call fs%get_cfl(time%dt,time%cfl)
			call time%adjust_dt()
			call time%increment()
			
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

				! ! Build mid-time scalar
				! fm%SC=0.5_WP*(fm%SC+fm%SCold)

				! Form velocity gradient
				call fs%get_gradu(gradu)

				! ! ============= SCALAR SOLVER =======================  
            	! ! Reset interpolation metrics to QUICK scheme
				! call fm%metric_reset()

				! ! Calculate explicit SC prior to checking bounds
				! pre_check: block
				!    ! Explicit calculation of resSC from scalar equation
				!    call fm%get_drhoSCdt(resSC,fs%U,fs%V,fs%W)
				!    ! Assemble explicit residual
				!    resSC=-2.0_WP*(fm%rho*fm%SC-fm%rho*fm%SCold)+time%dt*resSC
				!    ! Apply this residual
				!    SC_=2.0_WP*fm%SC-fm%SCold+resSC
				! end block pre_check
				
				! ! Check boundedess of explicit SC calculation
				! call fm%metric_modification(SC=SC_,SCmin=0.0_WP)
	
				! ! Calculate explicit SC post checking bounds
				! post_check: block
				!    ! Explicit calculation of resSC from scalar equation
				!    call fm%get_drhoSCdt(resSC,fs%U,fs%V,fs%W)
				!    ! Assemble explicit residual
				!    resSC=-2.0_WP*(fm%rho*fm%SC-fm%rho*fm%SCold)+time%dt*resSC
				! end block post_check
				
				! ! Add FENE source terms
				! fene: block
				!    ! Calculate CgradU terms
				!    call fm%get_CgradU(gradu)    
				!    ! Calculate T terms
				!    call fm%get_stressTensor(lambda,Lmax,visc_p)     
				!    ! Add source terms to calculated residual
				!    resSC=resSC+(fm%CgradU-(1.00_WP/visc_p)*fm%T)*time%dt
				! end block fene
	
				! ! Form implicit residual
				! call fm%solve_implicit(time%dt,resSC,fs%U,fs%V,fs%W)
	
				! ! Apply this residual
				! fm%SC=2.0_WP*fm%SC-fm%SCold+resSC
	
				! ! Apply other boundary conditions on the resulting field
				! call fm%apply_bcond(time%t,time%dt)
				! ! ===================================================
				
				! ============= VELOCITY SOLVER ======================
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

				! ! Add in polymer stress
				! polymer: block
				! 	integer :: i,j,k
				! 	real (WP) :: H
				! 	! Calculate updated elastic tensor terms
				! 	call fm%get_stressTensor(lambda,Lmax,visc_p)
				! 	! Get its divergence 
				! 	call fm%get_divT(fs) 
				! 	! Add divT to momentum equation for G.ge.0
				! 	H=0.0_WP
				! 	do k=fs%cfg%kmin_,fs%cfg%kmax_
                !   		do j=fs%cfg%jmin_,fs%cfg%jmax_
                !     		do i=fs%cfg%imin_,fs%cfg%imax_
				! 				! Use level set array to set Heavside value
				! 				if (vf%G(i,j,k).ge.0.0_WP) then
				! 					H=1.0_WP
				! 				else
				! 					H=0.0_WP
				! 				end  if 
				! 				if (fs%umask(i,j,k).eq.0) resU=resU+H*fm%divT(i,j,k,1)*time%dt !> x face/U velocity
				! 				if (fs%vmask(i,j,k).eq.0) resV=resV+H*fm%divT(i,j,k,2)*time%dt !> y face/V velocity
				! 				if (fs%wmask(i,j,k).eq.0) resW=resW+H*fm%divT(i,j,k,3)*time%dt !> z face/W velocity
				! 			end do
				! 		end do
				! 	end do
			 	! end block polymer
				
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
				call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
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
			call vf%get_max()
			call mfile%write()
			call cflfile%write()
			
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