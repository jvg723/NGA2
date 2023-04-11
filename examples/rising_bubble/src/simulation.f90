!> Various definitions and tools for running an NGA2 simulation
module simulation
	use precision,         only: WP
	use geometry,          only: cfg
	!use hypre_str_class,   only: hypre_str
	use hypre_uns_class,   only: hypre_uns
	use ddadi_class,       only: ddadi
	use tpns_class,        only: tpns
	use vfs_class,         only: vfs
	use fene_class,        only: fene
	use timetracker_class, only: timetracker
	use ensight_class,     only: ensight
   	use surfmesh_class,    only: surfmesh
	use event_class,       only: event
	use monitor_class,     only: monitor
	implicit none
	private
	
	!> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
	!type(hypre_str),   public :: ps
   	type(hypre_uns),   public :: ps
	type(ddadi),       public :: vs,ss
	type(tpns),        public :: fs
	type(vfs),         public :: vf
	type(fene),        public :: fm
	type(timetracker), public :: time
	
	!> Ensight postprocessing
   	type(surfmesh) :: smesh
	type(ensight)  :: ens_out
	type(event)    :: ens_evt
	
	!> Simulation monitor file
	type(monitor) :: mfile,cflfile
	
	public :: simulation_init,simulation_run,simulation_final
	
	!> Private work arrays
	real(WP), dimension(:,:,:), 	allocatable :: resU,resV,resW
	real(WP), dimension(:,:,:), 	allocatable :: Ui,Vi,Wi
	real(WP), dimension(:,:,:,:),   allocatable :: resSC,SC_,trSC_
    real(WP), dimension(:,:,:,:,:), allocatable :: gradu
	real(WP), dimension(:,:,:,:),   allocatable :: fR,CgradU
	real(WP), dimension(:,:,:),   	allocatable :: SRmag
	real(WP), dimension(:,:,:),   	allocatable :: visc_p,visc_norm
   	real(WP), dimension(:,:,:,:), 	allocatable :: SR
	
	!> Problem definition
	real(WP), dimension(3) :: center
	real(WP) :: radius

	!> Fluid viscosity (solvent,polymer inital,liquid total,gas)
	real(WP) :: visc_s,visc_p0,visc_l,visc_g

	!> Shear thinning parameters
	real(WP) :: n,alpha,visc_0,visc_inf

	!> Artifical diffusivity for conformation tensor
	real(WP) :: stress_diff

	!> FENE-P model parameters
	real(WP) :: Lmax,lambda
	
contains


   !> Function that defines a level set function for a rising bubble problem
	function levelset_rising_bubble(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G
		! Create the bubble
	   	G=-radius+sqrt(sum((xyz-center)**2))
	end function levelset_rising_bubble
	
	
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
			allocate(gradu(3,3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(SRmag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        	allocate(SR   (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			! Scalar solver
			allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6))
			allocate(SC_   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6)) !< Temp SC array for checking bquick bound
			allocate(trSC_ (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1)) !< Temp trSC array for checking bquick bound
			allocate(fR    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6)) !< Array to hold relaxation function for FENE
            allocate(CgradU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,6)) !< Sum of distortion terms (CdotU and (CdotU)T)
			! Rate depdent polymer viscosity
			allocate(visc_p(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			! Array for normalized viscosity
			allocate(visc_norm(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
			use vfs_class, only: elvira,VFhi,VFlo
			integer :: i,j,k,n,si,sj,sk
			real(WP), dimension(3,8) :: cube_vertex
			real(WP), dimension(3) :: v_cent,a_cent
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
			! Create a VOF solver
			vf=vfs(cfg=cfg,reconstruction_method=elvira,name='VOF')
			! Initialize a bubble
			center=[0.0_WP,0.01_WP,0.0_WP]
			radius=0.005_WP
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
						vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
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
		   	!use hypre_str_class, only: pcg_pfmg,pcg_smg
        	use hypre_uns_class, only: gmres_amg,pcg_amg
         	! Create flow solver
			fs=tpns(cfg=cfg,name='Two-phase NS')
			! Assign constant viscosity to each phase
			call param_read('Solvent dynamic viscosity',visc_s) !; fs%visc_l=visc_l
			call param_read('Gas dynamic viscosity',visc_g); fs%visc_g=visc_g
			! Assign constant density to each phase
		   	call param_read('Liquid density',fs%rho_l)
			call param_read('Gas density',fs%rho_g)
			! Read in surface tension coefficient
			call param_read('Surface tension coefficient',fs%sigma)
			! Assign acceleration of gravity
		   	call param_read('Gravity',fs%gravity)
			! Configure pressure solver
			!ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_smg,nst=7)
			!ps%maxlevel=16
			ps=hypre_uns(cfg=cfg,name='Pressure',method=pcg_amg,nst=7)
			call param_read('Pressure iteration',ps%maxit)
			call param_read('Pressure tolerance',ps%rcvg)
			! Configure implicit velocity solver
			vs=ddadi(cfg=cfg,name='Velocity',nst=7)
			! Setup the solver
			call fs%setup(pressure_solver=ps,implicit_solver=vs)
			! Zero initial field
			fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
			! Calculate cell-centered velocities and divergence
			call fs%interp_vel(Ui,Vi,Wi)
			call fs%get_div()
	   	end block create_and_initialize_flow_solver

	   	! Set inital viscosity of liquid phase
		solvent_viscosity: block
			! Carreau model parameters	
			call param_read('Power law constant',n)
			call param_read('Shear rate parameter',alpha)
			!> Carreau model
			! ! No polymer viscosity
			! visc_p=0.0_WP
			! ! No infintie shear rate viscosity
			! visc_inf=0.0_WP
			! ! Zero shear rate viscosity
			! visc_0=visc_s; fs%visc_l=visc_0
			! visc_norm=fs%visc_l/visc_0
			!> Hyribd model
			! Zero shear rate viscosity
			call param_read('Zero shear rate viscosity',visc_0)
			! Inital Polymer viscosity
			visc_p0=visc_0-visc_s
			! Initial polymer viscosity at rest
			visc_p=visc_p0
			! Initial liquid phase viscosity
			fs%visc_l=visc_s+visc_p0
			visc_norm=fs%visc_l/(visc_s+visc_p0)
	 	end block solvent_viscosity
		
		! Create a FENE model 
		create_fene: block 
			use hypre_uns_class,   only: gmres
			use multiscalar_class, only: quick,bquick
			use fene_class,        only: FENECR
			integer :: i,j,k
			! Create FENE model solver
			fm=fene(cfg=cfg,model=FENECR,scheme=bquick,name='FENE')
			! Assign aritifical stress diffusivisty
			call param_read('Stress diffusivisty',stress_diff)
			fm%diff=stress_diff
			! Assign density value of 1
			fm%rho=1.00_WP
			! Maximum extensibility of polymer chain
			call param_read('Maximum extension of polymer chain',Lmax)
			! Relaxation time for polymer
			call param_read('Polymer Relaxation Time',lambda)
			! Configure the scalar solver
			! Configure the scalar solver
			ss=ddadi(cfg=cfg,name='Scalar',nst=13)
			! Setup the solver
			call fm%setup(implicit_solver=ss)
			! Intalize conformation tensor to identity matrix
			fm%SC(:,:,:,1)=1.00_WP !Cxx
			fm%SC(:,:,:,2)=0.00_WP !Cyx
			fm%SC(:,:,:,3)=0.00_WP !Czx
			fm%SC(:,:,:,4)=1.00_WP !Cyy
			fm%SC(:,:,:,5)=0.00_WP !Czy
			fm%SC(:,:,:,6)=1.00_WP !Czz
			! Calculate the relaxation function
			call fm%get_relaxationFunction(fR,Lmax)
			! Build stress tensor
			do k=fs%cfg%kmino_,fs%cfg%kmaxo_
				do j=fs%cfg%jmino_,fs%cfg%jmaxo_
					do i=fs%cfg%imino_,fs%cfg%imaxo_
						fm%T(i,j,k,1)=(visc_p(i,j,k)/lambda)*fR(i,j,k,1)
						fm%T(i,j,k,2)=(visc_p(i,j,k)/lambda)*fR(i,j,k,2)
						fm%T(i,j,k,3)=(visc_p(i,j,k)/lambda)*fR(i,j,k,3)
						fm%T(i,j,k,4)=(visc_p(i,j,k)/lambda)*fR(i,j,k,4)
						fm%T(i,j,k,5)=(visc_p(i,j,k)/lambda)*fR(i,j,k,5)
						fm%T(i,j,k,6)=(visc_p(i,j,k)/lambda)*fR(i,j,k,6)
					end do
				end do
			end do
			! Get its divergence 
			call fm%get_divT(fs)  
	  	end block create_fene 
	   

	   	! Create surfmesh object for interface polygon output
      	create_smesh: block
			smesh=surfmesh(nvar=0,name='plic')
			call vf%update_surfmesh(smesh)
      	end block create_smesh


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
           call ens_out%add_surface('plic',smesh)
		   call ens_out%add_scalar('SRmag',SRmag)
           call ens_out%add_scalar('normvisc_l',visc_norm)
		   call ens_out%add_scalar('trC',fm%trC)
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
	   
		! print *, 'epsilon',epsilon(0.0_WP)
		! print *, '0.1+epsilon',0.1_WP+epsilon(0.0_WP)
	   
	end subroutine simulation_init
	
	
	!> Perform an NGA2 simulation
	subroutine simulation_run
		implicit none
		
		! Perform time integration
	   	do while (.not.time%done())
			
			! Increment time
		   	call fs%get_cfl(time%dt,time%cfl)
		   	call time%adjust_dt()
			call time%increment()

			! Calculate SR
			call fs%get_strainrate(SR=SR)

			! Model shear thinning viscosity
			update_viscosity: block
			   integer :: i,j,k
			   do k=fs%cfg%kmino_,fs%cfg%kmaxo_
				  do j=fs%cfg%jmino_,fs%cfg%jmaxo_
					 do i=fs%cfg%imino_,fs%cfg%imaxo_
						! Strain rate magnitude
						SRmag(i,j,k)=sqrt(2.00_WP*SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
						! ! Carreau Model
						! fs%visc_l(i,j,k)=visc_inf+(visc_0-visc_inf)*(1.00_WP+(alpha*SRmag(i,j,k))**2)**((n-1.00_WP)/2.00_WP)
						! visc_norm(i,j,k)=fs%visc_l(i,j,k)/visc_0
						! Hybird model
						visc_p(i,j,k)=visc_p0*(1.00_WP+(alpha*SRmag(i,j,k))**2)**((n-1.00_WP)/2.00_WP)
						fs%visc_l(i,j,k)=visc_s+visc_p(i,j,k)
						visc_norm(i,j,k)=fs%visc_l(i,j,k)/(visc_s+visc_p0)
					 end do
				  end do
			   end do
			   call fs%cfg%sync(fs%visc_l)
			end block update_viscosity

			! Remember old scalars
			fm%SCold=fm%SC
			
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
				   ! Temp trSC array for Bquick  
				   trSC_(:,:,:,1)=SC_(:,:,:,1)+SC_(:,:,:,2)+SC_(:,:,:,3)
				end block pre_check
				
				! Check boundedess of explicit SC calculation
				call fm%metric_modification(SC=trSC_,SCmin=3.0_WP,SCmax=100.0_WP)
	
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
				   call fm%get_CgradU(gradu,CgradU)    
				   ! Calculate the relaxation function
				   call fm%get_relaxationFunction(fR,Lmax)     
				   ! Add source terms to calculated residual
				   resSC=resSC+(CgradU-(fR/lambda))*time%dt
				end block fene
	
				! Form implicit residual
				call fm%solve_implicit(time%dt,resSC,fs%U,fs%V,fs%W)
	
				! Apply this residual
				fm%SC=2.0_WP*fm%SC-fm%SCold+resSC
	
				! Apply other boundary conditions on the resulting field
				call fm%apply_bcond(time%t,time%dt)
				! ===================================================
				
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
				! 	! Calculate the relaxation function
				! 	call fm%get_relaxationFunction(fR,Lmax)
				! 	! Build stress tensor
				! 	do k=fs%cfg%kmino_,fs%cfg%kmaxo_
				! 		do j=fs%cfg%jmino_,fs%cfg%jmaxo_
				! 		    do i=fs%cfg%imino_,fs%cfg%imaxo_
				! 				fm%T(i,j,k,1)=(visc_p(i,j,k)/lambda)*fR(i,j,k,1)
				! 				fm%T(i,j,k,2)=(visc_p(i,j,k)/lambda)*fR(i,j,k,2)
				! 				fm%T(i,j,k,3)=(visc_p(i,j,k)/lambda)*fR(i,j,k,3)
				! 				fm%T(i,j,k,4)=(visc_p(i,j,k)/lambda)*fR(i,j,k,4)
				! 				fm%T(i,j,k,5)=(visc_p(i,j,k)/lambda)*fR(i,j,k,5)
				! 				fm%T(i,j,k,6)=(visc_p(i,j,k)/lambda)*fR(i,j,k,6)
				! 		    end do
				! 		end do
				! 	 end do    
				! 	! Get its divergence 
				! 	call fm%get_divT(fs) 
				! 	! Add divT to momentum equation 
				! 	do k=fs%cfg%kmin_,fs%cfg%kmax_
                !   		do j=fs%cfg%jmin_,fs%cfg%jmax_
                !     		do i=fs%cfg%imin_,fs%cfg%imax_
				! 				! Use volume fraction to apply divergence of polymer stress
				! 					! if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+vf%VF(i,j,k)*fm%divT(i,j,k,1)*time%dt !> x face/U velocity
				! 					! if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+vf%VF(i,j,k)*fm%divT(i,j,k,2)*time%dt !> y face/V velocity
				! 					! if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+vf%VF(i,j,k)*fm%divT(i,j,k,3)*time%dt !> z face/W velocity
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
				! ====================================================
				
				! Increment sub-iteration counter
			   	time%it=time%it+1
				
			end do
			
			! Recompute interpolated velocity and divergence
		   	call fs%interp_vel(Ui,Vi,Wi)
			call fs%get_div()
			
			! Output to ensight
		   	if (ens_evt%occurs()) then
				call vf%update_surfmesh(smesh)
				call ens_out%write_data(time%t)
         	end if
			
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
		deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu,resSC,visc_p,fR,CgradU)
	   
	end subroutine simulation_final
	
	
end module simulation