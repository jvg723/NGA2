!> Various definitions and tools for running an NGA2 simulation
module simulation
	use precision,         only: WP
	use geometry,          only: cfg
	use hypre_str_class,   only: hypre_str
	use ddadi_class,       only: ddadi
	use tpns_class,        only: tpns
	use vfs_class,         only: vfs
   	! use tpscalar_class,    only: tpscalar
	use tpfene_class,      only: tpfene
	use timetracker_class, only: timetracker
	use ensight_class,     only: ensight
	use surfmesh_class,    only: surfmesh
	use event_class,       only: event
	use monitor_class,     only: monitor
	implicit none
	private
	
	!> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
	type(hypre_str),   public :: ps
	type(ddadi),       public :: vs,ss
	type(tpns),        public :: fs
	type(vfs),         public :: vf
	type(tpfene),      public :: nn
   	! type(tpscalar),    public :: sc
	type(timetracker), public :: time
	
	!> Ensight postprocessing
	type(surfmesh) :: smesh
	type(ensight)  :: ens_out
	type(event)    :: ens_evt
	
	!> Simulation monitor file
	type(monitor) :: mfile,cflfile,scfile
	
	public :: simulation_init,simulation_run,simulation_final
	
	!> Private work arrays
	real(WP), dimension(:,:,:,:), 	allocatable :: resSC,SCtmp
	real(WP), dimension(:,:,:,:,:), allocatable :: gradu 
   	real(WP), dimension(:,:,:), 	allocatable :: resU,resV,resW
	real(WP), dimension(:,:,:), 	allocatable :: Ui,Vi,Wi
	
	!> Problem definition
	real(WP), dimension(3) :: center
	real(WP) :: radius,depth

contains


   !> Function that defines a level set function for a falling drop problem
	function levelset_falling_drop(xyz,t) result(G)
		implicit none
		real(WP), dimension(3),intent(in) :: xyz
		real(WP), intent(in) :: t
		real(WP) :: G
		! Create the droplet
	   G=radius-sqrt(sum((xyz-center)**2))
	   ! Add the pool
	   G=max(G,depth-xyz(2))
	end function levelset_falling_drop
	
	
	!> Initialization of problem solver
	subroutine simulation_init
		use param, only: param_read
		implicit none
		
		
		! Allocate work arrays
	   allocate_work_arrays: block
         	allocate(resSC(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
			allocate(SCtmp(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
			allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
		   	allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
			allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
			real(WP) :: vol,area
			integer, parameter :: amr_ref_lvl=4
			! Create a VOF solver
         	call vf%initialize(cfg=cfg,reconstruction_method=lvira,name='VOF',store_detailed_flux=.true.)
		   	! Initialize to a droplet and a pool
		   	center=[0.0_WP,0.05_WP,0.0_WP]
		   	radius=0.01_WP
		   	depth =0.025_WP
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
							call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_falling_drop,0.0_WP,amr_ref_lvl)
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
			use hypre_str_class, only: pcg_pfmg2
         	use mathtools,       only: Pi
			! Create flow solver
			fs=tpns(cfg=cfg,name='Two-phase NS')
			! Assign constant viscosity to each phase
			call param_read('Liquid dynamic viscosity',fs%visc_l)
			call param_read('Gas dynamic viscosity',fs%visc_g)
			! Assign constant density to each phase
		   	call param_read('Liquid density',fs%rho_l)
			call param_read('Gas density',fs%rho_g)
			! Read in surface tension coefficient
			call param_read('Surface tension coefficient',fs%sigma)
			call param_read('Static contact angle',fs%contact_angle)
         	fs%contact_angle=fs%contact_angle*Pi/180.0_WP
			! Assign acceleration of gravity
			call param_read('Gravity',fs%gravity)
			! Configure pressure solver
			ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         	ps%maxlevel=10
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
	   
      
      	! Create a two phase fene solver
      	create_fene: block
		  	use fene_class,		only: fenep
      		integer :: i,j,k
      		! Create fene solver
      		call nn%tpfene_initialize(cfg=cfg,model=fenep,name='fenep')
			! Maximum extensibility of polymer chain
			call param_read('Maximum polymer extensibility',nn%Lmax)
			! Relaxation time for polymer
			call param_read('Polymer relaxation time',nn%trelax)
			! Polymer viscosity at zero strain rate
			call param_read('Polymer viscosity',nn%visc)
      		! Assign zero diffusivity
      		nn%diff=0.0_WP
      		! Setup without an implicit solver
      		call nn%setup()
      		! Initialize scalar fields
      		do k=cfg%kmino_,cfg%kmaxo_
      		   	do j=cfg%jmino_,cfg%jmaxo_
      		   	   	do i=cfg%imino_,cfg%imaxo_
      		   	    	! Liquid scalar
      		   	    	if (vf%VF(i,j,k).gt.0.0_WP) then
      		   	     	   	! We are in the liquid
      		   	     	   	if (cfg%ym(j).gt.depth) then
      		   	     	    	! We are above the pool (trC=15)
								nn%SC(i,j,k,1)=5.0_WP !< Cxx
								nn%SC(i,j,k,2)=0.0_WP !< Cxy
								nn%SC(i,j,k,3)=0.0_WP !< Cxx
								nn%SC(i,j,k,4)=5.0_WP !< Cyy
								nn%SC(i,j,k,5)=0.0_WP !< Cyz
								nn%SC(i,j,k,6)=5.0_WP !< Czz
      		   	     	   	else
      		   	     	    	! We are in the pool (trC=3)
								nn%SC(i,j,:k,1)=1.0_WP !< Cxx
								nn%SC(i,j,:k,2)=0.0_WP !< Cxy
								nn%SC(i,j,:k,3)=0.0_WP !< Cxx
								nn%SC(i,j,:k,4)=1.0_WP !< Cyy
								nn%SC(i,j,:k,5)=0.0_WP !< Cyz
								nn%SC(i,j,:k,6)=1.0_WP !< Czz
      		   	     	   	end if
      		   	     	end if
      		   	     	! ! Gas scalar
						! if (vf%VF(i,j,k).lt.1.0_WP) then
						! 	! We are in the gas
						! 	nn%SC(i,j,k,:)=0.0_WP
						! end if
      		   	   	end do
      		   	end do
      		end do
      	end block create_fene
      
      
	   	! Create surfmesh object for interface polygon output
      	create_smesh: block
      	   smesh=surfmesh(nvar=0,name='plic')
      	   call vf%update_surfmesh(smesh)
      	end block create_smesh


	   	! Add Ensight output
	   	create_ensight: block
        	integer :: nsc
			! Create Ensight output from cfg
			ens_out=ensight(cfg=cfg,name='FallingDrop')
			! Create event for Ensight output
			ens_evt=event(time=time,name='Ensight output')
			call param_read('Ensight output period',ens_evt%tper)
			! Add variables to output
			call ens_out%add_vector('velocity',Ui,Vi,Wi)
			call ens_out%add_scalar('VOF',vf%VF)
			call ens_out%add_scalar('pressure',fs%P)
			call ens_out%add_scalar('curvature',vf%curv)
         	call ens_out%add_surface('plic',smesh)
        	do nsc=1,nn%nscalar
            	call ens_out%add_scalar(trim(nn%SCname(nsc)),nn%SC(:,:,:,nsc))
        	end do
		   ! Output to ensight
		   if (ens_evt%occurs()) call ens_out%write_data(time%t)
		end block create_ensight
	   
	   
	   	! Create a monitor file
	   	create_monitor: block
        	integer :: nsc
		   	! Prepare some info about fields
			call fs%get_cfl(time%dt,time%cfl)
			call fs%get_max()
			call vf%get_max()
         	call nn%get_max(VF=vf%VF)
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
        	! Create scalar monitor
        	scfile=monitor(nn%cfg%amRoot,'scalar')
        	call scfile%add_column(time%n,'Timestep number')
        	call scfile%add_column(time%t,'Time')
        	do nsc=1,nn%nscalar
				call scfile%add_column(nn%SCmin(nsc),trim(nn%SCname(nsc))//'_min')
				call scfile%add_column(nn%SCmax(nsc),trim(nn%SCname(nsc))//'_max')
				call scfile%add_column(nn%SCint(nsc),trim(nn%SCname(nsc))//'_int')
         	end do
         	call scfile%write()
	   end block create_monitor
	   
	   
	end subroutine simulation_init
	
	
	!> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
	subroutine simulation_run
    	use tpns_class, only: static_contact,arithmetic_visc
		implicit none
		
		! Perform time integration
	   	do while (.not.time%done())
			
			! Increment time
			call fs%get_cfl(time%dt,time%cfl)
			call time%adjust_dt()
			call time%increment()
			
			! Remember old VOF
			vf%VFold=vf%VF
         
         	! Remember old SC
         	nn%SCold=nn%SC
         
			! Remember old velocity
			fs%Uold=fs%U
			fs%Vold=fs%V
			fs%Wold=fs%W
		   
			! Prepare old staggered density (at n)
			call fs%get_olddensity(vf=vf)
			
			! VOF solver step
		   	call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)

			! Calculate grad(U)
			call fs%get_gradU(gradU)
			
			! ============= SCALAR SOLVER =======================

			! Explicit calculation of dSC/dt from scalar equation
			call nn%get_dSCdt(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,VFold=vf%VFold,VF=vf%VF,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)

			! Add viscoleastic source terms
            viscoelastic_src: block
               use fene_class, only: fenep,lptt,eptt
               integer :: n
               ! Streching and distortion term
               call nn%get_CgradU(gradU,SCtmp)
               do n=1,6
                  resSC(:,:,:,n)=resSC(:,:,:,n)+vf%VF(:,:,:)*SCtmp(:,:,:,n)
               end do
               ! Relaxation term
               call nn%get_relax(SCtmp,time%dt)
               do n=1,6
                  resSC(:,:,:,n)=resSC(:,:,:,n)+vf%VF(:,:,:)*SCtmp(:,:,:,n)
               end do
            end block viscoelastic_src

    		! Now transport our phase-specific scalars
    		advance_scalar: block
    			integer :: nsc
    			real(WP) :: p,q
    			! Advance scalar fields
    			do nsc=1,nn%nscalar
    			   p=real(nn%phase(nsc),WP); q=1.0_WP-2.0_WP*p
    			   where (nn%mask.eq.0.and.vf%VF.ne.p) nn%SC(:,:,:,nsc)=((p+q*vf%VFold)*nn%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/(p+q*vf%VF)
    			   where (vf%VF.eq.p) nn%SC(:,:,:,nsc)=0.0_WP
    			end do
    		end block advance_scalar

			! Apply boundary conditions
			call nn%apply_bcond(time%t,time%dt)
			! ===================================================
         
			! Prepare new staggered viscosity (at n+1)
		   	call fs%get_viscosity(vf=vf,strat=arithmetic_visc)
         
		   	! Perform sub-iterations
		   	do while (time%it.le.time%itmax)
            
				! Build mid-time velocity
			   	fs%U=0.5_WP*(fs%U+fs%Uold)
				fs%V=0.5_WP*(fs%V+fs%Vold)
				fs%W=0.5_WP*(fs%W+fs%Wold)
				
				! ============= VELOCITY SOLVER ======================
				! Preliminary mass and momentum transport step at the interface
			   	call fs%prepare_advection_upwind(dt=time%dt)
				
				! Explicit calculation of drho*u/dt from NS
			   	call fs%get_dmomdt(resU,resV,resW)
				
				! Add momentum source terms
			   	call fs%addsrc_gravity(resU,resV,resW)

				! Add polymer stress term
					polymer_stress: block
					use fene_class, only: fenep
					integer :: i,j,k,n
					real(WP), dimension(:,:,:), allocatable :: Txy,Tyz,Tzx
					real(WP), dimension(:,:,:,:), allocatable :: stress
					real(WP) :: coeff
					! Allocate work arrays
					allocate(stress(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:6))
					allocate(Txy   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
					allocate(Tyz   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
					allocate(Tzx   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
					! Build liquid stress tensor
					select case (nn%model)
					case (fenep)
						call nn%get_relax(stress,time%dt)
						do n=1,6
							stress(:,:,:,n)=-nn%visc*vf%VF*stress(:,:,:,n)          
						end do
				   	end select
				   	! Interpolate tensor components to cell edges
				   	do k=cfg%kmin_,cfg%kmax_+1
						do j=cfg%jmin_,cfg%jmax_+1
							do i=cfg%imin_,cfg%imax_+1
								Txy(i,j,k)=sum(fs%itp_xy(:,:,i,j,k)*stress(i-1:i,j-1:j,k,2))
								Tyz(i,j,k)=sum(fs%itp_yz(:,:,i,j,k)*stress(i,j-1:j,k-1:k,5))
								Tzx(i,j,k)=sum(fs%itp_xz(:,:,i,j,k)*stress(i-1:i,j,k-1:k,3))
							end do
						end do
				   	end do
				   	! Add divergence of stress to residual
				   	do k=fs%cfg%kmin_,fs%cfg%kmax_
					  	do j=fs%cfg%jmin_,fs%cfg%jmax_
							do i=fs%cfg%imin_,fs%cfg%imax_
								if (fs%umask(i,j,k).eq.0) resU(i,j,k)=resU(i,j,k)+sum(fs%divu_x(:,i,j,k)*stress(i-1:i,j,k,1))&
								&                                                +sum(fs%divu_y(:,i,j,k)*Txy(i,j:j+1,k))     &
								&                                                +sum(fs%divu_z(:,i,j,k)*Tzx(i,j,k:k+1))
								if (fs%vmask(i,j,k).eq.0) resV(i,j,k)=resV(i,j,k)+sum(fs%divv_x(:,i,j,k)*Txy(i:i+1,j,k))     &
								&                                                +sum(fs%divv_y(:,i,j,k)*stress(i,j-1:j,k,4))&
								&                                                +sum(fs%divv_z(:,i,j,k)*Tyz(i,j,k:k+1))
								if (fs%wmask(i,j,k).eq.0) resW(i,j,k)=resW(i,j,k)+sum(fs%divw_x(:,i,j,k)*Tzx(i:i+1,j,k))     &
								&                                                +sum(fs%divw_y(:,i,j,k)*Tyz(i,j:j+1,k))     &                  
								&                                                +sum(fs%divw_z(:,i,j,k)*stress(i,j,k-1:k,6))        
						 	end do
					  	end do
				   	end do
				   	! Clean up
				   	deallocate(stress,Txy,Tyz,Tzx)
				end block polymer_stress
				
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
         	call nn%get_max(VF=vf%VF)
			call mfile%write()
			call cflfile%write()
         call scfile%write()
			
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
	   deallocate(resSC,resU,resV,resW,Ui,Vi,Wi)
	   
	end subroutine simulation_final
	
	
end module simulation