!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use ibconfig_class, only: ibconfig
   use precision,    only: WP
   implicit none
   private
   
   !> Two configs
   !> Block 1 (b1) -> Turbulent annular pipe, single phase
   type(ibconfig), target, public :: cfg1
   !> Block 2 (b2) -> Swirl atomizer, two-phase wtih VOF
   type(config),   target, public :: cfg2

   !> Pipe diameter 
   real(WP) :: Dout,Din

   public :: geometry_init,Din,Dout

contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Create config for block 1
      create_block1: block
         use sgrid_class, only: cartesian
         use parallel, only: group
         use ibconfig_class, only: bigot,sharp
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx,radius
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         ! Read in grid and geometry definition
         call param_read('1 Pipe length',Lx)
         call param_read('1 Outer Pipe diameter',Dout)
         call param_read('1 Inner Pipe diameter',Din)
         call param_read('1 ny',ny); allocate(y(ny+1))
         call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 nz',nz); allocate(z(nz+1))
         ! Domain lengths
         dx=Lx/real(nx,WP)
         no=6
         if (ny.gt.1) then
            Ly=Dout+real(2*no,WP)*Dout/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=Dout+real(2*no,WP)*Dout/real(ny-2*no,WP)
         else
            Lz=dx
         end if
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.5_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='annular_pipe')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=ibconfig(grp=group,decomp=partition,grid=grid)
         ! Create IB walls for this config
         !> Create IB field
         do k=cfg1%kmino_,cfg1%kmaxo_
            do j=cfg1%jmino_,cfg1%jmaxo_
               do i=cfg1%imino_,cfg1%imaxo_
                  radius=sqrt(cfg1%ym(j)**2+cfg1%zm(k)**2)
                  cfg1%Gib(i,j,k)=max(0.5_WP*Din-radius,radius-0.5_WP*Dout)
               end do
            end do
         end do
         !> Get normal vector
         call cfg1%calculate_normal()
         !> Get VF field
         call cfg1%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_block1
      
      ! Create config for block 2
      create_block2: block
         use sgrid_class, only: cartesian
         use parallel, only: group
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         ! Read in grid and geometry definition
         call param_read('2 Lx',Lx); call param_read('2 nx',nx); allocate(x(nx+1))
         call param_read('2 Ly',Ly); call param_read('2 ny',ny); allocate(y(ny+1))
         call param_read('2 Lz',Lz); call param_read('2 nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='swirl_atomizer')
         ! Read in partition
         call param_read('2 Partition',partition,short='p')
         ! Create partitioned grid
         cfg2=config(grp=group,decomp=partition,grid=grid)
         ! Create masks for this config
         cfg2%VF=1.0_WP
      end block create_block2
   

   end subroutine geometry_init
   
   
end module geometry
