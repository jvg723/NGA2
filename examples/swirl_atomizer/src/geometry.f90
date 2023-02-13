!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Two configs
   !> Block 1 (b1) -> Turbulent annular pipe, single phase
   type(config), target, public :: cfg1
   !> Block 2 (b2) -> Swirl atomizer, two-phase wtih VOF
   type(config), target, public :: cfg2

   !> Pipe diameter (1=inner 2=outer)
   real(WP) :: D1,D2

   public :: geometry_init,D1,D2,get_VF

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
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         ! Read in grid and geometry definition
         call param_read('1 Pipe length',Lx)
         call param_read('1 Outer Pipe diameter',D2)
         call param_read('1 Inner Pipe diameter',D1)
         call param_read('1 ny',ny); allocate(y(ny+1))
         call param_read('1 nx',nx); allocate(x(nx+1))
         call param_read('1 nz',nz); allocate(z(nz+1))
         ! Domain lengths
         dx=Lx/real(nx,WP)
         no=6
         if (ny.gt.1) then
            Ly=D2+real(2*no,WP)*D2/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=D2+real(2*no,WP)*D2/real(ny-2*no,WP)
         else
            Lz=dx
         end if
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
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='annular_pipe')
         ! Read in partition
         call param_read('1 Partition',partition,short='p')
         ! Create partitioned grid
         cfg1=config(grp=group,decomp=partition,grid=grid)
         ! Create masks for this config
         do k=cfg1%kmin_,cfg1%kmax_
            do j=cfg1%jmin_,cfg1%jmax_
               do i=cfg1%imin_,cfg1%imax_
                  cfg1%VF(i,j,k)=max(get_VF(i,j,k,'SC'),epsilon(1.0_WP))
               end do
            end do
         end do
         call cfg1%sync(cfg1%VF)
         call cfg1%calc_fluid_vol()
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

   !> Get volume fraction for direct forcing
   function get_VF(cfg,i,j,k,dir) result(VF)
      use config_class, only: config
      implicit none
      class(config), intent(in) :: cfg
      integer, intent(in)       :: i,j,k
      character(len=*)          :: dir
      real(WP)                  :: VF
      real(WP)                  :: r,eta,lam,delta,VFx,VFy,VFz
      real(WP), dimension(3)    :: norm
      select case(trim(dir))
      case('U','u')
         delta=(cfg%dxm(i)*cfg%dy(j)*cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=cfg%ym(j)/r; norm(3)=cfg%zm(k)/r
      case('V','v')
         delta=(cfg%dx(i)*cfg%dym(j)*cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(cfg%y(j)**2+cfg%zm(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=cfg%y(j)/r; norm(3)=cfg%zm(k)/r
      case('W','w')
         delta=(cfg%dx(i)*cfg%dy(j)*cfg%dzm(k))**(1.0_WP/3.0_WP)
         r=sqrt(cfg%ym(j)**2+cfg%z(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=cfg%ym(j)/r; norm(3)=cfg%z(k)/r
      case default
         delta=(cfg%dx(i)*cfg%dy(j)*cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
         norm(1)=0.0_WP; norm(2)=cfg%ym(j)/r; norm(3)=cfg%zm(k)/r
      end select
      lam=sum(abs(norm)); eta=0.065_WP*(1.0_WP-lam**2)+0.39_WP
      if (r.ge.0.5_WP*D1) then 
         VF=0.5_WP*(1.0_WP-tanh((r-0.5_WP*D2)/(sqrt(2.0_WP)*lam*eta*delta+epsilon(1.0_WP))))
      else
         VF=0.0_WP
      end if
   end function get_VF
   
   
end module geometry
