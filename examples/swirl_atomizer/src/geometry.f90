!> Various definitions and tools for initializing NGA2 config
module geometry
   use mpi_f08,        only: MPI_Group
   use config_class,   only: config
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Two configs
   !> Block 1 (b1) -> Turbulent annular pipe, single phase
   type(ibconfig), target, public :: cfg1
   !> Block 2 (b2) -> Swirl atomizer, two-phase wtih VOF+particle tracking
   type(config),   target, public :: cfg2

   !> Public declarations
   public :: geometry_init,Din,Dout
   public :: group1,isInGrp1
   public :: group2,isInGrp2

   !> Pipe diameter 
   real(WP) :: Dout,Din

   !> Two groups and partitions, along with logicals
   integer, dimension(3) :: partition1,partition2
   logical :: isInGrp1,isInGrp2
   type(MPI_Group) :: group1,group2
   

contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use parallel,    only: comm,group,nproc,rank
      use mpi_f08,     only: MPI_Group,MPI_Group_range_incl
      implicit none
      type(sgrid) :: grid
      
      mpi_groups: block
         integer, dimension(3,1) :: grange
         integer :: ierr
         ! Read in partition
         call param_read('1 Partition',partition1,short='p')
         call param_read('2 Partition',partition2,short='p')
         ! Create an MPI group along with logical for the turbulent pipe on the lowest ranks
         grange(:,1)=[0,product(partition1)-1,1]
         call MPI_Group_range_incl(group,1,grange,group1,ierr)
         isInGrp1=.false.; if (rank.le.product(partition1)-1) isInGrp1=.true.
         ! Create an MPI group along with logical for the swirl atomizer on the highest ranks
         grange(:,1)=[nproc-product(partition2),nproc-1,1]
         call MPI_Group_range_incl(group,1,grange,group2,ierr)
         isInGrp2=.false.; if (rank.ge.nproc-product(partition2)) isInGrp2=.true.
      end block mpi_groups

      ! Initialize grid 1
      if (isInGrp1) then
         ! Create config for block 1
         create_block1: block
            use sgrid_class, only: cartesian
            use parallel, only: group
            use ibconfig_class, only: bigot,sharp
            integer :: i,j,k,nx,ny,nz,no
            real(WP) :: Lx,Ly,Lz,dx,radius
            real(WP), dimension(:), allocatable :: x,y,z               
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
            ! Create partitioned grid
            cfg1=ibconfig(grp=group1,decomp=partition1,grid=grid)         
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
      end if
      
      ! Initialize grid 2
      if (isInGrp2) then
         ! Create config for block 2
         create_block2: block
            use sgrid_class, only: cartesian
            use parallel, only: group
            integer :: i,j,k,nx,ny,nz
            real(WP) :: Lx,Ly,Lz
            real(WP), dimension(:), allocatable :: x,y,z
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
            grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='swirl_atomizer')
            ! Create partitioned grid
            cfg2=config(grp=group2,decomp=partition2,grid=grid)
            ! Create masks for this config
            cfg2%VF=1.0_WP
         end block create_block2
      end if
   
   end subroutine geometry_init
   
   
end module geometry
