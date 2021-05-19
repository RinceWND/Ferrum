! parallel_mpi.f

! subroutines for communicating between processors using MPI

!______________________________________________________________________
!
      subroutine initialize_mpi
!----------------------------------------------------------------------
!  Set up MPI execution environment and define the POM communicator
!______________________________________________________________________

      use glob_domain, only: error_status, is_master, master_task
     &                     , my_task     , POM_COMM , POM_COMM_COARSE
      use mpi        , only: MPI_COMM_WORLD

      implicit none

      integer ierr

! initiate MPI environment
      call mpi_init(ierr)

! determine processor rank
      call mpi_comm_rank(MPI_COMM_WORLD,my_task,ierr)

      POM_COMM        = MPI_COMM_WORLD
      POM_COMM_COARSE = MPI_COMM_WORLD
      master_task  = 0
      error_status = 0

      is_master = ( my_task == master_task )

      return

      end

!______________________________________________________________________
!
      subroutine finalize_mpi
!----------------------------------------------------------------------
!  Terminate the MPI execution environment
!______________________________________________________________________

      use glob_domain, only: is_master
      use mpi        , only: mpi_finalize

      implicit none

      integer ierr

! terminate MPI environment : TODO: Check returned error code?
      call mpi_finalize(ierr)
      print *, "Finalized with ", ierr
! stop every process except master
      if ( .not.is_master ) stop


      end

!______________________________________________________________________
!
      subroutine distribute_mpi
!----------------------------------------------------------------------
!  Distribute the model domain across processors
!______________________________________________________________________

      use glob_domain
      use mpi        , only: mpi_comm_size

      implicit none

      integer i,j,ierr,nproc,nproc_x,nproc_y

! determine the number of processors
      call mpi_comm_size(POM_COMM,nproc,ierr)

! check number of processors
      if ( nproc /= n_proc ) then
        error_status = 1
        if ( is_master ) print '(a,i3,''/='',i3//a)'
     &  , 'Incompatible number of processors',nproc,n_proc
     &  , 'POM terminated with error'
        call finalize_mpi
        stop
      end if

! determine the number of processors in x
      if ( mod(im_global-2,im_local-2) == 0 ) then
        nproc_x = (im_global-2)/(im_local-2)
      else
        nproc_x = (im_global-2)/(im_local-2) + 1
      end if

! determine the number of processors in y
      if ( mod(jm_global-2,jm_local-2) == 0 ) then
        nproc_y = (jm_global-2)/(jm_local-2)
      else
        nproc_y = (jm_global-2)/(jm_local-2) + 1
      end if

! check local size
      if ( nproc_x*nproc_y > n_proc ) then
        error_status = 1
        if ( is_master ) print '(a//a)'
     &  , 'im_local or jm_local is too low'
     &  , 'POM terminated with error'
        call finalize_mpi
        stop
      end if

! detemine global and local indices
      im = im_local
      do i=1,im_local
        i_global(i) = 0
      end do
      do i=1,im
        i_global(i) = i + mod(my_task,nproc_x)*(im-2)
        if ( i_global(i) > im_global ) then
          im = i-1
          i_global(i) = 0
          cycle
        end if
      end do
      imm1 = im-1
      imm2 = im-2

      jm = jm_local
      do j=1,jm_local
        j_global(j) = 0
      end do
      do j=1,jm
        j_global(j) = j + (my_task/nproc_x)*(jm-2)
        if ( j_global(j) > jm_global ) then
          jm = j-1
          j_global(j) = 0
          cycle
        end if
      end do
      jmm1 = jm-1
      jmm2 = jm-2

      kbm1 = kb-1
      kbm2 = kb-2
! detemine global and local indices for coarse grids
!      im_coarse=im_local_coarse
!      do i=1,im_local_coarse
!        i_global_coarse(i)=0
!      end do
!      do i=1,im_coarse
!        i_global_coarse(i)=i+mod(my_task,nproc_x)*(im_coarse-2)
!        if(i_global_coarse(i).gt.im_global_coarse) then
!          im_coarse=i-1
!          i_global_coarse(i)=0
!          cycle
!        end if
!      end do
!      imm1=im-1
!      imm2=im-2

!      jm_coarse=jm_local_coarse
!      do j=1,jm_local_coarse
!        j_global_coarse(j)=0
!      end do
!      do j=1,jm_coarse
!        j_global_coarse(j)=j+(my_task/nproc_x)*(jm_coarse-2)
!        if(j_global_coarse(j).gt.jm_global_coarse) then
!          jm_coarse=j-1
!          j_global_coarse(j)=0
!          cycle
!        end if
!      end do
!!      jmm1=jm-1
!!      jmm2=jm-2

!!      kbm1=kb-1
!!      kbm2=kb-2

! determine the neighbors (tasks)
      n_east  = my_task + 1
      n_west  = my_task - 1
      n_north = my_task + nproc_x
      n_south = my_task - nproc_x

      if ( mod(n_east  ,nproc_x) == 0 ) n_east = -1
      if ( mod(n_west+1,nproc_x) == 0 ) n_west = -1
      if (  n_north         /nproc_x == nproc_y ) n_north = -1
      if ( (n_south+nproc_x)/nproc_x == 0       ) n_south = -1

      return
      end
!______________________________________________________________________
!fhx: a new distribute mpi for wind and assim data to do the interpolation.2010/12/06
!
      subroutine distribute_mpi_coarse
!----------------------------------------------------------------------
!  Distribute the wind/assim data domain across processors
!______________________________________________________________________

      use glob_domain
      use mpi        , only: mpi_comm_size

      implicit none

      integer i,j,ierr,nproc,nproc_x,nproc_y

! determine the number of processors
      call mpi_comm_size(POM_COMM_coarse,nproc,ierr)

! check number of processors
      if ( nproc /= n_proc ) then
        error_status = 1
        if ( is_master ) print '(a//a)'
     &  , 'Incompatible number of processors for coarse grid'
     &  , 'POM terminated with error'
        call finalize_mpi
        stop
      end if

! determine the number of processors in x
      if ( mod(im_global_coarse-2,im_local_coarse-2) == 0 ) then
        nproc_x = (im_global_coarse-2)/(im_local_coarse-2)
      else
        nproc_x = (im_global_coarse-2)/(im_local_coarse-2) + 1
      end if

! determine the number of processors in y
      if ( mod(jm_global_coarse-2,jm_local_coarse-2) == 0 ) then
        nproc_y = (jm_global_coarse-2)/(jm_local_coarse-2)
      else
        nproc_y = (jm_global_coarse-2)/(jm_local_coarse-2) + 1
      end if

! check local size
      if ( nproc_x*nproc_y > n_proc ) then
        error_status = 1
        if ( is_master ) print '(a//a)'
     &  , 'im_local_coarse or jm_local_coarse is too low'
     &  , 'POM terminated with error'
        call finalize_mpi
        stop
      end if


! detemine global and local indices for coarse grids
      im_coarse = im_local_coarse
      do i=1,im_local_coarse
        i_global_coarse(i) = 0
      end do

      do i=1,im_coarse
        i_global_coarse(i) = i + mod(my_task,nproc_x)*(im_coarse-2)
        if ( i_global_coarse(i) > im_global_coarse ) then
          im_coarse = i-1
          i_global_coarse(i) = 0
          cycle
        end if
      end do

!      if(n_west.eq.-1) then   
!      do i=1,im_coarse
!        i_global_coarse(i)=i+mod(my_task,nproc_x)*(im_coarse-2)
!        if(i_global_coarse(i).gt.im_global_coarse) then
!          im_coarse=i-1
!          i_global_coarse(i)=0
!          cycle
!        end if
!      end do
!      end if
        
       
!      imm1=im-1
!      imm2=im-2

      jm_coarse = jm_local_coarse
      do j=1,jm_local_coarse
        j_global_coarse(j) = 0
      end do
      do j=1,jm_coarse
        j_global_coarse(j) = j + (my_task/nproc_x)*(jm_coarse-2)
        if ( j_global_coarse(j) > jm_global_coarse ) then
          jm_coarse = j-1
          j_global_coarse(j) = 0
          cycle
        end if
      end do
!      jmm1=jm-1
!      jmm2=jm-2

!      kbm1=kb-1
!      kbm2=kb-2

      return
      end

!______________________________________________________________________
!
      subroutine sum0i_mpi(work,to)
!----------------------------------------------------------------------
!  Send integer sum of WORK to node TO
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INTEGER, MPI_SUM

      implicit none

      integer, intent(in)    :: to
      integer, intent(inout) :: work

      integer ierr, tmp

! sum data
      call mpi_reduce(work,tmp,1,MPI_INTEGER,MPI_SUM,to,POM_COMM,ierr)
      work=tmp

      return

      end

!______________________________________________________________________
!
      subroutine sum0d_mpi(work,to)
!----------------------------------------------------------------------
!  Send real sum of WORK to node TO
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL, MPI_SUM

      implicit none

      integer , intent(in)    :: to
      real(rk), intent(inout) :: work

      integer  ierr, MPI_RK
      real(rk) tmp
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if
! sum data
      call mpi_reduce(work,tmp,1,MPI_RK,MPI_SUM,to,POM_COMM,ierr)
      work=tmp

      return

      end
!______________________________________________________________________
!
      subroutine max0d_mpi(work,to)
!----------------------------------------------------------------------
!  Send real maximum of WORK to node TO
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_MAX, MPI_REAL

      implicit none

      integer , intent(in)    :: to
      real(rk), intent(inout) :: work

      integer  ierr, MPI_RK
      real(rk) tmp
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if
! get max
      call mpi_reduce(work,tmp,1,MPI_RK,MPI_MAX,to,POM_COMM,ierr)
      work=tmp

      return

      end
!______________________________________________________________________
!
      subroutine bcast0i_mpi(work,from)
!----------------------------------------------------------------------
!  Send integer WORK to all nodes from node FROM
!______________________________________________________________________

      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INTEGER

      implicit none

      integer, intent(in)    :: from
      integer, intent(inout) :: work

      integer ierr

! broadcast data
      call mpi_bcast(work,1,MPI_INTEGER,from,POM_COMM,ierr)

      return

      end
!______________________________________________________________________
!
      subroutine bcast0d_mpi(work,from)
!----------------------------------------------------------------------
!  Send real WORK to all nodes from node FROM
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL

      implicit none

      integer , intent(in)    :: from
      real(rk), intent(inout) :: work

      integer ierr, MPI_RK
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if
! broadcast data
      call mpi_bcast(work,1,MPI_RK,from,POM_COMM,ierr)

      return

      end

!______________________________________________________________________
!
      subroutine exchange2d_mpi(work,nx,ny)
!----------------------------------------------------------------------
!  Exchange ghost cells around 2D local grids
!  One band at a time
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: my_task, n_east, n_north
     &                     , n_south, n_west, POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL
     &                     , MPI_STATUS_SIZE

      implicit none

      integer , intent(in)    :: nx,ny
      real(rk), intent(inout) :: work(nx,ny)

      integer  i,j
      integer  ierr, MPI_RK
      integer  istatus(MPI_STATUS_SIZE)
      real(rk) send_east(ny) ,recv_west(ny)
      real(rk) send_west(ny) ,recv_east(ny)
      real(rk) send_north(nx),recv_south(nx)
      real(rk) send_south(nx),recv_north(nx)
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if

! send ghost cell data to the east
      if ( n_east /= -1 ) then
        do j=1,ny
          send_east(j) = work(nx-1,j)
        end do
        call mpi_send(send_east,ny,MPI_RK,n_east,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the west
      if ( n_west /= -1 ) then
        call mpi_recv(recv_west,ny,MPI_RK,n_west,n_west,
     &                POM_COMM,istatus,ierr)
        do j=1,ny
          work(1,j) = recv_west(j)
        end do
      end if

! send ghost cell data to the west
      if ( n_west /= -1 ) then
        do j=1,ny
          send_west(j) = work(2,j)
        end do
        call mpi_send(send_west,ny,MPI_RK,n_west,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the east
      if ( n_east /= -1 ) then
        call mpi_recv(recv_east,ny,MPI_RK,n_east,n_east,
     &                POM_COMM,istatus,ierr)
        do j=1,ny
          work(nx,j) = recv_east(j)
        end do
      end if

! send ghost cell data to the north
      if ( n_north /= -1 ) then
        do i=1,nx
          send_north(i) = work(i,ny-1)
        end do
        call mpi_send(send_north,nx,MPI_RK,n_north,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the south
      if ( n_south /= -1 ) then
        call mpi_recv(recv_south,nx,MPI_RK,n_south,n_south,
     &                POM_COMM,istatus,ierr)
        do i=1,nx
          work(i,1) = recv_south(i)
        end do
      end if

! send ghost cell data to the south
      if ( n_south /= -1 ) then
        do i=1,nx
          send_south(i) = work(i,2)
        end do
        call mpi_send(send_south,nx,MPI_RK,n_south,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the north
      if ( n_north /= -1 ) then
        call mpi_recv(recv_north,nx,MPI_RK,n_north,n_north,
     &                POM_COMM,istatus,ierr)
        do i=1,nx
          work(i,ny) = recv_north(i)
        end do
      end if

      return

      end

!______________________________________________________________________
!
      subroutine exchange3d_mpi(work,nx,ny,nz)
!----------------------------------------------------------------------
!  Exchange ghost cells around 3D local grids
!  One band at a time
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: my_task, n_east, n_north
     &                     , n_south, n_west, POM_COMM
      use mpi        , only: mpi_recv, mpi_send
     &                     , MPI_DOUBLE_PRECISION, MPI_REAL
     &                     , MPI_STATUS_SIZE

      implicit none

      integer , intent(in)    :: nx,ny,nz
      real(rk), intent(inout) :: work(nx,ny,nz)

      integer  i,j,k
      integer  ierr, MPI_RK
      integer  istatus(MPI_STATUS_SIZE)
      real(rk) send_east(ny*nz) ,recv_west(ny*nz)
      real(rk) send_west(ny*nz) ,recv_east(ny*nz)
      real(rk) send_north(nx*nz),recv_south(nx*nz)
      real(rk) send_south(nx*nz),recv_north(nx*nz)
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if

! send ghost cell data to the east
      if ( n_east /= -1 ) then
        do k=1,nz
          do j=1,ny
            i = j + (k-1)*ny
            send_east(i) = work(nx-1,j,k)
          end do
        end do
        call mpi_send(send_east,ny*nz,MPI_RK,n_east,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the west
      if ( n_west /= -1 ) then
        call mpi_recv(recv_west,ny*nz,MPI_RK,n_west,n_west,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
          do j=1,ny
            i = j + (k-1)*ny
            work(1,j,k) = recv_west(i)
          end do
        end do
      end if

! send ghost cell data to the west
      if ( n_west /= -1 ) then
        do k=1,nz
         do j=1,ny
            i = j + (k-1)*ny
            send_west(i) = work(2,j,k)
          end do
        end do
        call mpi_send(send_west,ny*nz,MPI_RK,n_west,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the east
      if ( n_east /= -1 ) then
        call mpi_recv(recv_east,ny*nz,MPI_RK,n_east,n_east,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
         do j=1,ny
            i = j + (k-1)*ny
            work(nx,j,k) = recv_east(i)
          end do
        end do
      end if

! send ghost cell data to the north
      if ( n_north /= -1 ) then
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            send_north(j) = work(i,ny-1,k)
          end do
        end do
        call mpi_send(send_north,nx*nz,MPI_RK,n_north,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the south
      if ( n_south /= -1 ) then
        call mpi_recv(recv_south,nx*nz,MPI_RK,n_south,n_south,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            work(i,1,k) = recv_south(j)
          end do
        end do
      end if

! send ghost cell data to the south
      if ( n_south /= -1 ) then
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            send_south(j) = work(i,2,k)
          end do
        end do
        call mpi_send(send_south,nx*nz,MPI_RK,n_south,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the north
      if ( n_north /= -1 ) then
        call mpi_recv(recv_north,nx*nz,MPI_RK,n_north,n_north,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            work(i,ny,k) = recv_north(j)
          end do
        end do
      end if

      return

      end

!______________________________________________________________________
!
      subroutine psum_mpi( work, nx, sum_out )
!----------------------------------------------------------------------
!  TODO: Check if this crap is even needed 
! ayumi 2010/6/1
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: is_master, n_proc, POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL

      implicit none

      integer , intent(in)  :: nx
      real(rk), intent(in)  :: work(nx)
      real(rk), intent(out) :: sum_out

      integer  i, j
      integer  ierr, MPI_RK
      real(rk) buf( nx, n_proc )
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if

      buf     = 0.
      sum_out = 0.
    
      call mpi_gather( work, nx, MPI_RK, 
     &                 buf , nx, MPI_RK,
     &                 0, POM_COMM, ierr )                      


      if ( is_master ) then
        do j = 1, n_proc
          do i = 1, nx
            sum_out = sum_out + buf(i,j)
          end do
        end do
      end if

      return

      end
!______________________________________________________________________
!fhx:Toni:npg
!
      subroutine order2d_mpi(work2,work4,nx,ny)
!----------------------------------------------------------------------
!  Convert a 2nd order 2D matrix to special 4th order 2D matrix
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: my_task, n_east, n_north
     &                     , n_south, n_west, POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL
     &                     , MPI_STATUS_SIZE

      implicit none

      integer , intent(in)  :: nx,ny
      real(rk), intent(in)  :: work2(nx,ny)
      real(rk), intent(out) :: work4(0:nx,0:ny)

      integer  i,j
      integer  ierr, MPI_RK
      integer  istatus(MPI_STATUS_SIZE)
      real(rk) send_east(ny) ,recv_west(ny)
      real(rk) send_north(nx),recv_south(nx)
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if

      work4 = 0.
      work4(1:nx,1:ny) = work2

! send ghost cell data to the east
      if ( n_east /= -1 ) then
        do j=1,ny
          send_east(j) = work2(nx-2,j)
        end do
        call mpi_send(send_east,ny,MPI_RK,n_east,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the west
      if ( n_west /= -1 ) then
        call mpi_recv(recv_west,ny,MPI_RK,n_west,n_west,
     &                POM_COMM,istatus,ierr)
        do j=1,ny
          work4(0,j) = recv_west(j)
        end do
      end if

! send ghost cell data to the north
      if ( n_north /= -1 ) then
        do i=1,nx
          send_north(i) = work2(i,ny-2)
        end do
        call mpi_send(send_north,nx,MPI_RK,n_north,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the south
      if ( n_south /= -1 ) then
        call mpi_recv(recv_south,nx,MPI_RK,n_south,n_south,
     &                POM_COMM,istatus,ierr)
        do i=1,nx
          work4(i,0) = recv_south(i)
        end do
      end if

      return

      end

!______________________________________________________________________
!fhx:Toni:npg
!
      subroutine order3d_mpi(work2,work4,nx,ny,nz)
!----------------------------------------------------------------------
!  Convert a 2nd order 3D matrix to special 4th order 3D matrix
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain, only: my_task, n_east, n_north
     &                     , n_south, n_west, POM_COMM
      use mpi        , only: MPI_DOUBLE_PRECISION, MPI_REAL
     &                     , MPI_STATUS_SIZE

      implicit none

      integer , intent(in)  :: nx,ny,nz
      real(rk), intent(in)  :: work2(nx,ny,nz)
      real(rk), intent(out) :: work4(0:nx,0:ny,nz)

      integer  i,j,k
      integer  ierr, MPI_RK
      integer  istatus(MPI_STATUS_SIZE)
      real(rk) send_east(ny*nz) ,recv_west(ny*nz)
      real(rk) send_north(nx*nz),recv_south(nx*nz)
      
      if ( rk == 8 ) then
        MPI_RK = MPI_DOUBLE_PRECISION
      else
        MPI_RK = MPI_REAL
      end if

      work4 = 0.
      work4(1:nx,1:ny,1:nz) = work2

! send ghost cell data to the east
      if ( n_east /= -1 ) then
        do k=1,nz
          do j=1,ny
            i = j + (k-1)*ny
            send_east(i) = work2(nx-2,j,k)
          end do
        end do
        call mpi_send(send_east,ny*nz,MPI_RK,n_east,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the west
      if ( n_west /= -1 ) then
        call mpi_recv(recv_west,ny*nz,MPI_RK,n_west,n_west,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
          do j=1,ny
            i = j + (k-1)*ny
            work4(0,j,k) = recv_west(i)
          end do
        end do
      end if

! send ghost cell data to the north
      if ( n_north /= -1 ) then
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            send_north(j) = work2(i,ny-2,k)
          end do
        end do
        call mpi_send(send_north,nx*nz,MPI_RK,n_north,my_task,
     &                POM_COMM,ierr)
      end if
! recieve ghost cell data from the south
      if ( n_south /= -1 ) then
        call mpi_recv(recv_south,nx*nz,MPI_RK,n_south,n_south,
     &                POM_COMM,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j = i + (k-1)*nx
            work4(i,0,k) = recv_south(j)
          end do
        end do
      end if

      return

      end
