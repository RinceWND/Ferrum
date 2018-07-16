! parallel_mpi.f

! subroutines for communicating between processors using MPI

!_______________________________________________________________________
      subroutine initialize_mpi
! set up MPI execution environment and define the POM communicator
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer ierr
! initiate MPI environment
      call mpi_init(ierr)
! determine processor rank
      call mpi_comm_rank(mpi_comm_world,my_task,ierr)
      pom_comm=mpi_comm_world
      pom_comm_coarse=mpi_comm_world
      master_task=0
      error_status=0
      return
      end

!_______________________________________________________________________
      subroutine finalize_mpi
! terminate the MPI execution environment
      implicit none
      include 'mpif.h'
      integer ierr
! terminate MPI environment
      call mpi_finalize(ierr)
      return
      end

!_______________________________________________________________________
      subroutine distribute_mpi
! distribute the model domain across processors
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer i,j,ierr,nproc,nproc_x,nproc_y

! determine the number of processors
      call mpi_comm_size(pom_comm,nproc,ierr)

! check number of processors
      if(nproc.ne.n_proc) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(a,i3''/=''i3//a)')
     $   'Incompatible number of processors',nproc,n_proc
     $  ,'POM terminated with error'
        call finalize_mpi
        stop
      end if

! determine the number of processors in x
      if(mod(im_global-2,im_local-2).eq.0) then
        nproc_x=(im_global-2)/(im_local-2)
      else
        nproc_x=(im_global-2)/(im_local-2) + 1
      end if

! determine the number of processors in y
      if(mod(jm_global-2,jm_local-2).eq.0) then
        nproc_y=(jm_global-2)/(jm_local-2)
      else
        nproc_y=(jm_global-2)/(jm_local-2) + 1
      end if

! check local size
      if(nproc_x*nproc_y.gt.n_proc) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(a//a)')
     $   'im_local or jm_local is too low','POM terminated with error'
        call finalize_mpi
        stop
      end if

! detemine global and local indices
      im=im_local
      do i=1,im_local
        i_global(i)=0
      end do
      do i=1,im
        i_global(i)=i+mod(my_task,nproc_x)*(im-2)
        if(i_global(i).gt.im_global) then
          im=i-1
          i_global(i)=0
          cycle
        end if
      end do
      imm1=im-1
      imm2=im-2

      jm=jm_local
      do j=1,jm_local
        j_global(j)=0
      end do
      do j=1,jm
        j_global(j)=j+(my_task/nproc_x)*(jm-2)
        if(j_global(j).gt.jm_global) then
          jm=j-1
          j_global(j)=0
          cycle
        end if
      end do
      jmm1=jm-1
      jmm2=jm-2

      kbm1=kb-1
      kbm2=kb-2
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
      n_east=my_task+1
      n_west=my_task-1
      n_north=my_task+nproc_x
      n_south=my_task-nproc_x

      if(mod(n_east,nproc_x).eq.0) n_east=-1
      if(mod(n_west+1,nproc_x).eq.0) n_west=-1
      if(n_north/nproc_x.eq.nproc_y) n_north=-1
      if((n_south+nproc_x)/nproc_x.eq.0) n_south=-1

      return
      end
!_______________________________________________________________________
!fhx: a new distribute mpi for wind and assim data to do the interpolation.2010/12/06
      subroutine distribute_mpi_coarse
! distribute the wind/assim data domain across processors
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer i,j,ierr,nproc,nproc_x,nproc_y

! determine the number of processors
      call mpi_comm_size(pom_comm_coarse,nproc,ierr)

! check number of processors
      if(nproc.ne.n_proc) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(a//a)')
     $'Incompatible number of processors C','POM terminated with error'
        call finalize_mpi
        stop
      end if

! determine the number of processors in x
      if(mod(im_global_coarse-2,im_local_coarse-2).eq.0) then
        nproc_x=(im_global_coarse-2)/(im_local_coarse-2)
      else
        nproc_x=(im_global_coarse-2)/(im_local_coarse-2) + 1
      end if

! determine the number of processors in y
      if(mod(jm_global_coarse-2,jm_local_coarse-2).eq.0) then
        nproc_y=(jm_global_coarse-2)/(jm_local_coarse-2)
      else
        nproc_y=(jm_global_coarse-2)/(jm_local_coarse-2) + 1
      end if

! check local size
      if(nproc_x*nproc_y.gt.n_proc) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(a//a)')
     $   'im_local_coarse or jm_local_coarse is too low',
     $   'POM terminated with error'
        call finalize_mpi
        stop
      end if


! detemine global and local indices for coarse grids
      im_coarse=im_local_coarse
      do i=1,im_local_coarse
        i_global_coarse(i)=0
      end do

      do i=1,im_coarse
        i_global_coarse(i)=i+mod(my_task,nproc_x)*(im_coarse-2)
        if(i_global_coarse(i).gt.im_global_coarse) then
          im_coarse=i-1
          i_global_coarse(i)=0
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

      jm_coarse=jm_local_coarse
      do j=1,jm_local_coarse
        j_global_coarse(j)=0
      end do
      do j=1,jm_coarse
        j_global_coarse(j)=j+(my_task/nproc_x)*(jm_coarse-2)
        if(j_global_coarse(j).gt.jm_global_coarse) then
          jm_coarse=j-1
          j_global_coarse(j)=0
          cycle
        end if
      end do
!      jmm1=jm-1
!      jmm2=jm-2

!      kbm1=kb-1
!      kbm2=kb-2

      return
      end

!_______________________________________________________________________
      subroutine sum0i_mpi(work,to)
! send sum of WORK to node TO
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer to
      integer work,tmp
      integer ierr

! sum data
      call mpi_reduce(work,tmp,1,mpi_integer,mpi_sum,to,pom_comm,ierr)
      work=tmp
      return
      end

!_______________________________________________________________________
      subroutine sum0d_mpi(work,to)
! send sum of WORK to node TO
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer to
      real(kind=rk) work,tmp
      integer ierr
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if
! sum data
      call mpi_reduce(work,tmp,1,mpi_rk,mpi_sum,to,pom_comm,ierr)
!      call mpi_reduce(work,tmp,1,mpi_double_precision
!     $     ,mpi_sum,to,pom_comm,ierr)
      work=tmp
      return
      end
!_______________________________________________________________________
      subroutine max0d_mpi(work,to)
! send sum of WORK to node TO
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer to
      real(kind=rk) work,tmp
      integer ierr
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if
! get max
      call mpi_reduce(work,tmp,1,mpi_rk,mpi_max,to,pom_comm,ierr)
      work=tmp
      return
      end
!_______________________________________________________________________
      subroutine min0d_all_mpi(work)
! send min of WORK to all nodes
      implicit none
      include 'mpif.h'
      include 'pom.h'
      real(kind=rk) work,tmp
      integer ierr
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if
! get max
      call mpi_allreduce(work,tmp,1,mpi_rk,mpi_min,pom_comm,ierr)
      work=tmp
      return
      end
!_______________________________________________________________________
      subroutine bcast0i_mpi(work,from)
! send WORK to all nodes from node FROM
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer from
      integer work
      integer ierr

! broadcast data
      call mpi_bcast(work,1,mpi_integer,from,pom_comm,ierr)
      return
      end
!_______________________________________________________________________
      subroutine bcast0d_mpi(work,from)
! send WORK to all nodes from node FROM
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer from
      real(kind=rk) work
      integer ierr
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if
! broadcast data
      call mpi_bcast(work,1,mpi_rk,from,pom_comm,ierr)
!      call mpi_bcast(work,1,mpi_double_precision,from,pom_comm,ierr)
      return
      end

!_______________________________________________________________________
      subroutine exchange2d_mpi(work,nx,ny)
! exchange ghost cells around 2d local grids
! one band at a time
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer nx,ny
      real(kind=rk) work(nx,ny)
      integer i,j
      integer ierr
      integer istatus(mpi_status_size)
      real(kind=rk) send_east(ny),recv_west(ny)
      real(kind=rk) send_west(ny),recv_east(ny)
      real(kind=rk) send_north(nx),recv_south(nx)
      real(kind=rk) send_south(nx),recv_north(nx)
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

! send ghost cell data to the east
      if(n_east.ne.-1) then
        do j=1,ny
          send_east(j)=work(nx-1,j)
        end do
        call mpi_send(send_east,ny,mpi_rk,n_east,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the west
      if(n_west.ne.-1) then
        call mpi_recv(recv_west,ny,mpi_rk,n_west,n_west,
     $                pom_comm,istatus,ierr)
        do j=1,ny
          work(1,j)=recv_west(j)
        end do
      end if

! send ghost cell data to the west
      if(n_west.ne.-1) then
        do j=1,ny
          send_west(j)=work(2,j)
        end do
        call mpi_send(send_west,ny,mpi_rk,n_west,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the east
      if(n_east.ne.-1) then
        call mpi_recv(recv_east,ny,mpi_rk,n_east,n_east,
     $                pom_comm,istatus,ierr)
        do j=1,ny
          work(nx,j)=recv_east(j)
        end do
      end if

! send ghost cell data to the north
      if(n_north.ne.-1) then
        do i=1,nx
          send_north(i)=work(i,ny-1)
        end do
        call mpi_send(send_north,nx,mpi_rk,n_north,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the south
      if(n_south.ne.-1) then
        call mpi_recv(recv_south,nx,mpi_rk,n_south,n_south,
     $                pom_comm,istatus,ierr)
        do i=1,nx
          work(i,1)=recv_south(i)
        end do
      end if

! send ghost cell data to the south
      if(n_south.ne.-1) then
        do i=1,nx
          send_south(i)=work(i,2)
        end do
        call mpi_send(send_south,nx,mpi_rk,n_south,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the north
      if(n_north.ne.-1) then
        call mpi_recv(recv_north,nx,mpi_rk,n_north,n_north,
     $                pom_comm,istatus,ierr)
        do i=1,nx
          work(i,ny)=recv_north(i)
        end do
      end if

      return
      end

!_______________________________________________________________________
      subroutine exchange3d_mpi(work,nx,ny,nz)
! exchange ghost cells around 3d local grids
! one band at a time
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer nx,ny,nz
      real(kind=rk) work(nx,ny,nz)
      integer i,j,k
      integer ierr
      integer istatus(mpi_status_size)
      real(kind=rk) send_east(ny*nz),recv_west(ny*nz)
      real(kind=rk) send_west(ny*nz),recv_east(ny*nz)
      real(kind=rk) send_north(nx*nz),recv_south(nx*nz)
      real(kind=rk) send_south(nx*nz),recv_north(nx*nz)
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

! send ghost cell data to the east
      if(n_east.ne.-1) then
        do k=1,nz
          do j=1,ny
            i=j+(k-1)*ny
            send_east(i)=work(nx-1,j,k)
          end do
        end do
        call mpi_send(send_east,ny*nz,mpi_rk,n_east,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the west
      if(n_west.ne.-1) then
        call mpi_recv(recv_west,ny*nz,mpi_rk,n_west,n_west,
     $                pom_comm,istatus,ierr)
        do k=1,nz
          do j=1,ny
            i=j+(k-1)*ny
            work(1,j,k)=recv_west(i)
          end do
        end do
      end if

! send ghost cell data to the west
      if(n_west.ne.-1) then
        do k=1,nz
         do j=1,ny
            i=j+(k-1)*ny
            send_west(i)=work(2,j,k)
          end do
        end do
        call mpi_send(send_west,ny*nz,mpi_rk,n_west,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the east
      if(n_east.ne.-1) then
        call mpi_recv(recv_east,ny*nz,mpi_rk,n_east,n_east,
     $                pom_comm,istatus,ierr)
        do k=1,nz
         do j=1,ny
            i=j+(k-1)*ny
            work(nx,j,k)=recv_east(i)
          end do
        end do
      end if

! send ghost cell data to the north
      if(n_north.ne.-1) then
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            send_north(j)=work(i,ny-1,k)
          end do
        end do
        call mpi_send(send_north,nx*nz,mpi_rk,n_north,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the south
      if(n_south.ne.-1) then
        call mpi_recv(recv_south,nx*nz,mpi_rk,n_south,n_south,
     $                pom_comm,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            work(i,1,k)=recv_south(j)
          end do
        end do
      end if

! send ghost cell data to the south
      if(n_south.ne.-1) then
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            send_south(j)=work(i,2,k)
          end do
        end do
        call mpi_send(send_south,nx*nz,mpi_rk,n_south,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the north
      if(n_north.ne.-1) then
        call mpi_recv(recv_north,nx*nz,mpi_rk,n_north,n_north,
     $                pom_comm,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            work(i,ny,k)=recv_north(j)
          end do
        end do
      end if

      return
      end

!_______________________________________________________________________
      subroutine psum_mpi( work, nx, sum_out )
! 
! ayumi 2010/6/1
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer nx
      real(kind=rk), intent(in) ::  work(nx)
      real(kind=rk), intent(out) ::  sum_out
      real(kind=rk) buf( nx, n_proc )
      integer i, j
      integer ierr
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

      buf = 0.0
      sum_out = 0.0
    
      call mpi_gather( work,nx, mpi_rk, 
     $                 buf, nx, mpi_rk,
     $                 0, pom_comm, ierr )                      


      if ( my_task == master_task ) then
         do j = 1, n_proc
            do i = 1, nx
               sum_out = sum_out + buf(i,j)
            enddo
         enddo
      endif


      return

      end
!_______________________________________________________________________
!fhx:Toni:npg
      subroutine order2d_mpi(work2,work4,nx,ny)
! convert a 2nd order 2d matrix to special 4th order 2d matrix
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer nx,ny
      real(kind=rk) work2(nx,ny),work4(0:nx,0:ny)
      integer i,j
      integer ierr
      integer istatus(mpi_status_size)
      real(kind=rk) send_east(ny),recv_west(ny)
      real(kind=rk) send_north(nx),recv_south(nx)
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

      work4=0.
      do i=1,nx
        do j=1,ny
          work4(i,j)=work2(i,j)
        end do
      end do

! send ghost cell data to the east
      if(n_east.ne.-1) then
        do j=1,ny
          send_east(j)=work2(nx-2,j)
        end do
        call mpi_send(send_east,ny,mpi_rk,n_east,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the west
      if(n_west.ne.-1) then
        call mpi_recv(recv_west,ny,mpi_rk,n_west,n_west,
     $                pom_comm,istatus,ierr)
        do j=1,ny
          work4(0,j)=recv_west(j)
        end do
      end if

! send ghost cell data to the north
      if(n_north.ne.-1) then
        do i=1,nx
          send_north(i)=work2(i,ny-2)
        end do
        call mpi_send(send_north,nx,mpi_rk,n_north,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the south
      if(n_south.ne.-1) then
        call mpi_recv(recv_south,nx,mpi_rk,n_south,n_south,
     $                pom_comm,istatus,ierr)
        do i=1,nx
          work4(i,0)=recv_south(i)
        end do
      end if

      return
      end

!_______________________________________________________________________
!fhx:Toni:npg
      subroutine order3d_mpi(work2,work4,nx,ny,nz)
! convert a 2nd order 3d matrix to special 4th order 3d matrix
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer nx,ny,nz
      real(kind=rk) work2(nx,ny,nz),work4(0:nx,0:ny,nz)
      integer i,j,k
      integer ierr
      integer istatus(mpi_status_size)
      real(kind=rk) send_east(ny*nz),recv_west(ny*nz)
      real(kind=rk) send_north(nx*nz),recv_south(nx*nz)
      integer mpi_rk
      
      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

      work4=0.
      do i=1,nx
        do j=1,ny
          do k=1,nz
            work4(i,j,k)=work2(i,j,k)
          end do
        end do
      end do

! send ghost cell data to the east
      if(n_east.ne.-1) then
        do k=1,nz
          do j=1,ny
            i=j+(k-1)*ny
            send_east(i)=work2(nx-2,j,k)
          end do
        end do
        call mpi_send(send_east,ny*nz,mpi_rk,n_east,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the west
      if(n_west.ne.-1) then
        call mpi_recv(recv_west,ny*nz,mpi_rk,n_west,n_west,
     $                pom_comm,istatus,ierr)
        do k=1,nz
          do j=1,ny
            i=j+(k-1)*ny
            work4(0,j,k)=recv_west(i)
          end do
        end do
      end if

! send ghost cell data to the north
      if(n_north.ne.-1) then
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            send_north(j)=work2(i,ny-2,k)
          end do
        end do
        call mpi_send(send_north,nx*nz,mpi_rk,n_north,my_task,
     $                pom_comm,ierr)
      end if
! recieve ghost cell data from the south
      if(n_south.ne.-1) then
        call mpi_recv(recv_south,nx*nz,mpi_rk,n_south,n_south,
     $                pom_comm,istatus,ierr)
        do k=1,nz
          do i=1,nx
            j=i+(k-1)*nx
            work4(i,0,k)=recv_south(j)
          end do
        end do
      end if

      return
      end


