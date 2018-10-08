!     eda:20110118:subroutine for drifter assmilation
c-------------------------------------------------------------------------
ceda: drifter assimilation modified from (in GFDL)
c     /home/xil/OIpsLag/gomc27_test87d_hcast2me_pseudo1.f
c     Drifter assimilation using the pseudo-Lagrangian Optimal
c     Interpolation (OI-psLag, see Ozgokmen et al, Assimilation of drifter 
c     observations in primitive equation models of midlatitude ocean 
c     circulation, JGR, 108(C7), 3238, 2003).
c-------------------------------------------------------------------------
!==================================================================
      subroutine assimdrf_main (d_in)

      use module_time
      implicit none
      include 'pom.h'
      
      type(date), intent(in) ::  d_in
      integer, dimension(4) :: itime1

      itime1(1)=d_in%year; itime1(2)=d_in%month
      itime1(3)=d_in%day; itime1(4)=d_in%hour

        call assimdrf_OIpsLag(itime1, 
     1   IM, JM, KB, u, v, east_e,north_e,
     2   zz, D, ub, vb, dz,dti, iint,my_task, master_task)


      end subroutine assimdrf_main
!==================================================================
      subroutine assimdrf_OIpsLag(itime1, 
     1     IM, JM, KB, u, v, alon, alat, zz, D,
     2     ub, vb, dz, dti, iint, my_task, master_task)

      implicit none
      include 'realkind'

      integer, parameter :: ndrfmax = 1000
      real(kind=rk), intent(in) :: dti
      integer, intent(in) :: iint,my_task, master_task
      integer, intent(in), dimension(4) :: itime1
      integer, intent(in) :: IM, JM, KB
      real(kind=rk), dimension (IM,JM), intent(in) :: alon, alat
      real(kind=rk), dimension (IM,JM), intent(in) :: D
      real(kind=rk), dimension (KB), intent(in) :: zz, dz
      real(kind=rk), dimension (IM,JM,KB), intent(inout) :: u,v, ub,vb


!local variable
      logical :: fexist, foundgridbox
      logical, external :: ingridbox

      character(*), parameter :: FileFormat =
     1     "(A7,1X,F10.6,1X,F10.6,1X,F10.3,1X,F10.3,1X,F10.3)"

      character :: infilee*18, FileNamee*30
      real(kind=rk), dimension( im, jm) :: dist2

      real(kind=rk), dimension( im, jm, kb ):: uold, vold

      integer :: i, j, k, n, ndrf
      integer :: IOStatus, iassimdrf
      integer, dimension (2) :: locind
      integer :: i0, j0, i1, j1

      real(kind=rk) :: gamma, usum, vsum,
     1     lat_in, lon_in, u_in, v_in, speed_in, 
     2     Lx, Ly, delu, delv, minlocal
      

      real(kind=rk),dimension (ndrfmax) :: ulag_o1,vlag_o1,ueu_b1,veu_b1

      real(kind=rk),dimension (ndrfmax) :: ulag_o, vlag_o, ueu_b, veu_b, 
     1     lat_drf, lon_drf, u_drf, v_drf

      character(len=20) id_in

      namelist/assimdrf_nml/ nassimdrf,nx, ny, beta, mindist,zd
      real(kind=rk) :: nassimdrf, nx, ny, beta, mindist,zd

         open(73,file='switch.nml',status='old')
         read(73, nml=assimdrf_nml)
         close(73)

      i1 = 1
      j1 = 1

      iassimdrf = max( int( nassimdrf * 24 * 3600 / dti ), 1 )
      if ( mod(iint-1, iassimdrf) /= 0 ) return

!-------------------------------------------------------------------
!eda:      Start drf assim
!-------------------------------------------------------------------
!
c     backup u and v

      uold = u
      vold = v

!eda: using d_in instead of itime1 in ori_code
         write( infilee, '( "uvlatlon",i4.4,3i2.2 )' ) 
     $        itime1(1), itime1(2), itime1(3), itime1(4)

         write(FileNamee,'(''in/assimdrf/'',a )') trim(infilee)

!eda: check if file exist or not
         inquire(file=trim(FileNamee),exist=fexist)   

      if(fexist) then

                if ( my_task == master_task ) then
                  write(*,'(/2a)') 
     $              "assimdrf_main : "
     $              , trim(infilee)
                endif

      ndrf = 0
      open(73, file=trim(FileNamee), status="old", iostat=IOStatus)

      if ( IOStatus == 0 ) then

         ! read one data entry
         read(73, FileFormat, iostat=IOStatus) 
     1        id_in, lat_in, lon_in, u_in, v_in, speed_in


         ! until the end of the file is reached
         do while ( IOStatus == 0 )

            ndrf = ndrf+1
            if ( ndrf > ndrfmax ) then
                write(*,'(a)') "Error: ndrf > ndrfmax."
                write(*,'(a)') "Use a larger ndrfmax."
               stop
            end if

            lat_drf(ndrf) = lat_in
            lon_drf(ndrf) = lon_in
            u_drf(ndrf) = u_in
            v_drf(ndrf) = v_in

            ! read the next data entry
            read(73, FileFormat, iostat=IOStatus) 
     1           id_in, lat_in, lon_in, u_in, v_in, speed_in

         end do

      end if
      close(73)

      do n = 1, ndrf

        ! find the gridbox that the drifter current position is in.
        ! (i1,j1) is the lower left point of the gridbox.

        ! find the nearest grid point (i0,j0)
        dist2 = 0.
        dist2 = (lon_drf(n) - alon(:,:))**2 + 
     1          (lat_drf(n) - alat(:,:))**2
        locind = minloc(dist2)
        i0 =  locind(1)
        j0 =  locind(2)

!eda: find local min dist (minlocal), if .gt. mindist,then return      

        minlocal=minval(dist2)
        if ( i0 == 1 .or. i0 == im .or. j0 == 1 .or. j0 == jm 
     &      .or. minlocal > mindist ) then

          foundgridbox = .false.

        else if ( ingridbox( alon(i0-1,j0-1), alon(i0-1,j0),
     1           alon(i0,j0), alon(i0,j0-1),
     2           alat(i0-1,j0-1), alat(i0-1,j0),
     3           alat(i0,j0), alat(i0,j0-1),
     4           lon_drf(n), lat_drf(n) ) ) then

          foundgridbox = .true.
          i1 = i0-1; j1 = j0-1

        else if ( ingridbox( alon(i0-1,j0), alon(i0-1,j0+1),
     1           alon(i0,j0+1), alon(i0,j0),
     2           alat(i0-1,j0), alat(i0-1,j0+1),
     3           alat(i0,j0+1), alat(i0,j0),
     4           lon_drf(n), lat_drf(n) ) ) then

          foundgridbox = .true.
          i1 = i0-1; j1 = j0

        else if ( ingridbox( alon(i0,j0), alon(i0,j0+1),
     1           alon(i0+1,j0+1), alon(i0+1,j0),
     2           alat(i0,j0), alat(i0,j0+1),
     3           alat(i0+1,j0+1), alat(i0+1,j0),
     4           lon_drf(n), lat_drf(n) ) ) then

          foundgridbox = .true.
          i1 = i0;   j1 = j0

        else if ( ingridbox( alon(i0,j0-1), alon(i0,j0),
     1           alon(i0+1,j0), alon(i0+1,j0-1),
     2           alat(i0,j0-1), alat(i0,j0),
     3           alat(i0+1,j0), alat(i0+1,j0-1),
     4           lon_drf(n), lat_drf(n) ) ) then

          foundgridbox = .true.
          i1 = i0;   j1 = j0-1

        else

          foundgridbox = .false.

        end if

        if ( foundgridbox ) then

c     Rotate input drifter velocities from west/east & south/north 
c     to curvillinear grid.
c     ulag, vlag: Lagrangian velocity

          call spher2curvi( u_drf(n), v_drf(n), 
     1      alon(i1,j1), alon(i1+1,j1), alat(i1,j1), 
     2      alat(i1+1,j1),ulag_o1(n), vlag_o1(n))

c     In the pseudo-Lagrangian assimilation, the model Eulerian velocity
c     is used instead of the model Lagrangian velocity. Eulerian velocity
c     at the drifter position is found by interpolating grid velocities.
c     ueu, veu: Eulerian velocity
c     Note that the positions of u and v are not in the center of the
c     grid cell.

          call BLINT( (alon(i1-1,j1)+alon(i1,j1))/2., 
     2      (alon(i1-1,j1+1)+alon(i1,j1+1))/2., 
     3      (alon(i1,j1+1)+alon(i1+1,j1+1))/2., 
     4      (alon(i1,j1)+alon(i1+1,j1))/2.,
     5      (alat(i1-1,j1)+alat(i1,j1))/2., 
     6      (alat(i1-1,j1+1)+alat(i1,j1+1))/2., 
     7      (alat(i1,j1+1)+alat(i1+1,j1+1))/2., 
     8      (alat(i1,j1)+alat(i1+1,j1))/2.,
     9      u(i1,j1,1), u(i1,j1+1,1), u(i1+1,j1+1,1), u(i1+1,j1,1),
     1      lon_drf(n), lat_drf(n), ueu_b1(n))

          call BLINT( (alon(i1,j1-1)+alon(i1,j1))/2., 
     2      (alon(i1,j1)+alon(i1,j1+1))/2., 
     3      (alon(i1+1,j1)+alon(i1+1,j1+1))/2., 
     4      (alon(i1+1,j1-1)+alon(i1+1,j1))/2.,
     5      (alat(i1,j1-1)+alat(i1,j1))/2., 
     6      (alat(i1,j1)+alat(i1,j1+1))/2., 
     7      (alat(i1+1,j1)+alat(i1+1,j1+1))/2., 
     8      (alat(i1+1,j1-1)+alat(i1+1,j1))/2.,
     9      v(i1,j1,1), v(i1,j1+1,1), v(i1+1,j1+1,1), v(i1+1,j1,1),
     1      lon_drf(n), lat_drf(n), veu_b1(n))


        else


        end if

      end do ! n

!eda:gather all of the data to master_task
           call gather( ueu_b1, ndrfmax, master_task, ueu_b)
           call gather( veu_b1, ndrfmax, master_task, veu_b)
           call gather(ulag_o1, ndrfmax, master_task,ulag_o)
           call gather(vlag_o1, ndrfmax, master_task,vlag_o)

!eda: broadcast to all processors
           call bcast(ueu_b,ndrfmax,master_task)
           call bcast(veu_b,ndrfmax,master_task)
           call bcast(ulag_o,ndrfmax,master_task)
           call bcast(vlag_o,ndrfmax,master_task)


      do i = 1, im
      do j = 1, jm

          Lx = (alon(i+1,j) - alon(i,j))*nx
          Ly = (alat(i,j+1) - alat(i,j))*ny

          usum = 0.
          vsum = 0.

          do n = 1, ndrf

            gamma = exp( -(
     1          ((lat_drf(n) - alat(i,j))/Ly)**2 +
     2          ((lon_drf(n) - alon(i,j))/Lx)**2 )/2. )

            usum = usum + gamma*(ulag_o(n) - ueu_b(n))
            vsum = vsum + gamma*(vlag_o(n) - veu_b(n))

          end do ! n

          u(i,j,1) = u(i,j,1) + beta*usum
          v(i,j,1) = v(i,j,1) + beta*vsum

      end do ! j
      end do ! i

      do i = 1, im
      do j = 1, jm

        do k = 2, kb-1

          u(i,j,k) = u(i,j,k) 
     1       + exp(zz(k)*D(i,j)/zd)*(u(i,j,1) - uold(i,j,1))
          v(i,j,k) = v(i,j,k) 
     1       + exp(zz(k)*D(i,j)/zd)*(v(i,j,1) - vold(i,j,1))

        end do

        delu = dot_product(-u(i,j,1:kb-1)+uold(i,j,1:kb-1),
     1       dz(1:kb-1))
        delv = dot_product(-v(i,j,1:kb-1)+vold(i,j,1:kb-1), 
     1       dz(1:kb-1))

        u(i,j,1:kb-1) = u(i,j,1:kb-1) + delu
        v(i,j,1:kb-1) = v(i,j,1:kb-1) + delv

      end do !j
      end do !i

      ub(:,:,1:kb-1) = u(:,:,1:kb-1)
      vb(:,:,1:kb-1) = v(:,:,1:kb-1)

!eda:
      call exchange3d_mpi(u(:,:,1:kb-1),im,jm,kb-1)
      call exchange3d_mpi(v(:,:,1:kb-1),im,jm,kb-1)
      call exchange3d_mpi(ub(:,:,1:kb-1),im,jm,kb-1)
      call exchange3d_mpi(vb(:,:,1:kb-1),im,jm,kb-1)
!
!eda:if no file, return
         else  !if fexist

                 if ( my_task == master_task ) then
                  write(*,'(/2a)') 
     $              "missing drf file at assimdrf_main : "
     $              , trim(infilee)
                 endif
            
         endif !if fexist  

      return

      end subroutine assimdrf_OIpsLag
!==================================================================
!
c
c-------------------------------------------------------------------
      SUBROUTINE BLINT(x1,x2,x3,x4,y1,y2,y3,y4,f1,f2,f3,f4,x,y,f)
C
C Bilinear interpolation subroutine.
C (Xi,Yi,fi) = data grid & values surounding model point (x,y)
C f = interpolated value at the model grid point.
C
      implicit none
      include 'realkind'

      real(kind=rk), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,
     &                             f1,f2,f3,f4,x,y
      real(kind=rk), intent(out) :: f
      real(kind=rk) :: a1,a2,a3,a4,b1,b2,b3,b4,A,B,C,t,s

      a1=x1-x2+x3-x4
      a2=-x1+x4
      a3=-x1+x2
      a4=x1-x
      b1=y1-y2+y3-y4
      b2=-y1+y4
      b3=-y1+y2
      b4=y1-y
      A=a3*b1-a1*b3
      B=b2*a3+b1*a4-a1*b4-a2*b3
      C=-a2*b4+a4*b2
      if(ABS(A*C).gt.0.002*B**2) then
         t=(-B-sqrt(B*B-4.*A*C))/(2.*A)
         else
         t=C/ABS(B)
      endif
!  10  CONTINUE
      A=a2*b1-a1*b2
      B=b3*a2+b1*a4-a1*b4-a3*b2
      C=-a3*b4+a4*b3
      if(ABS(A*C).gt.0.002*B**2) then
         s=(-B+sqrt(B*B-4.*A*C))/(2.*A)
         else
         s=-C/ABS(B)
      endif
!  20  CONTINUE
      f=f1*(1.-t)*(1.-s)+f2*t*(1.-s)+f3*s*t+f4*(1.-t)*s

      end subroutine BLINT
c
c-------------------------------------------------------------------
      subroutine spher2curvi(us, vs, lon1, lon2, lat1, lat2, uc, vc)
      ! Coordinate rotaion from spherical to curvi linear for u & v.
      ! lon1= alon(i-1,j);  lon2=alon(i,j)
      ! lat1= alat(i-1,j);  lat2=alat(i,j)

      implicit none
      include 'realkind'

      real(kind=rk), parameter :: radian = 3.14159265/180.

      real(kind=rk), intent(in) :: us, vs, lon1, lon2, lat1, lat2
      real(kind=rk), intent(out) :: uc, vc

      real(kind=rk) :: dxc, dyc, angc

      dxc = lon2 - lon1
      dyc = lat2 - lat1
      angc = atan(dyc/(dxc*cos(lat2*radian)))
      uc =   us*cos(angc) + vs*sin(angc)
      vc = - us*sin(angc) + vs*cos(angc)

      end subroutine spher2curvi

!==================================================================
      logical function ingridbox(x1,x2,x3,x4,y1,y2,y3,y4,x,y)

      implicit none
      include 'realkind'

      real(kind=rk), parameter :: pi = 3.14159

      real(kind=rk), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y
!      logical :: ingridbox
      real(kind=rk), dimension(4) :: xbox, ybox, angle
      integer :: ii

      xbox(1) = x1; xbox(2) = x2; xbox(3) = x3; xbox(4) = x4
      ybox(1) = y1; ybox(2) = y2; ybox(3) = y3; ybox(4) = y4

      do ii = 1, 4
         angle(ii) = atan2(ybox(ii) - y, xbox(ii) - x)
      end do

      if ( maxval(angle) - minval(angle) >= pi ) then
         ingridbox = .true.
      else
         ingridbox = .false.
      end if

c      write(*,*)xbox, ybox, angle, ingridbox

      end function ingridbox
!==================================================================
      subroutine bcast(data,elts,src)
      implicit none
      integer elts,src
      include 'mpif.h'

      integer ierr
      include 'pom.h'
      real(kind=rk) data(elts)
      integer mpi_rk

      if (rk==8) then
        mpi_rk = mpi_double_precision
      else
        mpi_rk = mpi_real
      end if

      call mpi_bcast(data,elts,mpi_rk,src,mpi_comm_world,ierr)

      end subroutine bcast
!==================================================================
      subroutine gather(work,nx,src,new)

      implicit none

      include 'mpif.h'
      include 'pom.h'
      integer nx,src
      integer ierr,i,j
      real(kind=rk), intent(in) ::  work(nx)
      real(kind=rk) ::  buf(nx,n_proc)
      real(kind=rk), intent(out) ::  new(nx)

      call MPI_GATHER(work,nx,MPI_REAL,
     &                 buf,nx,MPI_REAL,
     &                 src,MPI_COMM_WORLD,ierr)

!eda: select the good value after gather all of the data

       if ( my_task == src ) then
         do j = 1, n_proc
           do i = 1, nx     

        if ( abs(buf(i,j))> 1.e-5) then
            new(i)=buf(i,j)
        end if

           end do
         end do
       end if



      end subroutine gather
!==================================================================
