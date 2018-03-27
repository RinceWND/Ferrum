      module uvforce

      implicit none

      private

      public :: uvforce_init, uvforce_main

      include 'pom.h'

      logical, parameter :: read_transport = .false.
      integer, parameter :: n_bry = 3

!     days in month
      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,
     $                          31, 31, 30, 31, 30, 31/)


!      real(kind=rk), dimension( im_local ) ::
!     $     vabs_a, vabn_a, vabs_b, vabn_b
!
!      real(kind=rk), dimension( jm_local ) ::
!     $     uabe_a, uabw_a, uabe_b, uabw_b

      real(kind=rk), dimension( im_local, kb ) ::
     $     vbs_a, vbn_a, vbs_b, vbn_b

      real(kind=rk), dimension( jm_local, kb ) ::
     $     ube_a, ubw_a, ube_b, ubw_b

      real(kind=rk), dimension( 10 ) ::
     &     trans_a, trans_b, transport
      real(kind=rk), dimension( 2, 10 ) ::
     &     x_a, x_b, y_a, y_b, x, y

      integer :: mon_a, mon_b, sec_in_month, mid_in_month
      integer :: k
      character(len=7) :: infile_a, infile_b
      real(kind=rk) :: aa


      contains

!==============================================================
! Initialization variables for U & V boundary condition
!--------------------------------------------------------------
      subroutine uvforce_init( d_in )

      use module_time

      implicit none

!     intent(in)
      type(date), intent(in) :: d_in


!     check leap year
      if( ( mod( d_in%year, 4 ) .eq. 0
     $     .and. mod( d_in%year, 100 ) .ne. 0 )
     $     .or. mod( d_in%year, 400 ) .eq. 0 ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif

!     current time [sec] from the beginning of the month.

      sec_in_month = d_in%day * 24 * 3600
     $             + d_in%hour * 3600 + d_in%min * 60 + d_in%sec


!     mid-point [sec] in the month.

      mid_in_month = int( real( mday( d_in%month ) )/ 2.0
     $             * 24. * 3600. )


!     decide between which two months.

!     former half in the month.
      if ( sec_in_month .le. mid_in_month ) then
         mon_a = d_in%month - 1
!     latter half in the month
      else
         mon_a = d_in%month
      endif

      mon_b = mon_a + 1

      if ( mon_a ==  0 ) mon_a = 12
      if ( mon_b == 13 ) mon_b =  1


!     data open.
      if ( read_transport ) then

        call read_bc_transport_pnetcdf
     &      ( trans_a, y_a, x_a, "bc.trans.nc", mon_a )

        call read_bc_transport_pnetcdf
     &      ( trans_b, y_b, x_b, "bc.trans.nc", mon_b )

      else

        call read_bc_pnetcdf
     $      ( ube_a, ubw_a, vbs_a, vbn_a, "bc.nc", mon_a )
        call read_bc_pnetcdf
     $      ( ube_b, ubw_b, vbs_b, vbn_b, "bc.nc", mon_b )
!     TODO: Treat water transport increase properly in case of baroclinic velocities input
!           The below won't work for vertically non-homogeneous velocity field
        ube_a = sf_bf*ube_a
        ubw_a = sf_bf*ubw_a
        vbs_a = sf_bf*vbs_a
        vbn_a = sf_bf*vbn_a

      end if


      if ( my_task == master_task )
     $        write(*,'(/a/)') "---------- uvforce_init."


      return

      end subroutine uvforce_init
!--------------------------------------------------------------

!==============================================================
! Read & time-interpolation of U & V boundary condition
!--------------------------------------------------------------
      subroutine uvforce_main( d_in )

      use module_time

      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in



!     check leap year.
      if( ( mod( d_in%year, 4 ) .eq. 0
     $     .and. mod( d_in%year, 100 ) .ne. 0 )
     $     .or. mod( d_in%year, 400 ) .eq. 0 ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif


!     current time [sec] from the beginning of the month.

      sec_in_month = d_in%day * 24 * 3600
     $             + d_in%hour * 3600 + d_in%min * 60 + d_in%sec


!     mid-point [sec] in the month.

      mid_in_month = int( real( mday( d_in%month ) )/ 2. * 24.*3600. )


!     decide between which two months.

!     former half in the month.
      if ( sec_in_month .le. mid_in_month ) then
         mon_a = d_in%month - 1
         aa = real( sec_in_month )
     $      + real( mday( mon_a ) )/ 2. * 24. * 3600.

!     latter half in the month
      else
         mon_a = d_in%month
         aa = real( sec_in_month )
     $      - real( mday( mon_a ) )/ 2. * 24. * 3600.
      endif

      mon_b = mon_a + 1

      if ( mon_a ==  0 ) mon_a = 12
      if ( mon_b == 13 ) mon_b =  1

      aa = aa /
     $     ( real( mday( mon_a ) + mday( mon_b ) ) / 2. * 24. * 3600. )



!     data open.

      if ( sec_in_month > mid_in_month  .and.
     $     sec_in_month - mid_in_month <= int( dti ) ) then

!         uabe_a = uabe_b
!         uabw_a = uabw_b
!         vabs_a = vabs_b
!         vabn_a = vabn_b

         ube_a = ube_b
         ubw_a = ubw_b
         vbs_a = vbs_b
         vbn_a = vbn_b

        if ( read_transport ) then

          call read_bc_transport_pnetcdf
     &      ( n_bry, trans_b, y_b, x_b, "bc.trans.nc", mon_b )

        else

          write( infile_b, '( "bc",i2.2,".nc" )' ) mon_b
          call read_bc_pnetcdf
     $        ( ube_b, ubw_b, vbs_b, vbn_b, "bc.nc", mon_b )

          ube_b = sf_bf*ube_b
          ubw_b = sf_bf*ubw_b
          vbs_b = sf_bf*vbs_b
          vbn_b = sf_bf*vbn_b

        end if

      endif


!     time interpolation.

      if ( read_transport ) then

        transport = ( 1.0 - aa ) * trans_a + aa * trans_b
        y = ( 1.0 - aa ) * y_a + aa * y_b
        x = ( 1.0 - aa ) * x_a + aa * x_b

        call transport_to_velocities

!        print *, my_task, ": VOLTR [Sv]: ", transport(1:n_bry)
        if ( my_task == 1 ) then
          print *, my_task, ": TsuV [m/s]: "
     &         , minval(vabs), ";", maxval(vabs)
        end if
        if ( my_task == 3 ) then
          print *, my_task, ": TsgUSoy [m/s]: "
     &         , minval(uabe), ";", maxval(uabe)
        end if

      else

        ube = ( 1.0 - aa ) * ube_a + aa * ube_b
        ubw = ( 1.0 - aa ) * ubw_a + aa * ubw_b
        vbs = ( 1.0 - aa ) * vbs_a + aa * vbs_b
        vbn = ( 1.0 - aa ) * vbn_a + aa * vbn_b

  !     integrate vertically
        uabe = 0.
        uabw = 0.
        vabs = 0.
        vabn = 0.
        do k=1,kb
          uabe = uabe + ube(:,k)*dz(k)
          uabw = uabw + ubw(:,k)*dz(k)
          vabs = vabs + vbs(:,k)*dz(k)
          vabn = vabn + vbn(:,k)*dz(k)
        end do

      end if


      return

      end subroutine uvforce_main
!--------------------------------------------------------------

      subroutine transport_to_velocities

        implicit none

        include 'pom.h'

        integer i,j,n
        real(kind=rk) area
        integer x1,x2, y1,y2

        x1 = 0
        x2 = 0
        y1 = 0
        y2 = 0

        do n = 1, n_bry

! get boundary section area (the loop is convoluted; put only sections here, not volumes)

! -- north
          area = 0.

          if ( n_north == -1 ) then
            if ( y(1,n) == y(2,n)
     &                    .and. int(y(1,n)) == jm_global ) then

              x1 = minloc(i_global, 1, i_global >= x(1,n))
!              print *, my_task, " @N x1: ", x1
              x2 = maxloc(i_global, 1, i_global <= x(2,n))
!              print *, my_task, " @N x2: ", x2

              if ( x1 > 0 .and. x2 <= im ) then

                do i = x1, x2
                  area = area + 0.25
     &           * ( h(i,jm) + elf(i,jm) + h(i,jmm1) + elf(i,jmm1) )
     &           * ( dx(i,jm) + dx(i,jmm1) ) * dvm(i,jm)
                end do

              end if

            end if
          end if

          call sum0d_mpi( area, 0 )
          call bcast0d_mpi( area, 0 )

          if ( area > 0. ) then
            if ( n_north == -1 ) then
              if ( y(1,n) == y(2,n)
     &                    .and. int(y(1,n)) == jm_global ) then

                if ( x1 > 0 .and. x2 <= im ) then
                  vabn(x1:x2) = dvm(x1:x2,jm)*transport(n) / area
                end if
!              print *, my_task, " @N TR: ", transport(n)
!              print *, my_task, " @N mV: ", vabn(x1)

              end if
            end if
          end if

! -- east
          area = 0.

          if ( n_east == -1 ) then
            if ( x(1,n) == x(2,n)
     &                    .and. int(x(1,n)) == im_global ) then

              y1 = minloc(j_global, 1, j_global >= y(1,n))
!              print *, my_task, " @E y1: ", y1
              y2 = maxloc(j_global, 1, j_global <= y(2,n))
!              print *, my_task, " @E y2: ", y2

              if ( y1 > 0 .and. y2 <= jm ) then

                do j = y1, y2
                  area = area + 0.25
     &           * ( h(im,j) + elf(im,j) + h(imm1,j) + elf(imm1,j) )
     &           * ( dy(im,j) + dy(imm1,j) ) * dum(im,j)
                end do

              end if

            end if
          end if

          call sum0d_mpi( area, 0 )
          call bcast0d_mpi( area, 0 )

          if ( area > 0. ) then
            if ( n_east == -1 ) then
              if ( x(1,n) == x(2,n)
     &                    .and. int(x(1,n)) == im_global ) then

                if ( y1 > 0 .and. y2 <= jm ) then
                  uabe(y1:y2) = dum(im,y1:y2)*transport(n) / area
                end if
!              print *, my_task, " @E TR: ", transport(n)
!              print *, my_task, " @E mV: ", uabe(y1)

              end if
            end if
          end if

! -- south
          area = 0.

          if ( n_south == -1 ) then
            if ( y(1,n) == y(2,n) .and. int(y(1,n)) == 1 ) then

              x1 = minloc(i_global, 1, i_global >= x(1,n))
              x2 = maxloc(i_global, 1, i_global <= x(2,n))

              if ( x1 > 0 .and. x2 <= im ) then

                do i = x1, x2
                  area = area + 0.25
     &           * ( h(i,1) + elf(i,1) + h(i,2) + elf(i,2) )
     &           * ( dx(i,1) + dx(i,2) ) * dum(i,1)
                end do

              end if

            end if
          end if

          call sum0d_mpi( area, 0 )
          call bcast0d_mpi( area, 0 )

          if ( area > 0. ) then
            if ( n_south == -1 ) then
              if ( y(1,n) == y(2,n) .and. int(y(1,n)) == 1 ) then

                if ( x1 > 0 .and. x2 <= im ) then
                  vabs(x1:x2) = dvm(x1:x2,1)*transport(n) / area
                end if

              end if
            end if
          end if

! -- west
          area = 0.

          if ( n_west == -1 ) then
            if ( x(1,n) == x(2,n) .and. int(x(1,n)) == 1 ) then

              y1 = minloc(j_global, 1, j_global >= y(1,n))
!              print *, my_task, " @W y1: ", y1
              y2 = maxloc(j_global, 1, j_global <= y(2,n))
!              print *, my_task, " @W y2: ", y2

              if ( y1 > 0 .and. y2 <= jm ) then

                do j = y1, y2
                  area = area + 0.25
     &           * ( h(1,j) + elf(1,j) + h(2,j) + elf(2,j) )
     &           * ( dy(1,j) + dy(2,j) ) * dum(1,j)
                end do

              end if

            end if
          end if

          call sum0d_mpi( area, 0 )
          call bcast0d_mpi( area, 0 )

          if ( area > 0. ) then
            if ( n_west == -1 ) then
              if ( x(1,n) == x(2,n)
     &                    .and. int(x(1,n)) == im_global ) then

                if ( y1 > 0 .and. y2 <= jm ) then
                  uabw(y1:y2) = dum(1,y1:y2)*transport(n) / area
                end if
!              print *, my_task, " @W TR: ", transport(n)
!              print *, my_task, " @W mV: ", uabw(y1)

              end if
            end if
          end if

!        print *, my_task, "DONE!", n
        end do

      end subroutine transport_to_velocities
!--------------------------------------------------------------

      end module uvforce
