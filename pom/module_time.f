      module module_time

      implicit none

      type date
      integer :: year, month, day, hour, min, sec
      end type date

      interface operator(+)
      module procedure add_date
      end interface

      interface operator(-)
      module procedure dif_date
      end interface

      interface operator(>=)
      module procedure ge_date
      end interface

      interface operator(>)
      module procedure gt_date
      end interface

      contains


!______________________________________________________________________
!
      logical function ge_date( d1, d2 )
!----------------------------------------------------------------------
!  Returns whether d1 is equal to or earlier than d2.
!______________________________________________________________________
!
      type(date), intent(in) :: d1, d2


      ge_date = lge( date2str( d1 ), date2str( d2 ) )


      end ! function ge_date
!
!______________________________________________________________________
!
      logical function gt_date( d1, d2 )
!----------------------------------------------------------------------
!  Returns whether d1 is earlier than d2.
!______________________________________________________________________
!
      type(date), intent(in) :: d1, d2


      gt_date = lgt( date2str( d1 ), date2str( d2 ) )


      end ! function gt_date
!
!______________________________________________________________________
!
      type(date) function str2date( str ) result(d)
!----------------------------------------------------------------------
!  Converts datetime string (`YYYY-MM-DD*hh*ii*ss` format) to `date`
! type. `*` denotes any (ignored) character.
!______________________________________________________________________
!
      character(len=19) :: str

      read(str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                    
     $     d%year, d%month, d%day, d%hour, d%min, d%sec


      end ! function str2date
!
!______________________________________________________________________
!
      character(19) function date2str( d ) result(str)
!----------------------------------------------------------------------
!  Converts `date` type variable d to datetime string in
! `YYYY-MM-DD_hh_ii_ss` format.
!______________________________________________________________________
!
      type(date) :: d

      write(str, 
     $     '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,"_",i2.2,"_",i2.2)' )
     $     d%year, d%month, d%day, d%hour, d%min, d%sec


      end ! function date2str
!
!______________________________________________________________________
!
      type(date) function add_date( d, dt_in ) result(d_result)
!----------------------------------------------------------------------
!  Adds `dt_in` (s) to date `d`.
!______________________________________________________________________
!
      type(date), intent(in) :: d
      integer, intent(in) :: dt_in

      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,
     $                          31, 31, 30, 31, 30, 31/)


!     dt_in must be positive.
      if  ( dt_in < 0 ) then
        write(*,'(a,i0)')
     $  " Eror. module_time : dt_in must be positive. dt_in = ", dt_in
        stop
      end if

!     correct for leap year
      mday(2) = 28 + inc_leap(d%year)


!     seconds
      d_result%sec = d%sec + dt_in
      d_result%min = d_result%sec / 60
      d_result%sec = mod( d_result%sec, 60 )

!     minutes
      d_result%min  = d_result%min + d%min
      d_result%hour = d_result%min / 60
      d_result%min  = mod( d_result%min, 60 )

!     hour
      d_result%hour = d_result%hour + d%hour
      d_result%day  = d_result%hour / 24
      d_result%hour = mod( d_result%hour, 24 )

      d_result%day = d_result%day + d%day
      d_result%month = d%month
      d_result%year = d%year

      do while ( d_result%day > mday(d_result%month) )
        d_result%day   = d_result%day   - mday(d_result%month)
        d_result%month = d_result%month + 1
        if ( d_result%month == 13 ) then
          d_result%year  = d_result%year + 1
          d_result%month = 1
          mday(2) = 28 + inc_leap(d_result%year)
        end if
      end do


      end ! function add_date
!
!______________________________________________________________________
!
      integer function dif_date( d1, d2 )
!----------------------------------------------------------------------
!  Calculates time interval (s) from d1 to d2 ( d1 > d2 ).
!______________________________________________________________________
!
      type(date), intent(in) :: d1, d2

      integer :: mday(12) = (/ 31, 28, 31, 30, 31, 30,        
     $                         31, 31, 30, 31, 30, 31 /)

      integer :: day1, day2, sec1, sec2, iy


!     Number of days since 1 January d2%year.

      mday(2) = 28 + inc_leap(d1%year)
      day1 = sum( mday(1:d1%month) ) - mday(d1%month) + d1%day
      mday(2) = 28 + inc_leap(d2%year)
      day2 = sum( mday(1:d2%month) ) - mday(d2%month) + d2%day

      do iy = d2%year, d1%year - 1
        day1 = day1 + 365 + inc_leap(iy)
      end do

!     Number of seconds since 0h in each day.

      sec1 = 3600 * d1%hour + 60 * d1%min + d1%sec
      sec2 = 3600 * d2%hour + 60 * d2%min + d2%sec

!     sum diffs

      dif_date = 86400 * (day1 - day2) + (sec1 - sec2)


      end ! function dif_date
!
!______________________________________________________________________
!
      pure integer function inc_leap( year )
!----------------------------------------------------------------------
!  Returns a one day increment if `year` is leap year.
!______________________________________________________________________
!
      integer, intent(in) :: year


      if ( is_leap( year ) ) then
        inc_leap = 1
      else
        inc_leap = 0
      end if


      end ! function inc_leap
!
!______________________________________________________________________
!
      pure logical function is_leap( year )
!----------------------------------------------------------------------
!  Returns a one day increment if `year` is leap year.
!______________________________________________________________________
!
      integer, intent(in) :: year


      if ( ( mod(year, 4  ) == 0 .and. mod(year, 100) /= 0 ) .or.
     $       mod(year, 400) == 0 ) then
        is_leap = .true.
      else
        is_leap = .false.
      end if


      end ! function is_leap
!
!______________________________________________________________________
!
      integer function seconds_of_year( d )
!----------------------------------------------------------------------
!  Counts seconds since the start of the year for date `d`.
!______________________________________________________________________
!
        implicit none

        type(date), intent(in) :: d

        type(date) d0


        d0 = str2date("1979-01-01 00:00:00")
        d0%year = d%year

        seconds_of_year = d - d0


      end ! function seconds_of_year
!
!______________________________________________________________________
!
      integer function interval_of_year( d, i )
!----------------------------------------------------------------------
!  Counts the number of full time intervals `i` (s) for date `d`
! since the start of the year `d%year`.
!______________________________________________________________________
!
        implicit none

        type(date), intent(in) :: d
        integer   , intent(in) :: i

        type(date) d0


        d0 = str2date("1979-01-01 00:00:00")
        d0%year = d%year

        interval_of_year = ( d - d0 ) / i


      end ! function interval_of_year
!
!______________________________________________________________________
!
      real function chunk_of_year( d, i )
!----------------------------------------------------------------------
!  Counts the time since the start of the year for date `d`
! in chunks of `i`.
!______________________________________________________________________
!
        implicit none

        type(date), intent(in) :: d
        integer   , intent(in) :: i

        type(date) d0


        d0 = str2date("1979-01-01 00:00:00")
        d0%year = d%year

        chunk_of_year = ( d - d0 ) / real(abs(i))


      end ! function chunk_of_year
!
!______________________________________________________________________
!
      integer function max_chunks_in_year( y, i )
!----------------------------------------------------------------------
!  Gets the maximum time since the start of the year `y`
! in chunks of `i`.
!______________________________________________________________________
!
        implicit none

        integer, intent(in) :: y, i

        type(date) d0, d1


        d0 = str2date("1979-01-01 00:00:00")
        d1 = d0
        d0%year = y
        d1%year = y+1

        max_chunks_in_year = ( d1 - d0 ) / abs(i)


      end ! function max_chunks_in_year
!
!______________________________________________________________________
!
      end module
