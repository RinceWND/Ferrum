
      module mcsst

      use module_time

      implicit none

      private

      public :: mcsst_init,mcsst_main             

      include 'pom.h'

      real(kind=rk), parameter :: errval = -999.

      type(date) :: d1, d2, dstart, d_buf

      character(len=19) ::
     $     dstart_str = "1992-01-22_00:00:00"

      character(len=30) :: infile

      logical, dimension( im_local, jm_local ) ::
     $     mask_mcsst

      real(kind=rk), dimension( im_local, jm_local ) :: 
     $     mcsst_a, mcsst_b, mcsst_surf   

      integer ::  i, j, ndays

      logical :: lexist

      integer :: count_days            !fhx:debug

      contains

!==============================================================
!--------------------------------------------------------------
!fhx:mcsst:initialize mcsst
!--------------------------------------------------------------
      subroutine mcsst_init( d_in )

      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in
   
         dstart = str2date( dstart_str )

         if ( dstart > d_in ) then

            if ( my_task == master_task ) then
               write(*,'(/2a)') 
     $              "Before MCSST period. d_in : "
     $              , date2str(d_in)
               write(*,'(2a)') "no MCSST until "
     $              , dstart_str
            endif
            
            d1 = dstart
            d2 = dstart

         else

      !     decide the backward date.

      ! set initial d_buf just before or at the d_in date
 
            d_buf = dstart + 
     $    max(( int(( d_in - dstart ) / ( 86400 ) ) ) 
     $    * ( 86400 ), 0 ) 

           lexist = .FALSE.
           ndays = 1

           do while ( .NOT.lexist )

           write( infile, 
     $              '( "./in/mcsst/mcsst_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day

            inquire(file=trim(infile),exist=lexist) 

               if ( lexist ) then
                  d1 = d_buf     
               else
!                  if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing MCSST data at mcsst_init : "
!     $              , trim(infile)
!                  endif                     
               d_buf = dstart + 
     $    ( int(( d_in - dstart ) / ( 86400 )-ndays ) ) * ( 86400 )
               ndays = ndays + 1
               endif
           end do

      !     decide the forward date.
             
            d_buf = d1 + 24 * 3600
            lexist = .FALSE.
            count_days = 0      !fhx:debug
            do while ( .NOT.lexist )
               write( infile, 
     $              '( "./in/mcsst/mcsst_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
              inquire(file=trim(infile),exist=lexist)
  
               if ( lexist ) then 
                   d2 = d_buf
               else
!                 if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing MCSST data at mcsst_init : "
!     $              , trim(infile)
!                 endif

               count_days = count_days+1     !fhx:debug
               if(count_days>12) then   !fhx:debug
                 calc_tsurf_mc = .FALSE.
                 go to 100
               endif
                 
                 d_buf = d_buf + 24 * 3600
               endif       
              
            end do

         endif

!     read data

         write( infile, '( "mcsst_",i4.4,2i2.2,".nc" )' ) 
     $        d1%year, d1%month, d1%day

         call read_mcsst_pnetcdf( mcsst_a, trim(infile) )    
       
         write( infile, '( "mcsst_",i4.4,2i2.2,".nc" )' ) 
     $        d2%year, d2%month, d2%day

         call read_mcsst_pnetcdf( mcsst_b, trim(infile) )

       
         if ( my_task == master_task ) then
            write(*,'(/4a)')  "MCSST backward & forward : "
     $           , date2str(d1)," ", date2str(d2)
            write(*,'(/a/)')  "---------- mcsst_init. "
         endif            

100     return
      end subroutine mcsst_init
!--------------------------------------------------------------
!fhx:mcsst:read and interpolate mcsst
!--------------------------------------------------------------
      subroutine mcsst_main( d_in )

      implicit none

      ! intent(in)
      type(date), intent(in) :: d_in
    
      ! local
      integer :: sec_d1din, sec_d1d2    
      real(kind=rk) :: aa

      if ( dstart > d_in ) return

      if ( d_in >= d2 ) then     !  renew the backward and forward data

      ! set initial d_buf just before or at the d_in date
 
            d_buf = dstart + 
     $    max(( int(( d_in - dstart ) / ( 86400 ) ) ) 
     $    * ( 86400 ), 0 ) 
     
            lexist = .FALSE.
            do while ( .NOT.lexist )
               write( infile, 
     $              '( "./in/mcsst/mcsst_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               inquire(file=trim(infile),exist=lexist)         
               if ( lexist ) then
                  d1 = d_buf
               else
!                  if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing MCSST data at mcsst_main : "
!     $              , trim(infile)
!                 endif                      
               d_buf = dstart + 
     $    ( int(( d_in - dstart ) / ( 86400 )-1 ) ) * ( 86400 ) 
               endif

            end do

      !     decide the forward date.

            d_buf = d1 + 24 * 3600
            lexist = .FALSE.            
            count_days = 0      !fhx:debug

            do while (  .NOT.lexist )
               write( infile, 
     $              '( "./in/mcsst/mcsst_",i4.4,2i2.2,".nc" )' ) 
     $              d_buf%year, d_buf%month, d_buf%day
               inquire(file=trim(infile),exist=lexist)         
               if ( lexist ) then 
                   d2 = d_buf
               else
!                 if ( my_task == master_task ) then
!                   write(*,'(/2a)') 
!     $              "missing MCSST data at mcsst_main : "
!     $              , trim(infile)
!                 endif
               count_days = count_days+1     !fhx:debug
               if(count_days>12) then   !fhx:debug
                 calc_tsurf_mc = .FALSE.
                 go to 200
               endif
                 d_buf = d_buf + 24 * 3600
               endif       
              
            end do

!     read data
     
         write( infile, '( "mcsst_",i4.4,2i2.2,".nc" )' ) 
     $        d1%year, d1%month, d1%day

         call read_mcsst_pnetcdf( mcsst_a, trim(infile) )    

         write( infile, '( "mcsst_",i4.4,2i2.2,".nc" )' ) 
     $        d2%year, d2%month, d2%day

         call read_mcsst_pnetcdf( mcsst_b, trim(infile) )

         endif    ! if ( d_in >= d2 ) then  

!interpolation in time
!     seconds from the previous data date.

      sec_d1din = d_in - d1 

!     time interval between the data.    
    
      sec_d1d2 = d2 -d1

      aa = real( sec_d1din ) / real( sec_d1d2 )

!     error mask
     
      where ( mcsst_a == errval .or. mcsst_b == errval ) 
         mask_mcsst = .false.
         mcsst_surf = 0.0
      elsewhere
         mask_mcsst = .true. 
         mcsst_surf = ( 1.0 - aa ) * mcsst_a + aa * mcsst_b
      endwhere

      tsurf = mcsst_surf

      if ( my_task == n_proc/2-1 ) then 
         write(*,'(/3(a,f12.5),a,2i4)') 
     $        "mcsst done: tb@k1 = ", tb(im/2,jm/2,1), 
     $        ", tsurf = ", tsurf(im/2,jm/2),
     $        ", mcsst_surf = ", mcsst_surf(im/2,jm/2),
     $        ", (i,j) = ", i_global(im/2),j_global(jm/2)
         
         write(*,'(/a/)')  "---------- mcsst_main. "

      endif

200   return

      end subroutine mcsst_main
!--------------------------------------------------------------
      end module mcsst
