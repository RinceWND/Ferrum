      module tsforce

      implicit none

      private

      public :: tsforce_init, tsforce_main, tsforce_tsflx

      include 'pom.h'
      
!     days in month
      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,               
     $                          31, 31, 30, 31, 30, 31/)
                             
!     coefficients for relaxation 
!     c1 = 1/30 [m/day]
!     real(kind=rk), parameter :: c1 = 3.858024691e-7 
!     c1 = 1 [m/day] !lyo:20110202:
      real(kind=rk), parameter :: c1 = 1.157407407e-5 !lyo:pac10:exp001->007;exp302:
      real(kind=rk) :: sstrelx, sssrelx


!     buffers for tsurf, ssurf etc
      real(kind=rk), dimension( im_local, jm_local ) ::
     $     tsurf_a, ssurf_a, tsurf_b, ssurf_b, uht, swr, tair, emp
  
      real(kind=rk), dimension( im_local, jm_local, kb ) ::
     $     tc_a, sc_a, tc_b, sc_b,
     $     tm_a, sm_a, tm_b, sm_b,
     $     tmean, smean, rm_a, rm_b
     
      real(kind=rk), dimension( im_local, jm_local ) ::
     $     uw_a, vw_a, uw_b, vw_b
    
      integer :: mon_a, mon_b, sec_in_month, mid_in_month
      integer :: i, j, k, nb
      character*13 :: infile_a, infile_b
      real(kind=rk) :: aa


      contains

!==============================================================
! Initialization variables for TS climatorology.
!--------------------------------------------------------------
      subroutine tsforce_init( d_in )

      use module_time

      implicit none

!     intent(in)
      type(date), intent(in) :: d_in 
  

!     initialize

      wtsurf = 0.0
      wssurf = 0.0
      nb = 0

!     Read backward and forward TS climatology in initial state

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
     $             * 24 * 3600 )


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
         call read_tsclim_monthly_pnetcdf
     $        ( tm_a, sm_a, tc_a, sc_a, "ts_clim.nc", mon_a )

         call read_tsclim_monthly_pnetcdf
     $        ( tm_b, sm_b, tc_b, sc_b, "ts_clim.nc", mon_b )
     
!         call read_tsclim_monthly_pnetcdf_obs
!     $        ( rm_a, tc_a, sc_a, "ts_clim_old.nc", mon_a )
!
!         call read_tsclim_monthly_pnetcdf_obs
!     $        ( rm_b, tc_b, sc_b, "ts_clim_old.nc", mon_b )
     
!         call read_wind_monthly_pnetcdf
!     $        ( uw_a, vw_a,       "sfrc_old.nc",    mon_a )
!
!         call read_wind_monthly_pnetcdf
!     $        ( uw_b, vw_b,       "sfrc_old.nc",    mon_b )

!         write( infile_a, '( "ssts",i2.2,".nc" )' ) mon_a
!         call read_ssts_monthly_pnetcdf
!     $        ( tsurf_a, ssurf_a, trim(infile_a) )

!         write( infile_b, '( "ssts",i2.2,".nc" )' ) mon_b
!         call read_ssts_monthly_pnetcdf
!     $        ( tsurf_b, ssurf_b, trim(infile_b) )

!          tsurf_a = tm_a(:,:,1)
!          tsurf_b = tm_b(:,:,1)
!          ssurf_a = sm_a(:,:,1)
!          ssurf_b = sm_b(:,:,1)
          tsurf_a = tc_a(:,:,1)
          tsurf_b = tc_b(:,:,1)
          ssurf_a = sc_a(:,:,1)
          ssurf_b = sc_b(:,:,1)


      if ( my_task == master_task ) 
     $        write(*,'(/a/)') "---------- tsforce_init."


      return

      end subroutine tsforce_init
!--------------------------------------------------------------

!==============================================================
! Read & time-interpolation of TS climatorology.
!--------------------------------------------------------------
      subroutine tsforce_main( d_in )

      use module_time

      implicit none
      
      integer :: ii, jj
      
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
      
      mid_in_month = int( real( mday( d_in%month ) )/ 2.0  
     $             * 24 * 3600 )


!     decide between which two months.
 
!     former half in the month.    
      if ( sec_in_month .le. mid_in_month ) then
         mon_a = d_in%month - 1
         aa = real( sec_in_month ) 
     $      + real( mday( mon_a ) )/ 2.0 * 24. * 3600.

!     latter half in the month
      else
         mon_a = d_in%month 
         aa = real( sec_in_month ) 
     $      - real( mday( mon_a ) )/ 2.0 * 24. * 3600.
      endif
      
      mon_b = mon_a + 1

      if ( mon_a ==  0 ) mon_a = 12
      if ( mon_b == 13 ) mon_b =  1

      aa = aa / 
     $     ( real( mday( mon_a ) + mday( mon_b ) ) / 2.0 * 24. * 3600. )            



!     data open.

      if ( sec_in_month > mid_in_month  .and.
     $     sec_in_month - mid_in_month <= int( dti ) ) then

         tm_a = tm_b
         sm_a = sm_b
         tc_a = tc_b
         sc_a = sc_b
         uw_a = uw_b
         vw_a = vw_b

         tsurf_a = tsurf_b
         ssurf_a = ssurf_b

!         write( infile_b, '( "tsclimib",i2.2,".nc" )' ) mon_b
         call read_tsclim_monthly_pnetcdf
     $        ( tm_b, sm_b, tc_b, sc_b, "ts_clim.nc", mon_b )
!         call read_tsclim_monthly_pnetcdf_obs
!     $        ( rm_b, tc_b, sc_b, "ts_clim_old.nc", mon_b )
!         call read_wind_monthly_pnetcdf
!     $        ( uw_b, vw_b,       "sfrc_old.nc",    mon_b )

!         write( infile_b, '( "ssts",i2.2,".nc" )' ) mon_b
!         call read_ssts_monthly_pnetcdf
!     $        ( tsurf_b, ssurf_b, "ts_clim.nc", mon_b )
     
         tsurf_b = tc_b(:,:,1)
         ssurf_b = sc_b(:,:,1)

      endif

     

!     time interpolation.

      tclim = ( 1.0 - aa ) * tc_a + aa * tc_b
      sclim = ( 1.0 - aa ) * sc_a + aa * sc_b
!      tmean = ( 1.0 - aa ) * tm_a + aa * tm_b
!      smean = ( 1.0 - aa ) * sm_a + aa * sm_b
      rmean = ( 1.0 - aa ) * rm_a + aa * rm_b
      
!      wusurf = ( 1.0 - aa ) * uw_a + aa * uw_b
!      wvsurf = ( 1.0 - aa ) * vw_a + aa * vw_b

      tsurf = ( 1.0 - aa ) * tsurf_a + aa * tsurf_b
      ssurf = ( 1.0 - aa ) * ssurf_a + aa * ssurf_b
      

!     calculation of rmean.
      
!      call dens( smean, tmean, rmean )

!     set boundary condition.

!----------------------------------------------------------------------!
!lyo:20110224:alu:stcc:
!     Use tob* and sob* instead - see below commented lines 
!     to be deleted in future - will get rid of tbe,sbe,tbw,sbw,tbn,sbn,
!     tbs & sbs througout the code
!     do k=1,kb
!        do j=1,jm
!           tbe( j, k ) = tclim( im, j, k ) * fsm( im, j )
!           sbe( j, k ) = sclim( im, j, k ) * fsm( im, j )
!           tbw( j, k ) = tclim(  1, j, k ) * fsm(  1, j )
!           sbw( j, k ) = sclim(  1, j, k ) * fsm(  1, j )
!         enddo
!        do i=1,im
!           tbn( i, k ) = tclim( i, jm, k ) * fsm( i, jm )
!           sbn( i, k ) = sclim( i, jm, k ) * fsm( i, jm )
!           tbs( i, k ) = tclim( i,  1, k ) * fsm( i,  1 )
!           sbs( i, k ) = sclim( i,  1, k ) * fsm( i,  1 )
!        enddo
!     enddo
!----------------------------------------------------------------------!
      do k=1,kb
         do j=1,jm
            do i=1,nfw
               tobw(i,j,k) = tclim( i, j, k) * fsm( i, j)
               sobw(i,j,k) = sclim( i, j, k) * fsm( i, j)
               enddo
            do i=1,nfe
               ii=im-i+1
               tobe(i,j,k) = tclim(ii, j, k) * fsm(ii, j)
               sobe(i,j,k) = sclim(ii, j, k) * fsm(ii, j)
               enddo
               tbw(j,k) = tobw(1,j,k); sbw(j,k) = sobw(1,j,k)
               tbe(j,k) = tobe(1,j,k); sbe(j,k) = sobe(1,j,k)
          enddo
         do i=1,im
            do j=1,nfs
               tobs(i,j,k) = tclim( i, j, k) * fsm( i, j)
               sobs(i,j,k) = sclim( i, j, k) * fsm( i, j)
               enddo
            do j=1,nfn
               jj=jm-j+1
               tobn(i,j,k) = tclim( i,jj, k) * fsm( i,jj)
               sobn(i,j,k) = sclim( i,jj, k) * fsm( i,jj)
               enddo
               tbs(i,k) = tobs(i,1,k); sbs(i,k) = sobs(i,1,k)
               tbn(i,k) = tobn(i,1,k); sbn(i,k) = sobn(i,1,k)
          enddo
      enddo
!----------------------------------------------------------------------!

      return

      end subroutine tsforce_main
!--------------------------------------------------------------

!--------------------------------------------------------------
! Calculate heat/salt fluxes wtsurf and wssurf.
!--------------------------------------------------------------
      subroutine tsforce_tsflx( d_in )

!     surface heat/salt fluxes.

      use module_time

      implicit none
      
      integer n
      
      ! intent(in)
      type(date), intent(in) :: d_in
      type(date) d_tmp
      real(kind=rk) sec_in_day

      logical :: lexist    

 

      sec_in_day = d_in%hour*3600 + d_in%min*60 + d_in%sec

      d_tmp = str2date("1979-01-01 00:00:00")
      d_tmp%year = d_in%year
      n = int(dif_date(d_in, d_tmp)/(86400.)*4.)+1
      
      do j=1,jm
         do i=1,im
            
            sstrelx = 1.0d0
!lyo:20110202:
!           sssrelx = ( 1.0d0 + tanh( 0.002d0 * (h(i,j)-1000.d0)))*0.5d0 
!           sssrelx = ( 1.0d0 + tanh(0.0005d0 * (h(i,j)-2500.d0)))*0.5d0 

            wtsurf( i, j ) 
     $           = c1 * sstrelx * ( tb( i, j, 1 ) - tsurf( i, j ) )
            wssurf( i, j ) 
     $           = c1 * sssrelx * ( sb( i, j, 1 ) - ssurf( i, j ) )
            
         enddo
      enddo

      if (n/=nb) then
        write(*,'(a,i5)') "Reading heat record ",n
        nb = n
        write( infile_b, '( a3,".",i4.4,".nc" )' )
     $        "hfl", d_in%year

        inquire(file='in/heat/'//trim(infile_b),
     $           exist=lexist)

        swrad  = 0.
        if(lexist) then
          call read_heat_pnetcdf(
     $                  uht,swr,tair,emp,infile_b,n)
!     $                 wtsurf(1:im,1:jm),swrad(1:im,1:jm),infile_b,n)
          wtsurf(1:im,1:jm) = wtsurf+uht+emp*(tair-t(1:im,1:jm,1))*3986.
          swrad(1:im,1:jm)  = swr
          wtsurf = wtsurf/(rhoref*3986.)
          swrad  = -swrad/(rhoref*3986.)
        else
          if ( my_task == master_task ) then
            write(*,'(/2a)') 
     $      "missing heat data at wind_main : "
     $              , trim(infile_b)
          endif
        endif

      end if

! TODO: interpolation

      return
      
      end subroutine tsforce_tsflx
!--------------------------------------------------------------


      end module tsforce
