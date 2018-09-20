      module tsforce

      use air        , only: wssurf, wtsurf
      use glob_const , only: rk
      use glob_domain, only: im, jm, kb
      use glob_grid  , only: fsm
      use glob_ocean , only: rho, s, sb, sclim, t, tb, tclim, u, v
      use wind

      implicit none

      private

      public :: tsforce_init, tsforce_main, tsforce_tsflx

      real(rk), external :: heatlat

      logical, parameter :: read_bulk     = .true.
     &                    , read_flux     = .false.
     &                    , relax_surface = .true.
     &                    , corr_surface  = .false.
!     days in month
      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,
     $                          31, 31, 30, 31, 30, 31/)

!     coefficients for relaxation
!     c1 = 1/30 [m/day]
      real(rk), parameter :: c1 = 3.858024691e-7
!     c1 = 1 [m/day] !lyo:20110202:
!      real(kind=rk), parameter :: c1 = 1.157407407e-5 !lyo:pac10:exp001->007;exp302:
      real(rk) :: sstrelx, sssrelx


!     buffers for tsurf, ssurf etc
      real(rk), allocatable, dimension(:,:) ::
     &     el_a, el_b
     &    ,cloud_a, cloud_b, cloud
     &    ,rain_a, rain_b, rain
     &    ,shum_a, shum_b, shum
     &    ,tair_a, tair_b, tair
!     &    ,tsurf_a, ssurf_a, tsurf_b, ssurf_b
     &    ,pres_a, pres_b, pres
     &    ,uwnd_a, vwnd_a, uwnd_b, vwnd_b, uwnd, vwnd
     &    ,tskin_a, tskin_b, tskin, sskin, emp
     &    ,sh,sh_a,sh_b, lh,lh_a,lh_b
     &    ,sr,sr_a,sr_b, lr,lr_a,lr_b
     &    ,dlr,dlr_a,dlr_b
     &    ,dtemp,dt_a,dt_b, dsalt,ds_a,ds_b

      real(rk), allocatable, dimension(:,:,:) ::
     $     tc_a, sc_a, tc_b, sc_b,
     $     tm_a, sm_a, tm_b, sm_b,
     $     tmean, smean

!      real(kind=rk), dimension( im_local, jm_local ) ::
!     $     uw_a, vw_a, uw_b, vw_b

      integer :: mon_a, mon_b, sec_in_month, mid_in_month
      integer :: i, j, nb, mb, db
      character(len=14) :: infile_b
      real(rk) aa, bb
      real(rk) days_in_year

      real(rk), parameter :: rho_fw = 1000.

      contains

!==============================================================
! Initialization variables for TS climatorology.
!--------------------------------------------------------------
      subroutine tsforce_init( d_in )

      use glob_domain, only: is_master
      use module_time

      implicit none

!     intent(in)
      type(date), intent(in) :: d_in
      type(date) d_tmp
      integer n


! Allocate arrays
      allocate(
     &  el_a(im,jm),el_b(im,jm),cloud_a(im,jm),cloud_b(im,jm)
     &, cloud(im,jm),rain_a(im,jm),rain_b(im,jm),rain(im,jm)
     &, shum_a(im,jm),shum_b(im,jm),shum(im,jm),tair_a(im,jm)
     &, tair_b(im,jm),tair(im,jm),pres_a(im,jm),pres_b(im,jm)
     &, pres(im,jm),uwnd_a(im,jm),vwnd_a(im,jm),uwnd_b(im,jm)
     &, vwnd_b(im,jm),uwnd(im,jm),vwnd(im,jm),tskin_a(im,jm)
     &, tskin_b(im,jm),tskin(im,jm),sskin(im,jm),emp(im,jm)
     &, sh(im,jm),sh_a(im,jm),sh_b(im,jm),lh(im,jm),lh_a(im,jm)
     &, lh_b(im,jm),sr(im,jm),sr_a(im,jm),sr_b(im,jm)
     &, lr(im,jm),lr_a(im,jm),lr_b(im,jm),dlr(im,jm)
     &, dlr_a(im,jm),dlr_b(im,jm),dtemp(im,jm),dt_a(im,jm)
     &, dt_b(im,jm),dsalt(im,jm),ds_a(im,jm),ds_b(im,jm)
     &, tc_a(im,jm,kb),tc_b(im,jm,kb),sc_a(im,jm,kb)
     &, sc_b(im,jm,kb),tm_a(im,jm,kb),tm_b(im,jm,kb)
     &, sm_a(im,jm,kb),sm_b(im,jm,kb),tmean(im,jm,kb)
     &, smean(im,jm,kb)
     & )


!     initialize

      wtsurf = 0.0
      wssurf = 0.0
      nb = 0
      mb = 0
      db = 0
      pres_a = 1013.
      pres_b = 1013.
      cloud_a = 0.
      cloud_b = 0.
      shum_a = 0.
      shum_b = 0.
      rain_a = 0.
      rain_b = 0.
      tair_a = 15.
      tair_b = 15.
      dtemp = 0.
      dt_a = 0.
      dt_b = 0.
      dsalt = 0.
      sr = 0.
      uwnd_a = 0.
      uwnd_b = 0.
      vwnd_a = 0.
      vwnd_b = 0.
      tskin_a = 15.
      tskin_b = 15.


!     Read backward and forward TS climatology in initial state

!     check leap year
      if(    ( mod( d_in%year,   4 ) == 0
     $   .and. mod( d_in%year, 100 ) /= 0 )
     $    .or. mod( d_in%year, 400 ) == 0   ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif

!     current time [sec] from the beginning of the month.

      sec_in_month = d_in%day * 24 * 3600
     $             + d_in%hour * 3600 + d_in%min * 60 + d_in%sec

      days_in_year = real( mday( d_in%month ) - 31 + d_in%day )

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
         call read_tsclim_monthly_pnetcdf
     $        ( tm_a, sm_a, tc_a, sc_a, el_a, "ts_clim.nc", mon_a )

         call read_tsclim_monthly_pnetcdf
     $        ( tm_b, sm_b, tc_b, sc_b, el_b, "ts_clim.nc", mon_b )

      if ( corr_surface ) then
        call read_surf_corr_pnetcdf( dt_a, ds_a, "dt.nc", mon_a )
        call read_surf_corr_pnetcdf( dt_b, ds_b, "dt.nc", mon_b )
      end if
!      call read_mean_ts_z_pnetcdf(tm_b, sm_b, 44, int(days_in_year))

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
!          tsurf_a = tc_a(:,:,1)
!          tsurf_b = tc_b(:,:,1)
!          ssurf_a = sc_a(:,:,1)
!          ssurf_b = sc_b(:,:,1)
      d_tmp = str2date("1979-01-01 00:00:00")
      d_tmp%year = d_in%year
      n = int((d_in-d_tmp)/86400.*4.)+1

      if ( read_bulk ) then

        write( infile_b, '( i4.4,".nc" )' ) d_in%year

        call read_bulk_pnetcdf( pres_a,  tair_a, shum_a
     &                        , rain_a, cloud_a, uwnd_a
     &                        , vwnd_a, tskin_a, infile_b, n)

      end if

      if ( read_flux ) then

        write( infile_b, '( i4.4 )' ) d_in%year

        call read_heat_pnetcdf( sh_a,  lh_a, sr_a
     &                        , lr_a, dlr_a, infile_b, n )

      end if

! Read next record
      n = n+1
      if (n > 4*(365+inc_leap(d_in%year))) then
        n = mod(n,4*(365+inc_leap(d_in%year)))
        d_tmp%year = d_tmp%year+1
      end if

      if ( read_bulk ) then

        write( infile_b, '( i4.4,".nc" )' ) d_in%year

        call read_bulk_pnetcdf( pres_b,  tair_b, shum_b
     &                        , rain_b, cloud_b, uwnd_b
     &                        , vwnd_b, tskin_b, infile_b, n)

      end if

      if ( read_flux ) then

        write( infile_b, '( i4.4 )' ) d_in%year

        call read_heat_pnetcdf(sh_b,lh_b,sr_b,lr_b,dlr_b,infile_b,n)

      end if

      nb = n

      if ( is_master ) print '(/a/)', "---------- tsforce_init."


      return

      end subroutine tsforce_init
!--------------------------------------------------------------

!==============================================================
! Read & time-interpolation of TS climatorology.
!--------------------------------------------------------------
      subroutine tsforce_main( d_in )

      use bry
      use glob_domain, only: n_east, n_north, n_south, n_west
      use model_run  , only: dti
      use glob_ocean , only: rmean, ssurf, tsurf
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

      days_in_year = real( mday( d_in%month ) - 31 + d_in%day )

!     mid-point [sec] in the month.

      mid_in_month = int( real( mday( d_in%month ) )/ 2.0
     $             * 24. * 3600. )


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
         el_a = el_b
         if ( corr_surface ) then
           dt_a = dt_b
           ds_a = ds_b
         end if
!         uw_a = uw_b
!         vw_a = vw_b

!         tsurf_a = tsurf_b
!         ssurf_a = ssurf_b

!         write( infile_b, '( "tsclimib",i2.2,".nc" )' ) mon_b
         call read_tsclim_monthly_pnetcdf
     $        ( tm_b, sm_b, tc_b, sc_b, el_b, "ts_clim.nc", mon_b )
         if ( corr_surface ) then
           call read_surf_corr_pnetcdf( dt_b, ds_b, "dt.nc", mon_b)
         end if
!         call read_tsclim_monthly_pnetcdf_obs
!     $        ( rm_b, tc_b, sc_b, "ts_clim_old.nc", mon_b )
!         call read_wind_monthly_pnetcdf
!     $        ( uw_b, vw_b,       "sfrc_old.nc",    mon_b )

!         write( infile_b, '( "ssts",i2.2,".nc" )' ) mon_b
!         call read_ssts_monthly_pnetcdf
!     $        ( tsurf_b, ssurf_b, "ts_clim.nc", mon_b )

!         tsurf_b = tc_b(:,:,1)
!         ssurf_b = sc_b(:,:,1)

      endif

!      if ( db /= int( days_in_year ) ) then
!        db = int( days_in_year )
!        tm_a = tm_b
!        sm_a = sm_b
!        call read_mean_ts_z_pnetcdf( tm_b, sm_b, 44, db+1 )
!      end if


!     time interpolation.

      tclim = ( 1. - aa ) * tc_a       + aa * tc_b
      sclim = ( 1. - aa ) * sc_a       + aa * sc_b

      if ( n_north == -1 ) then
        EL_bry % NTH(:,1)
     &      = ( 1. - aa ) * el_a(:,jm) + aa * el_b(:,jm)
      end if

      if ( n_east == -1 ) then
        EL_bry % EST(1,:)
     &      = ( 1. - aa ) * el_a(im,:) + aa * el_b(im,:)
      end if

      if ( n_south == -1 ) then
        EL_bry % STH(:,1)
     &      = ( 1. - aa ) * el_a(:, 1) + aa * el_b(:, 1)
      end if

      if ( n_west == -1 ) then
        EL_bry % WST(1,:)
     &      = ( 1. - aa ) * el_a( 1,:) + aa * el_b( 1,:)
      end if

      if ( corr_surface ) then
        dtemp=( 1. - aa ) * dt_a       + aa * dt_b
        dsalt=( 1. - aa ) * ds_a       + aa * ds_b
      end if
      bb = aa
!      bb = days_in_year - int( days_in_year )
      tmean = ( 1.0 - bb ) * tm_a + bb * tm_b
      smean = ( 1.0 - bb ) * sm_a + bb * sm_b

!      wusurf = ( 1.0 - aa ) * uw_a + aa * uw_b
!      wvsurf = ( 1.0 - aa ) * vw_a + aa * vw_b

      tsurf = t(:,:,1)+dtemp !( 1.0 - aa ) * tsurf_a + aa * tsurf_b
      ssurf = s(:,:,1)+dsalt !( 1.0 - aa ) * ssurf_a + aa * ssurf_b


!     calculation of rmean.

      call dens( smean, tmean, rmean )
!      if (106 > i_global(1) .and. 106 < i_global(im)
!     & .and. 214 > j_global(1) .and. 214 < j_global(jm) ) then
!        print *, rmean(106-i_global(1)+1,214-j_global(jm)+1,:)
!      end if

      return

      end subroutine tsforce_main
!--------------------------------------------------------------

!--------------------------------------------------------------
! Calculate heat/salt fluxes wtsurf and wssurf.
!--------------------------------------------------------------
      subroutine tsforce_tsflx( d_in )

!     surface heat/salt fluxes.

      use air        , only: e_atmos, swrad, wusurf, wvsurf
      use config     , only: calc_bulk, rhoref, sbias, sf_hf, tbias
      use glob_domain, only: my_task
      use glob_grid  , only: east_e, north_e
      use glob_misc  , only: hi, ice
      use module_time

      implicit none

      integer n

      ! intent(in)
      type(date), intent(in) :: d_in
      type(date) d_tmp
      real(kind=rk) bb
      integer sec_in_day
!      character(len=120) datestr

      logical :: lexist


      sec_in_day = d_in%hour*3600 + d_in%min*60 + d_in%sec
      bb = real( mod( sec_in_day, 6 * 3600 ) ) / ( 6. * 3600. )

      d_tmp = str2date("1979-01-01 00:00:00")
      d_tmp%year = d_in%year
      n = int((d_in-d_tmp+6*3600)/86400.*4.)+1
      if (n > 4*(365+inc_leap(d_in%year))) then
        n = mod(n,4*(365+inc_leap(d_in%year)))
        d_tmp%year = d_tmp%year+1
      end if

      if ( n/=nb ) then

        nb = n

        if ( read_bulk ) then

          tair_a = tair_b
          pres_a = pres_b
          shum_a = shum_b
          rain_a = rain_b
          cloud_a= cloud_b
          tskin_a= tskin_b
          uwnd_a = uwnd_b
          vwnd_a = vwnd_b

          write( infile_b, '( i4.4,".nc" )' ) d_tmp%year

          inquire(file='in/heat/'//trim(infile_b),
     $           exist=lexist)

          swrad  = 0.
          call read_bulk_pnetcdf(pres_b,tair_b,shum_b,rain_b
     $                     ,cloud_b,uwnd_b,vwnd_b,tskin_b,infile_b,n)

        end if

        if ( read_flux ) then

          sh_a = sh_b
          lh_a = lh_b
          sr_a = sr_b
          lr_a = lr_b
          dlr_a= dlr_b

          write( infile_b, '( i4.4 )' ) d_in%year

          call read_heat_pnetcdf(sh_b,lh_b,sr_b,lr_b,dlr_b,infile_b,n)

        end if

      end if

      tskin = (1.-bb)*tskin_a + bb*tskin_b
      sskin = (1.-bb)*sc_a(:,:,1) + bb*sc_b(:,:,1)

      if ( read_bulk ) then
        pres  = (1.-bb)*pres_a + bb*pres_b
        rain  = (1.-bb)*rain_a + bb*rain_b
      end if
! Simplified version of inverse barometer:
!   IB(mm) = -9.948 * ( ΔRdry(mbars) – 1013.3 )
! In POM it should be inversed in sign (it seems)
      e_atmos = 0.01 * ( pres - 1013. )

!      datestr = "dbg."//date2str(d_in)
!      call write_sflx(datestr, tair,shum,rain,tskin,pres,cloud )
      if ( read_flux ) then

        sh = (1.-bb)*sh_a + bb*sh_b
        lh = (1.-bb)*lh_a + bb*lh_b
        sr = (1.-bb)*sr_a + bb*sr_b
        lr = (1.-bb)*lr_a + bb*lr_b
        dlr= (1.-bb)*dlr_a+ bb*dlr_b
        swrad = -sr*fsm

      else
        dlr= -370.
      end if

      if ( calc_bulk ) then

        tair  = (1.-bb)*tair_a + bb*tair_b
        shum  = (1.-bb)*shum_a + bb*shum_b
        cloud = (1.-bb)*cloud_a + bb*cloud_b

        uwnd = ( 1.0 - bb ) * uwnd_a + bb * uwnd_b
        vwnd = ( 1.0 - bb ) * vwnd_a + bb * vwnd_b

        call bulk(im,jm,tbias,fsm,t(:,:,1),east_e,north_e,
     $              d_in%year,d_in%month,d_in%day,
     $              d_in%hour,d_in%min,
     $              wusurf,wvsurf,wtsurf,swrad,emp,
     $              uwnd,vwnd,u(:,:,1),v(:,:,1),rho(:,:,1),
     $              tair,shum,rain,cloud,pres,-dlr,my_task)

        wssurf = emp*(sb(:,:,1)+sbias)
! Relax to skin TS... Skin salinity is just climatology
        if ( relax_surface ) then
          wtsurf = wtsurf
     &           + c1 * ( tb(:,:,1) - tskin(:,:) - dtemp )
!     &                / max(-h(:,:)*zz(1), 1.)    !rwnd: (linear) prevention of overheating a thick layer
          wssurf = wssurf
     $           + c1 * ( sb(:,:,1) - sskin(:,:) )
        end if
        wtsurf = (1.-ice)*wtsurf - ice*2.13*(-4.-t(:,:,1))/hi/4.1876d6
!          write(*,*) my_task, "WT:", minval(wtsurf),maxval(wtsurf)
!          write(*,*) my_task, "SW:", minval(swrad),maxval(swrad)
!          write(*,*) my_task, "MP:", minval(emp),maxval(emp)
      else


        do j=1,jm
          do i=1,im

            sstrelx = 1.
            sssrelx = 1.
!lyo:20110202:
!           sssrelx = ( 1.0d0 + tanh( 0.002d0 * (h(i,j)-1000.d0)))*0.5d0
!           sssrelx = ( 1.0d0 + tanh(0.0005d0 * (h(i,j)-2500.d0)))*0.5d0

            swrad( i, j ) = sr( i, j )/3986./rhoref
            if ( read_bulk ) then
              emp( i , j ) = lh(i,j)/heatlat(t(i,j,1))/rhoref
     &                      -rain(i,j)/rho_fw
            end if
            wtsurf( i, j ) = ( lr(i,j)+sh(i,j)+lh(i,j) )
     &                       /3986./rhoref
     &             + emp(i,j)*(tair(i,j)-t(i,j,1))   ! Not needed, since it already is included in latent heat flux, right?
            wtsurf( i, j ) = sf_hf*wtsurf( i, j )
            wssurf( i, j ) = -emp( i, j ) * ( s( i, j, 1 )+sbias )

            if ( relax_surface ) then
              wtsurf(i,j) = wtsurf(i,j)
     &           + c1 * sstrelx * ( tb( i, j, 1 ) - tskin( i, j ) )
              wssurf(i,j) = wssurf(i,j)
     &           + c1 * sssrelx * ( sb( i, j, 1 ) - sskin( i, j ) )
            end if

          end do
        end do

      end if

      return

      end subroutine tsforce_tsflx
!--------------------------------------------------------------


      end module tsforce
