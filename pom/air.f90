!______________________________________________________________________
!
! Module `AIR` (air.f90)
!----------------------------------------------------------------------
!  Module for applying surface forcing.
!
!  Author  : RinceWND
!  Created : 2018-09-06
!______________________________________________________________________

module air

  use glob_const, only: rk

  implicit none

  public

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  real(rk), parameter :: rhoa = 1.22

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED ! Is set according to the external flag `use_air`.

  integer read_int   ! Reading interval in days

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

  logical          & ! Interpolate between records:
    INTERP_CLOUD   & !  total cloud cover
  , INTERP_HEAT    & !  all heat fluxes
  , INTERP_HUMID   & !  humidity
  , INTERP_PRES    & !  atm. pressure
  , INTERP_RAIN    & !  precipitation rate
  , INTERP_SST     & !  sea surface temperature
  , INTERP_STRESS  & !  momentum flux
  , INTERP_TAIR    & !  air temperature
  , INTERP_WIND    & !  wind
                     ! Read surface forcing:
  , READ_BULK      & !  bulk variables (even if not using bulk)
  , READ_HEAT      & !  heat fluxes
  , READ_STRESS    & !  momentum flux
  , READ_WIND      & !  wind (even if using wind stress)
                     ! Treat surface forcing differently:
  , USE_BULK       & !  heat flux is estimated using bulk formulas
  , USE_COARE      & !  use COARE algorithm for bulk formulations
  , USE_STRESS       !  momentum flux is read directly from files

  integer*1        &
    LWRAD_FORMULA    ! Bulk formula to estimate longwave radiation

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=256)  & ! Full paths (with filenames) to:
    bulk_path         & !  atmospheric parameters file
  , flux_path         & !  surface fluxes file
  , wind_path           !  wind file

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(len=32)  &
    dlrad_name       & ! Downward longwave radiation
  , lheat_name       & ! Latent heat flux
  ,  lrad_name       & ! Longwave net radiation
  , humid_name       & ! Relative humidity
  ,  pres_name       & ! Atm. pressure
  ,  rain_name       & ! Precipitation rate
  , sheat_name       & ! Sensible heat flux
  ,  srad_name       & ! Shortwave net radiation
  ,   sst_name       & ! Sea surface temperature
  ,  tair_name       & ! Air temperature
  ,  tcld_name       & ! Total cloud cover
  ,  ustr_name       & ! U-component of momentum flux
  ,  uwnd_name       & ! U-component of wind
  ,  vstr_name       & ! V-component of momentum flux
  ,  vwnd_name       & ! V-component of wind
  ,  wgst_name         ! Wind gustiness [NOT IMPLEMENTED]


!----------------------------------------------------------------------
! Air-Sea interaction related arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    e_atmos                  & ! atmospheric pressure
  , swrad                    & ! short wave radiation incident
                               !  on the ocean surface
  , uwsrf                    & ! wind speed in x-direction
  , vfluxb                   & ! volume flux through water column
                               !  surface at time n-1
  , vfluxf                   & ! volume flux through water column
                               !  surface at time n+1
  , vwsrf                    & ! wind speed in y-direction
  , wssurf                   & ! <ws(0)> salinity flux at the surface
  , wtsurf                   & ! <wt(0)> temperature flux
                               !  at the surface
  , wusurf                   & ! <wu(0)> momentum flux at the surface
  , wvsurf                     ! <wv(0)> momentum flux at the surface

!----------------------------------------------------------------------
! Private variables
!----------------------------------------------------------------------
  real(kind=rk)              &
       , private             &
       , allocatable         &
       , dimension(:,:,:) :: &
    cloud                    & ! Total cloud cover [fraction 0..1]
  , dlrad                    & ! Downward longwave radiation [W/m^2]
  , humid                    & ! Relative humidity [%]
                               !  or specific humidity [kg/??]
  , lheat                    & ! Latent heat flux [W/m^2]
  , lrad                     & ! Net longwave radiation [W/m^2]
  , pres                     & ! Atmospheric pressure [Pa]
  , rain                     & ! Precipetation rate [m/s]
  , sheat                    & ! Sensible heat flux [W/m^2]
  , srad                     & ! Net shortwave radiation [W/m^2]
  , sst                      & ! Sea surface temperature [K] [degC]
  , temp                     & ! Air temperature (@2m) [K] [degC]
  , ustr                     & ! U-component of momentum flux
                               !  [m^2/s^2] [N/m^2]
  , uwnd                     & ! U-component of wind velocity [m/s]
  , vstr                     & ! V-component of momentum flux
                               !  [m^2/s^2] [N/m^2]
  , vwnd                       ! V-component of wind velocity [m/s]
!______________________________________________________________________
!  Radiation and heatflux parameters are negative downward
! (ocean-ward).
!  Net longwave radiation should already include downward
! longwave radiation.


  contains

!______________________________________________________________________
!
    subroutine initialize_air( config_file )
!----------------------------------------------------------------------
!  Initialize air module.
!______________________________________________________________________

      use config     , only: use_air
      use glob_domain, only: is_master

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/air_nml/                                         &
        INTERP_CLOUD , INTERP_HEAT, INTERP_HUMID , INTERP_PRES  &
      , INTERP_RAIN  , INTERP_SST , INTERP_STRESS, INTERP_TAIR  &
      , INTERP_WIND  , READ_BULK  , READ_HEAT    , READ_STRESS  &
      , READ_WIND    , USE_BULK   , USE_COARE    , USE_STRESS   &
      , LWRAD_FORMULA, bulk_path  , flux_path    , wind_path    &
      , read_int

      namelist/air_vars_nml/                           &
        dlrad_name, lheat_name,  lrad_name, rain_name  &
      ,  pres_name, humid_name, sheat_name,  sst_name  &
      ,  srad_name,  tair_name,  tcld_name, ustr_name  &
      ,  uwnd_name,  vstr_name,  vwnd_name, wgst_name

      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_AIR ) DISABLED = .true.

! Allocate mandatory arrays even if not active
      if ( DISABLED ) then
        call allocate_air
        return
      end if

! Initialize variables with their defaults
      read_int = 3600 * 6 ! 86400 / 4 (4xDaily)

      bulk_path = "in/surf/"
      flux_path = "in/surf/"
      wind_path = "in/surf/"

      dlrad_name = "dlwr"
      lheat_name = "lht"
      lrad_name  = "lwr"
      rain_name  = "prate"
      pres_name  = "pres"
      humid_name = "rhum"
      sheat_name = "sht"
      sst_name   = "sst"
      srad_name  = "swr"
      tair_name  = "tair"
      tcld_name  = "cloud"
      ustr_name  = "uflx"
      uwnd_name  = "uwnd"
      vstr_name  = "vflx"
      vwnd_name  = "vwnd"
      wgst_name  = "gust"

      READ_WIND = .false.

! Override configuration
      open ( 73, file = config_file, status = 'old' )
      read ( 73, nml = air_nml )
      read ( 73, nml = air_vars_nml )
      close( 73 )

! Variables management
      pos = len(trim(bulk_path))
      if ( bulk_path(pos:pos) == "/" ) then
        bulk_path = trim(bulk_path)//"blk."
      end if

      pos = len(trim(flux_path))
      if ( flux_path(pos:pos) == "/" ) then
        flux_path = trim(flux_path)//"flx."
      end if

      pos = len(trim(wind_path))
      if ( wind_path(pos:pos) == "/" ) then
        wind_path = trim(wind_path)//"wnd."
      end if

      if ( is_master ) then
        print *, "Bulk data     : ", trim(bulk_path)
        print *, "Surface fluxes: ", trim(flux_path)
        print *, "Wind          : ", trim(wind_path)
      end if

! Allocate necessary arrays
      call allocate_air

      call msg_print("AIR MODULE INITIALIZED", 1, "")

    end subroutine
!______________________________________________________________________
!
    subroutine allocate_air
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________

!      use config     , only: calc_wind
      use glob_domain, only: im, jm

      implicit none

      integer*1 N ! Interpolation array extension size
                  ! The structure is following:
                  !   var at n   step: var(:,:,1)
                  !   var at n-1 step: var(:,:,2)
                  !   var at n+1 step: var(:,:,3)

! Allocate core arrays
      allocate(         &
        e_atmos(im,jm)  &
      , swrad(im,jm)    &
      , uwsrf(im,jm)    &
      , vfluxb(im,jm)   &
      , vfluxf(im,jm)   &
      , vwsrf(im,jm)    &
      , wssurf(im,jm)   &
      , wtsurf(im,jm)   &
      , wusurf(im,jm)   &
      , wvsurf(im,jm)   &
       )

! Initialize mandatory arrays
      e_atmos= 1013.
      swrad  = 0.
      uwsrf  = 0.
      vfluxb = 0.
      vfluxf = 0.
      vwsrf  = 0.
      wssurf = 0.
      wtsurf = 0.
      wusurf = 0.
      wvsurf = 0.

! Quit if the module is not used.
      if ( DISABLED ) return

! Allocate optional arrays
      if ( read_wind ) then
        N = 1
        if ( interp_wind ) N = 3
        allocate(         &
          uwnd(im,jm,N)   &
        , vwnd(im,jm,N)   &
         )
        uwnd = 0.
        vwnd = 0.
      end if

      if ( read_stress ) then
        N = 1
        if ( interp_stress ) N = 3
        allocate(         &
          ustr(im,jm,N)   &
        , vstr(im,jm,N)   &
         )
      end if

      if ( read_bulk ) then
        if ( interp_cloud ) then
          allocate( cloud(im,jm,3) )
        else
          allocate( cloud(im,jm,1) )
        end if
        if ( interp_humid ) then
          allocate( humid(im,jm,3) )
        else
          allocate( humid(im,jm,1) )
        end if
        if ( interp_pres ) then
          allocate( pres(im,jm,3) )
        else
          allocate( pres(im,jm,1) )
        end if
        if ( interp_rain ) then
          allocate( rain(im,jm,3) )
        else
          allocate( rain(im,jm,1) )
        end if
        if ( interp_tair ) then
          allocate( temp(im,jm,3) )
        else
          allocate( temp(im,jm,1) )
        end if
      end if

      if ( read_heat ) then
        N = 1
        if ( interp_heat ) N = 3
        allocate(         &
          sheat(im,jm,N)  &
        , srad(im,jm,N)   &
        , lheat(im,jm,N)  &
        , lrad(im,jm,N)   &
         )
      end if

    end subroutine ! allocate_air
!______________________________________________________________________
!
    subroutine air_init( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields before experiment's start.
!______________________________________________________________________

!      use glob_domain, only: is_master
!      use config     , only: calc_wind
      use module_time
      use glob_ocean , only: tb

      implicit none

      type(date), intent(in) :: d_in

      integer            max_in_prev, max_in_this, ncid
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

! Set ncid to -1 to open first file. Then set it to 0 to close previous file and open new one. TODO: Stupid, I know.
      ncid = -1

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      chunk     = chunk_of_year( d_in, read_int )
      record(1) = int(chunk)

      if ( chunk - record(1) < .5 ) then
        record(2) = record(1)
      else
        record(2) = record(1) + 1
      end if
      record(3) = record(2) + 1

      if ( record(2) == 0 ) then
        record(2) = max_in_prev + 1 ! TODO: [ NEEDS TESTING ]
        year(2) = d_in%year - 1
      elseif ( record(3) == max_in_this + 1 ) then
        record(3) = 1
        year(3) = d_in%year + 1
      end if

! Read wind if momentum flux is derived from wind
      if ( read_wind ) then
        if ( interp_wind ) then
          call read_var_nc( wind_path, uwnd_name, uwnd(:,:,2)  &
                          , year(2), record(2), ncid )
          call read_var_nc( wind_path, vwnd_name, vwnd(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            uwnd(:,:,2) = 0.
            vwnd(:,:,2) = 0.
          end if
        else
          call read_var_nc( wind_path, uwnd_name, uwnd(:,:,1)  &
                          , year(1), record(1), ncid )
          call read_var_nc( wind_path, vwnd_name, vwnd(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            uwnd(:,:,1) = 0.
            vwnd(:,:,1) = 0.
          end if
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Wind is set to zero.")
        end if
      end if

      ncid = 0

      if ( read_bulk ) then
        if ( interp_humid ) then
          call read_var_nc( bulk_path,humid_name, humid(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) humid(:,:,2) = 1.
        else
          call read_var_nc( bulk_path,humid_name, humid(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) humid(:,:,1) = 1.
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Humidity is set to 100%")
        end if
        if ( interp_rain ) then
          call read_var_nc( bulk_path, rain_name,  rain(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) rain(:,:,2) = 0.
        else
          call read_var_nc( bulk_path, rain_name,  rain(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) rain(:,:,1) = 0.
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Precip.rate is set to zero.")
        end if
        if ( interp_pres ) then
          call read_var_nc( bulk_path, pres_name,  pres(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) pres(:,:,2) = 1013.
        else
          call read_var_nc( bulk_path, pres_name,  pres(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) pres(:,:,1) = 1013.
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Atm.pressure is set to 1013 hPa.")
        end if
        if ( interp_sst ) then
          call read_var_nc( bulk_path,  sst_name,   sst(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) sst(:,:,2) = tb(:,:,1)
        else
          call read_var_nc( bulk_path,  sst_name,   sst(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) sst(:,:,1) = tb(:,:,1)
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "SST is set to surface level temp.")
        end if
        if ( interp_tair ) then
          call read_var_nc( bulk_path, tair_name,  temp(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) temp(:,:,2) = tb(:,:,1)
        else
          call read_var_nc( bulk_path, tair_name,  temp(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) temp(:,:,1) = tb(:,:,1)
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Tair is set to surface level temp.")
        end if
        if ( interp_cloud ) then
          call read_var_nc( bulk_path, tcld_name, cloud(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) cloud(:,:,2) = 0.
        else
          call read_var_nc( bulk_path, tcld_name, cloud(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) cloud(:,:,1) = 0.
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Cloud cover is set to zero.")
        end if
      end if

      ncid = 0

      if ( read_stress ) then
        if ( interp_stress ) then
          call read_var_nc( flux_path, ustr_name, ustr(:,:,2)  &
                          , year(2), record(2), ncid )
          call read_var_nc( flux_path, vstr_name, vstr(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            ustr(:,:,2) = 0.
            vstr(:,:,2) = 0.
          end if
        else
          call read_var_nc( flux_path, ustr_name, ustr(:,:,1)  &
                          , year(1), record(1), ncid )
          call read_var_nc( flux_path, vstr_name, vstr(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            ustr(:,:,1) = 0.
            vstr(:,:,1) = 0.
          end if
        end if
        if ( ncid == -1 ) then
          call msg_print("", 2, "Wind stress is set to zero.")
        end if
      end if

      ncid = 0

      if ( read_heat ) then
        if ( interp_heat ) then
          call read_var_nc( flux_path, dlrad_name, dlrad(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            dlrad(:,:,2) = 370.
            call msg_print("", 2, "Downlongrad. is set to 370 W/m^2.")
          end if
          call read_var_nc( flux_path, lheat_name, lheat(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            lheat(:,:,2) = 0.
            call msg_print("", 2, "Latent heat is set to zero.")
          end if
          call read_var_nc( flux_path,  lrad_name,  lrad(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            lrad(:,:,2) = 0.
            call msg_print("", 2, "Net longrad. is set to zero.")
          end if
          call read_var_nc( flux_path, sheat_name, sheat(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            sheat(:,:,2) = 0.
            call msg_print("", 2, "Sensible heat is set to zero.")
          end if
          call read_var_nc( flux_path,  srad_name,  srad(:,:,2)  &
                          , year(2), record(2), ncid )
          if ( ncid == -1 ) then
            srad(:,:,2) = 0.
            call msg_print("", 2, "Net shortrad. is set to zero.")
          end if
        else
          call read_var_nc( flux_path, dlrad_name, dlrad(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            dlrad(:,:,1) = 0.
            call msg_print("", 2, "Downlongrad. is set to 370 W/m^2.")
          end if
          call read_var_nc( flux_path, lheat_name, lheat(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            lheat(:,:,1) = 0.
            call msg_print("", 2, "Latent heat is set to zero.")
          end if
          call read_var_nc( flux_path,  lrad_name,  lrad(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            lrad(:,:,1) = 0.
            call msg_print("", 2, "Net longrad. is set to zero.")
          end if
          call read_var_nc( flux_path, sheat_name, sheat(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            sheat(:,:,1) = 0.
            call msg_print("", 2, "Sensible heat is set to zero.")
          end if
          call read_var_nc( flux_path,  srad_name,  srad(:,:,1)  &
                          , year(1), record(1), ncid )
          if ( ncid == -1 ) then
            srad(:,:,1) = 0.
            call msg_print("", 2, "Net shortrad. is set to zero.")
          end if
        end if
      end if


      call msg_print("AIR INITIALIZED", 2, "")


      return

    end subroutine air_init
!______________________________________________________________________
!
    subroutine air_step( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields during experiment.
!______________________________________________________________________

!      use glob_const , only: SEC2DAY
      use config     , only: rhow => rhoref
      use glob_domain, only: im, is_master, jm
      use module_time
      use glob_ocean , only: u, v
      use model_run  , only: dti, iint

      implicit none

      type(date), intent(in) :: d_in

      logical            ADVANCE_REC, ADVANCE_REC_INT
      integer            i, j, max_in_prev, max_in_this, ncid, secs
      integer                          &
      , dimension(3)  :: record, year
      real               a, cda, chunk, uvabs


! Quit if the module is not used.
      if ( DISABLED ) return

! Set ncid to -1 to open first file. Then set it to 0 to close previous file and open new one. TODO: Stupid, I know.
      ncid = -1

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      chunk     = chunk_of_year( d_in, read_int )
      record(1) = int(chunk)
      secs      = seconds_of_year(d_in) - (real(record(1))+.5)*read_int

      if ( chunk - record(1) < .5 ) then ! TODO: test (real - int) = ?
        record(2) = record(1)
        a = chunk - record(1) + .5
      else
        record(2) = record(1) + 1
        a = chunk - record(1) - .5
      end if
      record(3) = record(2) + 1

      if ( record(2) == 0 ) then
        record(2) = max_in_prev + 1 ! TODO: [ NEEDS TESTING ]
        year(2) = d_in%year - 1
      elseif ( record(3) == max_in_this + 1 ) then
        record(3) = 1
        year(3) = d_in%year + 1
      end if

      if ( iint == 1 .or. ( secs >= 0 .and. secs < dti ) ) then
        ADVANCE_REC_INT = .true.
        if ( interp_wind ) then
          uwnd(:,:,2) = uwnd(:,:,3)
          vwnd(:,:,2) = vwnd(:,:,3)
        end if
      else
        ADVANCE_REC_INT = .false.
      end if

      if ( (chunk-record(1))*read_int <= int(dti) ) then ! TODO: test this one as well.
        ADVANCE_REC = .true.
      else
        ADVANCE_REC = .false.
      end if

!      if ( is_master ) print *, "II", (chunk-record(1)) &
!                                    , (chunk-record(1)-.5) &
!                              , secs, a  &
!                              , ADVANCE_REC_INT, ADVANCE_REC

      call read_all( ADVANCE_REC_INT, 3, year, record )
      call read_all( ADVANCE_REC    , 1, year, record )

      if ( interp_wind ) then
        uwnd(:,:,1) = ( 1. - a ) * uwnd(:,:,2) + a * uwnd(:,:,3)
        vwnd(:,:,1) = ( 1. - a ) * vwnd(:,:,2) + a * vwnd(:,:,3)
      end if

! Calculate wind stress (m^2/s^2)
      do j=1,jm
        do i=1,im

          uwsrf(i,j) = uwnd(i,j,1)
          vwsrf(i,j) = vwnd(i,j,1)

          uvabs = sqrt( (uwnd(i,j,1)-u(i,j,1))**2  &
                      + (vwnd(i,j,1)-v(i,j,1))**2 )

          if (     uvabs <= 11. ) then
            cda = .0012
          elseif ( uvabs <= 19. ) then
            cda = .00049  + .000065 *uvabs
          elseif ( uvabs <= 100. ) then
            cda = .001364 + .0000234*uvabs - 2.31579e-7*uvabs**2
          else
            cda = .00138821
          end if

          wusurf(i,j) = -rhoa/rhow*cda*uvabs * (uwnd(i,j,1)-u(i,j,1))
          wvsurf(i,j) = -rhoa/rhow*cda*uvabs * (vwnd(i,j,1)-v(i,j,1))

        end do
      end do


    end subroutine air_step
!______________________________________________________________________
!
    subroutine calculate_fluxes
!----------------------------------------------------------------------
!  Generates and interpolates all surface fluxes.
!______________________________________________________________________

      implicit none


    end subroutine
!______________________________________________________________________
!
    subroutine read_all( execute, n, year, record )

      use glob_ocean, only: tb

      implicit none

      logical              , intent(in) :: execute
      integer              , intent(in) :: n
      integer, dimension(3), intent(in) :: record, year

      integer            ncid
      character(len=128) desc


      if ( .not. execute ) return

      if ( n == 3 ) then
        write(desc,'("Reading interp. forcing record #",i4," @ ",i4)') &
            record(3), year(3)
      else
        write(desc,'("Reading surface forcing record #",i4," @ ",i4)') &
            record(1), year(1)
      end if

      call msg_print("", 1, desc)

      ncid = -1

      if ( read_wind ) then

        if ((     interp_wind .and. n==3) .or.       &
            (.not.interp_wind .and. n==1)      ) then
          call read_var_nc( wind_path, uwnd_name, uwnd(:,:,n)  &
                          , year(n), record(n), ncid )
          call read_var_nc( wind_path, vwnd_name, vwnd(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            uwnd(:,:,n) = 0.
            vwnd(:,:,n) = 0.
            call msg_print("", 2, "Wind is set to zero.")
          end if
        end if

      end if

      ncid = 0

      if ( read_bulk ) then

        if ((     interp_humid .and. n==3) .or.       &
            (.not.interp_humid .and. n==1)      ) then
          call read_var_nc( bulk_path,humid_name, humid(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            humid(:,:,n) = 1.
            call msg_print("", 2, "Humidity is set to 100%")
          end if
        end if

        if ((     interp_rain .and. n==3) .or.       &
            (.not.interp_rain .and. n==1)      ) then
          call read_var_nc( bulk_path, rain_name,  rain(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            rain(:,:,n) = 0.
            call msg_print("", 2, "Precip.rate is set to zero.")
          end if
        end if

        if ((     interp_pres .and. n==3) .or.       &
            (.not.interp_pres .and. n==1)      ) then
          call read_var_nc( bulk_path, pres_name,  pres(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            pres(:,:,n) = 1013.
            call msg_print("", 2, "Atm.pressure is set to 1013 hPa.")
          end if
        end if

        if ((     interp_sst .and. n==3) .or.       &
            (.not.interp_sst .and. n==1)      ) then
          call read_var_nc( bulk_path,  sst_name,   sst(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            sst(:,:,n) = tb(:,:,1)
            call msg_print("", 2, "SST is set to surf. level temp.")
          end if
        end if

        if ((     interp_tair .and. n==3) .or.       &
            (.not.interp_tair .and. n==1)      ) then
          call read_var_nc( bulk_path, tair_name,  temp(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            temp(:,:,n) = tb(:,:,1)
            call msg_print("", 2, "Tair is set to surf. level temp.")
          end if
        end if

        if ((     interp_cloud .and. n==3) .or.       &
            (.not.interp_cloud .and. n==1)      ) then
          call read_var_nc( bulk_path, tcld_name, cloud(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            cloud(:,:,n) = 0.
            call msg_print("", 2, "Cloud cover is set to zero.")
          end if
        end if

      end if ! if READ_BULK

      ncid = 0

      if ( read_stress ) then

        if ((     interp_stress .and. n==3) .or.       &
            (.not.interp_stress .and. n==1)      ) then
          call read_var_nc( flux_path, ustr_name, ustr(:,:,n)  &
                          , year(n), record(n), ncid )
          call read_var_nc( flux_path, vstr_name, vstr(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            ustr(:,:,n) = 0.
            vstr(:,:,n) = 0.
            call msg_print("", 2, "Wind stress is set to zero.")
          end if
        end if

      end if

      ncid = 0

      if ( read_heat ) then

        if ((     interp_heat .and. n==3) .or.       &
            (.not.interp_heat .and. n==1)      ) then
          call read_var_nc( flux_path, dlrad_name, dlrad(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            dlrad(:,:,n) = 370.
            call msg_print("", 2, "Downlonrad. is set to 370 W/m^2.")
          end if
          call read_var_nc( flux_path, lheat_name, lheat(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            lheat(:,:,n) = 0.
            call msg_print("", 2, "Latent heat is set to zero.")
          end if
          call read_var_nc( flux_path,  lrad_name,  lrad(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            lrad(:,:,n) = 0.
            call msg_print("", 2, "Net longrad. is set to zero.")
          end if
          call read_var_nc( flux_path, sheat_name, sheat(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            sheat(:,:,n) = 0.
            call msg_print("", 2, "Sensible heat is set to zero.")
          end if
          call read_var_nc( flux_path,  srad_name,  srad(:,:,n)  &
                          , year(n), record(n), ncid )
          if ( ncid == -1 ) then
            srad(:,:,n) = 0.
            call msg_print("", 2, "Net shortrad. is set to zero.")
          end if
        end if
      end if ! if READ_HEAT


    end subroutine
!______________________________________________________________________
!
    pure character(len=256) function get_filename( path, year )
!----------------------------------------------------------------------
!  Costructs filename string in `<path>YYYY<FORMAT_EXT>` format.
!______________________________________________________________________

      implicit none

      character(len=*), intent(in) :: path
      integer         , intent(in) :: year


      write( get_filename, '( a, i4.4, a )' ) trim(path)      &
                                            , year            &
                                            , trim(FORMAT_EXT)

    end function

!= I/O SECTION ========================================================
!______________________________________________________________________
!
    subroutine read_var_nc( path, var_name, var   &
                          , year, record  , ncid )
!----------------------------------------------------------------------
!  Read a variable (NC format).
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_3d

      integer                   , intent(inout) :: ncid
      integer                   , intent(in   ) :: record, year
      real(rk), dimension(im,jm), intent(  out) :: var
      character(len=*)          , intent(in   ) :: path, var_name

      integer                  varid, status
      integer(MPI_OFFSET_KIND) start(4), edge(4)
      character(len=256)       filename, netcdf_file


      start = 1
      edge  = 1

      if ( ncid <= 0 ) then
! open netcdf file
        if ( ncid == 0 ) then
          call check( nf90mpi_close( ncid )              &
                    , 'nfmpi_close:'//trim(netcdf_file) )
        end if

        filename = get_filename( path, year )
        netcdf_file = trim(filename)
        status = nf90mpi_open( POM_COMM, netcdf_file, NF_NOWRITE   &
                             , MPI_INFO_NULL, ncid )
        if ( status /= NF_NOERR ) then
          call msg_print("", 2, "Failed reading `"//trim(var_name)  &
                              //"` from `"//trim(filename)//"`")
          ncid = -1
          return
        end if
      end if

! get variable
      call check( nf90mpi_inq_varid( ncid, var_name, varid )  &
                , 'nfmpi_inq_varid: '//trim(var_name) )

! set reading area
      start(1) = i_global(1)
      start(2) = j_global(1)
      start(3) = record
      edge(1) = im
      edge(2) = jm
      edge(3) =  1

! get data
      call check( get_var_real_3d                    &
                  ( ncid, varid, start, edge, var )  &
                , 'get_var_real: '//trim(var_name) )


      return

      end

!______________________________________________________________________
!
    subroutine check(status, routine)
!----------------------------------------------------------------------
!  Checks for NetCDF I/O error and exits with an error message if hits.
!______________________________________________________________________

      use glob_domain, only: error_status, is_master
      use pnetcdf    , only: nf90mpi_strerror, NF_NOERR

      implicit none

      integer         , intent(in) :: status
      character(len=*), intent(in) :: routine


      if ( status /= NF_NOERR ) then
        error_status = 1
        if ( is_master ) then
          print '(/a,a)', 'IO error at module `AIR`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
          stop
        end if
      end if


      return

    end



end module air
