!______________________________________________________________________
!
! Module `AIR` (air.f90)
!----------------------------------------------------------------------
!  Module for applying surface forcing.
!
!  Author  : RinceWND
!  Created : 2018-09-06
!______________________________________________________________________
!
module air

  use glob_const, only: DEG2RAD, PI, rk

  implicit none

  public

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  real(rk), parameter ::     & !
    CH     =     .66e-3      &  ! Sensible heat transfer coefficient (for stable case)
  , Cp     = 1004.8          &  ! Specific heat capacity of air
  , emiss  =     .97         &  ! Ocean emissivity
  , expsi  =     .622        &  ! ???
  , Ps     = 1013.           &  ! Surface air pressure
  , Rd     =  287.           &  ! Dry air gas constant
  , RHOa   =    1.22         &  ! Density of air
  , sigma  =    5.670367e-8  &  ! Stefan-Boltzmann constant (kg s^-3 K^-4)
  , solar  = 1350.              ! Solar constant (W/m^2)

  integer(1)                &
      , parameter           &
      , private    ::       &
    lwBERLIAND = 3          &
  , lwBIGNAMI  = 0          &
  , lwHERZFELD = 2          &
  , lwMAY      = 1
!

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED ! Is set according to the external flag `use_air`.

  integer, private :: read_int  ! interval for reading (days)
  real   , private :: a         ! time-interpolation factor

  character(10)                           &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

  logical          &
    CALC_SWR       & ! Calculate shortwave radiation internally (Reed)
                     ! Interpolate between records:
  , INTERP_CLOUD   & !  total cloud cover
  , INTERP_HEAT    & !  all heat fluxes
  , INTERP_HUMID   & !  humidity
  , INTERP_PRES    & !  atm. pressure
  , INTERP_RAIN    & !  precipitation rate
  , INTERP_SSS     & !  sea surface salinity
  , INTERP_SST     & !  sea surface temperature
  , INTERP_STRESS  & !  momentum flux
  , INTERP_TAIR    & !  air temperature
  , INTERP_WIND    & !  wind
                     ! Read surface forcing:
  , READ_BULK      & !  bulk variables (even if not using bulk)
  , READ_CLOUD     & !  total cloud cover to use in long- and shortwave radiation calculations
  , READ_HEAT      & !  heat fluxes
  , READ_STRESS    & !  momentum flux
  , READ_WIND      & !  wind (even if using wind stress)
  , TAPER_BRY      & ! Apply tapering along boundaries
                     ! Treat surface forcing differently:
  , USE_BULK       & !  heat flux is estimated using bulk formulas
  , USE_CALENDAR   & !  lookup for forcing file's record using model date (otherwise reading starts from the first record)
  , USE_COARE      & !  use COARE algorithm for bulk formulations
  , USE_DQDSST     & !  use heat flux sensitivity to restore temperature
  , USE_FLUXES     & !  momentum flux is read directly from files !TODO: if set to false, wtsurf crashes the calculations
  , USE_RAMP         !  use temporal linear damping on wind stress

  integer(1)        &
    LWRAD_FORMULA    ! Bulk formula to estimate longwave radiation

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(256)      & ! Full paths (with filenames) to:
    bulk_path         & !  atmospheric parameters file
  , flux_path         & !  surface fluxes file
  , wind_path           !  wind file

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(32)      &
    dlrad_name       & ! Downward longwave radiation
  , lheat_name       & ! Latent heat flux
  ,  lrad_name       & ! Longwave net radiation
  , humid_name       & ! Relative humidity
  ,  pres_name       & ! Atm. pressure
  ,  rain_name       & ! Precipitation rate
  , sheat_name       & ! Sensible heat flux
  ,  srad_name       & ! Shortwave net radiation
  ,   sss_name       & ! Sea surface salinity
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
  real(rk)                   &
       , allocatable         &
       , dimension(:,:)   :: &
    e_atmos                  & ! atmospheric pressure
  , swrad                    & ! short wave radiation incident
                               !  on the ocean surface
  , taper_mask               & ! mask for surface forcing tapering
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
  real(rk)                   &
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
  , rain                     & ! Precipitation rate [m/s]
  , sheat                    & ! Sensible heat flux [W/m^2]
  , srad                     & ! Net shortwave radiation [W/m^2]
  , sss                      & ! Sea surface salinity [psu]
  , sst                      & ! Sea surface temperature [K] [degC]
  , temp                     & ! Air temperature (@2m) [K] [degC]
  , ustr                     & ! U-component of momentum flux
                               !  [m^2/s^2] [N/m^2]
  , uwnd                     & ! U-component of wind velocity [m/s]
  , vstr                     & ! V-component of momentum flux
                               !  [m^2/s^2] [N/m^2]
  , vwnd                       ! V-component of wind velocity [m/s]
!______________________________________________________________________
!  Radiation and heat flux parameters are negative downward
! (ocean-ward).
!  Net longwave radiation should already include downward
! longwave radiation.


  contains

!______________________________________________________________________
!
    subroutine initialize_mod( config_file )
!----------------------------------------------------------------------
!  Initialize air module.
!----------------------------------------------------------------------
! called by: initialize_arrays [initialize.f]
!
! calls    : allocate_arrays   [air]
!______________________________________________________________________
!
      use config     , only: use_air
      use glob_domain, only: im,jm, is_master, n_east,n_north,n_south,n_west
      use glob_grid  , only: fsm

      implicit none

      character(*), intent(in) :: config_file

      integer pos

      namelist/air_nml/                                             &
        CALC_SWR     , INTERP_CLOUD , INTERP_HEAT  , INTERP_HUMID   &
      , INTERP_PRES  , INTERP_RAIN  , INTERP_SSS   , INTERP_SST     &
      , INTERP_STRESS, INTERP_TAIR  , INTERP_WIND  , READ_BULK      &
      , READ_CLOUD   , READ_HEAT    , READ_STRESS  , READ_WIND      &
      , TAPER_BRY    , USE_BULK     , USE_CALENDAR , USE_COARE      &
      , USE_DQDSST   , USE_FLUXES   , USE_RAMP     , LWRAD_FORMULA  &
      , bulk_path    , flux_path    , wind_path    , read_int

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
        call allocate_arrays
        return
      end if

! Initialize variables with their defaults
      read_int = 3600 * 3 ! 8xDaily

      CALC_SWR     = .true.
      READ_HEAT    = .false.
      READ_STRESS  = .false.
      READ_WIND    = .false.
      TAPER_BRY    = .false.
      READ_CLOUD   = .true.
      USE_BULK     = .false.
      USE_COARE    = .false.
      USE_FLUXES   = .false.
      USE_DQDSST   = .true.
      USE_CALENDAR = .true.
      USE_RAMP     = .false.

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
      sss_name   = ""
      sst_name   = "sst"
      srad_name  = "swr"
      tair_name  = "tair"
      tcld_name  = "cloud"
      ustr_name  = "uflx"
      uwnd_name  = "uwnd"
      vstr_name  = "vflx"
      vwnd_name  = "vwnd"
      wgst_name  = "gust"

      INTERP_HEAT   = .false.
      INTERP_PRES   = .false.
      INTERP_STRESS = .false.

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
        print *, "------------"
        print *, " calculate solar radiation: ", CALC_SWR
        print *, " read heat fluxes         : ", READ_HEAT
        print *, " read wind stress         : ", READ_STRESS
        print *, " read wind velocity       : ", READ_WIND
        print *, " taper boundary           : ", TAPER_BRY
        print *, " read cloud cover         : ", READ_CLOUD
        print *, " calculate bulk heat flux : ", USE_BULK
        print *, " use TOGA COARE for"
        print *, "  sensible and latent flux: ", USE_COARE
        print *, " heat fluxes enabled      : ", USE_FLUXES
        print *, " apply dQdSST correction  : ", USE_DQDSST
        print *, " use date of year to read : ", USE_CALENDAR
        print *, " ramp forcing             : ", USE_RAMP
      end if

! Allocate necessary arrays
      call allocate_arrays

      if ( n_east==-1 ) then
        where ( fsm(im,:) /= 0. )
          taper_mask(im  ,:) = 0.
          taper_mask(im-1,:) = 0.2
          taper_mask(im-2,:) = 0.5
          taper_mask(im-3,:) = 0.8
        end where
      end if

      if ( n_west==-1 ) then
        where ( fsm(1,:) /= 0. )
          taper_mask(1,:) = 0.
          taper_mask(2,:) = 0.
          taper_mask(3,:) = 0.2
          taper_mask(4,:) = 0.5
          taper_mask(5,:) = 0.8
        end where
      end if

      if ( n_north==-1 ) then
        where ( fsm(:,jm) /= 0. )
          taper_mask(:,jm  ) = 0.
          taper_mask(:,jm-1) = 0.2*taper_mask(:,jm-1)
          taper_mask(:,jm-2) = 0.5*taper_mask(:,jm-2)
          taper_mask(:,jm-3) = 0.8*taper_mask(:,jm-3)
        end where
      end if

      if ( n_south==-1 ) then
        where ( fsm(:,1) /= 0. )
          taper_mask(:,1) = 0.
          taper_mask(:,2) = 0.
          taper_mask(:,3) = 0.2*taper_mask(:,3)
          taper_mask(:,4) = 0.5*taper_mask(:,4)
          taper_mask(:,5) = 0.8*taper_mask(:,5)
        end where
      end if

      call msg_print("AIR MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_mod
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!----------------------------------------------------------------------
! called by: initialize_mod    [air]
!______________________________________________________________________
!
      use glob_domain, only: im, jm

      implicit none

      integer(1) N ! Interpolation array extension size
                   ! The structure is following:
                   !   var at n   step: var(:,:,1)
                   !   var at n-1 step: var(:,:,2)
                   !   var at n+1 step: var(:,:,3)


! Allocate core arrays
      allocate(            &
        e_atmos(im,jm)     &
      , swrad(im,jm)       &
      , taper_mask(im,jm)  &
      , uwsrf(im,jm)       &
      , vfluxb(im,jm)      &
      , vfluxf(im,jm)      &
      , vwsrf(im,jm)       &
      , wssurf(im,jm)      &
      , wtsurf(im,jm)      &
      , wusurf(im,jm)      &
      , wvsurf(im,jm)      &
       )

! Initialize mandatory arrays
      e_atmos    = 0.
      swrad      = 0.
      taper_mask = 1.
      uwsrf      = 0.
      vfluxb     = 0.
      vfluxf     = 0.
      vwsrf      = 0.
      wssurf     = 0.
      wtsurf     = 0.
      wusurf     = 0.
      wvsurf     = 0.

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
        pres = Ps
        if ( interp_sss ) then
          allocate( sss(im,jm,3) )
        else
          allocate( sss(im,jm,1) )
        end if
        if ( interp_sst ) then
          allocate( sst(im,jm,3) )
        else
          allocate( sst(im,jm,1) )
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

      N = 1
      if ( read_heat .or. use_fluxes ) then
        if ( interp_heat ) N = 3
        allocate(         &
          sheat(im,jm,N)  &
        , srad(im,jm,N)   &
        , lheat(im,jm,N)  &
        , lrad(im,jm,N)   &
         )
      end if
      allocate( dlrad(im,jm,N) )
      dlrad = 370.


    end ! subroutine allocate_arrays
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields before experiment's start.
!----------------------------------------------------------------------
! called by: modules_initial_step  [initialize.f]
!
! calls    : chunk_of_year         [module_time.f]
!            max_chunks_in_year    [module_time.f]
!            read_all              [air]
!______________________________________________________________________
!
!      use glob_domain, only: is_master
      use module_time
      use config     , only: nbct
      use seaice     , only: icec!, itsurf

      implicit none

      type(date), intent(in) :: d_in

      integer                       max_in_prev, max_in_this, ncid
      integer, dimension(3)      :: record, year
      real(rk)                      a, chunk


! Quit if the module is not used.
      if ( DISABLED ) return

! Set ncid to -1 to open first file. Then set it to 0 to close previous file and open new one. TODO: Stupid, I know.
      ncid = -1

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      if ( USE_CALENDAR ) then

        chunk     = chunk_of_year( d_in, read_int )
        record(1) = int(chunk)

        if ( chunk - record(1) < .5 ) then
          record(2) = record(1)
          a = chunk - record(1) + .5
        else
          record(2) = record(1) + 1
          a = chunk - record(1) - .5
        end if
        record(3) = record(2) + 1

        if ( record(2) == 0 ) then
          record(2) = max_in_prev ! TODO: [ NEEDS TESTING ]
          year(2) = d_in%year - 1
        elseif ( record(3) == max_in_this ) then
          record(3) = 1
          year(3) = d_in%year + 1
        end if

      else

        record = [ 1, 1, 2 ]
        a = 0._rk

      end if

! HARDCODED!!! TODO: Remove
!      record(1) = d_in%month
!      record(2) = d_in%month
!      record(3) = d_in%month+1

! Read wind if momentum flux is derived from wind
      call read_all( .true., 1, year, record )
      call read_all( .true., 2, year, record )
      call read_all( .true., 3, year, record )

      if ( READ_WIND ) then

        if ( interp_wind ) then
          call linint_vec( uwnd(:,:,2), vwnd(:,:,2)  &
                         , uwnd(:,:,3), vwnd(:,:,3)  &
                         , a                         &
                         , uwnd(:,:,1), vwnd(:,:,1)  &
                         )
        end if

! TODO: Are these surface arrays even needed?
        uwsrf = uwnd(:,:,1)
        vwsrf = vwnd(:,:,1)

! Calculate wind stress
        if ( USE_BULK ) then
          if ( READ_HEAT ) then
            dlrad(:,:,1) = ( 1. - a ) * dlrad(:,:,2) + a * dlrad(:,:,3)
          end if
          call wind_to_stress( uwsrf, vwsrf, wusurf, wvsurf, 1 )
          call calculate_fluxes
        end if

      end if

      if ( USE_FLUXES ) then
!        if ( INTERP_STRESS ) then
!          ustr(:,:,1) = ( 1. - a ) * ustr(:,:,2) + a * ustr(:,:,3)
!          vstr(:,:,1) = ( 1. - a ) * vstr(:,:,2) + a * vstr(:,:,3)
          call wind_to_stress( uwsrf, vwsrf, wusurf, wvsurf, 1 )
!        else
!          wusurf = ustr(:,:,1)
!          wvsurf = vstr(:,:,1)
!        end if
        if ( INTERP_HEAT ) then
          wtsurf = ( 1. - a ) * ( lrad(:,:,2)+sheat(:,:,2)+lheat(:,:,2) )  &
                 +        a   * ( lrad(:,:,3)+sheat(:,:,3)+lheat(:,:,3) )
          swrad  = ( 1. - a ) * srad(:,:,2) + a * srad(:,:,3)
        else
          wtsurf = lrad(:,:,1) + sheat(:,:,1) + lheat(:,:,1)
          swrad  = srad(:,:,1)
        end if
      end if

! If solar radiation does not penetrate water
      if ( nbct == 1 ) wtsurf = wtsurf + swrad

! Simple parameterisation
!  Assume that at least 10% solar penetration exists even at 10/10 sea ice concentration.
! Scale linearly the rest.
!  Assume that there is no thermal radiation from sea at 10/10 concentrations.
! To adjust subice sea temperature surface relaxation is being used at the moment.
!  Assume that there is no salt exchange at 10/10 concentrations.
      swrad  =  swrad*(1.-icec*0.9)
      wtsurf = wtsurf*(1.-icec)! + itsurf*icec
      wssurf = wssurf*(1.-icec)

      if (.false.) call river_flux

      if ( TAPER_BRY ) call taper_forcing

      call msg_print("AIR INITIALIZED", 2, "")


    end ! subroutine init
!
!______________________________________________________________________
!
    subroutine step( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields during experiment.
!----------------------------------------------------------------------
! called by: update_bc             [advance.f]
!
! calls    : calculate_fluxes      [air]
!            chunk_of_year         [module_time.f]
!            max_chunks_in_year    [module_time.f]
!            read_all              [air]
!            wind_to_stress        [air]
!______________________________________________________________________
!
      use clim       , only: relax_surface
      use config     , only: nbct
      use module_time
      use model_run  , only: dti, iint, sec_of_year
      use seaice     , only: icec!, itsurf

      implicit none

      type(date), intent(in) :: d_in

      logical            ADVANCE_REC, ADVANCE_REC_INT
      integer            max_in_prev, max_in_this, secs
      integer                          &
      , dimension(3)  :: record, year
      real(rk)           a, chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      secs = sec_of_year

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

      ADVANCE_REC     = .false.
      ADVANCE_REC_INT = .false.

! Decide on the record to read
      if ( USE_CALENDAR ) then

        chunk     = chunk_of_year( d_in, read_int )
        record(1) = int(chunk)

        if ( secs - int(real(record(1))*read_int) < dti ) then ! TODO: test this one as well.
          ADVANCE_REC = .true.
        else
          ADVANCE_REC = .false.
        end if

        secs = secs - int((real(record(1))+.5)*read_int)

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

        if ( secs >= 0 .and. secs < dti ) ADVANCE_REC_INT = .true.

      else

        record(1) = int( iint*dti / read_int ) + 1
        record(2) = record(1)
        record(3) = record(2) + 1

! TODO: right now it interpolates between records (edges of rec1-rec2 span). Implement interpolation between the centers of record spans.
        a = modulo( real(iint*dti), real(read_int) )
        if ( a < dti ) then
          if ( a >= 0. ) then
            ADVANCE_REC = .true.
          end if
        end if
        a = a / read_int
!        print *, "A: ", a

        ADVANCE_REC_INT = ADVANCE_REC

      end if

      if ( ADVANCE_REC_INT ) then
        if ( interp_wind ) then
          uwnd(:,:,2) = uwnd(:,:,3)
          vwnd(:,:,2) = vwnd(:,:,3)
        end if
        if ( .not.CALC_SWR .and. INTERP_HEAT ) then
          srad(:,:,2)  = srad(:,:,3)
        end if
        if ( USE_FLUXES ) then
          if ( INTERP_STRESS ) then
            ustr(:,:,2) = ustr(:,:,3)
            vstr(:,:,2) = vstr(:,:,3)
          end if
          if ( INTERP_HEAT ) then
            dlrad(:,:,2) = dlrad(:,:,3)
            lheat(:,:,2) = lheat(:,:,3)
            lrad(:,:,2)  = lrad(:,:,3)
            sheat(:,:,2) = sheat(:,:,3)
            srad(:,:,2)  = srad(:,:,3)
          end if
        end if
        if ( USE_BULK ) then
          if ( INTERP_TAIR  ) temp(:,:,2)  = temp(:,:,3)
          if ( INTERP_PRES  ) pres(:,:,2)  = pres(:,:,3)
          if ( INTERP_CLOUD ) cloud(:,:,2) = cloud(:,:,3)
          if ( INTERP_HUMID ) humid(:,:,2) = humid(:,:,3)
          if ( INTERP_SST   ) sst(:,:,2)   = sst(:,:,3)
        end if
      end if

!      record(1) = record(1) + 1

!      if ( is_master ) print *, "II", uwnd(50,50,2), uwnd(50,50,3) &
!                              , a, ADVANCE_REC_INT, ADVANCE_REC

      call read_all( ADVANCE_REC_INT, 3, year, record )
      call read_all( ADVANCE_REC    , 1, year, record )

      if ( READ_WIND ) then

        if ( interp_wind ) then
          call linint_vec( uwnd(:,:,2), vwnd(:,:,2)  &
                         , uwnd(:,:,3), vwnd(:,:,3)  &
                         , a                         &
                         , uwnd(:,:,1), vwnd(:,:,1)  &
                         )
        end if

! TODO: Are these surface arrays even needed?
        uwsrf = uwnd(:,:,1)
        vwsrf = vwnd(:,:,1)

        if ( .not.CALC_SWR ) then
          swrad = ( 1. - a ) * srad(:,:,2) + a * srad(:,:,3)
        end if

! Calculate wind stress
        if ( USE_BULK ) then
          call wind_to_stress( uwsrf, vwsrf, wusurf, wvsurf, 1 )
          call calculate_fluxes
        end if

      end if

      if ( USE_FLUXES ) then
!        if ( INTERP_STRESS ) then
!          wusurf = ( 1. - a ) * ustr(:,:,2) + a * ustr(:,:,3)
!          wvsurf = ( 1. - a ) * vstr(:,:,2) + a * vstr(:,:,3)
        call wind_to_stress( uwsrf, vwsrf, wusurf, wvsurf, 1 )
!        else
!          wusurf = ustr(:,:,1)
!          wvsurf = vstr(:,:,1)
!        end if
        if ( INTERP_HEAT ) then
          wtsurf = ( 1. - a ) * ( lrad(:,:,2)+sheat(:,:,2)+lheat(:,:,2) )  &
                 +        a   * ( lrad(:,:,3)+sheat(:,:,3)+lheat(:,:,3) )
          swrad  = ( 1. - a ) * srad(:,:,2) + a * srad(:,:,3)
        else
          wtsurf = lrad(:,:,1) + sheat(:,:,1) + lheat(:,:,1)
          swrad  = srad(:,:,1)
        end if
      end if

! If solar radiation does not penetrate water
      if ( nbct == 1 ) wtsurf = wtsurf + swrad

! Simple parameterisation
      swrad  =  swrad*(1.-icec)
      wtsurf = wtsurf*(1.-icec)! + itsurf*icec ! [TODO] itsurf is unstable after spinup for some reason
      wssurf = wssurf*(1.-icec)
!      print *, "WT:", minval(wtsurf), maxval(wtsurf)
!      print *, "WS:", minval(wssurf), maxval(wssurf)
!      print *, "IT:", minval(itsurf), maxval(itsurf)
!      print *, "SW:", minval(swrad) , maxval(swrad)
!      print *, "Ci:", minval(icec)  , maxval(icec)

      if ( TAPER_BRY ) call taper_forcing

      if ( .false. ) call river_flux

! Relax surface to climatology
      call relax_surface( wssurf, wtsurf, sss, sst )


    end ! subroutine step
!
!______________________________________________________________________
!
    subroutine taper_forcing
!----------------------------------------------------------------------
!  Primitive surface fluxes tapering at the border
!----------------------------------------------------------------------
! called by: step [air]
!______________________________________________________________________
!
      use model_run, only: iint, iend

      implicit none

      real(rk) ramp


      ramp = 1.
      if ( USE_RAMP ) ramp = iint/iend

      swrad  = swrad *taper_mask
      wssurf = wssurf*taper_mask
      wtsurf = wtsurf*taper_mask
      wusurf = ramp*wusurf*taper_mask
      wvsurf = ramp*wvsurf*taper_mask


    end ! subroutine
!
!______________________________________________________________________
!
    subroutine wind_to_stress( uwnd, vwnd, ustr, vstr, mode )
!----------------------------------------------------------------------
!  Calculates momentum flux (m^2/s^2) from wind velocity (m/s).
!----------------------------------------------------------------------
! called by: step  [air]
!______________________________________________________________________
!
      use glob_const , only: rhow => rhoref
      use glob_domain, only: im, jm
      use glob_ocean , only: t, u, v
      use seaice     , only: icec

      implicit none

      integer                   , intent(in   ) :: mode
      real(rk), dimension(im,jm), intent(in   ) :: uwnd, vwnd
      real(rk), dimension(im,jm), intent(  out) :: ustr, vstr

      integer  i, j
      real(rk) cda, uvabs


      cda = .00125 ! Default value and compiler warning eliminator

      do j = 1,jm
        do i = 1,im

          uvabs = sqrt( (uwnd(i,j)-u(i,j,1))**2  &
                      + (vwnd(i,j)-v(i,j,1))**2 )

          select case ( mode )

! mpiPOM formula with high wind-speed limit (XQ Yin, LY Oey, 2007)
            case ( 1 )

              if (     uvabs <=  11. ) then
                cda = .0012_rk
              elseif ( uvabs <=  19. ) then
                cda = .00049_rk  + uvabs* .000065_rk
              elseif ( uvabs <= 100. ) then
                cda = .001364_rk + uvabs*(.0000234_rk - uvabs*2.31579e-7_rk)
              else
                cda = .00138821_rk
              end if

! Hellerman and Rosenstein
            case ( 2 )

              cda = t(i,j,1)-temp(i,j,1)
              cda = .934e-3 + uvabs*(  .788e-4 - uvabs*.616e-6   &
                                     - .214e-5*cda )             &
                            +   cda*(  .868e-4 -   cda*.12e-5 )
              if ( cda < .00125 ) cda = .00125

          end select

!  Assume that free moving ice packs even at 10/10 concentrations pass at least half the momentum from atmosphere.
! TODO: If info about fast ice is available, damp any momentum flux over such areas to zero.
          ustr(i,j) = -rhoa/rhow*cda*uvabs * (uwnd(i,j)-u(i,j,1)) *(1.-icec(i,j)*0.5)
          vstr(i,j) = -rhoa/rhow*cda*uvabs * (vwnd(i,j)-v(i,j,1)) *(1.-icec(i,j)*0.5)

        end do
      end do


    end ! subroutine wind_to_stress
!
!______________________________________________________________________
!
    subroutine calculate_fluxes
!----------------------------------------------------------------------
!  Generates and interpolates all surface fluxes.
!----------------------------------------------------------------------
! called by: step         [air]
!
! calls    : bulk_fluxes  [air]
!______________________________________________________________________
!
      implicit none


      if ( READ_BULK ) then

        if ( interp_cloud ) then
          cloud(:,:,1) = ( 1. - a ) * cloud(:,:,2) + a * cloud(:,:,3)
        end if
        if ( interp_humid ) then
          humid(:,:,1) = ( 1. - a ) * humid(:,:,2) + a * humid(:,:,3)
        end if
        if ( interp_pres ) then
          pres(:,:,1) = ( 1. - a ) * pres(:,:,2) + a * pres(:,:,3)
        end if
        if ( interp_tair ) then
          temp(:,:,1) = ( 1. - a ) * temp(:,:,2) + a * temp(:,:,3)
        end if
!        print *, "CLOUD: ", maxval(cloud(:,:,1))
!        print *, "HUMID: ", maxval(humid(:,:,1))
!        print *, "PRESS: ", maxval(pres(:,:,1))
!        print *, "TEMPR: ", maxval(temp(:,:,1))

        e_atmos = 0.01 * ( pres(:,:,1) - 1013. )

! Calculate fluxes
!        if ( USE_BULK ) then
          call bulk_fluxes
!        end if

      end if


    end ! subroutine calculate_fluxes
!
!______________________________________________________________________
!
    subroutine bulk_fluxes
!----------------------------------------------------------------------
! Based on: http://pelagos.oc.phys.uoa.gr/mfstep/alermo/bulk_dis2.f
!
! Report: http://pelagos.oc.phys.uoa.gr/mfstep/alermo/airsea_report.pdf
!
!  This subroutine provides surface boundary conditions for momentum,
! heat and salt equations solved by the hydrodynamic model in cases
! when the atmospheric model provides cloud cover data instead of the
! net solar radiation flux and the downward longwave radiation flux as
! it is assumed in the first version of bulk code. The net solar
! radiation is calculated according to the Reed formula while the net
! longwave radiation can be calculated according to Bignami, May,
! Herzfeld or Berliand formula (see `LWRAD_FORMULA` flag below).
!
!  Momentum, heat and freshwater fluxes are calculated from the
! atmospheric parameters (wind velocity, air temperature, relative
! humidity, precipitation, cloud coverage) provided by the weather
! prediction model and the model's sst using proper air-sea bulk
! formulae (Castellari et al., 1998, Korres and Lascaratos 2003).
! ( Castellari et al., 1998. Journal of Marine Systems, 18, 89-114 ;
!   Korres and Lascaratos, 2003. Annales Geophysicae, 21, 205-220.)
!
!  All units are S.I. (M.K.S.)
!
!  The user has the option to calculate the net longwave radiation flux
! according to Bignami, May, (pseudo)Herzfeld, and Berliand formula.
! This is done through the `LWRAD_FORMULA` parameter (0 - 3)
!
!  This subroutine provides its output into the arrays:
!
! 1. WUSURF(): X-component of wind stress divided by (-1)*rho
! 2. WVSURF(): Y-component of wind stress divided by (-1)*rho
! 3. SWRAD() : Net solar radiation heat flux divided by (-1)*rho*Cpw
! 4. WTSURF(): Net upward heat flux divided by rho*Cpw
! 5. PME()   : precipitation minus evaporation rate
! ( RHO: sea water density, Cpw: specific heat capacity of seawater )
!
!
!  This subroutine needs the following input:
!
!
! Model related data
! 1. IM     : NUMBER OF GRID POINTS IN X
! 2. JM     : NUMBER OF GRID POINTS IN Y
! 3. TBIAS  : CONSTANT VALUE SUBTRACTED FROM MODEL'S TEMPERATURE FIELD
! 4. FSM()  : THE MODEL MASK (1.:SEA GRID POINT, 0.:LAND GRID POINT)
! 5. TSURF(): MODEL'S TEMPERATURE FIELD AT THE TOP VERTICAL LEVEL AND AT THE
!             CENTRAL TIME LEVEL
! 6. ALON() : LOGNITUDE OF GRID POINTS
! 7. ALAT() : LATITUDE OF GRID POINTS
! 9. IYR    : INTEGER VALUE CORRESPONDING TO CURRENT YEAR(for example iyr=2002)
! 10. IMO   : INTEGER VALUE CORRESPONDING TO CURRENT MONTH (1-->12)
! 11. IDAY  : INTEGER VALUE CORRESPONDING TO CURRENT DAY (1-->31)
! 12. IHOUR : INTEGER VALUE CORRESPONDING TO CURRENT HOUR (0-->23)
! 13: IMIN  : INEGER VALUE CORRESPONDING TO CURRENT MINUTES (0--59)
!
! ATMOSPHERIC DATA:
! 6.  UAIR() : X-COMPONENT OF AIR VELOCITY (in m/s) AT 10m ABOVE SEA SURFACE
! 7.  VAIR() : Y-COMPONENT OF AIR VELOCITY (in m/s) AT 10m ABOVE SEA SURFACE
! 8.  TAIR() : AIR TEMPERATURE (in deg Kelvin) AT 2m ABOVE SEA SURFACE
! 9.  RHUM() : RELATIVE HUMIDITY (%) AT 2m ABOVE SEA SURFACE
! 10. RAIN() : PRECIPITATION RATE (in m/s)
! 11. CLOUD()  : CLOUD COVERAGE IN TENTHS (0.-->1.)
! 12. PRES() : ATMOSPHERIC PRESSURE AT SURFACE (hPa)
!
!
! Important:
! SUBROUTINE BULK REQUIRES THAT THE ATMOSPHERIC DATA ARE ALREADY
! INTERPOLATED IN TIME AND SPACE (i.e. mapped onto the model grid and
! interpolated to the current time step. The user has to write his/her own
! subroutines in order to map the raw atmospheric data onto the model grid
! and interpolate them to the current time step)
!______________________________________________________________________
!
      use glob_const , only: C2K, DEG2RAD, rhoref, Rho_Cpw, rk
      use glob_domain, only: im, jm !    , my_task
      use glob_grid  , only: east_e, fsm, north_e
      use glob_ocean , only: rho, s, t, u, v
      use config     , only: tbias, sbias
      use model_run  , only: dtime

      implicit none

      real(rk) unow, vnow, tnow, pnow, precip, cld  &
             , sst_model, QBW, lwrd

!      real(rk), external :: cd, heatlat, esk

      integer  i, j
      real(rk) ce2, ch2, const               &
             , deltemp                       &
             , ea12, esatair, esatoce, evap  &
             , fe, fh, humnow                &
             , Qe, Qh, Qu                    &
             , rhom, rnow                    &
             , sol_net, sp, ss, sstk, stp    &
             , tnowk                         &
             , usrf, vsrf                    &
             , wair, wflux, wsatair, wsatoce
      real(rk), dimension(im,jm) :: pme ! Precipitation minus evaporation [m/s]


      const = expsi/Ps ! 0.622/1013.

      do j = 1,jm
      do i = 1,im

        if ( fsm(i,j) == 0. ) cycle

        unow      = uwnd(i,j,1)
        vnow      = vwnd(i,j,1)
        tnow      = temp(i,j,1)
        pnow      = pres(i,j,1)
        rnow      =  rho(i,j,1)*rhoref + 1000._rk
        humnow    = humid(i,j,1)
        precip    = rain(i,j,1)/1000._rk ! rwnd: precipitation rate from kg/(m2*s) to m/s
        cld       = maxval((/0._rk,cloud(i,j,1)/100._rk/)) ! rwnd: total cloud cover from % to tenths
        sst_model = t(i,j,1) + tbias
        lwrd      = dlrad(i,j,1)

        if ( i < im ) then
          usrf = .5_rk*(u(i+1,j,1)+u(i,j,1))
        else
          usrf = .5_rk*(u(i-1,j,1)+u(i,j,1))
        end if

        if ( j < jm ) then
          vsrf = .5_rk*(v(i,j+1,1)+v(i,j,1))
        else
          vsrf = .5_rk*(v(i,j-1,1)+v(i,j,1))
        end if

! SST_MODEL IS THE MODEL'S TEMPERATURE (in deg Celsius) AT THE TOP LEVEL
!
! --- compute wind speed magnitude for bulk coeff.
!
        SP = sqrt(unow*unow+vnow*vnow)

!
! --- SST data converted in Kelvin degrees
! --- TNOW is already in Kelvin
!
        sstk  = sst_model + C2K
        tnowk = tnow      + C2K
!
!
! ---calculates the Saturation Vapor Pressure at air temp. and at sea temp.
! ---esat(Ta) , esat(Ts)
!
        esatair = bucksat(tnow     ,pnow)
        esatoce = bucksat(sst_model,pnow)
!
! --- calculates the saturation mixing ratios at air temp. and sea temp.
! --- wsat(Ta) , wsat(Ts)
!
        wsatair = (expsi/pnow) * esatair
        wsatoce = (expsi/pnow) * esatoce
!
! --- calculates the mixing ratio of the air
! --- w(Ta)
!
        wair = .01_rk * humnow * wsatair
!            if ( isnan(wair) ) then
!              print *, "[ WAIR ]"
!              print *, "humnow", humnow
!              print *, "wsatair", wsatair
!              print *, "pnow", pnow
!              print *, "esatair", esatair
!              stop
!            end if
!
! --- calculates the density of  moist air
!
        rhom = 100._rk*(pnow/Rd)*(expsi*(1.+wair)/(tnowk*(expsi+wair)))
!            if ( isnan(rhom) ) then
!              print *, "[ RHOM ]"
!              print *, "pnow", pnow
!              print *, "Rd", Rd
!              print *, "expsi", expsi
!              print *, "wair", wair
!              print *, "tnowk", tnowk
!              stop
!            end if
!
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the net longwave radiation flux at the sea surface (QBW)
! according to Bignami (Bignami et al., 1995) or May formula (May,1986)
!
! Bignami et al., 1995: Longwave radiation budget in the Mediterranean
! Sea, J.Geophys Res., 100, 2501-2514.
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
        ea12 = .01_rk*humnow*esatair

        select case ( LWRAD_FORMULA )

          case ( lwBIGNAMI )
            QBW = .98_rk*sigma*sstk**4 - sigma*tnowk**4  &
                 *(  .653_rk + .00535_rk*ea12    )       &
                 *( 1._rk    + .1762_rk *cld*cld )

          case ( lwMAY )
            if ( cld < 0._rk ) cld = 0._rk
            QBW = ( 1._rk - .75_rk*(cld**3.4_rk) )                 &
                 *( sigma*(tnowk**4)*( .4_rk - .05_rk*sqrt(ea12) ) &
                  + 4._rk*sigma*(tnowk**3)*(sstk-tnowk) )

          case ( lwHERZFELD )
            QBW = ( sigma*.96_rk*( 1._rk-(.92e-5_rk*tnowk*tnowk) )*tnowk**4  &
                +4._rk*sigma*.96_rk*( C2K+temp(i,j,1)**3 )*(sstk-tnowk) )  &
                 *( 1._rk - .75_rk*cos(north_e(i,j)*DEG2RAD)*cld) ! cos(phi) here is an improvised `beta` coefficient as a function of latitude from Herzfeld

          case default  ! lwBERLIAND - Berliand (1952)
            QBW = .97_rk*sigma*                               &
                 (     tnowk**4 * (.39_rk-.05_rk*sqrt(ea12))  &
                  *( 1._rk - .6823_rk*cld*cld )               &
                  + 4._rk*tnowk**3 * (sstk-tnowk) )

        end select

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- calculate the term : ( Ts - Ta )
        deltemp = sstk - tnowk

        sol_net = 0._rk
        if ( CALC_SWR ) then
! --- Calculate net solar radiation flux according to Reed (Reed,1977) formula
!
! Reed, R.K.,1977: On estimating insolation over the ocean, J.Phys.
! Oceanogr. 17, 854-871.

          sol_net = sol_rad( cld, east_e(i,j)*DEG2RAD            &
                           , north_e(i,j)*DEG2RAD, dtime%year    &
                           , dtime%month, dtime%day, dtime%hour  &
                           , dtime%min )

! --- 1. Divide net solar radiation flux by rho*Cpw and reverse sign
          swrad(i,j) = -sol_net/rho_cpw
        end if
!

        if ( USE_COARE ) then

          call coare( (sqrt((unow-usrf)**2+(vnow-vsrf)**2)),      &
                      10._rk, temp(i,j,1), 2._rk, humnow, 2._rk,  &
                      pnow,sst_model,sol_net,lwrd,north_e(i,j),   &
                      600._rk, (3.6e6_rk*precip), 0._rk, 0._rk,   &
                      QH, QE, Evap, ch2, ce2 )
!            Evap = Evap*rho/3.6e6
!            if (my_task==1.and.i==50.and.j==50) then
!!            print *, "3.5 EVAP: ", Evap
!            write(61, '(6(f12.7,x),f12.7)') QE,QH,QBW,sol_net,
!     &                           (qbw+qh+qe-sol_net),lwrd,3.6e6*precip
!            end if

!            call coare30((unow-usurf(i,j)),(vnow-vsurf(i,j)),
!     &                    10._rk, tair(i,j), 2._rk, rhnow, 2._rk,
!     &                    pnow,sst_model,rnow,cld,precip*1000.,sol_net,
!     &                   -QBW, QH, QE, Evap )

!            if (my_task==1.and.i==50.and.j==50) then
!!            print *, "3.5 EVAP: ", Evap
!            write(62, '(6(f12.7,x),f12.7)') QE,QH,QBW,sol_net,
!     &                           (qbw+qh+qe-sol_net),lwrd,precip*1000.
!            end if
!            if (i==50.and.j==50) then
!            print *, "3.0 EVAP: ", Evap
!            print *, "3.0 QE:   ", QE
!            print *, "3.0 QH:   ", QH
!            print *, "RAIN===   ", precip
!            end if

        else
! Calculate turbulent exchange coefficients according to Kondo scheme
! ( Kondo, J., 1975. Boundary Layer Meteorology, 9, 91-112)

! ---- Kondo scheme
          ss = 0.
          fe = 0.
          fh = 0.
          if ( sp > .0 ) ss = deltemp/(sp**2.)
!
! --- calculate the Stability Parameter :
!
          stp = ss*abs(ss) / (abs(ss)+.01_rk)
!
! --- for stable condition :
          if ( ss < 0. ) then
            if ( (stp > -3.3) .and. (stp < 0.) ) then
              fh = .1_rk + .03_rk*stp + .9_rk*exp(4.8_rk*stp)
              fe = fh
            else
              if ( stp <= -3.3_rk ) then
                fh = 0.
                fe = fh
              end if
            end if
!
! --- for unstable condition :
          else
            fh = 1._rk + .63_rk*sqrt(stp)
            fe = fh
          end if

          ch2 = 1.3e-3_rk*fh
          ce2 = 1.5e-3_rk*fe

!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the sensible (QH) latent (QE) heat flux
!                    and the evaporation rate (EVAP)
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
          QH = rhom*cp*ch2*sp*deltemp
!
! --- calculates the term : esat(Ts)-r*esat(Ta)
!
          EVAP = esatoce - humnow*.01_rk*esatair
!
! --- calculate the term : Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
! --- Evaporation rate [kg/(m2*sec)]
!
          EVAP = rhom*ce2*sp*evap*const
!
! --- calculate the LATENT HEAT FLUX  QE in MKS ( watt/m*m )
! --- QE = L*rhom*Ce*|V|*[esat(Ts)-r*esat(Ta)]0.622/1013
!
          QE = evap*heat_latent(sst_model)

        end if
!
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the water flux (WFLUX) in m/sec
!---- -- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
        WFLUX =  evap/rhoref - precip

        pme(i,j)= -wflux
!            if ( isnan(pme(i,j)) ) then
!              print *, "[[[[[[]]]]]]"
!              print *, "pme: ", pme(i,j)
!              print *, "wflux: ", wflux
!              print *, "evap: ", evap
!              print *, "precip: ", precip
!              stop
!            end if
!
! Important note for Princeton Ocean Model users:
! THE SALT FLUX ( WSSURF() ) IN POM MAIN CODE (REQUIRED FOR PROFT) SHOULD
! BE CALCULATED AS:
!       do 3072 j = 1, jm
!       do 3072 i = 1, im
!  3072 WSSURF(I,J)=PME(I,J)*(VF(I,J,1)+SBIAS)

!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate  the net upward flux (QU) at the sea surface
!---- --- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
! --- calculates : Qu = Qb + QH + QE
!
        QU = qbw + qh + qe
!
! --- 1. Divide upward heat flux by rho*Cpw
!
        wtsurf(i,j) = QU/rho_cpw
!

!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Calculate the wind stress components (TAUX, TAUY)
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
! --- Calculate  the Drag Coefficient Cd accoriding to
!                         Hellerman & Rosenstein (1983)
!
!            cd1 = cd(sp,deltemp)

! --- Calculate  the wind stresses in MKS ( newton/m*m )
! --- taux= rhom*Cd*|V|u     tauy= rhom*Cd*|V|v
!
!            TauX = rhom*cd1*sp*unow
!            TauY = rhom*cd1*sp*vnow

! --- Reverse Sign and divide by sea water density
!            wusurf(i,j) = -taux/rho
!            wvsurf(i,j) = -tauy/rho

!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
! Multiply all fluxes by model mask
!---- ------ ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
!
        wusurf(i,j) = wusurf(i,j)*fsm(i,j)
        wvsurf(i,j) = wvsurf(i,j)*fsm(i,j)
        wtsurf(i,j) = wtsurf(i,j)*fsm(i,j)
        wssurf(i,j) = pme(i,j)*(s(i,j,1)+sbias)*fsm(i,j) ! sb? vf?
!        if ( abs(wssurf(i,j)).gt.1.e-2 ) then
!          print *, wssurf(i,j), pme(i,j), s(i,j,1)
!          stop
!        end if

        if ( USE_DQDSST ) then ! From Roms_tools (Penven, Pierrick, et al. "Software tools for pre-and post-processing of oceanic regional simulations." Environmental Modelling & Software 23.5 (2008): 660-662.)
          wtsurf(i,j) = wtsurf(i,j)                                 &
                      + ( 4._rk*sigma*tnowk**3                      &
                        + rhom*ch2*SP*Cp                            &
                        + rhom*ce2*SP*heat_latent(sst_model)        &
                         *2353._rk*log(10._rk)/tnowk**2             &
                         *( humnow*2.541e4_rk*exp(-5415._rk/tnowk)  & ! 2.541e4 instead of 2.541e6 is to convert `humnow` from % to fraction
                           *.620689655_rk )                         &
                        ) * (sst_model-sst(i,j,1)) / rho_cpw
        end if

      end do
      end do


    end ! subroutine bulk_fluxes
!
!______________________________________________________________________
!
    subroutine coare( u, zu, t, zt, rh, zq, P, ts, Rs, Rl, lat        &
                        ,zi, rain, cp, sigH, hsb, hlb, Evap, ch, ce )
!----------------------------------------------------------------------
! Input:
!
!     u = relative wind speed (m/s) at height zu(m)
!     t = bulk air temperature (degC) at height zt(m)
!    rh = relative humidity (%) at height zq(m)
!     P = surface air pressure (mb) (default = 1015)
!    ts = water temperature (degC)
!    Rs = downward shortwave radiation (W/m^2) (default = 150)
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude
!    zi = PBL height (m) (default = 600m)
!  rain = rain rate (mm/hr)
!    cp = phase speed of dominant waves (m/s)
!  sigH =  significant wave height (m)
!
! jcool=0 if surface temperature is true surface skin temperature,
!         otherwise ts is assumed the bulk temperature.
!
! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the
!  COARE algorithm, J. Climate, 16, 571-590.
!
      implicit none

      integer(1), parameter   :: jcool = 1

      real(rk), intent(in   ) :: u,zu, t,zt, rh,zq, P, ts       &
                               , Rs,Rl, lat, zi, rain, cp,sigH
      real(rk), intent(  out) :: Ce, Ch

      real(rk)  A, a1, a2, ad, Al, alq                             &
              , B, Bd, be, Bf, Beta, bigc                          &
              , CC, Cd, Cd10, cdhf, Ch10                           &
              , charn, charnC, charnS, charnW, cpa, cpv            &
              , cpw, cqhf, Ct, Ct10, cthf                          &
              , dels, dq, dqer, dt, dter                           &
              , Evap, fdg, gf, grav                                &
              , hbb, hlb, hlwebb, hsb, hsbb                        &
              , L, L10, Le                                  &
              , Pv, Q, qcol, qout       &
              , Qs, qsr                          &
              , Rgas, rhoa, rhodry, rhow           &
              , Ribcu, Ribu, Rnl, Rns, rr        &
              , ta, tau, tcw, tdk, tkt  &
              , tsr, tssr, tvsr                  &
              , u10, u10N, ug, umax     &
              , usr, ut, visa, visw, von, wbar, wetc        &
              , xlamx, zet, zetu, zo, zo10, zoq, zoS, zot, zot10   &
              , zref

      integer(1) nits, i


      logical waveage, seastate

      waveage  = .false.
      seastate = .false.

      L = 1.  ! compiler warning bypass

! convert rh to specific humidity (the functions below return g/kg)
      Qs = qsat26sea(ts,P)/1000._rk       ! surface water specific humidity [kg/kg]
      call qsat26air(t,P,rh, Q, Pv)       ! specific humidity of air [kg/kg]
      Q = Q/1000._rk

!-----------  set constants ----------------------------------------------
      zref = 10.   _rk
      Beta =  1.2  _rk                    ! Given as 1.25 in Fairall et al.(1996)
      von  =   .4  _rk                    ! Von Karman's "number"
      fdg  = 1.    _rk                    ! Turbulent Prandtl number
      tdk  = 273.16_rk
      grav = grv(lat)

!-----------  air constants ----------------------------------------------
      Rgas   = 287.1_rk                            ! Gas const. dry air [J/kg/K]
      Le     = ( 2.501_rk - .00237_rk*ts )*1.e6_rk ! Latent Heat of vaporization [J/kg]
      cpa    = 1004.67_rk                          ! Specific heat of dry air [J/kg/K] (Businger 1982)
      cpv    = cpa*( 1._rk + .84_rk*Q )            ! Moist air - currently not used (Businger 1982)
      rhoa   =  P    *100._rk                      &
              /( Rgas*(t+tdk) * (1._rk+.61_rk*Q) ) ! Moist air density [kg/m^3]
      rhodry = (P-Pv)*100._rk                      &
              /( Rgas*(t+tdk)                    ) ! Dry air density [kg/m^3]
      visa   = 1.326e-5_rk*( 1._rk+t*6.542e-3_rk   &
                            +    t*t*8.301e-6_rk   &
                            -  t*t*t*4.840e-9_rk ) !  Kinematic viscosity of dry air [m2/s], Andreas (1989).

!-----------  cool skin constants  ---------------------------------------
      Al   = 2.1e-5_rk*( ts+3.2_rk )**.79_rk       ! Water thermal expansion coef.
      be   =  .026_rk                              ! Salinity expansion coef.
      cpw  = 4000._rk                              ! Specific heat of water [J/kg/K]
      rhow = 1025._rk                              ! Seawater density
      visw = 1.e-6_rk                              ! Kinematic viscosity of water [m^2/s]
      tcw  =  .6  _rk                              ! Thermal conductivity of water [W/m/K]
      bigc = 16._rk*grav*cpw*(rhow*visw)**3 / (tcw**2 * rhoa**2)
      wetc =  .622_rk*Le*Qs / ( Rgas*(ts+tdk)**2 ) ! correction for dq;slope of sat. vap.

!-----------  net radiation fluxes ---------------------------------------
      Rns =  .945_rk * Rs                          ! albedo correction
      Rnl =  .97_rk*( 5.67e-8_rk*( ts-.3_rk*jcool+tdk )**4 - Rl) ! initial value

!----------------  begin bulk loop --------------------------------------------

!-----------  first guess ------------------------------------------------
!        du = u-us                          ! u is already relative to surface speed
      dt = ts-t-.0098_rk*zt
      dq = Qs-Q
      ta = t+tdk
      ug = .5_rk
      dter  = .3_rk
      dqer  = 0._rk
      ut    = sqrt( u**2 + ug**2 )
      u10   = ut * log(10._rk/1.e-4_rk) / log(zu/1.e-4_rk)
      usr   = .035_rk * u10
      zo10  = .011_rk*usr**2/grav + .11_rk*visa/usr
      Cd10  = ( von / log(10._rk/zo10) )**2
      Ch10  = .00115_rk
      Ct10  = Ch10 / sqrt(Cd10)
      zot10 = 10._rk / exp(von/Ct10)
      Cd    = ( von / log(zu/zo10) )**2
      Ct    =   von / log(zt/zot10)
      CC    =   von * Ct/Cd
      Ribcu = -250._rk*zu/(zi*Beta**3) ! 1./.004 = 250. I just like multiplication better
      Ribu  = -grav*zu/ta*( (dt-dter*jcool)+.61_rk*ta*dq )/ut**2
      if ( Ribu < 0. ) then
        zetu = CC*Ribu/( 1._rk +       Ribu/Ribcu )
      else
        zetu = CC*Ribu*( 1._rk + 3._rk*Ribu/CC    )
      end if
      L10 = zu/zetu                         ! Monin-Obukhov length
      gf  = ut/u
      usr = ut*von / ( log(zu/zo10) - psiu_40(zu/L10) )
      tsr = -(dt-     dter*jcool)*von*fdg                          &
                                 /(log(zt/zot10)-psit_26(zt/L10))
      qsr = -(dq-wetc*dter*jcool)*von*fdg                          &
                                 /(log(zq/zot10)-psit_26(zq/L10))
      tkt = .001_rk

!----------------------------------------------------------
!  The following gives the new formulation for the
!  Charnock variable
!----------------------------------------------------------
      charnC =   .011 _rk
      umax   = 22.    _rk ! 19.
      a1     =   .0016_rk !   .0017
      a2     =-  .0035_rk !-  .0050

      charnC = a1*u10 + a2
      if ( u10 > umax ) charnC = a1*umax+a2
      if ( charnC < .011 ) charnC = .011_rk

      A = .114_rk   ! wave-age dependent coefficients
      B = .622_rk

      Ad= .091_rk   ! Sea-state/wave-age dependent coefficients
      Bd= 2.  _rk

      charnW = A*(usr/cp)**B
      zoS    = sigH*Ad*(usr/cp)**Bd
      charnS = zoS*grav/usr/usr

      charn = .011_rk !*ones(N,1);
      if ( ut > 18. ) then
        charn = .018_rk
      elseif ( ut > 10. ) then
        charn = .011_rk + (ut-10._rk)/8._rk*.007_rk ! ~/8. was ~/(18.-10.) and ~*.007 was ~*(.018-0.011)
      end if
!        charnC = charn ! Fix?

      nits = 10   ! number of iterations

      if ( zetu > 50. ) nits = 1
!--------------  bulk loop --------------------------------------------------

      do i = 1, nits

        zet = von*grav*zu/ta*(tsr+.61_rk*ta*qsr)/(usr**2)
!          zet = von*grav*zu*(tsr*(1.+.61*Q))
        if (waveage) then
          if (seastate) then
            charn = charnS
          else
            charn = charnW
          end if
        else
          charn = charnC
        end if
        L  = zu/zet
        zo = charn*usr**2/grav + .11_rk*visa/usr      ! surface roughness
!          if (zo<1.d-10) zo = 1.d-10
        rr = zo*usr/visa
!        zoq= min(1.6e-4_rk, 5.8e-5_rk/rr**.72)        ! These thermal roughness lengths give Stanton and
!        zot= zoq                                      ! Dalton numbers that closely approximate COARE 3.0
        zot = min( 1.00e-4_rk/rr**.55_rk, 2.4e-4_rk/rr**1.2_rk )
        zoq = min( 1.15e-4_rk           , 5.5e-5_rk/rr** .6_rk )
        cdhf = von    /(log(zu/zo) - psiu_26(zu/L))
        cqhf = von*fdg/(log(zq/zoq)- psit_26(zq/L))
        cthf = von*fdg/(log(zt/zot)- psit_26(zt/L))
        usr  = ut*cdhf
        qsr  =-(dq-wetc*dter*jcool)*cqhf
        tsr  =-(dt-     dter*jcool)*cthf
        tvsr = tsr+.61_rk*ta*qsr
        tssr = tsr+.51_rk*ta*qsr
        Bf   =-grav/ta*usr*tvsr
        ug   = .2_rk
        if ( Bf > 0. ) ug = max(.2_rk, Beta*(Bf*zi)**.333_rk )
        ut = sqrt(u**2 + ug**2)
        gf = ut/u
        hsb=-rhoa*cpa*usr*tsr
        hlb=-rhoa*Le*usr*qsr
        qout = Rnl+hsb+hlb
        dels = Rns*( .065_rk + 11._rk*tkt - 6.6e-5_rk/tkt  &
                    *( 1._rk-exp(-1250._rk*tkt) ) )        ! 1./0.0008 = 1250.
        qcol = qout-dels
        alq  = Al*qcol + be*hlb*cpw/Le
        if ( alq > 0. ) then
          xlamx = 6._rk / ( 1._rk+(bigc*alq/usr**4)**.75_rk )**.333_rk
          tkt   = xlamx*visw/(sqrt(rhoa/rhow)*usr)
        else
          xlamx= 6._rk
          tkt = min(.01_rk, xlamx*visw/(sqrt(rhoa/rhow)*usr))
        end if
        dter = qcol*tkt/tcw
        dqer = wetc*dter
        Rnl  = .97_rk*( 5.67e-8_rk*(ts-dter*jcool+tdk)**4 - Rl ) ! update dter
        u10N = usr/von/gf*log(10._rk/zo)
        charnC = a1*u10N+a2
        if ( u10N > umax ) charnC = a1*umax+a2
        charnW = A*(usr/cp)**B
        zoS = sigH*Ad*(usr/cp)**Bd - .11_rk*visa/usr
        charnS = zoS*grav/usr/usr
        if ( charnC < .011 ) charnC = .011_rk

      end do

!----------------  compute fluxes  --------------------------------------------
      tau  =  rhoa*usr*usr/gf       ! wind stress [N/m^2]
      hsb  = -rhoa*cpa*usr*tsr      ! sensible heat flux [W/m^2]
      hlb  = -rhoa*Le*usr*qsr       ! latent heat flux [W/m^2]
      hbb  = -rhoa*cpa*usr*tvsr     ! buoyancy flux
      hsbb = -rhoa*cpa*usr*tssr     ! sonic heat flux
      wbar = 1.61_rk*hlb/Le/(1._rk+1.61_rk*Q)/rhoa + hsb/rhoa/cpa/ta
      hlwebb = hlb + rhoa*wbar*Q*Le
      Evap   = hlwebb/Le          ! evaporation rate [kg/m^2/s]
!        hlb = hlb + hlwebb ! ?????? Webb correction to latent heat flux already in ef via zoq/rr function so return hlwebb

!-----  compute transfer coeffs relative to ut @ meas. ht  --------------------
      Cd = tau/rhoa/ut/max(.1_rk,u)
      Ch =-usr*tsr/ut/(dt-dter*jcool)
      Ce =-usr*qsr/ut/(dq-dqer*jcool)


    end !subroutine coare
!
!______________________________________________________________________
!
    pure subroutine qsat26air(T,P,rh,q,em)
!----------------------------------------------------------------------
! computes saturation specific humidity [g/kg]
! given T [degC], rh [%], and P [mb]
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: T, P, rh
      real(rk), intent(out):: q, em


      em =  .01_rk* rh * bucksat(T,P)
      q  = 622._rk*em/( P - .378_rk*em )


    end !subroutine qsat26air
!
!______________________________________________________________________
!
    pure real(rk) function qsat26sea(T,P)
!----------------------------------------------------------------------
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: T, P
      real(rk) es


      es = .98_rk * bucksat(T,P) ! reduction at sea surface
      qsat26sea = 622._rk*es/( P - .378_rk*es )


    end !function qsat26sea
!
!______________________________________________________________________
!
    pure real(rk) function bucksat(T,P)
!----------------------------------------------------------------------
! computes saturation vapor pressure [mb]
! given T [degC] and P [mb]
!______________________________________________________________________
!
      implicit none

      real(rk), intent(in) :: T, P


      bucksat = 6.1121_rk * exp( 17.502_rk*T/(T+240.97_rk) )  &
                          * ( 1.0007_rk+3.46e-6_rk*P )


    end ! function bucksat
!
!______________________________________________________________________
!
    pure real(rk) function psit_26(zet)
!----------------------------------------------------------------------
! computes temperature structure function
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: zet

      real(rk) dzet, f, psik, psic, x


      dzet = min(50._rk, .35_rk*zet)         ! stable
      if ( zet < 0._rk ) then                  ! unstable
        x = (1._rk-15._rk  *zet)**.5_rk
        psik = 2._rk *log((1._rk+x     )/2._rk)
        x = (1._rk-34.15_rk*zet)**.3333_rk
        psic = 1.5_rk*log((1._rk+x+x**2)/3._rk)               &
             - sqrt(3._rk)*atan((1._rk+2._rk*x)/sqrt(3._rk))  &
             + 4._rk*atan(1._rk)/sqrt(3._rk)
        f = zet**2/(1._rk+zet**2)
        psit_26 = (1._rk-f)*psik + f*psic
      else
        psit_26 = -( (1._rk+.6667_rk*zet)**1.5_rk                  &
                + .6667_rk*(zet-14.28_rk)*exp(-dzet) + 8.525_rk )
      end if

    end ! function psit_26
!
!______________________________________________________________________
!
    pure real(rk) function grv(lat)
!----------------------------------------------------------------------
! computes g [m/sec^2] given lat in deg
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: lat

      real(rk), parameter  :: gamma = 9.7803267715_rk  &
                            , c1    =  .0052790414_rk  &
                            , c2    =  .0000232718_rk  &
                            , c3    =  .0000001262_rk  &
                            , c4    =  .0000000007_rk

      real(rk) phi, x


      phi = lat*pi/180._rk
      x   = sin(phi)**2
      grv = gamma*( 1._rk + x*c1 + x*x*c2 + x*x*x*c3 + x*x*x*x*c4 )


    end ! function grv
!
!______________________________________________________________________
!
    pure real(rk) function psiu_40(zet)
!----------------------------------------------------------------------
! computes velocity structure function
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: zet

      real(rk), parameter :: a = 1._rk        &
                           , b = 3._rk/4._rk  &
                           , c = 5._rk        &
                           , d =  .35_rk

      real(rk) dzet, f, psik, psic, x


      dzet = min(50._rk, .35_rk*zet)           ! stable
      psiu_40 = -(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
      if ( zet < 0._rk ) then                  ! unstable
        x = ( 1._rk-18._rk*zet)**.25_rk
        psik = 2._rk *log((1._rk+x     )/2._rk)  &
             +        log((1._rk+x*x   )/2._rk)  &
             - 2._rk*( atan(x)+atan(1._rk) )
        x = ( 1._rk-10._rk*zet)**.3333_rk
        psic = 1.5_rk*log((1._rk+x+x**2)/3._rk)               &
             - sqrt(3._rk)*atan((1._rk+2._rk*x)/sqrt(3._rk))  &
             + 4._rk*atan(1._rk)/sqrt(3._rk)
        f = zet**2 / ( 1._rk+zet**2 )
        psiu_40 = (1._rk-f)*psik + f*psic
      end if


    end !function psiu_40
!
!______________________________________________________________________
!
    pure real(rk) function psiu_26(zet)
!----------------------------------------------------------------------
! computes velocity structure function
!______________________________________________________________________
!
      implicit none

      real(rk), intent(in) :: zet

      real(rk), parameter :: a =  .7_rk       &
                           , b = 3._rk/4._rk  &
                           , c = 5._rk        &
                           , d =  .35_rk

      real(rk) dzet, f, psik, psic, x


      dzet = min(50._rk, .35_rk*zet)            ! stable
      if ( zet < 0._rk ) then                   ! unstable
        x = (1._rk-15._rk  *zet)**.25_rk
        psik = 2._rk *log((1._rk+x     )/2._rk)  &
             +        log((1._rk+x*x   )/2._rk)  &
             - 2._rk*( atan(x)+atan(1._rk) )
        x = (1._rk-10.15_rk*zet)**.3333_rk
        psic = 1.5_rk*log((1._rk+x+x**2)/3._rk)               &
             - sqrt(3._rk)*atan((1._rk+2._rk*x)/sqrt(3._rk))  &
             + 4._rk*atan(1._rk)/sqrt(3._rk)
        f = zet**2 / ( 1+zet**2 )
        psiu_26 = (1._rk-f)*psik + f*psic
      else
        psiu_26 =-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
      end if


    end !function psiu_26
!
!______________________________________________________________________
!
    pure real(rk) function RHcalc(T,P,Q)
!----------------------------------------------------------------------
!  Computes relative humidity given T [degC], P [hPa], and Q [
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in) :: T, P, Q

      real(rk) es, em


      es = bucksat(T,P)
      em = Q*P/(0.378_rk*Q+0.622_rk)
      RHcalc = 100._rk*em/es


    end ! function RHcalc
!!
!!______________________________________________________________________
!!
!    real(rk) function ESK(t)
!!----------------------------------------------------------------------
!!  Computes the saturation water vapor pressure from temperature (K)
!! (Lowe,1977; JAM,16,100,1977)
!!______________________________________________________________________
!!
!      implicit none
!
!      real(rk), intent(in) :: t
!
!      real(rk), dimension(7) :: a
!
!      data  a / 6984.505294              &
!              , -188.9039310             &
!              ,    2.133357675           &
!              ,   -1.288580973e-2_rk     &
!              ,    4.393587233e-5_rk     &
!              ,   -8.023923082e-8_rk     &
!              ,    6.136820929e-11_rk /
!
!
!      esk = a(1)+t*a(2)+t**2*a(3)+t**3*a(4)+t**4*a(5)+t**5*a(6)+t**6*a(7)
!
!
!    end ! function esk
!
!______________________________________________________________________
!
    real(rk) function sol_rad( acl,alon,alat,iyr,imt,idy,ihr,ime )
!--------------------------------------------------------------------
!  Computes shortwave solar radiation
!______________________________________________________________________
!
      implicit none


      real(rk), intent(in)  :: alat, alon, acl
      integer , intent(in)  :: iyr, imt, idy, ihr, ime

      dimension alpham(12),alb1(20),za(20),dza(19)

      integer imt1, iyr1, intT1, intT2, jab
      real(rk) :: albedo, alb1, alpha, alpham, aozone         &
                , bb1, bb2                                    &
                , capC, capG, capL, cosZen, DEC, dza, dZen    &
                , eclips, epsiln, g360, gha, gha360, RAD2DEG  &
                , SHA, SMLT, solar, SolAlt, SunBet, SunDec    &
                , tau, ThSun, TRM111, TRM112, TRM11, UT       &
                , qatten, qdiff, qdir, qtot, qzer             &
                , xl360, XLCT, yrdays, za, zen

      solar  = 1350.  _rk
      tau    =     .7 _rk
      aozone =     .09_rk
      yrdays =  365.  _rk
      alb1 = [ .719, .656, .603, .480, .385, .300, .250 &
             , .193, .164, .131, .103, .084, .071, .061 &
             , .054, .039, .036, .032, .031, .030       ]

      za = [ 90., 88., 86., 84., 82., 80., 78., 76., 74., 70. &
           , 66., 62., 58., 54., 50., 40., 30., 20., 10., 0.0 ]

      dza( 1: 8) =  2._rk
      dza( 9:14) =  4._rk
      dza(15:19) = 10._rk
!
! --- albedo monthly values from Payne (1972) as means of the values
! --- at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
! --- band of the Mediterranean Sea ) :
      alpham = [ .09_rk, .08_rk, .06_rk, .06_rk, .06_rk, .06_rk &
               , .06_rk, .06_rk, .06_rk, .07_rk, .09_rk, .10_rk ]

! --- albedo monthly values from Budyko (1963) at 50N
!      alpham = [ .16_rk, .12_rk, .09_rk, .07_rk, .07_rk, .06_rk &
!                ,.06_rk, .07_rk, .07_rk, .08_rk, .11_rk, .16_rk ]
!
!--------------------- calculations start -----------------------------
!
! --- sun hour angle :
!
      eclips = 23.439_rk*DEG2RAD
      RAD2DEG = 1._rk/DEG2RAD

      XLCT = ( 1._rk*ihr ) + ( 1._rk*ime / 60._rk )
      UT   = XLCT

      if ( imt > 2 ) then
        iyr1 = iyr
        imt1 = imt-3
      else
        iyr1 = iyr-1
        imt1 = imt+9
      end if

      intT1 = int(  30.6_rk * imt1      + .5_rk )
      intT2 = int( 365.25_rk*(iyr1-1976._rk)    )
      SMLT  = ( (UT/24._rk) + idy + intT1+intT2 - 8707.5_rk) / 36525._rk
      epsiln=  23.4393_rk -      .013_rk*SMLT
      capG  = 357.528_rk  + 35999.05_rk *SMLT

      if ( capG > 360._rk ) then
        g360 = capG - int( capG/360._rk )*360._rk
      else
        g360 = capG
      end if

      capC = 1.915_rk*sin(       g360*DEG2RAD )  &
           +  .02_rk *sin( 2._rk*g360*DEG2RAD )
      capL = 280.46_rk + 36000.77_rk*SMLT + capC

      if ( capL > 360._rk ) then
        xl360 = capL - int( capL/360._rk ) *360._rk
      else
        xl360 = capL
      end if

      alpha = xl360 - 2.466_rk*sin( 2._rk*xl360*DEG2RAD )  &
                    +  .053_rk*sin( 4._rk*xl360*DEG2RAD )
      gha = 15._rk*UT - 180._rk- capC + xl360 - alpha

      if ( gha > 360._rk ) then
        gha360 = gha - int( gha/360._rk )*360._rk
      else
        gha360 = gha
      end if

      DEC = atan( tan( epsiln*DEG2RAD ) * sin( alpha*DEG2RAD ) )*RAD2DEG

!     Calculate Solar Hour Angle
      ThSun = ( GHA360 + alon*RAD2DEG )*DEG2RAD
      SHA = GHA360 + ( alon*RAD2DEG )


! --- sun declination :
      SUNDEC = DEC * DEG2RAD

      TRM111 = sin( alat ) * sin( DEC*DEG2RAD )
      TRM112 =-cos( alat ) * cos( DEC*DEG2RAD )
      TRM11  = TRM111 - TRM112

! --- solar noon altitude in degrees :
      SolAlt = asin( TRM11 ) * RAD2DEG
      SunBet = SolAlt

!
! --- cosine of the solar zenith angle :
!
      coszen = sin(alat)*sin(sundec)+cos(alat)*cos(sundec)*cos(ThSun)

      if ( coszen <= 1.e-4_rk ) then
        coszen = 0._rk
        qatten = 0._rk
      else
        qatten = tau**(1._rk/coszen)
      end if
      qzer  = coszen * solar
      qdir  = qzer * qatten
      qdiff = ((1._rk-aozone)*qzer - qdir) * .5_rk
      qtot  =  qdir + qdiff
!
! --- ( radiation as from Reed(1977), Simpson and Paulson(1979) )
!
!
!-----------------------------------------------------------------------
! --- calculates the albedo as a function of the solar zenith angle :
! --- ( after Payne jas 1972 )
!-----------------------------------------------------------------------
!
! --- solar zenith angle in degrees :
!
      zen = RAD2DEG*acos(coszen)
!
      if     ( zen >= 74._rk ) then
        jab = int( .5 _rk*(90._rk-zen) +  1._rk )
      elseif ( zen >= 50._rk ) then
        jab = int( .23_rk*(74._rk-zen) +  9._rk )
      else
        jab = int( .1 _rk*(50._rk-zen) + 15._rk )
      end if
!
      dzen = ( za(jab) - zen ) / dza(jab)
!
      albedo=alb1(jab)+dzen*(alb1(jab+1)-alb1(jab))
!
! --- calculates SHORT WAVE FLUX ( watt/m*m )
! --- ( Rosati,Miyakoda 1988 ; eq. 3.8 )
!
!
      bb1 = .62 _rk     ! Reed
      bb2 = .636_rk     ! Isemer et al. (1989)
!
      sol_rad = qtot*(1._rk - bb1*acl + .0019_rk*sunbet)*(1._rk-albedo)
!        if (qshort.gt.qtot) qshort=qtot
      sol_rad = minval( [sol_rad, qtot] )


    end ! function sol_rad
!
!______________________________________________________________________
!
    pure real(rk) function heat_latent(t)
!----------------------------------------------------------------------
!  Calculates the Latent Heat of Vaporization ( J/kg ) as function of
! the temperature ( Celsius degrees )
! ( A. Gill  pag. 607 )
!
! --- Constant Latent Heat of Vaporization
!     L = 2.501e+6  (MKS)
!
      implicit none

      real(rk), intent(in)  :: t


      heat_latent = 2.5008e+6_rk - 2.3e+3_rk * t


    end ! function heat_latent
!
!______________________________________________________________________
!
    subroutine river_flux

      use glob_domain, only: im, jm, i_global, j_global
      use glob_grid  , only: art
      use glob_ocean , only: tsurf
      use model_run  , only: dtime

      implicit none

      real(rk), parameter, dimension(12) :: tumen_dis = (/ 30.2463253  &
                                                         , 44.2542105  &
                                                         ,105.8621385  &
                                                         ,156.2726807  &
                                                         ,357.9148493  &
                                                         ,526.118025   &
                                                         ,584.7622891  &
                                                         ,988.0466263  &
                                                         ,593.8361866  &
                                                         ,277.2579819  &
                                                         ,125.0181446  &
                                                         , 50.4105421  &
                                                        /)             &
! Below is amur discharge estimated from precipitation
!                                          , amur_dis  = (/  1982.  &
!                                                         ,  1324.  &
!                                                         ,  1058.  &
!                                                         ,  3231.  &
!                                                         , 14094.  &
!                                                         , 15948.  &
!                                                         , 15553.  &
!                                                         , 19291.  &
!                                                         , 20813.  &
!                                                         , 16596.  &
!                                                         ,  6162.  &
!                                                         ,  2441.  &
!                                                        /)
! Below is Amur discharge estimated from Mozherovsky
                                          , amur_dis  = (/  2480.  &
                                                         ,  1900.  &
                                                         ,  1620.  &
                                                         ,  3130.  &
                                                         , 15580.  &
                                                         , 17370.  &
                                                         , 15890.  &
                                                         , 18970.  &
                                                         , 22350.  &
                                                         , 17900.  &
                                                         ,  6850.  &
                                                         ,  2750.  &
                                                        /)
!      integer, parameter :: tumen_i = 7, tumen_j = 108
      integer, parameter :: tumen_i = 5, tumen_j = 59 ! jes403
!      integer, parameter :: amur_i = 31, amur_j  = 233 ! jes403
      integer, parameter :: amur_i = 36, amur_j  = 389 ! tat001
      integer i,j


!      if (     i_global(1) <= tumen_i .and. i_global(im) >= tumen_i  &
!         .and. j_global(1) <= tumen_j .and. j_global(jm) >= tumen_j ) then
!         i = tumen_i - i_global(1) + 1
!         j = tumen_j - j_global(1) + 1
!         vfluxf(i,j) = -tumen_dis( dtime%month ) / art(i,j)
!         tsurf(i,j) = max(temp(i,j,1),0.)
!!         wtsurf(i,j) = vfluxf(i,j)*max(temp(i,j,1),0.)
!         print *, "Tumen discharge: ", vfluxf(i,j), wtsurf(i,j),temp(i,j,1)
!      end if

      if (     i_global(1) <= amur_i .and. i_global(im) >= amur_i  &
         .and. j_global(1) <= amur_j .and. j_global(jm) >= amur_j ) then
         i = amur_i - i_global(1) + 1
         j = amur_j - j_global(1) + 1
         vfluxf(i,j) = -amur_dis( dtime%month ) / art(i,j)
         tsurf(i,j) = max(temp(i,j,1),0._rk)
!         wtsurf(i,j) = vfluxf(i,j)*max(temp(i,j,1),0.)
         print *, "Amur discharge: ", vfluxf(i,j), wtsurf(i,j),temp(i,j,1)
      end if


    end ! subroutine
!
!______________________________________________________________________
!
    subroutine linint_vec( u1, v1, u2, v2, a, u_int, v_int )

      use glob_domain, only: im, jm

      implicit none

      real(rk), dimension(im,jm), intent(in   ) :: u1, u2, v1, v2
      real(rk), dimension(im,jm), intent(  out) :: u_int, v_int
      real(rk)                  , intent(in   ) :: a
      real(rk), dimension(im,jm)                :: dir, spd


      u_int = ( 1. - a ) * u1 + a * u2
      v_int = ( 1. - a ) * v1 + a * v2
! Do not interpolate wind linearly with rough temporal resolution or quite variable wind direction.
! First, get the "correct" direction from linear interpolation, and then recalculate the absolute value.
! The obvious caveat is when w(n+1) and w(n) are perfectly in opposite to one another the resulted direction defaults to zero.
      dir = atan2(v_int,u_int)
!          where ( dir /= 0. ) ! TODO: Catch opposites with a condition, maybe?
      spd = ( 1. - a ) * sqrt( u1*u1 + v1*v1 )  &
           +       a   * sqrt( u2*u2 + v2*v2 )
      u_int = spd * cos( dir )
      v_int = spd * sin( dir )


    end ! subroutine linint_vec
!
!______________________________________________________________________
!
    subroutine read_all( execute, n, year, record )

      use glob_ocean, only: tb

      implicit none

      logical              , intent(in) :: execute
      integer              , intent(in) :: n
      integer, dimension(3), intent(in) :: record, year

      integer        ncid
      character(128) desc


      if ( .not. execute ) return

      if ( n >= 2 ) then
        write(desc,'("Reading interp. forcing record #",i4," @ ",i4)') &
            record(n), year(n)
      else
        write(desc,'("Reading surface forcing record #",i4," @ ",i4)') &
            record(1), year(1)
      end if

      call msg_print("", 1, desc)

! Read wind file
      if ( read_wind ) then

        ncid = file_open_nc( wind_path, year(n) )

        if ( ncid /= -1 ) then

          if ((     interp_wind .and. n>=2) .or.       &
              (.not.interp_wind .and. n==1)      ) then
            if ( read_var_nc( uwnd_name, uwnd(:,:,n), record(n), ncid ) /= 0  &
            .or. read_var_nc( vwnd_name, vwnd(:,:,n), record(n), ncid ) /= 0 ) then
              uwnd(:,:,n) = 0.
              vwnd(:,:,n) = 0.
              call msg_print("", 2, "Wind is defaulted to zero.")
            end if
          end if

          call check( file_close_nc( ncid ), "nf_close" )

        end if

      end if

! Read bulk file
      if ( read_bulk ) then

        ncid = file_open_nc( bulk_path, year(n) )

        if ( ncid /= -1 ) then

          if ((     interp_humid .and. n>=2) .or.       &
              (.not.interp_humid .and. n==1)      ) then
            if ( read_var_nc( humid_name, humid(:,:,n), record(n), ncid ) /= 0 ) then
              humid(:,:,n) = .85
              call msg_print("", 2, "Humidity is defaulted to 85%")
            end if
          end if

          if ((     interp_rain .and. n>=2) .or.       &
              (.not.interp_rain .and. n==1)      ) then
            if ( read_var_nc( rain_name,  rain(:,:,n), record(n), ncid ) /= 0 ) then
              rain(:,:,n) = 0.
              call msg_print("", 2, "Precip.rate is defaulted to zero.")
            end if
          end if

          if ((     interp_pres .and. n>=2) .or.       &
              (.not.interp_pres .and. n==1)      ) then
            if ( read_var_nc( pres_name,  pres(:,:,n), record(n), ncid ) /= 0 ) then
              pres(:,:,n) = 1013.
              call msg_print("", 2, "Atm.pressure is defaulted to 1013 hPa.")
            end if
          end if

          if ((     interp_sst .and. n>=2) .or.       &
              (.not.interp_sst .and. n==1)      ) then
            if ( read_var_nc( sst_name,   sst(:,:,n), record(n), ncid ) /= 0 ) then
              sst(:,:,n) = tb(:,:,1)
              call msg_print("", 2, "SST is defaulted to surf. level temp.")
            end if
          end if

          if ((     interp_tair .and. n>=2) .or.       &
              (.not.interp_tair .and. n==1)      ) then
            if ( read_var_nc( tair_name,  temp(:,:,n), record(n), ncid ) /= 0 ) then
              temp(:,:,n) = tb(:,:,1)
              call msg_print("", 2, "Tair is defaulted to surf. level temp.")
            end if
          end if

          if ( READ_CLOUD ) then
          if ((     interp_cloud .and. n>=2) .or.       &
              (.not.interp_cloud .and. n==1)      ) then
            if ( read_var_nc( tcld_name, cloud(:,:,n), record(n), ncid ) /= 0 ) then
              cloud(:,:,n) = 0.
              call msg_print("", 2, "Cloud cover is defaulted to zero.")
            end if
          end if
          else
            cloud(:,:,n) = 0.
          end if

          call check( file_close_nc( ncid ), "nf_close" )

        end if

      end if ! if READ_BULK

! Read flux file
      if ( read_stress .or. read_heat ) then

        ncid = file_open_nc( flux_path, year(n) )

        if ( read_stress .and. ncid /= -1 ) then

          if ((     interp_stress .and. n>=2) .or.       &
              (.not.interp_stress .and. n==1)      ) then
            if ( read_var_nc( ustr_name, ustr(:,:,n), record(n), ncid ) /= 0  &
            .or. read_var_nc( vstr_name, vstr(:,:,n), record(n), ncid ) /= 0 ) then
              ustr(:,:,n) = 0.
              vstr(:,:,n) = 0.
              call msg_print("", 2, "Wind stress is defaulted to zero.")
            end if
          end if

        end if

        if ( read_heat .and. ncid /= -1 ) then

          if ((     interp_heat .and. n>=2) .or.       &
              (.not.interp_heat .and. n==1)      ) then
            if ( read_var_nc( dlrad_name, dlrad(:,:,n), record(n), ncid ) /= 0 ) then
              dlrad(:,:,n) = 370.
              call msg_print("", 2, "Downlonrad. is set to 370 W/m^2.")
            end if
            if ( read_var_nc( lheat_name, lheat(:,:,n), record(n), ncid ) /= 0 ) then
              lheat(:,:,n) = 0.
              call msg_print("", 2, "Latent heat is set to zero.")
            end if
            if ( read_var_nc(  lrad_name,  lrad(:,:,n), record(n), ncid ) /= 0 ) then
              lrad(:,:,n) = 0.
              call msg_print("", 2, "Net longrad. is set to zero.")
            end if
            if ( read_var_nc( sheat_name, sheat(:,:,n), record(n), ncid ) /= 0 ) then
              sheat(:,:,n) = 0.
              call msg_print("", 2, "Sensible heat is set to zero.")
            end if
            if ( read_var_nc(  srad_name,  srad(:,:,n), record(n), ncid ) /= 0 ) then
              srad(:,:,n) = 0.
              call msg_print("", 2, "Net shortrad. is set to zero.")
            end if
          end if

        end if
        
        if ( ncid /= -1 ) then
          call check( file_close_nc( ncid ), "nf_close" )
        end if

      end if ! if READ_FLUX


    end ! subroutine read_all
!
!______________________________________________________________________
!
    pure character(256) function get_filename( path, year )
!----------------------------------------------------------------------
!  Constructs filename string in `<path>YYYY<FORMAT_EXT>` format.
!______________________________________________________________________
!
      implicit none

      character(*), intent(in) :: path
      integer     , intent(in) :: year


      if ( path(len(trim(path)):len(trim(path))) == "." ) then
        write( get_filename, '( a, i4.4, a )' ) trim(path)      &
                                               , year            &
                                               , trim(FORMAT_EXT)
      else
        get_filename = path
      end if


    end ! function get_filename
!
!
!= I/O SECTION ========================================================
!______________________________________________________________________
!
    integer(1) function read_var_nc( var_name, var, record, ncid )
!----------------------------------------------------------------------
!  Read a variable (NC format).
!______________________________________________________________________
!
      use glob_const , only: C2K, rhoref, rho_cpw, rk
      use glob_domain
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_3d

      integer                   , intent(in   ) :: ncid, record
      real(rk), dimension(im,jm), intent(  out) :: var
      character(*)              , intent(in   ) :: var_name

      integer                  status, varid
      integer(MPI_OFFSET_KIND) start(4), edge(4)
      character(64)            units


      if ( var_name == '' ) then
        read_var_nc = -1
        return
      end if

      read_var_nc = 0

      start = 1
      edge  = 1

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

! convert data if necessary
      status = nf90mpi_get_att( ncid, varid, "units", units )
      if (var_name=='precip') print *, "!!!!", units
      if ( status == NF_NOERR ) then
        select case ( trim(units) )
          case ( "K" )
            var = var - C2K
          case ( "Pa" )
            var = var/100.
          case ( "W m**-2", "W m-2", "W/m^2", "W/m2" )
            var = var/rho_cpw
          case ( "N m**-2", "N m-2", "N/m^2", "N/m2" )
            var = -var/rhoref
          case ( "mm/(3hours)", "mm/3hr", "mm/(3hr)")
            var = var/10800000._rk
        end select
      end if


    end ! subroutine read_var_nc
!
!___________________________________________________________________
!
    integer function file_open_nc( path, year )
!-------------------------------------------------------------------
!  Opens netcdf file for reading.
!___________________________________________________________________
!
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INFO_NULL
      use pnetcdf    , only: nf90mpi_open, NF_NOERR, NF_NOWRITE

      implicit none

      integer     , intent(in) :: year
      character(*), intent(in) :: path

      integer        status
      character(256) filename, netcdf_file


      filename = get_filename( path, year )
      netcdf_file = trim(filename)
      status = nf90mpi_open( POM_COMM, netcdf_file, NF_NOWRITE   &
                           , MPI_INFO_NULL, file_open_nc )
      if ( status /= NF_NOERR ) then
        call msg_print("", 2, "Failed to open `"//trim(filename)//"`")
        file_open_nc = -1
      end if


    end ! function file_open_nc
!
!___________________________________________________________________
!
    integer function file_close_nc( ncid )
!-------------------------------------------------------------------
!  Opens netcdf file for reading.
!___________________________________________________________________
!
      use pnetcdf, only: nf90mpi_close

      implicit none

      integer, intent(in) :: ncid


      file_close_nc = nf90mpi_close( ncid )


    end ! function file_close_nc
!
!______________________________________________________________________
!
    subroutine check(status, routine)
!----------------------------------------------------------------------
!  Checks for NetCDF I/O error and exits with an error message if hits.
!______________________________________________________________________
!
      use glob_domain, only: is_master
      use pnetcdf    , only: nf90mpi_strerror, NF_NOERR

      implicit none

      integer     , intent(in) :: status
      character(*), intent(in) :: routine


      if ( status /= NF_NOERR ) then
        if ( is_master ) then
          print '(/a,a)', 'IO error at module `AIR`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
!          stop
        end if
      end if


    end ! subroutine check


end module air
