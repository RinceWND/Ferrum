!______________________________________________________________________
!
! Module `TIDE` (tide.f90)
!----------------------------------------------------------------------
!  Module for applying tides.
!
!  Author  : RinceWND
!  Created : 2019-02-11
!______________________________________________________________________
!
module tide

  use module_time
  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  public

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  type(date) tide_offset

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED   & ! Is set according to the external flag `use_air`.
                    , TPXO_FILE    ! Original TPXO Netcdf file used.

  real   , private :: a         ! time-interpolation factor

  character(10)                           &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(PATH_LEN)  & ! Full paths (with filenames) to:
    tide_el_path       & !  atmospheric parameters file
  , tide_uv_path         !  surface fluxes file

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(len=VAR_LEN)  &
      cons_name           & ! Downward longwave radiation
  , el_amp_name           & ! Downward longwave radiation
  , el_pha_name           & ! Latent heat flux
  ,    lat_name           & ! Latitude at z-points
  ,  lat_u_name           & ! Latitude at u-points
  ,  lat_v_name           & ! Latitude at v-points
  ,    lon_name           & ! Longitude at z-points
  ,  lon_u_name           & ! Longitude at u-points
  ,  lon_v_name           & ! Longitude at v-points
  , ua_amp_name           & ! Longwave net radiation
  , ua_pha_name           & ! Relative humidity
  , va_amp_name           & ! Atm. pressure
  , va_pha_name             ! Precipitation rate


!----------------------------------------------------------------------
! Tide related arrays
!----------------------------------------------------------------------
  character(4)               &
       , allocatable         &
       , dimension(:)   ::   &
    constituent

  real(rk)                   &
       , allocatable         &
       , dimension(:,:) ::   &
    lat                      & ! Latitude at z-points
  , lat_u                    & ! Latitude at u-points
  , lat_v                    & ! Latitude at v-points
  , lon                      & ! Longitude at z-points
  , lon_u                    & ! Longitude at u-points
  , lon_v                    & ! Longitude at v-points
  , tide_el                  & ! Total cloud cover [fraction 0..1]
  , tide_el_b                &
  , tide_mask                &
  , tide_ua                  &
  , tide_ua_b                &
  , tide_va                  &
  , tide_va_b

  real(rk)                   &
       , allocatable         &
       , dimension(:,:,:) :: &
    el_amp                   & ! atmospheric pressure
  , el_pha                   & ! short wave radiation incident
  , ua_amp                   & ! wind speed in x-direction
  , ua_pha                   & ! volume flux through water column
  , va_amp                   & ! volume flux through water column
  , va_pha                     ! wind speed in y-direction

  integer                                 &
       , private                          &
       , parameter :: CONS_OVERALL = 104

  integer(8)           &
       , private   ::  &
    ie_lft             &
  , ie_rht             &
  , iu_lft             &
  , iu_rht             &
  , iv_lft             &
  , iv_rht             &
  , je_bot             &
  , je_top             &
  , ju_bot             &
  , ju_top             &
  , jv_bot             &
  , jv_top

  integer(8)           &
    ncons

  character(6)         &
       , dimension(12) :: active_cons

  type T_constituent
    character(6)              name   ! Name of constituent
    real(rk)                  freq   ! Constituent frequency
    real(rk)                  speed  ! Constituent angular speed
    integer                   k      ! Wave number
    real(rk), dimension(4) :: V      ! Astronomical constants
    real(rk), dimension(3) :: u      ! Astronomical constants
    real(rk), dimension(4) :: f      ! Reduction constants
  end type

  type(T_constituent), dimension(CONS_OVERALL) :: con


  contains

!______________________________________________________________________
!
    subroutine initialize_mod( config_file )
!----------------------------------------------------------------------
!  Initialize tide module.
!______________________________________________________________________
!
      use config     , only: use_tide
      use glob_domain, only: is_master

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/tide/                                        &
        active_cons, tide_el_path, tide_uv_path

      namelist/tide_vars/                                   &
        el_amp_name, el_pha_name,    lat_name,  lat_u_name  &
      ,  lat_v_name,    lon_name,  lon_u_name,  lon_v_name  &
      , ua_amp_name, ua_pha_name, va_amp_name, va_pha_name  &
      ,   cons_name


      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_TIDE ) DISABLED = .true.

      if ( DISABLED ) return

      TPXO_FILE = .true.

      tide_offset = str2date("1900-01-01_00:00:00")
      active_cons    = ""
      active_cons(1) = "m2"
      active_cons(2) = "s2"
      active_cons(3) = "n2"
      active_cons(4) = "k2"

! Initialize variables with their defaults
      tide_el_path = "in/tide/"
      tide_uv_path = "in/tide/"

      el_amp_name = "ha"
      el_pha_name = "hp"
         lat_name = "lat_z"
         lon_name = "lon_z"
       lat_u_name = "lat_u"
       lon_u_name = "lon_u"
       lat_v_name = "lat_v"
       lon_v_name = "lon_v"
      ua_amp_name = "ua"
      ua_pha_name = "up"
      va_amp_name = "va"
      va_pha_name = "vp"
        cons_name = "con"

! Override configuration
      open ( 73, file = config_file, status = 'old' )
      read ( 73, nml = tide )
      read ( 73, nml = tide_vars )
      close( 73 )

! Variables management
      pos = len(trim(tide_el_path))
      if ( tide_el_path(pos:pos) == "/" ) then
        tide_el_path = trim(tide_el_path)  &
                     //"h_tpxo9.v1"//FORMAT_EXT
      end if

      pos = len(trim(tide_uv_path))
      if ( tide_uv_path(pos:pos) == "/" ) then
        tide_uv_path = trim(tide_uv_path)  &
                     //"u_tpxo9.v1"//FORMAT_EXT
      end if

      if ( is_master ) then
        print *, "Tidal elevation: ", trim(tide_el_path)
        print *, "Tidal currents : ", trim(tide_uv_path)
      end if

! Allocate necessary arrays
      call allocate_arrays

! Generate tidal constituents:       name, period       , t,   [  a,  s,  h, p ]
!                                                              [ c1, c2, c3 ]
!                                                              [ c0, c1, c2, c3, pow ]
!  - semidiurnal:
      con( 1) = specify_constituent( "m2", 28.9841042_rk, 2, [   0.    , -2.    ,  2.    , 0.     ]      &
                                                           , [ - 2.14  ,  0.    ,  0.    ]               &
                                                           , [   1.0004, -0.0373,  0.0002,  0.    , 1. ] ) ! Principal lunar semidiurnal
      con( 2) = specify_constituent( "s2", 30.       _rk, 2, [   0.    ,  0.    ,  0.    , 0.     ]      &
                                                           , [   0.    ,  0.    ,  0.    ]               &
                                                           , [   1.    ,  0.    ,  0.    ,  0.    , 1. ] ) ! Principal solar semidiurnal
      con( 3) = specify_constituent( "n2", 28.4397295_rk, 2, [   0.    , -3.    ,  2.    , 1.     ]      &
                                                           , [ - 2.14  ,  0.    ,  0.    ]               &
                                                           , [   1.0004, -0.0373,  0.0002,  0.    , 1. ] ) ! Larger lunar elliptic semidiurnal
      con( 4) = specify_constituent( "k2", 30.0821373_rk, 2, [   0.    ,  0.    ,  2.    , 0.     ]      &
                                                           , [ -17.74  ,  0.68  , -0.04  ]               &
                                                           , [   1.0241,  0.2863,  0.0083, -0.0015, 1. ] ) ! Lunisolar declinational semidiurnal
!  - diurnal:
      con( 5) = specify_constituent( "o1", 13.9430356_rk, 1, [ 270.    , -2.    ,  1.    , 0.     ]      &
                                                           , [  10.8   , -1.34  ,  0.19  ]               &
                                                           , [   1.0089,  0.1871, -0.0147,  0.0014, 1. ] ) ! Lunar diurnal
      con( 6) = specify_constituent( "p1", 14.9589314_rk, 1, [ 270.    ,  0.    , -1.    , 0.     ]      &
                                                           , [   0.    ,  0.    ,  0.    ]               &
                                                           , [   1.    ,  0.    ,  0.    ,  0.    , 1. ] ) ! Solar diurnal
      con( 7) = specify_constituent( "q1", 13.3986609_rk, 1, [ 270.    , -3.    ,  1.    , 1.     ]      &
                                                           , [  10.8   , -1.34  ,  0.19  ]               &
                                                           , [   1.0089,  0.1871, -0.0147,  0.0014, 1. ] ) ! Larger lunar elliptic diurnal
      con( 8) = specify_constituent( "k1", 15.0410686_rk, 1, [  90.    ,  0.    ,  1.    , 0.     ]      &
                                                           , [ - 8.86  ,  0.68  , -0.07  ]               &
                                                           , [   1.006 ,  0.115 , -0.0088,  0.0006, 1. ] ) ! Lunisolar declinational diurnal
!  - shallow water:
      con( 9) = specify_constituent( "m4" ,  57.9682084_rk, 4, [   0.    , -4.    ,  4.    , 0.     ]      &
                                                             , [ - 4.28  ,  0.    ,  0.    ]               &
                                                             , [   1.0004, -0.0373,  0.0002,  0.    , 2. ] ) ! Shallow water overtides of principal lunar constituent, quarter-diurnal
      con(10) = specify_constituent( "m6" ,  86.9523127_rk, 6, [   0.    , -6.    ,  6.    , 0.     ]      &
                                                             , [ - 2.14  ,  0.    ,  0.    ]               &
                                                             , [   1.0004, -0.0373,  0.0002,  0.    , 1. ] ) ! Shallow water overtides of principal lunisolar constituent, sixth-diurnal
      con(11) = specify_constituent( "ms4",  58.9841042_rk, 4, [   0.    , -2.    ,  2.    , 0.     ]      &
                                                             , [  10.8   , -1.34  ,  0.19  ]               &
                                                             , [   1.0004, -0.0373,  0.0002,  0.    , 3. ] ) ! Shallow water lunisolar quarter-diurnal
      con(12) = specify_constituent( "m8" , 115.93642  _rk, 8, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water lunar eighth-diurnal
      con(13) = specify_constituent( "s8" , 120.       _rk, 8, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water solar eighth-diurnal
      con(14) = specify_constituent( "k3" ,  43.31412  _rk, 3, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water lunisolar declinational third-diurnal
      con(15) = specify_constituent( "m3" ,  43.4761563_rk, 3, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water lunar third-diurnal
!  - low-frequency:
      con(16) = specify_constituent( "msf", 1.0158958_rk, 2, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Lunisolar synodic fortnightly
      con(17) = specify_constituent( "mf" , 1.0980331_rk, 2, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Lunar fortnightly
      con(18) = specify_constituent( "mm" ,  .5443747_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Lunar monthly
      con(19) = specify_constituent( "msm",  .4716   _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Solar monthly
      con(20) = specify_constituent( "ssa",  .0821373_rk, 2, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Solar semiannual
      con(21) = specify_constituent( "sa" ,  .0410686_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Solar annual
!  - superimposed:
      con(22) = specify_constituent( "mu2" ,  27.9682084_rk, 1, [0.,0.,0.,0.]                 &
                                                              , [0.,0.,0.]                    &
                                                              , [1.0004,-0.0373,0.0002,0.,1.] ) ! Variational (M2,S2)
      con(23) = specify_constituent( "j1"  ,  15.5854433_rk, 1, [0.,0.,0.,0.]              &
                                                              , [0.,0.,0.]                 &
                                                              , [1.013,0.168,-0.017,0.,1.] ) ! Smaller lunar elliptic diurnal
      con(24) = specify_constituent( "nu2" ,  28.5125831_rk, 1, [0.,0.,0.,0.]                 &
                                                              , [0.,0.,0.]                    &
                                                              , [1.0004,-0.0373,0.0002,0.,1.] ) ! Larger lunar evectional
      con(25) = specify_constituent( "oo1" ,  16.1391017_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Lunar diurnal
      con(26) = specify_constituent( "2n2" ,  27.8953548_rk, 1, [0.,0.,0.,0.]                 &
                                                              , [0.,0.,0.]                    &
                                                              , [1.0004,-0.0373,0.0002,0.,1.] ) ! Lunar elliptical semidiurnal second-order
      con(27) = specify_constituent( "m1"  ,  14.4920521_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Smaller lunar elliptic diurnal ! or 14.496694
      con(28) = specify_constituent( "t2"  ,  29.9589333_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Larger solar elliptic
      con(29) = specify_constituent( "l2"  ,  29.5284789_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Smaller lunar elliptic semidiurnal
      con(30) = specify_constituent( "mn4" ,  57.4238337_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water quarter diurnal
      con(31) = specify_constituent( "pi1" ,  14.9178647_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(32) = specify_constituent( "2sm2",  31.0158958_rk, 1, [0.,0.,0.,0.]                 &
                                                              , [0.,0.,0.]                    &
                                                              , [1.0004,-0.0373,0.0002,0.,1.] ) ! Shallow water semidiurnal
      con(33) = specify_constituent( "phi1",  15.1232059_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(34) = specify_constituent( "2ms6",  87.9682084_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(35) = specify_constituent( "sn4" ,  58.4397295_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(36) = specify_constituent( "ksi1",  14.5695476_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(37) = specify_constituent( "mo3" ,  42.9271398_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(38) = specify_constituent( "2mn6",  96.407938 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(39) = specify_constituent( "mk3" ,  44.0251729_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water terdiurnal
      con(40) = specify_constituent( "msn6",  87.4238337_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(41) = specify_constituent( "2sm6",  88.9841042_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(42) = specify_constituent( "s9"  , 135.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallowwater overtides of the principal solar, ninth-diurnal
      con(43) = specify_constituent( "s7"  , 105.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallowwater overtides of the principal solar, seventh-diurnal
      con(44) = specify_constituent( "s6"  ,  90.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallowwater overtides of the principal solar, sixth-diurnal
      con(45) = specify_constituent( "s5"  ,  75.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallowwater overtides of the principal solar, fifth-diurnal
      con(46) = specify_constituent( "msk6",  89.066241 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(47) = specify_constituent( "2mk6",  88.0503457_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(48) = specify_constituent( "sk4" ,  60.0821373_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(49) = specify_constituent( "s4"  ,  60.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Shallow water overtides of principal solar constituent
      con(50) = specify_constituent( "mk4" ,  59.0662415_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(51) = specify_constituent( "sk3" ,  45.0410686_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(52) = specify_constituent( "so3" ,  43.943036 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(53) = specify_constituent( "kj2" ,  30.626512 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(54) = specify_constituent( "msn2",  30.5443747_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(55) = specify_constituent( "mks2",  29.0662415_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(56) = specify_constituent( "op2" ,  28.9019669_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(57) = specify_constituent( "mns2",  27.4238337_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(58) = specify_constituent( "oq2" ,  27.3416964_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(59) = specify_constituent( "so1" ,  16.0569644_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(60) = specify_constituent( "mp1" ,  14.0251729_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(61) = specify_constituent( "s3"  ,  45.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(62) = specify_constituent( "r2"  ,  30.0410667_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?

      con(63) = specify_constituent( "lamda2", 29.4556253_rk, 1, [0.,0.,0.,0.]                 &
                                                               , [0.,0.,0.]                    &
                                                               , [1.0004,-0.0373,0.0002,0.,1.] ) ! Smaller lunar evectional
      con(64) = specify_constituent( "theta1", 15.5125897_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(65) = specify_constituent( "psi1"  , 15.0821353_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(66) = specify_constituent( "s1"    , 15.       _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Solar diurnal
      con(67) = specify_constituent( "rho1"  , 13.4715145_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! Larger lunar evectional diurnal
      con(68) = specify_constituent( "sigma1", 12.9271398_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(69) = specify_constituent( "2q1"   , 12.8542862_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(70) = specify_constituent( "no1"   , 14.4966939_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(71) = get_constituent( "ksi1" )
      con(71) % name = "lp1"
      con(72) = get_constituent( "pi1" )
      con(72) % name = "tk1"
      con(73) = get_constituent( "psi1" )
      con(73) % name = "rp1"
      con(74) = get_constituent( "phi1" )
      con(74) % name = "kp1"
      con(75) = specify_constituent( "kq1", 16.6834764_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(76) = get_constituent( "mu2" )
      con(76) % name = "2ms2"
      con(77) = specify_constituent( "mq3", 42.3827651_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(78) = get_constituent( "mo3" )
      con(78) % name = "2mk3"
      con(79) = specify_constituent( "sp3"  ,  44.9589314_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(80) = specify_constituent( "2mns4",  56.407938 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(81) = specify_constituent( "3mk4" ,  56.8701754_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(82) = specify_constituent( "3ms4" ,  56.9523127_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(83) = specify_constituent( "2msk4",  57.8660711_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(84) = specify_constituent( "2mks4",  58.0503457_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(85) = specify_constituent( "3mn4" ,  58.5125831_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(86) = specify_constituent( "2smk4",  58.9019669_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(87) = specify_constituent( "2snm4",  59.4556253_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(88) = specify_constituent( "2msn4",  59.5284789_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(89) = specify_constituent( "2smn4",  60.5443747_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(90) = specify_constituent( "3sm4" ,  61.0158958_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(91) = specify_constituent( "2skm4",  61.0980331_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(92) = specify_constituent( "mnk6" ,  87.505971 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(93) = specify_constituent( "2mn8" , 114.8476674_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(94) = specify_constituent( "3mn8" , 115.3920422_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(95) = specify_constituent( "2msn8", 116.407938 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(96) = specify_constituent( "2mnk8", 116.4900753_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(97) = specify_constituent( "3ms8" , 116.9523127_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(98) = specify_constituent( "3mk8" , 117.0344499_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(99) = specify_constituent( "2smn8", 117.4238337_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(100)= specify_constituent( "msnk8", 117.505971 _rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(101)= specify_constituent( "2ms8" , 117.9682084_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(102)= specify_constituent( "2msk8", 118.0503457_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(103)= specify_constituent( "3sm8" , 118.9841042_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?
      con(104)= specify_constituent( "2smk8", 119.0662415_rk, 1, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] ) ! ?


      call msg_print("TIDE MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_air
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________
!
      use glob_domain, only: im, jm, POM_COMM, my_task
      use grid       , only:  east_e,  east_u,  east_v  &
                           , north_e, north_u, north_v
      use io
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      integer                                   file_id, varid, ndims, status
      integer(8)                                nx, ny
      integer(MPI_OFFSET_KIND), dimension(3) :: edge, start
      integer(MPI_OFFSET_KIND), dimension(:)  &
                              , allocatable  :: cons_dimlens
      real(rk) , allocatable, dimension(:,:) :: src_lat, src_lat_u, src_lat_v  &
                                              , src_lon, src_lon_u, src_lon_v  &
                                              , tmp


! Read files
      file_id = file_open( tide_el_path, NF90_NOWRITE )
      if ( file_id < 0 ) then
        call msg_print( "", 3, "Failed to open file `"//trim(tide_el_path)//"`" )
      end if

! If the original TPXO files are used, then interpolate into the model grid
      if ( TPXO_FILE ) then

        ndims = var_rank( file_id, el_amp_name )
        if     ( ndims > 3 ) then
          call msg_print( "", 1, "Input file variable rank is not 3." )
        elseif ( ndims < 3 ) then
          call msg_print( "TIDES DISABLED", 1, "Input file variable rank is less than 3." )
          DISABLED = .true.
          return
        end if
        allocate( cons_dimlens(ndims) )
        cons_dimlens = var_shape( file_id, el_amp_name )
        ncons = cons_dimlens(3)
        nx    = cons_dimlens(2)
        ny    = cons_dimlens(1)
        allocate( constituent(ncons) )

        start = 1
        edge  = [ ny, nx, ncons ]

        allocate( src_lat(nx,ny), src_lat_u(nx,ny), src_lat_v(nx,ny)    &
                , src_lon(nx,ny), src_lon_u(nx,ny), src_lon_v(nx,ny)    &
                , tmp(ny,nx) )

        status = var_read( file_id, trim(lat_name), tmp, start, edge )
        src_lat = transpose( tmp )
        status = var_read( file_id, trim(lon_name), tmp, start, edge )
        src_lon = transpose( tmp )

        file_id = file_close( file_id )


        file_id = file_open( tide_uv_path, NF90_NOWRITE )

        status = var_read( file_id, trim(lat_u_name), tmp, start, edge )
        src_lat_u = transpose( tmp )
        status = var_read( file_id, trim(lon_u_name), tmp, start, edge )
        src_lon_u = transpose( tmp )

        status = var_read( file_id, trim(lat_v_name), tmp, start, edge )
        src_lat_v = transpose( tmp )
        status = var_read( file_id, trim(lon_v_name), tmp, start, edge )
        src_lon_v = transpose( tmp )

        file_id = file_close( file_id )

  ! Crop domain and store these extents for further use
        ie_lft = maxloc(src_lon(:,1), 1, src_lon(:,1) < minval(east_e))
        ie_rht = minloc(src_lon(:,1), 1, src_lon(:,1) > maxval(east_e))
        iu_lft = maxloc(src_lon_u(:,1), 1, src_lon_u(:,1) < minval(east_u))
        iu_rht = minloc(src_lon_u(:,1), 1, src_lon_u(:,1) > maxval(east_u))
        iv_lft = maxloc(src_lon_v(:,1), 1, src_lon_v(:,1) < minval(east_v))
        iv_rht = minloc(src_lon_v(:,1), 1, src_lon_v(:,1) > maxval(east_v))
        ie_lft = max( ie_lft-1, 1_8 )
        ie_rht = min( ie_rht+1, nx  )
        iu_lft = max( iu_lft-1, 1_8 )
        iu_rht = min( iu_rht+1, nx  )
        iv_lft = max( iv_lft-1, 1_8 )
        iv_rht = min( iv_rht+1, nx  )
        if ( src_lat(1,1) < src_lat(1,ny) ) then
          je_bot = maxloc(src_lat(1,:), 1, src_lat(1,:) < minval(north_e))
          je_top = minloc(src_lat(1,:), 1, src_lat(1,:) > maxval(north_e))
          ju_bot = maxloc(src_lat_u(1,:), 1, src_lat_u(1,:) < minval(north_u))
          ju_top = minloc(src_lat_u(1,:), 1, src_lat_u(1,:) > maxval(north_u))
          jv_bot = maxloc(src_lat_v(1,:), 1, src_lat_v(1,:) < minval(north_v))
          jv_top = minloc(src_lat_v(1,:), 1, src_lat_v(1,:) > maxval(north_v))
          je_bot = max( je_bot-1, 1_8 )
          je_top = min( je_top+1, ny  )
          ju_bot = max( ju_bot-1, 1_8 )
          ju_top = min( ju_top+1, ny  )
          jv_bot = max( jv_bot-1, 1_8 )
          jv_top = min( jv_top+1, ny  )
        else
          print *, "[X] UNIMPLEMENTED : Inverted latitude"
          stop
        end if

  !      print *, "DST_lon: ", minval(east_e), ":", maxval(east_e)
  !      print *, "SRC_lon: ", src_lon(ie_lft,1), ":", src_lon(ie_rht,1)
  !      print *, "DST_lat: ", minval(north_e), ":", maxval(north_e)
  !      print *, "SRC_lat: ", src_lat(1,je_bot), ":", src_lat(1,je_top)

        deallocate( src_lat  , src_lon           &
                  , src_lat_u, src_lon_u         &
                  , src_lat_v, src_lon_v, tmp )

        allocate( lat( ie_rht-ie_lft+1, je_top-je_bot+1 )    &
                , lon( ie_rht-ie_lft+1, je_top-je_bot+1 )    &
                , lat_u( iu_rht-iu_lft+1, ju_top-ju_bot+1 )  &
                , lon_u( iu_rht-iu_lft+1, ju_top-ju_bot+1 )  &
                , lat_v( iv_rht-iv_lft+1, jv_top-jv_bot+1 )  &
                , lon_v( iv_rht-iv_lft+1, jv_top-jv_bot+1 ) )

!  If the tidal phases and amplitudes are already interpolated into the grid,
! read the number of constituents.
      else

        file_id = file_open( tide_el_path, NF90_NOWRITE )

! We presume that the processed tidal file is correct, having (im,jm,ncons) shape.
        allocate( cons_dimlens(3) )
        cons_dimlens = var_shape( file_id, el_amp_name )
        ncons = cons_dimlens(3)

        file_id = file_close( file_id )

      end if

! Allocate tidal arrays
      allocate(              &
        el_amp(im,jm,ncons)  &
      , el_pha(im,jm,ncons)  &
      , ua_amp(im,jm,ncons)  &
      , ua_pha(im,jm,ncons)  &
      , va_amp(im,jm,ncons)  &
      , va_pha(im,jm,ncons)  &
      , tide_el  (im,jm)     &
      , tide_el_b(im,jm)     &
      , tide_mask(im,jm)     &
      , tide_ua  (im,jm)     &
      , tide_ua_b(im,jm)     &
      , tide_va  (im,jm)     &
      , tide_va_b(im,jm)     &
       )

       el_amp  = 0.
       el_pha  = 0.
       ua_amp  = 0.
       ua_pha  = 0.
       va_amp  = 0.
       va_pha  = 0.
       tide_el = 0.
       tide_ua = 0.
       tide_va = 0.


    end ! subroutine allocate_mod
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  TODO: Fill up
!______________________________________________________________________
!
      use glob_domain, only: im, i_global, jm, j_global, my_task, POM_COMM
      use grid       , only: fsm, east_e, east_u, east_v, rot  &
                           , north_e, north_u, north_v
      use io
      use module_time
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      type(date), intent(in) :: d_in

      integer                                    file_id, varid, dimlen, i,j,k
      integer(MPI_OFFSET_KIND)                :: start(3), edge(3)
      real(rk), allocatable, dimension(:,:,:) :: tmp, var


! Quit if the module is not used.
      if ( DISABLED ) return

! Shrink free surface mask
      tide_mask = fsm
      allocate( tmp(im,jm,1) )
      do k = 1, 8
        tmp(:,:,1) = tide_mask
        do j = 2, jm-1
          do i = 2, im-1
            if ( (     tmp(i,j+1,1) < 1. .or. tmp(i+1,j,1) < 1.    &
                  .or. tmp(i,j-1,1) < 1. .or. tmp(i-1,j,1) < 1. )  &
                .and. tmp(i,j,1) == 1. ) then
              tide_mask(i,j) = 0. !real(k-1) / real(15-k)
            end if
          end do
        end do
        do j = 2, jm-1
          if ( tmp(1,j,1) == 1.                                    &
               .and. (     tmp(1,j+1,1) < 1. .or. tmp(2,j,1) < 1.  &
                      .or. tmp(1,j-1,1) < 1. ) ) then
            tide_mask(1,j) = 0. !real(k-1) / real(15-k)
          end if
          if ( tmp(im,j,1) == 1.                                       &
               .and. (     tmp(im,j+1,1) < 1. .or. tmp(im-1,j,1) < 1.  &
                      .or. tmp(im,j-1,1) < 1. ) ) then
            tide_mask(im,j) = 0. !real(k-1) / real(15-k)
          end if
        end do
        do i = 2, im-1
          if ( tmp(i,1,1) == 1.                                      &
               .and. (     tmp(i  ,2,1) < 1. .or. tmp(i+1,1,1) < 1.  &
                      .or. tmp(i-1,1,1) < 1. ) ) then
            tide_mask(i,1) = 0. !real(k-1) / real(15-k)
          end if
          if ( tmp(i,jm,1) == 1.                                       &
               .and. (     tmp(i,jm-1,1) < 1. .or. tmp(i+1,jm,1) < 1.  &
                      .or. tmp(i-1,jm,1) < 1. ) ) then
            tide_mask(i,jm) = 0. !real(k-1) / real(15-k)
          end if
        end do
      end do
      deallocate( tmp )

      if ( TPXO_FILE ) then

! Read rho-values (and constituent names)
        call check( nf90mpi_open( POM_COMM, tide_el_path, NF_NOWRITE   &
                                , MPI_INFO_NULL, file_id )                &
                  , "nf_open: "//tide_el_path )

        call check( nf90mpi_inq_varid( file_id, cons_name, varid )  &
                  , "nf_inq_varid: "//cons_name )
        call check( nf90mpi_get_var_all( file_id, varid, constituent )  &
                  , "nf_get_var: "//cons_name )
        print *, my_task, ":", constituent, shape(constituent)

! Read coordinates
        allocate( tmp( je_top-je_bot+1, ie_rht-ie_lft+1, 1 ) )

        start = (/ je_bot         , ie_lft         , 1_MPI_OFFSET_KIND /)
        edge  = (/ je_top-je_bot+1, ie_rht-ie_lft+1, 1_MPI_OFFSET_KIND /)

        call check( nf90mpi_inq_varid( file_id, lat_name, varid )  &
                  , "nf_inq_varid: "//lat_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lat_name )
        lat = transpose( tmp(:,:,1) )

        call check( nf90mpi_inq_varid( file_id, lon_name, varid )  &
                  , "nf_inq_varid: "//lon_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lon_name )
        lon = transpose( tmp(:,:,1) )
      
!      print *, "_DST_lon: ", minval(east_e), ":", maxval(east_e)
!      print *, "_SRC_lon: ", minval(lon), ":", maxval(lon)
!      print *, "_DST_lat: ", minval(north_e), ":", maxval(north_e)
!      print *, "_SRC_lat: ", minval(lat), ":", maxval(lat)

! Read amplitudes and phases
        deallocate( tmp )
        allocate( tmp( je_top-je_bot+1, ie_rht-ie_lft+1, ncons )    &
                , var( ie_rht-ie_lft+1, je_top-je_bot+1, ncons ) )
      
        start = (/ je_bot         , ie_lft         ,          1_MPI_OFFSET_KIND /)
        edge  = (/ je_top-je_bot+1, ie_rht-ie_lft+1, int(ncons,MPI_OFFSET_KIND) /)

        call check( nf90mpi_inq_varid( file_id, el_amp_name, varid )  &
                  , "nf_inq_varid: "//el_amp_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//el_amp_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( ie_rht-ie_lft+1, je_top-je_bot+1       &
                            , lon   , lat    , var(:,:,i)            &
                            , im             , jm                    &
                            , east_e, north_e, el_amp(:,:,i), 'V' )
        end do

        call check( nf90mpi_inq_varid( file_id, el_pha_name, varid )  &
                  , "nf_inq_varid: "//el_pha_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp          &
                                      , start, edge )  &
                  , "nf_get_var: "//el_pha_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( ie_rht-ie_lft+1, je_top-je_bot+1       &
                            , lon   , lat    , var(:,:,i)            &
                            , im             , jm                    &
                            , east_e, north_e, el_pha(:,:,i), 'P' )
        end do

        call check( nf90mpi_close( file_id ), "nf_close"//tide_uv_path )

! Read uv-values
        call check( nf90mpi_open( POM_COMM, tide_uv_path, NF_NOWRITE   &
                                , MPI_INFO_NULL, file_id )                &
                  , "nf_open: "//tide_uv_path )

! Read u-coordinates
        deallocate( tmp )
        allocate( tmp( ju_top-ju_bot+1, iu_rht-iu_lft+1, 1 ) )

        start = (/ ju_bot         , iu_lft         , 1_MPI_OFFSET_KIND /)
        edge  = (/ ju_top-ju_bot+1, iu_rht-iu_lft+1, 1_MPI_OFFSET_KIND /)

        call check( nf90mpi_inq_varid( file_id, lat_u_name, varid )  &
                  , "nf_inq_varid: "//lat_u_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lat_u_name )
        lat_u = transpose( tmp(:,:,1) )

        call check( nf90mpi_inq_varid( file_id, lon_u_name, varid )  &
                  , "nf_inq_varid: "//lon_u_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lon_u_name )
        lon_u = transpose( tmp(:,:,1) )

! Read u-amplitudes and phases
        deallocate( tmp, var )
        allocate( tmp( ju_top-ju_bot+1, iu_rht-iu_lft+1, ncons )    &
                , var( iu_rht-iu_lft+1, ju_top-ju_bot+1, ncons ) )
      
        start = (/ ju_bot         , iu_lft         ,          1_MPI_OFFSET_KIND /)
        edge  = (/ ju_top-ju_bot+1, iu_rht-iu_lft+1, int(ncons,MPI_OFFSET_KIND) /)

        call check( nf90mpi_inq_varid( file_id, ua_amp_name, varid )  &
                  , "nf_inq_varid: "//ua_amp_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//ua_amp_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( iu_rht-iu_lft+1, ju_top-ju_bot+1       &
                            , lon_u , lat_u  , var(:,:,i)            &
                            , im             , jm                    &
                            , east_u, north_u, ua_amp(:,:,i), 'V' )
        end do

        call check( nf90mpi_inq_varid( file_id, ua_pha_name, varid )  &
                  , "nf_inq_varid: "//ua_pha_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp          &
                                      , start, edge )  &
                  , "nf_get_var: "//ua_pha_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( iu_rht-iu_lft+1, ju_top-ju_bot+1       &
                            , lon_u , lat_u  , var(:,:,i)            &
                            , im             , jm                    &
                            , east_u, north_u, ua_pha(:,:,i), 'P' )
        end do

! Read v-coordinates
        deallocate( tmp )
        allocate( tmp( jv_top-jv_bot+1, iv_rht-iv_lft+1, 1 ) )

        start = (/ jv_bot         , iv_lft         , 1_MPI_OFFSET_KIND /)
        edge  = (/ jv_top-jv_bot+1, iv_rht-iv_lft+1, 1_MPI_OFFSET_KIND /)

        call check( nf90mpi_inq_varid( file_id, lat_v_name, varid )  &
                  , "nf_inq_varid: "//lat_v_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lat_v_name )
        lat_v = transpose( tmp(:,:,1) )

        call check( nf90mpi_inq_varid( file_id, lon_v_name, varid )  &
                  , "nf_inq_varid: "//lon_v_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//lon_v_name )
        lon_v = transpose( tmp(:,:,1) )

! Read v-amplitudes and phases
        deallocate( tmp, var )
        allocate( tmp( jv_top-jv_bot+1, iv_rht-iv_lft+1, ncons )    &
                , var( iv_rht-iv_lft+1, jv_top-jv_bot+1, ncons ) )

        start = (/ jv_bot         , iv_lft         ,          1_MPI_OFFSET_KIND /)
        edge  = (/ jv_top-jv_bot+1, iv_rht-iv_lft+1, int(ncons,MPI_OFFSET_KIND) /)

        call check( nf90mpi_inq_varid( file_id, va_amp_name, varid )  &
                  , "nf_inq_varid: "//va_amp_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp               &
                                      , start, edge )  &
                  , "nf_get_var: "//va_amp_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( iv_rht-iv_lft+1, jv_top-jv_bot+1       &
                            , lon_v , lat_v  , var(:,:,i)            &
                            , im             , jm                    &
                            , east_v, north_v, va_amp(:,:,i), 'V' )
        end do

        call check( nf90mpi_inq_varid( file_id, va_pha_name, varid )  &
                  , "nf_inq_varid: "//va_pha_name )
        call check( nf90mpi_get_var_all( file_id, varid, tmp          &
                                      , start, edge )  &
                  , "nf_get_var: "//va_pha_name )
! interpolate into the grid
        do i = 1, int(ncons)
          var(:,:,i) = transpose( tmp(:,:,i) )
          call interpolate_2D( iv_rht-iv_lft+1, jv_top-jv_bot+1       &
                            , lon_v , lat_v  , var(:,:,i)            &
                            , im             , jm                    &
                            , east_v, north_v, va_pha(:,:,i), 'P' )
        end do
      
        call check( nf90mpi_close( file_id ), "nf_close"//tide_uv_path )

        deallocate( tmp, var )

! Rotate uv-variables to conform to the grid using `tmp` and `var` as u and v
        allocate( tmp( im, jm, ncons )    &
                , var( im, jm, ncons ) )

! MUST NOT rotate phases, but the resulting tidal vectors... TODO: Does rotation of amplitude do the job?
        do i = 1, int(ncons)
          tmp(:,:,i) = ua_amp(:,:,i)*cos(rot) - va_amp(:,:,i)*sin(rot)
          var(:,:,i) = ua_amp(:,:,i)*sin(rot) + va_amp(:,:,i)*cos(rot)
        end do
        ua_amp = tmp/100.
        va_amp = var/100.

        ! do i = 1, int(ncons)
          ! tmp(:,:,i) = ua_pha(:,:,i)*cos(rot) - va_pha(:,:,i)*sin(rot)
          ! var(:,:,i) = ua_pha(:,:,i)*sin(rot) + va_pha(:,:,i)*cos(rot)
        ! end do
        ! where ( tmp >  180. ) tmp = tmp - 360.
        ! where ( tmp < -180. ) tmp = tmp + 360.
        ! where ( var >  180. ) var = var - 360.
        ! where ( var < -180. ) var = var + 360.
        ! ua_pha = tmp
        ! va_pha = var

        deallocate( tmp, var )

      else

        file_id = file_open( tide_el_path, NF90_NOWRITE )

        start = [ i_global(1), j_global(1),          1 ]
        edge  = [ im         , jm         , int(ncons) ]
        call check( var_read( file_id, el_amp_name, el_amp       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(el_amp_name) )
        call check( var_read( file_id, el_pha_name, el_pha       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(el_pha_name) )
        call check( var_read( file_id, ua_amp_name, ua_amp       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(ua_amp_name) )
        call check( var_read( file_id, ua_pha_name, ua_pha       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(ua_pha_name) )
        call check( var_read( file_id, va_amp_name, va_amp       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(va_amp_name) )
        call check( var_read( file_id, va_pha_name, va_pha       &
                            , start, edge )                      &
                  , "[tide]var_read:init "//trim(va_pha_name) )

        file_id = file_close( file_id )

      end if

      call step( d_in )

      call msg_print("TIDE INITIALIZED", 1, "")


    end ! subroutine init
!
!______________________________________________________________________
!
    subroutine step( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields during experiment.
!______________________________________________________________________
!
      use glob_const , only: DEG2RAD, GRAV
      use glob_domain, only: im, imm1, jm
      use grid       , only: fsm
      use module_time
      use model_run  , only: dti, iint, sec_of_year
!      use glob_ocean , only: d, el, elf

      implicit none

      type(date), intent(in) :: d_in

      integer                   i
      real(rk)                  f, h, N, p, s, tick, uv
      real(rk), dimension(4) :: tmp
      type(T_constituent)       this_con


! Quit if the module is not used.
      if ( DISABLED ) return

      ! tick = (d_in - tide_offset)/3600.
      tick = d_in % hour + d_in % min / 60._rk + d_in % sec / 3600._rk

      tide_el_b = tide_el
      tide_ua_b = tide_ua
      tide_va_b = tide_va
      tide_el = 0.
      tide_ua = 0.
      tide_va = 0.

      tmp = get_astronomy( d_in )
      h = tmp(1)
      s = tmp(2)
      p = tmp(3)
      N = tmp(4)

      do i = 1, int(ncons)
        this_con = get_constituent( constituent(i) )
        f  = this_con%f(1) + this_con%f(2)*cos(N) + this_con%f(3)*cos(2.*N) + this_con%f(4)*cos(3.*N)
        uv = this_con%V(1) + this_con%V(2)*s      + this_con%V(3)*h         + this_con%V(4)*p
        uv = uv            + this_con%u(1)*sin(N) + this_con%u(2)*sin(2.*N) + this_con%u(3)*sin(3.*N)
        tide_el(im,:) = tide_el(im,:) + f*el_amp(im,:,i)*cos( (this_con%speed*tick + uv - el_pha(im,:,i))*DEG2RAD )
        tide_ua(im,:) = tide_ua(im,:) + f*ua_amp(im,:,i)*cos( (this_con%speed*tick + uv - ua_pha(im,:,i))*DEG2RAD )
        tide_va(im,:) = tide_va(im,:) + f*va_amp(im,:,i)*cos( (this_con%speed*tick + uv - va_pha(im,:,i))*DEG2RAD )
        tide_el(2:imm1,jm) = tide_el(2:imm1,jm)                   &
                           + f*el_amp(2:imm1,jm,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - el_pha(2:imm1,jm,i) )*DEG2RAD )
        tide_ua(2:imm1,jm) = tide_ua(2:imm1,jm)                   &
                           + f*ua_amp(2:imm1,jm,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - ua_pha(2:imm1,jm,i) )*DEG2RAD )
        tide_va(2:imm1,jm) = tide_va(2:imm1,jm)                   &
                           + f*va_amp(2:imm1,jm,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - va_pha(2:imm1,jm,i) )*DEG2RAD )
        tide_el( 1,:) = tide_el( 1,:) + f*el_amp( 1,:,i)*cos( (this_con%speed*tick + uv - el_pha( 1,:,i))*DEG2RAD )
        tide_ua( 2,:) = tide_ua( 2,:) + f*ua_amp( 2,:,i)*cos( (this_con%speed*tick + uv - ua_pha( 2,:,i))*DEG2RAD )
        tide_va( 1,:) = tide_va( 1,:) + f*va_amp( 1,:,i)*cos( (this_con%speed*tick + uv - va_pha( 1,:,i))*DEG2RAD )
        tide_el(2:imm1, 1) = tide_el(2:imm1, 1)                   &
                           + f*el_amp(2:imm1, 1,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - el_pha(2:imm1, 1,i) )*DEG2RAD )
        tide_ua(2:imm1, 1) = tide_ua(2:imm1, 1)                   &
                           + f*ua_amp(2:imm1, 1,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - ua_pha(2:imm1, 1,i) )*DEG2RAD )
        tide_va(2:imm1, 2) = tide_va(2:imm1, 2)                   &
                           + f*va_amp(2:imm1, 2,i)                &
                            *cos( ( this_con%speed*tick + uv      &
                                  - va_pha(2:imm1, 2,i) )*DEG2RAD )
      end do
      !tide_el = tide_el*fsm!tide_mask
      !tide_ua = tide_ua*fsm!tide_mask
      !tide_va = tide_va*fsm!tide_mask

      ! where( tide_mask(2:im,:) == 1. .and. tide_mask(1:imm1,:) == 0. )
      ! tide_ua(2:im,:) = tide_ua(2:im,:)/2.
        ! tide_ua(2:im,:) = tide_ua(2:im,:)  &
                        ! + (.5*(elf(1:imm1,:)+elf(2:im,:))-el(2:im,:))*sqrt(grav/d(2:im,:))
      ! end where
      ! where( tide_mask(2:im,:) == 0. .and. tide_mask(1:imm1,:) == 1. )
      ! tide_ua(1:imm1,:) = tide_ua(1:imm1,:)/2.
        ! tide_ua(2:im,:) = tide_ua(2:im,:)  &
                        ! - (el(2:im,:)-tide_el(1:imm1,:))*sqrt(grav/d(1:imm1,:))
      ! end where
      ! call exchange2d_mpi(tide_ua,im,jm)

      ! where( tide_mask(:,2:jm) == 1. .and. tide_mask(:,1:jm) == 0. )
        ! tide_va(:,2:jm) = tide_va(:,2:jm)  &
                        ! + (el(:,1:jmm1)-tide_el(:,2:jm))*sqrt(grav/d(:,2:jm))
      ! end where
      ! where( tide_mask(:,2:jm) == 0. .and. tide_mask(:,1:jmm1) == 1. )
        ! tide_va(:,2:jm) = tide_va(:,2:jm)  &
                        ! - (el(:,2:jm)-tide_el(:,1:jmm1))*sqrt(grav/d(:,1:jmm1))
      ! end where
      ! call exchange2d_mpi(tide_va,im,jm)


    end ! subroutine step
!
!______________________________________________________________________
!
    type (T_constituent) function specify_constituent( name, speed, wavenumber, V_arg, u_arg, f_arg )
!----------------------------------------------------------------------
!  Defines tidal constituent
!
! Input:
!  name       - name of constituent
!  speed      - angular speed
!  wavenumber - tidal wavenumber !TODO: Not sure about this word. Maybe tidal constituent subscript?
!  V_arg      - constants of equilibrium argument (V + u): [ a, p, s, h ]
!  u_arg      - coefficients for (V + u) argument
!
!  The corresponding coefficients are derived from u formulation for
! each constituent as in Table II of "Руководство по обработке и
! предсказанию приливов", 1941 (in Russian).
!
!  In this application "u" and "f" arguments have forms of:
!               c1*sin(N) + c2*sin(2*N) + c3*sin(3*N), and
!          c0 + c1*cos(N) + c2*cos(2*N) + c3*cos(3*N) respectively.
!
!  E.g. for M2:
!
!   V = 2t+2h-2s (and V_0 = 2h_0-2s_0),
! so a (in degrees) = 0, p = 0, s = -2 and h = 2.
! Coefficient for t siginifies wavenumber;
!
!   u = 2ξ-2ν,
! where ξ = 11.87*sin(N) - 1.34*sin(2*N) + 0.19*sin(3*N),
! and   ν = 12.94*sin(N) - 1.34*sin(2*N) + 0.19*sin(3*N),
! hence c1 = -2.14, c2 = 0, c3 = 0.
!______________________________________________________________________
!
      implicit none

      character(*), intent(in) :: name
      real(rk)    , intent(in) :: speed
      integer     , intent(in) :: wavenumber
      real        , dimension(4)  &
                  , intent(in) :: V_arg, f_arg
      real        , dimension(3)  &
                  , intent(in) :: u_arg


      specify_constituent % name  = name
      specify_constituent % speed = max( speed, 1.e-4_rk )
      specify_constituent % freq  = 360. / specify_constituent % speed
      specify_constituent % k     = wavenumber

      specify_constituent % V     = real( V_arg, rk )
      specify_constituent % u     = real( u_arg, rk )
      specify_constituent % f     = real( f_arg, rk )


    end ! function
!
!______________________________________________________________________
!
    type (T_constituent) function get_constituent( name )

      implicit none

      character(*), intent(in) :: name

      integer i


      get_constituent = specify_constituent( "n/a", 0._rk, 0, [0.,0.,0.,0.], [0.,0.,0.], [1.,0.,0.,0.,1.] )

      do i =  1, CONS_OVERALL
        if ( name == trim( con(i)%name ) ) then
          get_constituent = con(i)
          return
        end if
      end do


    end ! function
!
!______________________________________________________________________
!
    function get_astronomy( cur_date )

      use glob_const , only: DEG2RAD
      use module_time

      implicit none

      real(rk), dimension(4) :: get_astronomy

      type(date), intent(in) :: cur_date

      integer, dimension(12), parameter ::                             &
               monthly_day_of_year = [   0,  31,  59,  90, 120, 151    &
                                     , 181, 212, 243, 273, 304, 334 ]
      integer  day_of_year   , year
      real(rk) time_remainder


      day_of_year    = seconds_of_year( cur_date )/86400 + 1
      time_remainder = seconds_of_year( cur_date )/86400.+ 1. - day_of_year

      year           = cur_date % year
      time_remainder = time_remainder + real(day_of_year+(year-1901)/4,rk)

      get_astronomy(1) = 277.025_rk + 129.3848 _rk*real(year-1900,rk)  &
                                    +  13.1764 _rk*time_remainder
      get_astronomy(2) = 280.19 _rk -    .23872_rk*real(year-1900,rk)  &
                                    +    .98565_rk*time_remainder
      get_astronomy(3) = 334.385_rk +  40.66249_rk*real(year-1900,rk)  &
                                    +    .1114 _rk*time_remainder
      get_astronomy(4) = 259.157_rk -  19.32818_rk*real(year-1900,rk)  &
                                    -    .05295_rk*time_remainder
      get_astronomy = modulo( get_astronomy, 360._rk )
      get_astronomy(4) = get_astronomy(4)*DEG2RAD


    end !
!
!______________________________________________________________________
!
    subroutine read_all( execute, n, year, record )

      implicit none

      logical              , intent(in) :: execute
      integer              , intent(in) :: n
      integer, dimension(3), intent(in) :: record, year

      integer            ncid
      character(len=128) desc


      if ( .not. execute ) return


    end ! subroutine read_all
!
!______________________________________________________________________
!
    subroutine define_tide_init( file_id, dimensions, type )
!----------------------------------------------------------------------
!  Defines tide module variables for initial output.
!______________________________________________________________________
!
      use io, only: dim_define, var_define

      implicit none

      integer              , intent(in) :: file_id   , type
      integer, dimension(:), intent(in) :: dimensions

      integer varid, dimid


      if ( DISABLED ) return

      dimid = dim_define( file_id, 'ncons', int(ncons) )
      varid = var_define( file_id, 'tide_u_p', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal u-current phase'                  &
                        , 'deg'                                    &
                        , 0, 0., 'east_u north_u' )
      varid = var_define( file_id, 'tide_v_p', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal v-current phase'                  &
                        , 'deg'                                    &
                        , 0, 0., 'east_v north_v' )
      varid = var_define( file_id, 'tide_e_p', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal elevation phase'                  &
                        , 'deg'                                    &
                        , 0, 0., 'east_e north_e' )

      varid = var_define( file_id, 'tide_u_a', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal u-current amplitude'              &
                        , 'deg'                                    &
                        , 0, 0., 'east_u north_u' )
      varid = var_define( file_id, 'tide_v_a', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal v-current amlitude'               &
                        , 'deg'                                    &
                        , 0, 0., 'east_v north_v' )
      varid = var_define( file_id, 'tide_e_a', type                &
                        , [ dimensions(1), dimensions(2), dimid ]  &
                        , 'tidal elevation amplitude'              &
                        , 'deg'                                    &
                        , 0, 0., 'east_e north_e' )


    end ! subroutine
!
!______________________________________________________________________
!
    subroutine define_tide( file_id, dimensions, type )
!----------------------------------------------------------------------
!  Defines tide module variables for output.
!______________________________________________________________________
!
      use io, only: var_define

      implicit none

      integer              , intent(in) :: file_id   , type
      integer, dimension(:), intent(in) :: dimensions

      integer varid


      if ( DISABLED ) return

      varid = var_define( file_id, 'tide_ua', type  &
                        , dimensions                &
                        , 'tide u-velocity'         &
                        , 'm/s'                     &
                        , 0, 0., 'east_u north_u' )
      varid = var_define( file_id, 'tide_va', type  &
                        , dimensions                &
                        , 'tide v-velocity'         &
                        , 'm/s'                     &
                        , 0, 0., 'east_v north_v' )
      varid = var_define( file_id, 'tide_el', type  &
                        , dimensions                &
                        , 'tide elevation'          &
                        , 'm/s'                     &
                        , 0, 0., 'east_e north_e' )


    end ! subroutine
!
!______________________________________________________________________
!
    subroutine write_tide_init( file_id, start, stride )
!----------------------------------------------------------------------
!  Writes tide module variables.
!______________________________________________________________________
!
      use io , only: var_write
      use mpi, only: MPI_OFFSET_KIND

      implicit none

      integer                               , intent(in) :: file_id
      integer(MPI_OFFSET_KIND), dimension(:), intent(in) :: start, stride


      if ( DISABLED ) return

      call var_write( file_id, "tide_u_p", ua_pha                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )
      call var_write( file_id, "tide_v_p", va_pha                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )
      call var_write( file_id, "tide_e_p", el_pha                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )

      call var_write( file_id, "tide_u_a", ua_amp                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )
      call var_write( file_id, "tide_v_a", va_amp                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )
      call var_write( file_id, "tide_e_a", el_amp                    &
                    , [ start(1) , start(2) , 1_MPI_OFFSET_KIND ]    &
                    , [ stride(1), stride(2), ncons             ] )


    end ! subroutine
!
!______________________________________________________________________
!
    subroutine write_tide( file_id, start, stride )
!----------------------------------------------------------------------
!  Writes tide module variables.
!______________________________________________________________________
!
      use io , only: var_write
      use mpi, only: MPI_OFFSET_KIND

      implicit none

      integer                               , intent(in) :: file_id
      integer(MPI_OFFSET_KIND), dimension(:), intent(in) :: start, stride


      if ( DISABLED ) return

      call var_write( file_id, "tide_ua", tide_ua, start, stride )
      call var_write( file_id, "tide_va", tide_va, start, stride )
      call var_write( file_id, "tide_el", tide_el, start, stride )


    end ! subroutine
!
!______________________________________________________________________
!
    pure character(len=256) function get_filename( path, year )
!----------------------------------------------------------------------
!  Costructs filename string in `<path>YYYY<FORMAT_EXT>` format.
!______________________________________________________________________
!
      implicit none

      character(len=*), intent(in) :: path
      integer         , intent(in) :: year


      write( get_filename, '( a, i4.4, a )' ) trim(path)      &
                                            , year            &
                                            , trim(FORMAT_EXT)

    end ! function get_filemname
!
!______________________________________________________________________
!
! INTERPOLATION ROUTINES (TODO: move to separate module)
!______________________________________________________________________
!
    subroutine interpolate_2D( src_x_n, src_y_n, src_x, src_y, src_raw   &
                             , dst_x_n, dst_y_n, dst_x, dst_y, dst, flag )
!----------------------------------------------------------------------
!  Linearly interpolates 2D fields.
!______________________________________________________________________
!
      use glob_domain!, only: my_task

      implicit none

      character                 , intent(in   ) :: flag
      integer(8)                , intent(in   ) :: src_x_n, src_y_n
      integer                   , intent(in   ) :: dst_x_n, dst_y_n
      real(rk), dimension(src_x_n,src_y_n)                            &
                                , intent(in   ) :: src_raw, src_x, src_y
      real(rk), dimension(dst_x_n,dst_y_n)                            &
                                , intent(in   ) ::      dst_x, dst_y
      real(rk), dimension(dst_x_n,dst_y_n)                            &
                                , intent(  out) :: dst

      integer                  i, j, num
      integer, dimension(2) :: x_bl, x_tl, x_tr, x_br  &
                             , y_bl, y_tl, y_tr, y_br, pos
      real(rk)                 wgt, val
      real(rk), dimension(src_x_n,src_y_n)      :: src , src_b
!      integer , dimension(src_x_n,src_y_n)      :: mask, mask_b


      src_b = src_raw
      num = count( src_b == 0._rk )

      do while ( num > 0 )

        src = src_b

        do j = 1, src_y_n
          do i = 1, src_x_n

            if ( src(i,j) == 0._rk ) then

              val = 0.
              wgt = 0.

              if ( j /= 1 ) then
                if ( src(i,j-1) /= 0._rk ) then
                  val = val + src(i,j-1)
                  wgt = wgt + 1.
                end if
              end if
              if ( j /= 1 .and. i /= src_x_n ) then
                if ( src(i+1,j-1) /= 0._rk ) then
                  val = val + .5*src(i+1,j-1)
                  wgt = wgt + .5
                end if
              end if
              if ( i /= src_x_n ) then
                if ( src(i+1,j) /= 0._rk ) then
                  val = val + src(i+1,j)
                  wgt = wgt + 1.
                end if
              end if
              if ( j /= src_y_n .and. i /= src_x_n ) then
                if ( src(i+1,j+1) /= 0._rk ) then
                  val = val + .5*src(i+1,j+1)
                  wgt = wgt + .5
                end if
              end if
              if ( j /= src_y_n ) then
                if ( src(i,j+1) /= 0._rk ) then
                  val = val + src(i,j+1)
                  wgt = wgt + 1.
                end if
              end if
              if ( j /= src_y_n .and. i /= 1 ) then
                if ( src(i-1,j+1) /= 0._rk ) then
                  val = val + .5*src(i-1,j+1)
                  wgt = wgt + .5
                end if
              end if
              if ( i /= 1 ) then
                if ( src(i-1,j) /= 0._rk ) then
                  val = val + src(i-1,j)
                  wgt = wgt + 1.
                end if
              end if
              if ( j /= 1 .and. i /= 1 ) then
                if ( src(i-1,j-1) /= 0._rk ) then
                  val = val + .5*src(i-1,j-1)
                  wgt = wgt + .5
                end if
              end if
              if ( wgt > 0. ) then
                num = num-1
                src_b(i,j) = val / wgt
              end if
            end if
          end do
        end do
        !num = count( src_b == 0._rk )
      end do
      do j = 1, dst_y_n
        do i = 1, dst_x_n
          pos = get_bottom_left( dst_x(i,j), dst_y(i,j), src_x, src_y, src_x_n, src_y_n )
          dst(i,j) = bilint( dst_x(i,j), dst_y(i,j)                                    &
                           , src_x(pos(1)  ,pos(2)  ), src_x(pos(1)  ,pos(2)+1)    &
                           , src_x(pos(1)+1,pos(2)+1), src_x(pos(1)+1,pos(2)  )    &
                           , src_y(pos(1)  ,pos(2)  ), src_y(pos(1)  ,pos(2)+1)    &
                           , src_y(pos(1)+1,pos(2)+1), src_y(pos(1)+1,pos(2)  )    &
                           , src(pos(1)  ,pos(2)  ), src(pos(1)  ,pos(2)+1)    &
                           , src(pos(1)+1,pos(2)+1), src(pos(1)+1,pos(2)  )    &
                           , 0._rk, flag )
          ! if ( i_global(1)<141.488 .and. i_global(im)>141.488 .and. j_global(1)<53.2246 .and. j_global(jm)>53.2246 ) then
          ! if ( dst_x(i,j)>141.48 .and. dst_x(i,j)<141.49 .and. dst_y(i,j)>53.22 .and. dst_y(i,j)<53.23 ) then
            ! print *, "DST: ", dst_y(i,j), dst_x(i,j), "||", i,j, "[",my_task,"]"
            ! print *, "x1:", src_x(pos(1),pos(2))
            ! print *, "x2:", src_x(pos(1),pos(2)+1)
            ! print *, "x3:", src_x(pos(1)+1,pos(2)+1)
            ! print *, "x4:", src_x(pos(1)+1,pos(2))
            ! print *, "y1:", src_y(pos(1),pos(2))
            ! print *, "y2:", src_y(pos(1),pos(2)+1)
            ! print *, "y3:", src_y(pos(1)+1,pos(2)+1)
            ! print *, "y4:", src_y(pos(1)+1,pos(2))
            ! print *, "f1:", src(pos(1),pos(2))
            ! print *, "f2:", src(pos(1),pos(2)+1)
            ! print *, "f3:", src(pos(1)+1,pos(2)+1)
            ! print *, "f4:", src(pos(1)+1,pos(2))
            ! print *, "res: ", dst(i,j)
            ! call finalize_mpi
            ! stop
          ! end if
          ! end if
        end do
      end do
      !call finalize_mpi
      ! stop


      return

    end ! subroutine interpolate_to_grid
!
!______________________________________________________________________
!
    function get_bottom_left( tgt_x, tgt_y, x, y, x_n, y_n ) result( pos )
!----------------------------------------------------------------------
!  Finds indicies of the bottom left corner cell in x and y
! curvilinear coordinates.
!  NOTE: geographical degrees are not converted to meters.
!______________________________________________________________________
!
      use glob_domain, only: my_task
      implicit none

      integer(8)                   , intent(in) :: x_n, y_n
      real(rk)                     , intent(in) :: tgt_x, tgt_y
      real(rk), dimension(x_n, y_n), intent(in) :: x, y

      integer, dimension(2) :: pos

!      integer                          clk0, clk1, clkrate
      real(rk), dimension(x_n, y_n) :: dist


!      call system_clock(clk0,clkrate)
!      clk0 = clk0/clkrate

      dist = sqrt( (tgt_x-x)**2 + (tgt_y-y)**2 )

!      call system_clock(clk1,clkrate)
!      clk1 = clk1/clkrate
!      print *, "<<<", ( clk1 - clk0 )
      
      pos = minloc(dist, x < tgt_x .and. y < tgt_y )
      
!      write(40+my_task,*) tgt_x,";",tgt_y,";",x(pos(1),pos(2)),";",y(pos(1),pos(2)),":",minval(dist),";",dist(pos(1),pos(2))


    end ! finction
!
!______________________________________________________________________
!
    real(rk) function bilint(  x,  y, x1, x2, x3, x4, y1, y2, y3, y4  &
                            , f1, f2, f3, f4, fill, flag )
!----------------------------------------------------------------------
!  Performs bilinear interpolation.
!______________________________________________________________________
!
      use glob_const, only: DEG2RAD, RAD2DEG

      implicit none

      real(rk), intent(in) :: x,y, x1,x2,x3,x4, y1,y2,y3,y4, f1,f2,f3,f4, fill
      character, intent(in):: flag

      integer  count
      real(rk) t,s, a1,a2,a3,a4, b1,b2,b3,b4, A,B,C, tot_w, w1,w2,w3,w4
      real(rk), parameter :: SMALL = 1.e-3_rk


      a1 =   x1 - x2 + x3 - x4
      a2 = - x1           + x4
      a3 = - x1 + x2
      a4 =   x1                - x

      b1 =   y1 - y2 + y3 - y4
      b2 = - y1           + y4
      b3 = - y1 + y2
      b4 =   y1                - y
      
!      print *, "aa: ", a1, a2, a3, a4
!      print *, "bb: ", b1, b2, b3, b4

      A = a3*b1         - a1*b3
      B = a4*b1 + a3*b2 - a2*b3 - a1*b4
      C =         a4*b2         - a2*b4
      
!      print *, "tABC: ", A, B, C

!      if ( abs( A ) > SMALL ) then
      if ( abs( A*C ) > .002*B*B ) then
        t = ( - B - sqrt( B*B - 4.*A*C ) ) / ( 2.*A )
      else
        t = C / abs( B )
      end if

      A = a2*b1 - a1*b2
      B = a4*b1 - a3*b2 + a2*b3 - a1*b4
      C =                 a4*b3 - a3*b4
      
!      print *, "sABC: ", A, B, C

!      if ( abs( A ) > SMALL ) then
      if ( abs( A*C ) > .002*B*B ) then
        s = ( - B + sqrt( B*B - 4.*A*C ) ) / ( 2.*A )
      else
        s = - C / abs( B )
      end if

      w1 = ( 1.-t ) * ( 1.-s )
      w2 =      t   * ( 1.-s )
      w3 =      t   *      s
      w4 = ( 1.-t ) *      s
!      print *, "WW: ", w1, w2, w3, w4
!      print *, "ts: ", t, s

     count = 0

     if ( w1 < 0. .and. w3 < 0. ) then
       count = 1
       w1 = 0.
       w2 = 1.
       w3 = 0.
       w4 = 0.
     end if
     if ( w3 < 0. .and. w4 < 0. ) then
       count = 2
       w3 = 0.
       w4 = 0.
     end if
     if ( w2 < 0. .and. w4 < 0. ) then
       count = 3
       w1 = 1.
       w2 = 0.
       w3 = 0.
       w4 = 0.
     end if
     if ( w2 < 0. .and. w3 < 0. ) then
       count = 4
       w2 = 0.
       w3 = 0.
     end if
     if ( w1 < 0. .and. w3 < 0. ) then
       count = 5
       w1 = 1.
       w2 = 0.
       w3 = 0.
       w4 = 1.
     end if
     if ( w1 < 0. .and. w2 < 0. ) then
       count = 6
       w1 = 0.
       w2 = 0.
     end if
     if ( w2 < 0. .and. w4 < 0. ) then
       count = 7
       w1 = 1.
       w2 = 0.
       w3 = 1.
       w4 = 0.
     end if
     if ( w1 < 0. .and. w4 < 0. ) then
       count = 8
       w1 = 0.
       w4 = 0.
     end if

!       if ( f1 == fill ) w1 = 0.
!       if ( f2 == fill ) w2 = 0.
!       if ( f3 == fill ) w3 = 0.
!       if ( f4 == fill ) w4 = 0.
!      print *, "FILL: ", (f1==fill), (f2==fill), (f3==fill), (f4==fill)

      ! tot_w = w1 + w2 + w3 + w4
      ! if ( tot_w > 0. .and. tot_w /= 1. ) then
        ! w1 = w1 / tot_w
        ! w2 = w2 / tot_w
        ! w3 = w3 / tot_w
        ! w4 = w4 / tot_w
      ! end if

      if ( flag == 'P' ) then
        ! bilint = 0.
        ! a2 = f1-f2
        ! a3 = f1-f3
        ! a4 = f1-f4
        ! b1 = f1
        ! b2 = f2
        ! b3 = f3
        ! b4 = f4
        ! if ( a2 <= -180. ) then
          ! b2 = f2-360.
        ! elseif ( a2 >= 180. ) then
          ! b2 = f2+360.
        ! end if
        ! if ( a3 <= -180. ) then
          ! b3 = f3-360.
        ! elseif ( a3 >= 180. ) then
          ! b3 = f3+360.
        ! end if
        ! if ( a4 <= -180. ) then
          ! b4 = f4-360.
        ! elseif ( a4 >= 180. ) then
          ! b4 = f4+360.
        ! end if
        b1 = cos(f1*DEG2RAD)
        b2 = cos(f2*DEG2RAD)
        b3 = cos(f3*DEG2RAD)
        b4 = cos(f4*DEG2RAD)
        bilint = w1*b1 + w2*b2 + w3*b3 + w4*b4
        b1 = sin(f1*DEG2RAD)
        b2 = sin(f2*DEG2RAD)
        b3 = sin(f3*DEG2RAD)
        b4 = sin(f4*DEG2RAD)
        bilint = atan2( w1*b1+w2*b2+w3*b3+w4*b4, bilint )*RAD2DEG
        ! if ( bilint > 360. ) then
          ! bilint = bilint - 360.
        ! elseif ( bilint < 0. ) then
          ! bilint = bilint + 360.
        ! end if
      else
        bilint = w1*f1 + w2*f2 + w3*f3 + w4*f4
      end if
!      print *, "INT:", f1,f2,f3,f4, "[", w1,w2,w3,w4 ,"] =>", bilint
!      stop

    end !
!
!______________________________________________________________________
!
    real(rk) function nbrint(  x,  y, x1, x2, x3, x4, y1, y2, y3, y4  &
                            , f1, f2, f3, f4, fill   )
!----------------------------------------------------------------------
!  Performs bilinear interpolation.
!______________________________________________________________________
!
      implicit none

      real(rk), intent(in) :: x,y, x1,x2,x3,x4, y1,y2,y3,y4, f1,f2,f3,f4, fill

      real(rk), dimension(4) :: d, f
      integer i


      f(1) = f1
      f(2) = f2
      f(3) = f3
      f(4) = f4

      d(1) = sqrt( (x-x1)**2 + (y-y1)**2 )
      d(2) = sqrt( (x-x2)**2 + (y-y2)**2 )
      d(3) = sqrt( (x-x3)**2 + (y-y3)**2 )
      d(4) = sqrt( (x-x4)**2 + (y-y4)**2 )

      do i = 1, 4
        if ( f(i) == fill ) d(i) = HUGE(0._rk)
      end do

      nbrint = f(minloc(d,1))


    end ! 

end module tide
