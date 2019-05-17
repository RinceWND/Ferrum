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
  real(rk), parameter ::     & !
    Cp     = 1005.             ! Specific heat capacity of air
    
  type(date) tide_offset

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED ! Is set according to the external flag `use_air`.

  real   , private :: a         ! time-interpolation factor

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=PATH_LEN)  & ! Full paths (with filenames) to:
    tide_el_path           & !  atmospheric parameters file
  , tide_uv_path             !  surface fluxes file

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

  integer              &
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

  type T_constituent
    character(6) name
    real(rk)     freq
    integer      k
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
        tide_el_path, tide_uv_path

      namelist/tide_vars/                                   &
        el_amp_name, el_pha_name,    lat_name,  lat_u_name  &
      ,  lat_v_name,    lon_name,  lon_u_name,  lon_v_name  &
      , ua_amp_name, ua_pha_name, va_amp_name, va_pha_name  &
      ,   cons_name


      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_TIDE ) DISABLED = .true.

      if ( DISABLED ) return

      tide_offset = str2date("1992-01-01_00:00:00")

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

! Generate tidal constituents
!  - semidiurnal:
      con( 1) = specify_constituent( "m2", 28.9841042_rk, 2 ) ! Principal lunar semidiurnal
      con( 2) = specify_constituent( "s2", 30.       _rk, 2 ) ! Principal solar semidiurnal
      con( 3) = specify_constituent( "n2", 28.4397295_rk, 2 ) ! Larger lunar elliptic semidiurnal
      con( 4) = specify_constituent( "k2", 30.0821373_rk, 2 ) ! Lunisolar declinational semidiurnal
!  - diurnal:
      con( 5) = specify_constituent( "o1", 13.9430356_rk, 1 ) ! Lunar diurnal
      con( 6) = specify_constituent( "p1", 14.9589314_rk, 1 ) ! Solar diurnal
      con( 7) = specify_constituent( "q1", 13.3986609_rk, 1 ) ! Larger lunar elliptic diurnal
      con( 8) = specify_constituent( "k1", 15.0410686_rk, 1 ) ! Lunisolar declinational diurnal
!  - shallow water:
      con( 9) = specify_constituent( "m4" ,  57.9682084_rk, 4 ) ! Shallow water overtides of principal lunar constituent, quarter-diurnal
      con(10) = specify_constituent( "m6" ,  86.9523127_rk, 6 ) ! Shallow water overtides of principal lunisolar constituent, sixth-diurnal
      con(11) = specify_constituent( "ms4",  58.9841042_rk, 4 ) ! Shallow water lunisolar quarter-diurnal
      con(12) = specify_constituent( "m8" , 115.93642  _rk, 8 ) ! Shallow water lunar eighth-diurnal
      con(13) = specify_constituent( "s8" , 120.       _rk, 8 ) ! Shallow water solar eighth-diurnal
      con(14) = specify_constituent( "k3" ,  43.31412  _rk, 3 ) ! Shallow water lunisolar declinational third-diurnal
      con(15) = specify_constituent( "m3" ,  43.4761563_rk, 3 ) ! Shallow water lunar third-diurnal
!  - low-frequency:
      con(16) = specify_constituent( "msf", 1.0158958_rk, 2 ) ! Lunisolar synodic fortnightly
      con(17) = specify_constituent( "mf" , 1.0980331_rk, 2 ) ! Lunar fortnightly
      con(18) = specify_constituent( "mm" ,  .5443747_rk, 1 ) ! Lunar monthly
      con(19) = specify_constituent( "msm",  .4716   _rk, 1 ) ! Solar monthly
      con(20) = specify_constituent( "ssa",  .0821373_rk, 2 ) ! Solar semiannual
      con(21) = specify_constituent( "sa" ,  .0410686_rk, 1 ) ! Solar annual
!  - superimposed:
      con(22) = specify_constituent( "mu2" ,  27.9682084_rk, 1 ) ! Variational (M2,S2)
      con(23) = specify_constituent( "j1"  ,  15.5854433_rk, 1 ) ! Smaller lunar elliptic diurnal
      con(24) = specify_constituent( "nu2" ,  28.5125831_rk, 1 ) ! Larger lunar evectional
      con(25) = specify_constituent( "oo1" ,  16.1391017_rk, 1 ) ! Lunar diurnal
      con(26) = specify_constituent( "2n2" ,  27.8953548_rk, 1 ) ! Lunar elliptical semidiurnal second-order
      con(27) = specify_constituent( "m1"  ,  14.4920521_rk, 1 ) ! Smaller lunar elliptic diurnal ! or 14.496694
      con(28) = specify_constituent( "t2"  ,  29.9589333_rk, 1 ) ! Larger solar elliptic
      con(29) = specify_constituent( "l2"  ,  29.5284789_rk, 1 ) ! Smaller lunar elliptic semidiurnal
      con(30) = specify_constituent( "mn4" ,  57.4238337_rk, 1 ) ! Shallow water quarter diurnal
      con(31) = specify_constituent( "pi1" ,  14.9178647_rk, 1 ) ! ?
      con(32) = specify_constituent( "2sm2",  31.0158958_rk, 1 ) ! Shallow water semidiurnal
      con(33) = specify_constituent( "phi1",  15.1232059_rk, 1 ) ! ?
      con(34) = specify_constituent( "2ms6",  87.9682084_rk, 1 ) ! ?
      con(35) = specify_constituent( "sn4" ,  58.4397295_rk, 1 ) ! ?
      con(36) = specify_constituent( "ksi1",  14.5695476_rk, 1 ) ! ?
      con(37) = specify_constituent( "mo3" ,  42.9271398_rk, 1 ) ! ?
      con(38) = specify_constituent( "2mn6",  96.407938 _rk, 1 ) ! ?
      con(39) = specify_constituent( "mk3" ,  44.0251729_rk, 1 ) ! Shallow water terdiurnal
      con(40) = specify_constituent( "msn6",  87.4238337_rk, 1 ) ! ?
      con(41) = specify_constituent( "2sm6",  88.9841042_rk, 1 ) ! ?
      con(42) = specify_constituent( "s9"  , 135.       _rk, 1 ) ! Shallowwater overtides of the principal solar, ninth-diurnal
      con(43) = specify_constituent( "s7"  , 105.       _rk, 1 ) ! Shallowwater overtides of the principal solar, seventh-diurnal
      con(44) = specify_constituent( "s6"  ,  90.       _rk, 1 ) ! Shallowwater overtides of the principal solar, sixth-diurnal
      con(45) = specify_constituent( "s5"  ,  75.       _rk, 1 ) ! Shallowwater overtides of the principal solar, fifth-diurnal
      con(46) = specify_constituent( "msk6",  89.066241 _rk, 1 ) ! ?
      con(47) = specify_constituent( "2mk6",  88.0503457_rk, 1 ) ! ?
      con(48) = specify_constituent( "sk4" ,  60.0821373_rk, 1 ) ! ?
      con(49) = specify_constituent( "s4"  ,  60.       _rk, 1 ) ! Shallow water overtides of principal solar constituent
      con(50) = specify_constituent( "mk4" ,  59.0662415_rk, 1 ) ! ?
      con(51) = specify_constituent( "sk3" ,  45.0410686_rk, 1 ) ! ?
      con(52) = specify_constituent( "so3" ,  43.943036 _rk, 1 ) ! ?
      con(53) = specify_constituent( "kj2" ,  30.626512 _rk, 1 ) ! ?
      con(54) = specify_constituent( "msn2",  30.5443747_rk, 1 ) ! ?
      con(55) = specify_constituent( "mks2",  29.0662415_rk, 1 ) ! ?
      con(56) = specify_constituent( "op2" ,  28.9019669_rk, 1 ) ! ?
      con(57) = specify_constituent( "mns2",  27.4238337_rk, 1 ) ! ?
      con(58) = specify_constituent( "oq2" ,  27.3416964_rk, 1 ) ! ?
      con(59) = specify_constituent( "so1" ,  16.0569644_rk, 1 ) ! ?
      con(60) = specify_constituent( "mp1" ,  14.0251729_rk, 1 ) ! ?
      con(61) = specify_constituent( "s3"  ,  45.       _rk, 1 ) ! ?
      con(62) = specify_constituent( "r2"  ,  30.0410667_rk, 1 ) ! ?

      con(63) = specify_constituent( "lamda2", 29.4556253_rk, 1 ) ! Smaller lunar evectional
      con(64) = specify_constituent( "theta1", 15.5125897_rk, 1 ) ! ?
      con(65) = specify_constituent( "psi1"  , 15.0821353_rk, 1 ) ! ?
      con(66) = specify_constituent( "s1"    , 15.       _rk, 1 ) ! Solar diurnal
      con(67) = specify_constituent( "rho1"  , 13.4715145_rk, 1 ) ! Larger lunar evectional diurnal
      con(68) = specify_constituent( "sigma1", 12.9271398_rk, 1 ) ! ?
      con(69) = specify_constituent( "2q1"   , 12.8542862_rk, 1 ) ! ?
      con(70) = specify_constituent( "no1"   , 14.4966939_rk, 1 ) ! ?
      con(71) = get_constituent( "ksi1" )
      con(71) % name = "lp1"
      con(72) = get_constituent( "pi1" )
      con(72) % name = "tk1"
      con(73) = get_constituent( "psi1" )
      con(73) % name = "rp1"
      con(74) = get_constituent( "phi1" )
      con(74) % name = "kp1"
      con(75) = specify_constituent( "kq1", 16.6834764_rk, 1 ) ! ?
      con(76) = get_constituent( "mu2" )
      con(76) % name = "2ms2"
      con(77) = specify_constituent( "mq3", 42.3827651_rk, 1 ) ! ?
      con(78) = get_constituent( "mo3" )
      con(78) % name = "2mk3"
      con(79) = specify_constituent( "sp3"  ,  44.9589314_rk, 1 ) ! ?
      con(80) = specify_constituent( "2mns4",  56.407938 _rk, 1 ) ! ?
      con(81) = specify_constituent( "3mk4" ,  56.8701754_rk, 1 ) ! ?
      con(82) = specify_constituent( "3ms4" ,  56.9523127_rk, 1 ) ! ?
      con(83) = specify_constituent( "2msk4",  57.8660711_rk, 1 ) ! ?
      con(84) = specify_constituent( "2mks4",  58.0503457_rk, 1 ) ! ?
      con(85) = specify_constituent( "3mn4" ,  58.5125831_rk, 1 ) ! ?
      con(86) = specify_constituent( "2smk4",  58.9019669_rk, 1 ) ! ?
      con(87) = specify_constituent( "2snm4",  59.4556253_rk, 1 ) ! ?
      con(88) = specify_constituent( "2msn4",  59.5284789_rk, 1 ) ! ?
      con(89) = specify_constituent( "2smn4",  60.5443747_rk, 1 ) ! ?
      con(90) = specify_constituent( "3sm4" ,  61.0158958_rk, 1 ) ! ?
      con(91) = specify_constituent( "2skm4",  61.0980331_rk, 1 ) ! ?
      con(92) = specify_constituent( "mnk6" ,  87.505971 _rk, 1 ) ! ?
      con(93) = specify_constituent( "2mn8" , 114.8476674_rk, 1 ) ! ?
      con(94) = specify_constituent( "3mn8" , 115.3920422_rk, 1 ) ! ?
      con(95) = specify_constituent( "2msn8", 116.407938 _rk, 1 ) ! ?
      con(96) = specify_constituent( "2mnk8", 116.4900753_rk, 1 ) ! ?
      con(97) = specify_constituent( "3ms8" , 116.9523127_rk, 1 ) ! ?
      con(98) = specify_constituent( "3mk8" , 117.0344499_rk, 1 ) ! ?
      con(99) = specify_constituent( "2smn8", 117.4238337_rk, 1 ) ! ?
      con(100)= specify_constituent( "msnk8", 117.505971 _rk, 1 ) ! ?
      con(101)= specify_constituent( "2ms8" , 117.9682084_rk, 1 ) ! ?
      con(102)= specify_constituent( "2msk8", 118.0503457_rk, 1 ) ! ?
      con(103)= specify_constituent( "3sm8" , 118.9841042_rk, 1 ) ! ?
      con(104)= specify_constituent( "2smk8", 119.0662415_rk, 1 ) ! ?


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
      use glob_grid  , only:  east_e,  east_u,  east_v  &
                           , north_e, north_u, north_v
      use mpi        , only: MPI_INFO_NULL
      use pnetcdf

      implicit none

      integer                                   ncid, varid, ndims
      integer(8)                                nx, ny
      integer  , allocatable, dimension(:)   :: cons_dimids, cons_dimlens
      real(rk) , allocatable, dimension(:,:) :: src_lat, src_lat_u, src_lat_v  &
                                              , src_lon, src_lon_u, src_lon_v  &
                                              , tmp


! Read files
      call check( nf90mpi_open( POM_COMM, tide_el_path, NF_NOWRITE   &
                              , MPI_INFO_NULL, ncid )                &
                , "nf_open: "//tide_el_path )

      call check( nf90mpi_inq_varid( ncid, el_amp_name, varid )  &
                , "nf_inq_varid for dimension extraction: "//el_amp_name )
      call check( nf90mpi_inquire_variable( ncid, varid, ndims = ndims ) &
                , "nf_inq_var_ndims for dimension extraction: "//el_amp_name )
      allocate( cons_dimids(ndims), cons_dimlens(ndims) )
      call check( nf90mpi_inquire_variable( ncid, varid, dimids = cons_dimids ) &
                , "nf_inq_var_dimids for dimension extraction: "//cons_name )
      call check( nf90mpi_inquire_dimension( ncid, cons_dimids(3), len = ncons ) &
                , "nf_inq_dim_len: number of constituents" )
      allocate( constituent(ncons) )
      call check( nf90mpi_inquire_dimension( ncid, cons_dimids(2), len = nx ) &
                , "nf_inq_dim_len: number of longitudes" )
      call check( nf90mpi_inquire_dimension( ncid, cons_dimids(1), len = ny ) &
                , "nf_inq_dim_len: number of latitudes" )

      allocate( src_lat(nx,ny), src_lat_u(nx,ny), src_lat_v(nx,ny)    &
              , src_lon(nx,ny), src_lon_u(nx,ny), src_lon_v(nx,ny)    &
              , tmp(ny,nx) )

      call check( nf90mpi_inq_varid( ncid, lat_name, varid )  &
                , "nf_inq_varid: "//lat_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lat_name )
      src_lat = transpose( tmp )
      call check( nf90mpi_inq_varid( ncid, lon_name, varid )  &
                , "nf_inq_varid: "//lon_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lon_name )
      src_lon = transpose( tmp )

      call check( nf90mpi_close( ncid ), "nf_close"//tide_el_path )


      call check( nf90mpi_open( POM_COMM, tide_uv_path, NF_NOWRITE   &
                              , MPI_INFO_NULL, ncid )                &
                , "nf_open: "//tide_uv_path )

      call check( nf90mpi_inq_varid( ncid, lat_u_name, varid )  &
                , "nf_inq_varid: "//lat_u_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lat_u_name )
      src_lat_u = transpose( tmp )
      call check( nf90mpi_inq_varid( ncid, lon_u_name, varid )  &
                , "nf_inq_varid: "//lon_u_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lon_u_name )
      src_lon_u = transpose( tmp )

      call check( nf90mpi_inq_varid( ncid, lat_v_name, varid )  &
                , "nf_inq_varid: "//lat_v_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lat_v_name )
      src_lat_v = transpose( tmp )
      call check( nf90mpi_inq_varid( ncid, lon_v_name, varid )  &
                , "nf_inq_varid: "//lon_v_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp )  &
                , "nf_get_var: "//lon_v_name )
      src_lon_v = transpose( tmp )

! Crop domain and store these extents for further use
      ie_lft = maxloc(src_lon(:,1), 1, src_lon(:,1) < minval(east_e))
      ie_rht = minloc(src_lon(:,1), 1, src_lon(:,1) > maxval(east_e))
      iu_lft = maxloc(src_lon_u(:,1), 1, src_lon_u(:,1) < minval(east_u))
      iu_rht = minloc(src_lon_u(:,1), 1, src_lon_u(:,1) > maxval(east_u))
      iv_lft = maxloc(src_lon_v(:,1), 1, src_lon_v(:,1) < minval(east_v))
      iv_rht = minloc(src_lon_v(:,1), 1, src_lon_v(:,1) > maxval(east_v))
      if ( src_lat(1,1) < src_lat(1,ny) ) then
        je_bot = maxloc(src_lat(1,:), 1, src_lat(1,:) < minval(north_e))
        je_top = minloc(src_lat(1,:), 1, src_lat(1,:) > maxval(north_e))
        ju_bot = maxloc(src_lat_u(1,:), 1, src_lat_u(1,:) < minval(north_u))
        ju_top = minloc(src_lat_u(1,:), 1, src_lat_u(1,:) > maxval(north_u))
        jv_bot = maxloc(src_lat_v(1,:), 1, src_lat_v(1,:) < minval(north_v))
        jv_top = minloc(src_lat_v(1,:), 1, src_lat_v(1,:) > maxval(north_v))
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

! Allocate tidal arrays
      allocate(              &
        el_amp(im,jm,ncons)  &
      , el_pha(im,jm,ncons)  &
      , ua_amp(im,jm,ncons)  &
      , ua_pha(im,jm,ncons)  &
      , va_amp(im,jm,ncons)  &
      , va_pha(im,jm,ncons)  &
      , tide_el(im,jm)       &
      , tide_el_b(im,jm)     &
      , tide_mask(im,jm)     &
      , tide_ua(im,jm)       &
      , tide_ua_b(im,jm)     &
      , tide_va(im,jm)       &
      , tide_va_b(im,jm)     &
       )

       el_amp = 0.
       el_pha = 0.
       ua_amp = 0.
       ua_pha = 0.
       va_amp = 0.
       va_pha = 0.
       tide_el = 0.
       tide_ua = 0.
       tide_va = 0.


    end ! subroutine allocate_mod
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields before experiment's start.
!______________________________________________________________________
!
      use glob_domain, only: im, jm, my_task, POM_COMM
      use glob_grid  , only: fsm, east_e, east_u, east_v, rot  &
                           , north_e, north_u, north_v
      use air, only:e_atmos
      use module_time
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      type(date), intent(in) :: d_in

      integer                                    ncid, varid, dimlen, i,j,k
      integer(MPI_OFFSET_KIND)                :: start(3), edge(3)
      real(rk), allocatable, dimension(:,:,:) :: tmp, var


! Quit if the module is not used.
      if ( DISABLED ) return

! Expand free surface mask
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
! Read rho-values (and constituent names)
      call check( nf90mpi_open( POM_COMM, tide_el_path, NF_NOWRITE   &
                              , MPI_INFO_NULL, ncid )                &
                , "nf_open: "//tide_el_path )

      call check( nf90mpi_inq_varid( ncid, cons_name, varid )  &
                , "nf_inq_varid: "//cons_name )
      call check( nf90mpi_get_var_all( ncid, varid, constituent )  &
                , "nf_get_var: "//cons_name )
      print *, my_task, ":", constituent, shape(constituent)

! Read coordinates
      allocate( tmp( je_top-je_bot+1, ie_rht-ie_lft+1, 1 ) )

      start = (/ je_bot         , ie_lft         , 1 /)
      edge  = (/ je_top-je_bot+1, ie_rht-ie_lft+1, 1 /)

      call check( nf90mpi_inq_varid( ncid, lat_name, varid )  &
                , "nf_inq_varid: "//lat_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
                                     , start, edge )  &
                , "nf_get_var: "//lat_name )
      lat = transpose( tmp(:,:,1) )

      call check( nf90mpi_inq_varid( ncid, lon_name, varid )  &
                , "nf_inq_varid: "//lon_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
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
      
      start = (/ je_bot         , ie_lft         , 1          /)
      edge  = (/ je_top-je_bot+1, ie_rht-ie_lft+1, int(ncons) /)

      call check( nf90mpi_inq_varid( ncid, el_amp_name, varid )  &
                , "nf_inq_varid: "//el_amp_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
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

      call check( nf90mpi_inq_varid( ncid, el_pha_name, varid )  &
                , "nf_inq_varid: "//el_pha_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp          &
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

      call check( nf90mpi_close( ncid ), "nf_close"//tide_uv_path )

! Read uv-values
      call check( nf90mpi_open( POM_COMM, tide_uv_path, NF_NOWRITE   &
                              , MPI_INFO_NULL, ncid )                &
                , "nf_open: "//tide_uv_path )

! Read u-coordinates
      deallocate( tmp )
      allocate( tmp( ju_top-ju_bot+1, iu_rht-iu_lft+1, 1 ) )

      start = (/ ju_bot         , iu_lft         , 1 /)
      edge  = (/ ju_top-ju_bot+1, iu_rht-iu_lft+1, 1 /)

      call check( nf90mpi_inq_varid( ncid, lat_u_name, varid )  &
                , "nf_inq_varid: "//lat_u_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
                                     , start, edge )  &
                , "nf_get_var: "//lat_u_name )
      lat_u = transpose( tmp(:,:,1) )

      call check( nf90mpi_inq_varid( ncid, lon_u_name, varid )  &
                , "nf_inq_varid: "//lon_u_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
                                     , start, edge )  &
                , "nf_get_var: "//lon_u_name )
      lon_u = transpose( tmp(:,:,1) )

! Read u-amplitudes and phases
      deallocate( tmp, var )
      allocate( tmp( ju_top-ju_bot+1, iu_rht-iu_lft+1, ncons )    &
              , var( iu_rht-iu_lft+1, ju_top-ju_bot+1, ncons ) )
      
      start = (/ ju_bot         , iu_lft         , 1          /)
      edge  = (/ ju_top-ju_bot+1, iu_rht-iu_lft+1, int(ncons) /)

      call check( nf90mpi_inq_varid( ncid, ua_amp_name, varid )  &
                , "nf_inq_varid: "//ua_amp_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
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

      call check( nf90mpi_inq_varid( ncid, ua_pha_name, varid )  &
                , "nf_inq_varid: "//ua_pha_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp          &
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

      start = (/ jv_bot         , iv_lft         , 1 /)
      edge  = (/ jv_top-jv_bot+1, iv_rht-iv_lft+1, 1 /)

      call check( nf90mpi_inq_varid( ncid, lat_v_name, varid )  &
                , "nf_inq_varid: "//lat_v_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
                                     , start, edge )  &
                , "nf_get_var: "//lat_v_name )
      lat_v = transpose( tmp(:,:,1) )

      call check( nf90mpi_inq_varid( ncid, lon_v_name, varid )  &
                , "nf_inq_varid: "//lon_v_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
                                     , start, edge )  &
                , "nf_get_var: "//lon_v_name )
      lon_v = transpose( tmp(:,:,1) )

! Read v-amplitudes and phases
      deallocate( tmp, var )
      allocate( tmp( jv_top-jv_bot+1, iv_rht-iv_lft+1, ncons )    &
              , var( iv_rht-iv_lft+1, jv_top-jv_bot+1, ncons ) )
      
      start = (/ jv_bot         , iv_lft         , 1          /)
      edge  = (/ jv_top-jv_bot+1, iv_rht-iv_lft+1, int(ncons) /)

      call check( nf90mpi_inq_varid( ncid, va_amp_name, varid )  &
                , "nf_inq_varid: "//va_amp_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp               &
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

      call check( nf90mpi_inq_varid( ncid, va_pha_name, varid )  &
                , "nf_inq_varid: "//va_pha_name )
      call check( nf90mpi_get_var_all( ncid, varid, tmp          &
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
      
      call check( nf90mpi_close( ncid ), "nf_close"//tide_uv_path )

      deallocate( tmp, var )

! Rotate uv-variables to conform to the grid using `tmp` and `var` as u and v
      allocate( tmp( im, jm, ncons )    &
              , var( im, jm, ncons ) )

      do i = 1, int(ncons)
        tmp(:,:,i) = ua_amp(:,:,i)*cos(rot) - va_amp(:,:,i)*sin(rot)
        var(:,:,i) = ua_amp(:,:,i)*sin(rot) + va_amp(:,:,i)*cos(rot)
      end do
      ua_amp = tmp/100.
      va_amp = var/100.

      do i = 1, int(ncons)
        tmp(:,:,i) = ua_pha(:,:,i)*cos(rot) - va_pha(:,:,i)*sin(rot)
        var(:,:,i) = ua_pha(:,:,i)*sin(rot) + va_pha(:,:,i)*cos(rot)
      end do
      where ( tmp > 360. ) tmp = tmp - 360.
      where ( tmp <   0. ) tmp = tmp + 360.
      where ( var > 360. ) var = var - 360.
      where ( var <   0. ) var = var + 360.
      ua_pha = tmp
      va_pha = var

      deallocate( tmp, var )

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
      use glob_const , only: DEG2RAD
!      use glob_domain, only: is_master
      use module_time
      use model_run  , only: dti, iint, sec_of_year
      use seaice     , only: icec

      implicit none

      type(date), intent(in) :: d_in

      integer             i, tick
      type(T_constituent) this_con


! Quit if the module is not used.
      if ( DISABLED ) return

      tick = d_in - tide_offset

      tide_el_b = tide_el
      tide_ua_b = tide_ua
      tide_va_b = tide_va
      tide_el = 0.
      tide_ua = 0.
      tide_va = 0.

      do i = 1, int(ncons)
        this_con = get_constituent( constituent(i) )
        tide_el   = tide_el + el_amp(:,:,i)*cos( this_con%freq*tick + el_pha(:,:,i)*DEG2RAD )
        tide_ua = tide_ua + ua_amp(:,:,i)*cos( this_con%freq*tick + ua_pha(:,:,i)*DEG2RAD )
        tide_va = tide_va + va_amp(:,:,i)*cos( this_con%freq*tick + va_pha(:,:,i)*DEG2RAD )
      end do
      tide_el = tide_el*tide_mask
      tide_ua = tide_ua*tide_mask
      tide_va = tide_va*tide_mask


    end ! subroutine step
!
!______________________________________________________________________
!
    type (T_constituent) function specify_constituent( name, speed, wavenumber )

      implicit none

      character(*), intent(in) :: name
      real(rk)    , intent(in) :: speed
      integer     , intent(in) :: wavenumber


      specify_constituent % name = name
      specify_constituent % freq = 360./speed
      specify_constituent % k    = wavenumber


    end ! function
!
!______________________________________________________________________
!
    type (T_constituent) function get_constituent( name )

      implicit none

      character(*), intent(in) :: name

      integer i

      get_constituent = specify_constituent( "n/a", 0._rk, 0 )
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
!
!= I/O SECTION ========================================================
!______________________________________________________________________
!
    subroutine read_var_nc( var_name, var, record, ncid )
!----------------------------------------------------------------------
!  Read a variable (NC format).
!______________________________________________________________________
!
      use glob_const , only: C2K, rk
      use glob_domain
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_3d

      integer                   , intent(inout) :: ncid
      integer                   , intent(in   ) :: record
      real(rk), dimension(im,jm), intent(  out) :: var
      character(len=*)          , intent(in   ) :: var_name

      integer                  status, varid
      integer(MPI_OFFSET_KIND) start(4), edge(4)
      character(len=64)        units


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
      if ( status == NF_NOERR ) then
        select case ( trim(units) )
          case ( "K" )
            var = var - C2K
          case ( "Pa" )
            var = var/100.
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

      integer         , intent(in) :: year
      character(len=*), intent(in) :: path

      integer            status
      character(len=256) filename, netcdf_file


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
      use glob_domain, only: error_status, is_master
      use pnetcdf    , only: nf90mpi_strerror, NF_NOERR

      implicit none

      integer         , intent(in) :: status
      character(len=*), intent(in) :: routine


      if ( status /= NF_NOERR ) then
        error_status = 1
        if ( is_master ) then
          print '(/a,a)', 'IO error at module `TIDE`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
!          call finalize_mpi
!          stop
        end if
        call finalize_mpi
        stop
      end if


    end ! subroutine check
!
! INTERPOLATION ROUTINES (TODO: move to separate module)
!______________________________________________________________________
!
    subroutine interpolate_2D( src_x_n, src_y_n, src_x, src_y, src    &
                             , dst_x_n, dst_y_n, dst_x, dst_y, dst, flag )
!----------------------------------------------------------------------
!  Linearly interpolates 2D fields.
!______________________________________________________________________
!
                             use glob_domain, only: my_task
      implicit none

      character                 , intent(in   ) :: flag
      integer                   , intent(in   ) :: dst_x_n, dst_y_n   &
                                                 , src_x_n, src_y_n
      real(rk), dimension(src_x_n,src_y_n)                            &
                                , intent(in   ) :: src, src_x, src_y
      real(rk), dimension(dst_x_n,dst_y_n)                            &
                                , intent(in   ) ::      dst_x, dst_y
      real(rk), dimension(dst_x_n,dst_y_n)                            &
                                , intent(  out) :: dst

      integer                  i, j
      integer, dimension(2) :: x_bl, x_tl, x_tr, x_br  &
                             , y_bl, y_tl, y_tr, y_br, pos


      do j = 1, dst_y_n
        do i = 1, dst_x_n
          pos = get_bottom_left( dst_x(i,j), dst_y(i,j), src_x, src_y, src_x_n, src_y_n )
!          print *, "DST: ", dst_y(i,j), dst_x(i,j)
!          print *, "Y_BL:", y_bl, src_y(y_bl(1),y_bl(2)), src_x(y_bl(1),y_bl(2))
!          print *, "X_BL:", x_bl, src_y(x_bl(1),x_bl(2)), src_x(x_bl(1),x_bl(2))
!          y_tl = minloc( src_y, src_y > dst_y(i,j) )
!          x_tl = minloc( src_x, src_x > dst_x(i,j) )
!          dst(i,j) = src(pos(1),pos(2))
!          write( 40+my_task, * ) dst_x(i,j),";",dst_y(i,j),";",src_x(pos(1),pos(2)),";",src_y(pos(1),pos(2)),";",dst(i,j)
          dst(i,j) = bilint( dst_x(i,j), dst_y(i,j)                                    &
                           , src_x(pos(1)  ,pos(2)  ), src_x(pos(1)  ,pos(2)+1)    &
                           , src_x(pos(1)+1,pos(2)+1), src_x(pos(1)+1,pos(2)  )    &
                           , src_y(pos(1)  ,pos(2)  ), src_y(pos(1)  ,pos(2)+1)    &
                           , src_y(pos(1)+1,pos(2)+1), src_y(pos(1)+1,pos(2)  )    &
                           , src(pos(1)  ,pos(2)  ), src(pos(1)  ,pos(2)+1)    &
                           , src(pos(1)+1,pos(2)+1), src(pos(1)+1,pos(2)  )    &
                           , 0._rk, flag )
        end do
      end do


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

      integer                      , intent(in) :: x_n, y_n
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

!      count = 0
!
!      if ( w1 < 0. .and. w3 < 0. ) then
!        count = 1
!        w1 = 0.
!        w2 = 1.
!        w3 = 0.
!        w4 = 0.
!      end if
!      if ( w3 < 0. .and. w4 < 0. ) then
!        count = 2
!        w3 = 0.
!        w4 = 0.
!      end if
!      if ( w2 < 0. .and. w4 < 0. ) then
!        count = 3
!        w1 = 1.
!        w2 = 0.
!        w3 = 0.
!        w4 = 0.
!      end if
!      if ( w2 < 0. .and. w3 < 0. ) then
!        count = 4
!        w2 = 0.
!        w3 = 0.
!      end if
!      if ( w1 < 0. .and. w3 < 0. ) then
!        count = 5
!        w1 = 1.
!        w2 = 0.
!        w3 = 0.
!        w4 = 1.
!      end if
!      if ( w1 < 0. .and. w2 < 0. ) then
!        count = 6
!        w1 = 0.
!        w2 = 0.
!      end if
!      if ( w2 < 0. .and. w4 < 0. ) then
!        count = 7
!        w1 = 1.
!        w2 = 0.
!        w3 = 1.
!        w4 = 0.
!      end if
!      if ( w1 < 0. .and. w4 < 0. ) then
!        count = 8
!        w1 = 0.
!        w4 = 0.
!      end if

      if ( f1 == fill ) w1 = 0.
      if ( f2 == fill ) w2 = 0.
      if ( f3 == fill ) w3 = 0.
      if ( f4 == fill ) w4 = 0.
!      print *, "FILL: ", (f1==fill), (f2==fill), (f3==fill), (f4==fill)

      tot_w = w1 + w2 + w3 + w4
      if ( tot_w > 0. .and. tot_w /= 1. ) then
        w1 = w1 / tot_w
        w2 = w2 / tot_w 
        w3 = w3 / tot_w
        w4 = w4 / tot_w
      end if

      if ( flag == 'P' ) then
        bilint = 0.
        a2 = f1-f2
        a3 = f1-f3
        a4 = f1-f4
        b1 = f1
        b2 = f2
        b3 = f3
        b4 = f4
        if ( a2 <= -180. ) then
          b2 = f2-360.
        elseif ( a2 >= 180. ) then
          b2 = f2+360.
        end if
        if ( a3 <= -180. ) then
          b3 = f3-360.
        elseif ( a3 >= 180. ) then
          b3 = f3+360.
        end if
        if ( a4 <= -180. ) then
          b4 = f4-360.
        elseif ( a4 >= 180. ) then
          b4 = f4+360.
        end if
        bilint = w1*b1 + w2*b2 + w3*b3 + w4*b4
        if ( bilint > 360. ) then
          bilint = bilint - 360.
        elseif ( bilint < 0. ) then
          bilint = bilint + 360.
        end if
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
