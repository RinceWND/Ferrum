!______________________________________________________________________
!
! Module `CLIM` (clim.f90)
!----------------------------------------------------------------------
!  Module for managing climatology and restoring.
!
!  Author  : RinceWND
!  Created : 2019-02-04
!______________________________________________________________________
!
module clim

  use glob_const , only: PATH_LEN, rk, VAR_LEN
  use glob_domain, only: im, jm, kb

  implicit none

  public

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  real(rk), private :: a          & ! time-interpolation factor
                     , deep_rel   & ! relaxation period for deep cells
                     , h_thres    & ! apply TS relaxation to cells deeper than set threshold
                     , surf_rel     ! relaxation period for surface

  integer(1)                            & ! Surface relaxation
       , private    :: RELAX_SURF_TEMP  & !  for temperature
                     , RELAX_SURF_SALT    !  for salinity
                                          ! 0 - no relaxation
                                          ! 1 - relax to climatology
                                          ! 2 - relax to sst/sss

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

  logical, private :: INTERP_CLIM  & ! Interpolation flag
                    , RELAX_TS       ! Deep cells relaxation flag
!----------------------------------------------------------------------
! Climatology time mode
!----------------------------------------------------------------------
  integer(1) CLIM_MODE
!______________________________________________________________________
!  0 = Yearly
!  1 = Seasonal
!  2 = Monthly
!  3 = Daily

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=PATH_LEN)  & ! Full paths (with filenames) to:
    ts_clim_path           & ! climatology file for T and S
  , ts_mean_path           & ! horizontally averaged (in z-coord.) TS
  , uv_clim_path             ! climatology file for U and V

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(len=VAR_LEN)  &
    el_clim_name          & ! Surface elevation
  ,  s_clim_name          & ! Salinity
  ,  s_mean_name          & ! Horizontally averaged salinity
  ,  t_clim_name          & ! Temperature
  ,  t_mean_name          & ! Horizontally averaged temperature
  ,  u_clim_name          & ! X-velocity
  , ua_clim_name          & ! Vertically integrated x-velocity
  ,  v_clim_name          & ! Y-velocity
  , va_clim_name            ! Vertically integrated y-velocity

!----------------------------------------------------------------------
! Climatological arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
     eclim                   &
  , s_relx                   & ! mask for salinity surface relaxation
  , t_relx                   & ! mask for temperature surface relaxation
  , uaclim                   &
  , vaclim
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:,:) :: &
    rmean                    &
  , sclim                    &
  , smean                    &
  , tclim                    &
  , tmean                    &
  , uclim                    &
  , vclim

!----------------------------------------------------------------------
! Private interpolation arrays
!----------------------------------------------------------------------
  real(kind=rk)                &
       , private               &
       , allocatable           &
       , dimension(:,:,:) ::   &
    el_int                     &
  , ua_int                     &
  , va_int
  real(kind=rk)                &
       , private               &
       , allocatable           &
       , dimension(:,:,:,:) :: &
    sc_int                     &
  , sm_int                     &
  , tc_int                     &
  , tm_int                     &
  , uc_int                     &
  , vc_int


  contains

!______________________________________________________________________
!
    subroutine initialize_mod( config_file )
!----------------------------------------------------------------------
!  Initialize clim module.
!______________________________________________________________________
!
      use glob_domain, only: is_master
      use glob_grid  , only: h

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/clim/                                                  &
        CLIM_MODE      , deep_rel       , h_thres     , INTERP_CLIM   &
      , RELAX_SURF_SALT, RELAX_SURF_TEMP, RELAX_TS    , surf_rel      &
      , ts_clim_path   , ts_mean_path   , uv_clim_path

      namelist/clim_vars/                        &
        el_clim_name, s_clim_name,  s_mean_name  &
      ,  t_clim_name, t_mean_name,  u_clim_name  &
      , ua_clim_name, v_clim_name, va_clim_name


! Initialize with default values
      INTERP_CLIM = .true.
      RELAX_TS    = .true.
      RELAX_SURF_SALT = 0
      RELAX_SURF_TEMP = 0

      h_thres  = 1000._rk
      deep_rel = 15._rk
      surf_rel =  2._rk

      CLIM_MODE = 2

      ts_clim_path = "in/clim/"
      ts_mean_path = "in/clim/"
      uv_clim_path = "in/clim/"

      el_clim_name = "eclim"
       s_clim_name = "sclim"
       s_mean_name = "smean"
       t_clim_name = "tclim"
       t_mean_name = "tmean"
       u_clim_name = "uclim"
      ua_clim_name = "uaclim"
       v_clim_name = "vclim"
      va_clim_name = "vaclim"

! Override configuration
      open ( 73, file = config_file, status = 'old' )
      read ( 73, nml = clim )
      read ( 73, nml = clim_vars )
      close( 73 )

! Variables management
      pos = len(trim(ts_clim_path))
      if ( ts_clim_path(pos:pos) == "/" ) then
        ts_clim_path = trim(ts_clim_path)//"clim."
      end if

      pos = len(trim(ts_mean_path))
      if ( ts_mean_path(pos:pos) == "/" ) then
        ts_mean_path = trim(ts_mean_path)//"mean."
      end if

      pos = len(trim(uv_clim_path))
      if ( uv_clim_path(pos:pos) == "/" ) then
        uv_clim_path = trim(uv_clim_path)//"clim."
      end if

      if ( is_master ) then
        print *, "Climatology TS: ", trim(ts_clim_path)
        print *, "Background TS : ", trim(ts_mean_path)
        print *, "Climatology UV: ", trim(uv_clim_path)
      end if

! Allocate necessary arrays
      call allocate_arrays

! Always allow surface relaxation only over deep ocean > 1000m (TODO: Is it really necessary?)
      if ( RELAX_SURF_SALT > 0 ) then
        s_relx = ( 1._rk + tanh( .002_rk*(h-1000._rk) ) )*.5_rk
      end if

      if ( RELAX_SURF_TEMP > 0 ) then
        t_relx = ( 1._rk + tanh( .002_rk*(h-1000._rk) ) )*.5_rk
      end if

! Convert relaxation periods to frequencies
      deep_rel = 1._rk / (deep_rel*86400._rk)
      surf_rel = 1._rk / (surf_rel*86400._rk)

      call msg_print("CLIM MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_mod
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________
!
      implicit none

      integer(1) N         ! Interpolation array extension size
                           ! The structure is following:
                           !   var at n-1 step: var(:,:,2)
                           !   var at n+1 step: var(:,:,3)

      N = 2
      if ( .not.INTERP_CLIM ) N = 1

! Allocate core arrays
      allocate(          &
        eclim(im,jm)     &
      , rmean(im,jm,kb)  &
      , sclim(im,jm,kb)  &
      , smean(im,jm,kb)  &
      , tclim(im,jm,kb)  &
      , tmean(im,jm,kb)  &
      , uclim(im,jm,kb)  &
      , uaclim(im,jm)    &
      , vclim(im,jm,kb)  &
      , vaclim(im,jm)    &
       )

! Initialize mandatory arrays
      eclim = 0.
      sclim = 0.
      smean = 0.
      tclim = 0.
      tmean = 0.
      uclim = 0.
      uaclim= 0.
      vclim = 0.
      vaclim= 0.

! Allocate interpolation arrays
      allocate(                 &
        el_int(im,jm,   2:N+1)  &
      , sc_int(im,jm,kb,2:N+1)  &
      , sm_int(im,jm,kb,2:N+1)  &
      , tc_int(im,jm,kb,2:N+1)  &
      , tm_int(im,jm,kb,2:N+1)  &
      , ua_int(im,jm,   2:N+1)  &
      , uc_int(im,jm,kb,2:N+1)  &
      , va_int(im,jm,   2:N+1)  &
      , vc_int(im,jm,kb,2:N+1)  &
       )

! Allocate depth-relaxation arrays
      if ( RELAX_SURF_SALT > 0 ) then
        allocate(        &
          s_relx(im,jm)  &
         )
      end if

      if ( RELAX_SURF_TEMP > 0 ) then
        allocate(        &
          t_relx(im,jm)  &
         )
      end if


    end ! subroutine allocate_arrays
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields before experiment's start.
!______________________________________________________________________
!
      use glob_domain, only: is_master
      use module_time
      use model_run  , only: sec_of_month, mid_in_month, mid_in_nbr

      implicit none

      type(date), intent(in) :: d_in

      logical                  fexist
      integer, dimension(3) :: record


! Determine the record to read
      record(1) = d_in%month

      if ( INTERP_CLIM ) then

        if ( sec_of_month <= mid_in_month ) then
          record(2) = d_in%month - 1
          if ( record(2) == 0 ) record(2) = 12
          a = real( sec_of_month + mid_in_nbr   )  &
             /real( mid_in_nbr   + mid_in_month )
        else
          record(2) = d_in%month
          a = real( sec_of_month - mid_in_month )  &
             /real( mid_in_nbr   + mid_in_month )
        end if

        record(3) = mod( record(2), 12 ) + 1

      else
        record(2) = record(1)
      end if

! Read climatology
      inquire ( file = trim(ts_clim_path), exist = fexist )
      call msg_print("", 6, "Read TS climatology:")
      if ( fexist ) then
        call read_clim_ts_pnetcdf( tc_int(:,:,:,2), sc_int(:,:,:,2)  &
                                 , ts_clim_path, record(2) )
        if ( INTERP_CLIM ) then
          call read_clim_ts_pnetcdf( tc_int(:,:,:,3), sc_int(:,:,:,3)  &
                                   , ts_clim_path, record(3) )
        end if
      else
        call msg_print("", 2, "FAILED...")
        tc_int = 15.
        sc_int = 33.
      end if

! Read background TS
      inquire ( file = trim(ts_mean_path), exist = fexist )
      call msg_print("", 6, "Read background TS:")
      if ( fexist ) then
        call read_mean_ts_pnetcdf( tm_int(:,:,:,2), sm_int(:,:,:,2)  &
                                 , ts_mean_path, record(2) )
        if ( INTERP_CLIM ) then
          call read_mean_ts_pnetcdf( tm_int(:,:,:,3), sm_int(:,:,:,3)  &
                                   , ts_mean_path, record(3) )
        end if
      else
        call msg_print("", 2, "FAILED...")
        tm_int = tc_int
        sm_int = sc_int
      end if

      if ( INTERP_CLIM ) then
! Interpolate TS and get background density
        sclim = ( 1. - a )*sc_int(:,:,:,2) + a*sc_int(:,:,:,3)
        tclim = ( 1. - a )*tc_int(:,:,:,2) + a*tc_int(:,:,:,3)
        smean = ( 1. - a )*sm_int(:,:,:,2) + a*sm_int(:,:,:,3)
        tmean = ( 1. - a )*tm_int(:,:,:,2) + a*tm_int(:,:,:,3)
      else
        sclim = sc_int(:,:,:,2)
        tclim = tc_int(:,:,:,2)
        smean = sm_int(:,:,:,2)
        tmean = tm_int(:,:,:,2)
      end if

      call dens( smean, tmean, rmean )

! TODO: Gather stats from all processors
      if ( is_master ) then
        print '(/a,f7.3,a,f7.3,a)',"Background temperature:     ("  &
                                 , minval(tmean), ":"              &
                                 , maxval(tmean), ")"
        print '(a,f7.3,a,f7.3,a)', "Background salinity:        ("  &
                                 , minval(smean), ":"              &
                                 , maxval(smean), ")"
        print '(a,f7.3,a,f7.3,a)', "Climatological temperature: ("  &
                                 , minval(tclim), ":"              &
                                 , maxval(tclim), ")"
        print '(a,f7.3,a,f7.3,a)', "Climatological salinity:    ("  &
                                 , minval(sclim), ":"              &
                                 , maxval(sclim), ")"
!        print '(a,f7.3,a,f7.3,a)', "Climatological x-current:   ("  &
!                                 , minval(uc_int), ":"              &
!                                 , maxval(uc_int), ")"
!        print '(a,f7.3,a,f7.3,a)', "Climatological y-current:   ("  &
!                                 , minval(vc_int), ":"              &
!                                 , maxval(vc_int), ")"
!        print '(a,f7.3,a,f7.3,a)', "Climatological 2D x-current:("  &
!                                 , minval(ua_int), ":"              &
!                                 , maxval(ua_int), ")"
!        print '(a,f7.3,a,f7.3,a)', "Climatological 2D y-current:("  &
!                                 , minval(va_int), ":"              &
!                                 , maxval(va_int), ")"
      end if

      call msg_print("CLIM INITIALIZED", 2, "")


    end ! subroutine init
!
!______________________________________________________________________
!
    subroutine step( d_in )
!----------------------------------------------------------------------
!  Reads climatology fields.
!______________________________________________________________________
!
!      use glob_const , only: SEC2DAY
!      use glob_domain, only: is_master
      use module_time
      use model_run  , only: dti, iint                               &
                           , mid_in_month, mid_in_nbr, sec_of_month

      implicit none

      type(date), intent(in) :: d_in

      logical                   ADVANCE_REC_INT, fexist
      integer, dimension(3)  :: record


! Determine the record to read
      record(1) = d_in%month

      if ( sec_of_month <= mid_in_month ) then
        record(2) = d_in%month - 1
        if ( record(2) == 0 ) record(2) = 12
        a = real( sec_of_month + mid_in_nbr   )  &
           /real( mid_in_nbr   + mid_in_month )
      else
        record(2) = d_in%month
        a = real( sec_of_month - mid_in_month )  &
           /real( mid_in_nbr   + mid_in_month )
      end if

      record(3) = mod( record(2), 12 ) + 1

      ADVANCE_REC_INT = .false.
      if ( INTERP_CLIM ) then
        if ( sec_of_month - mid_in_month <= dti .and.  &
             record(2) == record(1) ) then
          if ( iint > 1 ) ADVANCE_REC_INT = .true.
        end if
      else
        if ( sec_of_month <= dti ) ADVANCE_REC_INT = .true.
      end if

      if ( ADVANCE_REC_INT ) then

        if ( INTERP_CLIM ) then
          sc_int(:,:,:,2) = sc_int(:,:,:,3)
          tc_int(:,:,:,2) = tc_int(:,:,:,3)
          sm_int(:,:,:,2) = sm_int(:,:,:,3)
          tm_int(:,:,:,2) = tm_int(:,:,:,3)
        end if

        ! Read climatology
        inquire ( file = trim(ts_clim_path), exist = fexist )
        call msg_print("", 6, "Read TS climatology:")
        if ( fexist ) then
          if ( INTERP_CLIM ) then
            call read_clim_ts_pnetcdf( tc_int(:,:,:,3), sc_int(:,:,:,3)  &
                                     , ts_clim_path, record(3) )
          else
            call read_clim_ts_pnetcdf( tc_int(:,:,:,2), sc_int(:,:,:,2)  &
                                     , ts_clim_path, record(1) )
          end if
        else
          call msg_print("", 2, "FAILED...")
          tc_int(:,:,:,3) = 15.
          sc_int(:,:,:,3) = 33.
        end if

  ! Read background TS
        inquire ( file = trim(ts_mean_path), exist = fexist )
        call msg_print("", 6, "Read background TS:")
        if ( fexist ) then
          if ( INTERP_CLIM ) then
            call read_mean_ts_pnetcdf( tm_int(:,:,:,3), sm_int(:,:,:,3)  &
                                     , ts_mean_path, record(3) )
          else
            call read_mean_ts_pnetcdf( tm_int(:,:,:,2), sm_int(:,:,:,2)  &
                                     , ts_mean_path, record(1) )
          end if
        else
          call msg_print("", 2, "FAILED...")
          tm_int(:,:,:,3) = tc_int(:,:,:,3)
          sm_int(:,:,:,3) = sc_int(:,:,:,3)
        end if

      end if

! Interpolate TS and get background density
      if ( INTERP_CLIM ) then
        sclim = ( 1. - a )*sc_int(:,:,:,2) + a*sc_int(:,:,:,3)
        tclim = ( 1. - a )*tc_int(:,:,:,2) + a*tc_int(:,:,:,3)
        smean = ( 1. - a )*sm_int(:,:,:,2) + a*sm_int(:,:,:,3)
        tmean = ( 1. - a )*tm_int(:,:,:,2) + a*tm_int(:,:,:,3)
      else
        sclim = sc_int(:,:,:,2)
        tclim = tc_int(:,:,:,2)
        sclim = sm_int(:,:,:,2)
        tclim = tm_int(:,:,:,2)
      end if

      call dens( smean, tmean, rmean )


    end ! subroutine step
!
!______________________________________________________________________
!
    subroutine relax_to_clim( temp, salt )
!----------------------------------------------------------------------
!  Relaxes temperature and salinity to climatology.
!______________________________________________________________________
!
      use glob_ocean , only: hz

      implicit none

      real(rk), dimension(im,jm,kb), intent(inout) :: salt, temp


      if ( .not.RELAX_TS ) return

      where ( hz .gt. h_thres )
        temp = deep_rel * ( tclim - temp )
        salt = deep_rel * ( sclim - salt )
      end where


    end ! subroutine relax_ts
!
!______________________________________________________________________
!
    subroutine relax_surface( wssurf, wtsurf, sss, sst )
!----------------------------------------------------------------------
!  Relax surface fluxes to climatologies
!______________________________________________________________________
!
      use glob_ocean, only: sb, tb

      implicit none

      real(rk), dimension(im,jm), intent(in   ) :: sss, sst
      real(rk), dimension(im,jm), intent(inout) :: wssurf, wtsurf


      select case ( RELAX_SURF_SALT )
        case ( 0 )
          return
        case ( 1 )
          wssurf = wssurf + deep_rel*s_relx*( sb(:,:,1) - sclim(:,:,1) )
        case ( 2 )
          wssurf = wssurf + surf_rel * ( sb(:,:,1) - sss )
      end select

      select case ( RELAX_SURF_TEMP )
        case ( 0 )
          return
        case ( 1 )
          wtsurf = wtsurf + deep_rel*t_relx*( tb(:,:,1) - tclim(:,:,1) )
        case ( 2 )
          wtsurf = wtsurf + surf_rel * ( tb(:,:,1) - sst )
      end select


    end ! subroutine relax_surface
!
!______________________________________________________________________
!
    pure character(len=256) function get_filename( path, year )
!----------------------------------------------------------------------
!  Constructs filename string in `<path>YYYY<FORMAT_EXT>` format.
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
          print '(/a,a)', 'IO error at module `CLIM`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
          stop
        end if
      end if


    end ! subroutine check


end module clim
