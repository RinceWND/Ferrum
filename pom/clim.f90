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
  use glob_domain, only: im, is_master, jm, km, kmm1

  implicit none

  public

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical              NO_CLIM, NO_MEAN ! Flags for climatology and background files existence

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
  , uv_clim_path           & ! climatology file for U and V
  , wattype_path             ! water type spatial distribution file

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
  , va_clim_name          & ! Vertically integrated y-velocity
  , wattype_name            ! Water type id

!----------------------------------------------------------------------
! Climatological and additional arrays
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

  integer                    &
       , allocatable         &
       , dimension(:,:)   :: &
    wtype                      ! Water type (in ids of ntp)

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
      use config     , only: ntp
      use glob_domain, only: is_master
      use grid       , only: h
      use model_run  , only: dti

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/clim/                                                  &
        CLIM_MODE      , deep_rel       , h_thres     , INTERP_CLIM   &
      , RELAX_SURF_SALT, RELAX_SURF_TEMP, RELAX_TS    , surf_rel      &
      , ts_clim_path   , ts_mean_path   , uv_clim_path, wattype_path

      namelist/clim_vars/                        &
        el_clim_name, s_clim_name,  s_mean_name  &
      ,  t_clim_name, t_mean_name,  u_clim_name  &
      , ua_clim_name, v_clim_name, va_clim_name  &
      , wattype_name


! Initialize with default values
      INTERP_CLIM = .true.
      RELAX_TS    = .true.
      RELAX_SURF_SALT = 0
      RELAX_SURF_TEMP = 0

      NO_CLIM = .false.
      NO_MEAN = .false.

      h_thres  = 1000._rk
      deep_rel = 15._rk
      surf_rel =  2._rk

      CLIM_MODE = 2

      ts_clim_path = "in/clim/"
      ts_mean_path = "in/clim/"
      uv_clim_path = "in/clim/"
      wattype_path = ""

      el_clim_name = "eclim"
       s_clim_name = "sclim"
       s_mean_name = "smean"
       t_clim_name = "tclim"
       t_mean_name = "tmean"
       u_clim_name = "uclim"
      ua_clim_name = "uaclim"
       v_clim_name = "vclim"
      va_clim_name = "vaclim"
      wattype_name = "ntp"

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

      deep_rel = max( deep_rel, dti/86400._rk )
      surf_rel = max( surf_rel, dti/86400._rk )

      if ( is_master ) then
        print *, "Climatology TS: ", trim(ts_clim_path)
        print *, "Background TS : ", trim(ts_mean_path)
        print *, "Climatology UV: ", trim(uv_clim_path)
        print *, "-----"
        print '(a,f8.2)', "Deep (slow) relaxation period:    ", deep_rel
        print '(a,f8.2)', "Surface (fast) relaxation period: ", surf_rel
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

! Initialise water type
      wtype = ntp

! Convert relaxation periods to frequencies
      deep_rel = 1._rk / (deep_rel*86400._rk)
      surf_rel = 1._rk / (surf_rel*86400._rk)

      call msg_print("CLIM MODULE READY", 1, "")


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
      , rmean(im,jm,km)  &
      , sclim(im,jm,km)  &
      , smean(im,jm,km)  &
      , tclim(im,jm,km)  &
      , tmean(im,jm,km)  &
      , uclim(im,jm,km)  &
      , uaclim(im,jm)    &
      , vclim(im,jm,km)  &
      , vaclim(im,jm)    &
      , wtype(im,jm)     &
       )

! Initialize mandatory arrays
      eclim = 0.
      rmean = 0.
      sclim = 0.
      smean = 0.
      tclim = 0.
      tmean = 0.
      uclim = 0.
      uaclim= 0.
      vclim = 0.
      vaclim= 0.
      wtype = 0

! Allocate interpolation arrays
      allocate(                 &
        el_int(im,jm,   2:N+1)  &
      , sc_int(im,jm,km,2:N+1)  &
      , sm_int(im,jm,km,2:N+1)  &
      , tc_int(im,jm,km,2:N+1)  &
      , tm_int(im,jm,km,2:N+1)  &
      , ua_int(im,jm,   2:N+1)  &
      , uc_int(im,jm,km,2:N+1)  &
      , va_int(im,jm,   2:N+1)  &
      , vc_int(im,jm,km,2:N+1)  &
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
      use config     , only: ntp
      use glob_domain, only: is_master
      use io
      use module_time
      use model_run  , only: sec_of_month, mid_in_month, mid_in_nbr
      use glob_ocean , only: elb, sb, tb, ub, vb
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      type(date), intent(in) :: d_in

      integer, dimension(3) :: record
      integer                  status


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

      el_int(:,:  ,2) = elb
      uc_int(:,:,:,2) = ub
      vc_int(:,:,:,2) = vb

      if ( INTERP_CLIM ) then
        el_int(:,:  ,3) = elb
        uc_int(:,:,:,3) = ub
        vc_int(:,:,:,3) = vb
      end if

! Read climatology
      call msg_print("", 6, "Read TS climatology:")
      status = read_clim( trim(ts_clim_path), record(2)          &
                        , t_clim_name       , s_clim_name        &
                        , tc_int(:,:,:,2)   , sc_int(:,:,:,2) )
      if ( status /= 0 ) then
        NO_CLIM = .true.
        call msg_print("", 2, "CLIMATOLOGY TO BE DERIVED FROM INITIAL CONDITIONS")
        tc_int(:,:,:,2) = tb
        sc_int(:,:,:,2) = sb
        if ( INTERP_CLIM ) then
          tc_int(:,:,:,3) = tb
          sc_int(:,:,:,3) = sb
        end if
      elseif ( INTERP_CLIM ) then
        status = read_clim( trim(ts_clim_path), record(3)          &
                          , t_clim_name       , s_clim_name        &
                          , tc_int(:,:,:,3)   , sc_int(:,:,:,3) )
      end if

! Read background TS
      call msg_print("", 6, "Read background TS:")
      status = read_clim( trim(ts_mean_path), record(2)          &
                        , t_mean_name       , s_mean_name        &
                        , tm_int(:,:,:,2)   , sm_int(:,:,:,2) )
      if ( status /= 0 ) then
        NO_MEAN = .true.
        tm_int =  5._rk
        sm_int = 32._rk
      elseif ( INTERP_CLIM ) then
        status = read_clim( trim(ts_mean_path), record(3)          &
                          , t_mean_name       , s_mean_name        &
                          , tm_int(:,:,:,3)   , sm_int(:,:,:,3) )
      end if

! TODO: IS the below part actually necessary? We should be able to derive density at the beginning of a step.
      if ( INTERP_CLIM ) then
! Interpolate TS and get background density
        if ( .not.NO_CLIM ) then
          sclim = ( 1. - a )*sc_int(:,:,:,2) + a*sc_int(:,:,:,3)
          tclim = ( 1. - a )*tc_int(:,:,:,2) + a*tc_int(:,:,:,3)
        end if
        if ( .not.NO_MEAN ) then
          smean = ( 1. - a )*sm_int(:,:,:,2) + a*sm_int(:,:,:,3)
          tmean = ( 1. - a )*tm_int(:,:,:,2) + a*tm_int(:,:,:,3)
        end if
        eclim = ( 1. - a )*el_int(:,:,2) + a*el_int(:,:,3)
        uclim = ( 1. - a )*uc_int(:,:,:,2) + a*uc_int(:,:,:,3)
        vclim = ( 1. - a )*vc_int(:,:,:,2) + a*vc_int(:,:,:,3)
      else

        if ( .not.NO_CLIM ) then
          sclim = sc_int(:,:,:,2)
          tclim = tc_int(:,:,:,2)
        end if
        if ( .not.NO_MEAN ) then
          smean = sm_int(:,:,:,2)
          tmean = tm_int(:,:,:,2)
        end if

        eclim = el_int(:,:,2)
        uclim = uc_int(:,:,:,2)
        vclim = vc_int(:,:,:,2)

      end if

      call dens( smean, tmean, rmean )

! Read water type data if specifed
      if ( wattype_path /= "" ) then
        status = read_wtype( trim(wattype_path), record(1), wattype_name, wtype )
        where ( wtype == 0 )
          wtype = ntp
        end where
      else
        wtype = ntp
      end if

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

      logical                   ADVANCE_REC_INT
      integer, dimension(3)  :: record
      integer                   status


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
          if ( .not.NO_CLIM ) then
            sc_int(:,:,:,2) = sc_int(:,:,:,3)
            tc_int(:,:,:,2) = tc_int(:,:,:,3)
          end if
          if ( .not.NO_MEAN ) then
            sm_int(:,:,:,2) = sm_int(:,:,:,3)
            tm_int(:,:,:,2) = tm_int(:,:,:,3)
          end if
        end if

        if ( .not.NO_CLIM ) then
! Read climatology
          call msg_print("", 6, "Read TS climatology:")
          if ( INTERP_CLIM ) then
            status = read_clim( trim(ts_clim_path), record(3)          &
                              , t_clim_name       , s_clim_name        &
                              , tc_int(:,:,:,3)   , sc_int(:,:,:,3) )
          else
            status = read_clim( trim(ts_clim_path), record(1)          &
                              , t_clim_name       , s_clim_name        &
                              , tc_int(:,:,:,2)   , sc_int(:,:,:,2) )
          end if
        end if

        if ( .not.NO_MEAN ) then
! Read background TS
          call msg_print("", 6, "Read background TS:")
          if ( INTERP_CLIM ) then
            status = read_clim( trim(ts_mean_path), record(3)          &
                              , t_mean_name       , s_mean_name        &
                              , tm_int(:,:,:,3)   , sm_int(:,:,:,3) )
          else
            status = read_clim( trim(ts_mean_path), record(1)          &
                              , t_mean_name       , s_mean_name        &
                              , tm_int(:,:,:,2)   , sm_int(:,:,:,2) )
          end if
        end if

      end if

! Interpolate TS and get background density
      if ( INTERP_CLIM ) then
        if ( .not.NO_CLIM ) then
          sclim = ( 1. - a )*sc_int(:,:,:,2) + a*sc_int(:,:,:,3)
          tclim = ( 1. - a )*tc_int(:,:,:,2) + a*tc_int(:,:,:,3)
        end if
        if ( .not.NO_MEAN ) then
          smean = ( 1. - a )*sm_int(:,:,:,2) + a*sm_int(:,:,:,3)
          tmean = ( 1. - a )*tm_int(:,:,:,2) + a*tm_int(:,:,:,3)
        end if
      else
        if ( .not.NO_CLIM ) then
          sclim = sc_int(:,:,:,2)
          tclim = tc_int(:,:,:,2)
        end if
        if ( .not.NO_MEAN ) then
          smean = sm_int(:,:,:,2)
          tmean = tm_int(:,:,:,2)
        end if
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

      real(rk), dimension(im,jm,km), intent(inout) :: salt, temp


      if ( .not.RELAX_TS ) return

      where ( hz .gt. h_thres )
        temp = temp + deep_rel * ( tclim - temp )
        salt = salt + deep_rel * ( sclim - salt )
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
        case ( 1 )
          wssurf = wssurf + deep_rel*s_relx*( sb(:,:,1) - sclim(:,:,1) )
        case ( 2 )
          wssurf = wssurf + surf_rel * ( sb(:,:,1) - sss )
      end select

      select case ( RELAX_SURF_TEMP )
        case ( 0 )
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
    integer function read_clim( filepath , record       &
                              , temp_name, salt_name    &
                              , temp     , salt      )

      use glob_domain, only: i_global, j_global
      use grid       , only: fsm
      use io
      use mpi        , only: MPI_OFFSET_KIND
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      character(*), intent(in)          :: filepath, salt_name, temp_name
      integer     , intent(in)          :: record
      real(rk)    , dimension(im,jm,km)  &
                  , intent(out)         :: temp, salt

      character(PATH_LEN)                    :: tmp_str
      integer(MPI_OFFSET_KIND), dimension(4) :: edge, start
      integer file_id


      tmp_str = ""

! open file
      if ( is_master ) then
        write(tmp_str, '("Reading file ",a," @ ",i2)') filepath, record
        call msg_print("", 6, trim(tmp_str))
      end if

      file_id = file_open( filepath, NF90_NOWRITE )
      if ( file_id < 0 ) then
        read_clim = -1
        return
      end if

! define domain to read
      start = [ i_global(1), j_global(1),    1, record ]
      edge  = [ im         , jm         , kmm1,      1 ]

! get data
      read_clim = var_read( file_id, temp_name, temp, start, edge )
      read_clim = var_read( file_id, salt_name, salt, start, edge )

      where ( fsm == 0. )
        temp = 0.
        salt = 0.
      end where
      ! temp(:,:,kb) = temp(:,:,kmm1)
      ! salt(:,:,kb) = salt(:,:,kmm1)

! close file
      file_id = file_close( file_id )


    end ! function read_clim
!
!___________________________________________________________________
!
    integer function read_wtype( filepath, record, name, wtype )

      use config     , only: ntp
      use glob_domain, only: i_global, j_global
      use grid       , only: fsm
      use io
      use mpi        , only: MPI_OFFSET_KIND
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      character(*), intent(in)          :: filepath, name
      integer     , intent(in)          :: record
      integer     , dimension(im,jm)  &
                  , intent(out)         :: wtype

      character(PATH_LEN)                    :: tmp_str
      integer(MPI_OFFSET_KIND), dimension(3) :: edge, start
      integer file_id


      tmp_str = ""

! open file
      if ( is_master ) then
        write(tmp_str, '("Reading file ",a," @ ",i2)') filepath, record
        call msg_print("", 6, trim(tmp_str))
      end if

      file_id = file_open( filepath, NF90_NOWRITE )
      if ( file_id < 0 ) then
        read_wtype = -1
        return
      end if

! define domain to read
      start = [ i_global(1), j_global(1), record ]
      edge  = [ im         , jm         ,      1 ]

      if ( var_rank( file_id, name ) == 2 ) then
        start(3) = 1
      end if

! get data
      read_wtype = var_read( file_id, name, wtype, start, edge )

! close file
      file_id = file_close( file_id )


    end ! function read_wtype
!
!___________________________________________________________________
!

end module clim
