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

  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  public

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  integer, private :: read_int  ! interval for reading (days)
  real   , private :: a         ! time-interpolation factor

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

!----------------------------------------------------------------------
! Climatology time mode
!----------------------------------------------------------------------
  integer(1) CLIM_MODE
!______________________________________________________________________
!  0 = Yearly
!  1 = Seasonal
!  2 = Monthly

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
  , uaclim                   &
  , vaclim
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:,:) :: &
    sclim                    &
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

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/clim/                                         &
        CLIM_MODE, ts_clim_path, ts_mean_path, uv_clim_path

      namelist/clim_vars/                        &
        el_clim_name, s_clim_name,  s_mean_name  &
      ,  t_clim_name, t_mean_name,  u_clim_name  &
      , ua_clim_name, v_clim_name, va_clim_name


! Initialize with default values
      CLIM_MODE = 3

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

      call msg_print("CLIM MODULE INITIALIZED", 1, "")

    end ! subroutine initialize_air
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________
!
      use glob_domain, only: im, jm, kb

      implicit none

      integer*1 N ! Interpolation array extension size
                  ! The structure is following:
                  !   var at n-1 step: var(:,:,1)
                  !   var at n+1 step: var(:,:,2)


! Allocate core arrays
      allocate(          &
        eclim(im,jm)     &
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
      allocate(             &
        el_int(im,jm,   N)  &
      , sc_int(im,jm,kb,N)  &
      , sm_int(im,jm,kb,N)  &
      , tc_int(im,jm,kb,N)  &
      , tm_int(im,jm,kb,N)  &
      , ua_int(im,jm,   N)  &
      , uc_int(im,jm,kb,N)  &
      , va_int(im,jm,   N)  &
      , vc_int(im,jm,kb,N)  &
       )


    end ! subroutine allocate_arrays
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields before experiment's start.
!______________________________________________________________________
!
!      use glob_domain, only: is_master
!      use config     , only: calc_wind
      use module_time
!      use glob_ocean , only: tb

      implicit none

      type(date), intent(in) :: d_in
      
      integer, dimension(3) :: record, year

! Read wind if momentum flux is derived from wind
      call read_all( .true., 1, year, record )
      call read_all( .true., 2, year, record )

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
      use model_run  , only: dti, iint, sec_of_year

      implicit none

      type(date), intent(in) :: d_in

      logical            ADVANCE_REC, ADVANCE_REC_INT
      integer            max_in_prev, max_in_this, secs
      integer                          &
      , dimension(3)  :: record, year
      real               chunk



      call read_all( ADVANCE_REC_INT, 3, year, record )
      call read_all( ADVANCE_REC    , 1, year, record )

    end ! subroutine step
!
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

      if ( n >= 2 ) then
        write(desc,'("Reading interp. record #",i4," @ ",i4)') &
            record(n), year(n)
      else
        write(desc,'("Reading clim. record #",i4," @ ",i4)') &
            record(1), year(1)
      end if

      call msg_print("", 1, desc)


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
          print '(/a,a)', 'IO error at module `AIR`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
          stop
        end if
      end if


    end ! subroutine check


end module clim
