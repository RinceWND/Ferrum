!______________________________________________________________________
!
! Module `RIVER` (river.f90)
!----------------------------------------------------------------------
!  Module for applying river volume influx.
!
!  Author  : RinceWND
!  Created : 2020-07-25
!______________________________________________________________________
!
module river

  use module_time
  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  public

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED  ! Is set according to the external flag `use_river`.

  character(10)                           &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(PATH_LEN)  & ! Full paths (with filenames) to:
    river_path           !  river discharge file

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(VAR_LEN)  &
    discharge_name      ! River discharge

!----------------------------------------------------------------------
! River discharge related arrays
!----------------------------------------------------------------------
  real(rk)                     &
       , allocatable           &
       , dimension(:,:,:) ::   &
    discharge                  ! River discharge


  contains

!______________________________________________________________________
!
    subroutine initialize_mod( config_file )
!----------------------------------------------------------------------
!  Initialize tide module.
!______________________________________________________________________
!
      use config     , only: use_river
      use glob_domain, only: is_master

      implicit none

      character(*), intent(in) :: config_file

      integer pos

      namelist/river/                                       &
        river_path, discharge_name


      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_RIVER ) DISABLED = .true.

      if ( DISABLED ) return

! Initialize variables with their defaults
      river_path     = "rivers.nc"
      discharge_name = "disch"

! Override configuration
      open ( 73, file = config_file, status = 'old' )
      read ( 73, nml = river )
      close( 73 )

      if ( is_master ) then
        print *, "River discharge file: ", trim(river_path)
      end if

! Allocate necessary arrays
      call allocate_arrays

      call msg_print("RIVER MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_air
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________
!
      use glob_domain, only: im, jm
      use grid       , only: art
      use io
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf    , only: NF90_NOWRITE

      implicit none


      allocate( discharge(im,jm,12) )

      discharge = 0.


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

      integer                                    file_id
      integer(MPI_OFFSET_KIND), dimension(3)  :: start, edge


! Quit if the module is not used.
      if ( DISABLED ) return

! Read data
      start = [ i_global(1), j_global(1),  1 ]
      edge  = [          im,          jm, 12 ]

! Read files
      file_id = file_open( river_path, NF90_NOWRITE )
      if ( file_id < 0 ) then
        call msg_print( "", 3, "Failed to open file `"//trim(river_path)//"`" )
      else
        call check( var_read( file_id, trim(discharge_name)  &
                            , discharge, start, edge )       &
                  , "[river]: var_read" )
      end if

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
      use air , only: vfluxf
      use grid, only: art

      implicit none

      type(date), intent(in) :: d_in


! Quit if the module is not used.
      if ( DISABLED ) return

      vfluxf = -discharge(:,:,d_in%month) / art


    end ! subroutine step


end module river
