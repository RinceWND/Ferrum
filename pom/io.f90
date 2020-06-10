!______________________________________________________________________
!
! Module `IO` (io.f90)
!----------------------------------------------------------------------
!  Module for IO operations. Generally for output and restart IO.
! Other module-specific subroutines should be in their respective
! modules.
!  This module can be switched to another one that uses file format
! suitable for user. All the subroutines should have the same
! signature of course. Otherwise you need to edit the source code.
!
!  Author  : RinceWND
!  Created : 2018-09-07
!______________________________________________________________________
!
module io

  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  interface att_write
    module procedure att_write_chr
    module procedure att_write_int
    module procedure att_write_flt
  end interface

  interface var_read
    module procedure var_read_0d
    module procedure var_read_1d
    module procedure var_read_2d
    module procedure var_read_3d
  end interface

  interface var_write
    module procedure var_write_0d
    module procedure var_write_1d
    module procedure var_write_2d
    module procedure var_write_3d
  end interface


  public :: att_write , check     , dim_define , is_error             &
          , file_open , file_close, file_create, file_end_definition  &
          , var_define, var_rank  , var_read   , var_shape            &
          , var_write

  private

  character(len=10), parameter :: FORMAT_EXT = "nc" ! Output format extension


  contains
!
!______________________________________________________________________
!
    integer function var_define( ncid     , name , vartype, vdims    &
                               , long_name, units, nofill , fillval  &
                               , coords                              )
!----------------------------------------------------------------------
!  Defines variable in NetCDF file.
!______________________________________________________________________
!
      use pnetcdf, only: nf90mpi_def_var, nf90mpi_def_var_fill  &
                       , nf90mpi_put_att                        &
                       , NF_NOERR

      integer     , intent(in)   :: ncid       & ! ID of NetCDF file (in define mode)
                                  , vartype    & ! Variable's type
                                  , nofill       ! Do not write missing values
                                                 !  (useful for static masks)

      integer     , intent(in)                 &
                  , dimension(:) :: vdims        ! Variable's dimensions array

      character(*), intent(in)   :: name       & ! Variable's name (ID)
                                  , long_name  & ! Variable's human-readable name
                                  , units      & ! Variable's units
                                  , coords       ! Variable's coordinates (opt.)

      real        , intent(in)   :: fillval      ! Missing value


! define variable
      call check( nf90mpi_def_var( ncid, name, vartype  &
                                 , vdims, var_define )  &
                , "var_define:def_var `"//name//"`" )
! define fill value
      if ( nofill > -1 ) then
        call check( nf90mpi_def_var_fill( ncid, var_define   &
                                        , nofill, fillval )  &
                  , "var_define:def_var_fill `"//name//"`" )
      end if
! define attributes
      call check( nf90mpi_put_att( ncid, var_define, 'long_name'  &
     &                           , trim(long_name) )              &
                , "var_define:put_att("//long_name//") @ `"//name//"`" )
      call check( nf90mpi_put_att( ncid, var_define, 'units'  &
     &                           ,trim(units) )               &
                , "var_define:put_att("//units//") @ `"//name//"`" )
! add coordinates attribute, if necessary
      if ( coords /= "" ) then
        call check( nf90mpi_put_att( ncid, var_define, 'coordinates'  &
     &                             , trim(coords) )                   &
                  , "var_define:put_att("//coords//") @ `"//name//"`" )
      end if


    end ! function var_define
!
!___________________________________________________________________
!
    integer function dim_define( ncid, name, length )
!-------------------------------------------------------------------
!  Creates netcdf file.
!___________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_UNLIMITED, nf90mpi_def_dim

      implicit none

      integer     , intent(in) :: ncid
      integer(4)  , intent(in) :: length
      character(*), intent(in) :: name

      integer(MPI_OFFSET_KIND) dim_len


      dim_len = max( int(length        ,MPI_OFFSET_KIND)    &
                   , int(NF90_UNLIMITED,MPI_OFFSET_KIND) )
      call check( nf90mpi_def_dim( ncid, name, dim_len, dim_define )  &
                , 'dim_define:def_dim '//name )


    end ! function dim_define
!
!___________________________________________________________________
!
    integer function var_rank( ncid, name )

      use pnetcdf, only: nf90mpi_inq_varid, nf90mpi_inquire_variable

      implicit none

      integer     , intent(in) :: ncid
      character(*), intent(in) :: name

      integer varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , "inq_varid:var_rank "//name )
      call check( nf90mpi_inquire_variable( ncid, varid, ndims = var_rank ) &
                , "inquire_variable(ndims):var_rank "//name )


    end function
!
!___________________________________________________________________
!
    function var_shape( ncid, name )
!-------------------------------------------------------------------
!  Returns a variable's shape.
!___________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: nf90mpi_inq_varid          &
                       , nf90mpi_inquire_dimension  &
                       , nf90mpi_inquire_variable

      implicit none

      integer(MPI_OFFSET_KIND), dimension(:), allocatable :: var_shape

      integer     , intent(in) :: ncid
      character(*), intent(in) :: name

      integer, dimension(:), allocatable :: dimids
      integer                               i, ndims, varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , "inq_varid:var_shape "//name )
      call check( nf90mpi_inquire_variable( ncid, varid, ndims = ndims ) &
                , "inquire_variable(ndims):var_shape "//name )
      allocate( dimids(ndims), var_shape(ndims) )
      call check( nf90mpi_inquire_variable( ncid, varid, dimids = dimids ) &
                , "inquire_variable(dimids):var_shape "//name )
      do i = 1, ndims
        call check( nf90mpi_inquire_dimension( ncid, dimids(i), len = var_shape(i) ) &
                  , "inquire_dimension(len):var_shape "//name )
      end do


    end function var_shape
!
!___________________________________________________________________
!
    subroutine att_write_chr( ncid, varid, name, value )
!-------------------------------------------------------------------
!  Creates netcdf file.
!___________________________________________________________________
!
      use pnetcdf, only: NF_GLOBAL, nf90mpi_put_att

      implicit none

      integer     , intent(in) :: ncid, varid
      character(*), intent(in) :: name, value

      integer status


      if ( varid == -1 ) then
        status = nf90mpi_put_att( ncid, NF_GLOBAL     &
                                , name, trim(value) )
      else
        status = nf90mpi_put_att( ncid, varid         &
                                , name, trim(value) )
      end if
      call check( status, 'nf_put_att_chr: '//name)


    end ! function att_write_chr
!
!___________________________________________________________________
!
    subroutine att_write_int( ncid, varid, name, value )
!-------------------------------------------------------------------
!  Creates netcdf file.
!___________________________________________________________________
!
      use pnetcdf, only: NF_GLOBAL, nf90mpi_put_att

      implicit none

      integer     , intent(in) :: ncid, varid, value
      character(*), intent(in) :: name

      integer status


      if ( varid == -1 ) then
        status = nf90mpi_put_att( ncid, NF_GLOBAL, name, value )
      else
        status = nf90mpi_put_att( ncid, varid, name, value )
      end if
      call check( status, 'nf_put_att_int: '//name)


    end ! function att_write_int
!
!___________________________________________________________________
!
    subroutine att_write_flt( ncid, varid, name, value )
!-------------------------------------------------------------------
!  Creates netcdf file.
!___________________________________________________________________
!
      use pnetcdf, only: NF_GLOBAL, nf90mpi_put_att

      implicit none

      integer     , intent(in) :: ncid, varid
      character(*), intent(in) :: name
      real(4)     , intent(in) :: value

      integer status


      if ( varid == -1 ) then
        status = nf90mpi_put_att( ncid, NF_GLOBAL, name, value )
      else
        status = nf90mpi_put_att( ncid, varid, name, value )
      end if
      call check( status, 'nf_put_att_flt: '//name)


    end ! function att_write_flt
!
!___________________________________________________________________
!
    integer function file_create( path )
!-------------------------------------------------------------------
!  Creates netcdf file.
!___________________________________________________________________
!
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INFO_NULL
      use pnetcdf    , only: nf90mpi_create               &
                           , NF_64BIT_OFFSET, NF_CLOBBER

      implicit none

      character(*), intent(in) :: path


      call check( nf90mpi_create( POM_COMM, trim( path )        &
                                , NF_CLOBBER+NF_64BIT_OFFSET    &
                                , MPI_INFO_NULL, file_create )  &
                , 'nf_create: '//path )


    end ! function file_create
!
!___________________________________________________________________
!
    integer function file_open( path, mode )
!-------------------------------------------------------------------
!  Opens netcdf file for reading.
!___________________________________________________________________
!
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INFO_NULL
      use pnetcdf    , only: nf90mpi_open, NF90_NOERR

      implicit none

      character(*), intent(in) :: path
      integer     , intent(in) :: mode

      integer status


      file_open = 0

      status = nf90mpi_open( POM_COMM, trim(path), mode  &
                           , MPI_INFO_NULL, file_open )
      if ( is_error( status ) ) then
        call msg_print("", 2, "Failed to open `"//trim(path)//"`")
        file_open = status
      end if


    end ! function file_open
!
!___________________________________________________________________
!
    integer function file_close( ncid )
!-------------------------------------------------------------------
!  Opens netcdf file for reading.
!___________________________________________________________________
!
      use pnetcdf, only: nf90mpi_close

      implicit none

      integer, intent(in) :: ncid


      file_close = nf90mpi_close( ncid )


    end ! function file_close
!
!______________________________________________________________________
!
    subroutine file_end_definition( ncid )
!----------------------------------------------------------------------
!  Finalizes file definition.
!______________________________________________________________________
!
      use pnetcdf, only: nf90mpi_enddef

      implicit none

      integer, intent(in) :: ncid


      call check( nf90mpi_enddef( ncid )  &
                , "nf_enddef" )


    end ! function file_end_defintion
!
!______________________________________________________________________
!
    integer function var_read_0d( ncid, name, var, start )
!----------------------------------------------------------------------
!  Macro for reading scalar variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_get_att            &
                       , nf90mpi_get_var_all        &
                       , nf90mpi_inquire_attribute  &
                       , nf90mpi_inq_varid

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start

      integer  status, varid, vartype
      real(rk) add_offset, missing_value, scale_factor


      var_read_0d = nf90mpi_inq_varid( ncid, name, varid )

      call check( var_read_0d, "inq_varid `"//trim(name)//"`" )
      if ( var_read_0d /= NF90_NOERR ) return

      add_offset    = 0.
      scale_factor  = 1.
      missing_value = 9.96920996838687e+36

      status = nf90mpi_inquire_attribute( ncid, varid, "add_offset", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "add_offset", add_offset )
      end if

      status = nf90mpi_inquire_attribute( ncid, varid, "scale_factor", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "scale_factor", scale_factor )
      end if

      var_read_0d = nf90mpi_get_var_all( ncid, varid, var, start )

      if ( var_read_0d /= NF90_NOERR ) return

      status = nf90mpi_inquire_attribute( ncid, varid            &
                                        , "_FillValue", vartype )
      if ( status == NF90_NOERR ) then

        status = nf90mpi_get_att( ncid, varid, "_FillValue", missing_value )
        call check( status, 'Failed reading `_FillValue` of '//name )

        if ( var == missing_value ) then
          var = 0.
        else
          var = var*scale_factor + add_offset
        end if

        return

      end if

      var = var*scale_factor + add_offset


    end ! function var_read_0d
!______________________________________________________________________
!
    integer function var_read_1d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 1D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_get_att            &
                       , nf90mpi_get_var_all        &
                       , nf90mpi_inquire_attribute  &
                       , nf90mpi_inq_varid

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:)                    &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer  status, varid, vartype
      real(rk) add_offset, missing_value, scale_factor


      var_read_1d = nf90mpi_inq_varid( ncid, name, varid )

      call check( var_read_1d, "inq_varid `"//trim(name)//"`" )
      if ( var_read_1d /= NF90_NOERR ) return

      add_offset    = 0.
      scale_factor  = 1.
      missing_value = 9.96920996838687e+36

      status = nf90mpi_inquire_attribute( ncid, varid, "add_offset", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "add_offset", add_offset )
      end if

      status = nf90mpi_inquire_attribute( ncid, varid, "scale_factor", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "scale_factor", scale_factor )
      end if

      var_read_1d = nf90mpi_get_var_all( ncid, varid, var, start, stride )

      if ( var_read_1d /= NF90_NOERR ) return

      status = nf90mpi_inquire_attribute( ncid, varid            &
                                        , "_FillValue", vartype )
      if ( status == NF90_NOERR ) then

        status = nf90mpi_get_att( ncid, varid, "_FillValue", missing_value )
        call check( status, 'Failed reading `_FillValue` of '//name )

        where ( var == missing_value )
          var = 0.
        elsewhere
          var = var*scale_factor + add_offset
        end where

        return

      end if

      var = var*scale_factor + add_offset


    end ! function var_read_1d
!______________________________________________________________________
!
    integer function var_read_2d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 2D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_get_att            &
                       , nf90mpi_get_var_all        &
                       , nf90mpi_inquire_attribute  &
                       , nf90mpi_inq_varid

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:,:)                  &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer  status, varid, vartype
      real(rk) add_offset, missing_value, scale_factor


      var_read_2d = nf90mpi_inq_varid( ncid, name, varid )

      call check( var_read_2d, "inq_varid `"//trim(name)//"`" )
      if ( var_read_2d /= NF90_NOERR ) return

      add_offset    = 0.
      scale_factor  = 1.
      missing_value = 9.96920996838687e+36

      status = nf90mpi_inquire_attribute( ncid, varid, "add_offset", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "add_offset", add_offset )
      end if

      status = nf90mpi_inquire_attribute( ncid, varid, "scale_factor", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "scale_factor", scale_factor )
      end if

      var_read_2d = nf90mpi_get_var_all( ncid, varid, var, start, stride )

      if ( var_read_2d /= NF90_NOERR ) return

      status = nf90mpi_inquire_attribute( ncid, varid            &
                                        , "_FillValue", vartype )
      if ( status == NF90_NOERR ) then

        status = nf90mpi_get_att( ncid, varid, "_FillValue", missing_value )
        call check( status, 'Failed reading `_FillValue` of '//name )

        where ( var == missing_value )
          var = 0.
        elsewhere
          var = var*scale_factor + add_offset
        end where

        return

      end if

      var = var*scale_factor + add_offset


    end ! function var_read_2d
!______________________________________________________________________
!
    integer function var_read_3d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 3D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_get_att            &
                       , nf90mpi_get_var_all        &
                       , nf90mpi_inquire_attribute  &
                       , nf90mpi_inq_varid

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:,:,:)                &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer  status, varid, vartype
      real(rk) add_offset, missing_value, scale_factor


      var_read_3d = nf90mpi_inq_varid( ncid, name, varid )

      call check( var_read_3d, "inq_varid `"//trim(name)//"`" )
      if ( var_read_3d /= NF90_NOERR ) return

      add_offset    = 0.
      scale_factor  = 1.
      missing_value = 9.96920996838687e+36

      status = nf90mpi_inquire_attribute( ncid, varid, "add_offset", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "add_offset", add_offset )
      end if

      status = nf90mpi_inquire_attribute( ncid, varid, "scale_factor", vartype )
      if ( status == NF90_NOERR ) then
        status = nf90mpi_get_att( ncid, varid, "scale_factor", scale_factor )
      end if

      var_read_3d = nf90mpi_get_var_all( ncid, varid, var, start, stride )

      if ( var_read_3d /= NF90_NOERR ) return

      status = nf90mpi_inquire_attribute( ncid, varid            &
                                        , "_FillValue", vartype )
      if ( status == NF90_NOERR ) then

        status = nf90mpi_get_att( ncid, varid, "_FillValue", missing_value )
        call check( status, 'Failed reading `_FillValue` of '//name )

        where ( var == missing_value )
          var = 0.
        elsewhere
          var = var*scale_factor + add_offset
        end where

        return

      end if

      var = var*scale_factor + add_offset


    end ! function var_read_3d
!______________________________________________________________________
!
    subroutine var_write_0d( ncid, name, var, start )
!----------------------------------------------------------------------
!  Macro for writing scalar variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: nf90mpi_inq_varid          &
                       , nf90mpi_put_var_all

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start

      integer varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , 'nf_inq_varid: '//trim(name) )

      call check( nf90mpi_put_var_all( ncid, varid, var, start )  &
                    , "nf_put_var: "//trim(name) )


    end ! function var_write_0d
!______________________________________________________________________
!
    subroutine var_write_1d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 1D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_inq_varid          &
                       , nf90mpi_put_var_all

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:)                    &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , 'nf_inq_varid: '//trim(name) )

      call check( nf90mpi_put_var_all( ncid, varid, var, start, stride )  &
                    , "nf_put_var: "//trim(name) )


    end ! function var_write_1d
!______________________________________________________________________
!
    subroutine var_write_2d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 2D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_inq_varid          &
                       , nf90mpi_put_var_all

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:,:)                  &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , 'nf_inq_varid: '//trim(name) )

      call check( nf90mpi_put_var_all( ncid, varid, var, start, stride )  &
                , "nf_put_var: "//trim(name) )


    end ! function var_write_2d
!______________________________________________________________________
!
    subroutine var_write_3d( ncid, name, var, start, stride )
!----------------------------------------------------------------------
!  Macro for reading 3D variable.
!______________________________________________________________________
!
      use mpi    , only: MPI_OFFSET_KIND
      use pnetcdf, only: NF90_NOERR                 &
                       , nf90mpi_inq_varid          &
                       , nf90mpi_put_var_all

      implicit none

      integer     , intent(in   ) :: ncid
      character(*), intent(in   ) :: name
      real(rk)    , dimension(:,:,:)                &
                  , intent(inout) :: var
      integer(MPI_OFFSET_KIND)                      &
                  , dimension(:)                    &
                  , intent(in   ) :: start, stride

      integer varid


      call check( nf90mpi_inq_varid( ncid, name, varid )  &
                , 'nf_inq_varid: '//trim(name) )

      call check( nf90mpi_put_var_all( ncid, varid, var, start, stride )  &
                , "nf_put_var: "//trim(name) )


    end ! function var_write_3d
!
!______________________________________________________________________
!
    pure logical function is_error( status )

      use pnetcdf, only: NF90_NOERR

      implicit none

      integer, intent(in) :: status


      is_error = ( status /= NF90_NOERR )


    end ! function is_error
!
!______________________________________________________________________
!
    subroutine check(status, routine)
!----------------------------------------------------------------------
!  Checks for NetCDF I/O error and exits with an error message if hits.
!______________________________________________________________________
!
      use glob_domain, only: error_status, my_task
      use pnetcdf    , only: nf90mpi_strerror, NF_NOERR

      implicit none

      integer         , intent(in) :: status
      character(len=*), intent(in) :: routine


      if ( status /= NF_NOERR ) then
        !error_status = 1 ! Suppress for now... TODO: Do we even need this redundant error check?..
        print '(/a,i4,a,a)', 'IO error @ ', my_task, ': ', routine
        print *, nf90mpi_strerror(status)
      end if


    end ! subroutine check
!

end module io