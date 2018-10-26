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

module io

  use glob_const, only: PATH_LEN

  implicit none

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=PATH_LEN)  & ! Full paths (with filenames) to:
      bc_path         & !  boundary condition file
  , clim_path         & !  climatology file
  , grid_path         & !  grid file
  ,   ic_path         & !  ic file
  , mean_path         & !  mean ts file
  ,  out_path           ! Full path to output directory
!______________________________________________________________________
!  If `bc_path`, `bulk_path` or `flux_path` end with a `.` (dot)
! then these files will be read as <path>YYYY.<FORMAT_EXT>.
!  If any path ends with `/` then the default filename will be
! appended to the path.

  public :: check            , initialize_io            &
          , read_grid_pnetcdf, read_initial , out_init  &
          , grid_path

  private

  character(len=10), parameter :: FORMAT_EXT = "nc" ! Output format extension


  contains

!______________________________________________________________________
!
    subroutine initialize_io( config_nml )
!----------------------------------------------------------------------
!  Initialize IO module.
!______________________________________________________________________

      use glob_domain, only: is_master

      implicit none

      character(len=*), intent(in) :: config_nml

      namelist/input_files_nml/                                        &
          bc_path, clim_path, grid_path, ic_path, mean_path,  out_path

      integer pos

! Initialize with default values
      bc_path   = "in/bc/"
      clim_path = "in/clim/"
      grid_path = "in/grid/"
      ic_path   = "in/init/"
      mean_path = "in/clim/"
      out_path  = "out/"

! Override configuration
      open ( 73, file = config_nml, status = 'old' )
      read ( 73, nml = input_files_nml )
      close( 73 )

! Manage inmput_files values
      pos = len(trim(bc_path))
      if ( bc_path(pos:pos) == "/" ) then
        bc_path = trim(bc_path)//"bc."//FORMAT_EXT
      end if

      pos = len(trim(clim_path))
      if ( clim_path(pos:pos) == "/" ) then
        clim_path = trim(clim_path)//"clim."//FORMAT_EXT
      end if

      pos = len(trim(grid_path))
      if ( grid_path(pos:pos) == "/" ) then
        grid_path = trim(grid_path)//"grid."//FORMAT_EXT
      end if

      pos = len(trim(ic_path))
      if ( ic_path(pos:pos) == "/" ) then
        ic_path = trim(ic_path)//"init."//FORMAT_EXT
      end if

      pos = len(trim(mean_path))
      if ( mean_path(pos:pos) == "/" ) then
        mean_path = trim(mean_path)//"mean."//FORMAT_EXT
      end if

      pos = len(trim(out_path))
      if ( out_path(pos:pos) /= "/" ) then
        out_path = trim(out_path)//"/"
      end if

! Print config
      if ( is_master ) then
        print *, "Grid          : ", trim(grid_path)
        print *, "Climatology   : ", trim(clim_path)
        print *, "Mean TS       : ", trim(mean_path)
        print *, "Initial cond. : ", trim(ic_path)
        print *, "Boundary cond.: ", trim(bc_path)
        print *, "--------------|----------------------------"
        print *, "Output to     : ", trim(out_path)
        call msg_print("CORE I/O MODULE INITIALIZED", 1, "")
      end if

    end subroutine ! initialize_io
!______________________________________________________________________
!
    subroutine read_initial
!----------------------------------------------------------------------
!  Reads initial fields (T,S,u,v,el).
!______________________________________________________________________

      use model_run , only: dtime
      use glob_ocean, only: elb, rho, rmean, sb, sclim, tb, tclim

      implicit none

      logical fexist


! Check if the climatology file exists and read it.
      inquire ( file = trim(ic_path), exist = fexist )
      call msg_print("", 6, "Read initial conditions:")
      if ( fexist ) then
! TODO: Make it possible to read separate initial file, not just clim.
        call read_initial_ts_pnetcdf( tb, sb, elb, dtime%month )
      else
        call msg_print("", 2, "FAILED...")
        tb = 15.
        sb = 33.
      end if

! Calculate initial density.
      call dens(sb,tb,rho)

! Read z-horizontally averaged TS for background density field.
      inquire ( file = trim(mean_path), exist = fexist )
      call msg_print("", 6, "Read background TS:")
      if ( fexist ) then
        call read_mean_ts_pnetcdf( tclim, sclim, dtime%month )
      else
        call msg_print("", 2, "FAILED...")
        tclim = 0.
        sclim = 0.
      end if

! Calculate background density.
!  Should use T and S horizontally averaged in z-coordinates and only
! after that interpolated to sigma-coordinates.
      call dens(sclim,tclim,rmean)

! Read climatology.
      inquire ( file = trim(clim_path), exist = fexist )
      call msg_print("", 6, "Read climatology:")
      if ( fexist ) then
        call read_clim_ts_pnetcdf( tclim, sclim, dtime%month )
      else
        call msg_print("", 2, "FAILED...")
        tclim = 15.
        sclim = 33.
      end if


    end subroutine ! read_initial
!______________________________________________________________________
!
    subroutine out_init( out_file )
!----------------------------------------------------------------------
!  Wrapper for output procedure.
!______________________________________________________________________

      implicit none

      character(len=*), intent(in) :: out_file

      character(len=5), parameter :: pfx = "init."


      call write_output_init_pnetcdf(      trim(out_path)    &
                                    //          pfx          &
                                    //     trim(out_file)    &
                                    //"."//trim(FORMAT_EXT) )

    end subroutine ! out_init


!______________________________________________________________________
!
    integer function def_var_pnetcdf(                                &
                                ncid   , name     , nvdims , vdims   &
                              , vartype, long_name, units  , nofill  &
                              , fillval, coords   , lcoords )
!----------------------------------------------------------------------
!  Defines variable in NetCDF file.
!______________________________________________________________________

      use pnetcdf, only: nf90mpi_def_var, nf90mpi_def_var_fill  &
                       , nf90mpi_put_att                        &
                       , NF_NOERR

      integer                  &
      , intent(in) :: ncid     & ! ID of NetCDF file (in define mode)
                    , nvdims   & ! Number of dimensions
                    , vartype  & ! Variable's type
                    , nofill     ! Do not write missing values
                                 !  (useful for static masks)

      integer                       &
      , dimension(nvdims)           &
      , intent(in)        :: vdims    ! Variable's dimensions array

      character(len=*)               &
      , intent(in)     :: name       & ! Variable's name (ID)
                        , long_name  & ! Variable's human-readable name
                        , units      & ! Variable's units
                        , coords       ! Variable's coordinates (opt.)

      logical, intent(in) :: lcoords  ! Write coordinates or not
      real   , intent(in) :: fillval  ! Missing value

      integer varid ! Temporary variable

! define variable
      call check( nf90mpi_def_var                          &
                    ( ncid, name, vartype, vdims, varid )  &
                , 'nf_def_var @ '//name )
! define fill value
      if ( nofill > -1 ) then
        call check( nf90mpi_def_var_fill                &
                      ( ncid, varid, nofill, fillval )  &
                  , 'nf_def_var_fill @ '//name )
      end if
! define attributes
      call check( nf90mpi_put_att                                  &
                    ( ncid, varid, 'long_name', trim(long_name) )  &
                , 'nf_put_att : long_name @ '//name )

      call check( nf90mpi_put_att                          &
                    ( ncid, varid, 'units', trim(units) )  &
                , 'nf_put_att : units @ '//name )
! add coordinates attribute, if necessary
      if ( lcoords ) then
        call check( nf90mpi_put_att                                 &
                      ( ncid, varid, 'coordinates', trim(coords) )  &
                  , 'nf_put_att : coordinates @ '//name )
      end if

      def_var_pnetcdf = varid

      return

    end function ! def_var_pnetcdf
!______________________________________________________________________
!
      subroutine read_grid_pnetcdf( filepath )
!----------------------------------------------------------------------
!  Read grid data.
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain
      use glob_grid
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_1d  &
                         , get_var_real_2d  &
                         , get_var_real_3d

      character(len=120), intent(in) :: filepath

      integer                                                       &
             dx_varid,      dy_varid,  east_c_varid,  east_e_varid  &
      ,  east_u_varid,  east_v_varid,     fsm_varid,       h_varid  &
      , north_c_varid, north_e_varid, north_u_varid, north_v_varid  &
      ,     rot_varid,       z_varid,      zz_varid
      integer dimn, ncid

      integer(MPI_OFFSET_KIND) start(4),edge(4)


      start = 1
      edge  = 1

! open netcdf file
      if ( is_master ) print '(/''Reading file '',a)', trim(filepath)
      call check( nf90mpi_open                                       &
                  (POM_COMM,filepath,NF_NOWRITE,MPI_INFO_NULL,ncid)  &
                , "nf_open: "//trim(filepath) )

! get variables
      call check( nf90mpi_inq_varid(ncid,'z',z_varid)  &
                , 'nf_inq_varid: z' )
      call check( nfmpi_inq_varid(ncid,'zz',zz_varid)  &
                , 'nfmpi_inq_varid: zz' )
      call check( nfmpi_inq_varid(ncid,'dx',dx_varid)  &
                , 'nfmpi_inq_varid: dx' )
      call check( nfmpi_inq_varid(ncid,'dy',dy_varid)  &
                , 'nfmpi_inq_varid: dy' )
      call check( nfmpi_inq_varid(ncid,'lon_u',east_u_varid)  &
                , 'nfmpi_inq_varid: east_u' )
      call check( nfmpi_inq_varid(ncid,'lon_v',east_v_varid)  &
                , 'nfmpi_inq_varid: east_v' )
      call check( nfmpi_inq_varid(ncid,'lon_e',east_e_varid)  &
                , 'nfmpi_inq_varid: east_e' )
      call check( nfmpi_inq_varid(ncid,'lon_c',east_c_varid)  &
                , 'nfmpi_inq_varid: east_c' )
      call check( nfmpi_inq_varid(ncid,'lat_u',north_u_varid)  &
                , 'nfmpi_inq_varid: north_u' )
      call check( nfmpi_inq_varid(ncid,'lat_v',north_v_varid)  &
                , 'nfmpi_inq_varid: north_v' )
      call check( nfmpi_inq_varid(ncid,'lat_e',north_e_varid)  &
                , 'nfmpi_inq_varid: north_e' )
      call check( nfmpi_inq_varid(ncid,'lat_c',north_c_varid)  &
                , 'nfmpi_inq_varid: north_c' )
      call check( nfmpi_inq_varid(ncid,'rot',rot_varid)  &
                , 'nfmpi_inq_varid: rot' )
      call check( nfmpi_inq_varid(ncid,'h',h_varid)  &
                , 'nfmpi_inq_varid: h' )
      call check( nfmpi_inq_varid(ncid,'fsm',fsm_varid)  &
                , 'nfmpi_inq_varid: fsm' )
!      call check( nfmpi_inq_varid(ncid,'dum',dum_varid)  &
!                , 'nfmpi_inq_varid: dum' )
!      call check( nfmpi_inq_varid(ncid,'dvm',dvm_varid)  &
!                , 'nfmpi_inq_varid: dvm' )
! get data
      call check( nfmpi_inq_varndims(ncid, z_varid, dimn)  &
                , 'inq_vardimn: z' )

      if (dimn==1) then
        start(1)=1
        edge(1)=kb
        call check( get_var_real_1d(ncid, z_varid,start,edge, z)  &
                  , 'get_var_real: z' )
        call check( get_var_real_1d(ncid,zz_varid,start,edge,zz)  &
                  , 'get_var_real: zz' )
      else if (dimn==3) then
! If z and zz are 3-dimensional get just the (1,1) cell distribution
        start(1)=1
        start(2)=1
        start(3)=1
        edge(1)=1
        edge(2)=1
        edge(3)=kb
        call check( get_var_real_3d(ncid, z_varid,start,edge, z)  &
                  , 'get_var_real: z' )
        call check( get_var_real_3d(ncid,zz_varid,start,edge,zz)  &
                  , 'get_var_real: zz' )
      end if

      start(1)=i_global(1)
      start(2)=j_global(1)
      edge(1)=im
      edge(2)=jm
      call check( get_var_real_2d(ncid,dx_varid,start,edge,dx)  &
                , 'get_var_real: dx' )
      call check( get_var_real_2d(ncid,dy_varid,start,edge,dy)  &
                , 'get_var_real: dy' )
      call check( get_var_real_2d(ncid,east_u_varid,start,edge,east_u) &
                , 'get_var_real: east_u' )
      call check( get_var_real_2d(ncid,east_v_varid,start,edge,east_v) &
                , 'get_var_real: east_v' )
      call check( get_var_real_2d(ncid,east_e_varid,start,edge,east_e) &
                , 'get_var_real: east_e' )
      call check( get_var_real_2d(ncid,east_c_varid,start,edge,east_c) &
                , 'get_var_real: east_c' )
      call check( get_var_real_2d                          &
                  (ncid,north_u_varid,start,edge,north_u)  &
                , 'get_var_real: north_u' )
      call check( get_var_real_2d                          &
                  (ncid,north_v_varid,start,edge,north_v)  &
                , 'get_var_real: north_v' )
      call check( get_var_real_2d                          &
                  (ncid,north_e_varid,start,edge,north_e)  &
                , 'get_var_real: north_e' )
      call check( get_var_real_2d                          &
                  (ncid,north_c_varid,start,edge,north_c)  &
                , 'get_var_real: north_c' )
      call check( get_var_real_2d(ncid,rot_varid,start,edge,rot)  &
                , 'get_var_real: rot' )
      call check( get_var_real_2d(ncid,h_varid,start,edge,h)  &
                , 'get_var_real: h' )
      call check( get_var_real_2d(ncid,fsm_varid,start,edge,fsm)  &
                , 'get_var_real: fsm' )

! close file:
      call check( nf90mpi_close(ncid), 'nf_close: '//trim(filepath) )

      return
      end
!______________________________________________________________________
!
    subroutine write_output_init_pnetcdf( out_file )
!----------------------------------------------------------------------
!  Write initial state output file.
!______________________________________________________________________

      use air        , only: wssurf, wtsurf, wusurf, wvsurf
      use bry
      use config     , only: mode, title, use_air
      use glob_const , only: rk
      use glob_domain
      use glob_grid
      use glob_misc
      use glob_ocean
      use model_run  , only: time, time_start
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf    , only: nf90mpi_close     , nf90mpi_create     &
                           , nf90mpi_def_dim   , nf90mpi_enddef     &
                           , nf90mpi_get_att   , nf90mpi_inq_varid  &
                           , nf90mpi_open      , nf90mpi_put_att    &
                           , nfmpi_put_vara_all                     &
                           , nfmpi_put_vara_real_all                &
                           , NF_64BIT_OFFSET, NF_CLOBBER            &
                           , NF_FLOAT       , NF_GLOBAL             &
                           , NF_NOERR       , NF_NOWRITE            &
                           , NF_UNLIMITED   , NF_WRITE

      implicit none

      character(len=*), intent(in) :: out_file  ! Output filename

      integer time_dimid, x_dimid, y_dimid, z_dimid
      integer  aamfrz_varid,    dum_varid,    dvm_varid,     dx_varid &
            ,      dy_varid, east_c_varid, east_e_varid, east_u_varid &
            ,  east_v_varid,    ele_varid,    eln_varid,    els_varid &
            ,     elw_varid,    frz_varid,    fsm_varid,      h_varid &
            ,      kh_varid,     km_varid,north_c_varid,north_e_varid &
            , north_u_varid,north_v_varid,    rho_varid,    rot_varid &
            ,       s_varid,      t_varid,   time_varid, wssurf_varid &
            ,  wtsurf_varid, wusurf_varid, wvsurf_varid,      z_varid &
            ,      zz_varid

      integer, dimension(4)      :: vdims
      integer(MPI_OFFSET_KIND)                       &
             , dimension(4)      :: start(4),edge(4)

      integer ncid

      real(kind=4), dimension(im      ) :: out1x
      real(kind=4), dimension(   jm   ) :: out1y
      real(kind=4), dimension(      kb) :: out1z
      real(kind=4), dimension(im,jm   ) :: out2
      real(kind=4), dimension(im,jm,kb) :: out3

      out1x = 0.
      out1y = 0.
      out1z = 0.
      out2 = 0.
      out3 = 0.

      call msg_print("", 6, "Writing file `"//trim(out_file)//"`")

      call check( nf90mpi_create                  &
                    ( POM_COMM, trim(out_file)    &
                    , NF_CLOBBER+NF_64BIT_OFFSET  &
                    , MPI_INFO_NULL, ncid )       &
                , 'nf_create: '//out_file )

! define global attributes
      call check( nf90mpi_put_att                            &
                  ( ncid, NF_GLOBAL, 'title', trim(title) )  &
                , 'nf_put_att: title @ '//out_file )

      call check( nf90mpi_put_att                                  &
                  ( ncid, NF_GLOBAL, 'description', 'Init file' )  &
                , 'nf_put_att: description @ '//out_file )

! define dimensions
      call check( nf90mpi_def_dim                                   &
                  ( ncid, 'time', int(NF_UNLIMITED,8), time_dimid ) &
                , 'nf_def_dim: time @ '//out_file )
      call check( nf90mpi_def_dim                                   &
                  ( ncid,    'z', int(          kb,8),    z_dimid ) &
                , 'nf_def_dim: z @ '//out_file )
      call check( nf90mpi_def_dim                                   &
                  ( ncid,    'y', int(   jm_global,8),    y_dimid ) &
                , 'nf_def_dim: y @ '//out_file )
      call check( nf90mpi_def_dim                                   &
                  ( ncid,    'x', int(   im_global,8),    x_dimid ) &
                , 'nf_def_dim: x @ '//out_file )

! define variables and their attributes
      vdims(1) = time_dimid
      time_varid = def_var_pnetcdf(ncid,'time',1,vdims,NF_FLOAT  &
                              ,'time','days since '//time_start  &
                              ,-1,0.,' ',.false.)

! Skip baroclinic output if model run is barotropic
      if ( mode /= 2 ) then

        vdims(1) = z_dimid
        z_varid = def_var_pnetcdf(ncid,'sigma',1,vdims,NF_FLOAT  &
                            ,'sigma of cell face','sigma_level'  &
                            ,-1,0.,' ',.false.)
        call check( nf90mpi_put_att( ncid, z_varid               &
                                   , 'standard_name'             &
                                   , 'ocean_sigma_coordinate' )  &
                  , 'nf_put_att: stdname @ z @ '//out_file )
        call check( nf90mpi_put_att( ncid, z_varid                  &
                                   , 'formula_terms'                &
                                   , 'sigma: z eta: elb depth: h' ) &
                  , 'nf_put_att: formterms @ z @ '//out_file )

        zz_varid = def_var_pnetcdf(ncid,'zz',1,vdims,NF_FLOAT  &
                        ,'sigma of cell centre','sigma_level'  &
                        ,-1,0.,' ',.false.)
        call check( nf90mpi_put_att( ncid, zz_varid              &
                                   , 'standard_name'             &
                                   , 'ocean_sigma_coordinate' )  &
                  , 'nf_put_att: stdname @ zz @ '//out_file )
        call check( nf90mpi_put_att( ncid, zz_varid                  &
                                   , 'formula_terms'                 &
                                   , 'sigma: zz eta: elb depth: h' ) &
                  , 'nf_put_att: formterms @ zz @ '//out_file )

      end if


      vdims(1) = x_dimid
      els_varid = def_var_pnetcdf(ncid,'els',1,vdims,NF_FLOAT  &
                                ,'South elevation BC','meter'  &
                                ,-1,0.,'east_e',.true.)
      eln_varid = def_var_pnetcdf(ncid,'eln',1,vdims,NF_FLOAT  &
                                ,'North elevation BC','meter'  &
                                ,-1,0.,'east_e',.true.)
      vdims(1) = y_dimid
      ele_varid = def_var_pnetcdf(ncid,'ele',1,vdims,NF_FLOAT  &
                                 ,'East elevation BC','meter'  &
                                 ,-1,0.,'north_e',.true.)
      elw_varid = def_var_pnetcdf(ncid,'elw',1,vdims,NF_FLOAT  &
                                 ,'West elevation BC','meter'  &
                                 ,-1,0.,'north_e',.true.)


      vdims(1) = x_dimid
      vdims(2) = y_dimid
      dx_varid = def_var_pnetcdf(ncid,'dx',2,vdims,NF_FLOAT  &
                             ,'grid increment in x','meter'  &
                             ,-1,0.,'east_e north_e',.true.)
      dy_varid = def_var_pnetcdf(ncid,'dy',2,vdims,NF_FLOAT  &
                             ,'grid increment in y','meter'  &
                             ,-1,0.,'east_e north_e',.true.)
! TODO: Make units of easting and northing flexible (degree/meter)
      east_u_varid = def_var_pnetcdf(ncid,'east_u',2,vdims,NF_FLOAT  &
                                    ,'easting of u-points','degree'  &
                                    ,-1,0.,'east_u north_u',.true.)
      east_v_varid = def_var_pnetcdf(ncid,'east_v',2,vdims,NF_FLOAT  &
                                    ,'easting of v-points','degree'  &
                                    ,-1,0.,'east_v north_v',.true.)
      east_e_varid = def_var_pnetcdf(ncid,'east_e',2,vdims,NF_FLOAT  &
                            ,'easting of elevation points','degree'  &
                            ,-1,0.,'east_e north_e',.true.)
      east_c_varid = def_var_pnetcdf(ncid,'east_c',2,vdims,NF_FLOAT  &
                                ,'easting of cell corners','degree'  &
                                ,-1,0.,'east_c north_c',.true.)
      north_u_varid= def_var_pnetcdf(ncid,'north_u',2,vdims,NF_FLOAT &
                                    ,'northing of u-points','degree' &
                                    ,-1,0.,'east_u north_u',.true.)
      north_v_varid= def_var_pnetcdf(ncid,'north_v',2,vdims,NF_FLOAT &
                                    ,'northing of v-points','degree' &
                                    ,-1,0.,'east_v north_v',.true.)
      north_e_varid= def_var_pnetcdf(ncid,'north_e',2,vdims,NF_FLOAT &
                            ,'northing of elevation points','degree' &
                            ,-1,0.,'east_e north_e',.true.)
      north_c_varid= def_var_pnetcdf(ncid,'north_c',2,vdims,NF_FLOAT &
                                ,'northing of cell corners','degree' &
                                ,-1,0.,'east_c north_c',.true.)
      rot_varid = def_var_pnetcdf(ncid,'rot',2,vdims,NF_FLOAT  &
               ,'Rotation angle of x-axis wrt. east','degree'  &
               ,-1,0.,'east_e north_e',.true.)
      h_varid = def_var_pnetcdf(ncid,'h',2,vdims,NF_FLOAT  &
                       ,'undisturbed water depth','metre'  &
                       ,-1,0.,'east_e north_e',.true.)
      fsm_varid = def_var_pnetcdf(ncid,'fsm',2,vdims,NF_FLOAT  &
                         ,'free surface mask','dimensionless'  &
                         ,-1,0.,'east_e north_e',.true.)
      dum_varid = def_var_pnetcdf(ncid,'dum',2,vdims,NF_FLOAT  &
                           ,'u-velocity mask','dimensionless'  &
                           ,-1,0.,'east_u north_u',.true.)
      dvm_varid = def_var_pnetcdf(ncid,'dvm',2,vdims,NF_FLOAT  &
                           ,'v-velocity mask','dimensionless'  &
                           ,-1,0.,'east_v north_v',.true.)
      frz_varid = def_var_pnetcdf(ncid,'frz',2,vdims,NF_FLOAT  &
                                ,'frz coeffi','dimensionless'  &
                                ,-1,0.,'east_v north_v',.true.)
      aamfrz_varid = def_var_pnetcdf(ncid,'aamfrz',2,vdims,NF_FLOAT  &
                             ,'aamfrz coefficients','dimensionless'  &
                             ,-1,0.,'east_v north_v',.true.)

! If atmospheric forcing is enabled...
      if ( use_air ) then
        vdims(1) = x_dimid
        vdims(2) = y_dimid
        vdims(3) = time_dimid
        wusurf_varid=def_var_pnetcdf(ncid,'wusurf',3,vdims,NF_FLOAT  &
                                 ,'x-momentum flux','metre^2/sec^2'  &
                                 ,-1,0.,'east_u north_u',.true.)
        wvsurf_varid=def_var_pnetcdf(ncid,'wvsurf',3,vdims,NF_FLOAT  &
                                 ,'y-momentum flux','metre^2/sec^2'  &
                                 ,-1,0.,'east_v north_v',.true.)
        wtsurf_varid=def_var_pnetcdf(ncid,'wtsurf',3,vdims,NF_FLOAT  &
                                     ,'temperature flux','deg m/s'   &
                                     ,-1,0.,'east_e north_e',.true.)
        wssurf_varid=def_var_pnetcdf(ncid,'wssurf',3,vdims,NF_FLOAT  &
                                     ,'salinity flux','psu m/s'      &
                                     ,-1,0.,'east_e north_e',.true.)
      end if

! Skip baroclinic output...
      if ( mode /= 2 ) then

        vdims(1) = x_dimid
        vdims(2) = y_dimid
        vdims(3) = z_dimid
        vdims(4) = time_dimid
        t_varid = def_var_pnetcdf(ncid,'t',4,vdims,NF_FLOAT  &
                          ,'potential temperature','degC'    &
                          ,-1,0.,'east_e north_e zz',.true.)
        s_varid = def_var_pnetcdf(ncid,'s',4,vdims,NF_FLOAT  &
                          ,'salinity x rho / rhoref','PSS'   &
                          ,-1,0.,'east_e north_e zz',.true.)
        rho_varid = def_var_pnetcdf(ncid,'rho',4,vdims,NF_FLOAT  &
                       ,'(density-1000)/rhoref','dimensionless'  &
                       ,-1,0.,'east_e north_e zz',.true.)
        kh_varid = def_var_pnetcdf(ncid,'kh',4,vdims,NF_FLOAT  &
                         ,'vertical diffusivity','metre2/sec'  &
                         ,-1,0.,'east_e north_e zz',.true.)
        km_varid = def_var_pnetcdf(ncid,'km',4,vdims,NF_FLOAT  &
                           ,'vertical viscosity','metre2/sec'  &
                           ,-1,0.,'east_e north_e zz',.true.)

      end if

! end definitions
      call check( nf90mpi_enddef(ncid)                     &
                , 'nf_enddef: init file @ '//out_file )


! write data
      start(1) = 1
      edge(1) = 1
      out1z = real(time,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,time_varid,start,edge,out1z )     &
                , 'nf_put_vara_real:time @ '//out_file )

! Skip baroclinic output...
      if ( mode /= 2 ) then
        start(1) = 1
        edge(1) = kb
        out1z = real(z,4)
        call check( nfmpi_put_vara_real_all               &
                    ( ncid,z_varid,start,edge,out1z )      &
                  , 'nf_put_var_real: z @ '//out_file )
        out1z = real(zz,4)
        call check( nfmpi_put_vara_real_all                &
                    ( ncid,zz_varid,start,edge,out1z )      &
                  , 'nf_put_var_real: zz @ '//out_file )
      end if

! East
      if ( n_east == -1 ) then
        start(1) = j_global(1)
        edge(1) = jm
        out1y = real(el_bry%est(1,:),4)
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,ele_varid,start,edge,out1y )  &
                  , 'nf_put_var_real: ele @ '//out_file )
      else
        start(1) = 1
        edge(1) = 0
        out1y = 0.
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,ele_varid,start,edge,out1y )  &
                  , 'nf_put_var_real: ele (dummy) @ '//out_file )
      end if
! West
      if ( n_west == -1 ) then
        start(1) = j_global(1)
        edge(1) = jm
        out1y = real(el_bry%wst(1,:),4)
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,elw_varid,start,edge,out1y )  &
                  , 'nf_put_var_real: elw @ '//out_file )
      else
        start(1) = 1
        edge(1) = 0
        out1y = 0.
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,elw_varid,start,edge,out1y )  &
                  , 'nf_put_var_real: elw (dummy) @ '//out_file )
      end if
! South
      if ( n_south == -1 ) then
        start(1) = i_global(1)
        edge(1) = im
        out1x = real(el_bry%sth(:,1),4)
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,els_varid,start,edge,out1x )  &
                  , 'nf_put_var_real: els @ '//out_file )
      else
        start(1) = 1
        edge(1) = 0
        out1x = 0.
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,els_varid,start,edge,out1x )  &
                  , 'nf_put_var_real: els (dummy) @ '//out_file )
      end if
! North
      if ( n_north == -1 ) then
        start(1) = i_global(1)
        edge(1) = im
        out1x = real(el_bry%nth(:,1),4)
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,eln_varid,start,edge,out1x )  &
                  , 'nf_put_var_real: eln @ '//out_file )
      else
        start(1) = 1
        edge(1) = 0
        out1x = 0.
        call check( nfmpi_put_vara_real_all              &
                    ( ncid,eln_varid,start,edge,out1x )  &
                  , 'nf_put_var_real: eln @ (dummy)'//out_file )
      end if


      start(1) = i_global(1)
      start(2) = j_global(1)
      edge(1) = im
      edge(2) = jm
      out2 = real(dx,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,dx_varid,start,edge,out2 )      &
                , 'nf_put_var_real: dx @ '//out_file )
      out2 = real(dy,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,dy_varid,start,edge,out2 )      &
                , 'nf_put_var_real: dy @ '//out_file )
      out2 = real(east_u,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,east_u_varid,start,edge,out2 )      &
                , 'nf_put_var_real: east_u @ '//out_file )
      out2 = real(east_v,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,east_v_varid,start,edge,out2 )      &
                , 'nf_put_var_real: east_v @ '//out_file )
      out2 = real(east_e,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,east_e_varid,start,edge,out2 )      &
                , 'nf_put_var_real: east_e @ '//out_file )
      out2 = real(east_c,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,east_c_varid,start,edge,out2 )      &
                , 'nf_put_var_real: east_c @ '//out_file )
      out2 = real(north_u,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,north_u_varid,start,edge,out2 )      &
                , 'nf_put_var_real: north_u @ '//out_file )
      out2 = real(north_v,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,north_v_varid,start,edge,out2 )      &
                , 'nf_put_var_real: north_v @ '//out_file )
      out2 = real(north_e,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,north_e_varid,start,edge,out2 )      &
                , 'nf_put_var_real: north_e @ '//out_file )
      out2 = real(north_c,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,north_c_varid,start,edge,out2 )      &
                , 'nf_put_var_real: north_c @ '//out_file )
      out2 = real(rot,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,rot_varid,start,edge,out2 )      &
                , 'nf_put_var_real: rot @ '//out_file )
      out2 = real(h,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,h_varid,start,edge,out2 )      &
                , 'nf_put_var_real: h @ '//out_file )
      out2 = real(fsm,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,fsm_varid,start,edge,out2 )      &
                , 'nf_put_var_real: fsm @ '//out_file )
      out2 = real(dum,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,dum_varid,start,edge,out2 )      &
                , 'nf_put_var_real: dum @ '//out_file )
      out2 = real(dvm,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,dvm_varid,start,edge,out2 )      &
                , 'nf_put_var_real: dvm @ '//out_file )

      if ( USE_SPONGE ) then
        out2 = real(frz,4)
        call check( nfmpi_put_vara_real_all                 &
                    ( ncid,frz_varid,start,edge,out2 )      &
                  , 'nf_put_var_real: frz @ '//out_file )
        out2 = real(aamfrz,4)
        call check( nfmpi_put_vara_real_all                    &
                    ( ncid,aamfrz_varid,start,edge,out2 )      &
                  , 'nf_put_var_real: aamfrz @ '//out_file )
      end if

! Write only if atmospheric forcing is enabled.
      if ( use_air ) then

        start(1) = i_global(1)
        start(2) = j_global(1)
        start(3) = 1
        edge(1) = im
        edge(2) = jm
        edge(3) = 1
        out2 = real(wusurf,4)
        call check( nfmpi_put_vara_real_all                     &
                    ( ncid,wusurf_varid,start,edge,out2 )       &
                  , 'nf_put_vara_real: wusurf @ '//out_file )
        out2 = real(wvsurf,4)
        call check( nfmpi_put_vara_real_all                     &
                    ( ncid,wvsurf_varid,start,edge,out2 )       &
                  , 'nf_put_vara_real: wvsurf @ '//out_file )
        out2 = real(wtsurf,4)
        call check( nfmpi_put_vara_real_all                     &
                    ( ncid,wtsurf_varid,start,edge,out2 )       &
                  , 'nf_put_vara_real: wtsurf @ '//out_file )
        out2 = real(wssurf,4)
        call check( nfmpi_put_vara_real_all                     &
                    ( ncid,wssurf_varid,start,edge,out2 )       &
                  , 'nf_put_vara_real: wssurf @ '//out_file )

      end if

! Skip baroclinic output...
      if ( mode /= 2 ) then

        start(1) = i_global(1)
        start(2) = j_global(1)
        start(3) = 1
        start(4) = 1
        edge(1) = im
        edge(2) = jm
        edge(3) = kb
        edge(4) = 1
        out3 = real(tb,4)
        call check( nfmpi_put_vara_real_all               &
                    ( ncid,t_varid,start,edge,out3 )      &
                  , 'nf_put_vara_real: t @ '//out_file )
        out3 = real(sb,4)
        call check( nfmpi_put_vara_real_all               &
                    ( ncid,s_varid,start,edge,out3 )      &
                  , 'nf_put_vara_real: s @ '//out_file )
        out3 = real(rho,4)
        call check( nfmpi_put_vara_real_all                 &
                    ( ncid,rho_varid,start,edge,out3 )      &
                  , 'nf_put_vara_real: rho @ '//out_file )
        out3 = real(kh,4)
        call check( nfmpi_put_vara_real_all                &
                    ( ncid,kh_varid,start,edge,out3 )      &
                  , 'nf_put_vara_real: kh @ '//out_file )
        out3 = real(km,4)
        call check( nfmpi_put_vara_real_all                &
                    ( ncid,km_varid,start,edge,out3 )      &
                  , 'nf_put_vara_real: km @ '//out_file )

      end if

! close file
      call check( nf90mpi_close(ncid)                   &
                , 'nf_close: output:  @ '//out_file )

      return

    end subroutine
!______________________________________________________________________
!
    subroutine read_initial_ts_pnetcdf( temp, salt, ssh, n )
!----------------------------------------------------------------------
!  Reads initial temperature, salinity and elevation.
! TODO: Add u and v. Add check for single record to read if not climatology.
!______________________________________________________________________

      use glob_const , only: rk
      use glob_domain
!      use glob_grid  , only: fsm
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_2d &
                         , get_var_real_3d

      integer , intent(in)  :: n
      real(rk), intent(out) :: temp(im,jm,kb),salt(im,jm,kb),ssh(im,jm)

      integer                  ncid, status
      integer                  el_varid, sb_varid, tb_varid
      integer(MPI_OFFSET_KIND) start(4),edge(4)


      start = 1
      edge  = 1

!  Open netcdf file.
      if ( is_master )                                              &
           print '(''reading file `'',a,''`''/)', trim(ic_path)
      call check( nf90mpi_open(                        &
                    POM_COMM, ic_path, NF_NOWRITE  &
                  , MPI_INFO_NULL, ncid )              &
                , "nfmpi_open: "//ic_path )

!  Get variables. [ TODO: parameterize varnames ]
      call check( nf90mpi_inq_varid( ncid, 'tclim', tb_varid )  &
                , "nfmpi_inq_varid: t_init @ "//ic_path )
      call check( nf90mpi_inq_varid( ncid, 'sclim', sb_varid )  &
                , "nfmpi_inq_varid: s_init @ "//ic_path )

      status = nf90mpi_inq_varid( ncid, 'eclim', el_varid )
      if ( status == NF_NOERR ) then
! Get elevation.
        start(1) = i_global(1)
        start(2) = j_global(1)
        start(3) = n
        edge(1) = im
        edge(2) = jm
        edge(3) = 1
        call check( get_var_real_2d( ncid, el_varid, start,edge, ssh ) &
                  , "nfmpi_get_vara_real_all"//ic_path )
      else
        if ( status == -49 ) then
          if ( is_master ) then
            call msg_print("", 2                                       &
                        , "Missing elevation data - falling back to 0.")
          end if
          ssh = 0.
        else
          call check(status, 'nfmpi_inq_varid: el_init @ '//ic_path)
        end if
      end if

      start(1) = i_global(1)
      start(2) = j_global(1)
      start(3) = 1
      start(4) = n
      edge(1) = im
      edge(2) = jm
      edge(3) = kb
      edge(4) = 1
      call check( get_var_real_3d(ncid,tb_varid,start,edge,temp)  &
                , "nf_get_var : temp @ "//ic_path )
      call check( get_var_real_3d(ncid,sb_varid,start,edge,salt)  &
                , "nf_get_var : salt @ "//ic_path )

!! TODO: move out
!      do k=1,kb-1
!        where (fsm==0.)
!          temp(:,:,k) = 0.
!          salt(:,:,k) = 0.
!        end where
!      end do
!      temp(:,:,kb) = temp(:,:,kb-1)
!      salt(:,:,kb) = salt(:,:,kb-1)
!      elb = elb*fsm

!  Close file.
      call check( nf90mpi_close(ncid), "nf_close @ "//ic_path)


    end subroutine

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
          print '(/a,a)', 'IO error: ', routine
          print *, nf90mpi_strerror(status)
          stop
        end if
      end if


      return

    end


end module io
