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

  use glob_const, only: PATH_LEN, VAR_LEN

  implicit none

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=PATH_LEN)  & ! Full paths (with filenames) to:
    grid_path              & !  grid file
  ,   ic_path              & !  ic file
  ,  out_path                ! Full path to output directory
!______________________________________________________________________
!  If `bc_path`, `bulk_path` or `flux_path` end with a `.` (dot)
! then these files will be read as <path>YYYY.<FORMAT_EXT>.
!  If any path ends with `/` then the default filename will be
! appended to the path.

!----------------------------------------------------------------------
! Varnames configuration
!----------------------------------------------------------------------
  character(len=VAR_LEN)  & ! Input variables names:
      dx_name         & !  grid's x-spacing (m)
  ,   dy_name         & !  grid's y-spacing (m)
  ,   ec_name         & !  grid's easting for cell corners (deg)
  ,   ee_name         & !  grid's easting for cell centers (deg)
  ,   eu_name         & !  grid's easting for u-points (deg)
  ,   ev_name         & !  grid's easting for v-points (deg)
  ,  fsm_name         & !  free surface mask ({0,1})
  ,    h_name         & !  bathymetry (m)
  ,   nc_name         & !  grid's northing for cell corners (deg)
  ,   ne_name         & !  grid's northing for cell centers (deg)
  ,   nu_name         & !  grid's northing for u-points (deg)
  ,   nv_name         & !  grid's northing for v-points (deg)
  ,  rot_name         & !  cell's angle (rad)
  ,    z_name         & !  vertical grid's cell faces positions (0..-1)
  ,   zz_name           !  vertical grid's cell centers positions (0..-1-eps)

!  TODO: make `io` module independent and move all irrelevant routines
! to their respective modules.
  public :: check        , file_open_nc     , file_close_nc  &
          , initialize_io, read_grid_pnetcdf, read_initial   &
          , out_init     , grid_path        , out_debug

  private

  character(len=10), parameter :: FORMAT_EXT = "nc" ! Output format extension


  contains

!______________________________________________________________________
!
    subroutine initialize_io( config_nml )
!----------------------------------------------------------------------
!  Initializes IO module.
!______________________________________________________________________
!
      use glob_domain, only: is_master

      implicit none

      character(len=*), intent(in) :: config_nml

      namelist/input_files_nml/                                       &
          grid_path, ic_path, out_path
      namelist/input_vars_nml/                                        &
           dx_name, dy_name, ec_name, ee_name, eu_name, ev_name       &
        , fsm_name,  h_name, nc_name, ne_name, nu_name, nv_name       &
        , rot_name,  z_name, zz_name

      integer number_of_modules_to_output, pos


! Initialize with default values
      grid_path = "in/grid/"
      ic_path   = "in/init/"
      out_path  = "out/"

      h_name  = "h"
      z_name  = "z"
      zz_name = "zz"
      dx_name = "dx"
      dy_name = "dy"
      ec_name = "lon_c"
      ee_name = "lon_e"
      eu_name = "lon_u"
      ev_name = "lon_v"
      nc_name = "lat_c"
      ne_name = "lat_e"
      nu_name = "lat_u"
      nv_name = "lat_v"
      fsm_name= "fsm"
      rot_name= "rot"

! Override configuration
      open ( 73, file = config_nml, status = 'old' )
      read ( 73, nml = input_files_nml )
      read ( 73, nml = input_vars_nml )
      close( 73 )

! Manage inmput_files values
      pos = len(trim(grid_path))
      if ( grid_path(pos:pos) == "/" ) then
        grid_path = trim(grid_path)//"grid."//FORMAT_EXT
      end if

      pos = len(trim(out_path))
      if ( out_path(pos:pos) /= "/" ) then
        out_path = trim(out_path)//"/"
      end if

! Print config
      if ( is_master ) then
        print *, "--------------|----------------------------"
        print *, " Grid         : ", trim(grid_path)
        print *, "--------------|----------------------------"
        print *, " Output to    : ", trim(out_path)
        print *, "--------------|----------------------------"
        call msg_print("CORE I/O MODULE INITIALIZED", 1, "")
      end if


    end ! subroutine initialize_io
!
!______________________________________________________________________
!
    subroutine read_initial
!----------------------------------------------------------------------
!  Reads initial fields (T,S,u,v,el).
!______________________________________________________________________
!
      use model_run  , only: dtime
      use glob_domain, only: kbm1
      use glob_grid  , only: dz, fsm
      use glob_ocean , only: elb, rho, sb, tb, uab, ub, vab, vb

      implicit none

      logical fexist
      integer k


! Check if the climatology file exists and read it.
      inquire ( file = trim(ic_path), exist = fexist )
      call msg_print("", 6, "Read initial conditions:")
      if ( fexist ) then
! TODO: Make it possible to read separate initial file, not just clim.
        call read_initial_ts_pnetcdf( tb, sb, ub, vb, elb, 1 ) !dtime%month )
!        call read_ts_z_pnetcdf( tb, sb, 40, dtime%month, ic_path )
      else
        call msg_print("", 2, "FAILED...")
        tb = 15.
        sb = 33.
      end if

      uab = 0.
      vab = 0.

      do k = 1,kbm1
        sb(:,:,k) = sb(:,:,k) * fsm
        tb(:,:,k) = tb(:,:,k) * fsm
        uab(:,:)  = uab(:,:) + ub(:,:,k)*dz(k)
        vab(:,:)  = vab(:,:) + vb(:,:,k)*dz(k)
      end do

! Calculate initial density.
      call dens(sb,tb,rho)

!      uab = .6*dum
!      if (n_south==-1) then
!        dvm(:,1)=0.
!        dum(:,1)=0.
!      end if

!      if(n_east.eq.-1) then
!        do j=1,jm
!          do i=1,im
!           aam2d(i,j)=100.e0*(1.+10.e0*exp(-(im-i)/10.e0))
!          end do
!        end do
!      end if


    end ! subroutine read_initial
!
!______________________________________________________________________
!
    subroutine out_init( out_file )
!----------------------------------------------------------------------
!  Wrapper for output procedure.
!______________________________________________________________________
!
      implicit none

      character(len=*), intent(in) :: out_file

      character(len=5), parameter :: pfx = "init."


      call write_output_init_pnetcdf(      trim(out_path)    &
                                    //          pfx          &
                                    //     trim(out_file)    &
                                    //"."//trim(FORMAT_EXT) )

    end ! subroutine out_init
!
!______________________________________________________________________
!
    subroutine out_debug( out_file )
!----------------------------------------------------------------------
!  Wrapper for debug output procedure.
!______________________________________________________________________
!
      implicit none

      character(len=*), intent(in) :: out_file

      character(len=4), parameter :: pfx = "dbg."


      call write_debug_pnetcdf(      trim(out_path)    &
                              //          pfx          &
                              //     trim(out_file)    &
                              //"."//trim(FORMAT_EXT) )


    end ! subroutine out_debug
!
!______________________________________________________________________
!
    integer function def_var_pnetcdf(                                &
                                ncid   , name     , nvdims , vdims   &
                              , vartype, long_name, units  , nofill  &
                              , fillval, coords   , lcoords )
!----------------------------------------------------------------------
!  Defines variable in NetCDF file.
!______________________________________________________________________
!
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


    end ! function def_var_pnetcdf
!
!______________________________________________________________________
!
    subroutine read_grid_pnetcdf( filepath )
!----------------------------------------------------------------------
!  Read grid data.
!______________________________________________________________________
!
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
      call check( nf90mpi_inq_varid(ncid,z_name,z_varid)  &
                , 'nf_inq_varid: z' )
      call check( nfmpi_inq_varid(ncid,zz_name,zz_varid)  &
                , 'nfmpi_inq_varid: zz' )
      call check( nfmpi_inq_varid(ncid,dx_name,dx_varid)  &
                , 'nfmpi_inq_varid: dx' )
      call check( nfmpi_inq_varid(ncid,dy_name,dy_varid)  &
                , 'nfmpi_inq_varid: dy' )
      call check( nfmpi_inq_varid(ncid,eu_name,east_u_varid)  &
                , 'nfmpi_inq_varid: east_u' )
      call check( nfmpi_inq_varid(ncid,ev_name,east_v_varid)  &
                , 'nfmpi_inq_varid: east_v' )
      call check( nfmpi_inq_varid(ncid,ee_name,east_e_varid)  &
                , 'nfmpi_inq_varid: east_e' )
      call check( nfmpi_inq_varid(ncid,ec_name,east_c_varid)  &
                , 'nfmpi_inq_varid: east_c' )
      call check( nfmpi_inq_varid(ncid,nu_name,north_u_varid)  &
                , 'nfmpi_inq_varid: north_u' )
      call check( nfmpi_inq_varid(ncid,nv_name,north_v_varid)  &
                , 'nfmpi_inq_varid: north_v' )
      call check( nfmpi_inq_varid(ncid,ne_name,north_e_varid)  &
                , 'nfmpi_inq_varid: north_e' )
      call check( nfmpi_inq_varid(ncid,nc_name,north_c_varid)  &
                , 'nfmpi_inq_varid: north_c' )
      call check( nfmpi_inq_varid(ncid,rot_name,rot_varid)  &
                , 'nfmpi_inq_varid: rot' )
      call check( nfmpi_inq_varid(ncid,h_name,h_varid)  &
                , 'nfmpi_inq_varid: h' )
      call check( nfmpi_inq_varid(ncid,fsm_name,fsm_varid)  &
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


    end ! subroutine read_grid_pnetcdf
!
!______________________________________________________________________
!
    subroutine write_output_init_pnetcdf( out_file )
!----------------------------------------------------------------------
!  Write initial state output file.
!______________________________________________________________________
!
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
        print *, "!!!", shape(el_bry%est), shape(out1y)
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


    end ! subroutine write_output_init_pnetcdf
!
!______________________________________________________________________
!
    subroutine write_debug_pnetcdf( out_file )
!----------------------------------------------------------------------
!  Write initial state output file.
!______________________________________________________________________
!
!      use air        , only: wssurf, wtsurf, wusurf, wvsurf
      use bry        , only: el_bry
      use config     , only: mode, title, use_air
      use glob_const , only: rk
      use glob_domain
      use glob_grid
      use glob_misc
      use glob_ocean
      use tide
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

      integer time_dimid, x_dimid, y_dimid, z_dimid, cons_dimid
      integer  aam2d_varid,    aam_varid,  advua_varid, advva_varid  &
            ,   advx_varid,   advy_varid,  adx2d_varid, ady2d_varid  &
            ,  drhox_varid,  drhoy_varid,  drx2d_varid, dry2d_varid  &
            ,e_atmos_varid,    egb_varid,    egf_varid,    el_varid  &
            ,    elb_varid,    ele_varid,    eln_varid,   elf_varid  &
            ,    els_varid,    elw_varid,     et_varid,   etb_varid  &
            ,    etf_varid, fluxua_varid, fluxva_varid,    kh_varid  &
            ,     km_varid,     kq_varid,      l_varid,    q2_varid  &
            ,    q2b_varid,    q2l_varid,   q2lb_varid,   rho_varid  &
            ,      s_varid,     sb_varid,      t_varid,  time_varid  &
            ,     tb_varid,      u_varid,     ua_varid,   uab_varid  &
            ,    uaf_varid,     ub_varid,     uf_varid,     v_varid  &
            ,     va_varid,     vb_varid,    vab_varid,   vaf_varid  &
            ,     vf_varid,      w_varid,     wr_varid,     z_varid  &
            ,     zz_varid, el_amp_varid

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

      call msg_print("", 6, "Writing debug `"//trim(out_file)//"`")

      call check( nf90mpi_create                  &
                    ( POM_COMM, trim(out_file)    &
                    , NF_CLOBBER+NF_64BIT_OFFSET  &
                    , MPI_INFO_NULL, ncid )       &
                , 'nf_create: '//out_file )

! define global attributes
      call check( nf90mpi_put_att                            &
                  ( ncid, NF_GLOBAL, 'title', trim(title) )  &
                , 'nf_put_att: title @ '//out_file )

      call check( nf90mpi_put_att                                   &
                  ( ncid, NF_GLOBAL, 'description', 'Debug file' )  &
                , 'nf_put_att: description @ '//out_file )

! define dimensions
      call check( nf90mpi_def_dim                                   &
                  ( ncid, 'time', int(NF_UNLIMITED,8), time_dimid ) &
                , 'nf_def_dim: time @ '//out_file )
      call check( nf90mpi_def_dim                                   &
                  ( ncid,  'nct', int(       ncons,8), cons_dimid ) &
                , 'nf_def_dim: constituents @ '//out_file )
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
      aam2d_varid = def_var_pnetcdf(ncid,'aam2d',2,vdims,NF_FLOAT  &
                             ,'vertical average of aam',''  &
                             ,-1,0.,'east_e north_e',.true.)
      advua_varid = def_var_pnetcdf(ncid,'advua',2,vdims,NF_FLOAT      &
                   ,'sum of the 2nd, 3rd and 4th terms in eq (18)',''  &
                   ,-1,0.,'east_e north_e',.true.)
      advva_varid = def_var_pnetcdf(ncid,'advva',2,vdims,NF_FLOAT  &
                   ,'sum of the 2nd, 3rd and 4th terms in eq (19)',''  &
                   ,-1,0.,'east_u north_u',.true.)
      adx2d_varid = def_var_pnetcdf(ncid,'adx2d',2,vdims,NF_FLOAT   &
                                   ,'vertical integral of advx',''  &
                                   ,-1,0.,'east_v north_v',.true.)
      ady2d_varid = def_var_pnetcdf(ncid,'ady2d',2,vdims,NF_FLOAT   &
                                   ,'vertical integral of advy',''  &
                                   ,-1,0.,'east_e north_e',.true.)
      drx2d_varid = def_var_pnetcdf(ncid,'drx2d',2,vdims,NF_FLOAT    &
                                   ,'vertical integral of drhox',''  &
                                   ,-1,0.,'east_c north_c',.true.)
      dry2d_varid = def_var_pnetcdf(ncid,'dry2d',2,vdims,NF_FLOAT   &
                                   ,'vertical integral of drhoy','' &
                                   ,-1,0.,'east_u north_u',.true.)
      e_atmos_varid = def_var_pnetcdf(ncid,'e_atmos',2,vdims,NF_FLOAT &
                                     ,'atmospheric pressure',''       &
                                     ,-1,0.,'east_v north_v',.true.)
      egb_varid = def_var_pnetcdf(ncid,'egb',2,vdims,NF_FLOAT          &
         ,'surface elevation use for pressure gradient at time n-1','' &
         ,-1,0.,'east_e north_e',.true.)
      egf_varid = def_var_pnetcdf(ncid,'egf',2,vdims,NF_FLOAT &
         ,'surface elevation use for pressure gradient at time n+1','' &
         ,-1,0.,'east_c north_c',.true.)
      el_varid = def_var_pnetcdf(ncid,'el',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the external mode at time n',''  &
         ,-1,0.,'east_e north_e',.true.)
      elb_varid = def_var_pnetcdf(ncid,'elb',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the external mode at time n-1','' &
         ,-1,0.,'east_e north_e',.true.)
      elf_varid = def_var_pnetcdf(ncid,'elf',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the external mode at time n+1','' &
         ,-1,0.,'east_e north_e',.true.)
      et_varid = def_var_pnetcdf(ncid,'et',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the internal mode at time n',''  &
         ,-1,0.,'east_u north_u',.true.)
      etb_varid = def_var_pnetcdf(ncid,'etb',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the internal mode at time n-1','' &
         ,-1,0.,'east_v north_v',.true.)
      etf_varid = def_var_pnetcdf(ncid,'etf',2,vdims,NF_FLOAT  &
         ,'surface elevation used in the internal mode at time n+1','' &
         ,-1,0.,'east_v north_v',.true.)
      fluxua_varid = def_var_pnetcdf(ncid,'fluxua',2,vdims,NF_FLOAT  &
                             ,'fluxua','dimensionless'  &
                             ,-1,0.,'east_v north_v',.true.)
      fluxva_varid = def_var_pnetcdf(ncid,'fluxva',2,vdims,NF_FLOAT  &
                             ,'fluxva','dimensionless'  &
                             ,-1,0.,'east_v north_v',.true.)
      ua_varid = def_var_pnetcdf(ncid,'ua',2,vdims,NF_FLOAT         &
                                ,'vertical mean of u at time n',''  &
                                ,-1,0.,'east_v north_v',.true.)
      uab_varid = def_var_pnetcdf(ncid,'uab',2,vdims,NF_FLOAT          &
                                 ,'vertical mean of u at time n-1',''  &
                                 ,-1,0.,'east_v north_v',.true.)
      uaf_varid = def_var_pnetcdf(ncid,'uaf',2,vdims,NF_FLOAT          &
                                 ,'vertical mean of u at time n+1',''  &
                                 ,-1,0.,'east_v north_v',.true.)
      va_varid = def_var_pnetcdf(ncid,'va',2,vdims,NF_FLOAT         &
                                ,'vertical mean of v at time n',''  &
                                ,-1,0.,'east_v north_v',.true.)
      vab_varid = def_var_pnetcdf(ncid,'vab',2,vdims,NF_FLOAT          &
                                 ,'vertical mean of v at time n-1',''  &
                                 ,-1,0.,'east_v north_v',.true.)
      vaf_varid = def_var_pnetcdf(ncid,'vaf',2,vdims,NF_FLOAT          &
                                 ,'vertical mean of v at time n+1',''  &
                                 ,-1,0.,'east_v north_v',.true.)

      vdims(1) = x_dimid
      vdims(2) = y_dimid
      vdims(3) = cons_dimid
      el_amp_varid = def_var_pnetcdf(ncid,'el_amp',3,vdims,NF_FLOAT    &
                                 ,'elevation amplitude',''  &
                                 ,-1,0.,'east_e north_e nct',.true.)

      vdims(1) = x_dimid
      vdims(2) = y_dimid
      vdims(3) = z_dimid
      vdims(4) = time_dimid
      aam_varid = def_var_pnetcdf(ncid,'aam',4,vdims,NF_FLOAT   &
                          ,'horizontal kinematic viscosity',''  &
                          ,-1,0.,'east_e north_e zz',.true.)
      advx_varid = def_var_pnetcdf(ncid,'advx',4,vdims,NF_FLOAT     &
                  ,'x-horizontal advection and diffusion terms',''  &
                  ,-1,0.,'east_e north_e zz',.true.)
      advy_varid = def_var_pnetcdf(ncid,'advy',4,vdims,NF_FLOAT     &
                  ,'y-horizontal advection and diffusion terms',''  &
                  ,-1,0.,'east_e north_e zz',.true.)
      drhox_varid = def_var_pnetcdf(ncid,'drhox',4,vdims,NF_FLOAT     &
               ,'x-component of the internal baroclinic pressure',''  &
               ,-1,0.,'east_e north_e zz',.true.)
      drhoy_varid = def_var_pnetcdf(ncid,'drhoy',4,vdims,NF_FLOAT     &
               ,'y-component of the internal baroclinic pressure',''  &
               ,-1,0.,'east_e north_e zz',.true.)
      kq_varid = def_var_pnetcdf(ncid,'kq',4,vdims,NF_FLOAT     &
                ,'vertical turbulent something',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      l_varid = def_var_pnetcdf(ncid,'l',4,vdims,NF_FLOAT     &
                ,'turbulence length scale',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      q2b_varid = def_var_pnetcdf(ncid,'q2b',4,vdims,NF_FLOAT     &
                ,'twice the turbulent kinetic energy at time n-1',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      q2_varid = def_var_pnetcdf(ncid,'q2',4,vdims,NF_FLOAT     &
                ,'twice the turbulent kinetic energy at time n',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      q2lb_varid = def_var_pnetcdf(ncid,'q2lb',4,vdims,NF_FLOAT     &
                ,'q2 x l at time n-1',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      q2l_varid = def_var_pnetcdf(ncid,'q2l',4,vdims,NF_FLOAT     &
                ,'q2 x l at time n',''  &
                ,-1,0.,'east_e north_e zz',.true.)
      t_varid = def_var_pnetcdf(ncid,'t',4,vdims,NF_FLOAT  &
                        ,'potential temperature','degC'    &
                        ,-1,0.,'east_e north_e zz',.true.)
      tb_varid = def_var_pnetcdf(ncid,'tb',4,vdims,NF_FLOAT    &
                        ,'potential temperature at n-1','degC' &
                        ,-1,0.,'east_e north_e zz',.true.)
      s_varid = def_var_pnetcdf(ncid,'s',4,vdims,NF_FLOAT  &
                        ,'salinity x rho / rhoref','PSS'   &
                        ,-1,0.,'east_e north_e zz',.true.)
      sb_varid = def_var_pnetcdf(ncid,'sb',4,vdims,NF_FLOAT      &
                        ,'salinity x rho / rhoref at n-1','PSS'  &
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
      u_varid = def_var_pnetcdf(ncid,'u',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      ub_varid = def_var_pnetcdf(ncid,'ub',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n-1',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      uf_varid = def_var_pnetcdf(ncid,'uf',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n-1',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      v_varid = def_var_pnetcdf(ncid,'v',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      vb_varid = def_var_pnetcdf(ncid,'vb',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n-1',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      vf_varid = def_var_pnetcdf(ncid,'vf',4,vdims,NF_FLOAT         &
                        ,'horizontal velocity in x at time n-1',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      w_varid = def_var_pnetcdf(ncid,'w',4,vdims,NF_FLOAT         &
                        ,'sigma coordinate vertical velocity',''  &
                        ,-1,0.,'east_e north_e zz',.true.)
      wr_varid = def_var_pnetcdf(ncid,'wr',4,vdims,NF_FLOAT          &
                        ,'real (z coordinate) vertical velocity',''  &
                        ,-1,0.,'east_e north_e zz',.true.)

! end definitions
      call check( nf90mpi_enddef(ncid)                     &
                , 'nf_enddef: debug file @ '//out_file )


! write data
      start(1) = 1
      edge(1) = 1
      out1z = real(time,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,time_varid,start,edge,out1z )     &
                , 'nf_put_vara_real:time @ '//out_file )

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

      call exchange2d_mpi(aam2d,im,jm)
      call exchange2d_mpi(advua,im,jm)
      call exchange2d_mpi(advva,im,jm)
      call exchange2d_mpi(adx2d,im,jm)
      call exchange2d_mpi(ady2d,im,jm)
      call exchange2d_mpi(drx2d,im,jm)
      call exchange2d_mpi(dry2d,im,jm)
      call exchange2d_mpi(egb,im,jm)
      call exchange2d_mpi(egf,im,jm)
      call exchange2d_mpi(el,im,jm)
      call exchange2d_mpi(elb,im,jm)
      call exchange2d_mpi(elf,im,jm)
      call exchange2d_mpi(et,im,jm)
      call exchange2d_mpi(etb,im,jm)
      call exchange2d_mpi(etb,im,jm)
      call exchange2d_mpi(fluxua,im,jm)
      call exchange2d_mpi(fluxva,im,jm)
      call exchange2d_mpi(ua,im,jm)
      call exchange2d_mpi(uab,im,jm)
      call exchange2d_mpi(uaf,im,jm)
      call exchange2d_mpi(va,im,jm)
      call exchange2d_mpi(vab,im,jm)
      call exchange2d_mpi(vaf,im,jm)

      call exchange3d_mpi(aam,im,jm,kb)
      call exchange3d_mpi(advx,im,jm,kb)
      call exchange3d_mpi(advy,im,jm,kb)
      call exchange3d_mpi(drhox,im,jm,kb)
      call exchange3d_mpi(drhoy,im,jm,kb)
      call exchange3d_mpi(kq,im,jm,kb)
      call exchange3d_mpi(l,im,jm,kb)
      call exchange3d_mpi(q2b,im,jm,kb)
      call exchange3d_mpi(q2,im,jm,kb)
      call exchange3d_mpi(q2l,im,jm,kb)
      call exchange3d_mpi(q2lb,im,jm,kb)
      call exchange3d_mpi(t,im,jm,kb)
      call exchange3d_mpi(tb,im,jm,kb)
      call exchange3d_mpi(s,im,jm,kb)
      call exchange3d_mpi(sb,im,jm,kb)
      call exchange3d_mpi(rho,im,jm,kb)
      call exchange3d_mpi(kh,im,jm,kb)
      call exchange3d_mpi(km,im,jm,kb)
      call exchange3d_mpi(u,im,jm,kb)
      call exchange3d_mpi(ub,im,jm,kb)
      call exchange3d_mpi(uf,im,jm,kb)
      call exchange3d_mpi(v,im,jm,kb)
      call exchange3d_mpi(vb,im,jm,kb)
      call exchange3d_mpi(vf,im,jm,kb)
      call exchange3d_mpi(w,im,jm,kb)
      call exchange3d_mpi(wr,im,jm,kb)

      call exchange3d_mpi(el_amp,im,jm,ncons)

      start(1) = i_global(1)
      start(2) = j_global(1)
      edge(1) = im
      edge(2) = jm
      out2 = real(aam2d,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,aam2d_varid,start,edge,out2 )      &
                , 'nf_put_var_real: aam2d @ '//out_file )
      out2 = real(advua,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,advua_varid,start,edge,out2 )      &
                , 'nf_put_var_real: advua @ '//out_file )
      out2 = real(advva,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,advva_varid,start,edge,out2 )      &
                , 'nf_put_var_real: advva @ '//out_file )
      out2 = real(adx2d,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,adx2d_varid,start,edge,out2 )      &
                , 'nf_put_var_real: adx2d @ '//out_file )
      out2 = real(ady2d,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,ady2d_varid,start,edge,out2 )      &
                , 'nf_put_var_real: ady2d @ '//out_file )
      out2 = real(drx2d,4)
      call check( nfmpi_put_vara_real_all                    &
                  ( ncid,drx2d_varid,start,edge,out2 )      &
                , 'nf_put_var_real: drx2d @ '//out_file )
      out2 = real(dry2d,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,dry2d_varid,start,edge,out2 )      &
                , 'nf_put_var_real: dry2d @ '//out_file )
      out2 = real(egb,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,egb_varid,start,edge,out2 )      &
                , 'nf_put_var_real: egb @ '//out_file )
      out2 = real(egf,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,egf_varid,start,edge,out2 )      &
                , 'nf_put_var_real: egf @ '//out_file )
      out2 = real(el,4)
      call check( nfmpi_put_vara_real_all                     &
                  ( ncid,el_varid,start,edge,out2 )      &
                , 'nf_put_var_real: el @ '//out_file )
      out2 = real(elb,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,elb_varid,start,edge,out2 )      &
                , 'nf_put_var_real: elb @ '//out_file )
      out2 = real(elf,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,elf_varid,start,edge,out2 )    &
                , 'nf_put_var_real: elf @ '//out_file )
      out2 = real(et,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,et_varid,start,edge,out2 )     &
                , 'nf_put_var_real: et @ '//out_file )
      out2 = real(etb,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,etb_varid,start,edge,out2 )      &
                , 'nf_put_var_real: etb @ '//out_file )
      out2 = real(etf,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,etf_varid,start,edge,out2 )      &
                , 'nf_put_var_real: etf @ '//out_file )
      out2 = real(fluxua,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,fluxua_varid,start,edge,out2 )    &
                , 'nf_put_var_real: fluxua @ '//out_file )
      out2 = real(fluxva,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,fluxva_varid,start,edge,out2 )    &
                , 'nf_put_var_real: fluxva @ '//out_file )
      out2 = real(ua,4)
      call check( nfmpi_put_vara_real_all              &
                  ( ncid,ua_varid,start,edge,out2 )    &
                , 'nf_put_var_real: ua @ '//out_file )
      out2 = real(uab,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,uab_varid,start,edge,out2 )    &
                , 'nf_put_var_real: uab @ '//out_file )
      out2 = real(uaf,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,uaf_varid,start,edge,out2 )    &
                , 'nf_put_var_real: uaf @ '//out_file )
      out2 = real(va,4)
      call check( nfmpi_put_vara_real_all              &
                  ( ncid,va_varid,start,edge,out2 )    &
                , 'nf_put_var_real: va @ '//out_file )
      out2 = real(vab,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,vab_varid,start,edge,out2 )    &
                , 'nf_put_var_real: vab @ '//out_file )
      out2 = real(vaf,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,vaf_varid,start,edge,out2 )    &
                , 'nf_put_var_real: vaf @ '//out_file )

      start(1) = i_global(1)
      start(2) = j_global(1)
      start(3) = 1
      start(4) = 1
      edge(1) = im
      edge(2) = jm
      edge(3) = ncons
      edge(4) = 1
      out3 = real(el_amp,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,el_amp_varid,start,edge,out3 )    &
                , 'nf_put_var_real: el_amp @ '//out_file )

      start(1) = i_global(1)
      start(2) = j_global(1)
      start(3) = 1
      start(4) = 1
      edge(1) = im
      edge(2) = jm
      edge(3) = kb
      edge(4) = 1
      out3 = real(aam,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,aam_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: aam @ '//out_file )
      out3 = real(advx,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,advx_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: advx @ '//out_file )
      out3 = real(advy,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,advy_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: advy @ '//out_file )
      out3 = real(drhox,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,drhox_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: drhox @ '//out_file )
      out3 = real(drhoy,4)
      call check( nfmpi_put_vara_real_all                  &
                  ( ncid,drhoy_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: drhoy @ '//out_file )
      out3 = real(kq,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,kq_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: kq @ '//out_file )
      out3 = real(l,4)
      call check( nfmpi_put_vara_real_all              &
                  ( ncid,l_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: l @ '//out_file )
      out3 = real(q2b,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,q2b_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: q2b @ '//out_file )
      out3 = real(q2,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,q2_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: q2 @ '//out_file )
      out3 = real(q2lb,4)
      call check( nfmpi_put_vara_real_all                 &
                  ( ncid,q2lb_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: q2lb @ '//out_file )
      out3 = real(q2l,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,q2l_varid,start,edge,out3 )     &
                , 'nf_put_vara_real: q2l @ '//out_file )
      out3 = real(t,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,t_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: t @ '//out_file )
      out3 = real(tb,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,tb_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: tb @ '//out_file )
      out3 = real(s,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,s_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: s @ '//out_file )
      out3 = real(sb,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,sb_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: sb @ '//out_file )
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
      out3 = real(ub,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,ub_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: ub @ '//out_file )
      out3 = real(u,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,u_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: u @ '//out_file )
      out3 = real(uf,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,uf_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: uf @ '//out_file )
      out3 = real(vb,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,vb_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: vb @ '//out_file )
      out3 = real(v,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,v_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: v @ '//out_file )
      out3 = real(vf,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,vf_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: vf @ '//out_file )
      out3 = real(w,4)
      call check( nfmpi_put_vara_real_all               &
                  ( ncid,w_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: w @ '//out_file )
      out3 = real(wr,4)
      call check( nfmpi_put_vara_real_all                &
                  ( ncid,wr_varid,start,edge,out3 )      &
                , 'nf_put_vara_real: wr @ '//out_file )

! close file
      call check( nf90mpi_close(ncid)                   &
                , 'nf_close: output:  @ '//out_file )


    end ! subroutine write_debug_pnetcdf
!
!______________________________________________________________________
!
    subroutine read_initial_ts_pnetcdf( temp, salt, u, v, ssh, n )
!----------------------------------------------------------------------
!  Reads initial temperature, salinity and elevation.
! TODO: Add u and v. Add check for single record to read if not climatology.
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain
!      use glob_grid  , only: fsm
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_2d &
                         , get_var_real_3d

      integer , intent(in)  :: n
      real(rk), intent(out) :: temp(im,jm,kb),salt(im,jm,kb)        &
                             , u(im,jm,kb),v(im,jm,kb), ssh(im,jm)

      integer                  ncid, status
      integer                  el_varid,sb_varid,tb_varid,u_varid,v_varid
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
      call check( nf90mpi_inq_varid( ncid, 't_init', tb_varid )  &
                , "nfmpi_inq_varid: t_init @ "//ic_path )
      call check( nf90mpi_inq_varid( ncid, 's_init', sb_varid )  &
                , "nfmpi_inq_varid: s_init @ "//ic_path )

      status = nf90mpi_inq_varid( ncid, 'einit', el_varid )
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

      status = nf90mpi_inq_varid( ncid, 'u_init', u_varid )
      if ( status == NF_NOERR ) then
        call check( get_var_real_3d(ncid,u_varid,start,edge,u)  &
                , "nf_get_var : u @ "//ic_path )
        u = u/100.
      else
        u = 0.
      end if
      status = nf90mpi_inq_varid( ncid, 'v_init', v_varid )
      if ( status == NF_NOERR ) then
        call check( get_var_real_3d(ncid,v_varid,start,edge,v)  &
                , "nf_get_var : v @ "//ic_path )
        v = v/100.
      else
        v = 0.
      end if

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


    end ! subroutine read_initial_ts_pnetcdf
!
!______________________________________________________________________
!
    subroutine read_ts_z_pnetcdf( temp, salt, ks, n, path )
!----------------------------------------------------------------------
!  Reads initial temperature, salinity and elevation.
! TODO: Add u and v. Add check for single record to read if not climatology.
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain
      use glob_grid  , only: zz, h
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_2d &
                         , get_var_real_3d

      integer                      , intent(in   ) :: ks, n
      character(len=*)             , intent(in   ) :: path
      real(rk), dimension(im,jm,kb), intent(  out) :: temp,salt

      integer                  ncid
      integer                  z_varid, sb_varid, tb_varid
      integer(MPI_OFFSET_KIND) start(4),edge(4)
      real(rk)                 tz(im,jm,ks),sz(im,jm,ks),z_z(ks)


      start = 1
      edge  = 1

!  Open netcdf file.
      if ( is_master )                                              &
           print '(''reading file `'',a,''`''/)', trim(ic_path)
      call check( nf90mpi_open(                        &
                    POM_COMM, path, NF_NOWRITE  &
                  , MPI_INFO_NULL, ncid )              &
                , "nfmpi_open: "//path )

!  Get variables. [ TODO: parameterize varnames ]
      call check( nf90mpi_inq_varid( ncid, 't', tb_varid )  &
                , "nfmpi_inq_varid: t_init @ "//ic_path )
      call check( nf90mpi_inq_varid( ncid, 's', sb_varid )  &
                , "nfmpi_inq_varid: s_init @ "//ic_path )
      call check( nf90mpi_inq_varid( ncid, 'z', z_varid )  &
                , "nfmpi_inq_varid: z_init @ "//ic_path )

      start(1) = 1
      edge(1) = ks
      call check( nf90mpi_get_var_all(ncid,z_varid,z_z,start,edge)  &
                , "nf_get_var : z @ "//ic_path )

      start(1) = i_global(1)
      start(2) = j_global(1)
      start(3) = 1
      start(4) = n
      edge(1) = im
      edge(2) = jm
      edge(3) = ks
      edge(4) = 1
      call check( get_var_real_3d(ncid,tb_varid,start,edge,tz)  &
                , "nf_get_var : temp @ "//ic_path )
      call check( get_var_real_3d(ncid,sb_varid,start,edge,sz)  &
                , "nf_get_var : salt @ "//ic_path )

!  Close file.
      call check( nf90mpi_close(ncid), "nf_close @ "//ic_path)

      call ztosig(z_z,tz,zz,h,temp,im,jm,ks,kb,n_west,n_east,n_south,n_north)
      call ztosig(z_z,sz,zz,h,salt,im,jm,ks,kb,n_west,n_east,n_south,n_north)


    end ! subroutine read_ts_z_pnetcdf
!
!___________________________________________________________________
!
    integer function file_open_nc( path )
!-------------------------------------------------------------------
!  Opens netcdf file for reading.
!___________________________________________________________________
!
      use glob_domain, only: POM_COMM
      use mpi        , only: MPI_INFO_NULL
      use pnetcdf    , only: nf90mpi_open, NF_NOERR, NF_NOWRITE

      implicit none

      character(len=*), intent(in) :: path

      integer                         status


      status = nf90mpi_open( POM_COMM, trim(path), NF_NOWRITE   &
                           , MPI_INFO_NULL, file_open_nc )
      if ( status /= NF_NOERR ) then
        call msg_print("", 2, "Failed to open `"//trim(path)//"`")
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
          print '(/a,a)', 'IO error: ', routine
          print *, nf90mpi_strerror(status)
          stop
        end if
      end if


    end ! subroutine check
!

end module io
