!______________________________________________________________________
!
! Module `GRID` (grid.f90)
!----------------------------------------------------------------------
!  Module for managing grid configuration.
!
!  Author  : RinceWND
!  Created : 2019-04-11
!______________________________________________________________________
!
module grid

  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  public

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(PATH_LEN)  &
    grid_path            ! Full path to grid file
!______________________________________________________________________
!  If `bc_path`, `bulk_path` or `flux_path` end with a `.` (dot)
! then these files will be read as <path>YYYY.<FORMAT_EXT>.
!  If any path ends with `/` then the default filename will be
! appended to the path.

!----------------------------------------------------------------------
! Varnames configuration
!----------------------------------------------------------------------
  character(VAR_LEN)  & ! Input variables names:
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

!----------------------------------------------------------------------
! Vertical discretization
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:)     :: &
    dz               & ! z(k)-z(k+1)
  , dzz              & ! zz(k)-zz(k+1)
  , z                & ! sigma coordinate from z=0 (surface) to z=-1 (bottom)
  , zz                 ! sigma coordinate, intermediate between z

!----------------------------------------------------------------------
! Horizontal discretization
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:,:)   :: &
    art              & ! cell area centered on T grid points
  , aru              & ! cell area centered on U grid points
  , arv              & ! cell area centered on V grid points
  , cor              & ! coriolis parameter
  , dum              & ! mask for u velocity
  , dvm              & ! mask for v velocity
  , dx               & ! grid spacing in x
  , dy               & ! grid spacing in y
  , east_c           & ! horizontal coordinate of cell corner points in x
  , east_e           & ! horizontal coordinate of elevation points in x
  , east_u           & ! horizontal coordinate of U points in x
  , east_v           & ! horizontal coordinate of V points in x
  , fsm              & ! mask for scalar variables
  , h                & ! bottom depth
  , north_c          & ! horizontal coordinate of cell corner points in y
  , north_e          & ! horizontal coordinate of elevation points in y
  , north_u          & ! horizontal coordinate of U points in y
  , north_v          & ! horizontal coordinate of V points in y
  , rot                ! rotation angle

  contains
!______________________________________________________________________
!
    subroutine initialize_mod( config_nml )
!----------------------------------------------------------------------
!  Initializes module.
!______________________________________________________________________
!
      use glob_domain, only: is_master

      implicit none

      character(*), intent(in) :: config_nml

      namelist/grid/                                             &
          grid_path                                              &
        ,  dx_name, dy_name, ec_name, ee_name, eu_name, ev_name  &
        , fsm_name,  h_name, nc_name, ne_name, nu_name, nv_name  &
        , rot_name,  z_name, zz_name

      integer pos


! Initialize with default values
      grid_path = "in/grid/"

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
      read ( 73, nml = grid )
      close( 73 )

! Manage inmput_files values
      pos = len(trim(grid_path))
      if ( grid_path(pos:pos) == "/" ) then
        grid_path = trim(grid_path)//"grid.nc"
      end if

! Print config
      if ( is_master ) then
        print *, "--------------|----------------------------"
        print *, " Grid         : ", trim(grid_path)
        print *, "--------------|----------------------------"
        call msg_print("GRID MODULE INITIALIZED", 1, "")
      end if


    end ! subroutine initialize_io
!
!______________________________________________________________________
!
    subroutine read_grid
!----------------------------------------------------------------------
!  Read grid data.
!______________________________________________________________________
!
      use mpi        , only: MPI_OFFSET_KIND
      use config     , only: n1d
      use glob_const , only: DEG2RAD, Ohm, rk, SMALL
      use glob_domain
      use io         , only: check, file_open, file_close, var_read
      use glob_domain, only: im, jm, kb, n_south, n_west

      implicit none

      integer                                   i, j, k, ncid
      integer(MPI_OFFSET_KIND), dimension(4) :: start, stride


      start  = 1
      stride = 1

      if ( is_master ) print '(/''Reading grid '',a)', trim( grid_path )
! open netcdf file
      ncid = file_open( grid_path )

      start(1)  =  1
      stride(1) = kb
      call check( var_read( ncid,  z_name,  z, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(z_name) )
      call check( var_read( ncid, zz_name, zz, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(zz_name) )

      start(1)  = i_global(1)
      start(2)  = j_global(1)
      stride(1) = im
      stride(2) = jm
      call check( var_read( ncid,  dx_name,      dx, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(dx_name) )
      call check( var_read( ncid,  dy_name,      dy, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(dy_name) )
      call check( var_read( ncid,  eu_name,  east_u, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(eu_name) )
      call check( var_read( ncid,  ev_name,  east_v, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(ev_name) )
      call check( var_read( ncid,  ee_name,  east_e, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(ee_name) )
      call check( var_read( ncid,  nu_name, north_u, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(nu_name) )
      call check( var_read( ncid,  nv_name, north_v, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(nv_name) )
      call check( var_read( ncid,  ne_name, north_e, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(ne_name) )
      call check( var_read( ncid, rot_name,     rot, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(rot_name) )
      call check( var_read( ncid,   h_name,       h, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(h_name) )
      call check( var_read( ncid, fsm_name,     fsm, start, stride )  &
                , '[grid]:read_grid:var_read `'//trim(fsm_name) )

! close file:
      ncid = file_close( ncid )


! derive "u" and "v" masks from "rho" mask
      dvm = fsm
      dum = fsm
      do j=1,jm-1
        do i=1,im
          if (fsm(i,j)==0..and.fsm(i,j+1)/=0.) dvm(i,j+1) = 0.
        end do
      end do
      do j=1,jm
        do i=1,im-1
          if (fsm(i,j)==0..and.fsm(i+1,j)/=0.) dum(i+1,j) = 0.
        end do
      end do
      call exchange2d_mpi(dum,im,jm)
      call exchange2d_mpi(dvm,im,jm)

! modify SMALL to zero distances to 1 meter
      where( dx < SMALL ) dx = 1.
      where( dy < SMALL ) dy = 1.

!fhx:debug:close north boundary for NWATL
!       if(n_north.eq.-1) then
!          fsm(:,jm)=0.0
!          h(:,jm)=1.0
!       endif
!
!lyo:scs1d:beg:Artificially increase dx & dy for 1d-simulations
      if (n1d.ne.0) then
        do j=1,jm; do i=1,im
          dx(i,j)=1.e9; dy(i,j)=1.e9
        enddo; enddo
      endif
!lyo:scs1d:end:
!
! derive vertical grid variables
!      z(kb) = -1. ! :DIRTY HACK!!! (for stupid grid)
!      zz(kb-1) = .5*(z(kb-1)+z(kb)) ! :
!      zz(kb) = 2.*zz(kb-1)-zz(kb-2) ! :
!!      if (my_task==1) fsm(:,jm-8:jm) = 0.
!      h = 1500.
      do k=1,kb-1
        dz(k) = z(k)- z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb) = dz(kb-1) !=0. !lyo:20110202 =0 is dangerous but checked
      dzz(kb)=dzz(kb-1) !=0. !thro' code - and found to be ok.

      call print_grid_vert

! set up Coriolis parameter
      cor = 2.*Ohm*sin(north_e*DEG2RAD)

! inertial period for temporal filter
!      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0

! calculate areas of "rho" cells
      art = dx*dy

! calculate areas of "u" and "v" cells
      do j=2,jm
        do i=2,im
          aru(i,j)=.25*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
      call exchange2d_mpi(aru,im,jm)
      call exchange2d_mpi(arv,im,jm)

      if (n_west == -1) then
        do j=1,jm
          aru(1,j)=aru(2,j)
          arv(1,j)=arv(2,j)
        end do
      end if

      if (n_south == -1) then
        do i=1,im
          aru(i,1)=aru(i,2)
          arv(i,1)=arv(i,2)
        end do
      end if


    end ! subroutine read_grid
!
!______________________________________________________________________
!
    subroutine out_init_define( ncid, x_dimid, y_dimid, z_dimid )
!----------------------------------------------------------------------
!  Defines variables for output.
!______________________________________________________________________
!
      use io     , only: att_write, var_define
      use pnetcdf, only: NF90_FLOAT

      implicit none

      integer, intent(in) :: ncid, x_dimid, y_dimid, z_dimid

      integer               :: varid
      integer, dimension(2) :: vdims

      vdims = 0

      vdims(1) = z_dimid
      varid = var_define( ncid, 'sigma', 1, vdims, NF90_FLOAT  &
                        , "sigma of cell face", "sigma_level"  &
                        , -1, 0., ' ', .false. )
      call att_write( ncid, varid, 'standard_name'             &
                                 , "ocean_sigma_coordinate" )
      call att_write( ncid, varid, 'formula_terms'                 &
                                 , "sigma: z eta: elb depth: h" )

      varid = var_define( ncid, 'zz', 1, vdims, NF90_FLOAT       &
                        , "sigma of cell centre", "sigma_level"  &
                        , -1, 0., ' ', .false. )
      call att_write( ncid, varid, 'standard_name'             &
                                 , "ocean_sigma_coordinate" )
      call att_write( ncid, varid, 'formula_terms'                  &
                                 , "sigma: zz eta: elb depth: h" )

      vdims(1) = x_dimid
      vdims(2) = y_dimid
      varid = var_define( ncid, 'dx', 2, vdims, NF90_FLOAT    &
                        , "grid increment in x", "metre"      &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'dy', 2, vdims, NF90_FLOAT    &
                        , "grid increment in y", "metre"      &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'east_u', 2, vdims, NF90_FLOAT  &
                        , "easting of u-points", "degree"       &
                        , -1, 0., "north_u east_u", .true. )
      varid = var_define( ncid, 'east_v', 2, vdims, NF90_FLOAT  &
                        , "easting of v-points", "degree"       &
                        , -1, 0., "north_v east_v", .true. )
      varid = var_define( ncid, 'east_e', 2, vdims, NF90_FLOAT     &
                        , "easting of elevation points", "degree"  &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'east_c', 2, vdims, NF90_FLOAT  &
                        , "easting of cell corners", "degree"   &
                        , -1, 0., "east_c north_c", .true. )
      varid = var_define( ncid, 'north_u', 2, vdims, NF90_FLOAT  &
                        , "northing of u-points", "degree"       &
                        , -1, 0., "north_u east_u", .true. )
      varid = var_define( ncid, 'north_v', 2, vdims, NF90_FLOAT  &
                        , "northing of v-points", "degree"       &
                        , -1, 0., "north_v east_v", .true. )
      varid = var_define( ncid, 'north_e', 2, vdims, NF90_FLOAT     &
                        , "northing of elevation points", "degree"  &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'north_c', 2, vdims, NF90_FLOAT  &
                        , "northing of cell corners", "degree"   &
                        , -1, 0., "east_c north_c", .true. )
      varid = var_define( ncid, 'rot', 2, vdims, NF90_FLOAT               &
                        , "rotation angle of x-axis wrt. east", "radian"  &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'h', 2, vdims, NF90_FLOAT     &
                        , "undisturbed water depth", "metre"  &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'fsm', 2, vdims, NF90_FLOAT     &
                        , "free surface mask", "dimensionless"  &
                        , -1, 0., "north_e east_e", .true. )
      varid = var_define( ncid, 'dum', 2, vdims, NF90_FLOAT   &
                        , "u-velocity mask", "dimensionless"  &
                        , -1, 0., "north_u east_u", .true. )
      varid = var_define( ncid, 'dvm', 2, vdims, NF90_FLOAT   &
                        , "v-velocity mask", "dimensionless"  &
                        , -1, 0., "north_v east_v", .true. )


    end ! subroutine out_init_define
!
!______________________________________________________________________
!
    subroutine out_write( ncid )
!----------------------------------------------------------------------
!  Writes to file.
!______________________________________________________________________
!
      use glob_domain, only: im, i_global, jm, j_global, kb
      use io         , only: var_write
      use mpi        , only: MPI_OFFSET_KIND

      implicit none

      integer, intent(in) :: ncid

      integer(MPI_OFFSET_KIND), dimension(4) :: start, stride


      start  = [  1, 1, 1, 1 ]
      stride = [ kb, 1, 1, 1 ]
      call var_write( ncid, 'sigma',  z, start, stride )
      call var_write( ncid,    'zz', zz, start, stride )

      start  = [ i_global(1), j_global(1), 1, 1 ]
      stride = [ im         , jm         , 1, 1 ]
      call var_write( ncid, 'dum'    , dum    , start, stride )
      call var_write( ncid, 'dvm'    , dvm    , start, stride )
      call var_write( ncid, 'dx'     , dx     , start, stride )
      call var_write( ncid, 'dy'     , dy     , start, stride )
      call var_write( ncid, 'east_c' , east_c , start, stride )
      call var_write( ncid, 'east_e' , east_e , start, stride )
      call var_write( ncid, 'east_u' , east_u , start, stride )
      call var_write( ncid, 'east_v' , east_v , start, stride )
      call var_write( ncid, 'fsm'    , fsm    , start, stride )
      call var_write( ncid, 'h'      , h      , start, stride )
      call var_write( ncid, 'north_c', north_c, start, stride )
      call var_write( ncid, 'north_e', north_e, start, stride )
      call var_write( ncid, 'north_u', north_u, start, stride )
      call var_write( ncid, 'north_v', north_v, start, stride )
      call var_write( ncid, 'rot'    , rot    , start, stride )


    end ! subroutine out_write


end module grid