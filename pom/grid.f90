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

  use glob_const , only: PATH_LEN, rk, VAR_LEN
  use glob_domain, only: im, jm, kb

  implicit none

  public

  integer(1), parameter ::  &
    vSIGMA        = 0       &  ! Native sigma coordinates
  , vGEOPOTENTIAL = 1          ! Z-coordinates
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

! Vertical discretisation switch
  integer(1) v_type

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
! Model discretisation
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:,:)   :: &
    art              & ! cell area centered on T grid points
  , aru              & ! cell area centered on U grid points
  , arv              & ! cell area centered on V grid points
  , cor              & ! coriolis parameter
  , dx               & ! grid spacing in x
  , dy               & ! grid spacing in y
  , east_c           & ! horizontal coordinate of cell corner points in x
  , east_e           & ! horizontal coordinate of elevation points in x
  , east_u           & ! horizontal coordinate of U points in x
  , east_v           & ! horizontal coordinate of V points in x
  , h                & ! bottom depth
  , norm_h           & ! sigma-coordinate normalization values (for averages correction)
  , north_c          & ! horizontal coordinate of cell corner points in y
  , north_e          & ! horizontal coordinate of elevation points in y
  , north_u          & ! horizontal coordinate of U points in y
  , north_v          & ! horizontal coordinate of V points in y
  , rot                ! rotation angle

  real(rk)                     &
       , allocatable           &
       , dimension(:,:,:)   :: &
    dum              & ! mask for u-velocity
  , dvm              & ! mask for v-velocity
  , dwm              & ! mask for w-velocity
  , dz               & ! sigma-"layer" thickness (z(k)-z(k+1))
  , dzz              & ! thickness of two halves of neigbour sigma-"layers" (zz(k)-zz(k+1))
  , fsm              & ! mask for scalar variables
  , z                & ! sigma coordinate from z=0 (surface) to z=-1 (bottom)
  , zz                 ! sigma coordinate, intermediate between z

  contains
!______________________________________________________________________
!
    subroutine initialize_mod( config_nml )
!----------------------------------------------------------------------
!  Initializes module.
!______________________________________________________________________
!
      implicit none

      character(*), intent(in) :: config_nml

      namelist/grid/                                              &
          grid_path, v_type                                       &
        ,  dx_name , dy_name, ec_name, ee_name, eu_name, ev_name  &
        , fsm_name ,  h_name, nc_name, ne_name, nu_name, nv_name  &
        , rot_name ,  z_name, zz_name

      integer pos

! Create module arrays
      call allocate_arrays

! Initialize with default values
      grid_path = "in/grid/"

      v_type = 0

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

! Manage input_files values
      pos = len(trim(grid_path))
      if ( grid_path(pos:pos) == "/" ) then
        grid_path = trim(grid_path)//"grid.nc"
      end if

! Print notification
      call msg_print("GRID MODULE READY", 1, "")


    end ! subroutine initialize_io
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!----------------------------------------------------------------------
! called by: initialize_mod    [air]
!______________________________________________________________________
!
      implicit none


! Allocate core arrays
      allocate(         &
        art    (im,jm)  &
      , aru    (im,jm)  &
      , arv    (im,jm)  &
      , cor    (im,jm)  &
      , dx     (im,jm)  &
      , dy     (im,jm)  &
      , east_c (im,jm)  &
      , east_e (im,jm)  &
      , east_u (im,jm)  &
      , east_v (im,jm)  &
      , h      (im,jm)  &
      , norm_h (im,jm)  &
      , north_c(im,jm)  &
      , north_e(im,jm)  &
      , north_u(im,jm)  &
      , north_v(im,jm)  &
      , rot    (im,jm)  &
      )

      allocate(         &
        dz (im,jm,kb)   &
      , dzz(im,jm,kb)   &
      , dum(im,jm,kb)   &
      , dvm(im,jm,kb)   &
      , dwm(im,jm,kb)   &
      , fsm(im,jm,kb)   &
      , z  (im,jm,kb)   &
      , zz (im,jm,kb)   &
      )


    end subroutine allocate_arrays
!
!______________________________________________________________________
!
    subroutine read_grid
!----------------------------------------------------------------------
!  Read grid data.
!______________________________________________________________________
!
!      use mpi        , only: MPI_OFFSET_KIND
      use config     , only: n1d
      use glob_const , only: DEG2RAD, Ohm, rk, SMALL
      use io         , only: check   , file_close, file_open  &
                           , is_error, var_rank  , var_read
      use glob_domain, only: i_global, imm1  , is_master  &
                           , j_global, jmm1  , kb, kbm1   &
                           , n_south , n_west
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      integer                                   file_id, i, j, k, z_rnk
      !integer(MPI_OFFSET_KIND), dimension(4) :: start, stride
      integer(8), dimension(4) :: start, stride


      start  = [ i_global(1), j_global(1),  1, 1 ]
      stride = [ im         , jm         , kb, 1 ]

      if ( is_master ) then
        call msg_print("", 6, "Reading grid `"//trim( grid_path )//"`")
      end if
! open netcdf file
      file_id = file_open( grid_path, NF90_NOWRITE )
      if ( is_error( file_id ) ) then
        call msg_print( "GRID CANNOT BE READ", 3, "Make sure the file exists and valid." )
      end if

      z_rnk = var_rank( file_id, z_name )

      if ( z_rnk == 3 ) then
        call check( var_read( file_id,  z_name,  z, start, stride )  &
                  , '[grid]:read_grid:var_read `'//trim(z_name) )
        call check( var_read( file_id, zz_name, zz, start, stride )  &
                  , '[grid]:read_grid:var_read `'//trim(zz_name) )
        call check( var_read( file_id, fsm_name, fsm                 &
                            , start(1:3), stride(1:3) )              &
                  , '[grid]:read_grid:var_read `'//trim(fsm_name) )
      else
        call check( var_read( file_id,  z_name,  z(1,1,:)           &
                            , start(3:3), stride(3:3) )             &
                  , '[grid]:read_grid:var_read `'//trim(z_name) )
        call check( var_read( file_id, zz_name, zz(1,1,:)           &
                            , start(3:3), stride(3:3) )             &
                  , '[grid]:read_grid:var_read `'//trim(zz_name) )
        call check( var_read( file_id, fsm_name, fsm(:,:,1)         &
                            , start(1:3), stride(1:3) )             &
                  , '[grid]:read_grid:var_read `'//trim(fsm_name) )
        do k = 1, kb
          z (:,:,k) = z (1,1,k)
          zz(:,:,k) = zz(1,1,k)
          fsm(:,:,k) = fsm(:,:,1)
        end do
      end if

      norm_h = 1.
      do k = kbm1, 1, -1
        where( fsm(:,:,k) == 0. ) norm_h = -z(:,:,k)
      end do
      where( norm_h == 0. ) norm_h = 1.
      norm_h = 1./norm_h

      call check( var_read( file_id,  dx_name, dx                 &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(dx_name) )
      call check( var_read( file_id,  dy_name, dy                 &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(dy_name) )
      call check( var_read( file_id,  eu_name, east_u             &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(eu_name) )
      call check( var_read( file_id,  ev_name, east_v             &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(ev_name) )
      call check( var_read( file_id,  ee_name, east_e             &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(ee_name) )
      call check( var_read( file_id,  nu_name, north_u            &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(nu_name) )
      call check( var_read( file_id,  nv_name, north_v            &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(nv_name) )
      call check( var_read( file_id,  ne_name, north_e            &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(ne_name) )
      call check( var_read( file_id, rot_name, rot                &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(rot_name) )
      call check( var_read( file_id,   h_name, h                  &
                          , start(1:3), stride(1:3) )             &
                , '[grid]:read_grid:var_read `'//trim(h_name) )

! close file:
      file_id = file_close( file_id )


! derive u- and v-mask from rho-mask
      dum(2:im,:,:) = fsm(2:im,:,:)*fsm(1:imm1,:,:)
      dvm(:,2:jm,:) = fsm(:,2:jm,:)*fsm(:,1:jmm1,:)
      call exchange3d_mpi(dum,im,jm,kb)
      call exchange3d_mpi(dvm,im,jm,kb)
      dwm(:,:,2:kb) = fsm(:,:,2:kb)*fsm(:,:,1:kbm1)
      
      if ( n_west  == -1 ) dum(1,:,:) = dum(2,:,:)
      if ( n_south == -1 ) dvm(:,1,:) = dvm(:,2,:)

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
      if ( n1d /= 0 ) then
        dx = 1.e9
        dy = 1.e9
      end if
!lyo:scs1d:end:
!
! derive vertical grid variables
!      z(kb) = -1. ! :DIRTY HACK!!! (for stupid grid)
!      zz(kb-1) = .5*(z(kb-1)+z(kb)) ! :
!      zz(kb) = 2.*zz(kb-1)-zz(kb-2) ! :
!!      if (my_task==1) fsm(:,jm-8:jm) = 0.
!      h = 1500.
      dz (:,:,1:kb-1) = z (:,:,1:kb-1) - z (:,:,2:kb)
      dzz(:,:,1:kb-1) = zz(:,:,1:kb-1) - zz(:,:,2:kb)
      dz (:,:,kb) = dz (:,:,kb-1)
      dzz(:,:,kb) = dzz(:,:,kb-1)

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

      if ( n_west == -1 ) then
        aru(1,:) = aru(2,:)
        arv(1,:) = arv(2,:)
      end if

      if ( n_south == -1 ) then
        aru(:,1) = aru(:,2)
        arv(:,1) = arv(:,2)
      end if

      call print_config


    end ! subroutine read_grid
!
!______________________________________________________________________
!
    subroutine print_config
!----------------------------------------------------------------------
!  Prints module configuration.
!______________________________________________________________________
!
      use glob_domain, only: is_master

      implicit none

      integer k


      if ( .not.is_master ) return

! Print initial (core) summary
      print '(/a)', "---- VERTICAL DISCRETISATION --------------------"
      print '(/3x,a,9x,a,8x,a,8x,a,8x,a)', 'k','z','zz','dz','dzz'
      do k = 1, kb
        print '(1x,i5,4f10.3)', k, z(2,2,k), zz(2,2,k)  &
                                 ,dz(2,2,k),dzz(2,2,k)
      end do
      
      print '(/a)', "---- HORIZONTAL STATS (only master atm) ---------"
      print '(/1x,a,"[",f8.3,":",f8.3,"]")', "Latitude   : "  &
          , minval(north_e), maxval(north_e)
      print '( 1x,a,"[",f8.3,":",f8.3,"]")', "Longitude  : "  &
          , minval(east_e) , maxval(east_e)
      print '(/1x,a,"[",f10.2,":",f10.2,"]")', "X cell size: "  &
          , minval(dx), maxval(dx)
      print '( 1x,a,"[",f10.2,":",f10.2,"]")', "Y cell size: "  &
          , minval(dy), maxval(dy)
      print '( 1x,a,"[",f10.0,":",f10.0,"]")', "Cell area  : "  &
          , minval(art), maxval(art)


    end ! subroutine print_config


end module grid
