! initialize.f

!  Initializes POM: defines constants, reads initial values, grid and
! initial conditions.

!______________________________________________________________________
!
      subroutine initialize
!----------------------------------------------------------------------
!  Initializes POM.
!______________________________________________________________________
!
      use config
      use glob_const , only: initialize_constants,rhoref
      use glob_domain
      use grid       , only: east_e, north_e, h
      use grid       , only: initialize_grid_module => initialize_mod
      use glob_ocean , only: aamfac, d, dt, el, et
     &                     , initialize_ocean => allocate_arrays
      use glob_out   , only: out_record, iprint
     &                     , initialize_mean_arrays => allocate_arrays
      use model_run  , only: dtime, time, time0

      implicit none

      integer, parameter :: TEST_CASE = 0 !1


! initialize the MPI execution environment and create communicator for
!internal POM communications
      call initialize_mpi

! Initialize various constants
      call initialize_constants

! Read domain distribution (TODO: get rid of hardcoded parameters)
      call read_domain_dist( 'domain.nml' )

! Read configuration files
      call read_config( 'config.nml' )
!______________________________________________________________________
!  `read_config` initializes all parameters to their default values
! prior to reading them.

! distribute the model domain across processors
      call distribute_mpi

! distribute the coarse wind/assim data domain across processors
      call distribute_mpi_coarse   !fhx:20101206

! decide on interpolation feasibility (sets calc_interp flag)
      call check_interpolation

! Initialize general IO and config
      call initialize_grid_module( 'config.nml' )

! Initialize grid
      call initialize_grid( TEST_CASE )

! Allocate and initialize arrays
      call initialize_ocean

! The model is initialized. Below is some updates. TODO: Separate?
      call msg_print("MODEL CORE INITIALIZED", 1, "")
      ! TODO: Revised up to here. Code below is "untested".

! Initialize modules
      call initialize_modules


      if ( .not.do_restart ) then
! read initial and lateral boundary conditions
        call initial_conditions

! Initialize time
        call initialize_time

! read previous (initial) record [TODO: NOTE! This step should be done before update_initial for clim to get rmean. But what happens to boundary?]
        call modules_initial_step( dtime )

! update initial conditions and set the remaining conditions
        call update_initial
        call create_initial("out/init."//trim(netcdf_file)//".nc")

! read laterel boundary conditions
       !call lateral_boundary_conditions !ayumi:20100407

! read M2 tidal amplitude & phase
        if ( use_tide ) call read_tide  ! TODO: move to initialize_modules?

      else

! read restart data from a previous run TODO: Shoul be in init module as a variation of read_initial
        call read_restart
! update time and water column
        !time = time0

        d  = h+el
        dt = h+et

        if ( append_output ) out_record = int(time/iprint)

! Initialize time
        call initialize_time
! initialise modules
        call modules_initial_step( dtime )

      end if


! calculate the bottom friction coefficient
      call bottom_friction

!lyo:pac10:Keep this to increase aam near (lono,lato)
!     xs=1.5; ys=1.5; fak=0.5; !lyo:pac10:Specified in pom.h
!     lono=-74.59; lato=36.70; !Correspond to lyo:coarse io=234; jo=185
      if ( lono/=999. .and. lato/=999. ) then
!set lono.OR.lato=999. to skip
        call incmix(aamfac,im,jm,1,lono,lato,east_e,north_e,xs,ys,fak)
      end if

! Output routines
      call initialize_mean_arrays

! write grid (NO initial conditions)
      if ( output_flag == 1 ) then
        call create_output("out/"//trim(netcdf_file)//".nc")
      end if

! check for errors
      call   sum0i_mpi(error_status,master_task)
      call bcast0i_mpi(error_status,master_task)
      if ( error_status /= 0 ) then
        if ( is_master ) print '(/a/a)'
     &                       , "Something went wrong during init"
     &                       , "POM terminated with error"
        call finalize_mpi
        stop
      end if

      call msg_print("INITIALIZATION COMPLETE", 1, "")


      end ! subroutine initialize
!
!______________________________________________________________________
!
      subroutine initialize_grid( TEST_CASE )
!----------------------------------------------------------------------
!  Reads of generates grid.
!______________________________________________________________________
!
      use config     , only: n1d
      use glob_const , only: DEG2RAD, Ohm, rk, SMALL
      use glob_domain, only: im, imm1, jm, jmm1, kb, n_west, n_south
      use grid
      use glob_ocean , only: d, dt, el, et

      implicit none

      integer, intent(in) :: TEST_CASE

      integer i,j    ,ierr


      if ( TEST_CASE == 0 ) then
! read in grid data
        call read_grid
      else
        call make_grid(TEST_CASE)
      end if

      call check_cfl_min


      end ! subroutine initialize_grid
!
!______________________________________________________________________
!
      subroutine initial_conditions
!----------------------------------------------------------------------
!  Set up initial and lateral boundary conditions.
!______________________________________________________________________
!
      use bry        , only: initial_conditions_boundary
      use config     , only: initial_file
      use glob_const , only: rk
      use glob_domain
      use grid       , only: dz, fsm
      use io
      use mpi        , only: MPI_OFFSET_KIND
      use glob_ocean , only: elb, rho, s, sb, t, tb, uab, ub, vab, vb
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      logical                      fexist
      integer                      file_id, i, j, k, record, status
      integer(MPI_OFFSET_KIND)
     &       , dimension(4)     :: edge, start
      character(32) :: el_name, s_name, t_name, u_name, v_name


      el_name = "elev"
       s_name = "mean_s" !"salt"
       t_name = "mean_t" !"temp"
       u_name = "u"
       v_name = "v"
      record = 1

! Initialize main parameters with "common" values
      tb  = 15.
      sb  = 33.
      elb =  0.
      ub  =  0.
      vb  =  0.

!      record = dtime % month
! Set computational-thread-specific domain to read in parallel
      start = [ i_global(1), j_global(1),  1, record ]
      edge  = [ im         , jm         , kb,      1 ]

! Read initial file
      call msg_print("", 6, "Read initial conditions:")
      file_id = file_open( trim(initial_file), NF90_NOWRITE )
      if ( .not.is_error(file_id) ) then
        status = var_read( file_id,  t_name,  tb, start, edge )
        status = var_read( file_id,  s_name,  sb, start, edge )
        status = var_read( file_id,  u_name,  ub, start, edge )
        status = var_read( file_id,  v_name,  vb, start, edge )
        status = var_read( file_id, el_name, elb, start, edge )
!        call read_initial_ts_pnetcdf( tb, sb, ub, vb, elb, 1 ) !dtime%month )
!        call read_ts_z_pnetcdf( tb, sb, 40, dtime%month, ic_path )
      end if

! Derive barotropic velocities
      uab = 0.
      vab = 0.
      do k = 1, kbm1
        sb(:,:,k) = sb(:,:,k) * fsm
        tb(:,:,k) = tb(:,:,k) * fsm
        uab(:,:)  = uab(:,:) + ub(:,:,k)*dz(k)
        vab(:,:)  = vab(:,:) + vb(:,:,k)*dz(k)
      end do

! Calculate initial density.
      call dens(sb,tb,rho)


      end
!
!______________________________________________________________________
!
      subroutine read_restart
!----------------------------------------------------------------------
!  Set up initial and lateral boundary conditions.
!______________________________________________________________________
!
      use air        , only: wusurf, wvsurf
      use config     , only: initial_file
      use glob_const , only: rk
      use glob_domain
      use io
      use mpi        , only: MPI_OFFSET_KIND
      use model_run  , only: time
      use glob_ocean
      use pnetcdf    , only: NF90_NOWRITE

      implicit none

      logical                      fexist
      integer                      file_id, i, j, k, record, status
      integer(MPI_OFFSET_KIND)
     &       , dimension(3)     :: edge, start


      start = [ i_global(1), j_global(1),  1 ]
      edge  = [ im         , jm         , kb ]

! Read initial file
      call msg_print("", 6, "Read restart file: `"
     &                      //trim(initial_file)//"`")
      file_id = file_open( trim(initial_file), NF90_NOWRITE )
      if ( .not.is_error(file_id) ) then
        status = var_read(file_id,'time'  ,time  , start(3:3) )
        status = var_read(file_id,'wusurf',wusurf, start(1:2),edge(1:2))
        status = var_read(file_id,'wvsurf',wvsurf, start(1:2),edge(1:2))
        status = var_read(file_id,'wubot' ,wubot , start(1:2),edge(1:2))
        status = var_read(file_id,'wvbot' ,wvbot , start(1:2),edge(1:2))
        status = var_read(file_id,'aam2d' ,aam2d , start(1:2),edge(1:2))
        status = var_read(file_id,'ua'    ,ua    , start(1:2),edge(1:2))
        status = var_read(file_id,'uab'   ,uab   , start(1:2),edge(1:2))
        status = var_read(file_id,'va'    ,va    , start(1:2),edge(1:2))
        status = var_read(file_id,'vab'   ,vab   , start(1:2),edge(1:2))
        status = var_read(file_id,'el'    ,el    , start(1:2),edge(1:2))
        status = var_read(file_id,'elb'   ,elb   , start(1:2),edge(1:2))
        status = var_read(file_id,'et'    ,et    , start(1:2),edge(1:2))
        status = var_read(file_id,'etb'   ,etb   , start(1:2),edge(1:2))
        status = var_read(file_id,'egb'   ,egb   , start(1:2),edge(1:2))
        status = var_read(file_id,'utb'   ,utb   , start(1:2),edge(1:2))
        status = var_read(file_id,'vtb'   ,vtb   , start(1:2),edge(1:2))
        status = var_read(file_id,'adx2d' ,adx2d , start(1:2),edge(1:2))
        status = var_read(file_id,'ady2d' ,ady2d , start(1:2),edge(1:2))
        status = var_read(file_id,'advua' ,advua , start(1:2),edge(1:2))
        status = var_read(file_id,'advva' ,advva , start(1:2),edge(1:2))
        status = var_read(file_id,'u'     ,u     , start     ,edge     )
        status = var_read(file_id,'ub'    ,ub    , start     ,edge     )
        status = var_read(file_id,'v'     ,v     , start     ,edge     )
        status = var_read(file_id,'vb'    ,vb    , start     ,edge     )
        status = var_read(file_id,'w'     ,w     , start     ,edge     )
        status = var_read(file_id,'t'     ,t     , start     ,edge     )
        status = var_read(file_id,'tb'    ,tb    , start     ,edge     )
        status = var_read(file_id,'s'     ,s     , start     ,edge     )
        status = var_read(file_id,'sb'    ,sb    , start     ,edge     )
        status = var_read(file_id,'rho'   ,rho   , start     ,edge     )
        status = var_read(file_id,'km'    ,km    , start     ,edge     )
        status = var_read(file_id,'kh'    ,kh    , start     ,edge     )
        status = var_read(file_id,'kq'    ,kq    , start     ,edge     )
        status = var_read(file_id,'l'     ,l     , start     ,edge     )
        status = var_read(file_id,'q2'    ,q2    , start     ,edge     )
        status = var_read(file_id,'q2b'   ,q2b   , start     ,edge     )
        status = var_read(file_id,'aam'   ,aam   , start     ,edge     )
        status = var_read(file_id,'q2l'   ,q2l   , start     ,edge     )
        status = var_read(file_id,'q2lb'  ,q2lb  , start     ,edge     )
      else
        call msg_print("", 2, "FAILED...")
        call finalize_mpi
        stop
      end if

! Close file
      file_id = file_close( file_id )


      end
!______________________________________________________________________
!
      subroutine initialize_time
!----------------------------------------------------------------------
!  Initialize `dtime` variable
!______________________________________________________________________
!
        use config     , only: do_restart
        use model_run  , only: dtime, time, time0, time_start
        use module_time

        implicit none


! Initialize time
        dtime = str2date( time_start ) + int(time*86400.)
        time0 = time


      end subroutine
!______________________________________________________________________
!
      subroutine create_output( out_file )
!----------------------------------------------------------------------
!  Create output file.
!______________________________________________________________________
!
        use air
        use config     , only: use_ice, mode, title
        use glob_const , only: MODE_BAROTROPIC, rk
        use glob_domain
        use grid
        use io
        use glob_ocean
        use glob_out
        use pnetcdf    , only: NF90_FLOAT
        use seaice     , only: icec, iceu, icev
        use tide       , only: define_tide, define_tide_init
     &                       , write_tide_init
        use model_run  , only: time, time_start
        use mpi        , only: MPI_OFFSET_KIND

        implicit none

        character(*), intent(in) :: out_file

!  Output floating point precision.
! Set to 4 to reduce output file's size. TODO: Unused yet
        integer, parameter :: wk = 4

        character(len=120) str_tmp!, netcdf_out_file
        integer time_dimid, x_dimid, y_dimid, z_dimid
        integer file_id, varid, status

        integer(MPI_OFFSET_KIND), dimension(4) :: start, edge


!     create file
        if ( is_master )
     &    call msg_print("", 6, "Creating `"//trim(out_file)//"`")

        file_id = file_create( out_file )

! define global attributes
        call att_write( file_id, -1, 'title'      , trim(title)   )
        call att_write( file_id, -1, 'description', 'Output file' )

! define dimensions
        time_dimid = dim_define( file_id, 'time',         0 )
        z_dimid    = dim_define( file_id, 'z'   ,        kb )
        y_dimid    = dim_define( file_id, 'y'   , jm_global )
        x_dimid    = dim_define( file_id, 'x'   , im_global )

! define variables and their attributes
        str_tmp  = 'days since '//time_start
        varid = var_define( file_id, 'time'
     &                    , NF90_FLOAT, [ time_dimid ]
     &                    , 'time'
     &                    , "days since "//time_start
     &                    , -1, 0., '' )

        varid = var_define( file_id, 'z'
     &                    , NF90_FLOAT, [ z_dimid ]
     &                    , 'sigma of cell face'
     &                    , 'sigma_level'
     &                    , -1, 0., '' )
        call att_write( file_id, varid
     &                    , 'standard_name'
     &                    , 'ocean_sigma_coordinate' )
        call att_write( file_id, varid
     &                    , 'formula_terms'
     &                    , 'sigma: z eta: elb depth: h' )

        varid = var_define( file_id, 'zz'
     &                    , NF90_FLOAT, [ z_dimid ]
     &                    , 'sigma of cell centre'
     &                    , 'sigma_level'
     &                    , -1, 0., '' )
        call att_write( file_id, varid
     &                    , 'standard_name'
     &                    , 'ocean_sigma_coordinate' )
        call att_write( file_id, varid
     &                    , 'formula_terms'
     &                    , 'sigma: zz eta: elb depth: h' )

        varid = var_define( file_id, 'east_u'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of u-points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'east_v'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of v-points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'east_e'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of elevation points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'north_u'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of u-points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'north_v'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of v-points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'north_e'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of elevation points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'rot'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'Rotation angle of x-axis wrt. east'
     &                    , 'degree'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'h'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'undisturbed water depth'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'fsm'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'free surface mask'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'dum'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'u - free surface mask'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'dvm'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'v - free surface mask'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_v north_v' )
        call define_tide_init( file_id, [ x_dimid, y_dimid ]
     &                       , NF90_FLOAT )

        varid = var_define( file_id, 'uab', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'depth-averaged u'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'vab', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'depth-averaged v'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'elb', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'surface elevation'
     &                    , 'metre'
     &                    , 0, -999., 'east_e north_e' )
        varid = var_define( file_id, 'wusurf', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'x-momentum flux'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'wvsurf', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'y-momentum flux'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'wtsurf', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'temperature flux'
     &                    , 'deg m/s'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'swrad', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'upward net heat flux'
     &                    , 'K m/s'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'wssurf', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'salinity flux'
     &                    , 'psu m/s'
     &                    , -1, 0., 'east_e north_e' )
        if ( use_ice ) then
          varid = var_define( file_id, 'icec', NF90_FLOAT
     &                      , [ x_dimid, y_dimid, time_dimid ]
     &                      , 'sea ice concentration'
     &                      , 'fraction'
     &                      , 0, 0., 'east_e north_e' )
        end if

        call define_tide( file_id, [ x_dimid, y_dimid, time_dimid ]
     &                  , NF90_FLOAT )

        if ( mode /= MODE_BAROTROPIC ) then
          varid = var_define( file_id, 'u', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'x-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_u north_u zz' )
          varid = var_define( file_id, 'v', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'y-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_v north_v zz' )
          varid = var_define( file_id, 'w', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'z-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_e north_e z' )
          varid = var_define( file_id, 't', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'potential temperature'
     &                      , 'degC'
     &                      , 0, -999., 'east_e north_e zz' )
          varid = var_define( file_id, 's', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'salinity x rho / rhoref'
     &                      , 'PSS'
     &                      , 0, 0., 'east_e north_e zz' )
          varid = var_define( file_id, 'rho', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , '(density-1000)/rhoref'
     &                      , 'dimensionless'
     &                      , -1, 0., 'east_e north_e zz' )
        end if

! end definitions
        call file_end_definition( file_id )

! write static data
        start = [  1, 1,1,1 ]
        edge  = [ kb, 1,1,1 ]
        call var_write( file_id, "z" , z , start, edge )
        call var_write( file_id, "zz", zz, start, edge )

        start(1:2) = [ i_global(1), j_global(1) ]
        edge (1:2) = [  im        , jm          ]
        call var_write( file_id, "east_u" , east_u , start, edge )
        call var_write( file_id, "east_v" , east_v , start, edge )
        call var_write( file_id, "east_e" , east_e , start, edge )
        call var_write( file_id, "north_u", north_u, start, edge )
        call var_write( file_id, "north_v", north_v, start, edge )
        call var_write( file_id, "north_e", north_e, start, edge )
        call var_write( file_id, "rot"    , rot    , start, edge )
        call var_write( file_id, "h"      , h      , start, edge )
        call var_write( file_id, "fsm"    , fsm    , start, edge )
        call var_write( file_id, "dum"    , dum    , start, edge )
        call var_write( file_id, "dvm"    , dvm    , start, edge )
        call write_tide_init( file_id, start, edge )

        file_id = file_close( file_id )


      end ! subroutine create_output


!______________________________________________________________________
!
      subroutine write_output( out_file, record, write_means )
!----------------------------------------------------------------------
!  Write output file.
!______________________________________________________________________
!
        use air
        use config     , only: use_ice, mode, title ! TODO: `use_ice` move output sections to be perfomed by their own modules
        use glob_const , only: MODE_BAROTROPIC, rk
        use glob_domain
        use grid
        use io
        use glob_ocean
        use pnetcdf    , only: NF90_WRITE
        use seaice     , only: icec, iceu, icev
        use tide       , only: write_tide
        use model_run  , only: time, time_start
        use mpi        , only: MPI_OFFSET_KIND
        use glob_out

        implicit none

        logical     , intent(in) :: write_means
        character(*), intent(in) :: out_file
        integer     , intent(in) :: record

        integer                                   file_id , status
        integer(MPI_OFFSET_KIND), dimension(4) :: start   , edge
        character(PATH_LEN)                       tmp_str


        tmp_str = ""
!     create file
        if ( is_master ) then
          write( tmp_str, '("Writing `",a,"` @ ",i5)' )
     &           trim(out_file), record
          call msg_print("", 6, tmp_str)
        end if

        file_id = file_open( out_file, NF90_WRITE ) ! TODO: Eliminate pnetcdf dependency. Pass own type ids and the responsible module will handle possible conversion.

        start = [ record, 1,1,1 ]
        edge  = [      1, 1,1,1 ]
        call var_write( file_id, "time", time, start )

        start = [ i_global(1), j_global(1), record, 1 ]
        edge  = [  im        , jm         ,      1, 1 ]
        if ( write_means ) then
          call var_write( file_id, "uab"   , uab_mean   , start, edge )
          call var_write( file_id, "vab"   , vab_mean   , start, edge )
          call var_write( file_id, "elb"   , elb_mean   , start, edge )
          call var_write( file_id, "wusurf", wusurf_mean, start, edge )
          call var_write( file_id, "wvsurf", wvsurf_mean, start, edge )
          call var_write( file_id, "wtsurf", wtsurf_mean, start, edge )
          call var_write( file_id, "wssurf", wssurf_mean, start, edge )
          call var_write( file_id, "swrad" , swrad_mean , start, edge )
        else
          call var_write( file_id, "uab"   , uab   , start, edge )
          call var_write( file_id, "vab"   , vab   , start, edge )
          call var_write( file_id, "elb"   , elb   , start, edge )
          call var_write( file_id, "wusurf", wusurf, start, edge )
          call var_write( file_id, "wvsurf", wvsurf, start, edge )
          call var_write( file_id, "wtsurf", wtsurf, start, edge )
          call var_write( file_id, "wssurf", wssurf, start, edge )
          call var_write( file_id, "swrad" , swrad , start, edge )
        end if
        if ( use_ice ) then
          call var_write( file_id, "icec", icec  , start, edge )
        end if
        call write_tide( file_id, start, edge )

        if ( mode /= MODE_BAROTROPIC ) then
          start(3:4) = [  1, record ]
          edge (3:4) = [ kb,      1 ]
          if ( write_means ) then
            call var_write( file_id, "u"   , u_mean  , start, edge )
            call var_write( file_id, "v"   , v_mean  , start, edge )
            call var_write( file_id, "w"   , w_mean  , start, edge )
            call var_write( file_id, "s"   , s_mean  , start, edge )
            call var_write( file_id, "t"   , t_mean  , start, edge )
            call var_write( file_id, "rho" , rho_mean, start, edge )
          else
            call var_write( file_id, "u"   , u  , start, edge )
            call var_write( file_id, "v"   , v  , start, edge )
            call var_write( file_id, "w"   , w  , start, edge )
            call var_write( file_id, "s"   , s  , start, edge )
            call var_write( file_id, "t"   , t  , start, edge )
            call var_write( file_id, "rho" , rho, start, edge )
          end if
        end if

!     close file
        file_id = file_close( file_id )


      end ! subroutine write_output
!______________________________________________________________________
!
      subroutine create_initial( out_file )
!----------------------------------------------------------------------
!  Create initial conditions file.
!______________________________________________________________________
!
        use air
        use clim       , only: rmean
        use config     , only: use_ice, mode, title
        use glob_const , only: MODE_BAROTROPIC, rk
        use glob_domain
        use grid
        use io
        use glob_ocean
        use glob_out
        use pnetcdf    , only: NF90_FLOAT
        use seaice     , only: icec, iceu, icev
        use model_run  , only: time, time_start
        use mpi        , only: MPI_OFFSET_KIND

        implicit none

        character(*), intent(in) :: out_file

!  Output floating point precision.
! Set to 4 to reduce output file's size. TODO: Unused yet
        integer, parameter :: wk = 4

        character(len=120) str_tmp!, netcdf_out_file
        integer time_dimid, x_dimid, y_dimid, z_dimid
        integer file_id, varid, status

        integer(MPI_OFFSET_KIND), dimension(4) :: start, edge


!     create file
        if ( is_master )
     &    call msg_print("", 6, "Creating `"//trim(out_file)//"`")

        file_id = file_create( out_file )

! define global attributes
        call att_write( file_id, -1, 'title'      , trim(title)   )
        call att_write( file_id, -1, 'description', 'Initial file' )

! define dimensions
        time_dimid = dim_define( file_id, 'time',         1 )
        z_dimid    = dim_define( file_id, 'z'   ,        kb )
        y_dimid    = dim_define( file_id, 'y'   , jm_global )
        x_dimid    = dim_define( file_id, 'x'   , im_global )

! define variables and their attributes
        str_tmp  = 'days since '//time_start
        varid = var_define( file_id, 'time'
     &                    , NF90_FLOAT, [ time_dimid ]
     &                    , 'time'
     &                    , "days since "//time_start
     &                    , -1, 0., '' )

        varid = var_define( file_id, 'z'
     &                    , NF90_FLOAT, [ z_dimid ]
     &                    , 'sigma of cell face'
     &                    , 'sigma_level'
     &                    , -1, 0., '' )
        call att_write( file_id, varid
     &                    , 'standard_name'
     &                    , 'ocean_sigma_coordinate' )
        call att_write( file_id, varid
     &                    , 'formula_terms'
     &                    , 'sigma: z eta: elb depth: h' )

        varid = var_define( file_id, 'zz'
     &                    , NF90_FLOAT, [ z_dimid ]
     &                    , 'sigma of cell centre'
     &                    , 'sigma_level'
     &                    , -1, 0., '' )
        call att_write( file_id, varid
     &                    , 'standard_name'
     &                    , 'ocean_sigma_coordinate' )
        call att_write( file_id, varid
     &                    , 'formula_terms'
     &                    , 'sigma: zz eta: elb depth: h' )

        varid = var_define( file_id, 'east_u'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of u-points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'east_v'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of v-points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'east_e'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'easting of elevation points'
     &                    , 'degrees east'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'north_u'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of u-points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'north_v'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of v-points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'north_e'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'northing of elevation points'
     &                    , 'degrees north'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'rot'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'Rotation angle of x-axis wrt. east'
     &                    , 'degree'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'h'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'undisturbed water depth'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'fsm'
     &                    , NF90_FLOAT, [ x_dimid, y_dimid ]
     &                    , 'free surface mask'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_e north_e' )

        varid = var_define( file_id, 'zeta', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'surface elevation'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'uwnd', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'meridional wind'
     &                    , 'metre/sec'
     &                    , 0, -999., 'east_e north_e' )
        varid = var_define( file_id, 'vwnd', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'zonal wind'
     &                    , 'metre/sec'
     &                    , 0, -999., 'east_e north_e' )
        varid = var_define( file_id, 'swrad', NF90_FLOAT
     &                    , [ x_dimid, y_dimid, time_dimid ]
     &                    , 'net shortwave radiation upwards'
     &                    , 'K m/s'
     &                    , -1, 0., 'east_e north_e' )

        if ( mode /= MODE_BAROTROPIC ) then
          varid = var_define( file_id, 'u', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'x-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_u north_u zz' )
          varid = var_define( file_id, 'v', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'y-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_v north_v zz' )
          varid = var_define( file_id, 'w', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'z-velocity'
     &                      , 'metre/sec'
     &                      , -1, 0., 'east_e north_e z' )
          varid = var_define( file_id, 't', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'potential temperature'
     &                      , 'degC'
     &                      , 0, -999., 'east_e north_e zz' )
          varid = var_define( file_id, 's', NF90_FLOAT
     &                      , [ x_dimid, y_dimid
     &                        , z_dimid, time_dimid ]
     &                      , 'salinity x rho / rhoref'
     &                      , 'PSS'
     &                      , 0, 0., 'east_e north_e zz' )
        end if

        varid = var_define( file_id, 'rho', NF90_FLOAT
     &                    , [ x_dimid, y_dimid
     &                      , z_dimid, time_dimid ]
     &                    , '(density-1000)/rhoref'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'rmean', NF90_FLOAT
     &                    , [ x_dimid, y_dimid
     &                      , z_dimid, time_dimid ]
     &                    , 'Background (density-1000)/rhoref'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_e north_e zz' )

! end definitions
        call file_end_definition( file_id )

! write static data
        start = [  1, 1,1,1 ]
        edge  = [ kb, 1,1,1 ]
        call var_write( file_id, "z" , z , start, edge )
        call var_write( file_id, "zz", zz, start, edge )

        start(1:2) = [ i_global(1), j_global(1) ]
        edge (1:2) = [  im        , jm          ]
        call var_write( file_id, "east_u" , east_u , start, edge )
        call var_write( file_id, "east_v" , east_v , start, edge )
        call var_write( file_id, "east_e" , east_e , start, edge )
        call var_write( file_id, "north_u", north_u, start, edge )
        call var_write( file_id, "north_v", north_v, start, edge )
        call var_write( file_id, "north_e", north_e, start, edge )
        call var_write( file_id, "rot"    , rot    , start, edge )
        call var_write( file_id, "h"      , h      , start, edge )
        call var_write( file_id, "fsm"    , fsm    , start, edge )
        call var_write( file_id, "zeta"   , elb    , start, edge )
        call var_write( file_id, "uwnd"   , uwsrf  , start, edge )
        call var_write( file_id, "vwnd"   , vwsrf  , start, edge )
        call var_write( file_id, "swrad"  , swrad  , start, edge )

        edge(3) = kb

        if ( mode /= MODE_BAROTROPIC ) then
          call var_write( file_id, "u"      , ub     , start, edge )
          call var_write( file_id, "v"      , vb     , start, edge )
          call var_write( file_id, "t"      , tb     , start, edge )
          call var_write( file_id, "s"      , sb     , start, edge )
          call var_write( file_id, "w"      , w      , start, edge )
        end if

        call var_write( file_id, "rho"    , rho    , start, edge )
        call var_write( file_id, "rmean"  , rmean  , start, edge )

        file_id = file_close( file_id )


      end ! subroutine create_initial
!______________________________________________________________________
!
      subroutine create_restart( out_file )
!----------------------------------------------------------------------
!  Create restart file.
!______________________________________________________________________
!
        use air
        use config     , only: use_ice, mode, spinup, title
        use glob_const , only: MODE_BAROTROPIC, rk
        use glob_domain
        use grid
        use io
        use glob_ocean
        use glob_out
        use pnetcdf    , only: NF90_DOUBLE
        use seaice     , only: icec, iceu, icev
        use model_run  , only: time, time_start, time0
        use mpi        , only: MPI_OFFSET_KIND

        implicit none

        character(*), intent(in) :: out_file

        integer time_dimid, x_dimid, y_dimid, z_dimid
        integer file_id, varid, status

        integer(MPI_OFFSET_KIND), dimension(4) :: start, edge


!     create file
        if ( is_master )
     &    call msg_print("", 6, "Creating `"//trim(out_file)//"`")

        file_id = file_create( out_file )

! define global attributes
        call att_write( file_id, -1, 'title'      , trim(title)   )
        call att_write( file_id, -1, 'description', 'Restart file' )

! define dimensions
        time_dimid = dim_define( file_id, 'time',         1 )
        z_dimid    = dim_define( file_id, 'z'   ,        kb )
        y_dimid    = dim_define( file_id, 'y'   , jm_global )
        x_dimid    = dim_define( file_id, 'x'   , im_global )

! define variables and their attributes
        varid = var_define( file_id, 'time'
     &                    , NF90_DOUBLE, [ time_dimid ]
     &                    , 'time'
     &                    , "days since "//time_start
     &                    , -1, 0., '' )

        varid = var_define( file_id, 'wusurf'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'x-momentum flux at the surface'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'wvsurf'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'y-momentum flux at the surface'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'wubot'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'x-momentum flux at the bottom'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'wvbot'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'y-momentum flux at the bottom'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'aam2d'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical average of aam'
     &                    , 'metre^2/sec'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'ua'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical mean of u'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'uab'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical mean of u at time -dt'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'va'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical mean of v'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'vab'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical mean of v at time -dt'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'el'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'surface elevation in external mode'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'elb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'surface elevation in external mode at -dt'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'et'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'surface elevation in internal mode'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'etb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'surface elevation in internal mode at -dt'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'egb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'surface elevation for pres. grad. at -dt'
     &                    , 'metre'
     &                    , -1, 0., 'east_e north_e' )
        varid = var_define( file_id, 'utb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'ua time averaged over dti'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'vtb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'va time averaged over dti'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'adx2d'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical integral of advx'
     &                    , '-'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'ady2d'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'vertical integral of advy'
     &                    , '-'
     &                    , -1, 0., 'east_v north_v' )
        varid = var_define( file_id, 'advua'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'sum of 2nd, 3rd and 4th terms in eq (18)'
     &                    , '-'
     &                    , -1, 0., 'east_u north_u' )
        varid = var_define( file_id, 'advva'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid ]
     &                    , 'sum of 2nd, 3rd and 4th terms in eq (19)'
     &                    , '-'
     &                    , -1, 0., 'east_v north_v' )

        varid = var_define( file_id, 'u'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'x-velocity'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u zz' )
        varid = var_define( file_id, 'ub'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'x-velocity at time -dt'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_u north_u zz' )
        varid = var_define( file_id, 'v'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'y-velocity'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v zz' )
        varid = var_define( file_id, 'vb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'y-velocity at time -dt'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_v north_v zz' )
        varid = var_define( file_id, 'w'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'sigma-velocity'
     &                    , 'metre/sec'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 't'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'potential temperature'
     &                    , 'K'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'tb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'potential temperature at time -dt'
     &                    , 'K'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 's'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'salinity x rho / rhoref'
     &                    , 'PSS'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'sb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'salinity x rho / rhoref at time -dt'
     &                    , 'PSS'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'rho'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , '(density - 1000) / rhoref'
     &                    , 'dimensionless'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'km'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'vertical kinematic viscosity'
     &                    , 'metre^2/sec'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'kh'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'vertical diffusivity'
     &                    , 'metre^2/sec'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'kq'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'kq'
     &                    , 'metre^2/sec'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'l'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'turbulent length scale'
     &                    , '-'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'q2'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'twice the turbulent kinetic energy'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'q2b'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'twice the turbulent kinetic energy at -dt'
     &                    , 'metre^2/sec^2'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'aam'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'horizontal kinematic viscosity'
     &                    , 'metre^2/sec'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'q2l'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'q2 x l'
     &                    , 'metre^3/sec^2'
     &                    , -1, 0., 'east_e north_e zz' )
        varid = var_define( file_id, 'q2lb'
     &                    , NF90_DOUBLE, [ x_dimid, y_dimid, z_dimid ]
     &                    , 'q2 x l at time -dt'
     &                    , 'metre^3/sec^2'
     &                    , -1, 0., 'east_e north_e zz' )

! end definitions
        call file_end_definition( file_id )

! write data
        start = [ i_global(1), j_global(1),  1, 1 ]
        edge  = [ im         , jm         , kb, 1 ]

        if ( spinup ) then
          call var_write( file_id, "time", time0, start(4:4) )
        else
          call var_write( file_id, "time", time , start(4:4) )
        end if

        call var_write( file_id, "wusurf" , wusurf , start, edge )
        call var_write( file_id, "wvsurf" , wvsurf , start, edge )
        call var_write( file_id, "wubot"  , wubot  , start, edge )
        call var_write( file_id, "wvbot"  , wvbot  , start, edge )
        call var_write( file_id, "aam2d"  , aam2d  , start, edge )
        call var_write( file_id, "ua"     , ua     , start, edge )
        call var_write( file_id, "uab"    , uab    , start, edge )
        call var_write( file_id, "va"     , va     , start, edge )
        call var_write( file_id, "vab"    , vab    , start, edge )
        call var_write( file_id, "el"     , el     , start, edge )
        call var_write( file_id, "elb"    , elb    , start, edge )
        call var_write( file_id, "et"     , et     , start, edge )
        call var_write( file_id, "etb"    , etb    , start, edge )
        call var_write( file_id, "egb"    , egb    , start, edge )
        call var_write( file_id, "utb"    , utb    , start, edge )
        call var_write( file_id, "vtb"    , vtb    , start, edge )
        call var_write( file_id, "adx2d"  , adx2d  , start, edge )
        call var_write( file_id, "ady2d"  , ady2d  , start, edge )
        call var_write( file_id, "advua"  , advua  , start, edge )
        call var_write( file_id, "advva"  , advva  , start, edge )

        call var_write( file_id, "u"      , u      , start, edge )
        call var_write( file_id, "ub"     , ub     , start, edge )
        call var_write( file_id, "v"      , v      , start, edge )
        call var_write( file_id, "vb"     , vb     , start, edge )
        call var_write( file_id, "w"      , w      , start, edge )
        call var_write( file_id, "t"      , t      , start, edge )
        call var_write( file_id, "tb"     , tb     , start, edge )
        call var_write( file_id, "s"      , s      , start, edge )
        call var_write( file_id, "sb"     , sb     , start, edge )
        call var_write( file_id, "rho"    , rho    , start, edge )
        call var_write( file_id, "km"     , km     , start, edge )
        call var_write( file_id, "kh"     , kh     , start, edge )
        call var_write( file_id, "kq"     , kq     , start, edge )
        call var_write( file_id, "l"      , l      , start, edge )
        call var_write( file_id, "q2"     , q2     , start, edge )
        call var_write( file_id, "q2b"    , q2b    , start, edge )
        call var_write( file_id, "aam"    , aam    , start, edge )
        call var_write( file_id, "q2l"    , q2l    , start, edge )
        call var_write( file_id, "q2lb"   , q2lb   , start, edge )

        file_id = file_close( file_id )


      end ! subroutine create_restart

!
!______________________________________________________________________
!lyo:pac10:beg:Here thro *end: replaced subr.lateral_boundary_conditions
      subroutine lateral_boundary_conditions
!----------------------------------------------------------------------
!  Read lateral boundary conditions.
!  Transport at eastern boundary for PROFS in GOM.
! ayumi 2010/4/7
!______________________________________________________________________
!
      use bry        , only: periodic_bc
      use glob_const , only: rk
      use glob_domain, only: i_global, im, im_global
     &                     , j_global, jm, jm_global, master_task
      use grid       , only: cor

      implicit none

      logical :: here, judge_inout !lyo:scs1d:
      integer :: i,j,ic,jc         !lyo:scs1d:
      real(kind=rk)    :: corcon            !lyo:scs1d:
!      character(len=120) in_file        !eda:uvforce

!     call read_uabe_pnetcdf(uabe)
!eda:uvforce
!      if (.not. calc_uvforce) then
!        write(in_file,'(a)') "bc.nc"
!        call read_bc_pnetcdf(uabe, uabw, vabs, vabn, in_file, 1)
!      endif

! Periodic in "x" and/or "y"?  !lyo:20110224:alu:stcc:
!     iperx.ne.0 if x-periodic; ipery.ne.0 if y-periodic               !
!     iperx(y) < 0 if south/north (west/east) walls are free-slip      !
!     iperx=-1; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!     iperx= 0; ipery= 0 !--> x-periodic & free-slip S/N boundaries
!lyo:scs1d:moved to namelist: pom.nml
!lyo:scs1d:cannot be beta-plane if double-periodic (note cor(*,*)
!     was previously defined in "call read_grid")
      if ( periodic_bc%x .and. periodic_bc%y ) then
      ic = (im_global+1)/2; jc = (jm_global+1)/2
         here = judge_inout( ic, jc,
     $                       i_global(1), i_global(im),
     $                       j_global(1), j_global(jm) )
         if ( here ) then
            corcon=cor(ic-i_global(1)+1,jc-j_global(1)+1)
         endif
         call bcast0d_mpi(corcon,1,master_task)
         do j=1,jm; do i=1,im !f-plane:
            cor(i,j)=corcon
            enddo; enddo
         endif !if (iperx.eq.1 .and. ipery.eq.1) then


      end
!
!______________________________________________________________________
!
      subroutine bfrz(mw,me,ms,mn,nw,ne,ns,nn,im,jm,nu,frz)
!----------------------------------------------------------------------
!lyo:20110224:alu:stcc:
!lyo:Modified from subroutine assimfrs in
!     /wrk/newshoni/hunglu/model/sbPOM/stcc_ideal/
!     stcc_alu_30TSrelx_60aam_timescle1d/pom/pom.f
!     which was modified from:
!     /archive/lyo/arctic/beringsea/codes_and_runscripts/
!     bering07_hass50m_faccof1p2_dtfac.f; see also
!     /wrk/aden/lyo/crwu/assimssh/how_to_put_assimssh_in_pom.txt
!lyo:Modified from /home/lyo/gom/hindcast/fgrid/gomn117.f subr.NFRSASIM
!----------------------------------------------------------------------
!     calculate boundary flow-relaxation array "frz"
!----------------------------------------------------------------------
!     nu          = unit# (fort.nu) for ASCII printout
!     mw,me,ms,mn = #buffer grid-pnts near west, east, south & north
!                   boundaries where assimilation frz=1.0; recommended
!                   buffer 100~200km;  m?=0 for no buffer;
!                   As a precaution, program stops if 0<m?<=3 (or <0)
!     nw,ne,ns,nn = n_west,n_east,n_south,n_north
!                   is =-1 if "bfrz" is being called by processor that
!                   shares the west, east, south or north boundary
!                   Just set all to -1 for a single-processor run
!     frz         = 1 at boundaries and tapered (tanh) =0  interior
!
!                ... l.oey --- lyo@princeton.edu (Jan/2008)
!----------------------------------------------------------------------

      use glob_const , only: rk

      implicit none

      integer, intent(in) :: mw,me,ms,mn,im,jm,nu
      integer, intent(in) :: nw,ne,ns,nn
      real(kind=rk), dimension(im,jm), intent(out) :: frz
      integer :: mmid,i,j,ii,jj
      integer, parameter :: my_task=0
      real(kind=rk) :: c,tanhm


      frz(:,:)=0.0       !Initialize interior

!     West:
      if(nw.eq.-1) then
       if (mw.gt.3) then  !west buffer: needs at least 4 pts
         frz(1,:)=1.0; mmid=mw/2; c=5./float(mw-mmid)
         do i=2,mw
           ii=i
           tanhm=0.5*(1.-tanh(float(i-mmid)*c))
           frz(ii,:)=tanhm
         enddo
       elseif(mw.eq.0.or.mw.eq.1) then
!      do nothing, i.e. frz remains = 0 for mw=0 or 1
       else !mw=2 or 3 or <0
        write(nu,'(''Stopped in subr.bfrz, mw ='',i4)') mw
        write( *,'(''Stopped in bfrz. proc# ,mw ='',2i4)') my_task,mw
        stop
       endif
      endif
!
!     East:
      if(ne.eq.-1) then
       if (me.gt.3) then  !east buffer:
         frz(im,:)=1.0; mmid=me/2; c=5./float(me-mmid)
         do i=2,me
           ii=im-i+1
           tanhm=0.5*(1.-tanh(float(i-mmid)*c))
           frz(ii,:)=tanhm
         enddo
       elseif(me.eq.0.or.me.eq.1) then
!      do nothing, i.e. frz remains = 0 for me=0 or 1
       else !me=2 or 3 or <0
        write(nu,'(''Stopped in subr.bfrz, me ='',i4)') me
        write( *,'(''Stopped in bfrz. proc# ,me ='',2i4)') my_task,me
        stop
       endif
      endif
!
!     South:
      if(ns.eq.-1) then
       if (ms.gt.3) then  !south buffer:
         frz(:,1)=1.0; mmid=ms/2; c=5./float(ms-mmid)
         do j=2,ms
           jj=j
           tanhm=0.5*(1.-tanh(float(j-mmid)*c))
           do i=1,im !lyo:debug:lyo:20110224:fxu:stcc:delete "if (nw.eq.-1.."
           frz(i,jj)=max(tanhm,frz(i,jj)) !takes care of SW & SE corners
           enddo
         enddo
       elseif(ms.eq.0.or.ms.eq.1) then
!      do nothing, i.e. frz remains = 0 for ms=0 or 1
       else
        write(nu,'(''Stopped in subr.bfrz, ms ='',i4)') ms
        write( *,'(''Stopped in bfrz. proc# ,ms ='',2i4)') my_task,ms
        stop
       endif
      endif
!
!     North:
      if(nn.eq.-1) then
       if (mn.gt.3) then  !north buffer:
         frz(:,jm)=1.0; mmid=mn/2; c=5./float(mn-mmid)
         do j=2,mn
           jj=jm-j+1
           tanhm=0.5*(1.-tanh(float(j-mmid)*c))
           do i=1,im !lyo:debug:lyo:20110224:fxu:stcc:delete "if (nw.eq.-1.."
           frz(i,jj)=max(tanhm,frz(i,jj)) !takes care of NW & NE corners
           enddo
         enddo
       elseif(mn.eq.0.or.mn.eq.1) then
!      do nothing, i.e. frz remains = 0 for mn=0 or 1
       else
        write(nu,'(''Stopped in subr.bfrz, mn ='',i4)') mn
        write( *,'(''Stopped in bfrz. proc# ,mn ='',2i4)') my_task,mn
        stop
       endif
      endif


      end
!lyo:pac10:end:
!
!______________________________________________________________________
!
      subroutine read_tide
!----------------------------------------------------------------------
!fhx:tide:read tidal amplitude & phase at the eastern boundary for PROFS
!______________________________________________________________________
!
        use glob_const , only: rk


!      call read_tide_east_pnetcdf(ampe,phae)
!        call read_tide_east_pnetcdf(ampe,phae,amue,phue)

!      if(my_task.eq.0)print*,ampe(10,1),phae(10,1),amue(10,1),phue(10,1)
!      if(my_task.eq.0)print*,ampe(10,2),phae(10,2),amue(10,2),phue(10,2)


      end
!fhx:tide:read_tide end
!
!______________________________________________________________________
!
      subroutine update_initial
!----------------------------------------------------------------------
!  Updates the initial conditions and sets the remaining
! initial conditions.
!______________________________________________________________________
!
      use air        , only: vfluxf
      use config     , only: aam_init, npg  ,do_restart,use_tide
      use glob_const , only: rk, SMALL
      use glob_domain, only: im, is_master, jm, kb, kbm1  ,my_task
      use grid       , only: dz, h
      use model_run,only:iint, dtime
      use glob_ocean , only: aam, d, drhox, drhoy, drx2d, dry2d, dt
     &                     , el, elb, et, etb, etf
     &                     , kh, km, kq, l, q2, q2b, q2l, q2lb
     &                     , rho, s, sb, t, tb, u, ua, ub, uab
     &                     , v, va, vb, vab, w
      use tide, only: tide_advance=>step, tide_el, tide_ua, tide_va

      implicit none

      integer i,j,k


      if ( .not.do_restart .and. use_tide ) then
        call tide_advance( dtime )
        uab = tide_ua
        vab = tide_va
        elb = tide_el
        etb = elb
      end if
      ua = uab
      va = vab
      el = elb
      et = etb
      etf= et
      w(:,:,1) = vfluxf

      d  = h + el
      dt = h + et

      do k=1,kb
        l(:,:,k) = .1*dt
      end do

      q2b = SMALL
      q2lb= l*q2b
      kh  = l*sqrt(q2b)
      km  = kh
      kq  = kh
      aam = aam_init

      do k=1,kbm1
        do i=1,im
          do j=1,jm
            q2(i,j,k)=q2b(i,j,k)
            q2l(i,j,k)=q2lb(i,j,k)
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            u(i,j,k)=ub(i,j,k)
            v(i,j,k)=vb(i,j,k)
          end do
        end do
      end do

      if ( is_master ) then
        if ( npg>5 .and. npg<0 ) then
          print '(/"[!] Invalid value for Pressure Gradient Scheme")'
          print '("   Defaulting to McCalpin 4th order (npg=2)")'
        end if
      end if

      call dens(s,t,rho)

      call pgscheme(npg)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do


      end
!______________________________________________________________________
!
      subroutine initialize_modules
!----------------------------------------------------------------------
!  Initialize all neccessary modules.
!______________________________________________________________________
!
      use air   , only: initialize_air => initialize_mod
      use bry   , only: initialize_bry => initialize_mod
      use clim  , only: initialize_clm => initialize_mod
      use seaice, only: initialize_ice => initialize_mod
      use tide  , only: initialize_tide=> initialize_mod

      implicit none


      call initialize_clm( 'config.nml' )
      call initialize_bry( 'config.nml' )
      call initialize_air( 'config.nml' )
      call initialize_ice( 'config.nml' )
      call initialize_tide('config.nml' )


      end ! subroutine initialize_modules
!
!______________________________________________________________________
!
      subroutine bottom_friction
!----------------------------------------------------------------------
!  Calculates the bottom friction coefficient.
!______________________________________________________________________
!
      use config     , only: cbcmax, cbcmin, z0b
      use glob_const , only: Kappa, rk
      use glob_domain, only: im, jm, kbm1
      use grid       , only: h, zz
      use glob_ocean , only: cbc

      implicit none

      integer i,j


! calculate bottom friction
      do j=1,jm
        do i=1,im
!lyo:correct:cbc(i,j)=(Kappa/log((1.+zz(kbm1))*h(i,j)/z0b))**2 !lyo:bug:
          cbc(i,j)=(Kappa/log(1.+(1.0+zz(kbm1))*h(i,j)/z0b))**2
          cbc(i,j)=max(cbcmin,cbc(i,j))
! if the following is invoked, then it is probable that the wrong
! choice of z0b or vertical spacing has been made:
          cbc(i,j)=min(cbcmax,cbc(i,j))
        end do
      end do


      end ! subroutine bottom_friction
!
!______________________________________________________________________
!
      subroutine ztosig(zs,tb,zz,h,t,im,jm,ks,kb,
     $                                   n_west,n_east,n_south,n_north)
!----------------------------------------------------------------------
!  Interpolates vertically from z- to sigma-coordinates.
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer , intent(in ) :: im, jm, ks, kb
      real(rk), intent(in ) :: zs(ks), tb(im,jm,ks), zz(kb)
     &                       , h(im,jm)
      real(rk), intent(out) :: t(im,jm,kb)
      integer , intent(in ) :: n_west,n_east,n_south,n_north

      integer  i,j,k
      real(rk) tmax, tin(ks), tout(kb), zzh(kb)


      do j=2,jm-1
        do i=2,im-1
          if ( h(i,j) > 1. ) then
! special interp on z-lev for cases of no data because h smoothing
            do k=1,ks
              tin(k) = tb(i,j,k)
              if ( zs(k)<=h(i,j) .and. tin(k)<.01 ) then
                tmax = max(tb(i-1,j,k),tb(i+1,j,k),
     &                     tb(i,j-1,k),tb(i,j+1,k))
                tin(k) = tmax
              end if
              if ( tin(k)<.01 .and. k/=1 ) tin(k) = tin(k-1)
            end do

            do k=1,kb
              zzh(k) = -zz(k)*h(i,j)
            end do

! vertical spline interp
            call splinc(zs,tin,ks,2.e30_rk,2.e30_rk,zzh,tout,kb)

            do k=1,kb
              t(i,j,k) = tout(k)
            end do

          end if
        end do
      end do
      call exchange3d_mpi(t,im,jm,kb)

! boundaries
      do k=1,kb
        do j=1,jm
          if ( n_west == -1 ) t( 1,j,k) = t(   2,j,k)
          if ( n_east == -1 ) t(im,j,k) = t(im-1,j,k)
        end do
        do i=1,im
          if ( n_south == -1 ) t(i, 1,k) = t(i,   2,k)
          if ( n_north == -1 ) t(i,jm,k) = t(i,jm-1,k)
        end do
      end do


      end
!
!______________________________________________________________________
!
      subroutine splinc(x,y,n,yp1,ypn,xnew,ynew,m)
!----------------------------------------------------------------------
!  Interpolate using splines.
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer, parameter :: nmax = 300

      integer               , intent(in ) :: n, m
      real(rk)              , intent(in ) :: yp1, ypn
      real(rk), dimension(n), intent(in ) :: x,y
      real(rk), dimension(m), intent(out) :: xnew,ynew

      integer                      i, k
      real(rk)                     p, qn, sig, un
      real(rk), dimension(nmax) :: y2,u


      if ( yp1 > .99e30 ) then
        y2(1) = 0.
        u(1)  = 0.
      else
        y2(1) = -.5
        u(1)  = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if

      do i=2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.
        y2(i) = (sig-1.)/p
        u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if ( ypn > .99e30 ) then
        qn = 0.
        un = 0.
      else
        qn =  .5
        un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if

      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
      end do

      do i=1,m
        call splint(x,y,y2,n,xnew(i),ynew(i))
      end do


      end
!
!______________________________________________________________________
!
      subroutine splint(xa,ya,y2a,n,x,y)
!----------------------------------------------------------------------
!  Spline interpolation? [ DESCRIPTION NOT PROVIDED ]
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: error_status

      implicit none

      integer               , intent(in ) :: n
      real(rk), dimension(n), intent(in ) :: xa, ya, y2a
      real(rk)              , intent(out) :: x, y

      integer k, khi, klo
      real(rk) a, b, h


      klo = 1
      khi = n
      do while ( khi-klo > 1 )
        k = (khi+klo)/2
        if ( xa(k) > x ) then
          khi = k
        else
          klo = k
        end if
      end do
      h = xa(khi) - xa(klo)
      if ( h == 0. )  then
        error_status = 1
        print '(/a)', 'Error: bad xa input in splint'
      end if
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+
     &       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.


      end ! subroutine splint
!
!______________________________________________________________________
!fhx:incmix:add subroutine incmix
      subroutine incmix(aam,im,jm,kb,lono,lato,x,y,xs,ys,fac)
!----------------------------------------------------------------------
!  Increase 'aam' by '(1+fac)*' at (lono,lato), then taper off to '1'
! in gaussian-manner for (x-xo, y-yo) > xs and ys.
!
!     Inputs: aam,x,y,xs,ys & fac
!     Output: aam is modified
!______________________________________________________________________
!
      use glob_const , only: rk

      implicit none

      integer                      , intent(in   ) :: im, jm, kb
      real(rk)                     , intent(in   ) :: lono, lato
     &                                              , xs, ys, fac
      real(rk), dimension(im,jm,kb), intent(inout) :: aam
      real(rk), dimension(im,jm   ), intent(in   ) :: x, y

      integer  i,j,k
      real(rk) factor,expon


!      print*,'incmix:', lono,lato
      do k=1,kb
        do j=1,jm
          do i=1,im
            factor=0.0
!      expon=((x(i,j)-x(io,jo))/xs)**2+((y(i,j)-y(io,jo))/ys)**2
            expon = ((x(i,j)-lono)/xs)**2+((y(i,j)-lato)/ys)**2
            if ( expon <= 10. ) then
              factor = fac*exp(-expon)
            end if
            aam(i,j,k) = aam(i,j,k)*(1.+factor)
          end do
        end do
      end do


      end
!lyo:pac10:exp013:
!______________________________________________________________________
!
      pure logical function judge_inout( i_in, j_in,
     &                                   imin_in, imax_in,
     &                                   jmin_in, jmax_in )
!----------------------------------------------------------------------
!  If the processor has ( i_in, j_in ) in its local domain.
!______________________________________________________________________
!
      implicit none

      integer, intent(in) ::    i_in,    j_in
     &                     , imin_in, imax_in
     &                     , jmin_in, jmax_in


      if (       ( i_in >= imin_in .and. i_in <= imax_in )
     &     .and. ( j_in >= jmin_in .and. j_in <= jmax_in ) ) then

         judge_inout = .true.

      else

         judge_inout = .false.

      endif


      end ! function judge_inout
!
!______________________________________________________________________
!
      subroutine make_grid( TEST_CASE )
!----------------------------------------------------------------------
!  Generate vertical and horizontal grid, topography, areas and masks.
!______________________________________________________________________
!
      use glob_const , only: DEG2RAD, rk
      use glob_domain, only: i_global, im_global, is_master
     &                     , j_global, jm_global, kb
     &                     , n_east, n_north, n_south, n_west
      use grid
      use glob_ocean , only: d, dt, el, et

      implicit none

      integer, intent(in) :: TEST_CASE

      integer i,j,k

      if ( is_master ) then
        call msg_print("IDEALIZED CASE: "//char(40+TEST_CASE),1,"")
      end if 

! generate grid
      do k = 1,kb
        z(k)  = -real(k-1)/real(kb-1)
      end do
      do k = 1,kb-1
        zz(k) = .5*(z(k+1)+z(k))
      end do
      zz(kb) = 2.*zz(kb-1)-zz(kb-2)

      do k=1,kb-1
        dz(k) = z(k)- z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb) = dz(kb-1)
      dzz(kb)=dzz(kb-1)


      do j = 1,jm
        do i = 1,im
          east_e(i,j) = 100.
     &        + 180.*real(i_global(i)-1)/real(im_global-1)
          north_e(i,j)= -80.
     &        + 160.*real(j_global(j)-1)/real(jm_global-1)
        end do
      end do

      do j = 1,jm
        do i = 2,im
          east_u(i,j) = .5*( east_e(i,j) + east_e(i-1,j) )
          north_u(i,j)= .5*( north_e(i,j)+ north_e(i-1,j))
        end do
      end do
      call exchange2d_mpi( east_u, im,jm)
      call exchange2d_mpi(north_u, im,jm)
      if (n_west==-1) then
        east_u(1,:) = 2.*east_u(2,:) - east_u(3,:)
        north_u(1,:)= 2.*north_u(2,:)- north_u(3,:)
      end if
      do j = 2,jm
        do i = 1,im
          east_v(i,j) = .5*( east_e(i,j) + east_e(i,j-1) )
          north_v(i,j)= .5*( north_e(i,j)+ north_e(i,j-1))
        end do
      end do
      call exchange2d_mpi( east_v, im,jm)
      call exchange2d_mpi(north_v, im,jm)
      if (n_south==-1) then
      east_v(:,1) = 2.*east_v(:,2) - east_v(:,3)
      north_v(:,1)= 2.*north_v(:,2)- north_v(:,3)
      end if

      do j = 1,jm
        do i = 2,im
          dx(i,j) = 111800.*cos(north_e(i,j)*DEG2RAD)
     &                     *abs(east_u(i,j)-east_u(i-1,j))
        end do
      end do
      do j = 2,jm
        do i = 1,im
          dy(i,j) = 111800.*cos(north_e(i,j)*DEG2RAD)
     &                     *abs(north_v(i,j)-north_v(i,j-1))
        end do
      end do
      call exchange2d_mpi(dx, im,jm)
      call exchange2d_mpi(dy, im,jm)
      if (n_west==-1) then
        dx(1,:) = dx(2,:)
        dy(1,:) = dy(2,:)
      end if
      if (n_south==-1) then
        dx(:,1) = dx(:,2)
        dy(:,1) = dy(:,2)
      end if

      rot = 0.
      do j = 1,jm
        do i = 1,im
          h(i,j)=100.-1000.*(tanh(5.*(real(i_global(i))
     &                               -real(im_global)/6.)
     &                              /(real(im_global)/6.
     &                               -real(im_global)))
     &                      +tanh(5.*(real(j_global(j))
     &                               -real(jm_global)/6.)
     &                              /(real(jm_global)/6.
     &                               -real(jm_global))))
          !h(i,j)=h(i,j)*(1.-real(i_global(i))/real(im_global))
          h(i,j)=h(i,j)+1000.*(tanh(5.*(real(im_global
     &                                      -i_global(i))
     &                               +real(im_global)/6.)
     &                              /(real(im_global)/6.
     &                               +real(im_global)))
     &                      +tanh(5.*(real(jm_global-
     &                                     j_global(j))
     &                               +real(jm_global)/6.)
     &                              /(real(jm_global)/6.
     &                               +real(jm_global))))
        end do
      end do

      fsm = 1.
      if (n_west==-1) fsm( 1,:) = 0.
      if (n_east==-1) fsm(im,:) = 0.
      if (n_south==-1) fsm(:, 1) = 0.
      if (n_north==-1) fsm(:,jm) = 0.

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
!     The followings are read in read_grid_pnetcdf:
!     z,zz,dx,dy
!     east_u,east_v,east_e,east_c
!     north_u,north_v,north_e,north_c
!     rot,h,fsm,dum,dvm


      call print_config

! set up Coriolis parameter
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29e-5*sin(north_e(i,j)*DEG2RAD)
          end do
        end do

! inertial period for temporal filter
!      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0

! calculate areas of "t" and "s" cells
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

      if (n_west.eq.-1) then
        aru(1,:)=aru(2,:)
        arv(1,:)=arv(2,:)
      end if

      if (n_south.eq.-1) then
        aru(:,1)=aru(:,2)
        arv(:,1)=arv(:,2)
      end if

      d  = h + el
      dt = h + et


      end

!______________________________________________________________________
! fhx:interp_flag
! fhx:fgrid
      subroutine check_interpolation
!----------------------------------------------------------------------
!  Decides on wether to use interpolation or not.
!______________________________________________________________________

      use config     , only: calc_interp
      use glob_domain, only: im_global, im_global_coarse, is_master
     &                     , jm_global, jm_global_coarse
     &                     , x_division, y_division

      implicit none

      if (  ((im_global_coarse-2)*x_division == (im_global-2))
     $ .and.((jm_global_coarse-2)*y_division == (jm_global-2)) ) then

        calc_interp = .TRUE.

      else if ( (im_global_coarse == im_global)
     $     .and.(jm_global_coarse == jm_global) ) then

        calc_interp = .FALSE.

      else

        if ( is_master ) print '(a//a)'
     &    , 'Incompatible number of *_global and *_global_coarse'
     &    , 'POM terminated with error'
        stop

      end if


      end ! subroutine

!______________________________________________________________________
!
      subroutine check_cfl_min
!----------------------------------------------------------------------
!  Estimates minimum barotropic timestep (sec)
!______________________________________________________________________
!
        use glob_const , only: grav, rk, SMALL
        use glob_domain, only: im, jm, master_task, my_task, POM_COMM
        use grid       , only: dx, dy, fsm, h
        use model_run  , only: dte
        use mpi        , only: MPI_DOUBLE, MPI_MIN, MPI_REAL
     &                       , mpi_reduce

        implicit none

        real(rk), dimension(im,jm) :: cfl
        real(rk)                      cflmin, tmp
        integer                       ierr, MPI_RK
        character(len=128)            msg, desc


        if ( rk == 8 ) then
          MPI_RK = MPI_DOUBLE
        else
          MPI_RK = MPI_REAL
        end if

        cfl = 0.
        cflmin = 1.d10

        cfl = fsm*.5 / sqrt(1./dx**2+1./dy**2) / sqrt(grav*(h+SMALL))

        cflmin = minval(cfl, cfl>0.)

        call mpi_reduce( cflmin, tmp, 1, MPI_RK, mpi_min
     &                 , master_task, POM_COMM, ierr     )

        if ( my_task == master_task ) then
          if ( cflmin < dte ) then
            write(msg, '(a,f5.2,a)')
     &                 "Timestep (", dte, ") is too large!"
            write(desc, '(a,f5.2,a)')
     &              "You are strongly advised to make dte smaller than "
     &                                                       ,cflmin,"."
            call msg_print(msg, 2, desc)
          end if
        end if


      end ! subroutine check_cfl_min
!______________________________________________________________________
!
      subroutine modules_initial_step( dtime )

        use air        , only: air_init => init
        use bry        , only: bry_init => init
        use clim       , only: clm_init => init
        use seaice     , only: ice_init => init
        use tide       , only: tide_init=> init
        use module_time

        implicit none

        type(date), intent(in) :: dtime


        call clm_init( dtime )
        call air_init( dtime )
        call bry_init( dtime )
        call ice_init( dtime )
        call tide_init(dtime )


      end ! subroutine
