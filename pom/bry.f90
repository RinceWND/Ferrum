!______________________________________________________________________
!
! Module `BRY` (bry.f90)
!----------------------------------------------------------------------
!  Module for applying boundary conditions.
!  Adopted from original POM `bcond` subroutine plus I/O routines and
! interpolation in time.
!
!  Author  : RinceWND
!  Created : 2018-09-02
!______________________________________________________________________
!
module bry

  use config     , only: rk ! SO... If a module uses some other module's variable, another module can access the variable through this "medium" module?
  use glob_domain, only: im, imm1, imm2, jm, jmm1, jmm2, kb, kbm1

  implicit none

!----------------------------------------------------------------------
! PRIVATE variables
!----------------------------------------------------------------------
  logical           &
       , private :: &
    DISABLED        & ! flag for active (read from file) boundaries
  , hasEAST         & ! flag for eastern boundary
  , hasNORTH        & ! flag for northern boundary
  , hasSOUTH        & ! flag for southern boundary
  , hasWEST         & ! flag for western boundary
  , INTERP_BRY        ! flag for active boundary interpolation [ NOT IMPLEMENTED ]

  integer           &
       , private :: &
    i, j, k         & ! simple counters
  , read_int          ! interval for reading (days) TODO: though should be months too

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"
       
       
  integer                        &
       , parameter               &
       , private   :: bcNTH = 0  & ! id for northern boundary
                    , bcEST = 1  & ! id for eastern boundary
                    , bcSTH = 2  & ! id for southern boundary
                    , bcWST = 3    ! id for western boundary

!----------------------------------------------------------------------
! PUBLIC variables
!----------------------------------------------------------------------
  public

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  integer*2, parameter ::  & ! Boundary conditions named constants
    bc0GRADIENT       = 0  & !  zero-gradient (free-slip) condition
  , bc3POINTSMOOTH    = 1  & !  smoothing with 3 interior points
  , bcCLAMPED         = 2  & !  forced value condition
  , bcFLATHER         = 3  & !  Flather (tidal) condition
  , bcINOUTFLOW       = 4  & !  in-, outflow condition
  , bcRADIATION       = 5  & !  radiation condition
  , bcORLANSKI        = 6  & !  Orlanski internal normal velocity condition
  , bcRADIATION_ENH   = 7  & !  enhanced radiation condition
  , bcGENFLATHER      = 8    !  generalized Flather (Oddo, Pinardi, 2010)

!----------------------------------------------------------------------
! Configuration
!----------------------------------------------------------------------
  logical USE_SPONGE   ! Flag for sponge zones

  real(kind=rk)      &
    hmax             & ! maximal water depth
  , rfe              & ! flag for eastern open boundary
  , rfn              & ! flag for northern open boundary
  , rfs              & ! flag for southern open boundary
  , rfw                ! flag for western open boundary

!-------------------------------------------------------------------
! Boundary relaxation thickness
!-------------------------------------------------------------------
  integer    &
    NFE      & ! relax pnts next to bound e
  , NFN      & ! relax pnts next to bound n
  , NFS      & ! relax pnts next to bound s
  , NFW        ! relax pnts next to bound w

!----------------------------------------------------------------------
! Custom types definition (some types are private, but not vars)
!----------------------------------------------------------------------
  type BC_PERIODIC    ! Periodic flags. Override following BC values.
    logical x       & !  periodic in x-direction if set to True
          , y         !  periodic in y-direction if set to True
  end type

  type BC_TYPE_VAL     ! Value type for all boundaries.
    integer*2  EAST  & !
            , NORTH  & !  ( self-explanatory )
            , SOUTH  & !
            , WEST     !
  end type

  type BC_TYPE_DIR            ! Differentiate directions for vector BC.
    type(BC_TYPE_VAL) NORM  & !  normal BC
                    , TANG    !  tangential BC
  end type

  type BC_TYPE                   ! BC core parameters:
    type(BC_TYPE_VAL)    ZETA  & !  elevation
                    ,      TS  & !  temperature and salinity (TODO: separate them?)
                    , VELVERT  & !  vertical velocity
                    ,    TURB    !  turbulent parameters (should not be used, though)
    type(BC_TYPE_DIR) VEL2D    & !  external velocity
                    , VEL3D      !  internal velocity
  end type

  type BRY_1D                        ! 1D arrays for boundaries
    real(kind=rk)                  & !
        , allocatable              & !
        , dimension(:)     :: est  & !   east
                            , nth  & !   north
                            , sth  & !   south
                            , wst    !   west
  end type BRY_1D

  type BRY_2D                        ! 2D arrays for boundaries
    real(kind=rk)                  & !
        , allocatable              & !
        , dimension(:,:)   :: est  & ! ( I think you get the idea )
                            , nth  & !
                            , sth  & !
                            , wst    !
  end type BRY_2D

  type BRY_3D                        ! 3D arrays for boundaries
    real(kind=rk)                  & !
        , allocatable              & !
        , dimension(:,:,:) :: est  & !
                            , nth  & !
                            , sth  & !
                            , wst    !
  end type BRY_3D

!----------------------------------------------------------------------
! Custom types variables
!----------------------------------------------------------------------
! configuration vars
  type (BC_PERIODIC) periodic_bc
  type (BC_TYPE)     BC
! data vars
  type (bry_2d) el_bry, ua_bry, va_bry
  type (bry_3d) s_bry, t_bry, u_bry, v_bry

!----------------------------------------------------------------------
! Sponge layer arrays
!----------------------------------------------------------------------
  real(kind=rk)              &
       , allocatable         &
       , dimension(:,:)   :: &
    aamfrz      & ! sponge factor - increased aam at boundaries
  , frz           ! flow-relax-coefficient at boundaries

!----------------------------------------------------------------------
! Path and variable names configuration
!----------------------------------------------------------------------
  character(len=256)  & ! Full paths (with filenames) to:
    bry_path            !  boundary file with t,s,el,u,v for each bry

  character(len=64)   & ! BC file variable names:
    el_name           & !  elevation
  ,  s_name           & !  salinity
  ,  t_name           & !  temperature
  ,  u_name           & !  u-velocity (normal for east and west)
  ,  v_name             !  v-velocity (tangential for east and west)

  contains
!______________________________________________________________________
!
    subroutine allocate_boundary
!----------------------------------------------------------------------
!  Allocates boundary arrays.
!______________________________________________________________________
!
      implicit none


! Allocate arrays
      if ( hasEAST ) then
        allocate(               &
          el_bry%est(NFE,jm)    &
        , ua_bry%est(NFE,jm)    &
        , va_bry%est(NFE,jm)    &
        , S_bry%est(NFE,jm,kb)  &
        , T_bry%est(NFE,jm,kb)  &
        , U_bry%est(NFE,jm,kb)  &
        , V_bry%est(NFE,jm,kb)  &
         )
      end if

      if ( hasNORTH ) then
        allocate(               &
          el_bry%nth(im,NFN)    &
        , ua_bry%nth(im,NFN)    &
        , va_bry%nth(im,NFN)    &
        , S_bry%nth(im,NFN,kb)  &
        , T_bry%nth(im,NFN,kb)  &
        , U_bry%nth(im,NFN,kb)  &
        , V_bry%nth(im,NFN,kb)  &
         )
      end if

      if ( hasSOUTH ) then
        allocate(               &
          el_bry%sth(im,NFS)    &
        , ua_bry%sth(im,NFS)    &
        , va_bry%sth(im,NFS)    &
        , S_bry%sth(im,NFS,kb)  &
        , T_bry%sth(im,NFS,kb)  &
        , U_bry%sth(im,NFS,kb)  &
        , V_bry%sth(im,NFS,kb)  &
         )
      end if

      if ( hasWEST ) then
        allocate(               &
          el_bry%wst(NFW,jm)    &
        , ua_bry%wst(NFW,jm)    &
        , va_bry%wst(NFW,jm)    &
        , S_bry%wst(NFW,jm,kb)  &
        , T_bry%wst(NFW,jm,kb)  &
        , U_bry%wst(NFW,jm,kb)  &
        , V_bry%wst(NFW,jm,kb)  &
         )
      end if

      if ( USE_SPONGE ) then
        allocate(        &
          aamfrz(im,jm)  &
        ,    frz(im,jm)  &
         )
      end if

!      if ( DISABLED ) return
!      
!      allocate( bry_path )

    end ! subroutine allocate_boundary
!______________________________________________________________________
!
    subroutine initialize_boundary( conf )
!----------------------------------------------------------------------
!  Initialize arrays for boundary conditions.
!______________________________________________________________________
!
      use config     , only: use_bry
      use glob_domain, only: is_master, n_east, n_north, n_south, n_west

      implicit none

      character(len=*), intent(in) :: conf

      logical periodic_x, periodic_y
      integer pos
      integer(1)                                                             &
                el_east,         el_north,         el_south,         el_west &
      ,         ts_east,         ts_north,         ts_south,         ts_west &
      ,       turb_east,       turb_north,       turb_south,       turb_west &
      , vel2d_norm_east, vel2d_norm_north, vel2d_norm_south, vel2d_norm_west &
      , vel2d_tang_east, vel2d_tang_north, vel2d_tang_south, vel2d_tang_west &
      , vel3d_norm_east, vel3d_norm_north, vel3d_norm_south, vel3d_norm_west &
      , vel3d_tang_east, vel3d_tang_north, vel3d_tang_south, vel3d_tang_west &
      ,    velvert_east,    velvert_north,    velvert_south,    velvert_west

      namelist/bry_nml/                                                &
        bry_path, el_name, Hmax  , interp_bry, periodic_x, periodic_y  &
      , read_int, s_name , t_name, USE_SPONGE, u_name    , v_name

      namelist/bry_cond/                                                     &
                el_east,         el_north,         el_south,         el_west &
      ,         ts_east,         ts_north,         ts_south,         ts_west &
      ,       turb_east,       turb_north,       turb_south,       turb_west &
      , vel2d_norm_east, vel2d_norm_north, vel2d_norm_south, vel2d_norm_west &
      , vel2d_tang_east, vel2d_tang_north, vel2d_tang_south, vel2d_tang_west &
      , vel3d_norm_east, vel3d_norm_north, vel3d_norm_south, vel3d_norm_west &
      , vel3d_tang_east, vel3d_tang_north, vel3d_tang_south, vel3d_tang_west &
      ,    velvert_east,    velvert_north,    velvert_south,    velvert_west


! Do not read boundary fields by default
      DISABLED = .true.
! Check for active boundary
      if ( use_bry ) DISABLED = .false.

! Set max depth
      hmax = 8000.

! Set default relaxation thickness
      NFE = 1
      NFN = 1
      NFS = 1
      NFW = 1

! Set boundary flags
      hasEAST  = .false.
      hasNORTH = .false.
      hasSOUTH = .false.
      hasWEST  = .false.

      if ( n_east  == -1 ) hasEAST  = .true.
      if ( n_north == -1 ) hasNORTH = .true.
      if ( n_south == -1 ) hasSOUTH = .true.
      if ( n_west  == -1 ) hasWEST  = .true.

! Default configuration flags
      USE_SPONGE = .false.

! Default periodic
      periodic_x = .false.
      periodic_y = .false.

! Default path to active boundary file
      bry_path   = "in/bry/"
! Default variable names
      el_name= "zeta"
      s_name = "salt"
      t_name = "temp"
      u_name = "u"
      v_name = "v"
      
! [ NOT IMPLEMENTED ] Interpolate active boundary by default (ignored if DISABLED)
      interp_bry = .true.

! Default readint interval (daily)
      read_int = 86400

! Default BC types
!   elevation
      el_east  = bc0GRADIENT
      el_north = bc0GRADIENT
      el_south = bc0GRADIENT
      el_west  = bc0GRADIENT
!   external velocity
      vel2d_norm_east  = bcFLATHER ! bcCLAMPED if RFE is 0 in old way
      vel2d_norm_north = bcFLATHER
      vel2d_norm_south = bcFLATHER
      vel2d_norm_west  = bcFLATHER
      vel2d_tang_east  = bcCLAMPED
      vel2d_tang_north = bcCLAMPED
      vel2d_tang_south = bcCLAMPED
      vel2d_tang_west  = bcCLAMPED
!   internal velocity
      vel3d_norm_east  = bcRADIATION
      vel3d_norm_north = bcRADIATION
      vel3d_norm_south = bcRADIATION
      vel3d_norm_west  = bcRADIATION
      vel3d_tang_east  = bc3POINTSMOOTH
      vel3d_tang_north = bc3POINTSMOOTH
      vel3d_tang_south = bc3POINTSMOOTH
      vel3d_tang_west  = bc3POINTSMOOTH
!   temp and salt
      ts_east  = bcINOUTFLOW
      ts_north = bcINOUTFLOW
      ts_south = bcINOUTFLOW
      ts_west  = bcINOUTFLOW
!   vertical velocity (IGNORED anyway)
      velvert_east  = bcCLAMPED
      velvert_north = bcCLAMPED
      velvert_south = bcCLAMPED
      velvert_west  = bcCLAMPED
!   turbulent parameters
      turb_east  = bcINOUTFLOW
      turb_north = bcINOUTFLOW
      turb_south = bcINOUTFLOW
      turb_west  = bcINOUTFLOW

! Override configuration
      open ( 73, file = conf, status = 'old' )
      read ( 73, nml = bry_nml )
      read ( 73, nml = bry_cond )
      close( 73 )

! Manage derived type variables
      periodic_bc % x = periodic_x
      periodic_bc % y = periodic_y

!   elevation
      BC % zeta % EAST  = el_east
      BC % zeta % NORTH = el_north
      BC % zeta % SOUTH = el_south
      BC % zeta % WEST  = el_west
!   external velocity
      BC % vel2d % NORM % EAST  = vel2d_norm_east
      BC % vel2d % NORM % NORTH = vel2d_norm_north
      BC % vel2d % NORM % SOUTH = vel2d_norm_south
      BC % vel2d % NORM % WEST  = vel2d_norm_west
      BC % vel2d % TANG % EAST  = vel2d_tang_east
      BC % vel2d % TANG % NORTH = vel2d_tang_north
      BC % vel2d % TANG % SOUTH = vel2d_tang_south
      BC % vel2d % TANG % WEST  = vel2d_tang_west
!   internal velocity
      BC % vel3d % NORM % EAST  = vel3d_norm_east
      BC % vel3d % NORM % NORTH = vel3d_norm_north
      BC % vel3d % NORM % SOUTH = vel3d_norm_south
      BC % vel3d % NORM % WEST  = vel3d_norm_west
      BC % vel3d % TANG % EAST  = vel3d_tang_east
      BC % vel3d % TANG % NORTH = vel3d_tang_north
      BC % vel3d % TANG % SOUTH = vel3d_tang_south
      BC % vel3d % TANG % WEST  = vel3d_tang_east
!   temp and salt
      BC % ts % EAST  = ts_east
      BC % ts % NORTH = ts_north
      BC % ts % SOUTH = ts_south
      BC % ts % WEST  = ts_west
!   vertical velocity (IGNORED anyway)
      BC % velvert % EAST  = velvert_east
      BC % velvert % NORTH = velvert_north
      BC % velvert % SOUTH = velvert_south
      BC % velvert % WEST  = velvert_west
!   turbulent parameters
      BC % turb % EAST  = turb_east
      BC % turb % NORTH = turb_north
      BC % turb % SOUTH = turb_south
      BC % turb % WEST  = turb_west
      
      pos = len(trim(bry_path))
      if ( bry_path(pos:pos) == "/" ) then
        bry_path = trim(bry_path)//"bry."
      end if

! Print summary [ TODO: print all the BC settings ]
      if ( is_master ) then
        if ( DISABLED ) then
          print *, "Active boundary: [ DISABLED ]"
        else
          print *, "Active boundary: "
          print *, " Bry file path: ", trim(bry_path)
        end if
      end if

! Allocate arrays
      call allocate_boundary

      call msg_print("BRY MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_boundary
!______________________________________________________________________
!
    subroutine initial_conditions_boundary
!----------------------------------------------------------------------
!  Sets up initial conditions for boundary arrays.
!______________________________________________________________________
!
      use glob_const , only: SEC2DAY
      use glob_domain, only: im, jm, n_east, n_north, n_south, n_west
      use glob_grid  , only: fsm
      use glob_ocean , only: elb, sclim, tclim
      use model_run  , only: dti

      implicit none

      integer  ii, jj
      real(rk) rdisp


! Initial conditions ! TODO: Move to a separate subroutine?
      if ( hasEAST ) then

        EL_bry % EST(1,:) = elb(im,:)

        UA_bry % EST = 0.
        VA_bry % EST = 0.

        U_bry % EST = 0.
        V_bry % EST = 0.

        do k=1,kb
          do j=1,jm
            do i=1,NFE
              ii=im-i+1
              T_bry%EST(i,j,k) = tclim(ii,j,k) * fsm(ii,j)
              S_bry%EST(i,j,k) = sclim(ii,j,k) * fsm(ii,j)
            end do
          end do
        end do

      end if


      if ( hasWEST ) then

        EL_bry % WST(1,:) = elb(1,:)

        UA_bry % WST = 0.
        VA_bry % WST = 0.

        U_bry % WST = 0.
        V_bry % WST = 0.

        do k=1,kb
          do j=1,jm
            do i=1,NFW
              T_bry%WST(i,j,k) = tclim(i,j,k) * fsm(i,j)
              S_bry%WST(i,j,k) = sclim(i,j,k) * fsm(i,j)
            end do
          end do
        end do

      end if


      if ( hasNORTH ) then

        EL_bry % NTH(:,1) = elb(:,jm)

        UA_bry % NTH = 0.
        VA_bry % NTH = 0.

        U_bry % NTH = 0.
        V_bry % NTH = 0.

        do k=1,kb
          do i=1,im
            do j=1,NFN
              jj=jm-j+1
              T_bry%NTH(i,j,k) = tclim(i,jj,k) * fsm(i,jj)
              S_bry%NTH(i,j,k) = sclim(i,jj,k) * fsm(i,jj)
            end do
          end do
        end do

      end if


      if ( hasSOUTH ) then

        EL_bry % STH(:,1) = elb(:,1)

        UA_bry % STH = 0.
        VA_bry % STH = 0.

        U_bry % STH = 0.
        V_bry % STH = 0.

        do k=1,kb
          do i=1,im
            do j=1,NFS
              T_bry%STH(i,j,k) = tclim(i,j,k) * fsm(i,j)
              S_bry%STH(i,j,k) = sclim(i,j,k) * fsm(i,j)
            end do
          end do
        end do

      end if


      if ( USE_SPONGE ) then
!lyo:pac10:beg:
! set lateral boundary FRZ & SPONGE layer parameters !lyo:20110224:alu:stcc:
!      call bfrz(   nfw,   nfe,    nfs,    nfn,
!     $          n_west,n_east,n_south,n_north,
!     $              im,    jm,      6,    frz)       !FRZ:
!     For FRZ in bcond: T[n+1]={ (1-frz)*T[n] + frz*Tobs },
!     the inverse-relax-time = frz/dti,
!     so that since "bfrz" gives frz=1 @boundaries,
!     frz=1 (=dti/dti) -> inverse-relax-time = 1/dti @boundaries, and
!     frz=dti/86400.   -> inverse-relax-time = 1/86400 @boundaries etc.
        rdisp = dti*SEC2DAY
        frz = frz*rdisp !lyo:stcc:mar_:dec_:
!     rdisp=1.;         frz(:,:)=frz(:,:)*rdisp !lyo:stcc:
!     For aam is increased by (1.+aamfrz(i,j)) in subr lateral_viscosity
        call bfrz(     0,     0,      0,      0   &
                 ,n_west,n_east,n_south,n_north   &
                 ,    im,    jm,      6, aamfrz)       !SPONGE:
        rdisp  = 0.
        aamfrz = aamfrz*rdisp
!lyo:pac10:end:
      end if

!      if ( hasEAST ) then
!        do j=2,jmm1
!          UA_bry%EST(1,j) = uab(imm1,j)
!          EL_bry%EST(1,j) = EL_bry%EST(1,j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
!        end do
!        EL_bry%EST(1,:) = (EL_bry%EST(1,:)-EL_bry%EST(1,jmm1/2))*fsm(im,:)
!      end if
!      if ( hasWEST ) then
!        dum(1,:) = 0.
!        do j=2,jmm1
!          UA_bry%WST(1,j) = uab(2,j)
!! set geostrophically conditioned elevations at the boundaries
!          EL_bry%WST(1,j) = EL_bry%WST(1,j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
!        end do
!        EL_bry%WST(1,:) = (EL_bry%WST(1,:)-EL_bry%WST(1,jmm1/2))*fsm(2,:)
!      end if


    end ! subroutine initial_conditions_boundary
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads initial lateral boundary conditions.
!______________________________________________________________________
!
      use module_time

      type(date), intent(in) :: d_in

      integer            max_in_prev, max_in_this
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      chunk     = chunk_of_year( d_in, read_int )
      record(1) = int(chunk)

      if ( chunk - record(1) < .5 ) then
        record(2) = record(1)
      else
        record(2) = record(1) + 1
      end if
      record(3) = record(2) + 1

      if ( record(2) == 0 ) then
        record(2) = max_in_prev + 1 ! TODO: [ NEEDS TESTING ]
        year(2) = d_in%year - 1
      elseif ( record(3) == max_in_this + 1 ) then
        record(3) = 1
        year(3) = d_in%year + 1
      end if

      record(1) = record(1) + 1

! Read wind if momentum flux is derived from wind
      call read_all( .true., 1, year, record )

      call msg_print("BRY INITIALIZED", 2, "")


    end ! subroutine init
!______________________________________________________________________
!
    subroutine step( d_in )
!----------------------------------------------------------------------
!  Reads lateral boundary conditions.
!______________________________________________________________________
!
      use glob_grid  , only: dz
!      use glob_domain, only: is_master
      use model_run  , only: dti, secs => sec_of_year
      use module_time

      implicit none

      type(date), intent(in) :: d_in

      logical            ADVANCE_REC
      integer            max_in_prev, max_in_this
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      chunk     = chunk_of_year( d_in, read_int )
      record(1) = int(chunk)

      ADVANCE_REC = .false.
      if ( secs - record(1)*read_int < dti ) then ! TODO: test this one as well.
        ADVANCE_REC = .true.
      end if

!      if ( is_master ) print *, "III", secs-record(1)*read_int, dti, ADVANCE_REC

      record(1) = record(1) + 1

      call read_all( ADVANCE_REC, 1, year, record )

      if ( ADVANCE_REC ) then
! Derive barotropic velocities
        if ( hasEAST ) then
          UA_bry % EST = 0.
          VA_bry % EST = 0.
          do k = 1, kb
            UA_bry % EST = UA_bry % EST                &
                         +  U_bry % EST(:,:,k) * dz(k)
            VA_bry % EST = VA_bry % EST                &
                         +  V_bry % EST(:,:,k) * dz(k)
          end do
        end if
        if ( hasNORTH ) then
          UA_bry % NTH = 0.
          VA_bry % NTH = 0.
          do k = 1, kb
            UA_bry % NTH = UA_bry % NTH                &
                         +  U_bry % NTH(:,:,k) * dz(k)
            VA_bry % NTH = VA_bry % NTH                &
                         +  V_bry % NTH(:,:,k) * dz(k)
          end do
        end if
        if ( hasSOUTH ) then
          UA_bry % STH = 0.
          VA_bry % STH = 0.
          do k = 1, kb
            UA_bry % STH = UA_bry % STH                &
                         +  U_bry % STH(:,:,k) * dz(k)
            VA_bry % STH = VA_bry % STH                &
                         +  V_bry % STH(:,:,k) * dz(k)
          end do
        end if
        if ( hasWEST ) then
          UA_bry % WST = 0.
          VA_bry % WST = 0.
          do k = 1, kb
            UA_bry % WST = UA_bry % WST                &
                         +  U_bry % WST(:,:,k) * dz(k)
            VA_bry % WST = VA_bry % WST                &
                         +  V_bry % WST(:,:,k) * dz(k)
          end do
        end if
      end if

    end ! subroutine step
!______________________________________________________________________
!
    subroutine read_all( execute, n, year, record )

      use glob_domain, only: i_global, j_global
      use glob_grid  , only: fsm
      use glob_ocean , only: sclim, tclim
      use mpi        , only: MPI_OFFSET_KIND

      implicit none

      logical              , intent(in) :: execute
      integer              , intent(in) :: n
      integer, dimension(3), intent(in) :: record, year

      integer                  ii, ncid
      integer(MPI_OFFSET_KIND) start(4), edge(4)
      real(rk)                 dummy(1,1,1)
      character(len=128)       desc


      if ( .not. execute ) return

      if ( n >= 2 ) then
        write(desc,'("[interp_bry not implemented yet] #",i4," @ ",i4)') &
            record(n), year(n)
      else
        write(desc,'("Reading BC record #",i4," @ ",i4)') &
            record(1), year(1)
      end if

      call msg_print("", 1, desc)

      ncid = -1

! EAST
      if ( hasEAST ) then
! set reading bounds for 3D vars
        if ( NFE > 1 ) then
          start(1) = i_global(1) + im-NFE
          start(2) = j_global(1)
          start(3) = 1
          start(4) = record(n)
          edge(1) = NFE
          edge(2) = jm
          edge(3) = kb
          edge(4) = 1
        else
          start(1) = j_global(1)
          start(2) = 1
          start(3) = record(n)
          start(4) = 1
          edge(1)   = jm
          edge(2)   = kb
          edge(3:4) = 1
        end if
! Temperature
        call read_var_3d_nc( bry_path, "east_"//t_name, T_bry%EST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,jm
              do i=1,NFE
                ii=im-i+1
                T_bry%EST(i,j,k) = tclim(ii,j,k) * fsm(ii,j)
              end do
            end do
          end do
          call msg_print("", 2, "Temp@EAST defaulted to climatology")
        end if
! Salinity
        call read_var_3d_nc( bry_path, "east_"//s_name, S_bry%EST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,jm
              do i=1,NFE
                ii=im-i+1
                S_bry%EST(i,j,k) = sclim(ii,j,k) * fsm(ii,j)
              end do
            end do
          end do
          call msg_print("", 2, "Salt@EAST defaulted to climatology")
        end if
! Normal velocity
        call read_var_3d_nc( bry_path, "east_"//u_name, U_bry%EST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          U_bry%EST = 0.
          call msg_print("", 2, "Uvel@EAST defaulted to zero")
        end if
! Tangential velocity
        call read_var_3d_nc( bry_path, "east_"//v_name, V_bry%EST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          V_bry%EST = 0.
          call msg_print("", 2, "Vvel@EAST defaulted to zero")
        end if
! set reading bounds for 2D vars
        if ( NFE > 1 ) then
          start(1) = i_global(1) + im-NFE
          start(2) = j_global(1)
          start(3) = record(n)
          start(4) = 1
          edge(1)   = NFE
          edge(2)   = jm
          edge(3:4) = 1
        else
          start(1)   = j_global(1)
          start(2)   = record(n)
          start(3:4) = 1
          edge(1)   = jm
          edge(2:4) = 1
        end if
! Elevation
        call read_var_2d_nc( bry_path, "east_"//el_name, EL_bry%EST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          EL_bry%EST = 0.
          call msg_print("", 2, "Elev@EAST defaulted to zero")
        end if

      else

        start = 1
        edge = 0
! Temperature
        call read_var_3d_nc( bry_path, "east_"//t_name, dummy  &
                           , year(n), start, edge, ncid )
! Salinity
        call read_var_3d_nc( bry_path, "east_"//s_name, dummy  &
                           , year(n), start, edge, ncid )
! Normal velocity
        call read_var_3d_nc( bry_path, "east_"//u_name, dummy  &
                           , year(n), start, edge, ncid )
! Tangential velocity
        call read_var_3d_nc( bry_path, "east_"//v_name, dummy  &
                           , year(n), start, edge, ncid )
! Elevation
        call read_var_2d_nc( bry_path, "east_"//el_name, dummy(1,:,:)  &
                           , year(n), start, edge, ncid )

      end if

! NORTH
      if ( hasNORTH ) then
! set reading bounds for 3D vars
        if ( NFN > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1) + jm-NFN
          start(3) = 1
          start(4) = record(n)
          edge(1) = im
          edge(2) = NFN
          edge(3) = kb
          edge(4) = 1
        else
          start(1) = i_global(1)
          start(2) = 1
          start(3) = record(n)
          start(4) = 1
          edge(1)   = im
          edge(2)   = kb
          edge(3:4) = 1
        end if
! Temperature
        call read_var_3d_nc( bry_path, "north_"//t_name, T_bry%NTH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,NFN
              do i=1,im
                ii=jm-j+1
                T_bry%NTH(i,j,k) = tclim(i,ii,k) * fsm(i,ii)
              end do
            end do
          end do
          call msg_print("", 2, "Temp@NORTH defaulted to climatology")
        end if
! Salinity
        call read_var_3d_nc( bry_path, "north_"//s_name, S_bry%NTH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,NFN
              do i=1,im
                ii=jm-j+1
                S_bry%NTH(i,j,k) = sclim(i,ii,k) * fsm(i,ii)
              end do
            end do
          end do
          call msg_print("", 2, "Salt@NORTH defaulted to climatology")
        end if
! Normal velocity
        call read_var_3d_nc( bry_path, "north_"//u_name, U_bry%NTH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          U_bry%NTH = 0.
          call msg_print("", 2, "Uvel@NORTH defaulted to zero")
        end if
! Tangential velocity
        call read_var_3d_nc( bry_path, "north_"//v_name, V_bry%NTH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          V_bry%NTH = 0.
          call msg_print("", 2, "Vvel@NORTH defaulted to zero")
        end if
! set reading bounds for 2D vars
        if ( NFN > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1) + jm-NFN
          start(3) = record(n)
          start(4) = 1
          edge(1)   = im
          edge(2)   = NFN
          edge(3:4) = 1
        else
          start(1)   = i_global(1)
          start(2)   = record(n)
          start(3:4) = 1
          edge(1)   = im
          edge(2:4) = 1
        end if
! Elevation
        call read_var_2d_nc( bry_path, "north_"//el_name, EL_bry%NTH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          EL_bry%NTH = 0.
          call msg_print("", 2, "Elev@NORTH defaulted to zero")
        end if

      else

        start = 1
        edge = 0
! Temperature
        call read_var_3d_nc( bry_path, "north_"//t_name, dummy  &
                           , year(n), start, edge, ncid )
! Salinity
        call read_var_3d_nc( bry_path, "north_"//s_name, dummy  &
                           , year(n), start, edge, ncid )
! Normal velocity
        call read_var_3d_nc( bry_path, "north_"//u_name, dummy  &
                           , year(n), start, edge, ncid )
! Tangential velocity
        call read_var_3d_nc( bry_path, "north_"//v_name, dummy  &
                           , year(n), start, edge, ncid )
! Elevation
        call read_var_2d_nc( bry_path, "north_"//el_name, dummy(1,:,:)  &
                           , year(n), start, edge, ncid )

      end if

! SOUTH
      if ( hasSOUTH ) then
! set reading bounds for 3D vars
        if ( NFS > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = 1
          start(4) = record(n)
          edge(1) = im
          edge(2) = NFS
          edge(3) = kb
          edge(4) = 1
        else
          start(1) = i_global(1)
          start(2) = 1
          start(3) = record(n)
          start(4) = 1
          edge(1)   = im
          edge(2)   = kb
          edge(3:4) = 1
        end if
! Temperature
        call read_var_3d_nc( bry_path, "south_"//t_name, T_bry%STH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,NFS
              do i=1,im
                T_bry%STH(i,j,k) = tclim(i,j,k) * fsm(i,j)
              end do
            end do
          end do
          call msg_print("", 2, "Temp@SOUTH defaulted to climatology")
        end if
! Salinity
        call read_var_3d_nc( bry_path, "south_"//s_name, S_bry%STH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,NFS
              do i=1,im
                S_bry%STH(i,j,k) = sclim(i,j,k) * fsm(i,j)
              end do
            end do
          end do
          call msg_print("", 2, "Salt@SOUTH defaulted to climatology")
        end if
! Normal velocity
        call read_var_3d_nc( bry_path, "south_"//u_name, U_bry%STH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          U_bry%STH = 0.
          call msg_print("", 2, "Uvel@SOUTH defaulted to zero")
        end if
! Tangential velocity
        call read_var_3d_nc( bry_path, "south_"//v_name, V_bry%STH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          V_bry%STH = 0.
          call msg_print("", 2, "Vvel@SOUTH defaulted to zero")
        end if
! set reading bounds for 2D vars
        if ( NFS > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = record(n)
          start(4) = 1
          edge(1)   = im
          edge(2)   = NFS
          edge(3:4) = 1
        else
          start(1)   = i_global(1)
          start(2)   = record(n)
          start(3:4) = 1
          edge(1)   = im
          edge(2:4) = 1
        end if
! Elevation
        call read_var_2d_nc( bry_path, "south_"//el_name, EL_bry%STH  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          EL_bry%STH = 0.
          call msg_print("", 2, "Elev@SOUTH defaulted to zero")
        end if

      else

        start = 1
        edge = 0
! Temperature
        call read_var_3d_nc( bry_path, "south_"//t_name, dummy  &
                           , year(n), start, edge, ncid )
! Salinity
        call read_var_3d_nc( bry_path, "south_"//s_name, dummy  &
                           , year(n), start, edge, ncid )
! Normal velocity
        call read_var_3d_nc( bry_path, "south_"//u_name, dummy  &
                           , year(n), start, edge, ncid )
! Tangential velocity
        call read_var_3d_nc( bry_path, "south_"//v_name, dummy  &
                           , year(n), start, edge, ncid )
! Elevation
        call read_var_2d_nc( bry_path, "south_"//el_name, dummy(1,:,:)  &
                           , year(n), start, edge, ncid )

      end if

! WEST
      if ( hasWEST ) then
! set reading bounds for 3D vars
        if ( NFW > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = 1
          start(4) = record(n)
          edge(1) = NFW
          edge(2) = jm
          edge(3) = kb
          edge(4) = 1
        else
          start(1) = j_global(1)
          start(2) = 1
          start(3) = record(n)
          start(4) = 1
          edge(1)   = jm
          edge(2)   = kb
          edge(3:4) = 1
        end if
! Temperature
        call read_var_3d_nc( bry_path, "west_"//t_name, T_bry%WST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,jm
              do i=1,NFW
                T_bry%WST(i,j,k) = tclim(i,j,k) * fsm(i,j)
              end do
            end do
          end do
          call msg_print("", 2, "Temp@WEST defaulted to climatology")
        end if
! Salinity
        call read_var_3d_nc( bry_path, "west_"//s_name, S_bry%WST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          do k=1,kb
            do j=1,jm
              do i=1,NFW
                S_bry%WST(i,j,k) = sclim(i,j,k) * fsm(i,j)
              end do
            end do
          end do
          call msg_print("", 2, "Salt@WEST defaulted to climatology")
        end if
! Normal velocity
        call read_var_3d_nc( bry_path, "west_"//u_name, U_bry%WST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          U_bry%WST = 0.
          call msg_print("", 2, "Uvel@WEST defaulted to zero")
        end if
! Tangential velocity
        call read_var_3d_nc( bry_path, "west_"//v_name, V_bry%WST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          V_bry%WST = 0.
          call msg_print("", 2, "Vvel@WEST defaulted to zero")
        end if
! set reading bounds for 2D vars
        if ( NFW > 1 ) then
          start(1) = i_global(1)
          start(2) = j_global(1)
          start(3) = record(n)
          start(4) = 1
          edge(1)   = NFW
          edge(2)   = jm
          edge(3:4) = 1
        else
          start(1)   = j_global(1)
          start(2)   = record(n)
          start(3:4) = 1
          edge(1)   = jm
          edge(2:4) = 1
        end if
! Elevation
        call read_var_2d_nc( bry_path, "west_"//el_name, EL_bry%WST  &
                           , year(n), start, edge, ncid )
        if ( ncid == -1 ) then
          EL_bry%WST = 0.
          call msg_print("", 2, "Elev@WEST defaulted to zero")
        end if

      else

        start = 1
        edge = 0
! Temperature
        call read_var_3d_nc( bry_path, "west_"//t_name, dummy  &
                           , year(n), start, edge, ncid )
! Salinity
        call read_var_3d_nc( bry_path, "west_"//s_name, dummy  &
                           , year(n), start, edge, ncid )
! Normal velocity
        call read_var_3d_nc( bry_path, "west_"//u_name, dummy  &
                           , year(n), start, edge, ncid )
! Tangential velocity
        call read_var_3d_nc( bry_path, "west_"//v_name, dummy  &
                           , year(n), start, edge, ncid )
! Elevation
        call read_var_2d_nc( bry_path, "west_"//el_name, dummy(1,:,:)  &
                           , year(n), start, edge, ncid )

      end if


    end ! subroutine read_all
!______________________________________________________________________
!
    subroutine bc_zeta
!----------------------------------------------------------------------
!  Apply external (2D) elevation boundary conditions.
!______________________________________________________________________
!
      use glob_const, only: g => GRAV
      use glob_grid , only: dx, dy, fsm, h
      use glob_ocean, only: elb, elf, ua, uab, uaf, va, vab, vaf
      use model_run , only: dti2

      implicit none


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi2d_mpi(elf,im,jm)
!...or...
      else
! EAST
        if ( hasEAST ) then

          select case ( BC % zeta % east )

            case ( bcFLATHER ) ! Complementary BCs for 2d vel should be bcFLATHER for tangential velocities.
                               ! TODO: Normal, however, should be deduced from continuity equation...
                               ! TODO: CHECK if these boundaries have any effect on interior at all...
                               ! [ DO NOT USE ]
              elf(imm1,2:jmm1) = ( elb(imm1,2:jmm1)                 &
                                  *(uab(imm1,2:jmm1)-sqrt(g*h(imm1,2:jmm1)))     &
                                  +uab(imm1,2:jmm1)*h(imm1,2:jmm1)  &
                                  -uab(imm1,2:jmm1)*h(im,2:jmm1) ) &
                                 /(uab(imm1,2:jmm1)-sqrt(g*h(imm1,2:jmm1)))


            elf(im,:) = elf(imm1,:)

          end select

        end if
! WEST
        if ( hasWEST ) then

          select case ( BC % zeta % west )

            case ( bcFLATHER )
              elf( 2,2:jmm1) = EL_bry % WST(1,2:jmm1)       &
                              - ( uaf(3,2:jmm1)             &
                                 -UA_bry % WST(1,2:jmm1) )  &
                                /sqrt(g*h(2,2:jmm1))        &
                              + elb(2,2:jmm1)               &
                              + dti2*(uab(2,2:jmm1)-uab(3,2:jmm1)) &
                               /dx(2,2:jmm1)

            elf( 1,:) = elf(2   ,:)

          end select

        end if

      end if


! Apply periodic BC in y-dimension
      if ( periodic_bc % y ) then

        call yperi2d_mpi(elf,im,jm)

      else
! NORTH
        if ( hasNORTH ) then

          select case ( BC % zeta % north )

            case ( bcFLATHER )
              elf(2:imm1,jmm1) = EL_bry % NTH(2:imm1,1)     &
                              + ( vaf(2:imm1,jmm1)          &
                                 -VA_bry % NTH(2:imm1,1) )  &
                                /sqrt(g*h(2:imm1,jmm1))     &
                              + elb(2:imm1,jmm1)            &
                              + dti2*(vab(2:imm1,jmm1)-vab(2:imm1,jmm2)) &
                               /dy(2:imm1,jmm1)

            elf(:,jm) = elf(:,jmm1)

          end select

        end if
! SOUTH
        if ( hasSOUTH ) then

          select case ( BC % zeta % south )

            case ( bcFLATHER )
              elf(2:imm1, 2) = EL_bry % STH(2:imm1,1)       &
                              - ( vaf(2:imm1,3)             &
                                 -VA_bry % STH(2:imm1,1) )  &
                                /sqrt(g*h(2:imm1,2))        &
                              + elb(2:imm1,2)               &
                              + dti2*(vab(2:imm1,2)-vab(2:imm1,3)) &
                               /dy(2:imm1,2)

            elf(:, 1) = elf(:,   2)

          end select

        end if

      end if

! Apply rho-mask (TODO: needed?)
      elf = elf*fsm


    end ! subroutine bc_zeta
!______________________________________________________________________
!
    subroutine bc_vel_ext
!----------------------------------------------------------------------
!  Apply external (2D) velocity boundary conditions.
!______________________________________________________________________
!
      use glob_const , only: GRAV
      use glob_grid  , only: dum, dvm, h, dx, dy
      use glob_ocean , only: el, uaf, vaf
      use model_run  , only: ramp

      implicit none


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi2d_mpi(uaf,im,jm)
        call xperi2d_mpi(vaf,im,jm)
        if ( hasNORTH ) then
          if ( BC % VEL2D % TANG % NORTH == bc0GRADIENT ) then
            uaf(:,jm) = uaf(:,jmm1)
            dum(:,jm) = 1.
          end if
        end if
        if ( hasSOUTH ) then
          if ( BC % VEL2D % TANG % SOUTH == bc0GRADIENT ) then
            uaf(:,1) = uaf(:,2)
            dum(:,1) = 1.
          end if
        end if

      else
! EAST
        if ( hasEAST ) then

          select case ( BC % VEL2D % NORM % EAST )

            case ( bc0GRADIENT )
              uaf(im,2:jmm1) = uaf(imm1,2:jmm1)

            case ( bc3POINTSMOOTH )
              uaf(im,2:jmm1) = ( uaf(imm1,1:jmm2)       &
                               + uaf(imm1,2:jmm1)       &
                               + uaf(imm1,3:jm  ) )/3.

            case ( bcCLAMPED )
              uaf(im,2:jmm1) = UA_bry%EST(1,2:jmm1)

            case ( bcFLATHER )
              uaf(im,2:jmm1) = UA_bry%EST(1,2:jmm1)              &
                         + sqrt(GRAV/h(imm1,2:jmm1))             &
                          *(el(imm1,2:jmm1)-EL_bry%EST(1,2:jmm1))

            case ( bcGENFLATHER )
              uaf(im,2:jmm1) = ( UA_bry%EST(1,2:jmm1)                   &
                                *(h(imm1,2:jmm1)+el(imm1,2:jmm1))       &
                                +sqrt(GRAV/h(imm1,2:jmm1))              &
                                *(el(imm1,2:jmm1)-EL_bry%EST(1,2:jmm1)))&
                               /(h(im,2:jmm1)+EL_bry%EST(1,2:jmm1))

            case default
              uaf(im,2:jmm1) = 0.

          end select

          uaf(im,2:jmm1) = ramp*uaf(im,2:jmm1)

          select case ( BC % VEL2D % TANG % EAST )

            case ( bc0GRADIENT )
              vaf(im,2:jmm1) = vaf(imm1,2:jmm1)

            case ( bc3POINTSMOOTH )
              vaf(im,2:jmm1) = ( vaf(imm1,1:jmm2)       &
                               + vaf(imm1,2:jmm1)       &
                               + vaf(imm1,3:jm  ) )/3.

            case ( bcCLAMPED )
              vaf(im,2:jmm1) = VA_bry%EST(1,2:jmm1)

            case ( bcFLATHER ) ! ??? Zero out tangential velocity gradient at imm1 point
                               ! (So we assume {partial U^T} over {partial n} = {partial U^T_exterior) over {partial n} = 0)
!              vaf(im,2:jmm1) = vaf(imm2,2:jmm1)
              vaf(im,2:jmm1) = vaf(imm1,2:jmm1)                     &
                             + (vaf(imm1,2:jmm1)-vaf(imm2,2:jmm1))  !&
!                              *dx(imm1,2:jmm1)/dx(imm2,2:jmm1)

            case default
              vaf(im,2:jmm1) = 0.

          end select

        end if

! WEST
        if ( hasWEST ) then

          select case ( BC % VEL2D % NORM % WEST )

            case ( bc0GRADIENT )
              uaf(2,2:jmm1) = uaf(3,2:jmm1)

            case ( bc3POINTSMOOTH )
              uaf(2,2:jmm1) = ( uaf(3,1:jmm2)       &
                              + uaf(3,2:jmm1)       &
                              + uaf(3,3:jm  ) )/3.

            case ( bcCLAMPED )
              uaf(2,2:jmm1) = UA_bry%WST(1,2:jmm1)

            case ( bcFLATHER )
              uaf(2,2:jmm1) = UA_bry%WST(1,2:jmm1)              &
                           - sqrt(GRAV/h(2,2:jmm1))             &
                            *(el(2,2:jmm1)-EL_bry%WST(1,2:jmm1))

            case ( bcGENFLATHER )
              uaf(2,2:jmm1) = ( UA_bry%WST(1,2:jmm1)                 &
                               *(h(2,2:jmm1)+el(2,2:jmm1))           &
                               -sqrt(GRAV/h(2,2:jmm1))               &
                               *(el(2,2:jmm1)-EL_bry%WST(1,2:jmm1))) &
                              /(h(1,2:jmm1)+EL_bry%WST(1,2:jmm1))

            case default
              uaf(2,2:jmm1) = 0.

          end select

          uaf(2,2:jmm1) = ramp*uaf(2,2:jmm1)
          uaf(1,2:jmm1) =      uaf(2,2:jmm1)

          select case ( BC % VEL2D % TANG % WEST )

            case ( bc0GRADIENT )
              vaf(1,2:jmm1) = vaf(2,2:jmm1)

            case ( bc3POINTSMOOTH )
              vaf(1,2:jmm1) = ( vaf(2,1:jmm2)       &
                              + vaf(2,2:jmm1)       &
                              + vaf(2,3:jm  ) )/3.

            case ( bcCLAMPED )
              vaf(1,2:jmm1) = VA_bry%WST(1,2:jmm1)

            case ( bcFLATHER ) ! ??? (See east boundary)
!              vaf(1,2:jmm1) = vaf(3,2:jmm1)
              vaf(1,2:jmm1) = vaf(2,2:jmm1)                   &
                             + (vaf(2,2:jmm1)-vaf(3,2:jmm1))  !&
!                              *dx(2,2:jmm1)/dx(3,2:jmm1)

            case default
              vaf(1,2:jmm1) = 0.

          end select

        end if

      end if ! periodic_x end

! Apply periodic BC in y-dimension
      if ( periodic_bc % y ) then

        call yperi2d_mpi(uaf,im,jm)
        call yperi2d_mpi(vaf,im,jm)
        if ( hasEAST ) then
          if ( BC % VEL2D % TANG % EAST == bc0GRADIENT ) then
            vaf(im,:) = vaf(imm1,:)
            dvm(im,:) = 1.0
          end if
        end if
        if ( hasWEST ) then
          if ( BC % VEL2D % TANG % WEST == bc0GRADIENT ) then
            vaf(1,:) = vaf(2,:)
            dvm(1,:) = 1.0
          end if
        end if

      else

! NORTH
        if ( hasNORTH ) then

          select case ( BC % VEL2D % NORM % NORTH )

            case ( bc0GRADIENT )
              vaf(2:imm1,jm) = vaf(2:imm1,jmm1)

            case ( bc3POINTSMOOTH )
              vaf(2:imm1,jm) = ( vaf(1:imm2,jmm1)       &
                               + vaf(2:imm1,jmm1)       &
                               + vaf(3:im  ,jmm1) )/3.

            case ( bcCLAMPED )
              vaf(2:imm1,jm) = VA_bry%NTH(2:imm1,1)

            case ( bcFLATHER )
              vaf(2:imm1,jm) = VA_bry%NTH(2:imm1,1)              &
                         + sqrt(GRAV/h(2:imm1,jmm1))             &
                          *(el(2:imm1,jmm1)-EL_bry%NTH(2:imm1,1))

            case ( bcGENFLATHER )
              vaf(2:imm1,jm) = ( VA_bry%NTH(2:imm1,1)                   &
                                *(h(2:imm1,jmm1)+el(2:imm1,jmm1))       &
                                +sqrt(GRAV/h(2:imm1,jmm1))              &
                                *(el(2:imm1,jmm1)-EL_bry%NTH(2:imm1,1)))&
                               /(h(2:imm1,jm)+EL_bry%NTH(2:imm1,1))

            case default
              vaf(2:imm1,jm) = 0.

          end select

          vaf(2:imm1,jm) = ramp*vaf(2:imm1,jm)

          select case ( BC % VEL2D % TANG % NORTH )

            case ( bc0GRADIENT )
              uaf(2:imm1,jm) = uaf(2:imm1,jmm1)

            case ( bc3POINTSMOOTH )
              uaf(2:imm1,jm) = ( uaf(1:imm2,jmm1)       &
                               + uaf(2:imm1,jmm1)       &
                               + uaf(3:im  ,jmm1) )/3.

            case ( bcCLAMPED )
              uaf(2:imm1,jm) = UA_bry%NTH(2:imm1,1)

            case ( bcFLATHER ) ! ??? (See east boundary)
!              uaf(2:imm1,jm) = uaf(2:imm1,jmm2)
              uaf(2:imm1,jm) = uaf(2:imm1,jmm2)                     &
                             + (uaf(2:imm1,jmm1)-uaf(2:imm1,jmm2))  !&
!                              *dy(2:imm1,jmm1)/dy(2:imm1,jmm2)

            case default
              uaf(2:imm1,jm) = 0.

          end select

        end if

! SOUTH
        if ( hasSOUTH ) then

          select case ( BC % VEL2D % NORM % SOUTH )

            case ( bc0GRADIENT )
              vaf(2:imm1,2) = vaf(2:imm1,3)

            case ( bc3POINTSMOOTH )
              vaf(2:imm1,2) = ( vaf(1:imm2,3)       &
                              + vaf(2:imm1,3)       &
                              + vaf(3:im  ,3) )/3.

            case ( bcCLAMPED )
              vaf(2:imm1,2) = VA_bry%STH(2:imm1,1)

            case ( bcFLATHER )
              vaf(2:imm1,2) = VA_bry%STH(2:imm1,1)              &
                           - sqrt(GRAV/h(2:imm1,2))             &
                            *(el(2:imm1,2)-EL_bry%STH(2:imm1,1))

            case ( bcGENFLATHER )
              vaf(2:imm1,2) = ( VA_bry%STH(2:imm1,1)                 &
                               *(h(2:imm1,2)+el(2:imm1,2))           &
                               -sqrt(GRAV/h(2:imm1,2))               &
                               *(el(2:imm1,2)-EL_bry%STH(2:imm1,1))) &
                              /(h(2:imm1,1)+EL_bry%STH(2:imm1,1))

            case default
              vaf(2:imm1,2) = 0.

          end select

          vaf(2:imm1,2) = ramp*vaf(2:imm1,2)
          vaf(2:imm1,1) = vaf(2:imm1,2)

          select case ( BC % VEL2D % TANG % SOUTH )

            case ( bc0GRADIENT )
              uaf(2:imm1,1) = uaf(2:imm1,2)

            case ( bc3POINTSMOOTH )
              uaf(2:imm1,1) = ( uaf(1:imm2,2)       &
                              + uaf(2:imm1,2)       &
                              + uaf(3:im  ,2) )/3.

            case ( bcCLAMPED )
              uaf(2:imm1,1) = UA_bry%STH(2:imm1,1)

            case ( bcFLATHER ) ! ??? (See east boundary)
!              uaf(2:imm1,1) = uaf(2:imm1,3)
              uaf(2:imm1,1) = uaf(2:imm1,2)                   &
                             + (uaf(2:imm1,2)-uaf(2:imm1,3))  !&
!                              *dy(2:imm1,2)/dy(2:imm1,3)

            case default
              uaf(2:imm1,1) = 0.

          end select

        end if

      end if ! periodic_y end


! Apply u- and v-masks
      uaf = uaf*dum
      vaf = vaf*dvm


    end ! subroutine bc_vel_ext
!______________________________________________________________________
!
    subroutine bc_vel_int
!----------------------------------------------------------------------
!  Apply internal (3D) velocity boundary conditions.
!______________________________________________________________________
!
      use glob_const , only: SMALL
      use glob_grid  , only: dum, dvm, h
      use glob_ocean , only: u, ub, uf, v, vb, vf, wubot, wvbot

      implicit none

      real(rk), dimension(im,2) :: grdx
      real(rk), dimension(2,jm) :: grdy
      real(rk)                  :: cff, cx, cy, dvdt, dvdx, dvdy


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi2d_mpi(wubot,im,jm)
        call xperi2d_mpi(wvbot,im,jm)
        call xperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        if ( hasNORTH ) then
          if ( BC % VEL3D % TANG % NORTH == bc0GRADIENT ) then
            wubot(:,jm) = wubot(:,jmm1)
            do k=1,kbm1
              uf(:,jm,k) = uf(:,jmm1,k)
            end do
          end if
        end if
        if ( hasSOUTH ) then
          if ( BC % VEL3D % TANG % SOUTH == bc0GRADIENT ) then
            wubot(:,1) = wubot(:,2)
            do k=1,kbm1
              uf(:,1,k) = uf(:,2,k)
            end do
          end if
        end if

      else
! EAST
        if ( hasEAST ) then

          select case ( BC % VEL3D % NORM % EAST )

            case ( bc0GRADIENT )
              uf(im,2:jmm1,1:kbm1) = u(imm1,2:jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(im,2:jmm1,1:kbm1) = ( u(imm1,1:jmm2,1:kbm1)       &
                                     + u(imm1,2:jmm1,1:kbm1)       &
                                     + u(imm1,3:jm  ,1:kbm1) )/3.

            case ( bcCLAMPED )
              uf(im,2:jmm1,1:kbm1) = U_bry%EST(1,2:jmm1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = sqrt( h(im,j) / hmax )
                  uf(im,j,k) = .25 * (     cff  * ( u(imm1,j-1,k)     &
                                 + 2.*u(imm1,j,k) + u(imm1,j+1,k) )   &
                                     + (1.-cff) * ( u(im  ,j-1,k)     &
                                 + 2.*u(im  ,j,k) + u(im  ,j+1,k) ) )
                end do
              end do

            case ( bcORLANSKI )
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = uf(im-1,j,k) + ub(im-1,j,k) - 2.*u(im-2,j,k)
                  if ( abs(cff) < .01 ) cff = sign(.01_rk,cff)
                  cff = ( ub(im-1,j,k)-uf(im-1,j,k) )/cff
                  if ( cff > 1. ) cff = 1.
                  if ( cff < 0. ) cff = 0.
                  uf(im,j,k) = ( (1.-cff)*ub(im  ,j,k)    &
                                + 2.*cff * u(im-1,j,k) )  &
                               / (1.+cff)
                end do
              end do

              case ( bcRADIATION_ENH )
              do k = 1,kbm1
                do j = 1,jmm1
                  grdy(1,j) = u(imm1,j+1,k)-u(imm1,j,k)
                  grdy(2,j) = u(im  ,j+1,k)-u(im  ,j,k)
                end do
                do j = 2,jmm1
                  dvdt =  u(imm1,j,k)-uf(imm1,j,k)
                  dvdx = uf(imm1,j,k)-uf(imm2,j,k)
                  if ( dvdt*dvdx < 0. ) dvdt = 0.
                  if ( dvdt*(grdy(1,j-1)+grdy(1,j)) > 0. ) then
                    dvdy = grdy(1,j-1)
                  else
                    dvdy = grdy(1,j  )
                  end if
                  cff = max(dvdx*dvdx + dvdy*dvdy, small)
                  cx  = dvdt*dvdx
                  cy  = min(cff,max(dvdt*dvdy,-cff))
                  uf(im,j,k) = ( cff*u(im,j,k) + cx*uf(imm1,j,k)      &
                                - max(cy,0.)*grdy(2,j-1)              &
                                - min(cy,0.)*grdy(2,j  ) )            &
                               / (cff+cx)
                end do
              end do

            case default
              uf(im,:,:) = 0.

          end select

          select case ( BC % VEL3D % TANG % EAST )

            case ( bc0GRADIENT )
              vf(im,2:jmm1,1:kbm1) = v(imm1,2:jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              vf(im,2:jmm1,1:kbm1) = ( v(imm1,1:jmm2,1:kbm1)       &
                                     + v(imm1,2:jmm1,1:kbm1)       &
                                     + v(imm1,3:jm  ,1:kbm1) )/3.

            case ( bcCLAMPED )
              vf(im,2:jmm1,1:kbm1) = V_bry%EST(1,2:jmm1,1:kbm1)

            case ( bcRADIATION ) ! NOT CORRECT!
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = sqrt( h(im,j) / hmax )
                  vf(im,j,k) = .25 * (     cff  * ( v(imm1,j-1,k)     &
                                 - 2.*v(imm1,j,k) + v(imm1,j+1,k) )   &
                                     + (1.-cff) * ( v(im  ,j-1,k)     &
                                 - 2.*v(im  ,j,k) + v(im  ,j+1,k) ) )
                end do
              end do

            case default
              vf(im,:,:) = 0.

          end select

        end if

! WEST
        if ( hasWEST ) then

          select case ( BC % VEL3D % NORM % WEST )

            case ( bc0GRADIENT )
              uf(2,2:jmm1,1:kbm1) = u(3,2:jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(2,2:jmm1,1:kbm1) = ( u(2,1:jmm2,1:kbm1)       &
                                    + u(2,2:jmm1,1:kbm1)       &
                                    + u(2,3:jm  ,1:kbm1) )/3.

            case ( bcCLAMPED )
              uf(2,2:jmm1,1:kbm1) = U_bry%WST(1,2:jmm1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = sqrt( h(1,j) / hmax )
                  uf(2,j,k) = .25 * (     cff  * ( u(3,j-1,k)     &
                                   + 2.*u(3,j,k) + u(3,j+1,k) )   &
                                    + (1.-cff) * ( u(2,j-1,k)     &
                                   + 2.*u(2,j,k) + u(2,j+1,k) ) )
                  uf(1,j,k) = uf(2,j,k)
                end do
              end do

            case ( bcORLANSKI )
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = uf(3,j,k) + ub(3,j,k) - 2.*u(4,j,k)
                  if ( abs(cff) < .01 ) cff = sign(.01_rk,cff)
                  cff = ( ub(3,j,k)-uf(3,j,k) )/cff
                  if ( cff > 1. ) cff = 1.
                  if ( cff < 0. ) cff = 0.
                  uf(2,j,k) = ( (1.-cff)*ub(2,j,k)    &
                               + 2.*cff * u(3,j,k) )  &
                              / (1.+cff)
                  uf(1,j,k) = uf(2,j,k)
                end do
              end do

            case ( bcRADIATION_ENH )
              do k = 1,kbm1
                do j = 1,jmm1
                  grdy(1,j) = u(2,j+1,k)-u(2,j,k)
                  grdy(2,j) = u(3,j+1,k)-u(3,j,k)
                end do
                do j = 2,jmm1
                  dvdt =  u(2,j,k)-uf(2,j,k)
                  dvdx = uf(2,j,k)-uf(3,j,k)
                  if ( dvdt*dvdx < 0. ) dvdt = 0.
                  if ( dvdt*(grdy(2,j-1)+grdy(2,j)) > 0. ) then
                    dvdy = grdy(2,j-1)
                  else
                    dvdy = grdy(2,j  )
                  end if
                  cff = max(dvdx*dvdx + dvdy*dvdy, small)
                  cx  = dvdt*dvdx
                  cy  = min(cff,max(dvdt*dvdy,-cff))
                  uf(2,j,k) = ( cff*u(2,j,k) + cx*uf(3,j,k)           &
                               - max(cy,0.)*grdy(1,j-1)               &
                               - min(cy,0.)*grdy(1,j  ) )             &
                              / (cff+cx)
                  uf(1,j,k) = uf(2,j,k)
                end do
              end do

            case default
              uf(1:2,:,:) = 0.

          end select

          select case ( BC % VEL3D % TANG % WEST )

            case ( bc0GRADIENT )
              vf(1,2:jmm1,1:kbm1) = v(2,2:jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              vf(1,2:jmm1,1:kbm1) = ( v(2,1:jmm2,1:kbm1)       &
                                    + v(2,2:jmm1,1:kbm1)       &
                                    + v(2,3:jm  ,1:kbm1) )/3.

            case ( bcCLAMPED )
              vf(1,2:jmm1,1:kbm1) = V_bry%WST(1,2:jmm1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do j = 2,jmm1
                  cff = sqrt( h(1,j) / hmax )
                  vf(1,j,k) = .25 * (     cff  * ( v(2,j-1,k)     &
                                   + 2.*u(2,j,k) + v(2,j+1,k) )   &
                                    + (1.-cff) * ( v(1,j-1,k)     &
                                   + 2.*u(1,j,k) + v(1,j+1,k) ) )
                end do
              end do

            case default
              vf(1,:,:) = 0.

          end select

        end if

      end if ! periodic_x end

! Apply periodic BC in y-dimension
      if ( periodic_bc % y ) then

        call yperi2d_mpi(wubot,im,jm)
        call yperi2d_mpi(wvbot,im,jm)
        call yperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        if ( hasEAST ) then
          if ( BC % VEL3D % TANG % EAST == bc0GRADIENT ) then
            wvbot(im,:) = wvbot(imm1,:)
            do k=1,kbm1
              vf(im,:,k) = vf(imm1,:,k)
            end do
          end if
        end if
        if ( hasWEST ) then
          if ( BC % VEL3D % TANG % WEST == bc0GRADIENT ) then
            wvbot(1,:) = wvbot(2,:)
            do k=1,kbm1
              vf(1,:,k) = vf(2,:,k)
            end do
          end if
        end if

      else
! NORTH
        if ( hasNORTH ) then

          select case ( BC % VEL3D % NORM % NORTH )

            case ( bc0GRADIENT )
              vf(2:imm1,jm,1:kbm1) = v(2:imm1,jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              vf(2:imm1,jm,1:kbm1) = ( v(1:imm2,jmm1,1:kbm1)       &
                                     + v(2:imm1,jmm1,1:kbm1)       &
                                     + v(3:im  ,jmm1,1:kbm1) )/3.

            case ( bcCLAMPED )
              vf(2:imm1,jm,1:kbm1) = V_bry%NTH(2:imm1,1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = sqrt( h(i,jm) / hmax )
                  vf(i,jm,k) = .25 * (     cff  * ( v(i-1,jmm1,k)     &
                                 + 2.*v(i,jmm1,k) + v(i+1,jmm1,k) )   &
                                     + (1.-cff) * ( v(i-1,jm  ,k)     &
                                 + 2.*v(i,jm  ,k) + v(i+1,jm  ,k) ) )
                end do
              end do

            case ( bcORLANSKI )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = vf(i,jm-1,k) + vb(i,jm-1,k) - 2.*v(i,jm-2,k)
                  if ( abs(cff) < .01 ) cff = sign(.01_rk,cff)
                  cff = ( vb(i,jm-1,k)-vf(i,jm-1,k) )/cff
                  if ( cff > 1. ) cff = 1.
                  if ( cff < 0. ) cff = 0.
                  vf(i,jm,k) = ( (1.-cff)*vb(i,jm  ,k)    &
                                + 2.*cff * v(i,jm-1,k) )  &
                               / (1.+cff)
                end do
              end do

            case ( bcRADIATION_ENH )
              do k = 1,kbm1
                do i = 2,im
                  grdx(i,1) = v(i,jmm1,k)-v(i-1,jmm1,k)
                  grdx(i,2) = v(i,jm  ,k)-v(i-1,jm  ,k)
                end do
                do i = 2,imm1
                  dvdt =  v(i,jmm1,k)-vf(i,jmm1,k)
                  dvdy = vf(i,jmm1,k)-vf(i,jmm2,k)
                  if ( dvdt*dvdy < 0. ) dvdt = 0.
                  if ( dvdt*(grdx(i,1)+grdx(i+1,1)) > 0. ) then
                    dvdx = grdx(i  ,1)
                  else
                    dvdx = grdx(i+1,1)
                  end if
                  cff = max(dvdx*dvdx + dvdy*dvdy, small)
                  cx  = min(cff,max(dvdt*dvdx,-cff))
                  cy  = dvdt*dvdy
                  vf(i,jm,k) = ( cff*v(i,jm,k) + cy*vf(i,jmm1,k)      &
                                - max(cx,0.)*grdx(i  ,2)              &
                                - min(cx,0.)*grdx(i+1,2) )            &
                               / (cff+cy)
                end do
              end do

            case default
              vf(:,jm,:) = 0.

          end select

          select case ( BC % VEL3D % TANG % NORTH )

            case ( bc0GRADIENT )
              uf(2:imm1,jm,1:kbm1) = u(2:imm1,jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(2:imm1,jm,1:kbm1) = ( u(1:imm2,jmm1,1:kbm1)       &
                                     + u(2:imm1,jmm1,1:kbm1)       &
                                     + u(3:im  ,jmm1,1:kbm1) )/3.

            case ( bcCLAMPED )
              uf(2:imm1,jm,1:kbm1) = U_bry%NTH(2:imm1,1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = sqrt( h(i,jm) / hmax )
                  uf(i,jm,k) = .25 * (     cff  * ( u(i-1,jmm1,k)     &
                                 + 2.*u(i,jmm1,k) + u(i+1,jmm1,k) )   &
                                     + (1.-cff) * ( u(i-1,jm  ,k)     &
                                 + 2.*u(i,jm  ,k) + u(i+1,jm  ,k) ) )
                end do
              end do

            case default
              uf(:,jm,:) = 0.

          end select

        end if

! SOUTH
        if ( hasSOUTH ) then

          select case ( BC % VEL3D % NORM % SOUTH )

            case ( bc0GRADIENT )
              vf(2:imm1,2,1:kbm1) = v(2:imm1,3,1:kbm1)

            case ( bc3POINTSMOOTH )
              vf(2:imm1,2,1:kbm1) = ( v(1:imm2,3,1:kbm1)       &
                                    + v(2:imm1,3,1:kbm1)       &
                                    + v(3:im  ,3,1:kbm1) )/3.

            case ( bcCLAMPED )
              vf(2:imm1,2,1:kbm1) = V_bry%STH(2:imm1,1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = sqrt( h(i,1) / hmax )
                  vf(i,2,k) = .25 * (     cff  * ( v(i-1,3,k)     &
                                   + 2.*v(i,3,k) + v(i+1,3,k) )   &
                                    + (1.-cff) * ( v(i-1,2,k)     &
                                   + 2.*v(i,2,k) + v(i+1,2,k) ) )
                end do
              end do

            case ( bcORLANSKI )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = vf(i,3,k) + vb(i,3,k) - 2.*v(i,4,k)
                  if ( abs(cff) < 0. ) cff = sign(.01_rk,cff)
                  cff = ( vb(i,3,k)-vf(i,3,k) )/cff
                  if ( cff > 1. ) cff = 1.
                  if ( cff < 0. ) cff = 0.
                  vf(i,2,k) = ( (1.-cff)*vb(i,2,k)    &
                               + 2.*cff * v(i,3,k) )  &
                              / (1.+cff)
                  vf(i,1,k) = vf(i,2,k)
                end do
              end do

            case ( bcRADIATION_ENH )
              do k = 1,kbm1
                do i = 2,im
                  grdx(i,1) = v(i,2,k)-v(i-1,2,k)
                  grdx(i,2) = v(i,3,k)-v(i-1,3,k)
                end do
                do i = 2,imm1
                  dvdt =  v(i,3,k)-vf(i,3,k)
                  dvdy = vf(i,3,k)-vf(i,4,k)
                  if ( dvdt*dvdy < 0. ) dvdt = 0.
                  if ( dvdt*(grdx(i,2)+grdx(i+1,2)) > 0. ) then
                    dvdx = grdx(i  ,2)
                  else
                    dvdx = grdx(i+1,2)
                  end if
                  cff = max(dvdx*dvdx + dvdy*dvdy, small)
                  cx  = min(cff,max(dvdt*dvdx,-cff))
                  cy  = dvdt*dvdy
                  vf(i,2,k) = ( cff*v(i,2,k) + cy*vf(i,3,k)           &
                               - max(cx,0.)*grdx(i  ,1)               &
                               - min(cx,0.)*grdx(i+1,1) )             &
                              / (cff+cy)
                  vf(i,1,k) = vf(i,2,k)
                end do
              end do

            case default
              vf(:,1:2,:) = 0.

          end select

          select case ( BC % VEL3D % TANG % SOUTH )

            case ( bc0GRADIENT )
              uf(2:imm1,1,1:kbm1) = u(2:imm1,2,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(2:imm1,1,1:kbm1) = ( u(1:imm2,2,1:kbm1)       &
                                    + u(2:imm1,2,1:kbm1)       &
                                    + u(3:im  ,2,1:kbm1) )/3.

            case ( bcCLAMPED )
              uf(2:imm1,1,1:kbm1) = U_bry%STH(2:imm1,1,1:kbm1)

            case ( bcRADIATION )
              do k = 1,kbm1
                do i = 2,imm1
                  cff = sqrt( h(i,1) / hmax )
                  uf(i,1,k) = .25 * (     cff  * ( u(i-1,2,k)     &
                                   + 2.*u(i,2,k) + u(i+1,2,k) )   &
                                    + (1.-cff) * ( u(i-1,1,k)     &
                                   + 2.*u(i,1,k) + u(i+1,1,k) ) )
                end do
              end do

            case default
              uf(:,1,:) = 0.

          end select

        end if

      end if ! periodic_y end


! Apply u- nd v-masks
      do k = 1,kbm1
        uf(:,:,k) = uf(:,:,k)*dum
        vf(:,:,k) = vf(:,:,k)*dvm
      end do


    end ! subroutine bc_vel_int
!______________________________________________________________________
!
    subroutine bc_ts
!----------------------------------------------------------------------
!  Apply temperature and salinity boundary conditions.
!
!   uf : temperature
!   vf : salinity
!______________________________________________________________________
!
      use glob_grid  , only: dx, dy, fsm, zz
      use glob_ocean , only: dt, s, t, u, uf, v, vf, w
      use model_run  , only: dti

      implicit none

      integer  ii, jj
      real(rk) u1, wm


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

      else
! EAST
        if ( hasEAST ) then

          select case ( BC % TS % EAST )

            case ( bc0GRADIENT )
              uf(im,:,1:kbm1) = t(imm1,:,1:kbm1)
              vf(im,:,1:kbm1) = s(imm1,:,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(im,2:jmm1,1:kbm1) = ( t(imm1,1:jmm2,1:kbm1)       &
                                     + t(imm1,2:jmm1,1:kbm1)       &
                                     + t(imm1,3:jm  ,1:kbm1) )/3.
              uf(im,1     ,1:kbm1) = ( t(imm1,1     ,1:kbm1)       &
                                     + t(imm1,2     ,1:kbm1) )*.5
              uf(im,  jm  ,1:kbm1) = ( t(imm1,  jmm1,1:kbm1)       &
                                     + t(imm1,  jm  ,1:kbm1) )*.5
              vf(im,2:jmm1,1:kbm1) = ( s(imm1,1:jmm2,1:kbm1)       &
                                     + s(imm1,2:jmm1,1:kbm1)       &
                                     + s(imm1,3:jm  ,1:kbm1) )/3.
              vf(im,1     ,1:kbm1) = ( s(imm1,1     ,1:kbm1)       &
                                     + s(imm1,2     ,1:kbm1) )*.5
              vf(im,  jm  ,1:kbm1) = ( s(imm1,  jmm1,1:kbm1)       &
                                     + s(imm1,  jm  ,1:kbm1) )*.5

            case ( bcCLAMPED )
              uf(im,:,1:kbm1) = T_bry%EST(1,:,1:kbm1)
              vf(im,:,1:kbm1) = S_bry%EST(1,:,1:kbm1)

            case ( bcINOUTFLOW )
              do k = 1,kbm1
                do j = 1,jm
                  u1 = 2.*u(im,j,k)*dti / (dx(im,j)+dx(imm1,j))
                  if ( u1 <= 0. ) then
                    uf(im,j,k) = t(im,j,k)                        &
                               - u1*(T_bry%EST(1,j,k)-t(im,j,k))
                    vf(im,j,k) = s(im,j,k)                        &
                               - u1*(S_bry%EST(1,j,k)-s(im,j,k))
                  else
                    uf(im,j,k) = t(im,j,k)                        &
                               - u1*(t(im,j,k)-t(imm1,j,k))
                    vf(im,j,k) = s(im,j,k)                        &
                               - u1*(s(im,j,k)-s(imm1,j,k))
                    if ( k/=1 .and. k/=kbm1 ) then
                      wm = .5 * (w(imm1,j,k)+w(imm1,j,k+1))*dti   &
                              / ( (zz(k-1)-zz(k+1))*dt(imm1,j) )
                      uf(im,j,k) = uf(im,j,k)                     &
                           - wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                      vf(im,j,k) = vf(im,j,k)                     &
                           - wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
                    end if
                  end if
                end do
              end do

              if ( NFE > 3 ) then
                do k = 1,kbm1
                  do j = 1,jm
                    do i = 1,NFE
                      ii = im-i+1
                      uf(ii,j,k) =    uf(ii,j,k)*(1.-frz(ii,j))  &
                                 + (T_bry%EST(i,j,k)*frz(ii,j))
                      vf(ii,j,k) =    vf(ii,j,k)*(1.-frz(ii,j))  &
                                 + (S_bry%EST(i,j,k)*frz(ii,j))
                    end do
                  end do
                end do
              end if

          end select

        end if

! WEST
        if ( hasWEST ) then

          select case ( BC % TS % WEST )

            case ( bc0GRADIENT )
              uf(1,:,1:kbm1) = t(2,:,1:kbm1)
              vf(1,:,1:kbm1) = s(2,:,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(1,2:jmm1,1:kbm1) = ( t(2,1:jmm2,1:kbm1)       &
                                    + t(2,2:jmm1,1:kbm1)       &
                                    + t(2,3:jm  ,1:kbm1) )/3.
              uf(1,1     ,1:kbm1) = ( t(2,1     ,1:kbm1)       &
                                    + t(2,2     ,1:kbm1) )*.5
              uf(1,  jm  ,1:kbm1) = ( t(2,  jmm1,1:kbm1)       &
                                    + t(2,  jm  ,1:kbm1) )*.5
              vf(1,2:jmm1,1:kbm1) = ( s(2,1:jmm2,1:kbm1)       &
                                    + s(2,2:jmm1,1:kbm1)       &
                                    + s(2,3:jm  ,1:kbm1) )/3.
              vf(1,1     ,1:kbm1) = ( s(2,1     ,1:kbm1)       &
                                    + s(2,2     ,1:kbm1) )*.5
              vf(1,  jm  ,1:kbm1) = ( s(2,  jmm1,1:kbm1)       &
                                    + s(2,  jm  ,1:kbm1) )*.5

            case ( bcCLAMPED )
              uf(1,:,1:kbm1) = T_bry%WST(1,:,1:kbm1)
              vf(1,:,1:kbm1) = S_bry%WST(1,:,1:kbm1)

            case ( bcINOUTFLOW )
              do k = 1,kbm1
                do j = 1,jm
                  u1 = 2.*u(2,j,k)*dti / (dx(1,j)+dx(2,j))
                  if ( u1 >= 0. ) then
                    uf(1,j,k) = t(1,j,k)                        &
                              - u1*(t(1,j,k)-T_bry%WST(1,j,k))
                    vf(1,j,k) = s(1,j,k)                        &
                              - u1*(s(1,j,k)-S_bry%WST(1,j,k))
                  else
                    uf(1,j,k) = t(1,j,k)                        &
                              - u1*(t(2,j,k)-t(1,j,k))
                    vf(1,j,k) = s(1,j,k)                        &
                              - u1*(s(2,j,k)-s(1,j,k))
                    if ( k/=1 .and. k/=kbm1 ) then
                      wm = .5 * ( w(2,j,k)+w(2,j,k+1) )*dti     &
                              / ( (zz(k-1)-zz(k+1))*dt(2,j) )
                      uf(1,j,k) = uf(1,j,k)                     &
                           - wm*(t(2,j,k-1)-t(2,j,k+1))
                      vf(1,j,k) = vf(1,j,k)                     &
                           - wm*(s(2,j,k-1)-s(2,j,k+1))
                    end if
                  end if
                end do
              end do

              if ( NFW > 3 ) then
                do k = 1,kbm1
                  do j = 1,jm
                    do i = 1,NFW
                      uf(i,j,k) =     uf(i,j,k)*(1.-frz(i,j))  &
                                + (T_bry%WST(i,j,k)*frz(i,j))
                      vf(i,j,k) =     vf(i,j,k)*(1.-frz(i,j))  &
                                + (S_bry%WST(i,j,k)*frz(i,j))
                    end do
                  end do
                end do
              end if

          end select

        end if

      end if ! periodic_x end

! Apply periodic BC in y-dimension
      if ( periodic_bc % y ) then

        call yperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

      else
! NORTH
        if ( hasNORTH ) then

          select case ( BC % TS % NORTH )

            case ( bc0GRADIENT )
              uf(:,jm,1:kbm1) = t(:,jmm1,1:kbm1)
              vf(:,jm,1:kbm1) = s(:,jmm1,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(2:imm1,jm,1:kbm1) = ( t(1:imm2,jmm1,1:kbm1)       &
                                     + t(2:imm1,jmm1,1:kbm1)       &
                                     + t(3:im  ,jmm1,1:kbm1) )/3.
              uf(1     ,jm,1:kbm1) = ( t(1     ,jmm1,1:kbm1)       &
                                     + t(2     ,jmm1,1:kbm1) )*.5
              uf(  im  ,jm,1:kbm1) = ( t(imm1  ,jmm1,1:kbm1)       &
                                     + t(im    ,jmm1,1:kbm1) )*.5
              vf(2:imm1,jm,1:kbm1) = ( s(1:imm2,jmm1,1:kbm1)       &
                                     + s(2:imm1,jmm1,1:kbm1)       &
                                     + s(3:im  ,jmm1,1:kbm1) )/3.
              vf(1     ,jm,1:kbm1) = ( s(1     ,jmm1,1:kbm1)       &
                                     + s(2     ,jmm1,1:kbm1) )*.5
              vf(  im  ,jm,1:kbm1) = ( s(imm1  ,jmm1,1:kbm1)       &
                                     + s(im    ,jmm1,1:kbm1) )*.5

            case ( bcCLAMPED )
              uf(:,jm,1:kbm1) = T_bry%NTH(:,1,1:kbm1)
              vf(:,jm,1:kbm1) = S_bry%NTH(:,1,1:kbm1)

            case ( bcINOUTFLOW )
              do k = 1,kbm1
                do i = 1,im
                  u1 = 2.*v(i,jm,k)*dti / (dy(i,jm)+dy(i,jmm1))
                  if ( u1 <= 0. ) then
                    uf(i,jm,k) = t(i,jm,k)                        &
                               - u1*(T_bry%NTH(i,1,k)-t(i,jm,k))
                    vf(i,jm,k) = s(i,jm,k)                        &
                               - u1*(S_bry%NTH(i,1,k)-s(i,jm,k))
                  else
                    uf(i,jm,k) = t(i,jm,k)                        &
                               - u1*(t(i,jm,k)-t(i,jmm1,k))
                    vf(i,jm,k) = s(i,jm,k)                        &
                               - u1*(s(i,jm,k)-s(i,jmm1,k))
                    if ( k/=1 .and. k/=kbm1 ) then
                      wm = .5 * (w(i,jmm1,k)+w(i,jmm1,k+1))*dti   &
                              / ( (zz(k-1)-zz(k+1))*dt(i,jmm1) )
                      uf(i,jm,k) = uf(i,jm,k)                     &
                                 - wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                      vf(i,jm,k) = vf(i,jm,k)                     &
                                 - wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
                    end if
                  end if
                end do
              end do

              if ( NFN > 3 ) then
                do k = 1,kbm1
                  do i = 1,im
                    do j = 1,NFN
                      jj = jm-j+1
                      uf(i,jj,k) =    uf(i,jj,k)*(1.-frz(i,jj))  &
                                 + (T_bry%NTH(i,j,k)*frz(i,jj))
                      vf(i,jj,k) =    vf(i,jj,k)*(1.-frz(i,jj))  &
                                 + (S_bry%NTH(i,j,k)*frz(i,jj))
                    end do
                  end do
                end do
              end if

          end select

        end if

! SOUTH
        if ( hasSOUTH ) then

          select case ( BC % TS % SOUTH )

            case ( bc0GRADIENT )
              uf(:,1,1:kbm1) = t(:,2,1:kbm1)
              vf(:,1,1:kbm1) = s(:,2,1:kbm1)

            case ( bc3POINTSMOOTH )
              uf(2:imm1,1,1:kbm1) = ( t(1:imm2,2,1:kbm1)       &
                                    + t(2:imm1,2,1:kbm1)       &
                                    + t(3:im  ,2,1:kbm1) )/3.
              uf(1     ,1,1:kbm1) = ( t(1     ,2,1:kbm1)       &
                                    + t(2     ,2,1:kbm1) )*.5
              uf(  im  ,1,1:kbm1) = ( t(imm1  ,2,1:kbm1)       &
                                    + t(im    ,2,1:kbm1) )*.5
              vf(2:imm1,1,1:kbm1) = ( s(1:imm2,2,1:kbm1)       &
                                    + s(2:imm1,2,1:kbm1)       &
                                    + s(3:im  ,2,1:kbm1) )/3.
              vf(1     ,1,1:kbm1) = ( s(1     ,2,1:kbm1)       &
                                    + s(2     ,2,1:kbm1) )*.5
              vf(  im  ,1,1:kbm1) = ( s(imm1  ,2,1:kbm1)       &
                                    + s(im    ,2,1:kbm1) )*.5

            case ( bcCLAMPED )
              uf(:,1,1:kbm1) = T_bry%STH(:,1,1:kbm1)
              vf(:,1,1:kbm1) = S_bry%STH(:,1,1:kbm1)

            case ( bcINOUTFLOW )
              do k = 1,kbm1
                do i = 1,im
                  u1 = 2.*v(i,2,k)*dti / (dy(i,1)+dy(i,2))
                  if ( u1 >= 0. ) then
                    uf(i,1,k) = t(i,1,k)                        &
                              - u1*(t(i,1,k)-T_bry%STH(i,1,k))
                    vf(i,1,k) = s(i,1,k)                        &
                              - u1*(s(i,1,k)-S_bry%STH(i,1,k))
                  else
                    uf(i,1,k) = t(i,1,k)                        &
                              - u1*(t(i,2,k)-t(i,1,k))
                    vf(i,1,k) = s(i,1,k)                        &
                              - u1*(s(i,2,k)-s(i,1,k))
                    if ( k/=1 .and. k/=kbm1 ) then
                      wm = .5 * ( w(i,2,k)+w(i,2,k+1) )*dti     &
                              / ( (zz(k-1)-zz(k+1))*dt(i,2) )
                      uf(i,1,k) = uf(i,1,k)                     &
                                - wm*(t(i,2,k-1)-t(i,2,k+1))
                      vf(i,1,k) = vf(i,1,k)                     &
                                - wm*(s(i,2,k-1)-s(i,2,k+1))
                    end if
                  end if
                end do
              end do

              if ( NFS > 3 ) then
                do k = 1,kbm1
                  do i = 1,im
                    do j = 1,NFS
                      uf(i,j,k) =    uf(i,j,k)*(1.-frz(i,j))  &
                                 + (T_bry%STH(i,j,k)*frz(i,j))
                      vf(i,j,k) =    vf(i,j,k)*(1.-frz(i,j))  &
                                 + (S_bry%STH(i,j,k)*frz(i,j))
                    end do
                  end do
                end do
              end if

          end select

        end if

      end if ! periodic_y end

! Allpy rho-mask
      do k = 1,kbm1
        uf(:,:,k) = uf(:,:,k)*fsm
        vf(:,:,k) = vf(:,:,k)*fsm
      end do


    end ! subroutine bc_ts
!______________________________________________________________________
!
    subroutine bc_vel_vert
!----------------------------------------------------------------------
!  Apply vertical velocity boundary conditions.
!______________________________________________________________________

      use glob_grid  , only: fsm
      use glob_ocean , only: w

      implicit none


! Apply rho-mask
      do k = 1,kbm1
        w(:,:,k) = w(:,:,k)*fsm
      end do


    end ! subroutine bc_vel_vert
!______________________________________________________________________
!
    subroutine bc_turb
!----------------------------------------------------------------------
!  Apply turbulent boundary conditions.
!
!   uf : q2
!   vf : q2l
!______________________________________________________________________
!
      use glob_const , only: SMALL
      use glob_grid  , only: dx, dy, fsm
      use glob_ocean , only: kh, km, kq, l, q2, q2l, u, uf, v, vf
      use model_run  , only: dti

      implicit none

      real(rk) u1


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(kh(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(km(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi(kq(:,:,1:kbm1),im,jm,kbm1)
        call xperi3d_mpi( l(:,:,1:kbm1),im,jm,kbm1)

      else
! EAST
        if ( hasEAST ) then

          do k = 1,kb
            do j = 1,jm
              u1 = 2.*u(im,j,k)*dti / (dx(im,j)+dx(imm1,j))
              if ( u1 <= 0. ) then
                uf(im,j,k) = q2(im,j,k)  - u1*(SMALL-q2(im,j,k))
                vf(im,j,k) = q2l(im,j,k) - u1*(SMALL-q2l(im,j,k))
              else
                uf(im,j,k) = q2(im,j,k)                      &
                           - u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k) = q2l(im,j,k)                     &
                           - u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end do
          end do

        end if

! WEST
        if ( hasWEST ) then

          do k = 1,kb
            do j = 1,jm
              u1 = 2.*u(2,j,k)*dti / (dx(1,j)+dx(2,j))
              if ( u1 >= 0. ) then
                uf(1,j,k) = q2(1,j,k)  - u1*(q2(1,j,k)-SMALL)
                vf(1,j,k) = q2l(1,j,k) - u1*(q2l(1,j,k)-SMALL)
              else
                uf(1,j,k) = q2(1,j,k)  - u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k) = q2l(1,j,k) - u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end do
          end do

        end if

      end if


! Apply periodic BC in y-dimension
      if ( periodic_bc % y ) then

        call yperi3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(kh(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(km(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi(kq(:,:,1:kbm1),im,jm,kbm1)
        call yperi3d_mpi( l(:,:,1:kbm1),im,jm,kbm1)

      else
! NORTH
        if ( hasNORTH ) then

          do k = 1,kb
            do i = 1,im
              u1 = 2.*v(i,jm,k)*dti / (dy(i,jm)+dy(i,jmm1))
              if ( u1 <= 0. ) then
                uf(i,jm,k) = q2(i,jm,k)  - u1*(SMALL-q2(i,jm,k))
                vf(i,jm,k) = q2l(i,jm,k) - u1*(SMALL-q2l(i,jm,k))
              else
                uf(i,jm,k) = q2(i,jm,k)                      &
                           - u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k) = q2l(i,jm,k)                     &
                           - u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end do
          end do

        end if

! SOUTH
        if ( hasSOUTH ) then

          do k=1,kb
            do i=1,im
              u1 = 2.*v(i,2,k)*dti / (dy(i,1)+dy(i,2))
              if ( u1 >= 0. ) then
                uf(i,1,k) = q2(i,1,k)  - u1*(q2(i,1,k)-SMALL)
                vf(i,1,k) = q2l(i,1,k) - u1*(q2l(i,1,k)-SMALL)
              else
                uf(i,1,k) = q2(i,1,k)  - u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k) = q2l(i,1,k) - u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end do
          end do

        end if

      end if


! Apply rho-mask
      do k=1,kb
        uf(:,:,k) = uf(:,:,k)*fsm + 1.e-10
        vf(:,:,k) = vf(:,:,k)*fsm + 1.e-10
      end do


    end ! subroutine bc_turb
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

    end ! function get_filename
!
!= I/O SECTION ========================================================
!______________________________________________________________________
!
    subroutine read_var_3d_nc( path, var_name, var          &
                             , year, start   , edge, ncid )
!----------------------------------------------------------------------
!  Read a variable (NC format).
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_3d

      integer                   , intent(inout) :: ncid
      integer                   , intent(in   ) :: year
      integer(MPI_OFFSET_KIND)                                 &
              , dimension(4)    , intent(in   ) :: start, edge
      real(rk), dimension(:,:,:), intent(  out) :: var
      character(len=*)          , intent(in   ) :: path, var_name

      integer                  varid, status
      character(len=256)       filename, netcdf_file, units


      if ( ncid <= 0 ) then
! open netcdf file
        if ( ncid == 0 ) then
          call check( nf90mpi_close( ncid )              &
                    , 'nfmpi_close:'//trim(netcdf_file) )
        end if

        filename = get_filename( path, year )
        netcdf_file = trim(filename)
        status = nf90mpi_open( POM_COMM, netcdf_file, NF_NOWRITE   &
                             , MPI_INFO_NULL, ncid )
        if ( status /= NF_NOERR ) then
          call msg_print("", 2, "Failed reading `"//trim(var_name)  &
                              //"` from `"//trim(filename)//"`")
          ncid = -1
          return
        end if
      end if

! get variable
      call check( nf90mpi_inq_varid( ncid, var_name, varid )  &
                , 'nfmpi_inq_varid: '//trim(var_name) )

! get data
      call check( get_var_real_3d                    &
                  ( ncid, varid, start, edge, var )  &
                , 'get_var_real: '//trim(var_name) )

! convert data if necessary
      status = nf90mpi_get_att( ncid, varid, "units", units )
      if ( status == NF_NOERR ) then
        select case ( trim(units) )
          case ( "cm", "cm/s", "cm s^-1" )
            var = var/100.
          case ( "mm", "mm/s", "mm s^-1" )
            var = var/1000.
        end select
      end if


    end ! subroutine read_var_3d_nc
!______________________________________________________________________
!
    subroutine read_var_2d_nc( path, var_name, var         &
                             , year, start   , edge, ncid )
!----------------------------------------------------------------------
!  Read a variable (NC format).
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain
      use mpi        , only: MPI_INFO_NULL, MPI_OFFSET_KIND
      use pnetcdf

      implicit none

      integer, external :: get_var_real_2d

      integer                   , intent(inout) :: ncid
      integer                   , intent(in   ) :: year
      integer(MPI_OFFSET_KIND)                                 &
              , dimension(4)    , intent(in   ) :: start, edge
      real(rk), dimension(:,:)  , intent(  out) :: var
      character(len=*)          , intent(in   ) :: path, var_name

      integer                  varid, status
      character(len=256)       filename, netcdf_file, units


      if ( ncid <= 0 ) then
! open netcdf file
        if ( ncid == 0 ) then
          call check( nf90mpi_close( ncid )              &
                    , 'nfmpi_close:'//trim(netcdf_file) )
        end if

        filename = get_filename( path, year )
        netcdf_file = trim(filename)
        status = nf90mpi_open( POM_COMM, netcdf_file, NF_NOWRITE   &
                             , MPI_INFO_NULL, ncid )
        if ( status /= NF_NOERR ) then
          call msg_print("", 2, "Failed reading `"//trim(var_name)  &
                              //"` from `"//trim(filename)//"`")
          ncid = -1
          return
        end if
      end if

! get variable
      call check( nf90mpi_inq_varid( ncid, var_name, varid )  &
                , 'nfmpi_inq_varid: '//trim(var_name) )

! get data
      call check( get_var_real_2d                    &
                  ( ncid, varid, start, edge, var )  &
                , 'get_var_real: '//trim(var_name) )

! convert data if necessary
      status = nf90mpi_get_att( ncid, varid, "units", units )
      if ( status == NF_NOERR ) then
        select case ( trim(units) )
          case ( "cm", "cm/s", "cm s^-1" )
            var = var/100.
          case ( "mm", "mm/s", "mm s^-1" )
            var = var/1000.
        end select
      end if


    end ! subroutine read_var_2d_nc
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
          print '(/a,a)', 'IO error at module `BRY`: ', routine
          print '("[",i4,"] ",a)', status, nf90mpi_strerror(status)
          stop
        end if
      end if


    end ! subroutine check


end module