!______________________________________________________________________
!
! Module `BRY` (bry.f90)
!----------------------------------------------------------------------
!  Module for applying boundary conditions.
!  Adopted from original POM `bcond` subroutine.
!
!  Author  : RinceWND
!  Created : 2018-09-02
!______________________________________________________________________

module bry

  use config     , only: rk ! SO... If a module uses some other module's variable, another module can access the variable through this "medium" module?
  use glob_domain, only: im, imm1, imm2, jm, jmm1, jmm2, kb, kbm1

  implicit none

!----------------------------------------------------------------------
! PRIVATE variables
!----------------------------------------------------------------------
  logical           &
       , private :: &
    hasEAST         & ! flag for eastern boundary
  , hasNORTH        & ! flag for northern boundary
  , hasSOUTH        & ! flag for southern boundary
  , hasWEST           ! flag for western boundary

  integer           &
       , private :: &
    i, j, k           ! simple counters

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
  , bcRADIATION       = 5    !  radiation condition

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


  contains

!______________________________________________________________________
!
    subroutine allocate_boundary
!----------------------------------------------------------------------
!  Allocates boundary arrays.
!______________________________________________________________________

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

    end subroutine allocate_boundary
!______________________________________________________________________
!
    subroutine initialize_boundary( config )
!----------------------------------------------------------------------
!  Initialize arrays for boundary conditions.
!______________________________________________________________________

      use glob_domain, only: n_east, n_north, n_south, n_west

      implicit none

      character(len=*), intent(in) :: config

      logical periodic_x, periodic_y

      namelist/bry_nml/                           &
        Hmax, periodic_x, periodic_y, USE_SPONGE

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


! Default BC types
!   elevation
      BC % zeta % EAST  = bc0GRADIENT
      BC % zeta % NORTH = bc0GRADIENT
      BC % zeta % SOUTH = bc0GRADIENT
      BC % zeta % WEST  = bc0GRADIENT
!   external velocity
      BC % vel2d % NORM % EAST  = bcFLATHER ! bcCLAMPED if RFE is 0 in old way
      BC % vel2d % NORM % NORTH = bcFLATHER
      BC % vel2d % NORM % SOUTH = bcFLATHER
      BC % vel2d % NORM % WEST  = bcFLATHER
      BC % vel2d % TANG % EAST  = bcCLAMPED
      BC % vel2d % TANG % NORTH = bcCLAMPED
      BC % vel2d % TANG % SOUTH = bcCLAMPED
      BC % vel2d % TANG % WEST  = bcCLAMPED
!   internal velocity
      BC % vel3d % NORM % EAST  = bcRADIATION
      BC % vel3d % NORM % NORTH = bcRADIATION
      BC % vel3d % NORM % SOUTH = bcRADIATION
      BC % vel3d % NORM % WEST  = bcRADIATION
      BC % vel3d % NORM % EAST  = bc3POINTSMOOTH
      BC % vel3d % NORM % NORTH = bc3POINTSMOOTH
      BC % vel3d % NORM % SOUTH = bc3POINTSMOOTH
      BC % vel3d % NORM % WEST  = bc3POINTSMOOTH
!   temp and salt
      BC % ts % EAST  = bcINOUTFLOW
      BC % ts % NORTH = bcINOUTFLOW
      BC % ts % SOUTH = bcINOUTFLOW
      BC % ts % WEST  = bcINOUTFLOW
!   vertical velocity (IGNORED anyway)
      BC % ts % EAST  = bcCLAMPED
      BC % ts % NORTH = bcCLAMPED
      BC % ts % SOUTH = bcCLAMPED
      BC % ts % WEST  = bcCLAMPED
!   turbulent parameters
      BC % ts % EAST  = bcINOUTFLOW
      BC % ts % NORTH = bcINOUTFLOW
      BC % ts % SOUTH = bcINOUTFLOW
      BC % ts % WEST  = bcINOUTFLOW

! Override configuration
      open ( 73, file = config, status = 'old' )
      read ( 73, nml = bry_nml )
      close( 73 )

! Manage derived type variables
      periodic_bc % x = periodic_x
      periodic_bc % y = periodic_y

! Allocate arrays
      call allocate_boundary


    end subroutine
!______________________________________________________________________
!
    subroutine initial_conditions_boundary
!----------------------------------------------------------------------
!  Sets up initial conditions for boundary arrays.
!______________________________________________________________________

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


    end subroutine initial_conditions_boundary

!______________________________________________________________________
!
    subroutine bc_zeta
!----------------------------------------------------------------------
!  Apply external (2D) elevation boundary conditions.
!______________________________________________________________________

      use glob_grid , only: fsm
      use glob_ocean, only: elf

      implicit none


! Apply periodic BC in x-dimension
      if ( periodic_bc % x ) then

        call xperi2d_mpi(elf,im,jm)
!...or...
      else
! EAST
        if ( hasEAST ) then

          select case ( BC % zeta % east )

            case default
              elf(im,:) = elf(imm1,:)

          end select

        end if
! WEST
        if ( hasWEST ) then

          select case ( BC % zeta % west )

            case ( bc0GRADIENT )
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

            case ( bc0GRADIENT )
              elf(:,jm) = elf(:,jmm1)

          end select

        end if
! SOUTH
        if ( hasSOUTH ) then

          select case ( BC % zeta % south )

            case ( bc0GRADIENT )
              elf(:, 1) = elf(:,   2)

          end select

        end if

      end if

! Apply rho-mask (TODO: needed?)
      elf = elf*fsm


      return

    end subroutine bc_zeta

!______________________________________________________________________
!
    subroutine bc_vel_ext
!----------------------------------------------------------------------
!  Apply external (2D) velocity boundary conditions.
!______________________________________________________________________

      use glob_const , only: GRAV
      use glob_grid  , only: dum, dvm, h
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

            case ( bcFLATHER ) ! NOT CORRECT!
              vaf(im,2:jmm1) = VA_bry%EST(1,2:jmm1)              &
                         - sqrt(GRAV/h(imm1,2:jmm1))             &
                          *(el(imm1,2:jmm1)-EL_bry%EST(1,2:jmm1))

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

          end select

          uaf(2,2:jmm1) = ramp*uaf(2,2:jmm1)
          uaf(1,2:jmm1) =      uaf(2,2:jmm1)

          select case ( BC % VEL2D % TANG % WEST )

            case ( bc0GRADIENT )
              vaf(2,2:jmm1) = vaf(3,2:jmm1)

            case ( bc3POINTSMOOTH )
              vaf(2,2:jmm1) = ( vaf(3,1:jmm2)       &
                              + vaf(3,2:jmm1)       &
                              + vaf(3,3:jm  ) )/3.

            case ( bcCLAMPED )
              vaf(2,2:jmm1) = VA_bry%WST(1,2:jmm1)

            case ( bcFLATHER ) ! NOT CORRECT!
              vaf(2,2:jmm1) = VA_bry%WST(1,2:jmm1)              &
                           + sqrt(GRAV/h(2,2:jmm1))             &
                            *(el(2,2:jmm1)-EL_bry%WST(1,2:jmm1))

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

            case ( bcFLATHER ) ! NOT CORRECT!
              uaf(2:imm1,jm) = UA_bry%NTH(2:imm1,1)              &
                         - sqrt(GRAV/h(2:imm1,jmm1))             &
                          *(el(2:imm1,jmm1)-EL_bry%NTH(2:imm1,1))

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

          end select

          vaf(2:imm1,2) = ramp*vaf(2:imm1,2)
          vaf(2:imm1,1) = vaf(2:imm1,2)

          select case ( BC % VEL2D % TANG % SOUTH )

            case ( bc0GRADIENT )
              uaf(2:imm1,2) = uaf(2:imm1,3)

            case ( bc3POINTSMOOTH )
              uaf(2:imm1,2) = ( uaf(1:imm2,3)       &
                              + uaf(2:imm1,3)       &
                              + uaf(3:im  ,3) )/3.

            case ( bcCLAMPED )
              uaf(2:imm1,2) = UA_bry%STH(2:imm1,1)

            case ( bcFLATHER ) ! NOT CORRECT!
              uaf(2:imm1,2) = UA_bry%STH(2:imm1,1)              &
                           + sqrt(GRAV/h(2:imm1,2))             &
                            *(el(2:imm1,2)-EL_bry%STH(2:imm1,1))

          end select

        end if

      end if ! periodic_y end


! Apply u- and v-masks
      uaf = uaf*dum
      vaf = vaf*dvm


      return

    end subroutine bc_vel_ext

!______________________________________________________________________
!
    subroutine bc_vel_int
!----------------------------------------------------------------------
!  Apply internal (3D) velocity boundary conditions.
!______________________________________________________________________

      use glob_grid  , only: dum, dvm, h
      use glob_ocean , only: u, uf, v, vf, wubot, wvbot

      implicit none

      real(rk) ga


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
                  ga = sqrt( h(im,j) / hmax )
                  uf(im,j,k) = .25 * (     ga  * ( u(imm1,j-1,k)     &
                                + 2.*u(imm1,j,k) + u(imm1,j+1,k) )   &
                                     + (1.-ga) * ( u(im  ,j-1,k)     &
                                + 2.*u(im  ,j,k) + u(im  ,j+1,k) ) )
                end do
              end do

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
                  ga = sqrt( h(im,j) / hmax )
                  vf(im,j,k) = .25 * (     ga  * ( v(imm1,j-1,k)     &
                                - 2.*v(imm1,j,k) + v(imm1,j+1,k) )   &
                                     + (1.-ga) * ( v(im  ,j-1,k)     &
                                - 2.*v(im  ,j,k) + v(im  ,j+1,k) ) )
                end do
              end do

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
                  ga = sqrt( h(1,j) / hmax )
                  uf(2,j,k) = .25 * (     ga  * ( u(3,j-1,k)     &
                                  + 2.*u(3,j,k) + u(3,j+1,k) )   &
                                    + (1.-ga) * ( u(2,j-1,k)     &
                                  + 2.*u(2,j,k) + u(2,j+1,k) ) )
                  uf(1,j,k) = uf(2,j,k)
                end do
              end do

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
                  ga = sqrt( h(1,j) / hmax )
                  vf(1,j,k) = .25 * (     ga  * ( v(2,j-1,k)     &
                                  + 2.*u(2,j,k) + v(2,j+1,k) )   &
                                    + (1.-ga) * ( v(1,j-1,k)     &
                                  + 2.*u(1,j,k) + v(1,j+1,k) ) )
                end do
              end do

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
                  ga = sqrt( h(i,jm) / hmax )
                  vf(i,jm,k) = .25 * (     ga  * ( v(i-1,jmm1,k)     &
                                + 2.*v(i,jmm1,k) + v(i+1,jmm1,k) )   &
                                     + (1.-ga) * ( v(i-1,jm  ,k)     &
                                + 2.*v(i,jm  ,k) + v(i+1,jm  ,k) ) )
                end do
              end do

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
                  ga = sqrt( h(i,jm) / hmax )
                  uf(i,jm,k) = .25 * (     ga  * ( u(i-1,jmm1,k)     &
                                + 2.*u(i,jmm1,k) + u(i+1,jmm1,k) )   &
                                     + (1.-ga) * ( u(i-1,jm  ,k)     &
                                + 2.*u(i,jm  ,k) + u(i+1,jm  ,k) ) )
                end do
              end do

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
                  ga = sqrt( h(i,1) / hmax )
                  vf(i,2,k) = .25 * (     ga  * ( v(i-1,3,k)     &
                                  + 2.*v(i,3,k) + v(i+1,3,k) )   &
                                    + (1.-ga) * ( v(i-1,2,k)     &
                                  + 2.*v(i,2,k) + v(i+1,2,k) ) )
                end do
              end do

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
                  ga = sqrt( h(i,1) / hmax )
                  uf(i,1,k) = .25 * (     ga  * ( u(i-1,2,k)     &
                                  + 2.*u(i,2,k) + u(i+1,2,k) )   &
                                    + (1.-ga) * ( u(i-1,1,k)     &
                                  + 2.*u(i,1,k) + u(i+1,1,k) ) )
                end do
              end do

          end select

        end if

      end if ! periodic_y end


! Apply u- nd v-masks
      do k = 1,kbm1
        uf(:,:,k) = uf(:,:,k)*dum
        vf(:,:,k) = vf(:,:,k)*dvm
      end do


      return

    end subroutine bc_vel_int


!______________________________________________________________________
!
    subroutine bc_ts
!----------------------------------------------------------------------
!  Apply temperature and salinity boundary conditions.
!
!   uf : temperature
!   vf : salinity
!______________________________________________________________________

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


      return

    end subroutine bc_ts


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


      return

    end subroutine bc_vel_vert

!______________________________________________________________________
!
    subroutine bc_turb
!----------------------------------------------------------------------
!  Apply turbulent boundary conditions.
!
!   uf : q2
!   vf : q2l
!______________________________________________________________________

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
        uf(:,:,k) = uf(:,:,k)*fsm
        vf(:,:,k) = vf(:,:,k)*fsm
      end do


      return

    end subroutine bc_turb


end module
