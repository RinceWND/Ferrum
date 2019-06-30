!______________________________________________________________________
!
! Module `SEAICE` (seaice.f90)
!----------------------------------------------------------------------
!  Module for adjusting surface forcing with sea ice concentration.
! And more?
!
!  Author  : RinceWND
!  Created : 2019-02-06
!______________________________________________________________________
!
module seaice

  use glob_const, only: PATH_LEN, rk, VAR_LEN

  implicit none

  private

  public :: initialize_mod, init, step, ice_advance, advect_ice

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  real(rk), parameter ::     &  !
    Cp     = 2113.           &  ! Specific heat capacity of "pure" ice
  , Cdwi   = 5.1e-3          &  ! Water-ice drag coefficent
!  , Cdwi   = 3.1e-3          &  ! Water-ice drag coefficent
  , Ch     =     .004        &  ! Heat transfer coefficient between ocean and sea ice (very rough approximation)
  , RHOiref=  920.              ! Reference density of sea ice (typically, 720-940 kg m^-3)

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED ! Is set according to the external flag `use_ice`.

  integer, private :: read_int     ! interval for reading (days)
  real   , private :: a         &  ! time-interpolation factor
                    , relx         ! relaxation period

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

  logical          &
    interp_ice       ! Interpolate in time

!----------------------------------------------------------------------
! Paths configuration
!----------------------------------------------------------------------
  character(len=PATH_LEN)  & ! Full paths (with filenames) to:
    ice_path                 !  sea ice concetration

!----------------------------------------------------------------------
! Input variables' names
!----------------------------------------------------------------------
  character(len=VAR_LEN)  &
    icec_name             & ! Ice concentration
  , iceh_name               ! Ice thickness estimate

!----------------------------------------------------------------------
! Sea ice related arrays
!----------------------------------------------------------------------
  real(rk)                   &
       , public              &
       , allocatable         &
       , dimension(:,:)   :: &
    icec                     & ! ice concentration (fraction)
  , icec_b                   & ! ice concentration (fraction) at time n-1
  , icec_f                   & ! ice concentration (fraction) at time n+1
  , iceh                     & ! ice thickness (m)
  , iceu                     & ! ice x-velocity (m/s)
  , iceu_b                   & ! ice x-velocity (m/s) at time n-1
  , iceu_f                   & ! ice x-velocity (m/s) at time n+1
  , icev                     & ! ice y-velocity (m/s)
  , icev_b                   & ! ice y-velocity (m/s) at time n-1
  , icev_f                   & ! ice y-velocity (m/s) at time n+1
  , tau_ia_u                 & ! u-component of ice-air stress (N/m^2)
  , tau_ia_v                 & ! v-component of ice-air stress (N/m^2)
  , tau_iw_u                 & ! u-component of ice-water stress (N/m^2)
  , tau_iw_v                 & ! v-component of ice-water stress (N/m^2)
  , itsurf                     ! ice-water heat flux

!----------------------------------------------------------------------
! Private variables
!----------------------------------------------------------------------
  real(rk)                   &
       , allocatable         &
       , dimension(:,:,:) :: &
    icec_int                 & !
  , iceh_int                   !


  contains
!______________________________________________________________________
!
    subroutine initialize_mod( config_file )
!----------------------------------------------------------------------
!  Initialize ice module.
!______________________________________________________________________
!
      use config     , only: use_ice
      use glob_domain, only: is_master

      implicit none

      character(len=*), intent(in) :: config_file

      integer pos

      namelist/ice/           &
        ice_path, read_int

      namelist/ice_vars/      &
        icec_name, iceh_name


      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_ICE ) DISABLED = .true.

! Initialize variables with their defaults
      read_int = 3600 * 24 ! 86400 (Daily)
      relx     = 1./60.!(.125*86400.) !1./(5.*86400.) ! 5-day relaxation

      ice_path = "in/surf/"

      icec_name = "icec"
      iceh_name = "iceh"
      
      interp_ice = .true.

! Override configuration
      open ( 73, file = config_file, status = 'old' )
      read ( 73, nml = ice )
      read ( 73, nml = ice_vars )
      close( 73 )

! Variables management
      pos = len(trim(ice_path))
      if ( ice_path(pos:pos) == "/" ) then
        ice_path = trim(ice_path)//"ice."
      end if

      if ( is_master ) then
        print *, "Sea ice data: ", trim(ice_path)
      end if

! Allocate necessary arrays
      call allocate_arrays

      call msg_print("SEAICE MODULE INITIALIZED", 1, "")


    end ! subroutine initialize_mod
!
!______________________________________________________________________
!
    subroutine allocate_arrays
!----------------------------------------------------------------------
!  Allocate necessary variables.
!______________________________________________________________________
!
!      use config     , only: calc_wind
      use glob_domain, only: im, jm

      implicit none

      integer(1) N ! Interpolation array extension size
                   ! The structure is following:
                   !   var at n   step: var(:,:,1)
                   !   var at n-1 step: var(:,:,2)
                   !   var at n+1 step: var(:,:,3)


! Allocate core arrays
      allocate(        &
        icec(im,jm)    &
      , icec_b(im,jm)  &
      , icec_f(im,jm)  &
      , iceh(im,jm)    &
      , iceu(im,jm)    &
      , iceu_b(im,jm)  &
      , iceu_f(im,jm)  &
      , icev(im,jm)    &
      , icev_b(im,jm)  &
      , icev_f(im,jm)  &
      , tau_ia_u(im,jm)&
      , tau_ia_v(im,jm)&
      , tau_iw_u(im,jm)&
      , tau_iw_v(im,jm)&
      , itsurf(im,jm)  &
       )

! Initialize mandatory arrays
      icec   = 0.
      icec_b = 0.
      icec_f = 0.
      iceh   =  .35
      iceu   = 0.
      iceu_b = 0.
      iceu_f = 0.
      icev   = 0.
      icev_b = 0.
      icev_f = 0.
      itsurf = 0.
      tau_ia_u = 0.
      tau_ia_v = 0.
      tau_iw_u = 0.
      tau_iw_v = 0.

! Allocate optional arrays
      N = 1
      if ( interp_ice ) N = 2

      allocate(           &
        icec_int(im,jm,2:N+1) &
      , iceh_int(im,jm,2:N+1) &
       )
      icec_int = 0.
      iceh_int = 0.


    end ! subroutine allocate_air
!
!______________________________________________________________________
!
    subroutine init( d_in )
!----------------------------------------------------------------------
!  Reads ice fields before experiment's start.
!______________________________________________________________________
!
!      use glob_domain, only: is_master
      use module_time
      use grid       , only: fsm
      use glob_ocean , only: u, v

      implicit none

      type(date), intent(in) :: d_in

      integer            max_in_prev, max_in_this, ncid
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      write(44,'(a4,20a12)') "i","c","cf","u","uf","v","vf","ui","uif","vi","vif"  &
                           ,"tauau","tauav","tauwu","tauwv","h","uw","vw","el","el+","el-"
! Set ncid to -1 to open first file. Then set it to 0 to close previous file and open new one. TODO: Stupid, I know.
      ncid = -1

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
        record(2) = max_in_prev ! TODO: [ NEEDS TESTING ]
        year(2) = d_in%year - 1
      elseif ( record(3) == max_in_this ) then
        record(3) = 1
        year(3) = d_in%year + 1
      end if

      record(1) = record(1) + 1

! Read ice fields
      if ( interp_ice ) then
        call read_all( .true., 2, year, record )
        call read_all( .true., 3, year, record )
        icec_int(:,:,2) = fsm*icec_int(:,:,2)
        icec_int(:,:,3) = fsm*icec_int(:,:,3)
        icec = ( 1. - a ) * icec_int(:,:,2) + a * icec_int(:,:,3)
        iceh = ( 1. - a ) * iceh_int(:,:,2) + a * iceh_int(:,:,3)
      else
        call read_all( .true., 1, year, record )
        icec = fsm*icec_int(:,:,1)
      end if

      iceu = .23*u(:,:,1)
      icev = .23*v(:,:,1)

      call msg_print("SEAICE INITIALIZED", 2, "")


    end ! subroutine init
!
!______________________________________________________________________
!
    subroutine step( d_in, uwnd, vwnd, u, v, wtsurf, swrad, tair )
!----------------------------------------------------------------------
!  Reads forcing fields during experiment.
!______________________________________________________________________
!
      use glob_const , only: rhoref
      use config     , only: tbias
      use glob_domain, only: im, jm!, is_master
      use grid       , only: fsm
      use module_time
      use model_run  , only: dti, sec_of_year
      use glob_ocean , only: t

      implicit none

      type(date), intent(in) :: d_in

      real(rk), dimension(im,jm), intent(in) ::  &
                         swrad, tair, wtsurf, u, uwnd, v, vwnd

      logical            ADVANCE_REC, ADVANCE_REC_INT
      integer            max_in_prev, max_in_this, secs
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      secs = sec_of_year

      year = d_in%year

      max_in_this = max_chunks_in_year( d_in%year  , read_int )
      max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
      chunk     = chunk_of_year( d_in, read_int )
      record(1) = int(chunk)

      if ( secs - int(real(record(1))*read_int) < dti ) then ! TODO: test this one as well.
        ADVANCE_REC = .true.
      else
        ADVANCE_REC = .false.
      end if

      secs = secs - int((real(record(1))+.5)*read_int)

      if ( chunk - record(1) < .5 ) then ! TODO: test (real - int) = ?
        record(2) = record(1)
        a = chunk - record(1) + .5
      else
        record(2) = record(1) + 1
        a = chunk - record(1) - .5
      end if
      record(3) = record(2) + 1

      if ( record(2) == 0 ) then
        record(2) = max_in_prev + 1 ! TODO: [ NEEDS TESTING ]
        year(2) = d_in%year - 1
      elseif ( record(3) == max_in_this + 1 ) then
        record(3) = 1
        year(3) = d_in%year + 1
      end if

      if ( secs >= 0 .and. secs < dti ) then
        ADVANCE_REC_INT = .true.
        if ( interp_ice ) then
          icec_int(:,:,2) = icec_int(:,:,3)
          iceh_int(:,:,2) = iceh_int(:,:,3)
        end if
      else
        ADVANCE_REC_INT = .false.
      end if

      record(1) = record(1) + 1

      if ( interp_ice ) then
        call read_all( ADVANCE_REC_INT, 3, year, record )
        icec_int(:,:,3) = fsm*icec_int(:,:,3)
        icec = icec + relx*( ( 1. - a ) * icec_int(:,:,2) + a * icec_int(:,:,3) - icec )
        iceh = ( 1. - a ) * iceh_int(:,:,2) + a * iceh_int(:,:,3)
      else
        call read_all( ADVANCE_REC    , 1, year, record )
        icec_int(:,:,1) = fsm*icec_int(:,:,1)
        icec = icec + relx*( icec_int(:,:,1) - icec )
      end if

!      call advect_ice( uwnd, vwnd, u, v, wtsurf, swrad, tair )
      itsurf = rhoiref/rhoref * Ch * ( t(:,:,1) + tbias + 1.82 )
!      print *, minval(itsurf), maxval(itsurf) ,":"
!      stop


    end ! subroutine step

!==============================================================
! Advance ice in time
!--------------------------------------------------------------
    subroutine ice_advance( wusurf, wvsurf, wtsurf, swrad, tair )

      use glob_const , only: grav, rhoref, SMALL
      use glob_domain, only: im, imm1, jm, jmm1, kb            &
                           , n_east, n_north, n_south, n_west
      use grid       , only: art, cor, dum, dvm, dx, dy, fsm
      use glob_ocean , only: el, t, u, v
      use model_run  , only: dte

      implicit none

      real(rk), dimension(im,jm), intent(in) :: wusurf, wvsurf  &
                                              , wtsurf, swrad, tair

      real(rk), dimension(im,jm) ::             &
        fx, fy, divu, pice, delx, dely          &
      , rhoi, duvi, uidx, vidy, fluxcx, fluxcy  &
      , ht

      real(rk) eeta, tmp, lhf, ti, tai, wtice, dh
      integer i,j, cnum
      real(rk), parameter :: turb_coef = 0.9

        ! FREEZ'T
      lhf = 323500. ! Latent heat of fusion [J/kg]
      ht  = (1.-icec)*(wtsurf-swrad)*4.1876d6-icec*(-4.+1.82)*2.13!/iceh ! Net heat flux [W/m2]
      do j = 1, jm
        do i = 1, im
          ti = .5*(tair(i,j)-1.82)
          tai= .5*(tair(i,j)-ti  )
          wtice = 0.
          if ( iceh(i,j) > 0.01 ) then
            wtice = 2.1/(900.*2.113+17.2e3*4./(tair(i,j)**2))*(tair(i,j)-2.*tai+ti)/(iceh(i,j)**2)
            wtice = wtice - 2.1/(2.113+17.2e3)*(ti*ti-2.*ti+1.82*1.82)/(iceh(i,j)**2)
          end if
          dh = 0.
          if ( icec(i,j) > 0. ) then
            dh = wtice*dte/lhf/930./(art(i,j)*icec(i,j))
          end if
          dh = dh + ht(i,j)*dte/lhf/930./(art(i,j)*(1-icec(i,j))*turb_coef)
          if ( dh > 0. ) then
            if ( t(i,j,1) < -1.82 ) then
              iceh(i,j) = iceh(i,j) + dh
            end if
          else
            if ( iceh(i,j) < -dh ) then
              iceh(i,j) = 0.
              icec(i,j) = 0.
            else
              iceh(i,j) = iceh(i,j) + dh
            end if
          end if
          if ( ht(i,j) > 0. .and. t(i,j,1) <= -1.82 .or. iceh(i,j) > 0. ) then
            iceh(i,j) = iceh(i,j) + dh
!            t(i,j,1) = -1.82
          end if
          if ( iceh(i,j) < 0. ) iceh(i,j) = 0.
!          icec(i,j) = icec(i,j) + exp(-icec(i,j)*8.)*.0147
          if ( ht(i,j) > 0. .and. t(i,j,1)<=-1.82 ) then
            icec(i,j) = icec(i,j) + (1-icec(i,j))*turb_coef
          elseif ( ht(i,j) < 0. .and. iceh(i,j) > 0. ) then
            if ( iceh(i,j) > 0.05 ) then
              icec(i,j) = icec(i,j)
            else
              icec(i,j) = icec(i,j) - 0.01
            end if
          end if
          if ( abs(iceu(i,j)) < SMALL .and. abs(icev(i,j)) < SMALL ) then
            iceu(i,j) = u(i,j,1)
            icev(i,j) = v(i,j,1)
          end if
        end do
      end do
      t(:,:,kb) = iceh

      eeta = 1.e2 !1.01e-7 ! 1010 cm2/s? ! The source claims the coefficient equals to 10^10 cm2/s! This gives unreallistic Infinities.

      delx = 0.
      dely = 0.
      duvi = 0.
      uidx = 0.
      vidy = 0.
      fluxcx = 0.
      fluxcy = 0.

!        u(:,:,1) = 0.
!        v(:,:,1) =  .2
!        el = 0.

! Calculate sea surface elevation gradient
      delx(2:im,:) = dum(2:im,:)*2.*(el(2:im,:)-el(1:imm1,:))  &
                         /(dx(2:im,:)+dx(1:imm1,:))
      dely(:,2:jm) = dvm(:,2:jm)*2.*(el(:,2:jm)-el(:,1:jmm1))  &
                         /(dy(:,2:jm)+dy(:,1:jmm1))
      call exchange2d_mpi(delx,im,jm)
      call exchange2d_mpi(dely,im,jm)

! Estimate sea ice density to a unit area
!      rhoi = 900.*max(iceh,.05_rk)*max(icec,.1_rk)
      rhoi = 900.

! Get wind stress over ice-free water (convert it from m^2/s^2 to N/m^2)
      tau_ia_u = -wusurf*rhoref
      tau_ia_v = -wvsurf*rhoref

! Calculate water-ice stress
      duvi=abs(sqrt((iceu-u(1:im,1:jm,1))**2+(icev-v(1:im,1:jm,1))**2))

      tau_iw_u= fsm*5.5e-3*rhoref*(iceu-u(1:im,1:jm,1))*duvi
      tau_iw_v= fsm*5.5e-3*rhoref*(icev-v(1:im,1:jm,1))*duvi

! Compute ice concentration fluxes
      do j=1,jm
        do i=2,imm1
          if (iceu(i,j)<0.) then
            fluxcx(i,j) = icec(i  ,j)*iceu(i,j)*dy(i,j)
          else
            fluxcx(i,j) = icec(i-1,j)*iceu(i,j)*dy(i,j)
          end if
        end do
      end do
      do j=2,jmm1
        do i=1,im
          if (icev(i,j)<0.) then
            fluxcy(i,j) = icec(i,j  )*icev(i,j)*dx(i,j)
          else
            fluxcy(i,j) = icec(i,j-1)*icev(i,j)*dx(i,j)
          end if
        end do
      end do
      call exchange2d_mpi(fluxcx,im,jm)
      call exchange2d_mpi(fluxcy,im,jm)

! Apply ice concentration boundary conditions
      if (n_north==-1) then
        do i=1,im
          if (icev(i,jm)>0.) then
            fluxcy(i,jm) = icec(i,jm)*icev(i,jm)*dx(i,jm)
          else
!            fluxcy(i,jm) = cibn(i)  *icev(i,jm)*dx(i,jm)
            fluxcy(i,jm) = ( fluxcy(i,jm-1)  &
                           + fluxcy(i,jm-2)  &
                           + fluxcy(i,jm-3) ) / 3.
          end if
        end do
      end if
      if (n_east==-1) then
        do j=1,jm
          if (iceu(im,j)>0.) then
            fluxcx(im,j) = icec(im,j)*iceu(im,j)*dy(im,j)
          else
!            fluxcx(im,j) = cibw(j)  *iceu(im,j)*dy(im,j)
            fluxcx(im,j) = ( fluxcx(im-1,j)  &
                           + fluxcx(im-2,j)  &
                           + fluxcx(im-3,j) ) / 3.
          end if
        end do
      end if
      if (n_south==-1) then
        do i=1,im
          if (icev(i,1)<0.) then
            fluxcy(i,1) = icec(i,1)*icev(i,1)*dx(i,1)
          else
!            fluxcy(i,1) = cibs(i) *icev(i,1)*dx(i,1)
            fluxcy(i,1) = ( fluxcy(i,2)  &
                          + fluxcy(i,3)  &
                          + fluxcy(i,4) ) / 3.
          end if
        end do
      end if
      if (n_west==-1) then
        do j=1,jm
          if (iceu(1,j)<0.) then
            fluxcx(1,j) = icec(1,j)*iceu(1,j)*dy(1,j)
          else
!            fluxcx(1,j) = cibw(j) *iceu(1,j)*dy(1,j)
            fluxcx(1,j) = ( fluxcx(2,j)  &
                          + fluxcx(3,j)  &
                          + fluxcx(4,j) ) / 3.
          end if
        end do
      end if

! Calculate velocity gradients
      uidx(2:im,:) = (iceu(2:im,:)-iceu(1:imm1,:))/dx(2:im,:)
      vidy(:,2:jm) = (icev(:,2:jm)-icev(:,1:jmm1))/dy(:,2:jm)
      call exchange2d_mpi(uidx,im,jm)
      call exchange2d_mpi(vidy,im,jm)
! Apply boundaries
      if (n_west==-1) uidx(1,:) = 0. !uidx(2,:)
      if (n_south==-1) vidy(:,1) = 0. !vidy(:,2)

! Derive divergency
      divu = uidx+vidy

! Get internal stress
      pice = 0.
      where (divu<0.) pice = -10.*divu

! Start ice concentration advection
      if (.false.) then
        call advtC(icec,icec_f)
      else
        do j=1,jmm1
          do i=1,imm1
            tmp = (fluxcx(i,j)-fluxcx(i+1,j)  &
                  +fluxcy(i,j)-fluxcy(i,j+1))
            icec_f(i,j) = icec(i,j)+dte*tmp/art(i,j)
            ! if there will be no ice in the cell, correct output fluxes (velocities) only
            if (icec_f(i,j) < 0.) then
              icec_f(i,j) = 0.
              cnum = 0
              tmp = 0.
              if (fluxcx(i  ,j)<0.) then
                iceu(i  ,j) = 0.
                cnum = cnum+1
                tmp = tmp+fluxcx(i,j)
              end if
              if (fluxcx(i+1,j)>0.) then
                iceu(i+1,j) = 0.
                cnum = cnum+1
                tmp = tmp+fluxcx(i+1,j)
              end if
              if (fluxcy(i,j  )<0.) then
                icev(i,j  ) = 0.
                cnum = cnum+1
                tmp = tmp+fluxcy(i,j)
              end if
              if (fluxcy(i,j+1)>0.) then
                icev(i,j+1) = 0.
                cnum = cnum+1
                tmp = tmp+fluxcy(i,j+1)
              end if

              if (fluxcx(i  ,j)<0.) then
                fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
              end if
              if (fluxcx(i+1,j)>0.) then
                fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
              end if
              if (fluxcy(i,j  )<0.) then
                fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
              end if
              if (fluxcy(i,j+1)>0.) then
                fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
              end if
            else
            ! if the cell will be water-free, correct input fluxes (velocities) only
              if (icec_f(i,j) > 1.) then
                icec_f(i,j) = 1.
                cnum = 0
                tmp = 0.
                if (fluxcx(i  ,j)>0.) then
                  iceu(i  ,j) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcx(i,j)
                end if
                if (fluxcx(i+1,j)<0.) then
                  iceu(i+1,j) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcx(i+1,j)
                end if
                if (fluxcy(i,j  )>0.) then
                  icev(i,j  ) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcy(i,j)
                end if
                if (fluxcy(i,j+1)<0.) then
                  icev(i,j+1) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcy(i,j+1)
                end if

                if (fluxcx(i  ,j)>0.) then
                  fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
                end if
                if (fluxcx(i+1,j)<0.) then
                  fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
                end if
                if (fluxcy(i,j  )>0.) then
                  fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
                end if
                if (fluxcy(i,j+1)<0.) then
                  fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
                end if
              end if
            end if
          end do
        end do
        icec_f = icec_f*fsm
      end if
      call exchange2d_mpi(icec_f,im,jm)

      fx = 0.
      fx(2:im,2:jm) = 2.*eeta                                 &
                       *((uidx(2:im,2:jm)-uidx(1:imm1,2:jm))  &
                         /(dx(2:im,2:jm)+dx(1:imm1,2:jm))     &
                        +(vidy(2:im,2:jm)-vidy(2:im,1:jmm1))  &
                         /(dy(2:im,2:jm)+dy(2:im,1:jmm1)))
      fy = fx
      fx(2:im,2:jm) = fx(2:im,2:jm) +  &
              eeta*((divu(2:im,2:jm)-divu(1:imm1,2:jm))/dy(2:im,2:jm))
      fx(2:im,2:jm) = fx(2:im,2:jm) -                        &
                        (pice(2:im,2:jm)-pice(1:imm1,2:jm))  &
                        /dy(2:im,2:jm)

      fy(2:im,2:jm) = fy(2:im,2:jm) +  &
              eeta*((divu(2:im,2:jm)-divu(2:im,1:jmm1))/dx(2:im,2:jm))
      fy(2:im,2:jm) = fy(2:im,2:jm) -                        &
                        (pice(2:im,2:jm)-pice(2:im,1:jmm1))  &
                        /dx(2:im,2:jm)
      call exchange2d_mpi(fx,im,jm)
      call exchange2d_mpi(fy,im,jm)

      do j=1,jmm1
        do i=1,imm1
          iceu_f(i,j) = 0.
          icev_f(i,j) = 0.
          if (icec_f(i,j)>0.) then
            iceu_f(i,j) = iceu(i,j) + dte*                      &
                        ( cor(i,j)*(icev(i,j)+icev(i,j+1))      &
                         - grav*delx(i,j)                       &
                         + (tau_ia_u(i,j)-tau_iw_u(i,j))/rhoi(i,j)  &
                         + fx(i,j) )
            if (icec_f(i+1,j)==0.) then
              iceu_f(i+1,j) = iceu(i+1,j) + dte*                           &
                        ( cor(i+1,j)*(icev(i+1,j)+icev(i+1,j+1))          &
                         - grav*delx(i+1,j)                           &
                         + (tau_ia_u(i+1,j)-tau_iw_u(i+1,j))/rhoi(i+1,j)  &
                         + fx(i+1,j) )
            end if
            icev_f(i,j) = icev(i,j) + dte*                           &
                        (-cor(i,j)*(iceu(i,j)+iceu(i+1,j))          &
                         - grav*dely(i,j)                       &
                         + (tau_ia_v(i,j)-tau_iw_v(i,j))/rhoi(i,j)  &
                         + fy(i,j) )
            if (icec_f(i,j+1)==0.) then
              icev_f(i,j+1) = icev(i,j+1) + dte*                           &
                        (-cor(i,j+1)*(iceu(i,j+1)+iceu(i+1,j+1))          &
                         - grav*dely(i,j+1)                           &
                         + (tau_ia_v(i,j+1)-tau_iw_v(i,j+1))/rhoi(i,j+1)  &
                         + fy(i,j+1) )
            end if
          end if
        end do
      end do
      call exchange2d_mpi(icec_f,im,jm)
      call exchange2d_mpi(iceu_f,im,jm)
      call exchange2d_mpi(icev_f,im,jm)
      if (n_east==-1) then
        iceu_f(im,:) = iceu_f(imm1,:)
        icev_f(im,:) = icev_f(imm1,:)
      end if
      if (n_north==-1) then
        iceu_f(:,jm) = iceu_f(:,jmm1)
        icev_f(:,jm) = icev_f(:,jmm1)
      end if

      iceu_f = iceu_f*dum
      icev_f = icev_f*dvm
      
!      do j = 1, jm
!        do i = 1, im
!          if ( iceu_f(i,j) == iceu_f(i,j)+1. .or. &
!               abs(iceu_f(i,j)) > 1. ) then
!            print *, "Infinite ICE_UF @ ", i,j, iceu_f(i,j)
!            print *, "ICE_U: ", iceu(i,j)
!            print *, "DTE  : ", dte
!            print *, "COR  : ", cor(i,j)
!            print *, "ICE_V: ", icev(i,j), icev(i,j+1)
!            print *, "GRAV : ", grav
!            print *, "DELX : ", delx(i,j)
!            print *, "TAi_U: ", tau_ia_u(i,j)
!            print *, "TAw_U: ", tau_iw_u(i,j)
!            print *, "RHOi : ", rhoi(i,j)
!            print *, "FX   : ", fx(i,j)
!            print *, "U    : ", u(i,j,1)
!            print *, "DUVI : ", duvi(i,j)
!            stop
!          end if
!          if ( icev_f(i,j) == icev_f(i,j)+1. .or. &
!               abs(icev_f(i,j)) > 1. ) then
!            print *, "Infinite ICE_VF @ ", i,j, icev_f(i,j)
!            print *, "ICE_V: ", icev(i,j)
!            print *, "DTE  : ", dte
!            print *, "COR  : ", cor(i,j)
!            print *, "ICE_U: ", iceu(i,j), iceu(i,j+1)
!            print *, "GRAV : ", grav
!            print *, "DELY : ", dely(i,j)
!            print *, "TAi_V: ", tau_ia_v(i,j)
!            print *, "TAw_V: ", tau_iw_v(i,j)
!            print *, "RHOi : ", rhoi(i,j)
!            print *, "FY   : ", fy(i,j)
!            print *, "V    : ", v(i,j,1)
!            print *, "DUVI : ", duvi(i,j)
!            stop
!          end if
!        end do
!      end do

      iceu_b = iceu
      iceu   = iceu_f
      icev_b = icev
      icev   = icev_f
      icec_b = icec
      icec   = icec_f


    end ! subroutine ice_advance
!
!
!
    subroutine advect_ice( dt, uwnd, vwnd, u, v, wtsurf, swrad, tair )

      use config     , only: smoth
      use glob_const , only: grav, rhoref, SMALL
      use glob_domain, only: im, imm1, imm2, jm, jmm1, jmm2, kb  &
                           , n_east, n_north, n_south, n_west    &
                           , i_global, j_global, my_task
      use grid       , only: art, cor, dum, dvm, dx, dy, fsm
      use glob_ocean , only: el, t, uf,vf
      use model_run  , only: dtime, iint

      implicit none

      real(rk)                  , intent(in) :: dt
      real(rk), dimension(im,jm), intent(in) :: u, v, uwnd, vwnd  &
                                              , wtsurf, swrad, tair

      real(rk), dimension(im,jm) ::             &
        fx, fy, divu, pice, delx, dely          &
      , rhoi, duvi, uidx, vidy, fluxcx, fluxcy  &
      , ht

      real(rk) eeta, tmp, lhf, ti, tai, wtice, dh
      integer i,j, cnum


      eeta = 1.e4 !5.e5 !1.e2 !1.01e-7 ! 1010 cm2/s? ! The source claims the coefficient equals to 10^10 cm2/s! This gives unreallistic Infinities.

      delx = 0.
      dely = 0.
      duvi = 0.
      uidx = 0.
      vidy = 0.
      fluxcx = 0.
      fluxcy = 0.

!        u(:,:,1) = 0.
!        v(:,:,1) =  .2
!        el = 0.
      iceh = 1.

! Calculate sea surface elevation gradient
      delx(2:im,:) = dum(2:im,:)*2.*(el(2:im,:)-el(1:imm1,:))  &
                         /(dx(2:im,:)+dx(1:imm1,:))
      dely(:,2:jm) = dvm(:,2:jm)*2.*(el(:,2:jm)-el(:,1:jmm1))  &
                         /(dy(:,2:jm)+dy(:,1:jmm1))
      call exchange2d_mpi(delx,im,jm)
      call exchange2d_mpi(dely,im,jm)

! Estimate sea ice density to a unit area
      rhoi = rhoiref*max(iceh,.05_rk)!*max(icec,.1_rk)
!      rhoi = 900.

! Compute ice concentration fluxes
      do j=1,jm
        do i=2,imm1
          if (iceu(i,j)<0.) then
            fluxcx(i,j) = .5*icec(i  ,j)*iceu(i,j)*(dy(i,j)+dy(i-1,j))
          else
            fluxcx(i,j) = .5*icec(i-1,j)*iceu(i,j)*(dy(i,j)+dy(i-1,j))
          end if
        end do
      end do
      do j=2,jmm1
        do i=1,im
          if (icev(i,j)<0.) then
            fluxcy(i,j) = .5*icec(i,j  )*icev(i,j)*(dx(i,j)+dx(i,j-1))
          else
            fluxcy(i,j) = .5*icec(i,j-1)*icev(i,j)*(dx(i,j)+dx(i,j-1))
          end if
        end do
      end do
      call exchange2d_mpi(fluxcx,im,jm)
      call exchange2d_mpi(fluxcy,im,jm)

! Apply ice concentration boundary conditions
      if (n_north==-1) then
        do i=2,im
          if (icev(i,jm)>0.) then
            fluxcy(i,jm) = .5*icec(i,jm)*icev(i,jm)*(dx(i,jm)+dx(i-1,jm))
          else
!            fluxcy(i,jm) = cibn(i)  *icev(i,jm)*dx(i,jm)
            fluxcy(i,jm) = fluxcy(i,jmm1)
!            ( fluxcy(i,jm-1)  &
!                           + fluxcy(i,jm-2)  &
!                           + fluxcy(i,jm-3) ) / 3.
          end if
        end do
      end if
      if (n_east==-1) then
        do j=2,jm
          if (iceu(im,j)>0.) then
            fluxcx(im,j) = .5*icec(im,j)*iceu(im,j)*(dy(im,j)+dy(im,j-1))
          else
!            fluxcx(im,j) = cibw(j)  *iceu(im,j)*dy(im,j)
            fluxcx(im,j) = fluxcx(im-1,j)
!            ( fluxcx(im-1,j)  &
!                           + fluxcx(im-2,j)  &
!                           + fluxcx(im-3,j) ) / 3.
          end if
        end do
      end if
      if (n_south==-1) then
        do i=2,im
          if (icev(i,1)<0.) then
            fluxcy(i,1) = .5*icec(i,1)*icev(i,1)*(dx(i,1)+dx(i-1,2))
          else
!            fluxcy(i,1) = cibs(i) *icev(i,1)*dx(i,1)
            fluxcy(i,1) = fluxcy(i,2)
!            ( fluxcy(i,2)  &
!                          + fluxcy(i,3)  &
!                          + fluxcy(i,4) ) / 3.
          end if
        end do
      end if
      if (n_west==-1) then
        do j=2,jm
          if (iceu(1,j)<0.) then
            fluxcx(1,j) = .5*icec(1,j)*iceu(1,j)*(dy(1,j)+dy(1,j-1))
          else
!            fluxcx(1,j) = cibw(j) *iceu(1,j)*dy(1,j)
            fluxcx(1,j) = fluxcx(2,j)
!            ( fluxcx(2,j)  &
!                          + fluxcx(3,j)  &
!                          + fluxcx(4,j) ) / 3.
          end if
        end do
      end if
! TODO: is another exchange necessary?
      call exchange2d_mpi(fluxcx,im,jm)
      call exchange2d_mpi(fluxcy,im,jm)

! Start ice concentration advection
      if (.false.) then
        call advtC(icec,icec_f,fluxcx,fluxcy)
      else
        do j=1,jmm1
          do i=1,imm1
            tmp = (fluxcx(i,j)-fluxcx(i+1,j)  &
                  +fluxcy(i,j)-fluxcy(i,j+1))
            icec_f(i,j) = icec(i,j)+dt*tmp/art(i,j)
            ! if there will be no ice in the cell, correct output fluxes (velocities) only
            if (icec_f(i,j) < 0.) then
              icec_f(i,j) = 0.
!              cnum = 0
!              tmp = 0.
!              if (fluxcx(i  ,j)<0.) then
!                iceu(i  ,j) = 0.
!                cnum = cnum+1
!                tmp = tmp+fluxcx(i,j)
!              end if
!              if (fluxcx(i+1,j)>0.) then
!                iceu(i+1,j) = 0.
!                cnum = cnum+1
!                tmp = tmp+fluxcx(i+1,j)
!              end if
!              if (fluxcy(i,j  )<0.) then
!                icev(i,j  ) = 0.
!                cnum = cnum+1
!                tmp = tmp+fluxcy(i,j)
!              end if
!              if (fluxcy(i,j+1)>0.) then
!                icev(i,j+1) = 0.
!                cnum = cnum+1
!                tmp = tmp+fluxcy(i,j+1)
!              end if

!              if (fluxcx(i  ,j)<0.) then
!                fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
!              end if
!              if (fluxcx(i+1,j)>0.) then
!                fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
!              end if
!              if (fluxcy(i,j  )<0.) then
!                fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
!              end if
!              if (fluxcy(i,j+1)>0.) then
!                fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
!              end if
            else
            ! if the cell will be water-free, correct input fluxes (velocities) only
              if (icec_f(i,j) > 1.) then
                icec_f(i,j) = 1.
!                cnum = 0
!                tmp = 0.
!                if (fluxcx(i  ,j)>0.) then
!                  iceu(i  ,j) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcx(i,j)
!                end if
!                if (fluxcx(i+1,j)<0.) then
!                  iceu(i+1,j) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcx(i+1,j)
!                end if
!                if (fluxcy(i,j  )>0.) then
!                  icev(i,j  ) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcy(i,j)
!                end if
!                if (fluxcy(i,j+1)<0.) then
!                  icev(i,j+1) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcy(i,j+1)
!                end if
!
!                if (fluxcx(i  ,j)>0.) then
!                  fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
!                end if
!                if (fluxcx(i+1,j)<0.) then
!                  fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
!                end if
!                if (fluxcy(i,j  )>0.) then
!                  fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
!                end if
!                if (fluxcy(i,j+1)<0.) then
!                  fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
!                end if
              end if
            end if
          end do
        end do
        icec_f = icec_f*fsm
      end if
      call exchange2d_mpi(fluxcx,im,jm)
      call exchange2d_mpi(fluxcy,im,jm)
      call exchange2d_mpi(icec_f,im,jm)
      call exchange2d_mpi(iceu,im,jm)
      call exchange2d_mpi(icev,im,jm)

! Get wind stress over ice-free water (convert it from m^2/s^2 to N/m^2)
      duvi = sqrt((iceu-uwnd)**2+(icev-vwnd)**2)

      tau_ia_u = .0012*duvi*(uwnd-iceu)
      tau_ia_v = .0012*duvi*(vwnd-iceu)
!      tau_ia_u = -sign( 12.*(1.-exp(-.002*duvi*abs(iceu-uwnd))), iceu-uwnd )
!      tau_ia_v = -sign( 12.*(1.-exp(-.002*duvi*abs(icev-vwnd))), icev-vwnd )

! Calculate water-ice stress
      duvi = sqrt((iceu-u)**2+(icev-v)**2)

      tau_iw_u = Cdwi*rhoref*(iceu-u)*duvi
      tau_iw_v = Cdwi*rhoref*(icev-v)*duvi

! Calculate velocity gradients
      uidx(2:im,:) = (iceu(2:im,:)-iceu(1:imm1,:))/dx(2:im,:)
      vidy(:,2:jm) = (icev(:,2:jm)-icev(:,1:jmm1))/dy(:,2:jm)
      call exchange2d_mpi(uidx,im,jm)
      call exchange2d_mpi(vidy,im,jm)
! Apply boundaries
      if (n_west==-1) uidx(1,:) = 0. !uidx(2,:)
      if (n_south==-1) vidy(:,1) = 0. !vidy(:,2)

! Derive divergency
      divu = uidx+vidy

! Get internal stress
      pice = 0.
      where (divu<0.) pice = -10.*eeta*divu

      if (.false.) then
      fx = 0.
      fx(2:im,2:jm) = 2.*eeta                                 &
                       *((uidx(2:im,2:jm)-uidx(1:imm1,2:jm))  &
                         /(dx(2:im,2:jm)+dx(1:imm1,2:jm))     &
                        +(vidy(2:im,2:jm)-vidy(2:im,1:jmm1))  &
                         /(dy(2:im,2:jm)+dy(2:im,1:jmm1)))
      fy = fx
      fx(2:im,2:jm) = fx(2:im,2:jm) +  &
              eeta*((divu(2:im,2:jm)-divu(1:imm1,2:jm))/dy(2:im,2:jm))
      fx(2:im,2:jm) = fx(2:im,2:jm) -                        &
                        (pice(2:im,2:jm)-pice(1:imm1,2:jm))  &
                        /dy(2:im,2:jm)

      fy(2:im,2:jm) = fy(2:im,2:jm) +  &
              eeta*((divu(2:im,2:jm)-divu(2:im,1:jmm1))/dx(2:im,2:jm))
      fy(2:im,2:jm) = fy(2:im,2:jm) -                        &
                        (pice(2:im,2:jm)-pice(2:im,1:jmm1))  &
                        /dx(2:im,2:jm)
      else
        fx = 0.
        fx(2:imm1,:) = eeta/(dx(2:imm1,:)**2)  &
                       *(iceu(1:imm2,:)+iceu(3:im,:)  &
                         -2.*iceu(2:imm1,:))
        call exchange2d_mpi(fx,im,jm)
        fx(2:im,:) = fx(2:im,:) + eeta/dx(2:im,:)  &
                     *( divu(2:im,:)-divu(1:imm1,:) )
        fx(2:im,:) = fx(2:im,:) -                        &
                        ( pice(2:im,:)-pice(1:imm1,:) )  &
                        /dx(2:im,:)
        fy = 0.
        fy(:,2:jmm1) = eeta/(dy(:,2:jmm1)**2)  &
                       *(icev(:,1:jmm2)+icev(:,3:jm)  &
                         -2.*icev(:,2:jmm1))
        call exchange2d_mpi(fy,im,jm)
        fy(:,2:jm) = fy(:,2:jm) + eeta/dy(:,2:jm) &
                     *( divu(:,2:jm)-divu(:,1:jmm1) )
        fy(:,2:jm) = fy(:,2:jm) -                        &
                        ( pice(:,2:jm)-pice(:,1:jmm1) )  &
                        /dy(:,2:jm)
      end if
      call exchange2d_mpi(fx,im,jm)
      call exchange2d_mpi(fy,im,jm)
      if ( n_west == -1 ) fx( 1,:) = fx(   2,:)
      if ( n_east == -1 ) fx(im,:) = fx(imm1,:)
      if ( n_south == -1 ) fy(:, 1) = fy(:,   2)
      if ( n_north == -1 ) fy(:,jm) = fy(:,jmm1)

      do j=1,jmm1
        do i=1,imm1
          iceu_f(i,j) = 0.
          icev_f(i,j) = 0.
          if (icec_f(i,j)>0.) then
            iceu_f(i,j) = iceu(i,j) + dt*                      &
                        ( .5*cor(i,j)*(icev(i,j)+icev(i,j+1))      &
                         - grav*delx(i,j)                       &
                         + (tau_ia_u(i,j)-tau_iw_u(i,j))/rhoi(i,j)  &
                         + fx(i,j)                                &
                        )
!            if (icec_f(i+1,j)==0.) then
!              iceu_f(i+1,j) = iceu(i+1,j) + dt*                           &
!                        ( .5*cor(i,j)*(icev(i+1,j)+icev(i+1,j+1))          &
!                         - grav*delx(i+1,j)                           &
!                         + (tau_ia_u(i,j)-tau_iw_u(i,j))/rhoi(i,j)  &
!                         + fx(i+1,j)  &
!                        )
!            end if
            icev_f(i,j) = icev(i,j) + dt*                           &
                        (-.5*cor(i,j)*(iceu(i,j)+iceu(i+1,j))          &
                         - grav*dely(i,j)                       &
                         + (tau_ia_v(i,j)-tau_iw_v(i,j))/rhoi(i,j)  &
                         + fy(i,j)   &
                        )
!            if (icec_f(i,j+1)==0.) then
!              icev_f(i,j+1) = icev(i,j+1) + dt*                           &
!                        (-.5*cor(i,j+1)*(iceu(i,j+1)+iceu(i+1,j+1))          &
!                         - grav*dely(i,j+1)                           &
!                         + (tau_ia_v(i,j)-tau_iw_v(i,j))/rhoi(i,j)  &
!                         + fy(i,j+1)   &
!                        )
!            end if
          end if
        end do
      end do
      call exchange2d_mpi(icec_f,im,jm)
      call exchange2d_mpi(iceu_f,im,jm)
      call exchange2d_mpi(icev_f,im,jm)
      if (n_east==-1) then
        iceu_f(im,:) = iceu_f(imm1,:)
        icev_f(im,:) = icev_f(imm1,:)
      end if
      if (n_north==-1) then
        iceu_f(:,jm) = iceu_f(:,jmm1)
        icev_f(:,jm) = icev_f(:,jmm1)
      end if

      where ( icec_f < 1.e-6 )
        icec_f = 0.
      end where
      where ( icec_f(1:imm1,:) < 1.e-6 .and. icec_f(2:im,:) < 1.e-6 )
        iceu_f(1:imm1,:) = 0.
      end where
      where ( icec_f(:,1:jmm1) < 1.e-6 .and. icec_f(:,2:jm) < 1.e-6 )
        icev_f(:,1:jmm1) = 0.
      end where

      iceu_f = iceu_f*dum
      icev_f = icev_f*dvm

!      if ( my_task==4 ) then
!        i = 25 - i_global(1) + 1
!        j = 269 - j_global(1) + 1
!        print *, "== i,j == ", i_global(i), j_global(j)
!        print *, "u   : ", u(i,j), uf(i,j,1)
!        print *, "v   : ", v(i,j), vf(i,j,1)
!        print *, "ui  : ", iceu(i,j), iceu_f(i,j)
!        print *, "vi  : ", icev(i,j), icev_f(i,j)
!        print *, "tai : ", tau_ia_u(i,j), tau_ia_v(i,j)
!        print *, "twi : ", tau_iw_u(i,j), tau_iw_v(i,j)
!        print *, "Fi  : ", fx(i+1,j), fy(i,j+1)
!        write(44,'(i4,20f12.7)') iint,icec(i,j),icec_f(i,j),u(i,j),uf(i,j,1),v(i,j),vf(i,j,1)  &
!                     , iceu(i,j),iceu_f(i,j),icev(i,j),icev_f(i,j)  &
!                     , tau_ia_u(i,j),tau_ia_v(i,j),tau_iw_u(i,j),tau_iw_v(i,j)  &
!                     , iceh(i,j),uwnd(i,j),vwnd(i,j),el(i,j),el(i,j+1)  &
!                     , el(i+1,j)
!      end if

      do j = 1, jm
        do i = 1, im
          if ( iceu_f(i,j) == iceu_f(i,j)+1. .or. &
               abs(iceu_f(i,j)) > 100. ) then
            print *, "Infinite ICE_UF @ ", i_global(i),j_global(j), iceu_f(i,j), my_task
            print *, "ICE_U: ", iceu(i,j)
            print *, "DTE  : ", dt
            print *, "COR  : ", cor(i,j)
            print *, "ICE_V: ", icev(i,j), icev(i,j+1)
            print *, "GRAV : ", grav
            print *, "DELX : ", delx(i,j)
            print *, "TAi_U: ", tau_ia_u(i,j)
            print *, "TAw_U: ", tau_iw_u(i,j)
            print *, "RHOi : ", rhoi(i,j)
            print *, "FX   : ", fx(i,j)
            print *, "U    : ", u(i,j)
            print *, "DUVI : ", duvi(i,j)
            stop
          end if
          if ( icev_f(i,j) == icev_f(i,j)+1. .or. &
               abs(icev_f(i,j)) > 100. ) then
            print *, "Infinite ICE_VF @ ", i_global(i),j_global(j), icev_f(i,j), my_task
            print *, "ICE_V: ", icev(i,j)
            print *, "DTE  : ", dt
            print *, "COR  : ", cor(i,j)
            print *, "ICE_U: ", iceu(i,j), iceu(i,j+1)
            print *, "GRAV : ", grav
            print *, "DELY : ", dely(i,j)
            print *, "TAi_V: ", tau_ia_v(i,j)
            print *, "TAw_V: ", tau_iw_v(i,j)
            print *, "RHOi : ", rhoi(i,j)
            print *, "FY   : ", fy(i,j)
            print *, "V    : ", v(i,j)
            print *, "DUVI : ", duvi(i,j)
            stop
          end if
        end do
      end do

      iceu = iceu + .5*smoth*( iceu_b + iceu_f -2.*iceu )
      icev = icev + .5*smoth*( icev_b + icev_f -2.*icev )
      icec = icec + .5*smoth*( icec_b + icec_f -2.*icec )

      iceu_b = iceu
      iceu   = iceu_f
      icev_b = icev
      icev   = icev_f
      icec_b = icec
      icec   = icec_f

    end ! subroutine advect_ice
!
!______________________________________________________________________
!
    subroutine read_all( execute, n, year, record )

      implicit none

      logical              , intent(in) :: execute
      integer              , intent(in) :: n
      integer, dimension(3), intent(in) :: record, year

      integer            ncid
      character(len=128) desc


      if ( .not. execute ) return

      if ( n >= 2 ) then
        write(desc,'("Reading interp. sea ice record #",i4," @ ",i4)') &
            record(n), year(n)
      else
        write(desc,'("Reading sea ice record #",i4," @ ",i4)') &
            record(1), year(1)
      end if

      call msg_print("", 1, desc)

! Read wind file
      ncid = file_open_nc( ice_path, year(n) )

      if ( ncid /= -1 ) then

        if ((     interp_ice .and. n>=2) .or.       &
            (.not.interp_ice .and. n==1)      ) then
          call read_var_nc( icec_name, icec_int(:,:,n), record(n), ncid )
          if ( ncid == -1 ) then
            icec_int(:,:,n) = 0.
            call msg_print("", 2, "Sec ice concentration is set to zero.")
          end if
        end if

        call check( file_close_nc( ncid ), "nf_close" )

      end if


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
      use glob_const , only: rk
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
          case ( "%" )
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

end ! module seaice