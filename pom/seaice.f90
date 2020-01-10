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

  public :: initialize_mod, init, step!, ice_advance

!----------------------------------------------------------------------
! Constants
!----------------------------------------------------------------------
  real(rk), parameter ::     &  !
    Cp     = 2113.           &  ! Specific heat capacity of "pure" ice
  , Ch     =     .004        &  ! Heat transfer coefficient between ocean and sea ice (very rough approximation)
  , RHOi   =  920.              ! Reference density of sea ice (typically, 720-940 kg m^-3)

!----------------------------------------------------------------------
! Configuration variables
!----------------------------------------------------------------------
  logical, private :: DISABLED ! Is set according to the external flag `use_ice`.

  integer, private :: read_int  ! interval for reading (days)
  real   , private :: a         ! time-interpolation factor

  character(len=10)                       &
       , parameter                        &
       , private   :: FORMAT_EXT = ".nc"

  logical          &
    interp_ice     & ! Interpolate in time
  , use_calendar     ! Use calendar to offset read record number from the start of a year

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
  , iceh                     & ! ice thickness (m)
  , iceu                     & ! ice x-velocity (m/s)
  , icev                     & ! ice y-velocity (m/s)
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

      namelist/ice/                       &
        ice_path, read_int, use_calendar

      namelist/ice_vars/      &
        icec_name, iceh_name


      DISABLED = .false.

! Configure module availability first
      if ( .not. USE_ICE ) DISABLED = .true.

! Initialize variables with their defaults
      read_int = 3600 * 24 ! 86400 (Daily)

      ice_path = "in/surf/"

      icec_name = "icec"
      iceh_name = "iceh"

      interp_ice   = .true.
      USE_CALENDAR = .true.

! Override configuration
      open  ( 73, file = config_file, status = 'old' )
      read  ( 73, nml = ice )
      rewind( 73 )
      read  ( 73, nml = ice_vars )
      close ( 73 )

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
      , iceh(im,jm)    &
      , iceu(im,jm)    &
      , icev(im,jm)    &
      , itsurf(im,jm)  &
       )

! Initialize mandatory arrays
      icec   = 0.
      iceh   =  .35
      iceu   = 0.
      icev   = 0.
      itsurf = 0.

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

      implicit none

      type(date), intent(in) :: d_in

      integer            max_in_prev, max_in_this, ncid
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

! Set ncid to -1 to open first file. Then set it to 0 to close previous file and open new one. TODO: Stupid, I know.
      ncid = -1

      year = d_in%year
      
! Decide on the record to read
      if ( USE_CALENDAR ) then

        max_in_this = max_chunks_in_year( d_in%year  , read_int )
        max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

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

      else

        record = [ 1, 1, 2 ]
        a = 0._rk

      end if

! Read ice fields
      if ( interp_ice ) then
        call read_all( .true., 2, year, record )
        call read_all( .true., 3, year, record )
      else
        call read_all( .true., 1, year, record )
      end if

      call msg_print("SEAICE INITIALIZED", 2, "")


    end ! subroutine init
!
!______________________________________________________________________
!
    subroutine step( d_in )
!----------------------------------------------------------------------
!  Reads forcing fields during experiment.
!______________________________________________________________________
!
      use glob_const , only: rhoref
      use config     , only: tbias
!      use glob_domain, only: is_master
      use module_time
      use model_run  , only: dti, iint, sec_of_year
      use glob_ocean , only: t

      implicit none

      type(date), intent(in) :: d_in

      logical            ADVANCE_REC, ADVANCE_REC_INT
      integer            max_in_prev, max_in_this, secs
      integer                          &
      , dimension(3)  :: record, year
      real               chunk


! Quit if the module is not used.
      if ( DISABLED ) return

      secs = sec_of_year

      year = d_in%year

      ADVANCE_REC     = .false.
      ADVANCE_REC_INT = .false.

! Decide on the record to read
      if ( USE_CALENDAR ) then

        max_in_this = max_chunks_in_year( d_in%year  , read_int )
        max_in_prev = max_chunks_in_year( d_in%year-1, read_int )

! Decide on the record to read
        chunk     = chunk_of_year( d_in, read_int )
        record(1) = int(chunk)

        if ( secs - int(real(record(1))*read_int) < dti ) then ! TODO: test this one as well.
          ADVANCE_REC = .true.
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
        end if

        record(1) = record(1) + 1

      else

        record(1) = int( iint*dti / read_int ) + 1
        record(2) = record(1)
        record(3) = record(2) + 1

! TODO: right now it interpolates between records (edges of rec1-rec2 span). Implement interpolation between the centers of record spans.
        a = modulo( real(iint*dti), real(read_int) )
        if ( a < dti ) then
          if ( a >= 0. ) then
            ADVANCE_REC = .true.
          end if
        end if
        a = a / read_int
!        print *, "A: ", a

        ADVANCE_REC_INT = ADVANCE_REC

      end if

      if ( interp_ice ) then
        call read_all( ADVANCE_REC_INT, 3, year, record )
        icec = ( 1. - a ) * icec_int(:,:,2) + a * icec_int(:,:,3)
        iceh = ( 1. - a ) * iceh_int(:,:,2) + a * iceh_int(:,:,3)
      else
        call read_all( ADVANCE_REC    , 1, year, record )
      end if

      itsurf = rhoi/rhoref * Ch * ( t(:,:,1) + tbias + 1.82 )
!      print *, minval(itsurf), maxval(itsurf) ,":"
!      stop


    end ! subroutine step

!==============================================================
! Advance ice in time
!--------------------------------------------------------------
!      subroutine ice_advance
!
!        implicit none
!
!        real(kind=rk), dimension(im,jm) ::
!     &                 fx, fy, divu, pice, delx, dely
!     &                ,rhoi, duvi, uidx, vidy, fluxcx, fluxcy
!     &                ,tauiau, tauiav, ht
!        real(kind=rk) eeta, tmp, lhf
!        integer i,j, cnum
!
!        ! FREEZ'T
!        lhf = 323500. ! Latent heat of fusion [J/kg]
!        ht  = (1.-ice)*(wtsurf-swrad)*4.1876d6-ice*(-4.+1.82)*2.13/hi ! Net heat flux [W/m2]
!        do j = 1, jm
!          do i = 1, im
!            if ( ht(i,j) > 0. .and. t(i,j,1) < -1.82 ) then
!              ice(i,j) = ice(i,j) + (1.-ice(i,j))*.47
!              hi(i,j)  = (1.-ice(i,j))*(ht(i,j)*dte/lhf/930.)
!     &                   +   ice(i,j) * hi(i,j)
!            end if
!          end do
!        end do
!
!        eeta = 1.e2 !1.01e-7 ! 1010 cm2/s? ! The source claims the coefficient equals to 10^10 cm2/s! This gives unreallistic Infinities.
!
!        delx = 0.
!        dely = 0.
!        duvi = 0.
!        uidx = 0.
!        vidy = 0.
!        fluxcx = 0.
!        fluxcy = 0.
!
!!        u(:,:,1) = 0.
!!        v(:,:,1) =  .2
!!        el = 0.
!
!! Calculate sea surface elevation gradient
!        delx(2:im,:) = dum(2:im,:)*2.*(el(2:im,:)-el(1:imm1,:))
!     &                     /(dx(2:im,:)+dx(1:imm1,:))
!        dely(:,2:jm) = dvm(:,2:jm)*2.*(el(:,2:jm)-el(:,1:jmm1))
!     &                     /(dy(:,2:jm)+dy(:,1:jmm1))
!        call exchange2d_mpi(delx,im,jm)
!        call exchange2d_mpi(dely,im,jm)
!
!! Get wind stress over ice-free water (convert it from m^2/s^2 to N/m^2)
!        tauiau = -wusurf*rhoref
!        tauiav = -wvsurf*rhoref
!
!! Estimate sea ice density to a unit area
!        rhoi = 900.*hi !*max(ice,.1)
!
!! Calculate water-ice stress
!        duvi=abs(sqrt((ui-u(1:im,1:jm,1))**2+(vi-v(1:im,1:jm,1))**2))
!
!        tauiwu= fsm*5.5e-3*rhoref*(ui-u(1:im,1:jm,1))*duvi
!        tauiwv= fsm*5.5e-3*rhoref*(vi-v(1:im,1:jm,1))*duvi
!
!! Compute ice concentration fluxes
!        do j=1,jm
!          do i=2,imm1
!            if (ui(i,j)<0.) then
!              fluxcx(i,j) = ice(i  ,j)*ui(i,j)*dy(i,j)
!            else
!              fluxcx(i,j) = ice(i-1,j)*ui(i,j)*dy(i,j)
!            end if
!          end do
!        end do
!        do j=2,jmm1
!          do i=1,im
!            if (vi(i,j)<0.) then
!              fluxcy(i,j) = ice(i,j  )*vi(i,j)*dx(i,j)
!            else
!              fluxcy(i,j) = ice(i,j-1)*vi(i,j)*dx(i,j)
!            end if
!          end do
!        end do
!        call exchange2d_mpi(fluxcx,im,jm)
!        call exchange2d_mpi(fluxcy,im,jm)
!
!! Apply ice concentration boundary conditions
!        if (n_north==-1) then
!          do i=1,im
!            if (vi(i,jm)>0.) then
!              fluxcy(i,jm) = ice(i,jm)*vi(i,jm)*dx(i,jm)
!            else
!              fluxcy(i,jm) = cibn(i)  *vi(i,jm)*dx(i,jm)
!            end if
!          end do
!        end if
!        if (n_east==-1) then
!          do j=1,jm
!            if (ui(im,j)>0.) then
!              fluxcx(im,j) = ice(im,j)*ui(im,j)*dy(im,j)
!            else
!              fluxcx(im,j) = cibw(j)  *ui(im,j)*dy(im,j)
!            end if
!          end do
!        end if
!        if (n_south==-1) then
!          do i=1,im
!            if (vi(i,1)<0.) then
!              fluxcy(i,1) = ice(i,1)*vi(i,1)*dx(i,1)
!            else
!              fluxcy(i,1) = cibs(i) *vi(i,1)*dx(i,1)
!            end if
!          end do
!        end if
!        if (n_west==-1) then
!          do j=1,jm
!            if (ui(1,j)<0.) then
!              fluxcx(1,j) = ice(1,j)*ui(1,j)*dy(1,j)
!            else
!              fluxcx(1,j) = cibw(j) *ui(1,j)*dy(1,j)
!            end if
!          end do
!        end if
!
!! Calculate velocity gradients
!        uidx(2:im,:) = (ui(2:im,:)-ui(1:imm1,:))/dx(2:im,:)
!        vidy(:,2:jm) = (vi(:,2:jm)-vi(:,1:jmm1))/dy(:,2:jm)
!        call exchange2d_mpi(uidx,im,jm)
!        call exchange2d_mpi(vidy,im,jm)
!! Apply boundaries
!        if (n_west==-1) uidx(1,:) = 0. !uidx(2,:)
!        if (n_south==-1) vidy(:,1) = 0. !vidy(:,2)
!
!! Derive divergency
!        divu = uidx+vidy
!
!! Get internal stress
!        pice = 0.
!        where (divu<0.) pice = -10.*divu
!
!! Start ice concentration advection
!        if (.false.) then
!          call advtC(ice,icf)
!        else
!          do j=1,jmm1
!            do i=1,imm1
!              tmp = (fluxcx(i,j)-fluxcx(i+1,j)
!     &              +fluxcy(i,j)-fluxcy(i,j+1))
!              icf(i,j) = ice(i,j)+dte*tmp/art(i,j)
!              ! if there will be no ice in the cell, correct output fluxes (velocities) only
!              if (icf(i,j) < 0.) then
!                icf(i,j) = 0.
!                cnum = 0
!                tmp = 0.
!                if (fluxcx(i  ,j)<0.) then
!                  ui(i  ,j) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcx(i,j)
!                end if
!                if (fluxcx(i+1,j)>0.) then
!                  ui(i+1,j) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcx(i+1,j)
!                end if
!                if (fluxcy(i,j  )<0.) then
!                  vi(i,j  ) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcy(i,j)
!                end if
!                if (fluxcy(i,j+1)>0.) then
!                  vi(i,j+1) = 0.
!                  cnum = cnum+1
!                  tmp = tmp+fluxcy(i,j+1)
!                end if
!
!                if (fluxcx(i  ,j)<0.) then
!                  fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
!                end if
!                if (fluxcx(i+1,j)>0.) then
!                  fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
!                end if
!                if (fluxcy(i,j  )<0.) then
!                  fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
!                end if
!                if (fluxcy(i,j+1)>0.) then
!                  fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
!                end if
!              else
!              ! if the cell will be water-free, correct input fluxes (velocities) only
!                if (icf(i,j) > 1.) then
!                  icf(i,j) = 1.
!                  cnum = 0
!                  tmp = 0.
!                  if (fluxcx(i  ,j)>0.) then
!                    ui(i  ,j) = 0.
!                    cnum = cnum+1
!                    tmp = tmp+fluxcx(i,j)
!                  end if
!                  if (fluxcx(i+1,j)<0.) then
!                    ui(i+1,j) = 0.
!                    cnum = cnum+1
!                    tmp = tmp+fluxcx(i+1,j)
!                  end if
!                  if (fluxcy(i,j  )>0.) then
!                    vi(i,j  ) = 0.
!                    cnum = cnum+1
!                    tmp = tmp+fluxcy(i,j)
!                  end if
!                  if (fluxcy(i,j+1)<0.) then
!                    vi(i,j+1) = 0.
!                    cnum = cnum+1
!                    tmp = tmp+fluxcy(i,j+1)
!                  end if
!
!                  if (fluxcx(i  ,j)>0.) then
!                    fluxcx(i  ,j)=fluxcx(i  ,j)-tmp/float(cnum)
!                  end if
!                  if (fluxcx(i+1,j)<0.) then
!                    fluxcx(i+1,j)=fluxcx(i+1,j)+tmp/float(cnum)
!                  end if
!                  if (fluxcy(i,j  )>0.) then
!                    fluxcy(i,j  )=fluxcy(i,j  )-tmp/float(cnum)
!                  end if
!                  if (fluxcy(i,j+1)<0.) then
!                    fluxcy(i,j+1)=fluxcy(i,j+1)+tmp/float(cnum)
!                  end if
!                end if
!              end if
!            end do
!          end do
!          icf = icf*fsm
!        end if
!        call exchange2d_mpi(icf,im,jm)
!
!        fx = 0.
!        fx(2:im,2:jm) = 2.*eeta
!     &                   *((uidx(2:im,2:jm)-uidx(1:imm1,2:jm))
!     &                     /(dx(2:im,2:jm)+dx(1:imm1,2:jm))
!     &                    +(vidy(2:im,2:jm)-vidy(2:im,1:jmm1))
!     &                     /(dy(2:im,2:jm)+dy(2:im,1:jmm1)))
!        fy = fx
!        fx(2:im,2:jm) = fx(2:im,2:jm) +
!     &          eeta*((divu(2:im,2:jm)-divu(1:imm1,2:jm))/dy(2:im,2:jm))
!        fx(2:im,2:jm) = fx(2:im,2:jm) -
!     &                    (pice(2:im,2:jm)-pice(1:imm1,2:jm))
!     &                    /dy(2:im,2:jm)
!
!        fy(2:im,2:jm) = fy(2:im,2:jm) +
!     &          eeta*((divu(2:im,2:jm)-divu(2:im,1:jmm1))/dx(2:im,2:jm))
!        fy(2:im,2:jm) = fy(2:im,2:jm) -
!     &                    (pice(2:im,2:jm)-pice(2:im,1:jmm1))
!     &                    /dx(2:im,2:jm)
!        call exchange2d_mpi(fx,im,jm)
!        call exchange2d_mpi(fy,im,jm)
!
!        do j=1,jmm1
!          do i=1,imm1
!            uif(i,j) = 0.
!            vif(i,j) = 0.
!            if (icf(i,j)>0.) then
!              uif(i,j) = ui(i,j) + dte*
!     &                    ( cor(i,j)*(vi(i,j)+vi(i,j+1))
!     &                     - grav*delx(i,j)
!     &                     + (tauiau(i,j)-tauiwu(i,j))/rhoi(i,j)
!     &                     + fx(i,j) )
!              if (icf(i+1,j)==0.) then
!                uif(i+1,j) = ui(i+1,j) + dte*
!     &                    ( cor(i+1,j)*(vi(i+1,j)+vi(i+1,j+1))
!     &                     - grav*delx(i+1,j)
!     &                     + (tauiau(i+1,j)-tauiwu(i+1,j))/rhoi(i+1,j)
!     &                     + fx(i+1,j) )
!              end if
!              vif(i,j) = vi(i,j) + dte*
!     &                    (-cor(i,j)*(ui(i,j)+ui(i+1,j))
!     &                     - grav*dely(i,j)
!     &                     + (tauiav(i,j)-tauiwv(i,j))/rhoi(i,j)
!     &                     + fy(i,j) )
!              if (icf(i,j+1)==0.) then
!                vif(i,j+1) = vi(i,j+1) + dte*
!     &                    (-cor(i,j+1)*(ui(i,j+1)+ui(i+1,j+1))
!     &                     - grav*dely(i,j+1)
!     &                     + (tauiav(i,j+1)-tauiwv(i,j+1))/rhoi(i,j+1)
!     &                     + fy(i,j+1) )
!              end if
!            end if
!          end do
!        end do
!        call exchange2d_mpi(icf,im,jm)
!        call exchange2d_mpi(uif,im,jm)
!        call exchange2d_mpi(vif,im,jm)
!        if (n_east==-1) then
!          uif(im,:) = uif(imm1,:)
!          vif(im,:) = vif(imm1,:)
!        end if
!        if (n_north==-1) then
!          uif(:,jm) = uif(:,jmm1)
!          vif(:,jm) = vif(:,jmm1)
!        end if
!
!        uif = uif*dum
!        vif = vif*dvm
!
!        uib = ui
!        ui  = uif
!        vib = vi
!        vi  = vif
!        icb = ice
!        ice = icf
!
!      end ! subroutine ice_advance
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
!  Constructs filename string in `<path>YYYY<FORMAT_EXT>` format.
!______________________________________________________________________
!
      implicit none

      character(len=*), intent(in) :: path
      integer         , intent(in) :: year


      if ( path(len(trim(path)):len(trim(path))) == "." ) then
        write( get_filename, '( a, i4.4, a )' ) trim(path)      &
                                               , year            &
                                               , trim(FORMAT_EXT)
      else
        get_filename = path
      end if

    end ! function get_filename
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