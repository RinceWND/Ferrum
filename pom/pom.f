! pom.f

! main program

      program pom

      use config
      use glob_domain, only: is_master
      use io
      use model_run  , only: dti, dtime, iend, iint

      use module_time
      use seaice

      implicit none

      real(8), external :: realtime

      integer mb
      real(8) tick0, tick1


! initialize model
      call initialize

      mb = dtime%month

      tick0 = realtime()
      if ( is_master ) then
        call msg_print("BEGIN NUMERICAL EXPERIMENT", 1, "")
        print '(" = ",a)', date2str( dtime )
      end if

! main loop
      do iint = 1,iend

!     advance model
        call advance
!        call ice_advance

!     drifter data assimilation  !eda:
!        if ( calc_assimdrf ) call assimdrf_main( dtime )

!     satellite SSHA data assimilation
!        if ( calc_assim ) call assim_main( dtime )

!     increment time
        dtime = dtime + int( dti )

        if ( is_master ) print '(" = ",a)', date2str( dtime )

! write output
        call output_results( dtime, monthly_flag, mb )

! write SURF output
!        call write_output_surf !fhx:20110131:

! write restart
        call write_restart( dtime )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TODO: CHECK PARALLEL RUN WITH DIFFERENT PROCs COUNT!!!
!       Uninitialized values too!                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call out_debug("ber."//date2str(dtime))
!      stop

      end do

      tick1 = realtime()
      call msg_print("NUMERICAL EXPERIMENT FINISHED", 1, "")

! finalize mpi
      call finalize_mpi

      if (is_master) print '(a,f7.2,a)',"Done in ",(tick1-tick0)/60.
     &                                 ," min."


      end program

!______________________________________________________________________
      real(kind=8) function realtime()
      call system_clock(i,j,k)
      realtime=dble(i)/dble(j)
      return
      end

!______________________________________________________________________
      subroutine output_results( d_in, output_mode, mb )

      use config     , only: calc_assim , netcdf_file
     &                     , output_flag, output_means
      use glob_domain, only: im, jm, km
      use glob_ocean , only: ssurf, tsurf
      use glob_out
      use model_run  , only: iint
      use module_time
!      use assim, only : assim_store_ssha

      implicit none

      type(date), intent(in) :: d_in
      integer, intent(in) :: output_mode
      integer, intent(inout) :: mb

      if ( netcdf_file /= 'nonetcdf' ) then

        if(( output_mode==0 .and. mod(iint,iprint)==0 ) .or.
     &     ( output_mode==1 .and. mb /= d_in%month )) then

          mb = d_in%month

          uab_mean    = uab_mean    / real ( num )
          vab_mean    = vab_mean    / real ( num )
          elb_mean    = elb_mean    / real ( num )
          wusurf_mean = wusurf_mean / real ( num )
          wvsurf_mean = wvsurf_mean / real ( num )
          wtsurf_mean = wtsurf_mean / real ( num )
          wssurf_mean = wssurf_mean / real ( num )
          swrad_mean  = swrad_mean  / real ( num )
          u_mean      = u_mean      / real ( num )
          v_mean      = v_mean      / real ( num )
          w_mean      = w_mean      / real ( num )
          t_mean      = t_mean      / real ( num )
          s_mean      = s_mean      / real ( num )
          rho_mean    = rho_mean    / real ( num )
          kh_mean     = kh_mean     / real ( num )
          km_mean     = km_mean     / real ( num )
          aam_mean    = aam_mean    / real ( num )



!     store ssha for data assimilation to t(k=kb)
!          if ( calc_assim )
!     &        call assim_store_ssha( t_mean(:,:,kb), 't(k=kb)', d_in )

!     store tsurf & ssurf in rho(k=kb) & s(k=kb) respectively
!     note that these will be satellite &/or monthly climatology

          s_mean(:,:,km)   = ssurf(:,:) !lyo:exp301!lyo:exp302:store ssurf!lyonew:
          rho_mean(:,:,km) = tsurf(:,:) !lyo:exp301!lyo:exp302:store tsurf!lyonew:

!     fill up ghost cells before output
          call exchange2d_mpi( uab_mean, im, jm )
          call exchange2d_mpi( vab_mean, im, jm )
          call exchange2d_mpi( elb_mean, im, jm )
          call exchange2d_mpi( wusurf_mean, im, jm )
          call exchange2d_mpi( wvsurf_mean, im, jm )
          call exchange2d_mpi( wtsurf_mean, im, jm )
          call exchange2d_mpi( wssurf_mean, im, jm )
          call exchange2d_mpi( swrad_mean, im, jm )
          call exchange3d_mpi( u_mean, im, jm, km )
          call exchange3d_mpi( v_mean, im, jm, km )
          call exchange3d_mpi( w_mean, im, jm, km )
          call exchange3d_mpi( t_mean, im, jm, km )
          call exchange3d_mpi( s_mean, im, jm, km )
          call exchange3d_mpi( rho_mean, im, jm, km )
          call exchange3d_mpi( kh_mean, im, jm, km )
          call exchange3d_mpi( km_mean, im, jm, km )
          call exchange3d_mpi( aam_mean, im, jm, km )


          if ( output_flag == 1 ) then

            call write_output(
     &        "out/"//trim(netcdf_file)//".nc"
     &       , out_record, output_means )

          else if (output_flag == 0) then

            call create_output(
     &        "out/"//trim(netcdf_file)//"."//
     &        date2str(d_in)//".nc" )
            call write_output(
     &        "out/"//trim(netcdf_file)//"."//
     &        date2str(d_in)//".nc", 1, output_means )

          end if

          uab_mean    = 0.
          vab_mean    = 0.
          elb_mean    = 0.
          wusurf_mean = 0.
          wvsurf_mean = 0.
          wtsurf_mean = 0.
          wssurf_mean = 0.
          swrad_mean  = 0.
          u_mean      = 0.
          v_mean      = 0.
          w_mean      = 0.
          t_mean      = 0.
          s_mean      = 0.
          rho_mean    = 0.
          kh_mean     = 0.
          km_mean     = 0.
          aam_mean    = 0.

          num = 0

          out_record = out_record + 1

        end if
      end if

      return
      end !subroutine output_results
!______________________________________________________________________
!
      subroutine write_output_surf !fhx:20110131: new subr.
!----------------------------------------------------------------------
!  Prepares output vars and calls output procedure.
!______________________________________________________________________

      use config     , only: netcdf_file, SURF_flag
      use glob_domain, only: im, jm
      use glob_out
      use model_run  , only: iint
      use module_time

      implicit none

      if (      netcdf_file /= 'nonetcdf'
     &    .and. mod(iint,iprints) == 0   ) then

        usrf_mean  = usrf_mean  / real ( nums )
        vsrf_mean  = vsrf_mean  / real ( nums )
        elsrf_mean = elsrf_mean / real ( nums )
        uwsrf_mean = uwsrf_mean / real ( nums )
        vwsrf_mean = vwsrf_mean / real ( nums )

!     fill up ghost cells before output
        call exchange2d_mpi( usrf_mean, im, jm )
        call exchange2d_mpi( vsrf_mean, im, jm )
        call exchange2d_mpi( elsrf_mean, im, jm )
        call exchange2d_mpi( uwsrf_mean, im, jm )
        call exchange2d_mpi( vwsrf_mean, im, jm )

!        if ( SURF_flag == 1 )
!     $    call write_SURF_pnetcdf(
!     $        "out/SRF."//trim(netcdf_file)//".nc")

        usrf_mean  = 0.0
        vsrf_mean  = 0.0
        elsrf_mean = 0.0
        uwsrf_mean = 0.0
        vwsrf_mean = 0.0

        nums = 0

      end if

      return
      end
!______________________________________________________________________
!
      subroutine write_restart( d_in )
!----------------------------------------------------------------------
!  Exchanges all output variables and writes restart file.
!______________________________________________________________________

      use air        , only: wusurf, wvsurf
      use config     , only: netcdf_file
      use glob_const , only: rk
      use glob_domain, only: im, jm, km
      use glob_ocean , only: aam, aam2d, advua, advva, adx2d, ady2d
     &                     , egb, el, elb, et, etb, kh, kmt, kq, l
     &                     , q2, q2b, q2l, q2lb, rho, s, sb, t, tb
     &                     , u, ua, uab, ub, utb, v, va, vab, vb, vtb
     &                     , w, wubot, wvbot
      use glob_out   , only: irestart
      use model_run  , only: dti, iend, iint
      use module_time

      implicit none

      type(date), intent(in) :: d_in

!lyo:exp302:The following changes will create a restart 1day after the
!     initial date_start0 to be used for next-day's ncast&fcast
      real(rk), parameter :: write_rst1d = 1.0
      integer irestart1d, ncount ! each run for next-day
      data ncount/0/             ! fcast

      irestart1d = nint(write_rst1d*86400./dti)
      if (     (mod(iint,irestart)   == 0                  )
     &    .or. (mod(iint,irestart1d) == 0 .and. ncount == 0)
     &    .or. ( iint == iend ) ) then
        ncount = 1

! fill up ghost cells
! They have to be filled in creating a restart file,
! which is independent of the number of nodes.

        call exchange2d_mpi(advua,im,jm)
        call exchange2d_mpi(advva,im,jm)
        call exchange2d_mpi(adx2d,im,jm)
        call exchange2d_mpi(ady2d,im,jm)
        call exchange2d_mpi(utb  ,im,jm)
        call exchange2d_mpi(vtb  ,im,jm)
!        call exchange2d_mpi(q2(:,:,1),im,jm)
!        call exchange2d_mpi(q2(:,:,kb),im,jm)
!        call exchange2d_mpi(q2b(:,:,1),im,jm)
!        call exchange2d_mpi(q2b(:,:,kb),im,jm)
!        call exchange2d_mpi(q2l(:,:,1),im,jm)
!        call exchange2d_mpi(q2l(:,:,kb),im,jm)
!        call exchange2d_mpi(q2lb(:,:,1),im,jm)
!        call exchange2d_mpi(q2lb(:,:,kb),im,jm)
!        call exchange2d_mpi(kq(:,:,1),im,jm)
!        call exchange2d_mpi(kq(:,:,kb),im,jm)
!        call exchange2d_mpi(t(:,:,1),im,jm)
!        call exchange2d_mpi(t(:,:,kb),im,jm)
!        call exchange2d_mpi(tb(:,:,1),im,jm)
!        call exchange2d_mpi(tb(:,:,kb),im,jm)
!        call exchange3d_mpi(ub,im,jm,kb)
!        call exchange3d_mpi(vb,im,jm,kb)

        call exchange2d_mpi(wusurf,im,jm)
        call exchange2d_mpi(wvsurf,im,jm)
        call exchange2d_mpi(wubot ,im,jm)
        call exchange2d_mpi(wvbot ,im,jm)
        call exchange2d_mpi(aam2d ,im,jm)
        call exchange2d_mpi(ua    ,im,jm)
        call exchange2d_mpi(uab   ,im,jm)
        call exchange2d_mpi(va    ,im,jm)
        call exchange2d_mpi(vab   ,im,jm)
        call exchange2d_mpi(el    ,im,jm)
        call exchange2d_mpi(elb   ,im,jm)
        call exchange2d_mpi(et    ,im,jm)
        call exchange2d_mpi(etb   ,im,jm)
        call exchange2d_mpi(egb   ,im,jm)

        call exchange3d_mpi(u   ,im,jm,km)
        call exchange3d_mpi(v   ,im,jm,km)
        call exchange3d_mpi(ub  ,im,jm,km)
        call exchange3d_mpi(vb  ,im,jm,km)
        call exchange3d_mpi(w   ,im,jm,km)
        call exchange3d_mpi(t   ,im,jm,km)
        call exchange3d_mpi(tb  ,im,jm,km)
        call exchange3d_mpi(s   ,im,jm,km)
        call exchange3d_mpi(sb  ,im,jm,km)
        call exchange3d_mpi(rho ,im,jm,km)
        call exchange3d_mpi(kmt ,im,jm,km)
        call exchange3d_mpi(kh  ,im,jm,km)
        call exchange3d_mpi(kq  ,im,jm,km)
        call exchange3d_mpi(l   ,im,jm,km)
        call exchange3d_mpi(q2  ,im,jm,km)
        call exchange3d_mpi(q2b ,im,jm,km)
        call exchange3d_mpi(aam ,im,jm,km)
        call exchange3d_mpi(q2l ,im,jm,km)
        call exchange3d_mpi(q2lb,im,jm,km)


        call create_restart(
     $        "out/"//date2str(d_in)//"."//
     $        trim(netcdf_file)//".rst" )

      end if


      return
      end
!______________________________________________________________________
!
      subroutine msg_print( msg, status, desc )
!----------------------------------------------------------------------
!  Prints a message. TODO: move to a separate module since other
! modules depend on this routine.
!______________________________________________________________________

        use glob_domain, only: is_master

        implicit none

        character(*), intent(in) :: msg, desc
        integer     , intent(in) :: status

        integer                 , parameter :: line_len = 56
        integer(2), dimension(6), parameter :: chr =
!     &                          (/ 9675  ! 1: white circle
     &                          (/ 79_2    ! 1: latin capital letter o
     &                           , 33_2    ! 2: exclamation mark
!     &                           , 10799 ! 3: vector or cross product "X"
     &                           , 88_2    ! 3: latin capital letter x
     &                           , 63_2    ! 4: question mark
!     &                           , 8270  ! 5: low asterisk
     &                           , 42_2    ! 5: asterisk
     &                           , 32_2    ! 6: space
     &                          /)
        integer               i
        character(64)         fmt01!, fmt02
        character(line_len-4) line


! If fatal error prepare to terminate execution
        if ( status == 3 ) then
          call finalize_mpi
! If not fatal leave just master to print the message
        elseif ( .not.is_master ) then
          return
        end if

        line = ""

! Print simple message if only `DESC` is present.
        if ( msg == '' .and. desc /= '' ) then
          print '("[",a,"] ",a)', char(chr(status)), desc
          if ( status == 3 ) stop
          return
        end if

        write( fmt01, '("(",i2,"a)")') line_len
        i = line_len-3 - len(msg)
        line(max(1,i):) = msg

        print '(/)'
        print fmt01, ("_",i=1,line_len)

        print '("[",a,"] ",a/)', char(chr(status)), line

        if ( desc == "" ) return

        print fmt01, ("-",i=1,line_len)
        print fmt01, (desc(i:i),i=1,len(desc))
        print fmt01, ("-",i=1,line_len)

        if ( status == 3 ) stop


      end subroutine
