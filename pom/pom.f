! pom.f

! main program

      program pom
      
      use module_time
      use river
      use tsforce
      use wind
      use assim
      use interp
      use mcsst        !fhx:mcsst
      use uvforce      !eda:uvforce
      use seaice

      implicit none
      include 'pom.h'

      logical spinup
      namelist/misc_nml/ spinup

      type(date) :: dtime

      open(73, file='switch.nml',status='old')
      read(73, nml=misc_nml)
      close(73)
      
! initialize model
      call initialize

! starting date and time
! read date from restart file name !fhx:
!     dtime = str2date( time_start(1:19) )
      dtime = str2date(read_rst_file(1:13)//":"//
     &read_rst_file(15:16)//":"//read_rst_file(18:19) )

      if ( calc_uvforce ) call uvforce_init(dtime) !eda:uvforce 
      call tsforce_init( dtime )
      if ( calc_tsurf_mc ) call mcsst_init(dtime) !fhx:mcsst
      if ( calc_interp ) call interp_init !fhx:interp_flag
      call wind_init( dtime )
      call river_init( dtime )
      call assim_init( dtime )
      if ( calc_ice ) call ice_init( dtime )
    
      if(my_task == master_task) then
        write(*,'(a)') 'End of initialization'
        write(*,*) 
        write(*,'("d = ",a)')  date2str( dtime )
      endif

      if(nread_rst.eq.0) !call write_output_init_pnetcdf !lyo:20110224:alu:stcc:!lyo:pac10:add this write_output_init*
     &    call write_output_init_pnetcdf( 
     &        "out/init."//trim(netcdf_file)//".nc")


! main loop
      do iint=1,iend

        
!     external forcings
        if (iint==1 .or. .not.spinup) then
          if ( calc_uvforce ) call uvforce_main(dtime)  !eda:uvforce
          call tsforce_main( dtime )
          if ( calc_tsurf_mc )call mcsst_main(dtime)  !fhx:mcsst
          if ( calc_tsforce ) call tsforce_tsflx( dtime )
          if ( calc_wind )    call wind_main( dtime )
          if ( calc_river )   call river_main( dtime, .false. )
          if ( calc_ice )     call ice_main( dtime )
        end if

       
!     advance model
!       call advance( dtime )    !lyo:???
        call advance    
!        call ice_advance

!     drifter data assimilation  !eda:

        if ( calc_assimdrf ) call assimdrf_main( dtime )

!     satellite SSHA data assimilation
        if ( calc_assim ) call assim_main( dtime )

       
        dtime = dtime + int( dti )

        if ( my_task.eq.master_task ) 
     &       write(*,'("d = ",a)')  date2str( dtime )

! write output
      call write_output( dtime )

! write SURF output
      call write_output_surf !fhx:20110131:

! write restart
      call write_restart( dtime )

!      call write_debug_pnetcdf("dbg."//date2str(dtime))

      end do


! finalize mpi
      call finalize_mpi

      stop
      end

!_______________________________________________________________________
      real(kind=8) function realtime()
      call system_clock(i,j,k)
      realtime=i/dble(j)
      return
      end


!_______________________________________________________________________
      subroutine write_output( d_in )

      use module_time
      use assim, only : assim_store_ssha

      implicit none
      include 'pom.h'

      type(date), intent(in) :: d_in

      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprint).eq.0) then

         
         uab_mean    = uab_mean    / real ( num )
         vab_mean    = vab_mean    / real ( num )
         elb_mean    = elb_mean    / real ( num )
         wusurf_mean = wusurf_mean / real ( num )
         wvsurf_mean = wvsurf_mean / real ( num )
         wtsurf_mean = wtsurf_mean / real ( num )
         wssurf_mean = wssurf_mean / real ( num )
         u_mean      = u_mean      / real ( num )
         v_mean      = v_mean      / real ( num )
         w_mean      = w_mean      / real ( num )
         t_mean      = t_mean      / real ( num )
         s_mean      = s_mean      / real ( num )
         rho_mean    = rho_mean    / real ( num )
         kh_mean     = kh_mean     / real ( num )
         km_mean     = km_mean     / real ( num )



!     store ssha for data assimilation to t(k=kb)

         if ( calc_assim ) 
     &        call assim_store_ssha( t_mean(:,:,kb), 't(k=kb)', d_in )

!     store tsurf & ssurf in rho(k=kb) & s(k=kb) respectively
!     note that these will be satellite &/or monthly climatology

      s_mean(:,:,kb)   = ssurf(:,:) !lyo:exp301!lyo:exp302:store ssurf!lyonew:
      rho_mean(:,:,kb) = tsurf(:,:) !lyo:exp301!lyo:exp302:store tsurf!lyonew:

!     fill up ghost cells before output
         call exchange2d_mpi( uab_mean, im, jm )
         call exchange2d_mpi( vab_mean, im, jm )
         call exchange2d_mpi( elb_mean, im, jm )
         call exchange2d_mpi( wusurf_mean, im, jm )
         call exchange2d_mpi( wvsurf_mean, im, jm )
         call exchange2d_mpi( wtsurf_mean, im, jm )
         call exchange2d_mpi( wssurf_mean, im, jm )
         call exchange3d_mpi( u_mean, im, jm, kb )
         call exchange3d_mpi( v_mean, im, jm, kb )
         call exchange3d_mpi( w_mean, im, jm, kb )
         call exchange3d_mpi( t_mean, im, jm, kb )
         call exchange3d_mpi( s_mean, im, jm, kb )
         call exchange3d_mpi( rho_mean, im, jm, kb )
         call exchange3d_mpi( kh_mean, im, jm, kb )
         call exchange3d_mpi( km_mean, im, jm, kb )





!         if ( my_task == 41 ) 
!     $        print*, im/2,jm/2,rot(im/2,jm/2),
!     $        uab_mean(im/2,jm/2),vab_mean(im/2,jm/2)
!
!     u,v winds ---> zonal and meridional winds
!         do j = 1, jm
!            do i = 1, im
!               u_tmp = uab_mean(i,j)
!               v_tmp = vab_mean(i,j)
!               uab_mean(i,j) 
!     $              = u_tmp * cos( rot(i,j) * deg2rad )
!     $              - v_tmp * sin( rot(i,j) * deg2rad )
!               vab_mean(i,j) 
!     $              = u_tmp * sin( rot(i,j) * deg2rad )
!     $              + v_tmp * cos( rot(i,j) * deg2rad )
!            enddo
!         enddo
!
!         if ( my_task == 41 ) 
!     $        print*, im/2,jm/2,
!     $        cos(rot(im/2,jm/2)*deg2rad),
!     $        uab_mean(im/2,jm/2),vab_mean(im/2,jm/2)
!
!         do j = 1, jm
!            do i = 1, im
!               u_tmp = wusurf_mean(i,j)
!               v_tmp = wvsurf_mean(i,j)
!               wusurf_mean(i,j) 
!     $              = u_tmp * cos( rot(i,j) * deg2rad )
!     $              - v_tmp * sin( rot(i,j) * deg2rad )
!               wvsurf_mean(i,j) 
!     $              = u_tmp * sin( rot(i,j) * deg2rad )
!     $              + v_tmp * cos( rot(i,j) * deg2rad )
!            enddo
!         enddo
!
!         do k=1,kbm1
!            do j = 1, jm
!               do i = 1, im
!                  u_tmp = u_mean(i,j,k)
!                  v_tmp = v_mean(i,j,k)
!                  u_mean(i,j,k) 
!     $                 = u_tmp * cos( rot(i,j) * deg2rad )
!     $                 - v_tmp * sin( rot(i,j) * deg2rad )
!                  v_mean(i,j,k) 
!     $                 = u_tmp * sin( rot(i,j) * deg2rad )
!     $                 + v_tmp * cos( rot(i,j) * deg2rad )
!               enddo
!            enddo
!         enddo
         
       if (output_flag == 1) then

         call write_output_pnetcdf( 
     $        "out/"//trim(netcdf_file)//".nc")

       else if (output_flag == 0) then

         call write_output_pnetcdf0( 
     $        "out/"//trim(netcdf_file)//"."//
     $        date2str(d_in)//".nc" )

       end if

         uab_mean    = 0.0
         vab_mean    = 0.0
         elb_mean    = 0.0
         wusurf_mean = 0.0
         wvsurf_mean = 0.0
         wtsurf_mean = 0.0
         wssurf_mean = 0.0
         u_mean      = 0.0
         v_mean      = 0.0
         w_mean      = 0.0
         t_mean      = 0.0
         s_mean      = 0.0
         rho_mean    = 0.0
         kh_mean     = 0.0
         km_mean     = 0.0
         
         num = 0

      endif

      return
      end
!-------------------------------------------------------

!
!_______________________________________________________________________
      subroutine write_output_surf !fhx:20110131: new subr.
      use module_time

      implicit none
      include 'pom.h'

      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprints).eq.0) then

         
         usrf_mean    = usrf_mean    / real ( nums )
         vsrf_mean    = vsrf_mean    / real ( nums )
         elsrf_mean   = elsrf_mean   / real ( nums )
         uwsrf_mean = uwsrf_mean / real ( nums )
         vwsrf_mean = vwsrf_mean / real ( nums )


!     fill up ghost cells before output
         call exchange2d_mpi( usrf_mean, im, jm )
         call exchange2d_mpi( vsrf_mean, im, jm )
         call exchange2d_mpi( elsrf_mean, im, jm )
         call exchange2d_mpi( uwsrf_mean, im, jm )
         call exchange2d_mpi( vwsrf_mean, im, jm )
       
       if (SURF_flag==1)
     $    call write_SURF_pnetcdf( 
     $        "out/SRF."//trim(netcdf_file)//".nc")
       
         usrf_mean    = 0.0
         vsrf_mean    = 0.0
         elsrf_mean   = 0.0
         uwsrf_mean = 0.0
         vwsrf_mean = 0.0
         
         nums = 0

      endif

      return
      end
!-------------------------------------------------------

!_______________________________________________________________________
      subroutine write_restart( d_in )

      use module_time

      implicit none
      include 'pom.h'

      type(date), intent(in) :: d_in

!lyo:exp302:The following changes will create a restart 1day after the 
!     initial date_start0 to be used for next-day's ncast&fcast
      real(kind=rk), parameter :: write_rst1d = 1.0 !lyo:exp302:daily analysis restart 
      integer irestart1d, ncount           !           each run for next-day 
      data ncount/0/                       !           fcast

      irestart1d=nint(write_rst1d*86400./dti)                    !lyo:exp302:
      if(     (mod(iint,irestart)   == 0                  )        !lyo:exp302:
     $   .or. (mod(iint,irestart1d) == 0 .and. ncount.eq.0) ) then !lyo:exp302:
         ncount=1                                                  !lyo:exp302:
!     if( mod(iint,irestart) == 0 ) then                           !lyo:exp302:

! fill up ghost cells 
! They have to be filled in creating a restart file,
! which is independent of the number of nodes.

        call exchange2d_mpi(advua,im,jm)
        call exchange2d_mpi(advva,im,jm)
        call exchange2d_mpi(adx2d,im,jm)
        call exchange2d_mpi(ady2d,im,jm)
        call exchange2d_mpi(utb,im,jm)
        call exchange2d_mpi(vtb,im,jm)
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
        call exchange2d_mpi(wubot,im,jm)
        call exchange2d_mpi(wvbot,im,jm)
        call exchange2d_mpi(aam2d,im,jm)
        call exchange2d_mpi(ua,im,jm)
        call exchange2d_mpi(uab,im,jm)
        call exchange2d_mpi(va,im,jm)
        call exchange2d_mpi(vab,im,jm)
        call exchange2d_mpi(el,im,jm)
        call exchange2d_mpi(elb,im,jm)
        call exchange2d_mpi(et,im,jm)
        call exchange2d_mpi(etb,im,jm)
        call exchange2d_mpi(egb,im,jm)

        call exchange3d_mpi(u,im,jm,kb)
        call exchange3d_mpi(v,im,jm,kb)
        call exchange3d_mpi(ub,im,jm,kb)
        call exchange3d_mpi(vb,im,jm,kb)
        call exchange3d_mpi(w,im,jm,kb)
        call exchange3d_mpi(t,im,jm,kb)
        call exchange3d_mpi(tb,im,jm,kb)
        call exchange3d_mpi(s,im,jm,kb)
        call exchange3d_mpi(sb,im,jm,kb)
        call exchange3d_mpi(rho,im,jm,kb)
        call exchange3d_mpi(km,im,jm,kb)
        call exchange3d_mpi(kh,im,jm,kb)
        call exchange3d_mpi(kq,im,jm,kb)
        call exchange3d_mpi(l,im,jm,kb)
        call exchange3d_mpi(q2,im,jm,kb)
        call exchange3d_mpi(q2b,im,jm,kb)
        call exchange3d_mpi(aam,im,jm,kb)
        call exchange3d_mpi(q2l,im,jm,kb)
        call exchange3d_mpi(q2lb,im,jm,kb)


         call write_restart_pnetcdf( 
     $        "out/"//date2str(d_in)//"."//
     $        trim(netcdf_file)//".rst" )
         
      endif


      return
      end
!-------------------------------------------------------
