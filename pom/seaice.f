      module seaice

      implicit none

      private

      public :: ice_init, ice_main, ice_advance

      include 'pom.h'


      real(kind=rk), dimension( im_local, jm_local ) ::
     &     ice_a, ice_b

    
      integer :: day_a, day_b
      character(len=16) :: infile_a, infile_b
      real(kind=rk) :: aa


      contains

!==============================================================
! Initialization variables for U & V boundary condition 
!--------------------------------------------------------------
      subroutine ice_init( d_in )

        use module_time

        implicit none

        type(date), intent(in) :: d_in 
        type(date) d_2
        real(kind=rk) secs


        d_2 = str2date("1979-01-01 00:00:00") ! TODO: Fix ice date generation
        d_2%year = d_in%year
        secs = d_in-d_2

        day_b = 2*(int(secs/172800.)+1)
        day_a = day_b-2

        d_2 = d_2 + day_a*86400

        aa = (day_b-day_a)/(day_b+day_a)


!     data open.
        write( infile_a, '( "ice.",i4,2(i2.2),".nc" )' )
     &                                      d_2%year, d_2%month, d_2%day  
        call read_ice_pnetcdf( ice_a, infile_a )

        d_2 = d_2 + 2*86400

        write( infile_b, '( "ice.",i4,2(i2.2),".nc" )' )
     &                                      d_2%year, d_2%month, d_2%day  
        call read_ice_pnetcdf( ice_b, infile_b )
     
        if ( my_task == master_task ) 
     $        write(*,'(/a/)') "---------- ice_init."

        ice = ( 1.0 - aa ) * ice_a + aa * ice_b

        return

      end subroutine ice_init
!--------------------------------------------------------------

!==============================================================
! Read and interpolate in time
!--------------------------------------------------------------
      subroutine ice_main( d_in )

        use module_time

        implicit none

        type(date), intent(in) :: d_in
        type(date) d_2
        real(kind=rk) secs
        real(kind=rk), dimension(im_local,jm_local) :: ci


        d_2 = str2date("1979-01-01 00:00:00")
        d_2%year = d_in%year
        secs = d_in-d_2

        day_b = 2*(int(secs/172800.)+1)

        d_2 = d_2 + day_b*86400

        aa = (day_b-day_a)/(day_b+day_a)

!     data open.

        if ( day_b-day_a > 4 ) then

          day_a = day_b-2

          ice_a = ice_b
          aa = 0.

          write( infile_b, '( "ice.",i4,2(i2.2),".nc" )' )
     &                                      d_2%year, d_2%month, d_2%day  
          call read_ice_pnetcdf( ice_b, infile_b )

        endif
     

!     time interpolation.

        ci = ( 1.0 - aa ) * ice_a + aa * ice_b

        cibn = ci(:,jm)
        cibe = ci(im,:)
        cibs = ci(:, 1)
        cibw = ci( 1,:)

      return

      end subroutine ice_main
!--------------------------------------------------------------
      
!==============================================================
! Advance ice in time
!--------------------------------------------------------------
      subroutine ice_advance

        implicit none

        real(kind=rk), dimension(im,jm) ::
     &                 fx, fy, divu, pice, delx, dely
     &                ,rhoi, duvi, uidx, vidy, fluxcx, fluxcy
     &                ,tauiau, tauiav
        real(kind=rk) eeta, tmp
        integer i,j, cnum

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
        delx(2:im,:) = dum(2:im,:)*2.*(el(2:im,:)-el(1:imm1,:))
     &                     /(dx(2:im,:)+dx(1:imm1,:))
        dely(:,2:jm) = dvm(:,2:jm)*2.*(el(:,2:jm)-el(:,1:jmm1))
     &                     /(dy(:,2:jm)+dy(:,1:jmm1))
        call exchange2d_mpi(delx,im,jm)
        call exchange2d_mpi(dely,im,jm)

! Get wind stress over ice-free water (convert it from m^2/s^2 to N/m^2)
        tauiau = -wusurf*rhoref
        tauiav = -wvsurf*rhoref

! Estimate sea ice density to a unit area
        rhoi = 900.*hi !*max(ice,.1)

! Calculate water-ice stress
        duvi=abs(sqrt((ui-u(1:im,1:jm,1))**2+(vi-v(1:im,1:jm,1))**2))

        tauiwu= fsm*5.5e-3*rhoref*(ui-u(1:im,1:jm,1))*duvi
        tauiwv= fsm*5.5e-3*rhoref*(vi-v(1:im,1:jm,1))*duvi

! Compute ice concentration fluxes
        do j=1,jm
          do i=2,imm1
            if (ui(i,j)<0.) then
              fluxcx(i,j) = ice(i  ,j)*ui(i,j)*dy(i,j)
            else
              fluxcx(i,j) = ice(i-1,j)*ui(i,j)*dy(i,j)
            end if
          end do
        end do
        do j=2,jmm1
          do i=1,im
            if (vi(i,j)<0.) then
              fluxcy(i,j) = ice(i,j  )*vi(i,j)*dx(i,j)
            else
              fluxcy(i,j) = ice(i,j-1)*vi(i,j)*dx(i,j)
            end if
          end do
        end do
        call exchange2d_mpi(fluxcx,im,jm)
        call exchange2d_mpi(fluxcy,im,jm)

! Apply ice concentration boundary conditions
        if (n_north==-1) then
          do i=1,im
            if (vi(i,jm)>0.) then
              fluxcy(i,jm) = ice(i,jm)*vi(i,jm)*dx(i,jm)
            else
              fluxcy(i,jm) = cibn(i)  *vi(i,jm)*dx(i,jm)
            end if
          end do
        end if
        if (n_east==-1) then
          do j=1,jm
            if (ui(im,j)>0.) then
              fluxcx(im,j) = ice(im,j)*ui(im,j)*dy(im,j)
            else
              fluxcx(im,j) = cibw(j)  *ui(im,j)*dy(im,j)
            end if
          end do
        end if
        if (n_south==-1) then
          do i=1,im
            if (vi(i,1)<0.) then
              fluxcy(i,1) = ice(i,1)*vi(i,1)*dx(i,1)
            else
              fluxcy(i,1) = cibs(i) *vi(i,1)*dx(i,1)
            end if
          end do
        end if
        if (n_west==-1) then
          do j=1,jm
            if (ui(1,j)<0.) then
              fluxcx(1,j) = ice(1,j)*ui(1,j)*dy(1,j)
            else
              fluxcx(1,j) = cibw(j) *ui(1,j)*dy(1,j)
            end if
          end do
        end if

! Calculate velocity gradients
        uidx(2:im,:) = (ui(2:im,:)-ui(1:imm1,:))/dx(2:im,:)
        vidy(:,2:jm) = (vi(:,2:jm)-vi(:,1:jmm1))/dy(:,2:jm)
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
          call advtC(ice,icf)
        else
          do j=1,jmm1
            do i=1,imm1
              tmp = (fluxcx(i,j)-fluxcx(i+1,j)
     &              +fluxcy(i,j)-fluxcy(i,j+1))
              icf(i,j) = ice(i,j)+dte*tmp/art(i,j)
              ! if there will be no ice in the cell, correct output fluxes (velocities) only
              if (icf(i,j) < 0.) then
                icf(i,j) = 0.
                cnum = 0
                tmp = 0.
                if (fluxcx(i  ,j)<0.) then
                  ui(i  ,j) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcx(i,j)
                end if
                if (fluxcx(i+1,j)>0.) then
                  ui(i+1,j) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcx(i+1,j)
                end if
                if (fluxcy(i,j  )<0.) then
                  vi(i,j  ) = 0.
                  cnum = cnum+1
                  tmp = tmp+fluxcy(i,j)
                end if
                if (fluxcy(i,j+1)>0.) then
                  vi(i,j+1) = 0.
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
                if (icf(i,j) > 1.) then
                  icf(i,j) = 1.
                  cnum = 0
                  tmp = 0.
                  if (fluxcx(i  ,j)>0.) then
                    ui(i  ,j) = 0.
                    cnum = cnum+1
                    tmp = tmp+fluxcx(i,j)
                  end if
                  if (fluxcx(i+1,j)<0.) then
                    ui(i+1,j) = 0.
                    cnum = cnum+1
                    tmp = tmp+fluxcx(i+1,j)
                  end if
                  if (fluxcy(i,j  )>0.) then
                    vi(i,j  ) = 0.
                    cnum = cnum+1
                    tmp = tmp+fluxcy(i,j)
                  end if
                  if (fluxcy(i,j+1)<0.) then
                    vi(i,j+1) = 0.
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
          icf = icf*fsm
        end if
        call exchange2d_mpi(icf,im,jm)

        fx = 0.
        fx(2:im,2:jm) = 2.*eeta
     &                   *((uidx(2:im,2:jm)-uidx(1:imm1,2:jm))
     &                     /(dx(2:im,2:jm)+dx(1:imm1,2:jm))
     &                    +(vidy(2:im,2:jm)-vidy(2:im,1:jmm1))
     &                     /(dy(2:im,2:jm)+dy(2:im,1:jmm1)))
        fy = fx
        fx(2:im,2:jm) = fx(2:im,2:jm) +
     &          eeta*((divu(2:im,2:jm)-divu(1:imm1,2:jm))/dy(2:im,2:jm))
        fx(2:im,2:jm) = fx(2:im,2:jm) -
     &                    (pice(2:im,2:jm)-pice(1:imm1,2:jm))
     &                    /dy(2:im,2:jm)

        fy(2:im,2:jm) = fy(2:im,2:jm) +
     &          eeta*((divu(2:im,2:jm)-divu(2:im,1:jmm1))/dx(2:im,2:jm))
        fy(2:im,2:jm) = fy(2:im,2:jm) -
     &                    (pice(2:im,2:jm)-pice(2:im,1:jmm1))
     &                    /dx(2:im,2:jm)
        call exchange2d_mpi(fx,im,jm)
        call exchange2d_mpi(fy,im,jm)

        do j=1,jmm1
          do i=1,imm1
            uif(i,j) = 0.
            vif(i,j) = 0.
            if (icf(i,j)>0.) then
              uif(i,j) = ui(i,j) + dte*
     &                    ( cor(i,j)*(vi(i,j)+vi(i,j+1))
     &                     - grav*delx(i,j)
     &                     + (tauiau(i,j)-tauiwu(i,j))/rhoi(i,j)
     &                     + fx(i,j) )
              if (icf(i+1,j)==0.) then
                uif(i+1,j) = ui(i+1,j) + dte*
     &                    ( cor(i+1,j)*(vi(i+1,j)+vi(i+1,j+1))
     &                     - grav*delx(i+1,j)
     &                     + (tauiau(i+1,j)-tauiwu(i+1,j))/rhoi(i+1,j)
     &                     + fx(i+1,j) )
              end if
              vif(i,j) = vi(i,j) + dte*
     &                    (-cor(i,j)*(ui(i,j)+ui(i+1,j))
     &                     - grav*dely(i,j)
     &                     + (tauiav(i,j)-tauiwv(i,j))/rhoi(i,j)
     &                     + fy(i,j) )
              if (icf(i,j+1)==0.) then
                vif(i,j+1) = vi(i,j+1) + dte*
     &                    (-cor(i,j+1)*(ui(i,j+1)+ui(i+1,j+1))
     &                     - grav*dely(i,j+1)
     &                     + (tauiav(i,j+1)-tauiwv(i,j+1))/rhoi(i,j+1)
     &                     + fy(i,j+1) )
              end if
            end if
          end do
        end do
        call exchange2d_mpi(icf,im,jm)
        call exchange2d_mpi(uif,im,jm)
        call exchange2d_mpi(vif,im,jm)
        if (n_east==-1) then
          uif(im,:) = uif(imm1,:)
          vif(im,:) = vif(imm1,:)
        end if
        if (n_north==-1) then
          uif(:,jm) = uif(:,jmm1)
          vif(:,jm) = vif(:,jmm1)
        end if

        uif = uif*dum
        vif = vif*dvm

        uib = ui
        ui  = uif
        vib = vi
        vi  = vif
        icb = ice
        ice = icf
        
      end subroutine ice_advance

      end module seaice