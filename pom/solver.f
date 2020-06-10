! solver.f

! main subroutines for POM

!______________________________________________________________________
!
      subroutine advave(curv2d)
!----------------------------------------------------------------------
!  Calculate horizontal advection and diffusion.
!----------------------------------------------------------------------
! called by: mode_interaction [advance.f]
!            mode_external    [advande.f]
!
! calls    : exchange2d_mpi   [parallel_mpi.f]
!______________________________________________________________________
!
      use config     , only: mode
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1
     &                     , n_south, n_west
      use grid       , only: aru, arv, dx, dy
      use glob_ocean , only: aam2d, advua, advva, cbc, d
     &                     , fluxua, fluxva, tps, ua, uab, va, vab
     &                     , wubot, wvbot

      implicit none


      real(rk), dimension(im,jm), intent(out) :: curv2d

      integer i,j

! u-advection and diffusion

! advective fluxes
      advua = 0.

      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do

! add viscous fluxes
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          tps(i,j)=.25*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
      call exchange2d_mpi(fluxua,im,jm)

      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do

! v-advection and diffusion
      advva = 0.

! advective fluxes
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do

      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do

! add viscous fluxes
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
      call exchange2d_mpi(fluxva,im,jm)

      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do

      if ( mode == 2 ) then

        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.25*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do
        call exchange2d_mpi(wubot,im,jm)

        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.25*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do
        call exchange2d_mpi(wvbot,im,jm)

        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do
        call exchange2d_mpi(curv2d,im,jm)

        do j=2,jmm1
          if(n_west.eq.-1) then
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          else
          do i=2,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          end if
        end do

        do i=2,imm1
          if(n_south.eq.-1) then
          do j=3,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          else
          do j=2,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          end if
        end do

      endif


      end ! subroutine advave
!
!______________________________________________________________________
!
      subroutine advct(xflux,yflux,curv)
!----------------------------------------------------------------------
!  Calculate the horizontal portions of momentum advection well in
! advance of their use in advu and advv so that their vertical
! integrals (created in the main program) may be used in the
! external (2-D) mode calculation.
!----------------------------------------------------------------------
! called by: lateral_viscosity [advance.f]
!
! calls    : exchange3d_mpi
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
     &                     , n_south, n_west
      use grid       , only: aru, arv, dx, dy
      use glob_ocean , only: aam, advx, advy, dt, u, ub, v, vb

      implicit none

      real(rk), dimension(im,jm,kb), intent(out) :: xflux, yflux, curv

      integer  i,j,k
      real(rk) dtaam


      curv  = 0.
      advx  = 0.
      xflux = 0.
      yflux = 0.

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25*( (v(i,j+1,k)+v(i,j,k))
     $                        *(dy(i+1,j)-dy(i-1,j))
     $                       -(u(i+1,j,k)+u(i,j,k))
     $                        *(dx(i,j+1)-dx(i,j-1)) )
     $                      /( dx(i,j)*dy(i,j) ) ! TODO: use `art` maybe
          end do
        end do
      end do
      call exchange3d_mpi(curv(:,:,1:kbm1),im,jm,kbm1)

! calculate x-component of velocity advection

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125*( (dt(i+1,j)+dt(i  ,j))*u(i+1,j,k)
     $                         +(dt(i  ,j)+dt(i-1,j))*u(i  ,j,k) )
     $                        *( u(i+1,j,k)+u(i,j,k) )
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125*( (dt(i  ,j)+dt(i  ,j-1))*v(i  ,j,k)
     $                         +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k) )
     $                        *( u(i,j,k)+u(i,j-1,k) )
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k) = xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam = .25*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $                 *( aam(i,j  ,k)+aam(i-1,j  ,k)
     $                   +aam(i,j-1,k)+aam(i-1,j-1,k) )
            yflux(i,j,k) = yflux(i,j,k)
     $                    -dtaam*( (ub(i,j,k)-ub(i,j-1,k))
     $                             /( dy(i,j  )+dy(i-1,j)
     $                               +dy(i,j-1)+dy(i-1,j-1) )
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                             /( dx(i,j  )+dx(i-1,j  )
     $                               +dx(i,j-1)+dx(i-1,j-1) ) )

            xflux(i,j,k) = dy(i,j)*xflux(i,j,k)
            yflux(i,j,k) = .25*( dx(i,j  )+dx(i-1,j  )
     $                          +dx(i,j-1)+dx(i-1,j-1) )*yflux(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k) = xflux(i,j  ,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i  ,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          if(n_west.eq.-1) then
            do i=3,imm1
              advx(i,j,k) = advx(i,j,k)
     $                     -aru(i,j)*.25
     $                      *( curv(i  ,j,k)*dt(i  ,j)
     $                         *(v(i  ,j+1,k)+v(i  ,j,k))
     $                        +curv(i-1,j,k)*dt(i-1,j)
     $                         *(v(i-1,j+1,k)+v(i-1,j,k)) )
            end do
          else
            do i=2,imm1
              advx(i,j,k) = advx(i,j,k)
     $                     -aru(i,j)*.25
     $                      *( curv(i  ,j,k)*dt(i  ,j)
     $                         *(v(i  ,j+1,k)+v(i  ,j,k))
     $                        +curv(i-1,j,k)*dt(i-1,j)
     $                         *(v(i-1,j+1,k)+v(i-1,j,k)) )
            end do
          end if
        end do
      end do

! calculate y-component of velocity advection

      advy  = 0.
      xflux = 0.
      yflux = 0.

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k) = .125*( (dt(i,j  )+dt(i-1,j  ))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k) )
     $                          *( v(i,j,k)+v(i-1,j,k) )
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k) = .125*( (dt(i,j+1)+dt(i,j  ))*v(i,j+1,k)
     $                           +(dt(i,j  )+dt(i,j-1))*v(i,j  ,k) )
     $                          *( v(i,j+1,k)+v(i,j,k) )
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam = .25*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $                 *( aam(i,j  ,k)+aam(i-1,j  ,k)
     $                   +aam(i,j-1,k)+aam(i-1,j-1,k) )
            xflux(i,j,k) = xflux(i,j,k)
     $                    -dtaam*( (ub(i,j,k)-ub(i,j-1,k))
     $                             /( dy(i,j  )+dy(i-1,j  )
     $                               +dy(i,j-1)+dy(i-1,j-1) )
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                             /( dx(i,j  )+dx(i-1,j  )
     $                               +dx(i,j-1)+dx(i-1,j-1) ) )
            yflux(i,j,k) = yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)

            xflux(i,j,k) = .25*( dy(i,j  )+dy(i-1,j  )
     $                          +dy(i,j-1)+dy(i-1,j-1) )*xflux(i,j,k)
            yflux(i,j,k) = dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k) = xflux(i+1,j,k)-xflux(i,j  ,k)
     $                   +yflux(i  ,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do

      do k=1,kbm1
        do i=2,imm1
          if(n_south.eq.-1) then
            do j=3,jmm1
              advy(i,j,k) = advy(i,j,k)
     $                     +arv(i,j)*.25
     $                      *( curv(i,j  ,k)*dt(i,j  )
     $                         *(u(i+1,j  ,k)+u(i,j  ,k))
     $                        +curv(i,j-1,k)*dt(i,j-1)
     $                         *(u(i+1,j-1,k)+u(i,j-1,k)) )
            end do
          else
            do j=2,jmm1
              advy(i,j,k) = advy(i,j,k)
     $                     +arv(i,j)*.25
     $                      *( curv(i,j  ,k)*dt(i,j  )
     $                         *(u(i+1,j  ,k)+u(i,j  ,k))
     $                        +curv(i,j-1,k)*dt(i,j-1)
     $                         *(u(i+1,j-1,k)+u(i,j-1,k)) )
            end do
          end if
        end do
      end do


      end ! subroutine advct
!
!______________________________________________________________________
!
      subroutine advq(qb,q,qf,xflux,yflux)
!----------------------------------------------------------------------
!  Calculates horizontal advection and diffusion, and vertical
! advection for turbulent quantities.
!----------------------------------------------------------------------
! called by: mode_internal [advance.f] (x2)
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: art, dum, dvm, dx, dy, dz, h
      use glob_ocean , only: aam, dt, etb, etf, u, v, w
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(in ) :: qb, q
      real(rk), dimension(im,jm,kb), intent(out) :: qf, xflux, yflux

      integer i,j,k


! do horizontal advection
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! do horizontal diffusion
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.25*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.25*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do

! do vertical advection, add flux terms, then step forward in time
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do


      end ! subroutine advq
!
!______________________________________________________________________
!
      subroutine advt1(fb,f,fclim,ff,xflux,yflux, var)
!----------------------------------------------------------------------
!  Integrate conservative scalar equations.
!  This is centred scheme, as originally provided in POM (previously
! called advt.)
!----------------------------------------------------------------------
! called by: mode_internal [advance.f] (x2)
!______________________________________________________________________
!
      use config     , only: tprni
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: art, dum, dvm, dx, dy, dz, h, zz
      use glob_ocean , only: aam, dt, etb, etf, tsurf, u, v, w, zflux
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(in   ) :: fclim
      real(rk), dimension(im,jm,kb), intent(inout) :: fb, f
      real(rk), dimension(im,jm,kb), intent(  out) :: ff,xflux,yflux
      character(len=1)             , intent(in   ) :: var

      integer       i,j,k
      real(kind=rk) relax !lyo:relax


      do j=1,jm
        do i=1,im
          f(i,j,kb)=f(i,j,kbm1)
          fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do

! do advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do

! add diffusive fluxes
      fb = fb - fclim

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.5*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.5*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do

      fb = fb + fclim

! do vertical advection
      do j=2,jmm1
        do i=2,imm1
!          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
!     for rivers 2010/5/08 ayumi
           if ( var == 'T' ) zflux(i,j,1)=tsurf(i,j)*w(i,j,1)*art(i,j)
           if ( var == 'S' ) zflux(i,j,1)=0.
           zflux(i,j,kb)=0.
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do

! add net horizontal fluxes and then step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
!lyo:relax to fclim in deep w/scales ~ 365days & 1000m. 2010/9/23
!fhx:test scale ~180 days & 500m. 2010/10/29
!exp301:!exp302:
!b420120325:relax=1.586e-8*(1.e0-exp(zz(k)*h(i,j)*5.e-4)) !730 days, 2000m
!           relax=3.171e-8*(1.e0-exp(zz(k)*h(i,j)*1.e-3)) !365 days, 1000m
            relax=6.43e-8*(1.-exp(zz(k)*h(i,j)*2.e-3)) !180 days, 500m
!           relax=0.0 !lyo:pac10:debug:
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
     $                 -relax*fclim(i,j,k)*dt(i,j)*art(i,j) !lyo:relax
            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)
     $                 -dti2*ff(i,j,k))
!    $                 /((h(i,j)+etf(i,j))*art(i,j))        !lyo:relax
     $                 /((h(i,j)+etf(i,j))*art(i,j)*(1.+relax*dti2))
          end do
        end do
      end do


      end ! subroutine advt1
!
!______________________________________________________________________
!
      subroutine advt2(fb,fclim,ff,xflux,yflux,var)
!----------------------------------------------------------------------
!  Integrate conservative scalar equations.
!  This is a first-order upstream scheme, which reduces implicit
! diffusion using the Smolarkiewicz iterative upstream scheme with an
! antidiffusive velocity.
!  It is based on the subroutines of Gianmaria Sannino
! (Inter-university Computing Consortium, Rome, Italy) and
! Vincenzo Artale (Italian National Agency for New Technology
! and Environment, Rome, Italy)
!----------------------------------------------------------------------
! called_by: mode_internal  [advance.f]
!
! calls    : bcast0d_mpi    [parallel_mpi.f]
!            exchange3d_mpi [parallel_mpi.f]
!            max0d_mpi      [parallel_mpi.f]
!            smol_adif      [solver.f]
!______________________________________________________________________
!
      use config     , only: nitera, tprni
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1, master_task
      use grid       , only: art, dum, dvm, dx, dy, dz, h
      use glob_ocean , only: aam, dt, etb, etf, tsurf, u, v, w, zflux
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(in   ) :: fclim
      real(rk), dimension(im,jm,kb), intent(inout) :: fb
      real(rk), dimension(im,jm,kb), intent(  out) :: ff,xflux,yflux
      character(len=1)             , intent(in   ) :: var

      integer                          i,j,k,itera
      real(rk), dimension(im,jm,kb) :: fbmem,xmassflux,ymassflux,zwflux
      real(rk), dimension(im,jm)    :: eta
      real(rk)                         eps, epsval  ! rwnd: iteration check


      select case ( var )
        case ('T')
          eps = 0.001
        case ('S')
          eps = 0.0001
        case default
          eps = 0.00001
      end select ! rwnd: iteration check

! calculate horizontal mass fluxes
      xmassflux = 0.
      ymassflux = 0.

      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=.25*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do

        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=.25*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do

! start Smolarkiewicz scheme
      do itera=1,nitera

! upwind advection scheme
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))

              yflux(i,j,k)=0.5
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do

        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.
!            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
!     for rivers 2010/5/08 ayumi
            if (itera == 1 ) then
              if ( var == 'T' )
     $              zflux(i,j,1)=tsurf(i,j)*w(i,j,1)*art(i,j)
              if ( var == 'S' )
     $              zflux(i,j,1)=0.
            end if
            zflux(i,j,kb)=0.
          end do
        end do

        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do

! add net advective fluxes and step forward in time
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
        ! next line added on 22-Jul-2009 by Raffaele Bernardello
        call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)

! calculate antidiffusion velocity
        call smol_adif(xmassflux,ymassflux,zwflux,ff)

        epsval = maxval(abs(ff(:,:,1:kbm1)-fbmem(:,:,1:kbm1))) ! rwnd: iteration check

        call max0d_mpi(epsval,master_task) ! rwnd: iteration check
        call bcast0d_mpi(epsval,master_task) ! rwnd: iteration check

        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do

        if (epsval < eps) exit ! rwnd: iteration check

! end of Smolarkiewicz scheme
      end do

!      write(*,*) my_task, ": ", var, ": nitera = ", itera  ! rwnd: iteration check

! add horizontal diffusive fluxes
      fb = fb - fclim

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.5/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.5/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do

      fb = fb + fclim

! add net horizontal fluxes and step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do


      end ! subroutine advt2
!
!______________________________________________________________________
!
      subroutine advtC(cbm,cf,ui,vi)
!----------------------------------------------------------------------
!  Integrate conservative scalar equations.
!  This is the 2-D modifiction of `advt2` for sea-ice concentration.
!----------------------------------------------------------------------
! called_by: ice_advance    [seaice.f]
!
! calls    : exchange2d_mpi [parallel_mpi.f]
!            smol_adifC     [solver.f]
!______________________________________________________________________
!
      use config     , only: nitera
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1
      use grid       , only: art, dx, dy
      use glob_ocean , only: etb, etf
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm), intent(in)  :: cbm, ui, vi
      real(rk), dimension(im,jm), intent(out) :: cf

      real(rk), dimension(im,jm) :: cbmem, eta, xflux, yflux
     &                            , xmassflux , ymassflux
!      real(kind=rk) eps, epsval  ! rwnd: iteration check
      integer i,j,itera


! calculate horizontal mass fluxes
      xmassflux = 0.
      ymassflux = 0.

      do j=2,jmm1
        do i=2,im
          xmassflux(i,j)=.5*(dy(i-1,j)+dy(i,j))*ui(i,j)
        end do
      end do

      do j=2,jm
        do i=2,imm1
          ymassflux(i,j)=.5*(dx(i,j-1)+dx(i,j))*vi(i,j)
        end do
      end do

      do j=1,jm
        do i=1,im
          eta(i,j)=etb(i,j)
        end do
      end do

      do j=1,jm
        do i=1,im
          cbmem(i,j)=cbm(i,j)
        end do
      end do

! start Smolarkiewicz scheme
      do itera=1,nitera

! upwind advection scheme
        do j=2,jm
          do i=2,im
            xflux(i,j) = .5
     $                      *((xmassflux(i,j)+abs(xmassflux(i,j)))
     $                       *cbmem(i-1,j)+
     $                        (xmassflux(i,j)-abs(xmassflux(i,j)))
     $                       *cbmem(i,j))

            yflux(i,j) = .5
     $                      *((ymassflux(i,j)+abs(ymassflux(i,j)))
     $                       *cbmem(i,j-1)+
     $                        (ymassflux(i,j)-abs(ymassflux(i,j)))
     $                       *cbmem(i,j))
          end do
        end do

! add net advective fluxes and step forward in time
        do j=2,jmm1
          do i=2,imm1
            cf(i,j) = xflux(i+1,j)-xflux(i,j)
     $               +yflux(i,j+1)-yflux(i,j)
            cf(i,j) = (cbmem(i,j)*art(i,j)
     $                   -dti2*cf(i,j))/art(i,j)
          end do
        end do
        ! next line added on 22-Jul-2009 by Raffaele Bernardello
        call exchange2d_mpi(cf,im,jm)

! calculate antidiffusion velocity
        call smol_adifC(xmassflux,ymassflux,cf)

!        epsval = maxval(abs(ff(:,:,1:kbm1)-fbmem(:,:,1:kbm1))) ! rwnd: iteration check
!
!        call max0d_mpi(epsval,master_task) ! rwnd: iteration check
!        call bcast0d_mpi(epsval,my_task) ! rwnd: iteration check

        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            cbmem(i,j)=cf(i,j)
          end do
        end do

!        if (epsval < eps) exit ! rwnd: iteration check

! end of Smolarkiewicz scheme
      end do

!      write(*,*) my_task, ": ", var, ": nitera = ", itera  ! rwnd: iteration check

! add horizontal diffusive fluxes
!      do j=2,jm
!        do i=2,im
!          xmassflux(i,j) = .5*(aam(i,j,1)+aam(i-1,j,1))
!          ymassflux(i,j) = .5*(aam(i,j,1)+aam(i,j-1,1))
!        end do
!      end do

!      do j=2,jm
!        do i=2,im
!         xflux(i,j) = -10.e2*tprni !-xmassflux(i,j)*tprni
!     $                   *(cb(i,j)-cb(i-1,j))*dum(i,j)
!     $                   *(dy(i,j)+dy(i-1,j))*0.5/(dx(i,j)+dx(i-1,j))
!         yflux(i,j) = -10.e2*tprni !-ymassflux(i,j)*tprni
!     $                   *(cb(i,j)-cb(i,j-1))*dvm(i,j)
!     $                   *(dx(i,j)+dx(i,j-1))*0.5/(dy(i,j)+dy(i,j-1))
!        end do
!      end do
!
!! add net horizontal fluxes and step forward in time
!      do j=2,jmm1
!        do i=2,imm1
!          cf(i,j) = cf(i,j)-dti2*(xflux(i+1,j)-xflux(i,j)
!     $                           +yflux(i,j+1)-yflux(i,j))
!     $                          /art(i,j)
!        end do
!      end do


      end ! subroutine advtC
!
!______________________________________________________________________
!
      subroutine advu
!----------------------------------------------------------------------
!  Do horizontal and vertical advection of u-momentum, and includes
! coriolis, surface slope and baroclinic terms.
!----------------------------------------------------------------------
! called by: mode_internal [advance.f]
!______________________________________________________________________
!
      use air        , only: e_atmos
      use glob_const , only: grav, rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: aru, cor, dy, dz, h
      use glob_ocean , only: advx, drhox, dt, egb, egf, etb, etf
     &                     , u, ub, uf, v, w
      use model_run  , only: dti2

      implicit none

      integer i,j,k


! do vertical advection
      uf = 0.

      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.25
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.125*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do

!  step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do


      end ! subroutine advu
!
!______________________________________________________________________
!
      subroutine advv
!----------------------------------------------------------------------
!  Do horizontal and vertical advection of v-momentum, and includes
! coriolis, surface slope and baroclinic terms.
!----------------------------------------------------------------------
! called by: mode_internal [advance.f]
!______________________________________________________________________
!
      use air        , only: e_atmos
      use glob_const , only: grav, rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: arv, cor, dx, dz, h
      use glob_ocean , only: advy, drhoy, dt, egb, egf, etb, etf
     &                     , u, v, vb, vf, w
      use model_run  , only: dti2

      implicit none

      integer i,j,k


! do vertical advection
      vf = 0.

      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.25
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.125*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do

! step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do


      end ! subroutine advv
!
!______________________________________________________________________
!
      subroutine baropg
!----------------------------------------------------------------------
!  Calculate  baroclinic pressure gradient.
!----------------------------------------------------------------------
! called by: pgscheme [advance.f]
!______________________________________________________________________
!
      use clim       , only: rmean
      use glob_const , only: grav, rk
      use glob_domain, only: imm1, jmm1, kb, kbm1
      use grid       , only: dum, dvm, dx, dy, zz
      use glob_ocean , only: drhox, drhoy, dt, rho
      use model_run  , only: ramp

      implicit none

      integer i,j,k


      rho = rho - rmean

! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.25*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.25*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.25*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.25*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

      rho = rho + rmean


      end ! subroutine baropg
!
!______________________________________________________________________
!fhx:Toni:npg
      subroutine baropg_mcc
!----------------------------------------------------------------------
!  Calculate baroclinic pressure gradient.
!  4th order correction terms, following McCalpin.
!----------------------------------------------------------------------
! called by: pgscheme [advance.f]
!______________________________________________________________________
!
      use clim       , only: rmean
      use glob_const , only: grav, rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1    ,my_task!REM:
     &                     , n_south, n_west
      use grid       , only: dum, dvm, dx, dy, dzz, zz
      use glob_ocean , only: d, drhox, drhoy, dt, density => rho
      use model_run  , only: ramp

      implicit none

      integer                              i     , j   , k
      real(rk), dimension(  im,  jm   ) :: d4    , ddx
      real(rk), dimension(  im,  jm,kb) :: drho  , rhou, rho
      real(rk), dimension(0:im,0:jm   ) :: d4th
      real(rk), dimension(0:im,0:jm,kb) :: rho4th


      rho = density - rmean

! convert a 2nd order matrices to special 4th order
! special 4th order case
      call order2d_mpi(d  ,d4th  ,im,jm)
      call order3d_mpi(rho,rho4th,im,jm,kb)

! compute terms correct to 4th order
      ddx  = 0.
      d4   = 0.
      rhou = 0.
      drho = 0.

! compute DRHO, RHOU, DDX and D4
      do j=1,jm
        do i=2,im
          do k=1,kbm1 ! TODO: Non-optimal loop nesting. Conduct performance comparisons.
            drho(i,j,k)=   (rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
          end do
          ddx(i,j)=   (d(i,j)-d(i-1,j))*dum(i,j)
          d4(i,j) =.5*(d(i,j)+d(i-1,j))*dum(i,j)
        end do
      end do

      if(n_west.eq.-1) then
        do j=1,jm
          do i=3,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i  ,j,k))-
     $                             2.*(rho(i  ,j,k)-rho(i-1,j,k))+
     $                     dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i  ,j,k)-rho(i+1,j,k))+
     $                     dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i  ,j))-
     $                        2.*(d(i  ,j)-d(i-1,j))+
     $                dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i  ,j)-d(i+1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
          end do
        end do
      else
        do j=1,jm
          do i=2,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                   (dum(i+1,j)*(rho(i+1,j,k)-rho   (i  ,j,k))-
     $                            2.*(rho(i  ,j,k)-rho   (i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i  ,j,k)-rho   (i+1,j,k))+
     $                     dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d   (i  ,j))-
     $                        2.*(d(i  ,j)-d   (i-1,j))+
     $                dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i  ,j)-d   (i+1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
          end do
        end do
      end if
! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do


      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                   +grav*0.5*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25*(dt(i,j)+dt(i-1,j))
     $                      *drhox(i,j,k)*dum(i,j)
     $                      *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

!!      print *, "BAROPG_MCC" !REM:
!!      if ( my_task==1 ) then
!!        print '(*(2(f15.7,";"),f15.7/))', (drhox(60,70,k),drho(60,70,k)
!!     &                            ,rhou(60,70,k), k=1,kbm1)
!!      end if
!!      call finalize_mpi
!!      stop

! compute terms correct to 4th order
      do i=1,im
        do j=1,jm
          ddx(i,j)=0.
          d4(i,j)=0.
        end do
      end do
      do k=1,kb
        do j=1,jm
          do i=1,im
            rhou(i,j,k)=0.
            drho(i,j,k)=0.
          end do
        end do
      end do

! compute DRHO, RHOU, DDX and D4
      do j=2,jm
        do i=1,im
          do k=1,kbm1
            drho(i,j,k)=   (rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
          end do
          ddx(i,j)=   (d(i,j)-d(i,j-1))*dvm(i,j)
          d4(i,j) =.5*(d(i,j)+d(i,j-1))*dvm(i,j)
        end do
      end do

      if(n_south.eq.-1) then
        do j=3,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j  ,k))-
     $                             2.*(rho(i,j  ,k)-rho(i,j-1,k))+
     $                     dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j  ,k)-rho(i,j+1,k))+
     $                     dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j  ))-
     $                        2.*(d(i,j  )-d(i,j-1))+
     $                dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j  )-d(i,j+1))+
     $               dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
          end do
        end do
      else
        do j=2,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho   (i,j  ,k))-
     $                             2.*(rho(i,j  ,k)-rho   (i,j-1,k))+
     $                     dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j  ,k)-rho   (i,j+1,k))+
     $                     dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dvm(i,j+1)*(d(i,j+1)-d   (i,j  ))-
     $                        2.*(d(i,j  )-d   (i,j-1))+
     $                dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j  )-d   (i,j+1))+
     $               dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
          end do
        end do
      end if

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                   +grav*0.5*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do


      end ! subroutine baropg_mcc
!
!______________________________________________________________________
!
      subroutine baropg_lin
!----------------------------------------------------------------------
!  Finite volume baroclinic pressure gradient scheme following Lin(*).
! TODO: fix and fill up the resource.
!----------------------------------------------------------------------
! called by: pgscheme [advance.f]
!______________________________________________________________________
!
      use glob_const , only: grav, rk ,rhoref!REM:
      use glob_domain, only: im, jm, kb, kbm1   ,my_task!REM:
      use grid       , only: aru,arv,dum, dvm, dx, dy, dz, dzz, h, z,zz
      use glob_ocean , only: drhox, drhoy, d, dt, el, et, density => rho
      use model_run  , only: ramp
      use clim, only: rmean

      implicit none

      integer       i,j,k
      !real(16) p(im,jm,kb),fx(im,jm,kb),fc(im,jm,kb)
      real(rk), dimension(im,jm,kb) :: p, fx, fc, rho
      real(rk) dh,cff,cff1


      rho = density - rmean

      p(:,:,1) = 0.
      do k = 1,kbm1
!        p(:,:,k+1) = p(:,:,k) + dt*dz(k)*(density(:,:,k))
        p(:,:,k+1) = d*dz(k)*(rho(:,:,k))
        fx(:,:,k) = .5*d*dz(k)*(p(:,:,k)+p(:,:,k+1))
      end do

      cff = .5*grav
      cff1= grav!/rhoref
!
!  Calculate pressure gradient in the XI-direction (m4/s2).
!
      fc(:,:,1) = 0.
      do k = 1,kbm1
        do j = 1,jm
          do i = 2,im
            if (dum(i,j)/=0.) then
              dh = z(k+1)*(d(i,j)-d(i-1,j))+el(i,j)-el(i-1,j)
              fc(i,j,k+1) = .5*dh*(p(i,j,k+1)+p(i-1,j,k+1))
              drhox(i,j,k) = cff*dz(k)*(dt(i-1,j)+dt(i,j))
     &                                 *( et(i-1,j)-et(i,j) )
     &                            + cff1*(fx(i-1,j,k)-
     &                                    fx(i  ,j,k)+
     &                                    fc(i,j,k  )-
     &                                    fc(i,j,k+1))
            end if
          end do
        end do
      end do
      do k = 1, kbm1
        drhox(2:im,:,k) = .5*(dy(2:im,:)+dy(1:im-1,:))*drhox(2:im,:,k)
      end do
!!      print *, "BAROPG_LIN" !REM:
!!      if ( my_task==1 ) then
!!        print *, "Di  :", d(25,25), h(25,25), el(25,25)
!!        print *, "Di-1:", d(24,25), h(24,25), el(24,25)
!!        print *, "dy:", dy(25,25), dy(24,25)
!!        do k = 1, kbm1
!!          fx(25,25,k) = z(k+1)*(h(25,25)-h(24,25))
!!        end do
!!        print '(*(3(f15.7,";"),f15.7/))', (fx(60,70,k),fc(60,70,k)
!!     &                            ,drhox(60,70,k),fx(60,70,k), k=1,kbm1)
!!      end if
!!      call finalize_mpi
!!      stop
!
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!
      fc(:,:,1) = 0.
      do k = 1,kbm1
        do j = 2,jm
          do i = 1,im
            if (dvm(i,j)/=0.) then
              dh = z(k+1)*(d(i,j)-d(i,j-1))+el(i,j)-el(i,j-1)
              fc(i,j,k+1) = .5*dh*(p(i,j,k+1)+p(i,j-1,k+1))
              drhoy(i,j,k) = cff*dz(k)*(dt(i,j-1)+dt(i,j))
     &                                 *( et(i,j-1)-et(i,j) )
     &                            + cff1*(fx(i,j-1,k)-
     &                                    fx(i,j  ,k)+
     &                                    fc(i,j,k  )-
     &                                    fc(i,j,k+1))
            end if
          end do
        end do
      end do
      do k = 1, kbm1
        drhoy(:,2:jm,k) = .5*(dx(:,2:jm)+dx(:,1:jm-1))*drhoy(:,2:jm,k)
      end do

      drhox = ramp*drhox
      drhoy = ramp*drhoy

!      if (ramp > 0) then
!        write(40+my_task,*) iint, drhox(50,50,:),drhoy(50,50,:)
!        write(50+my_task,*) iint, et(50,50), et(49,50), et(50,49),
!     &                   h(50,50), h(49,50), h(50,49), dz,
!     &                   dx(50,50), dx(49,50), dx(50,49),
!     &                   dy(50,50), dy(49,50), dy(50,49),
!     &               rho(50,50,:), "|", rho(49,50,:), "|", rho(50,49,:)
!      end if
!      if (iint >= 1) then
!        if (my_task==0) then
!          print *, iint, "z:     ", z(2)
!          print *, iint, "dz:    ", dz(2)
!          print *, iint, "hi-1:  ", h(49,50)
!          print *, iint, "h:     ", h(50,50)
!          print *, iint, "dy-1:  ", dy(49,50)
!          print *, iint, "dy:    ", dy(50,50)
!          print *, iint, "hzi-1: ", dz(2)*h(49,50)
!          print *, iint, "hzi:   ", dz(2)*h(50,50)
!          print *, iint, "rhoi-1:", rho(49,50,2)
!          print *, iint, "rho:   ", rho(50,50,2)
!          print *, iint, "pi-1:  ", p(49,50,2)
!          print *, iint, "p:     ", p(50,50,2)
!          print *, iint, "fxi-1: ", fx(49,50,2)
!          print *, iint, "fx:    ", fx(50,50,2)
!          print *, iint, "fc:    ", fc(50,50,2)
!          print *, z(3)*(dt(50,50)-dt(49,50))+et(50,50)-et(49,50)
!          print *, iint, "drhox: ", drhox(50,50,2)
!        end if
!        call finalize_mpi
!      end if
!      if (iint > 10) call finalize_mpi

      ! rho = rho+rmean


      end ! subroutine baropg_lin
!
!______________________________________________________________________
!
      subroutine baropg_song_std
!----------------------------------------------------------------------
!  Baroclinic pressure gradient scheme following Song(*).
! TODO: fix and fill up the resource.
!----------------------------------------------------------------------
! called by: pgscheme [advance.f]
!______________________________________________________________________
!
      use glob_const , only: grav, rhoref, rk
      use glob_domain, only: im, jm, kb, kbm1   ,my_task!:REM
      use grid       , only: dum, dvm, dx, dy, dz, z, zz  ,h!:REM
      use glob_ocean , only: d, drhox, drhoy, dt, density => rho  ,el!:REM
      use model_run  , only: ramp
      use clim,only:rmean

      implicit none

      integer       i,j,k
      real(rk) phix(2:im),phie(im)
      real(rk) fac,fac1,fac2,fac3,cff1,cff2,cff3,cff4,gamma
      real(rk) rho(im,jm,kb)


      rho = density - rmean

      fac  = 10000.      /rhoref
      fac1 =     .5 *grav!/rhoref
      fac2 = 1000.  *grav/rhoref
      fac3 =     .25*grav!/rhoref

      do j = 1,jm
        do i = 2,im
          cff1 = -zz(1)*(d(i,j)+d(i-1,j))
          phix(i) = fac1*(rho(i,j,1)-rho(i-1,j,1))*cff1
!          phix(i) = phix(i) + fac*(e_atmos(i,j)-e_atmos(i-1,j))
          phix(i) = phix(i)+                                              &
     &            (fac2+fac1*(rho(i,j,1)+rho(i-1,j,1)))*
     &            (el(i,j)-el(i-1,j))
          drhox(i,j,1) = -.25*dz(1)*(dt(i,j)+dt(i-1,j))*
     &                      phix(i)*(dy(i,j)+dy(i-1,j))*dum(i,j)
         end do
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
        do k = 2,kbm1
          do i = 2,im
            cff1 = 1./(d(i  ,j)*(zz(k-1)-zz(k))*
     &                 d(i-1,j)*(zz(k-1)-zz(k)))
            cff2 = (d(i,j)-d(i-1,j))*(zz(k)+zz(k-1))
     &            + 2.*(el(i,j)-el(i-1,j))
            cff3 = (zz(k-1)-zz(k))*(d(i,j)-d(i-1,j))
            gamma = .125*cff1*cff2*cff3

            cff1 = (1.+gamma)*(rho(i,j,k-1)-rho(i-1,j,k-1))+
     &             (1.-gamma)*(rho(i,j,k  )-rho(i-1,j,k  ))
            cff2 = rho(i,j,k-1)+rho(i-1,j,k-1)-                           &
     &             rho(i,j,k  )-rho(i-1,j,k  )
            cff3 = (d(i,j)+d(i-1,j))*(zz(k-1)-zz(k))
            cff4 = (1.+gamma)*
     &               (zz(k-1)*(d(i,j)-d(i-1,j))+el(i,j)-el(i-1,j))+
     &             (1.-gamma)*
     &               (zz(k  )*(d(i,j)-d(i-1,j))+el(i,j)-el(i-1,j))
            phix(i) = phix(i)+                                            &
     &                fac3*(cff1*cff3-cff2*cff4)
!
!            cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+                           &
!     &           rho(i,j,k  )-rho(i-1,j,k  )
!            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
!     &           rho(i,j,k  )-rho(i-1,j,k  )
!            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
!     &           z_r(i,j,k  )-z_r(i-1,j,k  )
!            cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+                           &
!     &           z_r(i,j,k  )-z_r(i-1,j,k  )
!            phix(i)=phix(i)+                                            &
!     &              fac3*(cff1*cff3-cff2*cff4)
            drhox(i,j,k) = -.25*dz(k)*(dt(i,j)+dt(i-1,j))*
     &                        phix(i)*(dy(i,j)+dy(i-1,j))*dum(i,j)
          end do
        end do
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
        if (j>=2) then
          do i = 1,im
            cff1 = (z(1)-zz(1))*(d(i,j)+d(i,j-1))
            phie(i) = fac1*(rho(i,j,1)-rho(i,j-1,1))*cff1
!            phie(i) = phie(i) + fac*(e_atmos(i,j)-e_atmos(i,j-1))
            phie(i) = phie(i)+                                            &
     &              (fac2+fac1*(rho(i,j,1)+rho(i,j-1,1)))*
     &              (z(1)*(d(i,j)-d(i,j-1)))! + el(i,j)-el(i,j-1))
            drhoy(i,j,1) = -.25*dz(1)*(dt(i,j)+dt(i,j-1))*
     &                        phie(i)*(dy(i,j)+dy(i,j-1))*dvm(i,j)
          end do
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
          do k = 2,kbm1
            do i = 1,im
              cff1 = 1./(d(i,j  )*(z(k-1)-z(k))*
     &                   d(i,j-1)*(z(k-1)-z(k)))
              cff2 = (d(i,j)-d(i,j-1))*(zz(k)+zz(k-1))
     &              + 2.*(el(i,j)-el(i,j-1))
              cff3 = (zz(k-1)-zz(k))*(d(i,j)-d(i,j-1))
              gamma = .125*cff1*cff2*cff3

              cff1 = (1.+gamma)*(rho(i,j,k-1)-rho(i,j-1,k-1))+
     &               (1.-gamma)*(rho(i,j,k  )-rho(i,j-1,k  ))
              cff2 = rho(i,j,k-1)+rho(i,j-1,k-1)-                         &
     &               rho(i,j,k  )-rho(i,j-1,k  )
              cff3 = (d(i,j)+d(i,j-1))*(zz(k-1)-zz(k))
              cff4 = (1.+gamma)*
     &                (zz(k-1)*(d(i,j)-d(i,j-1))+el(i,j)-el(i,j-1))+
     &               (1.-gamma)*
     &                (zz(k  )*(d(i,j)-d(i,j-1))+el(i,j)-el(i,j-1))
              phie(i) = phie(i)+                                          &
     &                  fac3*(cff1*cff3-cff2*cff4)
!
!              cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+                         &
!     &             rho(i,j,k  )-rho(i,j-1,k  )
!              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
!     &             rho(i,j,k  )-rho(i,j-1,k  )
!              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
!     &             z_r(i,j,k  )-z_r(i,j-1,k  )
!              cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+                         &
!     &             z_r(i,j,k  )-z_r(i,j-1,k  )
!              phie(i)=phie(i)+                                          &
!     &                fac3*(cff1*cff3-cff2*cff4)
              drhoy(i,j,k) = -.25*dz(k)*(dt(i,j)+dt(i,j-1))*
     &                          phie(i)*(dx(i,j)+dx(i,j-1))*dvm(i,j)
!              if (isnan(drhoy(i,j,k))) write(*,*) my_task,"::",i,j,k
            end do
          end do
        end if
      end do

!!      print *, "BAROPG_SONG" !REM:
!!      if ( my_task==1 ) then
!!        print *, "Di  :", d(25,25), h(25,25), el(25,25)
!!        print *, "Di-1:", d(24,25), h(24,25), el(24,25)
!!        print *, "dy:", dy(25,25), dy(24,25)
!!        print '(*(1(f15.7,";"),f15.7/))', (drhox(60,70,k)
!!     &                            ,drhox(60,70,k), k=1,kbm1)
!!      end if
!!      call finalize_mpi
!!      stop

      drhox = ramp*drhox
      drhoy = ramp*drhoy

!      if (ramp > 0) then
!        write(40+my_task,*) iint, drhox(50,50,:)
!        write(50+my_task,*) iint, el(50,50), rho(50,50,:)
!      end if
!      if (iint > 10) call finalize_mpi

!      rho = rho+rmean


      end ! subroutine baropg_song_std
!
!______________________________________________________________________
!
      subroutine baropg_shch
!----------------------------------------------------------------------
!  Baroclinic pressure gradient scheme following Shchepetkin(*).
! TODO: fix and fill up the resource.
!----------------------------------------------------------------------
! called by: pgscheme [advance.f]
!______________________________________________________________________
!
      use glob_const , only: grav, rhoref, rk
      use glob_domain, only: im, imm1, jm, jmm1, kb  ,my_task!REM:
      use grid       , only: dum, dvm, dx, dy, dz, z, zz
      use glob_ocean , only: d, drhox, drhoy, dt, density => rho, el
      use model_run  , only: ramp
      use clim       , only: rmean

      implicit none

      integer                     i,j,k
      real(rk), parameter :: OneFifth   = .2
     &                      ,OneTwelfth = 1./12.
     &                      ,eps        = 1.e-10
      real(rk)               GRho, GRho0, HalfGRho
      real(rk)               fac,cff,cff1,cff2
      real(rk), dimension(im,jm,kb):: p, rho
      real(rk), dimension(im,kb+1) :: idR, idZ
      real(rk), dimension(im,jm)   :: fc, aux, idRx, idZx


      rho = density - rmean

!
!-----------------------------------------------------------------------
!  Preliminary step (same for XI- and ETA-components:
!-----------------------------------------------------------------------
!
      GRho  = grav!/rhoref
      GRho0 = 1000.*GRho/rhoref
      HalfGRho = .5*GRho
      fac = 100./rhoref

      do j = 1,jm
        do k = 2,kb
          do i = 1,im
            idR(i,k) = rho(i,j,k-1) - rho(i,j,k)
            idZ(i,k) = (zz(k-1)-zz(k))*d(i,j)
          end do
        end do
        do i = 1,im
          idR(i,1) = idR(i,2)
          idZ(i,1) = idZ(i,2)
          idR(i,kb+1) = idR(i,kb)
          idZ(i,kb+1) = idZ(i,kb)
        end do
        do k=1,kb
          do i=1,im
            cff = 2.*idR(i,k)*idR(i,k+1)
            if (cff > eps) then
              idR(i,k) = cff/(idR(i,k)+idR(i,k+1))
            else
              idR(i,k) = 0.
            end if
            idZ(i,k) = 2.*idZ(i,k)*idZ(i,k+1)/(idZ(i,k)+idZ(i,k+1))
          end do
        end do
        do i=1,im
          cff1 = 1./(d(i,j)*(zz(1)-zz(2)))
          cff2 = .5*(rho(i,j,1)-rho(i,j,2))*d(i,j)*(-zz(1))*cff1
          p(i,j,1) = Grho0*  el(i,j)
     &             + Grho *(rho(i,j,1)+cff2)*d(i,j)*(-zz(1))
        end do
        do k = 2,kb
          do i = 1,im
            p(i,j,k) = p(i,j,k-1) +
     &                HalfGRho*((rho(i,j,k-1)+rho(i,j,k))*
     &                          (d(i,j)*(zz(k-1)-zz(k)))-
     &                          OneFifth*
     &                          ((idR(i,k-1)-idR(i,k))*
     &                           (d(i,j)*(zz(k-1)-zz(k))-
     &                            OneTwelfth*
     &                            (idZ(i,k-1)+idZ(i,k)))-
     &                           (idZ(i,k-1)-idZ(i,k))*
     &                           (rho(i,j,k-1)-rho(i,j,k)-
     &                            OneTwelfth*
     &                            (idR(i,k-1)+idR(i,k)))))
          end do
        end do
      end do
!
!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      do k = 1,kb
        do j = 1,jm
          do i = 2,im
            aux(i,j) = zz(k)*(d(i,j)-d(i-1,j))+el(i,j)-el(i-1,j)
            fc(i,j) = rho(i,j,k)-rho(i-1,j,k)
          end do
        end do
        do j = 1,jm
          do i = 1,imm1
            cff = 2.*aux(i,j)*aux(i+1,j)
            if (cff > eps) then
              cff1 = 1./(aux(i,j)+aux(i+1,j))
              idZx(i,j) = cff*cff1
            else
              idZx(i,j) = 0.
            end if
            cff1 = 2.*fc(i,j)*fc(i+1,j)
            if (cff1 > eps) then
              cff2 = 1./(fc(i,j)+fc(i+1,j))
              idRx(i,j) = cff1*cff2
            else
              idRx(i,j) = 0.
            end if
          end do
        end do

        do j = 1,jm
          do i = 2,im
            if (k==1) then
              drhox(i,j,k) = 0.
            else
              drhox(i,j,k) = drhox(i,j,k-1)
            end if
            drhox(i,j,k) = drhox(i,j,k)+
     &                    .25*(dy(i,j)+dy(i-1,j))*
     &                    (dz(k)*(dt(i,j)+dt(i-1,j)))*
     &                    (p(i-1,j,k)-p(i,j,k)-
     &                     HalfGRho*
     &                     ((rho(i,j,k)+rho(i-1,j,k))*
     &                      (zz(k)*(d(i,j)-d(i-1,j))
     &                             +el(i,j)-el(i-1,j))
     &                      -OneFifth*
     &                       ((idRx(i,j)-idRx(i-1,j))*
     &                        (zz(k)*(d(i,j)-d(i-1,j))
     &                               +el(i,j)-el(i-1,j)
     &                        -OneTwelfth*
     &                         (idZx(i,j)+idZx(i-1,j)))-
     &                        (idZx(i,j)-idZx(i-1,j))*
     &                        (rho(i,j,k)-rho(i-1,j,k)
     &                        -OneTwelfth*
     &                         (idRx(i,j)+idRx(i-1,j))))))
!     &                     /(.5*(dt(i,j)+dt(i-1,j)))
          end do
        end do
      end do
!!      do k=1,kb-1
!!        do j=2,jmm1
!!          do i=2,imm1
!!            drhox(i,j,k)=.25*(dt(i,j)+dt(i-1,j))
!!     $                      *drhox(i,j,k)*dum(i,j)
!!     $                      *(dy(i,j)+dy(i-1,j))
!!          end do
!!        end do
!!      end do
      print *, "BAROPG_SHCH" !REM:
      if ( my_task==1 ) then
!!        print *, "Di  :", d(25,25), h(25,25), el(25,25)
!!        print *, "Di-1:", d(24,25), h(24,25), el(24,25)
!!        print *, "dy:", dy(25,25), dy(24,25)
        print '(*(2(f18.7,";"),f18.7/))', (drhox(60,70,k)
     &                    ,drhox(60,70,k)/dt(60,70),dt(60,70), k=1,kb-1)
      end if
!!      call finalize_mpi
!!      stop
!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      do k = 1,kb
        do j = 2,jm
          do i = 1,im
            aux(i,j) = zz(k)*(d(i,j)-d(i,j-1))
            fc(i,j) = rho(i,j,k)-rho(i,j-1,k)
          end do
        end do
        do j = 1,jmm1
          do i = 1,im
            cff = 2.*aux(i,j)*aux(i,j+1)
            if (cff > eps) then
              cff1 = 1./(aux(i,j)+aux(i,j+1))
              idZx(i,j) = cff*cff1
            else
              idZx(i,j) = 0.
            end if
            cff1 = 2.*fc(i,j)*fc(i,j+1)
            if (cff1 > eps) then
              cff2 = 1./(fc(i,j)+fc(i,j+1))
              idRx(i,j) = cff1*cff2
            else
              idRx(i,j) = 0.
            end if
          end do
        end do

        do j = 2,jm
          do i = 1,im
            if (k==1) then
              drhoy(i,j,k) = 0.
            else
              drhoy(i,j,k) = drhoy(i,j,k-1)
            end if
            drhoy(i,j,k) = drhoy(i,j,k)+
     &                    .25*(dx(i,j)+dx(i,j-1))*
     &                    (dz(k)*(dt(i,j)+dt(i,j-1)))*
     &                    (p(i,j-1,k)-p(i,j,k)-
     &                     HalfGRho*
     &                     ((rho(i,j,k)+rho(i,j-1,k))*
     &                      (zz(k)*(d(i,j)-d(i,j-1))
     &                             +el(i,j)-el(i,j-1))
     &                      -OneFifth*
     &                       ((idRx(i,j)-idRx(i,j-1))*
     &                        (zz(k)*(d(i,j)-d(i,j-1))
     &                               +el(i,j)-el(i,j-1)
     &                        -OneTwelfth*
     &                         (idZx(i,j)+idZx(i,j-1)))-
     &                        (idZx(i,j)-idZx(i,j-1))*
     &                        (rho(i,j,k)-rho(i,j-1,k)
     &                        -OneTwelfth*
     &                         (idRx(i,j)+idRx(i,j-1))))))
!     &                     /(.5*(dt(i,j)+dt(i,j-1)))
          end do
        end do
      end do
!


      drhox = ramp*drhox
      drhoy = ramp*drhoy


      end ! subroutine baropg_shch
!
!______________________________________________________________________
!
      subroutine xperi2d_mpi(wrk,nx,ny) !lyo:20110224:alu:stcc:
!----------------------------------------------------------------------
!  Doing periodic bc in x.
!  Pass from east to west and also pass from west to east.
!----------------------------------------------------------------------
! called by: bc_vel_ext [bry.f90] (x2)
!            bc_vel_int [bry.f90] (x2)
!            bc_zeta    [bry.f90]
!
! calls    : mpi_recv   [mpi]
!            mpi_send   [mpi]
! TODO: move to `parallel_mpi.f`
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: im_global, im_local, my_task
     &                     , n_east, n_west, POM_COMM
      use mpi        , only: MPI_REAL, MPI_STATUS_SIZE
     &                     , mpi_recv, mpi_send

      implicit none

      integer , intent(in   ) :: nx,ny
      real(rk), intent(inout) :: wrk(nx,ny)

      integer  j,ierr,nproc_x
      integer  dest_task,sour_task
      integer  istatus(mpi_status_size)
      real(rk) sendbuf(ny),recvbuf(ny)


! determine the number of processors in x
      if(mod(im_global-2,im_local-2).eq.0) then
        nproc_x=(im_global-2)/(im_local-2)
      else
        nproc_x=(im_global-2)/(im_local-2) + 1
      end if


!lyo:scs1d:
      if (nproc_x.eq.1) then
        do j=1,ny
        wrk(nx,j)=wrk(3,j); wrk(1,j)=wrk(nx-2,j); wrk(2,j)=wrk(nx-1,j)
        enddo
      else
!
! The most east sudomains
      if(n_east.eq.-1) then
        dest_task=my_task-nproc_x+1
        sour_task=my_task-nproc_x+1

       ! first time to send
         do j=1,ny
           sendbuf(j)=wrk(nx-2,j)
         end do
! TODO: change MPI_REAL to MPI_RK and move the whole subroutine to parallel_mpi.f
         call mpi_send(sendbuf,ny,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)


       !first time to recieve
         call mpi_recv(recvbuf,ny,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do j=1,ny
          wrk(nx,j)=recvbuf(j)
         end do

       ! second time to send
         do j=1,ny
           sendbuf(j)=wrk(nx-1,j)
         end do
         call mpi_send(sendbuf,ny,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

      endif !if(n_east.eq.-1)

! The most west sudomains
      if(n_west.eq.-1) then
        sour_task=my_task+nproc_x-1
        dest_task=my_task+nproc_x-1

        ! first time to recieve
         call mpi_recv(recvbuf,ny,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do j=1,ny
           wrk(1,j)=recvbuf(j)
         end do


        ! first time to send
         do j=1,ny
           sendbuf(j)=wrk(3,j)
         end do
         call mpi_send(sendbuf,ny,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! second time to recieve
         call mpi_recv(recvbuf,ny,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)


         do j=1,ny
           wrk(2,j)=recvbuf(j)
         end do

      endif !if(n_west.eq.-1)

      endif !if (nproc_x.eq.1) then !lyo:scs1d:


      end ! subroutine xperi2d_mpi
!
!______________________________________________________________________
!
      subroutine xperi3d_mpi(wrk,nx,ny,nz) !lyo:20110224:alu:stcc:
!----------------------------------------------------------------------
!  Doing periodic bc in x.
!  Pass from east to west and also pass from west to east.
!----------------------------------------------------------------------
! called by: bc_ts      [bry.f90] (x2)
!            bc_turb    [bry.f90] (x6)
!            bc_vel_int [bry.f90] (x2)
!            bc_zeta    [bry.f90]
!
! calls    : mpi_recv   [mpi]
!            mpi_send   [mpi]
! TODO: move to `parallel_mpi.f`
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: im_global, im_local, my_task
     &                     , n_east, n_west, POM_COMM
      use mpi        , only: MPI_REAL, MPI_STATUS_SIZE
     &                     , mpi_recv, mpi_send

      implicit none

      integer , intent(in   ) :: nx,ny,nz
      real(rk), intent(inout) :: wrk(nx,ny,nz)

      integer  i,j,k
      integer  ierr,nproc_x
      integer  dest_task,sour_task
      integer  istatus(mpi_status_size)
      real(rk) sendbuf(ny*nz),recvbuf(ny*nz)


! determine the number of processors in x
      if(mod(im_global-2,im_local-2).eq.0) then
        nproc_x=(im_global-2)/(im_local-2)
      else
        nproc_x=(im_global-2)/(im_local-2) + 1
      end if

!lyo:scs1d:
      if (nproc_x.eq.1) then
        do k=1,nz; do j=1,ny
        wrk(nx,j,k)=wrk(3,j,k);
        wrk(1,j,k)=wrk(nx-2,j,k); wrk(2,j,k)=wrk(nx-1,j,k);
        enddo; enddo
      else
!
!The most east sudomains
      if(n_east.eq.-1) then
        dest_task=my_task-nproc_x+1
        sour_task=my_task-nproc_x+1

        ! first time to send
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
           sendbuf(i)=wrk(nx-2,j,k)
          end do
         end do
         call mpi_send(sendbuf,ny*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! first time to recieve
         call mpi_recv(recvbuf,ny*nz,MPI_REAL,sour_task,sour_task,
     $                 POM_COMM,istatus,ierr)
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
           wrk(nx,j,k)=recvbuf(i)
          end do
         end do

        ! second time to send
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
           sendbuf(i)=wrk(nx-1,j,k)
          end do
         end do
         call mpi_send(sendbuf,ny*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

      endif!if(n_east.eq.-1)

!The most west sudomains
      if(n_west.eq.-1) then
       sour_task=my_task+nproc_x-1
       dest_task=my_task+nproc_x-1

        ! first time to recieve
         call mpi_recv(recvbuf,ny*nz,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
            wrk(1,j,k)=recvbuf(i)
          end do
         end do


        ! first time to send
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
           sendbuf(i)=wrk(3,j,k)
          end do
         end do
         call mpi_send(sendbuf,ny*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! second time to recieve
         call mpi_recv(recvbuf,ny*nz,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do k=1,nz
          do j=1,ny
           i=j+(k-1)*ny
            wrk(2,j,k)=recvbuf(i)
          end do
         end do

      endif!if(n_west.eq.-1)

      endif !if (nproc_x.eq.1) then !lyo:scs1d:


      end ! subroutine xperi3d_mpi
!
!_______________________________________________________________________
!lyo:scs1d:add yperi*:ipery:beg:
      subroutine yperi2d_mpi(wrk,nx,ny)
!----------------------------------------------------------------------
!  Doing periodic bc in y.
!  Pass from north to south and also pass from south to north.
!----------------------------------------------------------------------
! called by: bc_vel_ext [bry.f90] (x2)
!            bc_vel_int [bry.f90] (x2)
!            bc_zeta    [bry.f90]
!
! calls    : mpi_recv   [mpi]
!            mpi_send   [mpi]
! TODO: move to `parallel_mpi.f`
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: jm_global, jm_local, my_task
     &                     , n_north, n_south, POM_COMM
      use mpi        , only: MPI_REAL, MPI_STATUS_SIZE
     &                     , mpi_recv, mpi_send

      implicit none

      integer , intent(in   ) :: nx, ny
      real(rk), intent(inout) :: wrk(nx,ny)

      integer  i,ierr,nproc_y
      integer  dest_task,sour_task
      integer  istatus(mpi_status_size)
      real(rk) sendbuf(nx),recvbuf(nx) !lyo:scs1d:


! determine the number of processors in y
      if(mod(jm_global-2,jm_local-2).eq.0) then
        nproc_y=(jm_global-2)/(jm_local-2)
      else
        nproc_y=(jm_global-2)/(jm_local-2) + 1
      end if

!lyo:scs1d:
      if (nproc_y.eq.1) then
        do i=1,nx
        wrk(i,ny)=wrk(i,3); wrk(i,1)=wrk(i,ny-2); wrk(i,2)=wrk(i,ny-1)
        enddo
      else
!
!The most north sudomains
      if(n_north.eq.-1) then
        dest_task=my_task-nproc_y+1
        sour_task=my_task-nproc_y+1

       ! first time to send
         do i=1,nx
           sendbuf(i)=wrk(i,ny-2)
         end do
         call mpi_send(sendbuf,nx,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)


       !first time to recieve
         call mpi_recv(recvbuf,nx,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do i=1,nx
          wrk(i,ny)=recvbuf(i)
         end do

       ! second time to send
         do i=1,nx
           sendbuf(i)=wrk(i,ny-1)
         end do
         call mpi_send(sendbuf,nx,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

      endif !if(n_north.eq.-1)

!The most south sudomains
      if(n_south.eq.-1) then
        sour_task=my_task+nproc_y-1
        dest_task=my_task+nproc_y-1

        ! first time to recieve
         call mpi_recv(recvbuf,nx,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do i=1,nx
           wrk(i,1)=recvbuf(i)
         end do


        ! first time to send
         do i=1,nx
           sendbuf(i)=wrk(i,3)
         end do
         call mpi_send(sendbuf,nx,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! second time to recieve
         call mpi_recv(recvbuf,nx,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)


         do i=1,nx
           wrk(i,2)=recvbuf(i)
         end do

      endif !if(n_south.eq.-1)

      endif !if (nproc_y.eq.1) then !lyo:scs1d:


      end ! subroutine yperi2d_mpi
!
!______________________________________________________________________
!
      subroutine yperi3d_mpi(wrk,nx,ny,nz)
!----------------------------------------------------------------------
!  Doing periodic bc in y.
!  Pass from north to south and also pass from south to north.
!----------------------------------------------------------------------
! called by: bc_ts      [bry.f90] (x2)
!            bc_turb    [bry.f90] (x6)
!            bc_vel_int [bry.f90] (x2)
!            bc_zeta    [bry.f90]
!
! calls    : mpi_recv   [mpi]
!            mpi_send   [mpi]
! TODO: move to `parallel_mpi.f`
!______________________________________________________________________
!
      use glob_const , only: rk
      use glob_domain, only: jm_global, jm_local, my_task
     &                     , n_north, n_south, POM_COMM
      use mpi        , only: MPI_REAL, MPI_STATUS_SIZE
     &                     , mpi_recv, mpi_send

      implicit none

      integer , intent(in   ) :: nx,ny,nz
      real(rk), intent(inout) :: wrk(nx,ny,nz)

      integer  i,j,k
      integer  ierr,nproc_y
      integer  dest_task,sour_task
      integer  istatus(mpi_status_size)
      real(rk) sendbuf(nx*nz),recvbuf(nx*nz)


! determine the number of processors in y
      if(mod(jm_global-2,jm_local-2).eq.0) then
        nproc_y=(jm_global-2)/(jm_local-2)
      else
        nproc_y=(jm_global-2)/(jm_local-2) + 1
      end if

!lyo:scs1d:
      if (nproc_y.eq.1) then
        do k=1,nz; do i=1,nx
        wrk(i,ny,k)=wrk(i,3,k);
        wrk(i,1,k)=wrk(i,ny-2,k); wrk(i,2,k)=wrk(i,ny-1,k)
        enddo; enddo
      else
!
!The most north sudomains
      if(n_north.eq.-1) then
        dest_task=my_task-nproc_y+1
        sour_task=my_task-nproc_y+1

        ! first time to send
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
           sendbuf(j)=wrk(i,ny-2,k)
          end do
         end do
         call mpi_send(sendbuf,nx*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! first time to recieve
         call mpi_recv(recvbuf,nx*nz,MPI_REAL,sour_task,sour_task,
     $                 POM_COMM,istatus,ierr)
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
           wrk(i,ny,k)=recvbuf(j)
          end do
         end do

        ! second time to send
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
           sendbuf(j)=wrk(i,ny-1,k)
          end do
         end do
         call mpi_send(sendbuf,nx*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

      endif!if(n_north.eq.-1)

!The most south sudomains
      if(n_south.eq.-1) then
       sour_task=my_task+nproc_y-1
       dest_task=my_task+nproc_y-1

        ! first time to recieve
         call mpi_recv(recvbuf,nx*nz,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
            wrk(i,1,k)=recvbuf(j)
          end do
         end do


        ! first time to send
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
           sendbuf(j)=wrk(i,3,k)
          end do
         end do
         call mpi_send(sendbuf,nx*nz,MPI_REAL,dest_task,my_task,
     $                   POM_COMM,ierr)

        ! second time to recieve
         call mpi_recv(recvbuf,nx*nz,MPI_REAL,sour_task,sour_task,
     $                POM_COMM,istatus,ierr)
         do k=1,nz
          do i=1,nx
           j=i+(k-1)*nx
            wrk(i,2,k)=recvbuf(j)
          end do
         end do

      endif!if(n_south.eq.-1)

      endif !if (nproc_y.eq.1) then !lyo:scs1d:


      end ! subroutine yperi3d_mpi
!
!______________________________________________________________________
!lyo:scs1d:add yperi*:ipery:end:
!lyo:pac10:beg:
!______________________________________________________________________
!
      subroutine bcond_nwatl(idx, ampe, phae, amue, phue)
!----------------------------------------------------------------------
! boundary conditions for the Mid-Atlantic Bight/NWAtl model.
! only the eastern boundary is open.
! annual mean uabe is specified (2-d velocity b.c.).
! radiation b.c. for 3-d uv,and advection b.c. for 3-d ts.
!----------------------------------------------------------------------
! called by: [NO CALLS]
!
! calls    : bcast0d_mpi [parallel_mpi.f]
!            psum_mpi    [parallel_mpi.f]
!______________________________________________________________________
!
      use bry
      use config     , only: calc_tide, ntide
      use glob_const , only: grav, pi, rk, small
      use glob_domain
      use grid       , only: dx, dy, dum, dvm, fsm, h
      use glob_ocean , only: elf, s, t, u, uaf, uf, vaf, vf, w
      use model_run  , only: dti, ramp, time
!      use river      , only: totq

      implicit none

      integer                   , intent(in) :: idx
      real(rk), dimension(im,jm), intent(in) :: ampe, amue, phae, phue

      integer i,j,k
      real(kind=rk) ga,u1, totq
!      real(kind=rk) hmax
      real(kind=rk) dum_area(jm), dum_flow(jm)
      real(kind=rk) sum_area, sum_flow, mean_uabe, uriv
      real(kind=rk) rdisp,rad_param,ome(6),ramt(6),phi0(6)   !fhx:tide


! dummy totq for compilation
      totq = 1.e6

      if(idx.eq.1) then

! external (2-D) elevation boundary conditions

         if(n_east.eq.-1) then
            do j=1,jm
               elf(im,j)=elf(imm1,j)
            end do
         endif

         do j=1,jm
            do i=1,im
               elf(i,j)=elf(i,j)*fsm(i,j)
            end do
         end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions


!       east

!     calculate net transport

         dum_area = 0d0
         dum_flow = 0d0


         if ( n_east == -1 ) then
         do j = 2,jmm1
          dum_area(j) = 0.25
     $       * ( h(im,j) + elf (im,j) + h(im-1,j) + elf(im-1,j) )
     $       * ( dy(im,j) + dy(im-1,j) ) * dum(im,j)

          dum_flow(j) = 0.25
     $       * ( h(im,j) + elf (im,j) + h(im-1,j) + elf(im-1,j) )
     $       * ( dy(im,j) + dy(im-1,j) ) * UA_bry%EST(1,j) * dum(im,j)
         end do
         endif


         call psum_mpi( dum_area, jm, sum_area )
         call psum_mpi( dum_flow, jm, sum_flow )


         if ( is_master ) then

!     mean u-velocity at the eastern boundary.

            mean_uabe = sum_flow / max( sum_area, small )

!     convert total river discharge [m3/s] to outgoing velocity [m/s]

            uriv = totq / max( sum_area, small )

         endif

!     broadcasting to all the nodes.

         call bcast0d_mpi( mean_uabe, 0 )
         call bcast0d_mpi( uriv, 0 )


!     adjust velocity at eastern boundary.

         if ( n_east == -1 ) then

         do j = 1,jm
           uaf(im,j) = ( UA_bry%EST(1,j) - mean_uabe+uriv ) * dum(im,j)
           vaf(im,j) = 0.0
         end do

!     tide  ................
!fhx:tide:!for tides coming in, use u = -sqrt(g/H)*elevation
          ramt = 0.
          ome  = 0.
          phi0 = 0.
       if(calc_tide) then

          rad_param=pi/180.0
          ome(1) = .08051140       !M2 in cycles per hour
          ome(2) = .08333333       !S2
          ome(5) = .04178075       !K1
          ome(6) = .03873066       !O1
          ome = ome*2.*pi*24.   !convert frequencies to radians/day
          rfe = 1.0
          ramt(1) = -0.5 !lyo:debug:1.0
          ramt(2) = 0.4
          ramt(5) = -0.35
          ramt(6) = 0.18
          phi0(1) = 0.   ! adjust phase
          phi0(2) = 106.*rad_param
          phi0(5) = 60.*rad_param
          phi0(6) = -40.*rad_param
       do j=1,jm
       rdisp=sqrt(grav/h(imm1,j))
       do k=1,ntide
       uaf(im,j)=uaf(im,j)+ramt(k)*(
     $   ramp*amue(j,k)*cos(ome(k)*time-phue(j,k)*rad_param+phi0(k))-
     $   rdisp*ampe(j,k)*cos(ome(k)*time-phae(j,k)*rad_param+phi0(k)))
       enddo

!      UAF(IM,J)=UAF(IM,J)+RAMT*( RAMP*AMUE(J,K)*COS(OME(K)*TIME-
!     1    PHUE(J,K)) - COVRHE(J)*AMPE(J,K)*COS(OME(K)*TIME-PHAE(J,K)) )

       enddo

       end if ! if(calc_tide)
!fhx:tide:end

         end if

!         if ( is_master ) print*, mean_uabe, uriv

         do j=1,jm
            do i=1,im
               uaf(i,j)=uaf(i,j)*dum(i,j)
               vaf(i,j)=vaf(i,j)*dvm(i,j)
            end do
         end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions


!     east
!     ratdiation boundary conditions.

         if(n_east.eq.-1) then

            hmax = 5300.

            do k=1,kbm1
               do j=2,jmm1
                  ga = sqrt( h(imm1,j) / hmax )
!                  ga = sqrt( h(im,j) / hmax )  !fhx:tide
                  uf(im,j,k)
     $                 = ga * ( 0.25 * u(imm1,j-1,k)
     $                 + 0.5 * u(imm1,j,k) + 0.25 * u(imm1,j+1,k) )
     $                 + ( 1.0 - ga ) * ( 0.25 * u(im,j-1,k)
     $                 + 0.5 * u(im,j,k) + 0.25 * u(im,j+1,k) )
               enddo
            enddo

         endif

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  uf(i,j,k)=uf(i,j,k)*dum(i,j)
                  vf(i,j,k)=vf(i,j,k)*dvm(i,j)
               end do
            end do
         end do

        return


      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)

!     east

         if(n_east.eq.-1) then

           do k=1,kbm1
             do j=1,jm
               u1=2.*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
               if(u1.le.0.) then
                 uf(im,j,k)=t(im,j,k)-u1*(T_bry%EST(1,j,k)-t(im,j,k))
                 vf(im,j,k)=s(im,j,k)-u1*(S_bry%EST(1,j,k)-s(im,j,k))
               else
                 uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
                 vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
!                     if(k.ne.1.and.k.ne.kbm1) then
!                        wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
!     $                       /((zz(k-1)-zz(k+1))*dt(imm1,j))
!                        uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
!                        vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
!                     endif
               end if
             enddo
           enddo

         end if


         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  uf(i,j,k)=uf(i,j,k)*fsm(i,j)
                  vf(i,j,k)=vf(i,j,k)*fsm(i,j)
               end do
            end do
         end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

!     east

         if(n_east.eq.-1) then
            do k=1,kb
               do j=1,jm
                  uf(im,j,k)=1.e-10
                  vf(im,j,k)=1.e-10
                  uf(1,j,k)=1.e-10
                  vf(1,j,k)=1.e-10
               enddo
            end do
         end if


        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      endif

      end ! subroutine bcond_nwatl
!
!______________________________________________________________________
!
      subroutine bcond0(idx)
!----------------------------------------------------------------------
!  Apply open boundary conditions.
!  Closed boundary conditions are automatically enabled through
! specification of the masks, dum, dvm and fsm, in which case the open
! boundary conditions, included below, will be overwritten.
!----------------------------------------------------------------------
! called by: [NO CALLS]
!______________________________________________________________________
!
      use bry
      use glob_const , only: grav, pi, rk, small
      use grid       , only: dum, dvm, dx, dy, fsm, h, zz
      use glob_domain
      use glob_ocean , only: dt, el, elf, q2, q2l, s, t
     &                     , u, uaf, uf, v, va, vaf, vf, w

      use model_run  , only: dte, dti, ramp, time
      implicit none

      integer idx
      integer i,j,k
      real(kind=rk) u1,wm,gae
      real(kind=rk) utide,etide

      if(idx.eq.1) then

! external (2-D) elevation boundary conditions
        do j=1,jm
          if(n_west.eq.-1) elf(1,j)=elf(2,j)
          if(n_east.eq.-1) elf(im,j)=elf(imm1,j)
        end do

        do i=1,im
          if(n_south.eq.-1) elf(i,1)=elf(i,2)
          if(n_north.eq.-1) elf(i,jm)=elf(i,jmm1)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions
        do j=2,jmm1
          ! west
          if(n_west.eq.-1) then
            uaf(2,j) = UA_bry%WST(1,j)
     &               - RFW*sqrt(grav/h(2,j))*(el(2,j)-EL_bry%WST(1,j))
            uaf(2,j) = ramp*uaf(2,j)
            uaf(1,j) = uaf(2,j)
            vaf(1,j) = 0.
          end if

          ! east
          if(n_east.eq.-1) then

        ! define an external tidal current
        utide=0.8*sin(2*pi*time*24./12.42)
        etide=1.1*cos(2*pi*time*24./12.42)

          uaf(im,j)=utide
     $               +rfe*sqrt(grav/h(imm1,j))
     $                         *(el(imm1,j)-etide)
          uaf(im,j)=ramp*uaf(im,j)
          gae=dte*sqrt(grav*h(imm1,j))/dx(imm1,j)
          vaf(im,j)=(va(im,j)+gae*vaf(imm1,j))/(1.+gae)
!          vaf(im,j)=0.e0


! ayumi 2010/4/7 ---------------------------------
          uaf(im,j) = UA_bry % EST(1,j)
          vaf(im,j) = 0.
!-------------------------------------------------

        end if

        end do

        do i=2,imm1
          ! south
          if ( n_south == -1 ) then
            vaf(i,2) = VA_bry%STH(i,1)
     &               - RFS*sqrt(grav/h(i,2))*(el(i,2)-EL_bry%STH(i,1))
            vaf(i,2) = ramp*vaf(i,2)
            vaf(i,1) = vaf(i,2)
            uaf(i,1) = 0.
          end if

          ! north
          if ( n_north == -1 ) then
            vaf(i,jm) = VA_bry%NTH(i,1)
     $         + RFN*sqrt(grav/h(i,jmm1))*(el(i,jmm1)-EL_bry%NTH(i,1))
            vaf(i,jm) = ramp*vaf(i,jm)
            uaf(i,jm) = 0.
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              u1=2.*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.) then
                uf(im,j,k)=t(im,j,k)-u1*(T_bry%EST(1,j,k)-t(im,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(S_bry%EST(1,j,k)-s(im,j,k))
              else
                uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(imm1,j))
                  uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                  vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
                endif
              end if
            end if

            ! west
            if(n_west.eq.-1) then
              u1=2.*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.) then
                uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-T_bry%WST(1,j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-T_bry%WST(1,j,k))
              else
                uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5*(w(2,j,k)+w(2,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(2,j))
                  uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                  vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
                end if
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.) then
                uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-T_bry%STH(i,1,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-S_bry%STH(i,1,k))
              else
                uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5*(w(i,2,k)+w(i,2,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,2))
                  uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                  vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
                end if
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.) then
                uf(i,jm,k)=t(i,jm,k)-u1*(T_bry%NTH(i,1,k)-t(i,jm,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(S_bry%NTH(i,1,k)-s(i,jm,k))
              else
                uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                  uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                  vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
                end if
              end if
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

        do k=1,kb
          do j=1,jm
            ! west
            if(n_west.eq.-1) then
              u1=2.*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.) then
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end if

            ! east
            if(n_east.eq.-1) then
              u1=2.*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.) then
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.) then
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.) then
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end if
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do

        return

      endif


      end ! subroutine bcond0
!
!______________________________________________________________________
!
      subroutine bcondorl(idx)
!----------------------------------------------------------------------
!  This is an optional subroutine replacing  bcond and using Orlanski's
! scheme (J. Comp. Phys. 21, 251-269, 1976), specialized for the
! seamount problem.
!----------------------------------------------------------------------
! called by: [NO CALLS]
!______________________________________________________________________
!
      use bry
      use glob_const , only: rk
      use glob_domain
      use grid       , only: dum, dvm, fsm
      use glob_ocean , only: elf, s, sb, t, tb
     &                     , u, ua, ub, uab, uaf, uf, vaf, vf, w

      implicit none

      integer idx
      real(kind=rk) cl,denom
      integer i,j,k


      if(idx.eq.1) then

! external (2-D) elevation boundary conditions
        do  j=1,jm
          if(n_west.eq.-1) elf(1,j)=elf(2,j)
          if(n_east.eq.-1) elf(im,j)=elf(imm1,j)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity  boundary conditions
        do j=2,jmm1
          ! east
          if(n_east.eq.-1) then
            denom=(uaf(im-1,j)+uab(im-1,j)-2.*ua(im-2,j))
            if(denom.eq.0.0)denom=0.01
            cl=(uab(im-1,j)-uaf(im-1,j))/denom
            if(cl.gt.1.) cl=1.
            if(cl.lt.0.) cl=0.
            uaf(im,j)=(uab(im,j)*(1.-cl)+2.*cl*ua(im-1,j))
     $                  /(1.+cl)
            vaf(im,j)=0.
          end if

          ! west
          if(n_west.eq.-1) then
            denom=(uaf(3,j)+uab(3,j)-2.*ua(4,j))
            if(denom.eq.0.0)denom=0.01
            cl=(uab(3,j)-uaf(3,j))/denom
            if(cl.gt.1.) cl=1.
            if(cl.lt.0.) cl=0.
            uaf(2,j)=(uab(2,j)*(1.-cl)+2.*cl*ua(3,j))
     $                 /(1.+cl)
            uaf(1,j)=uaf(2,j)
            vaf(1,j)=0.
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.*u(im-2,j,k))
              if(denom.eq.0.)denom=0.01
              cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              uf(im,j,k)=(ub(im,j,k)*(1.-cl)+2.*cl*u(im-1,j,k))
     $                    /(1.+cl)
              vf(im,j,k)=0.
            end if

            ! west
            if(n_west.eq.-1) then
              denom=(uf(3,j,k)+ub(3,j,k)-2.*u(4,j,k))
              if(denom.eq.0.)denom=0.01
              cl=(ub(3,j,k)-uf(3,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              uf(2,j,k)=(ub(2,j,k)*(1.-cl)+2.*cl*u(3,j,k))
     $                   /(1.+cl)
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=0.
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if (n_east == -1) then
              U_bry % EST(1,j,k) = ub(im,j,k)
              denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.*t(im-2,j,k))
              if(denom.eq.0.) denom=0.01
              cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              uf(im,j,k)=(tb(im,j,k)*(1.-cl)+2.*cl*t(im-1,j,k))
     $                    /(1.+cl)
              if(cl.eq.0..and.U_bry%EST(1,j,k)<=0.)
     &          uf(im,j,k)=T_bry%EST(1,j,k)

              denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.*s(im-2,j,k))
              if(denom.eq.0.) denom=0.01
              cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              vf(im,j,k)=(sb(im,j,k)*(1.-cl)+2.*cl*s(im-1,j,k))
     $                    /(1.+cl)
              if(cl.eq.0..and.U_bry%EST(1,j,k)<=0.)
     &          vf(im,j,k)=S_bry%EST(1,j,k)
            end if

            ! west
            if(n_west.eq.-1) then
              U_bry % WST(1,j,k) = ub(2,j,k)
              denom=(uf(2,j,k)+tb(2,j,k)-2.*t(3,j,k))
              if(denom.eq.0.) denom=0.01
              cl=(tb(2,j,k)-uf(2,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              uf(1,j,k)=(tb(1,j,k)*(1.-cl)+2.*cl*t(2,j,k))/(1.+cl)
              if(cl.eq.0..and.U_bry%WST(1,j,k)>=0.)
     &          uf(1,j,k)=T_bry%WST(1,j,k)

              denom=(vf(2,j,k)+sb(2,j,k)-2.*s(3,j,k))
              if(denom.eq.0.) denom=0.01
              cl=(sb(2,j,k)-vf(2,j,k))/denom
              if(cl.gt.1.) cl=1.
              if(cl.lt.0.) cl=0.
              vf(1,j,k)=(sb(1,j,k)*(1.-cl)+2.*cl*s(2,j,k))/(1.+cl)
              if(cl.eq.0..and.U_bry%WST(1,j,k)>=0.)
     &          vf(1,j,k)=S_bry%WST(1,j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions
        do k=1,kb
          do j=1,jm
            if(n_east.eq.-1) then
              uf(im,j,k)=1.e-10
              vf(im,j,k)=1.e-10
            end if
            if(n_west.eq.-1) then
              uf(1,j,k)=1.e-10
              vf(1,j,k)=1.e-10
            end if
          end do

          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      endif


      end ! subroutine bcondorl
!
!______________________________________________________________________
!
      subroutine dens(si,ti,rhoo)
!----------------------------------------------------------------------
!  Calculate (density-1000.)/rhoref.
! see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
!  NOTE: if pressure is not used in dens, buoyancy term (boygr) in
! subroutine profq must be changed (see note in subroutine profq.)
!----------------------------------------------------------------------
! called by: mode_internal [advance.f]
!            read_initial  [io.f90]    (x2)
!______________________________________________________________________
!
      use config     , only: sbias, tbias
      use glob_const , only: grav, rhoref, rk
      use glob_domain, only: im, jm, kb, kbm1
      use grid       , only: fsm, zz
      use glob_ocean , only: d

      implicit none

      real(rk), dimension(im,jm,kb), intent(in)  :: si, ti
      real(rk), dimension(im,jm,kb), intent(out) :: rhoo

      integer  i,j,k
      real(rk) cr,p,rhor,sr,tr,tr2,tr3,tr4


      do k=1,kbm1
        do j=1,jm
          do i=1,im

            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr

! approximate pressure in units of bars
            p=grav*rhoref*(-zz(k)*d(i,j))*1.e-5

            rhor=-0.157406+6.793952e-2*tr
     $           -9.095290e-3*tr2+1.001685e-4*tr3
     $           -1.120083e-6*tr4+6.536332e-9*tr4*tr

            rhor=rhor+(0.824493-4.0899e-3*tr
     $               +7.6438e-5*tr2-8.2467e-7*tr3
     $               +5.3875e-9*tr4)*sr
     $               +(-5.72466e-3+1.0227e-4*tr
     $               -1.6546e-6*tr2)*abs(sr)**1.5
     $               +4.8314e-4*sr*sr

            cr=1449.1+.0821*p+4.55*tr-.045*tr2
     $                 +1.34*(sr-35.)
            rhor=rhor+1.e5*p/(cr*cr)*(1.-2.*p/(cr*cr))

            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)

          end do
        end do
      end do


      end ! subroutine dens
!
!______________________________________________________________________
!
      subroutine profq(sm,sh,dh,cc)
!----------------------------------------------------------------------
!  Solve for q2 (twice the turbulent kinetic energy),
! q2l (q2 x turbulent length scale), km (vertical kinematic viscosity)
! and kh (vertical kinematic diffusivity), using a simplified version 
! of the level 2 1/2 model of Mellor and Yamada (1982).
!  In this version, the Craig-Banner sub-model whereby breaking
! wave tke is injected into the surface is included. However, we use an
! analytical solution to the near surface tke equation to solve for q2
! at the surface giving the same result as C-B diffusion. The new
! scheme is simpler and more robust than the latter scheme.
!----------------------------------------------------------------------
! called by: mode_internal  [advance.f]
!
! calls    : exchange2d_mpi [parallel_mpi.f]
!            exchange3d_mpi [parallel_mpi.f]
!______________________________________________________________________
!
      use air        , only: wusurf, wvsurf
      use config     , only: sbias, tbias, umol
      use glob_const , only: grav, kappa, rhoref, rk, small
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1!, kbm2
     &                     , n_east, n_north, n_south, n_west
      use grid       , only: dzz, dz, fsm, h, z, zz
      use glob_ocean , only: a, c, dtef, ee, etf, gg, kh, km, kq, l
     &                     , q2, q2b, q2lb, rho, s, t
     &                     , u, uf, v, vf, wubot, wvbot
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(out) :: cc, sh, sm
      real(rk), dimension(im,jm   ), intent(out) :: dh

      real(rk), dimension(im,jm,kb) :: boygr, gh, prod, stf
      real(rk), dimension(im,jm   ) :: l0
      real(rk)    a1,    a2,    b1,    b2,    c1
     &       , coef1, coef2, coef3, coef4, coef5
     &       ,const1,    e1,    e2,   ghc
     &       ,     p,   sef,    sp,    tp
     &       ,cbcnst, surfl,  shiw
     &       , utau2
      integer i, j, k, ki

      data a1,b1,a2,b2,c1/0.92,16.6,0.74,10.1,0.08/
      data e1/1.8/,e2/1.33/
      data sef/1./
      data cbcnst/100./surfl/2.e5/shiw/0.0/


      dh = h + etf

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.*umol)*.5
     $                                 /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.*umol)*.5
     $                               /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

! the following section solves the equation:
!     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

      ! surface and bottom boundary conditions
      const1=(16.6**(2./3.))*sef

      ! initialize fields that are not calculated on all boundaries
      ! but are later used there
      do i=1,im
        ee(i,jm,1)=0.
        gg(i,jm,1)=0.
        l0(i,jm)=0.
      end do
      do j=1,jm
        ee(im,j,1)=0.
        gg(im,j,1)=0.
        l0(im,j)=0.
      end do
      do i=1,im
        do j=1,jm
          do k=2,kbm1
            prod(i,j,k)=0.
          end do
        end do
      end do

      do j=1,jmm1
        do i=1,imm1
          utau2=sqrt((.5*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                           +(.5*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
          ! wave breaking energy- a variant of Craig & Banner (1994)
          ! see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2
          ! surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2/grav
          uf(i,j,kb)=sqrt((.5*(wubot(i,j)+wubot(i+1,j)))**2
     $                      +(.5*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
      call exchange2d_mpi(ee(:,:,1),im,jm)
      call exchange2d_mpi(gg(:,:,1),im,jm)
      call exchange2d_mpi(l0,im,jm)

      ! calculate speed of sound squared
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
            ! calculate pressure in units of decibars
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-4
            cc(i,j,k)=1449.1+.00821*p+4.55*tp -.045*tp**2
     $                                               +1.34*(sp-35.0)
            cc(i,j,k)=cc(i,j,k)/sqrt((1.-.01642*p/cc(i,j,k))
     $                                    *(1.-0.40*p/cc(i,j,k)**2))
          end do
        end do
      end do

      ! calculate buoyancy gradient
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                                               /(dzz(k-1)* h(i,j))
         ! note: comment out next line if dens does not include pressure
     $                    +(grav**2)*2./(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
!            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
! ayumi 2010/4/12
            l(i,j,k)=abs(q2lb(i,j,k)/(q2b(i,j,k)+small))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
!            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
! ayumi 2010/4/12
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)
     $           /(q2b(i,j,k)+small)
            gh(i,j,k)=min(gh(i,j,k),.028_rk)
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.
          gh(i,j,1)=0.
          gh(i,j,kb)=0.
        end do
      end do

! calculate production of turbulent kinetic energy:
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25*sef
     $                *((u(i,j,k)-u(i,j,k-1)+u(i+1,j,k)-u(i+1,j,k-1))**2
     $                +(v(i,j,k)-v(i,j,k-1)+v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                                            /(dzz(k-1)*dh(i,j))**2
                                  ! add shear due to internal wave field
     $                                      -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do

!   NOTE: Richardson # dep. dissipation correction (Mellor, 2001;
! Ezer, 2000), depends on ghc the critical number (empirical -6 to -2)
! to increase mixing.
      ghc=-6.0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.
            ! it is unclear yet if diss. corr. is needed when surf.
            ! waves are included.
!           if(gh(i,j,k).lt.0.e0)
!    $                     stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
!           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                                              /(b1*l(i,j,k)+small)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1./(a(i,j,k)+c(i,j,k)*(1.-ee(i,j,k-1))
     $                                    -(2.*dti2*dtef(i,j,k)+1.))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                                             -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do

! the following section solves the equation:
!     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
      do j=1,jm
        do i=1,im
          vf(i,j,1)=0.
          vf(i,j,kb)=0.
          ee(i,j,2)=0.
          gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
          vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
      end do
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)*(1.+e2*((1./abs(z(k)-z(1))
     $              +1./abs(z(k)-z(kb)))*l(i,j,k)/(dh(i,j)*kappa))**2)
          end do
        end do
      end do
      do k=3,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1./(a(i,j,k)+c(i,j,k)*(1.-ee(i,j,k-1))
     $                                         -(dti2*dtef(i,j,k)+1.))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                        +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do

      ! the following is to counter the problem of the ratio of two small
      ! numbers (l = q2l/q2) or one number becoming negative. Two
      ! options are included below. In this application, the second
      ! option, l was less noisy when uf or vf is small
      do k=2,kbm1
        do j=1,jm
          do i=1,im
!           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
!             uf(i,j,k)=small
!             vf(i,j,k)=0.1*dt(i,j)*small
!           end if
          uf(i,j,k)=abs(uf(i,j,k))
          vf(i,j,k)=abs(vf(i,j,k))
          end do
        end do
      end do

! the following section solves for km and kh
      coef4=18.*a1*a1+9.*a1*a2
      coef5=9.*a1*a2

      ! note that sm and sh limit to infinity when gh approaches 0.0288
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.-6.*a1/b1*stf(i,j,k))
            coef2=3.*a2*b2/stf(i,j,k)+18.*a1*a2
            coef3=a1*(1.-3.*c1-6.*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.-coef5*gh(i,j,k))
          end do
        end do
      end do

      ! there are 2 options for kq which, unlike km and kh, was not
      ! derived by Mellor and Yamada but was purely empirical based on
      ! neutral boundary layer data. The choice is whether or not it
      ! should be subject to the stability factor, sh. Generally,
      ! there is not a great difference in output
      do k=1,kb
        do j=1,jm
          do i=1,im
            prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(prod(i,j,k)*.41*sh(i,j,k)+kq(i,j,k))*.5
!            kq(i,j,k)=(prod(i,j,k)*.20+kq(i,j,k))*.5
            km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5
            kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5
          end do
        end do
      end do
      call exchange3d_mpi(km,im,jm,kb)
      call exchange3d_mpi(kh,im,jm,kb)

      ! cosmetics: make boundr. values as interior (even if not used,
      ! printout may show strange values)
      do k=1,kb
        do i=1,im
          if(n_north.eq.-1) then
            km(i,jm,k)=km(i,jmm1,k)*fsm(i,jm)
            kh(i,jm,k)=kh(i,jmm1,k)*fsm(i,jm)
          end if
          if(n_south.eq.-1) then
            km(i,1,k)=km(i,2,k)*fsm(i,1)
            kh(i,1,k)=kh(i,2,k)*fsm(i,1)
          end if
        end do
        do j=1,jm
          if(n_east.eq.-1) then
            km(im,j,k)=km(imm1,j,k)*fsm(im,j)
            kh(im,j,k)=kh(imm1,j,k)*fsm(im,j)
          end if
          if(n_west.eq.-1) then
            km(1,j,k)=km(2,j,k)*fsm(1,j)
            kh(1,j,k)=kh(2,j,k)*fsm(1,j)
          end if
        end do
      end do


      end ! subroutine profq
!
!______________________________________________________________________
!
      subroutine proft(f,wfsurf,fsurf,nbc,dh)
!----------------------------------------------------------------------
!  Solves for vertical diffusion of temperature and salinity using
! method described by Richmeyer and Morton (1967).
!  NOTE: wfsurf and swrad are negative values when water column is
! warming or salt is being added.
!----------------------------------------------------------------------
! called by: mode_internal [advance.f] (x2)
!______________________________________________________________________
!
      use air        , only: swrad
      use config     , only: ntp, umol
      use glob_const , only: rk
      use glob_domain, only: im, jm, kb, kbm1, kbm2, i_global, j_global
      use grid       , only: dz, dzz, h, z
      use glob_ocean , only: a, c, ee, etf, gg, kh
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(inout) :: f
      real(rk), dimension(im,jm   ), intent(in   ) :: fsurf, wfsurf
      real(rk), dimension(im,jm   ), intent(  out) :: dh
      integer                      , intent(in   ) :: nbc

      integer  i, j, k, ki
      real(rk) rad(im,jm,kb), r(5), ad1(5), ad2(5)

! irradiance parameters after Paulson and Simpson (1977)
!       ntp               1      2       3       4       5
!   Jerlov type           i      ia      ib      ii     iii
      data r   /         .58,    .62,    .67,    .77,    .78 /
      data ad1 /         .35,    .60,   1.  ,   1.5 ,   1.4  /
      data ad2 /       23.  ,  20.  ,  17.  ,  14.  ,   7.9  /


! surface boundary condition:
!       nbc   prescribed    prescribed   short wave
!             temperature      flux      penetration
!             or salinity               (temperature
!                                           only)
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! note that only 1 and 3 are allowed for salinity.

! the following section solves the equation
!     dti2*(kh*f')'-f=-fb
      dh = h + etf

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do


! calculate penetrative radiation. At the bottom any unattenuated
! radiation is deposited in the bottom layer
      rad = 0.

      if ( nbc == 2 .or. nbc == 4 ) then
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)*
     &                (     r(ntp) *exp(z(k)*dh(i,j)/ad1(ntp))
     &                 +(1.-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
      endif


      if ( nbc == 1 ) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.)
          end do
        end do


      elseif ( nbc == 2 ) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.)
          end do
        end do


      elseif ( nbc == 3 .or. nbc == 4) then

        ee(:,:,1) = 0.
        gg(:,:,1) = fsurf

      end if

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1./(a(i,j,k)+c(i,j,k)*(1.-ee(i,j,k-1))-1.)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do


! bottom adiabatic boundary condition
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.-ee(i,j,kbm2))-1.)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
!              if ( f(i,j,ki) == f(i,j,ki)+1 ) then
!                print *, "[ PROFT ]", i_global(i), j_global(j)
!                print *, " ee: ", ee(i,j,ki)
!                print *, " f+: ", f(i,j,ki+1)
!                print *, " gg: ", gg(i,j,ki)
!                print *, " sw: ", swrad(i,j)
!                print *, " wf: ", wfsurf(i,j)
!                stop
!              end if
          end do
        end do
      end do


      end ! subroutine proft
!
!______________________________________________________________________
!
      subroutine profu
!----------------------------------------------------------------------
!  Solves for vertical diffusion of x-momentum using method described
! by Richmeyer and Morton (1967).
!  NOTE: wusurf has the opposite sign to the wind speed.
!----------------------------------------------------------------------
! called by: mode_internal  [advance.f]
!
! calls    : exchange2d_mpi [parallel_mpi.f]
!______________________________________________________________________
!
      use air        , only: wusurf
      use config     , only: umol
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1, kbm2
      use grid       , only: dum, dz, dzz, h
      use glob_ocean , only: a, c, cbc, ee, etf, gg, km, tps
     &                     , ub, uf, vb, wubot
      use model_run  , only: dti2

      implicit none

      integer  i,j,k,ki
      real(rk) dh(im,jm)


! the following section solves the equation
!   dti2*(km*u')'-u=-ub
      dh = 1.

      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5
        end do
      end do

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1./(a(i,j,k)+c(i,j,k)*(1.-ee(i,j,k-1))-1.)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.25*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.
     $                    -(ee(i,j,kbm2)-1.)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wubot,im,jm)


      end ! subroutine profu
!
!______________________________________________________________________
!
      subroutine profv
!----------------------------------------------------------------------
!  Solves for vertical diffusion of x-momentum using method described
! by Richmeyer and Morton (1967).
!  NOTE: wvsurf has the opposite sign to the wind speed.
!----------------------------------------------------------------------
! called by: mode_internal  [advance.f]
!
! calls    : exchange2d_mpi [parallel_mpi.f]
!______________________________________________________________________
!
      use air        , only: wvsurf
      use config     , only: umol
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1, kbm2
      use grid       , only: dvm, dz, dzz, h
      use glob_ocean , only: a, c, cbc, ee, etf, gg, km, tps
     &                     , vb, ub, vf, wvbot
      use model_run  , only: dti2

      implicit none

      integer  i,j,k,ki
      real(rk) dh(im,jm)


! the following section solves the equation
!     dti2*(km*u')'-u=-ub
      dh = 1.

      do j=2,jm
        do i=2,im
          dh(i,j)=.5*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1./(a(i,j,k)+c(i,j,k)*(1.-ee(i,j,k-1))-1.)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.25*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.
     $                    -(ee(i,j,kbm2)-1.)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wvbot,im,jm)


      end ! subroutine profv
!
!______________________________________________________________________
!
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
!----------------------------------------------------------------------
!  Calculate the antidiffusive velocity used to reduce the numerical
! diffusion associated with the upstream differencing scheme.
!  This is based on a subroutine of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy).
!----------------------------------------------------------------------
! called by: advt2 [solver.f]
!______________________________________________________________________
!
      use config     , only: sw
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: aru, arv, dzz, fsm
      use glob_ocean , only: dt
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(inout) :: ff, xmassflux,
     &                                                ymassflux, zwflux

      integer                i,j,k
      real(rk)               mol,abs_1,abs_2
      real(rk)               udx,u2dt,vdy,v2dt,wdz,w2dt
      real(rk), parameter :: value_min = 1.e-9
     &                     , epsilon   = 1.e-14


! apply temperature and salinity mask
      do k=1,kb
        ff(:,:,k) = ff(:,:,k)*fsm
      end do

! recalculate mass fluxes with antidiffusion velocity
!      rewind(40+my_task)
!      write(40+my_task,*) dt, xmassflux
!      if ( is_master ) print *, "=== ",iext," ==="
!      print '(i1,x,a2,x,2(e14.7,x,2(i3,";")))', my_task, "dt"
!     &       , minval(dt, fsm/=0.), minloc(dt, fsm/=0.)
!     &       , maxval(dt), maxloc(dt)
!      print '(i1,x,a2,x,2(e15.8,x,3(i3,";")))', my_task, "xf"
!     &       , minval(xmassflux), minloc(xmassflux)
!     &       , maxval(xmassflux), maxloc(xmassflux)
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.
            else
              udx=abs(xmassflux(i,j,k))
!              write(40+my_task,*) i, j, k, aru(i,j), dt(i-1,j)
!     &             , dt(i,j), xmassflux(i,j,k), dti2, aru(i,j)
!              write(50+my_task,*)
!     &             dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.
!              write(60+my_task,*) aru(i,j)*(dt(i-1,j)+dt(i,j))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.
            end if
          end do
        end do
      end do

!      print '(i1,x,a2,x,2(e15.8,x,3(i3,";")))', my_task, "yf"
!     &       , minval(ymassflux), minloc(ymassflux)
!     &       , maxval(ymassflux), maxloc(ymassflux)
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.
            end if
          end do
        end do
      end do

!      print '(i1,x,a2,x,2(e15.8,x,3(i3,";")))', my_task, "zf"
!     &       , minval(zwflux), minloc(zwflux)
!     &       , maxval(zwflux), maxloc(zwflux)
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.
            end if
          end do
        end do
      end do


      end ! subroutine smol_adif
!
!______________________________________________________________________
!
      subroutine smol_adifC(xmassflux,ymassflux,cf)
!----------------------------------------------------------------------
! This is a 2D modification of `smol_adif` for sea-ice concentration.
!----------------------------------------------------------------------
! called by: advtC [solver.f]
!______________________________________________________________________
!
      use config     , only: sw
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb
      use grid       , only: aru, arv, fsm
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm), intent(inout) :: cf
     &                                           , xmassflux, ymassflux

      integer                i,j
      real(rk)               mol,abs_1,abs_2
      real(rk)               udx,u2dt,vdy,v2dt
      real(rk), parameter :: value_min = 1.e-9
     &                     , epsilon   = 1.e-14


! apply temperature and salinity mask
      cf = cf*fsm

! recalculate mass fluxes with antidiffusion velocity
      do j=2,jmm1
        do i=2,im
          if(cf(i,j).lt.value_min.or.
     $       cf(i-1,j).lt.value_min) then
            xmassflux(i,j) = 0.
          else
            udx=abs(xmassflux(i,j))
            u2dt=dti2*xmassflux(i,j)*xmassflux(i,j)*2.
     $          /aru(i,j)
            mol=(cf(i,j)-cf(i-1,j))
     $          /(cf(i-1,j)+cf(i,j)+epsilon)
            xmassflux(i,j)=(udx-u2dt)*mol*sw
            abs_1=abs(udx)
            abs_2=abs(u2dt)
            if(abs_1.lt.abs_2) xmassflux(i,j)=0.
          end if
        end do
      end do

      do j=2,jm
        do i=2,imm1
          if(cf(i,j).lt.value_min.or.
     $       cf(i,j-1).lt.value_min) then
            ymassflux(i,j)=0.
          else
            vdy=abs(ymassflux(i,j))
            v2dt=dti2*ymassflux(i,j)*ymassflux(i,j)*2.
     $          /arv(i,j)
            mol=(cf(i,j)-cf(i,j-1))
     $          /(cf(i,j-1)+cf(i,j)+epsilon)
            ymassflux(i,j)=(vdy-v2dt)*mol*sw
            abs_1=abs(vdy)
            abs_2=abs(v2dt)
            if(abs_1.lt.abs_2) ymassflux(i,j)=0.
          end if
        end do
      end do


      end ! subroutine smol_adifC
!
!______________________________________________________________________
!
      subroutine vertvl(xflux,yflux)
!----------------------------------------------------------------------
!  Calculates vertical velocity.
!----------------------------------------------------------------------
! called by: mode_internal [advance.f]
!______________________________________________________________________
!
      use air        , only: vfluxb, vfluxf
      use glob_const , only: rk
      use glob_domain, only: im, imm1, jm, jmm1, kb, kbm1
      use grid       , only: dx, dy, dz
      use glob_ocean , only: dt, etb, etf, u, v, w
      use model_run  , only: dti2

      implicit none

      real(rk), dimension(im,jm,kb), intent(out) :: xflux, yflux

      integer i,j,k


! reestablish boundary conditions
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do

! note: if one wishes to include freshwater flux, the surface velocity
! should be set to vflux(i,j). See also change made to 2-D volume
! conservation equation which calculates elf
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do


      end ! subroutine vertvl
!
!lyo:20110315:botwavedrag:add subr.botwavedrag & function fsinhinv
!______________________________________________________________________
!
       subroutine botwavedrag (im,jm,fsm,wusrf,wvsrf,kp,
     &                      wubot,wvbot,d,zzkbm1,z0b,cbcmin,cbcmax,cbc)
!----------------------------------------------------------------------!
!     This s.r.botwavedrag was modified from G.Mellor's s.r. botdrag   !
!     in /home/glm/MMS/pom08.f so that it now does not require calc.   !
!     wave inputs.  Given wind-stress and estimates of the peak wave   !
!     number "kp", the s.r. uses the Toba's 3/2 Law in Nielson's       !
!     formula to increase the bottom friction coefficient cbc due to   !
!     waves. This s.r. is intended for use in tropical cyclone         !
!     (hurricanes, typhoons) wind conditions, when kp ~ 2*pi/(200m)    !
!     e.g. Li et al (2002), and kp is input as a constant if wave is   !
!     not calculated in the POM model.                                 !
!                                                                      !
!     Ref: Jones & Toba, 2001: Wind Stress over the Ocean. Cambridge U !
!          Press, 307pp. See pages 46 & 116.                           !
!          Sig.Wave.Ht = Hs (=4*sqrt(ent/grav)) is used to express     !
!          Nielson's formula in terms of Hs.  Then Toba's Law eqn.4.39 !
!          is used to relate Hs to ts (=Significant wave period). Then !
!          eqn.(4.34) is used to relate omep (peak freq) to ts.        !
!          Finally, dispersion relation omep = sqrt(g*kp*tanh(kp*d))   !
!          is used to express uboscil = fn(kp,WindFricVel)             !
!                                                                      !
!          Nielson, P., 1992: Coastal bottom boundary layers and       !
!          sediment transport. World Scientific.                       !
!                                                                      !
!          Li et al, 2002: Observation of hurricane-generated ocean    !
!          swell refraction at the Gulf Stream north wall with the     !
!          RADARSAT-1 synthetic aperture radar. IEEE Transc. Geosci &  !
!          remote sensing, 40, 2131-2142.                              !
!                                                                      !
!     subr. is from /wrk/newshoni/lyo/apom/profs/workversions/         !
!                    lyo_botdrag2.f subr. botdrag2                     !
!                                                                      !
!           L.Oey (Jul/30/2010)(Mar/15/2011)                           !
!----------------------------------------------------------------------!
! Provides wave influenced bottom fricion. Derived from
! Nielson, P., 1992, Coastal bottom boundary layers and sediment
! transport. World Scientific.
! This aspect needs more basic research.
!
!  Input:
!    fsm(I,j) = 1 for water cell, = 0 for land cell
!    wusrf(i,j),wvsrf(i,j) = kinematic wind stress vector
!    kp      = peak frequency wave number (m-1), assumed =
!              the Significant wavenumber, a constant
!    wubot(i,j),wvbot(i,j) = bottom stress vector as determned
!      by bottom current (at zz(kb-1)) and cbc found here
!    zzkbm1 = zz(kb-1) (non-d.)
!    d(i,j) = water column depth (m)
!    z0b = roughness parameter (m)
!    cbcmin,cbcmax = limits on cbs (non-d.)
!
!  Output:
!    cbc(i,j) = bottom drag coefficient (non-d.)
!
!  Locals:
!    ts      = significant wave period (s), from kp via dispers relatn.
!----------------------------------------------------------------------
! called by: wind_main [wind.f]
!______________________________________________________________________
!
      use config    , only: rk
      use glob_const, only: KAPPA, GRAV!, PI

      implicit none

      integer                   , intent(in   ) :: im, jm
      real(rk)                  , intent(in   ) :: cbcmax, cbcmin
     &                                           , kp    , z0b
     &                                           , zzkbm1
      real(rk), dimension(im,jm), intent(in   ) :: d    , fsm
     &                                           , wubot, wvbot
     &                                           , wusrf, wvsrf
      real(rk), dimension(im,jm), intent(  out) :: cbc

      real(rk)               const   , uboscil , utau2
     &                     , utauwind, fsinhinv, z0a

      integer                i, j
      real(rk), parameter :: btoba    =  .062
     &                     , utau2min = 1.e-5


      const = 1.004849  !=btoba*sqrt(grav)*((pi/1.05)**(3/2))

      do j=1,jm
        do i=1,im
          if (fsm(i,j)==1.) then
!           uboscil=cp(i,j)*kp(i,j)*sqrt(2.*ent(i,j)/grav) !glm's orig
!    &               *fsinhinv(kp(i,j)*d(i,j))             !Nielson-fml
            utauwind = ( (wusrf(i,j)**2+wvsrf(i,j)**2)
     &                  /(grav*kp*tanh(kp*d(i,j)))    )**0.25
            uboscil  = const*utauwind*fsinhinv(kp*d(i,j))
            utau2    = sqrt(wubot(i,j)**2+wvbot(i,j)**2)+utau2min
            z0a      = z0b*(1.+0.05*uboscil**2/utau2)
            cbc(i,j) = ( kappa / log( 1.+(1.0+zzkbm1)*d(i,j)/z0a ) )**2
!     WRITE(10,'(''(1+zzkbm1)*d(i,j)='',1p1e13.5)')(1.0+zzkbm1)*d(i,j)
!     WRITE(10,'('' cbc = '',1p1e13.5)') cbc(i,j)
!     WRITE(10,'('' log = '',1p1e13.5)')log(1.+(1.0+zzkbm1)*d(i,j)/z0a)
            cbc(i,j) = max(cbcmin,cbc(i,j))
!     WRITE(10,'('' cbc = '',1p1e13.5)') cbc(i,j)
!
!     If the following is invoked, then it is probable that the wrong
!     choice of z0b or vertical spacing has been made:
!
            cbc(i,j) = min(cbcmax,cbc(i,j))
          end if
        end do
      end do

!
!     stop

!     WRITE(10,'('' z0a,z0b = '',1p2e13.5)') z0a,z0b
!
!     call exchange2d_mpi(cbc,im,jm) !lyo:don;t need?


      end ! subroutine botwavedrag
!
!lyo:20110315:botwavedrag:add function fsinhinv
!______________________________________________________________________
!
      function fsinhinv(x)
!----------------------------------------------------------------------
! called by: botwavedrag [solver.f]
!----------------------------------------------------------------------
      use config, only: rk

      implicit none

      real(rk) fsinhinv

      real(rk), intent(in) :: x

      integer                    ixd
      real(rk), dimension(15) :: fdat, xd

      data xd/0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7./
      data fdat/1.919,1.919,0.851,0.470,0.276,0.165,0.100,0.060,
     &          0.037,0.022,0.013,0.008,0.005,0.003,0.002/


      fsinhinv = 0.
      if (x < 7.) then
        ixd = int(2.*x)+1
        fsinhinv = fdat(ixd)+(fdat(ixd+1)-fdat(ixd))*(x-xd(ixd))*2.
      end if


      end ! function fsinhinv