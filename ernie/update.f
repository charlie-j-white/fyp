


     
      subroutine update(N, mach_inf, pres_rat, pres, u1, u2, u5, S, Sx,
     &                          dt, dx, cfl)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: mach_inf, pres_rat, dx, cfl
      double precision, intent(out) :: dt
      double precision, dimension(0:N+1), intent(out) :: pres, u1,u2,u5 
      double precision, dimension(0:N+1), intent(in) :: S, Sx

      integer :: ii
      double precision :: gam, pinf, pstag, rinf, rstag, ubc, mbc, pbc,
     &                  rbc, tau, cmax
      double precision, dimension(0:N+1) :: u1b, u2b, u5b,
     &                          f1, f2, f5, f1b, f2b, f5b, 
     &                          q1, q2, q5, q1b, q2b, q5b,
     &                          c 


      gam = 1.4d0





!
!
! calculate time step
!
!-----------------------------------------------------------------------
      do ii = 0,N+1
      c(ii) = dabs(dsqrt(gam*pres(ii)/u1(ii))) + dabs(u2(ii)/u1(ii))
      end do

      cmax = maxval(c)
      dt = cfl*dx/cmax
      tau = dt/dx
      !dt = 0.001



!      print*, dt






! initialise variables:
!       calculate f vector 0-->N+1
!       calculate q vector 0-->N+1
!-----------------------------------------------------------------------

      do ii = 0,N+1
             call fvector(u1(ii), u2(ii), u5(ii), 
     &                    f1(ii), f2(ii), f5(ii))
      end do
      
      do ii = 0,N+1
             call qvector(u1(ii), u2(ii), u5(ii), 
     &                    q1(ii), q2(ii), q5(ii), S(ii), Sx(ii))
      end do



! calculate predictor terms:
!       calculate u vector predictor terms 0-->N
!       calculate f vector predictor terms 0-->N
!       calculate q vector predictor terms 0-->N
!-----------------------------------------------------------------------

      do ii = 0,N
             u1b(ii) = u1(ii) - tau*(f1(ii+1)-f1(ii)) + dt*q1(ii)
             u2b(ii) = u2(ii) - tau*(f2(ii+1)-f2(ii)) + dt*q2(ii)
             u5b(ii) = u5(ii) - tau*(f5(ii+1)-f5(ii)) + dt*q5(ii)
      end do

      do ii = 0,N
             call fvector(u1b(ii), u2b(ii), u5b(ii), 
     &                    f1b(ii), f2b(ii), f5b(ii))
      end do

      do ii = 0,N
             call qvector(u1b(ii), u2b(ii), u5b(ii), 
     &                    q1b(ii), q2b(ii), q5b(ii), S(ii), Sx(ii))
      end do



! update:
!       calculate fully updated u vector 1-->N
!-----------------------------------------------------------------------

      do ii = 1,N
             u1(ii) = 0.5*(u1(ii)+u1b(ii))-0.5*tau*(f1b(ii)-f1b(ii-1))
     &                    + dt*q1b(ii)/2
             u2(ii) = 0.5*(u2(ii)+u2b(ii))-0.5*tau*(f2b(ii)-f2b(ii-1))
     &                    + dt*q2b(ii)/2
             u5(ii) = 0.5*(u5(ii)+u5b(ii))-0.5*tau*(f5b(ii)-f5b(ii-1))
     &                    + dt*q5b(ii)/2
             pres(ii) = (gam-1)*(u5(ii)-0.5d0*u2(ii)**2/u1(ii))
      end do



      end subroutine update
