!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     
      subroutine update(N, mach_inf, pres_rat, pres, u1, u2, u5)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: mach_inf, pres_rat
      double precision, dimension(0:N+1), intent(out) :: pres, u1,u2,u5 

      integer :: ii
      double precision :: gam, pinf, pstag, rinf, rstag, ubc, mbc, pbc,
     &                  rbc, tau
      double precision, dimension(0:N+1) :: u1b, u2b, u5b,
     &                          f1, f2, f5, f1b, f2b, f5b 


      

      tau = 0.5

! calculate f vector 0-->N+1
      do ii = 0,N+1
             call fvector(u1(ii), u2(ii), u5(ii), 
     &                    f1(ii), f2(ii), f5(ii))
      end do


! calculate u vector predictor terms  0-->N
      do ii = 0,N
             u1b(ii) = u1(ii) - tau*(f1(ii+1)-f1(ii))
             u2b(ii) = u2(ii) - tau*(f2(ii+1)-f2(ii))
             u5b(ii) = u5(ii) - tau*(f5(ii+1)-f5(ii))
      end do


! calculate f vector predictor terms 0-->N
      do ii = 0,N
             call fvector(u1b(ii), u2b(ii), u5b(ii), 
     &                    f1b(ii), f2b(ii), f5b(ii))
      end do


! calculate fully updated u vector 1-->N
      do ii = 1,N
             u1(ii) = 0.5*(u1(ii)+u1b(ii))-0.5*tau*(f1b(ii)-f1b(ii-1))
             u2(ii) = 0.5*(u2(ii)+u2b(ii))-0.5*tau*(f2b(ii)-f2b(ii-1))
             u5(ii) = 0.5*(u5(ii)+u5b(ii))-0.5*tau*(f5b(ii)-f5b(ii-1))
      end do



      end subroutine update
