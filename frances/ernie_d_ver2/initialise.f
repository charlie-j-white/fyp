



     
      subroutine initialise(N, mach_inf, pres_rat, pres, u1, u2, u5)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: mach_inf, pres_rat
      double precision, dimension(0:N+1), intent(out) :: pres, u1,u2,u5 

      integer :: ii
      double precision :: gam, pinf, pout, ptot, rinf, rout, vel


! calculate inflow and outflow values for pressure and density
      gam = 1.4
      pinf = 1.0/gam
      pout = pinf*pres_rat
      rinf = 1.0
      rout = rinf*pres_rat**(1/gam)
      ptot = pinf + 0.5*rinf*mach_inf**2

! linearly interpolate these values for the length of the nozzle
      do ii = 1,N

            pres(ii) = pinf + (pout-pinf)*(ii-1)/(N-1)
            u1(ii) = rinf + (rout-rinf)*(ii-1)/(N-1)
            vel = ((ptot-pres(ii))/(0.5*u1(ii)))**(0.5)
            u2(ii) = u1(ii)*vel
            u5(ii) = pres(ii)/(gam-1) + 0.5*u1(ii)*vel**2

      end do





      end subroutine initialise
