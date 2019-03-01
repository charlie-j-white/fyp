



     
      subroutine boundaries(N, mach_inf, pres_rat, pres, u1, u2, u5)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: mach_inf, pres_rat
      double precision, dimension(0:N+1), intent(out) :: pres, u1,u2,u5 

      integer :: ii
      double precision :: gam, pinf, pstag, rinf, rstag, ubc, mbc, pbc,
     &                  rbc


! useful values near the inflow
      gam = 1.4
      pinf = 1.0/gam
      pstag = pinf*(1+0.5*(gam-1)*mach_inf**2)**(gam/(gam-1))
      rinf = 1.0
      rstag = rinf*(1+0.5*(gam-1)*mach_inf**2)**(1.0/(gam-1))

      ubc = u2(1)/u1(1)
      mbc = ubc*ubc/(1 + 0.5*(gam-1)*mach_inf**2 - 0.5*(gam-1)*ubc**2)
      pbc = pstag/(1+0.5*(gam-1)*mbc**2)**(gam/(gam-1))
      rbc = rstag/(1+0.5*(gam-1)*mbc**2)**(1.0/(gam-1))

! halo cells at the inflow
      pres(0) = pbc
      u1(0) = rbc
      u2(0) = rbc*ubc
      u5(0) = pbc/(gam-1) + 0.5*rbc*ubc**2

! halo cells at the outflow
      pres(N+1) = pinf*pres_rat
      u1(N+1) = U1(N)
      u2(N+1) = U2(N)
      u5(N+1) = pres(N+1)/(gam-1) + 0.5*(u2(N+1)**2)/u1(N+1)



      end subroutine boundaries
