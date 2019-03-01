



     
      subroutine fvector(u1, u2, u5, f1, f2, f5)

      implicit none

      double precision, intent(in) :: u1, u2, u5
      double precision, intent(out) :: f1, f2, f5
      
      double precision :: pres, vel, gam

      gam = 1.40
      ! area terms S inlculded in the u variable, or cancelled


      

! clculate pressure
      pres = (gam-1)*(u5 - 0.5*(u2**2)/u1)

! calculate velocity
      vel = u2/u1


! calculate f vector
      f1 = vel*u1
      f2 = vel*u2 + pres
      f5 = vel*u5 + pres*vel


      end subroutine fvector



     
      subroutine jcost(N, J, pres, dx)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: dx
      double precision, intent(out) :: J
      double precision, dimension(0:N+1), intent(in) :: pres

      integer :: ii


      J = 0




!
! for the cost function, integrate vertical pressure (direction normal
! to the gross flow direction) on both upper and lower nozzle surfaces
!-----------------------------------------------------------------------

      do ii = 1,N
             J = J + 2*dx*pres(ii)
      end do



      end subroutine jcost




     
      subroutine qvector(u1, u2, u5, q1, q2, q5, S, Sx)

      implicit none

      double precision, intent(in) :: u1, u2, u5, S, Sx
      double precision, intent(out) :: q1, q2, q5
      
      double precision :: pres, gam

      gam = 1.40


      

! clculate pressure
      pres = (gam-1)*(u5 - 0.5*(u2**2)/u1)/S


! calculate q vector
      q1 = 0.0d0
      q2 = pres*Sx
      q5 = 0.0d0


      end subroutine qvector



     
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
     &                  rbc, tau, cmax, c
      double precision, dimension(0:N+1) :: u1b, u2b, u5b,
     &                          f1, f2, f5, f1b, f2b, f5b, 
     &                          q1, q2, q5, q1b, q2b, q5b


      gam = 1.4d0
      cmax = 0.0d0





!
!
! calculate time step
!
!-----------------------------------------------------------------------
      do ii = 0,N+1
      c = dabs(dsqrt(gam*pres(ii)/u1(ii))) + dabs(u2(ii)/u1(ii))
!      print*, c
      if (c .GE. cmax) then
              cmax = c
      end if
      end do

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

      subroutine main(area_rat, J, N)
      
      implicit none
      
      integer, intent(in) :: N
      double precision, intent(in) :: area_rat
      double precision, intent(out) :: J

      integer :: ii, jj, tmax
      double precision :: mach_inf, pres_rat, dx, dt, cfl
      double precision, dimension(0:N+1) :: S, Sx, 
     &          pres, u1, u2, u5 
      double precision, dimension(-1:N+1) :: xpos



! define primary flow variables
!-----------------------------------------------------------------------
      
      mach_inf = 0.3d0
      pres_rat = 0.8d0
      cfl = 0.9d0
      tmax = 100000



! allocate sizes for arrays 
!-----------------------------------------------------------------------

!      allocate(S(0:N+1))
!      allocate(Sx(0:N+1))
!      allocate(xpos(-1:N+1))
!      allocate(pres(0:N+1))
!      allocate(u1(0:N+1))
!      allocate(u2(0:N+1))
!      allocate(u5(0:N+1))



! create mesh 
!-----------------------------------------------------------------------
!
!             -1   0   1   2  ... N-2 N-1  N  N+1
!
! position         0                       L
! points       x   o   o   o  ...  o   o   o   x 
! cells          H   C   C    ...    C   C   H
!
!                0   1   2    ...   N-1  N  N+1
!

      call meshing(area_rat, N, S, Sx, dx, xpos)


! initialise variables 
!-----------------------------------------------------------------------

      call initialise(N, mach_inf, pres_rat, pres, u1, u2, u5)


! start time stepping
!-----------------------------------------------------------------------
      do ii = 1,tmax
             call boundaries(N, mach_inf, pres_rat, pres, u1, u2, u5)
             call update(N, mach_inf, pres_rat, pres, u1, u2, u5, 
     &                          S, Sx, dt, dx, cfl)
      end do



! calculate cost function
!-----------------------------------------------------------------------

      call jcost(N, J, pres, dx)




! correct variables for area
!-----------------------------------------------------------------------
!      do ii = 0,N+1
!             u1(ii) = u1(ii)/S(ii)
!             u2(ii) = u2(ii)/S(ii)
!             u5(ii) = u5(ii)/S(ii)
!      end do





! write info to file
!-----------------------------------------------------------------------
      open(100, file='solution.plt')
  
      write(100,*) 'VARIABLES = "x" "u1" "u2" "u5" "pres"'

101   format(5f12.8)

      do ii = 1,N

            write(100,101) xpos(ii),
     &                  0.5*(u1(ii)+u1(ii+1)),
     &                  0.5*(u2(ii)+u2(ii+1)),
     &                  0.5*(u5(ii)+u5(ii+1)),
     &                  0.5*(pres(ii)+pres(ii+1))
                                

      end do
        



! 
!-----------------------------------------------------------------------

      end subroutine main




     
      subroutine meshing(area_rat, N, S, Sx, dx, xpos)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: area_rat
      double precision, intent(out) :: dx
      double precision, dimension(0:N+1), intent(out) :: S, Sx
      double precision, dimension(-1:N+1), intent(out) :: xpos

      integer :: ii
      double precision :: L, x, r, pi







      L = 1.0
      r = area_rat
      pi = 3.1415
      dx = L/N


      do ii = 1,N

            xpos(ii) = ii*L/N
            x = (ii-0.5)*L/N 

            S(ii) = ((r+1)/(2*r))*(1+(r-1)/(r+1)*cos(2*pi*x/L))
            Sx(ii) = pi*((1-r)/(r*L))*sin(2*pi*x/L)

      end do

      S(0) = S(1) 
      S(N+1) = S(N)

      Sx(0) = Sx(1) 
      Sx(N+1) = Sx(N)


      xpos(-1) = -L/N
      xpos(0) = 0
      xpos(N+1) = 1+L/N











      end subroutine meshing




     
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
