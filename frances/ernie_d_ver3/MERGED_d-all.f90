!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7260) - 18 Jan 2019 10:11
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7260) - 18 Jan 2019 10:11
!
!  Differentiation of main in forward (tangent) mode:
!   variations   of useful results: j
!   with respect to varying inputs: area_rat
!   RW status of diff variables: j:out area_rat:in
!
SUBROUTINE MAIN_D(area_rat, area_ratd, j, jd)
  IMPLICIT NONE
!
!
!
!
! 
!-----------------------------------------------------------------------
!
  DOUBLE PRECISION, INTENT(IN) :: area_rat, j
  DOUBLE PRECISION, INTENT(IN) :: area_ratd, jd
!
  INTEGER :: n, ii, jj, tmax
  DOUBLE PRECISION :: mach_inf, pres_rat, dx, dt, cfl
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: s, sx, xpos, pres, u1, &
& u2, u5
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sd, sxd, presd, u1d, &
& u2d, u5d
!
!
!
! define primary flow variables
!-----------------------------------------------------------------------
  mach_inf = 0.3d0
  pres_rat = 0.8d0
  cfl = 0.9d0
  n = 1000
  tmax = 100000
!
!
!
! allocate sizes for arrays 
!-----------------------------------------------------------------------
!
  ALLOCATE(sd(0:n+1))
  ALLOCATE(s(0:n+1))
  ALLOCATE(sxd(0:n+1))
  ALLOCATE(sx(0:n+1))
  ALLOCATE(xpos(-1:n+1))
  ALLOCATE(presd(0:n+1))
  ALLOCATE(pres(0:n+1))
  ALLOCATE(u1d(0:n+1))
  ALLOCATE(u1(0:n+1))
  ALLOCATE(u2d(0:n+1))
  ALLOCATE(u2(0:n+1))
  ALLOCATE(u5d(0:n+1))
  ALLOCATE(u5(0:n+1))
!
!
!
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
!
  CALL MESHING_D(area_rat, area_ratd, n, s, sd, sx, sxd, dx, xpos)
!
!
! initialise variables 
!-----------------------------------------------------------------------
!
  CALL INITIALISE(n, mach_inf, pres_rat, pres, u1, u2, u5)
!
!
! start time stepping
!-----------------------------------------------------------------------
  DO ii=1,tmax
    CALL BOUNDARIES_D(n, mach_inf, pres_rat, pres, presd, u1, u1d, u2, &
&               u2d, u5, u5d)
    CALL UPDATE_D(n, mach_inf, pres_rat, pres, presd, u1, u1d, u2, u2d, &
&           u5, u5d, s, sd, sx, sxd, dt, dx, cfl)
  END DO
!
!
!
! calculate cost function
!-----------------------------------------------------------------------
!
  CALL JCOST_D(n, j, jd, pres, presd, dx)
!
  PRINT*, j
!
!
!
! correct variables for area
!-----------------------------------------------------------------------
!      do ii = 0,N+1
!             u1(ii) = u1(ii)/S(ii)
!             u2(ii) = u2(ii)/S(ii)
!             u5(ii) = u5(ii)/S(ii)
!      end do
!
!
!
!
!
! write info to file
!-----------------------------------------------------------------------
  OPEN(100, file='solution.plt') 
  WRITE(100, *) 'VARIABLES = "x" "u1" "u2" "u5" "pres"'
!
  DO ii=1,n
!
    WRITE(100, 101) xpos(ii), 0.5*(u1(ii)+u1(ii+1)), 0.5*(u2(ii)+u2(ii+1&
&   )), 0.5*(u5(ii)+u5(ii+1)), 0.5*(pres(ii)+pres(ii+1))
  END DO
!
 101 FORMAT(5f12.8)
END SUBROUTINE MAIN_D

!  Differentiation of boundaries in forward (tangent) mode:
!   variations   of useful results: u1 u2 u5 pres
!   with respect to varying inputs: u1 u2 u5 pres
!
!
!
!
SUBROUTINE BOUNDARIES_D(n, mach_inf, pres_rat, pres, presd, u1, u1d, u2&
& , u2d, u5, u5d)
  IMPLICIT NONE
!
!
!
!
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: mach_inf, pres_rat
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: pres, u1, u2, u5
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: presd, u1d, u2d, &
& u5d
!
  INTEGER :: ii
  DOUBLE PRECISION :: gam, pinf, pstag, rinf, rstag, ubc, mbc, pbc, rbc
  DOUBLE PRECISION :: ubcd, mbcd, pbcd, rbcd
  DOUBLE PRECISION :: pwx1
  DOUBLE PRECISION :: pwx1d
  DOUBLE PRECISION :: pwy1
  DOUBLE PRECISION :: pwr1
  DOUBLE PRECISION :: pwr1d
!
!
! useful values near the inflow
  gam = 1.4
  pinf = 1.0/gam
  pwx1 = 1 + 0.5*(gam-1)*mach_inf**2
  pwy1 = gam/(gam-1)
  pwr1 = pwx1**pwy1
  pstag = pinf*pwr1
  rinf = 1.0
  pwx1 = 1 + 0.5*(gam-1)*mach_inf**2
  pwy1 = 1.0/(gam-1)
  pwr1 = pwx1**pwy1
  rstag = rinf*pwr1
!
  ubcd = (u2d(1)*u1(1)-u2(1)*u1d(1))/u1(1)**2
  ubc = u2(1)/u1(1)
  mbcd = ((ubcd*ubc+ubc*ubcd)*(1+0.5*(gam-1)*mach_inf**2-0.5*(gam-1)*ubc&
&   **2)+ubc**3*0.5*(gam-1)*2*ubcd)/(1+0.5*(gam-1)*mach_inf**2-0.5*(gam-&
&   1)*ubc**2)**2
  mbc = ubc*ubc/(1+0.5*(gam-1)*mach_inf**2-0.5*(gam-1)*ubc**2)
  pwx1d = 0.5*(gam-1)*2*mbc*mbcd
  pwx1 = 1 + 0.5*(gam-1)*mbc**2
  pwy1 = gam/(gam-1)
  IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. pwy1 .EQ. INT(pwy1))) THEN
    pwr1d = pwy1*pwx1**(pwy1-1)*pwx1d
  ELSE IF (pwx1 .EQ. 0.0 .AND. pwy1 .EQ. 1.0) THEN
    pwr1d = pwx1d
  ELSE
    pwr1d = 0.0
  END IF
  pwr1 = pwx1**pwy1
  pbcd = -(pstag*pwr1d/pwr1**2)
  pbc = pstag/pwr1
  pwx1d = 0.5*(gam-1)*2*mbc*mbcd
  pwx1 = 1 + 0.5*(gam-1)*mbc**2
  pwy1 = 1.0/(gam-1)
  IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. pwy1 .EQ. INT(pwy1))) THEN
    pwr1d = pwy1*pwx1**(pwy1-1)*pwx1d
  ELSE IF (pwx1 .EQ. 0.0 .AND. pwy1 .EQ. 1.0) THEN
    pwr1d = pwx1d
  ELSE
    pwr1d = 0.0
  END IF
  pwr1 = pwx1**pwy1
  rbcd = -(rstag*pwr1d/pwr1**2)
  rbc = rstag/pwr1
!
! halo cells at the inflow
  presd(0) = pbcd
  pres(0) = pbc
  u1d(0) = rbcd
  u1(0) = rbc
  u2d(0) = rbcd*ubc + rbc*ubcd
  u2(0) = rbc*ubc
  u5d(0) = pbcd/(gam-1) + 0.5*(rbcd*ubc**2+rbc*2*ubc*ubcd)
  u5(0) = pbc/(gam-1) + 0.5*rbc*ubc**2
!
! halo cells at the outflow
  presd(n+1) = 0.D0
  pres(n+1) = pinf*pres_rat
  u1d(n+1) = u1d(n)
  u1(n+1) = u1(n)
  u2d(n+1) = u2d(n)
  u2(n+1) = u2(n)
  u5d(n+1) = presd(n+1)/(gam-1) + (0.5*2*u2(n+1)*u2d(n+1)*u1(n+1)-0.5*u2&
&   (n+1)**2*u1d(n+1))/u1(n+1)**2
  u5(n+1) = pres(n+1)/(gam-1) + 0.5*u2(n+1)**2/u1(n+1)
END SUBROUTINE BOUNDARIES_D

!  Differentiation of jcost in forward (tangent) mode:
!   variations   of useful results: j
!   with respect to varying inputs: pres
!
!
!
SUBROUTINE JCOST_D(n, j, jd, pres, presd, dx)
  IMPLICIT NONE
!
!
!
!
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx
  DOUBLE PRECISION, INTENT(OUT) :: j
  DOUBLE PRECISION, INTENT(OUT) :: jd
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(IN) :: pres
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(IN) :: presd
!
  INTEGER :: ii
!
!
  j = 0
  jd = 0.D0
!
!
!
!
!
! for the cost function, integrate vertical pressure (direction normal
! to the gross flow direction) on both upper and lower nozzle surfaces
!-----------------------------------------------------------------------
!
  DO ii=1,n
    jd = jd + 2*dx*presd(ii)
    j = j + 2*dx*pres(ii)
  END DO
END SUBROUTINE JCOST_D

!  Differentiation of meshing in forward (tangent) mode:
!   variations   of useful results: s sx
!   with respect to varying inputs: area_rat
!
!
!
!
SUBROUTINE MESHING_D(area_rat, area_ratd, n, s, sd, sx, sxd, dx, xpos)
  IMPLICIT NONE
!
!
!
!
!
!
!
!
!
!
!
!
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: area_rat
  DOUBLE PRECISION, INTENT(IN) :: area_ratd
  DOUBLE PRECISION, INTENT(OUT) :: dx
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: s, sx
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: sd, sxd
  DOUBLE PRECISION, DIMENSION(-1:n+1), INTENT(OUT) :: xpos
!
  INTEGER :: ii
  DOUBLE PRECISION :: l, x, r, pi
  DOUBLE PRECISION :: rd
  INTRINSIC COS
  INTRINSIC SIN
  DOUBLE PRECISION :: arg1
  INTEGER :: ii1
!
!
!
!
!
!
!
  l = 1.0
  rd = area_ratd
  r = area_rat
  pi = 3.1415
  dx = l/n
  DO ii1=0,n+1
    sd(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    sxd(ii1) = 0.D0
  END DO
!
!
  DO ii=1,n
!
    xpos(ii) = ii*l/n
    x = (ii-0.5)*l/n
!
    arg1 = 2*pi*x/l
    sd(ii) = (rd*2*r-(r+1)*2*rd)*(1+(r-1)/(r+1)*COS(arg1))/(2**2*r**2) +&
&     COS(arg1)*(rd*(r+1)-(r-1)*rd)/((r+1)*2*r)
    s(ii) = (r+1)/(2*r)*(1+(r-1)/(r+1)*COS(arg1))
    arg1 = 2*pi*x/l
    sxd(ii) = pi*SIN(arg1)*(-(rd*r*l)-(1-r)*l*rd)/(r**2*l**2)
    sx(ii) = pi*((1-r)/(r*l))*SIN(arg1)
  END DO
!
!
  sd(0) = sd(1)
  s(0) = s(1)
  sd(n+1) = sd(n)
  s(n+1) = s(n)
!
  sxd(0) = sxd(1)
  sx(0) = sx(1)
  sxd(n+1) = sxd(n)
  sx(n+1) = sx(n)
!
!
  xpos(-1) = -(l/n)
  xpos(0) = 0
  xpos(n+1) = 1 + l/n
END SUBROUTINE MESHING_D

!  Differentiation of update in forward (tangent) mode:
!   variations   of useful results: u1 u2 u5 pres
!   with respect to varying inputs: s u1 u2 u5 sx pres
!
!
!
SUBROUTINE UPDATE_D(n, mach_inf, pres_rat, pres, presd, u1, u1d, u2, u2d&
& , u5, u5d, s, sd, sx, sxd, dt, dx, cfl)
  IMPLICIT NONE
!
!
!
!
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: mach_inf, pres_rat, dx, cfl
  DOUBLE PRECISION, INTENT(OUT) :: dt
  DOUBLE PRECISION :: dtd
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: pres, u1, u2, u5
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(OUT) :: presd, u1d, u2d, &
& u5d
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(IN) :: s, sx
  DOUBLE PRECISION, DIMENSION(0:n+1), INTENT(IN) :: sd, sxd
!
  INTEGER :: ii
  DOUBLE PRECISION :: gam, pinf, pstag, rinf, rstag, ubc, mbc, pbc, rbc&
& , tau, cmax, c
  DOUBLE PRECISION :: taud, cmaxd, cd
  DOUBLE PRECISION, DIMENSION(0:n+1) :: u1b, u2b, u5b, f1, f2, f5, f1b, &
& f2b, f5b, q1, q2, q5, q1b, q2b, q5b
  DOUBLE PRECISION, DIMENSION(0:n+1) :: u1bd, u2bd, u5bd, f1d, f2d, f5d&
& , f1bd, f2bd, f5bd, q2d, q2bd
  INTRINSIC DSQRT
  INTRINSIC DABS
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: dabs0
  DOUBLE PRECISION :: dabs0d
  DOUBLE PRECISION :: dabs1
  DOUBLE PRECISION :: dabs1d
  DOUBLE PRECISION :: arg1
  DOUBLE PRECISION :: arg1d
  INTEGER :: ii1
!
!
  gam = 1.4d0
  cmax = 0.0d0
  cmaxd = 0.D0
!
!
!
!
!
!
!
! calculate time step
!
!-----------------------------------------------------------------------
  DO ii=0,n+1
    arg1d = (gam*presd(ii)*u1(ii)-gam*pres(ii)*u1d(ii))/u1(ii)**2
    arg1 = gam*pres(ii)/u1(ii)
    IF (arg1 .EQ. 0.0) THEN
      x1d = 0.D0
    ELSE
      x1d = arg1d/(2.D0*DSQRT(arg1))
    END IF
    x1 = DSQRT(arg1)
    IF (x1 .GE. 0.) THEN
      dabs0d = x1d
      dabs0 = x1
    ELSE
      dabs0d = -x1d
      dabs0 = -x1
    END IF
    IF (u2(ii)/u1(ii) .GE. 0.) THEN
      dabs1d = (u2d(ii)*u1(ii)-u2(ii)*u1d(ii))/u1(ii)**2
      dabs1 = u2(ii)/u1(ii)
    ELSE
      dabs1d = -((u2d(ii)*u1(ii)-u2(ii)*u1d(ii))/u1(ii)**2)
      dabs1 = -(u2(ii)/u1(ii))
    END IF
    cd = dabs0d + dabs1d
    c = dabs0 + dabs1
!      print*, c
    IF (c .GE. cmax) THEN
      cmaxd = cd
      cmax = c
    END IF
  END DO
!
  dtd = -(cfl*dx*cmaxd/cmax**2)
  dt = cfl*dx/cmax
!dt = 0.001
  taud = dtd/dx
  tau = dt/dx
  DO ii1=0,n+1
    f1d(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    f2d(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    f5d(ii1) = 0.D0
  END DO
!
!
!
!      print*, dt
!
!
!
!
!
!
! initialise variables:
!       calculate f vector 0-->N+1
!       calculate q vector 0-->N+1
!-----------------------------------------------------------------------
!
  DO ii=0,n+1
    CALL FVECTOR_D(u1(ii), u1d(ii), u2(ii), u2d(ii), u5(ii), u5d(ii), f1&
&            (ii), f1d(ii), f2(ii), f2d(ii), f5(ii), f5d(ii))
  END DO
  DO ii1=0,n+1
    q2d(ii1) = 0.D0
  END DO
  DO ii=0,n+1
    CALL QVECTOR_D(u1(ii), u1d(ii), u2(ii), u2d(ii), u5(ii), u5d(ii), q1&
&            (ii), q2(ii), q2d(ii), q5(ii), s(ii), sd(ii), sx(ii), sxd(&
&            ii))
  END DO
  DO ii1=0,n+1
    u1bd(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    u5bd(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    u2bd(ii1) = 0.D0
  END DO
!
!
!
! calculate predictor terms:
!       calculate u vector predictor terms 0-->N
!       calculate f vector predictor terms 0-->N
!       calculate q vector predictor terms 0-->N
!-----------------------------------------------------------------------
!
  DO ii=0,n
    u1bd(ii) = u1d(ii) - taud*(f1(ii+1)-f1(ii)) - tau*(f1d(ii+1)-f1d(ii)&
&     ) + q1(ii)*dtd
    u1b(ii) = u1(ii) - tau*(f1(ii+1)-f1(ii)) + dt*q1(ii)
    u2bd(ii) = u2d(ii) - taud*(f2(ii+1)-f2(ii)) - tau*(f2d(ii+1)-f2d(ii)&
&     ) + dtd*q2(ii) + dt*q2d(ii)
    u2b(ii) = u2(ii) - tau*(f2(ii+1)-f2(ii)) + dt*q2(ii)
    u5bd(ii) = u5d(ii) - taud*(f5(ii+1)-f5(ii)) - tau*(f5d(ii+1)-f5d(ii)&
&     ) + q5(ii)*dtd
    u5b(ii) = u5(ii) - tau*(f5(ii+1)-f5(ii)) + dt*q5(ii)
  END DO
  DO ii1=0,n+1
    f1bd(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    f5bd(ii1) = 0.D0
  END DO
  DO ii1=0,n+1
    f2bd(ii1) = 0.D0
  END DO
!
  DO ii=0,n
    CALL FVECTOR_D(u1b(ii), u1bd(ii), u2b(ii), u2bd(ii), u5b(ii), u5bd(&
&            ii), f1b(ii), f1bd(ii), f2b(ii), f2bd(ii), f5b(ii), f5bd(ii&
&            ))
  END DO
  DO ii1=0,n+1
    q2bd(ii1) = 0.D0
  END DO
!
  DO ii=0,n
    CALL QVECTOR_D(u1b(ii), u1bd(ii), u2b(ii), u2bd(ii), u5b(ii), u5bd(&
&            ii), q1b(ii), q2b(ii), q2bd(ii), q5b(ii), s(ii), sd(ii), sx&
&            (ii), sxd(ii))
  END DO
!
!
!
! update:
!       calculate fully updated u vector 1-->N
!-----------------------------------------------------------------------
!
  DO ii=1,n
    u1d(ii) = 0.5*(u1d(ii)+u1bd(ii)) - 0.5*(taud*(f1b(ii)-f1b(ii-1))+tau&
&     *(f1bd(ii)-f1bd(ii-1))) + q1b(ii)*dtd/2
    u1(ii) = 0.5*(u1(ii)+u1b(ii)) - 0.5*tau*(f1b(ii)-f1b(ii-1)) + dt*q1b&
&     (ii)/2
    u2d(ii) = 0.5*(u2d(ii)+u2bd(ii)) - 0.5*(taud*(f2b(ii)-f2b(ii-1))+tau&
&     *(f2bd(ii)-f2bd(ii-1))) + (dtd*q2b(ii)+dt*q2bd(ii))/2
    u2(ii) = 0.5*(u2(ii)+u2b(ii)) - 0.5*tau*(f2b(ii)-f2b(ii-1)) + dt*q2b&
&     (ii)/2
    u5d(ii) = 0.5*(u5d(ii)+u5bd(ii)) - 0.5*(taud*(f5b(ii)-f5b(ii-1))+tau&
&     *(f5bd(ii)-f5bd(ii-1))) + q5b(ii)*dtd/2
    u5(ii) = 0.5*(u5(ii)+u5b(ii)) - 0.5*tau*(f5b(ii)-f5b(ii-1)) + dt*q5b&
&     (ii)/2
    presd(ii) = (gam-1)*(u5d(ii)-(0.5d0*2*u2(ii)*u2d(ii)*u1(ii)-0.5d0*u2&
&     (ii)**2*u1d(ii))/u1(ii)**2)
    pres(ii) = (gam-1)*(u5(ii)-0.5d0*u2(ii)**2/u1(ii))
  END DO
END SUBROUTINE UPDATE_D

!  Differentiation of fvector in forward (tangent) mode:
!   variations   of useful results: f1 f2 f5
!   with respect to varying inputs: u1 u2 u5
!
!
!
!
SUBROUTINE FVECTOR_D(u1, u1d, u2, u2d, u5, u5d, f1, f1d, f2, f2d, f5, &
& f5d)
  IMPLICIT NONE
!
!
!
  DOUBLE PRECISION, INTENT(IN) :: u1, u2, u5
  DOUBLE PRECISION, INTENT(IN) :: u1d, u2d, u5d
  DOUBLE PRECISION, INTENT(OUT) :: f1, f2, f5
  DOUBLE PRECISION, INTENT(OUT) :: f1d, f2d, f5d
  DOUBLE PRECISION :: pres, vel, gam
  DOUBLE PRECISION :: presd, veld
!
! area terms S inlculded in the u variable, or cancelled
  gam = 1.40
!
!
!
! clculate pressure
  presd = (gam-1)*(u5d-(0.5*2*u2*u2d*u1-0.5*u2**2*u1d)/u1**2)
  pres = (gam-1)*(u5-0.5*u2**2/u1)
!
! calculate velocity
  veld = (u2d*u1-u2*u1d)/u1**2
  vel = u2/u1
!
!
! calculate f vector
  f1d = veld*u1 + vel*u1d
  f1 = vel*u1
  f2d = veld*u2 + vel*u2d + presd
  f2 = vel*u2 + pres
  f5d = veld*u5 + vel*u5d + presd*vel + pres*veld
  f5 = vel*u5 + pres*vel
END SUBROUTINE FVECTOR_D

!  Differentiation of qvector in forward (tangent) mode:
!   variations   of useful results: q2
!   with respect to varying inputs: s u1 u2 u5 sx
!
!
!
!
SUBROUTINE QVECTOR_D(u1, u1d, u2, u2d, u5, u5d, q1, q2, q2d, q5, s, sd, &
& sx, sxd)
  IMPLICIT NONE
!
!
!
  DOUBLE PRECISION, INTENT(IN) :: u1, u2, u5, s, sx
  DOUBLE PRECISION, INTENT(IN) :: u1d, u2d, u5d, sd, sxd
  DOUBLE PRECISION, INTENT(OUT) :: q1, q2, q5
  DOUBLE PRECISION, INTENT(OUT) :: q2d
  DOUBLE PRECISION :: pres, gam
  DOUBLE PRECISION :: presd
!
  gam = 1.40
!
!
!
! clculate pressure
  presd = ((gam-1)*(u5d-(0.5*2*u2*u2d*u1-0.5*u2**2*u1d)/u1**2)*s-(gam-1)&
&   *(u5-0.5*u2**2/u1)*sd)/s**2
  pres = (gam-1)*(u5-0.5*u2**2/u1)/s
!
!
! calculate q vector
  q1 = 0.0d0
  q2d = presd*sx + pres*sxd
  q2 = pres*sx
  q5 = 0.0d0
END SUBROUTINE QVECTOR_D

