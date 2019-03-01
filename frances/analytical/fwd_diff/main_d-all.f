C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.14 (r7260) - 18 Jan 2019 10:11
C
C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.14 (r7260) - 18 Jan 2019 10:11
C
C  Differentiation of update in forward (tangent) mode:
C   variations   of useful results: g
C   with respect to varying inputs: x1 x2
C   RW status of diff variables: g:out x1:in x2:in
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
C
      SUBROUTINE UPDATE_D(x1, x1d, x2, x2d, g, gd, gx, gy)
      IMPLICIT NONE
C
C
      DOUBLE PRECISION x1, x2, g, gx, gy
      DOUBLE PRECISION x1d, x2d, gd
      INTRINSIC SIN
      INTRINSIC COS
C
C
      gd = x1d*SIN(x1*x2) + x1*(x1d*x2+x1*x2d)*COS(x1*x2) + x2d
      g = x1*SIN(x1*x2) + x2
      gx = SIN(x1*x2) + x1*x2*COS(x1*x2)
      gy = x1*x1*COS(x1*x2) + 1
      END

