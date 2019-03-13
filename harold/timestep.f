!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine timestep(ni,ur,um,pres,dx,CFL,dt)
!
      IMPLICIT NONE
!
      INTEGER:: ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,pres,dx,dt
      DOUBLE PRECISION:: gam,c,CFL,u
      INTEGER:: i
!
      gam=1.4d0
      do i=1,ni
        c=dsqrt(gam*pres(i)/ur(i))
        u=um(i)/ur(i)
        dt(i)=CFL*dx(i)/(dabs(u)+dabs(c))
      enddo
!
      return 
      end subroutine timestep
!
!**********************************************************************
