!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine initialise(ni,pmachin,prat,ur,um,ue,pres)
!
      IMPLICIT NONE
!
      INTEGER:: ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,ue,pres
      DOUBLE PRECISION:: pmachin,gam,rinf,pinf,uinf,
     &                   prat,rout,uloc,ptot
      INTEGER:: i
!
      gam=1.4d0
      pinf=1.0d0/gam 
      rinf=1.0d0 
      rout=((pinf*prat/pinf)*(rinf**gam))**(1.0d0/gam)
      uinf=pmachin
      ptot=pinf+0.5*rinf*uinf*uinf
!
      do i=1,ni
        pres(i)=pinf+(pinf*prat-pinf)*dfloat(i-1)/dfloat(ni-1) 
        ur(i)=rinf+(rout-rinf)*dfloat(i-1)/dfloat(ni-1) 
        uloc=dsqrt((ptot-pres(i))/(0.5d0*ur(i)))
        um(i)=ur(i)*uloc
        ue(i)=pres(i)/(gam-1.0d0) + 0.5d0*ur(i)*uloc*uloc
      enddo
!
      return 
      end subroutine initialise
!
!**********************************************************************
