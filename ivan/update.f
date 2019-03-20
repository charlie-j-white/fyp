!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine update(ni,pres,area,dadx,res1,res2,res3,dt,
     &  ur,um,ue,ur0,um0,ue0)
!
      IMPLICIT NONE
!
      INTEGER:: ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: pres,dadx,area
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: res1,res2,res3,dt 
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,ue
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur0,um0,ue0
      DOUBLE PRECISION:: areac,dadxc,gam
      INTEGER:: i
!
      gam=1.4d0
      do i=1,ni
        areac=0.5d0*(area(i)+area(i-1))
        dadxc=0.5d0*(dadx(i)+dadx(i-1))
        ur(i)=ur0(i)-(dt(i)/areac)*res1(i)
        um(i)=um0(i)-(dt(i)/areac)*(res2(i)-pres(i)*dadxc) 
        ue(i)=ue0(i)-(dt(i)/areac)*res3(i)
        pres(i)=(gam-1.0d0)*(ue(i)-0.5d0*um(i)*um(i)/ur(i))
      enddo
      return
      end subroutine update
!
!**********************************************************************
