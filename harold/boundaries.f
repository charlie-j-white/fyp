!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!

      subroutine boundaries(ni,pmachin,prat,ur,um,ue,pres)
!
      IMPLICIT NONE
!
      INTEGER:: ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,ue,pres
      DOUBLE PRECISION:: pmachin,gam,pinf,uint,factor,g1,g2,pstag,
     &                   rstag,prat,M2,vel2,pm2,pr
!
!  SUBSONIC INFLOW AT LEFT BOUNDARY 
!
      gam=1.4d0
      g1=1.0d0/(gam-1.0d0)
      g2=gam/(gam-1.0d0)
      factor=(1.0d0 +0.5d0*(gam-1.0d0)*pmachin*pmachin)
      pstag=(1.0d0/gam)*(factor**g2)
      rstag=factor**g1
      M2=pmachin*pmachin
      pinf=1.0d0/gam 
!
      uint=um(1)/ur(1)
      vel2=uint*uint
      pm2=vel2/((gam-1.0d0)*(1.0d0/(gam-1.0d0)+0.5d0*M2-0.5d0*vel2))
      ur(0)=rstag/((1.0d0+0.5d0*(gam-1.0d0)*pm2)**g1)
      pr=pstag/((1.0d0+0.5d0*(gam-1.0d0)*pm2)**g2)
      ue(0)=2.5d0*pr+0.5d0*ur(0)*vel2
      um(0)=ur(0)*uint
      pres(0)=pr

      pres(-1)=pres(0)
      ur(-1)=ur(0)
      um(-1)=um(0)
      ue(-1)=ue(0)
!
!  SUBSONIC OUTFLOW AT RIGHT BOUNDARY
!
      pres(ni+1)=pinf*prat
      ur(ni+1)=ur(ni)
      um(ni+1)=um(ni)
      ue(ni+1)=pres(ni+1)/(gam-1.0d0) + 0.5d0*um(ni+1)*um(ni+1)/ur(ni+1)
      pres(ni+2)=pres(ni+1)
      ur(ni+2)=ur(ni+1)
      um(ni+2)=um(ni+1)
      ue(ni+2)=ue(ni+1)
!
      return
      end subroutine boundaries
!
!**********************************************************************
