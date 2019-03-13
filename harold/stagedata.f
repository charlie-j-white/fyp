!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine stagedata(nstages,alp,cflscale,cfl)
!
      IMPLICIT NONE
!
      DOUBLE PRECISION, DIMENSION(0:6):: alp
      DOUBLE PRECISION:: cflscale,cfl 
      INTEGER:: nstages
!
       alp(0)=1.0d0 
       if(nstages.eq.3) then
         alp(1)=0.25d0
         alp(2)=0.50d0
         alp(3)=1.0d0
         CFL=1.3d0*CFLscale
       elseif(nstages.eq.4) then
         alp(1)=0.11d0
         alp(2)=0.2766d0
         alp(3)=0.5d0
         alp(4)=1.0d0
         CFL=1.6d0*CFLscale
       elseif(nstages.eq.5) then
         alp(1)=0.059d0
         alp(2)=0.14d0
         alp(3)=0.273d0
         alp(4)=0.5d0
         alp(5)=1.0d0
         CFL=2.0d0*CFLscale
       elseif(nstages.eq.6) then
         alp(1)=0.037d0
         alp(2)=0.0851d0
         alp(3)=0.1521d0
         alp(4)=0.2562d0
         alp(5)=0.4512d0
         alp(6)=1.0d0
         CFL=2.5d0*CFLscale
      endif

       return
       end subroutine stagedata
!
!**********************************************************************
