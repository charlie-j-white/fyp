!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine nozzle_main(pmachin,arearat,prat,CFLscale,pk,
     &    nstages,maxit,ni,restol,J)
!
      IMPLICIT NONE
!
      INTEGER:: i,ni,nt,nstages,maxit,ns,INFO
      INTEGER, DIMENSION(3,3):: IPIV
      DOUBLE PRECISION:: J,jd,arearat,djda,ans
      DOUBLE PRECISION, DIMENSION(-1:ni+1):: 
     &             x,ur,um,ue,pres,ds,dx,res1,res2,res3,
     &             error1,error2,error3,dt,area,dadx,ur0,um0,ue0,
     &             vec0,vec1,
     &             res1d,res2d,res3d
      DOUBLE PRECISION, DIMENSION(0:6):: alp
      DOUBLE PRECISION, DIMENSION(1:3):: djdu
      DOUBLE PRECISION, DIMENSION(1,3):: adjoint
      DOUBLE PRECISION, DIMENSION(3,3):: fluxjac
      DOUBLE PRECISION, DIMENSION(-1:ni+1,3):: dRda,dRd1,dRd2,dRd3
      DOUBLE PRECISION:: CFL,pmachin,prat,pk,cflscale
      DOUBLE PRECISION:: rms,rms1,rms2,rms3,restol,resmax
!
!  READ INPUT DATA 
!
!      call readconfig(
!     &   pmachin,arearat,prat,CFLscale,pk,nstages,maxit,ni,restol)
!
!  POINTS GO FROM 0-NI
!  CELLS GO FROM 1-NI
!  TWO HALOS NEEDED AT EACH END 
!
      ni=ni-1
!
!      ALLOCATE(x(-1:ni+2))
!      ALLOCATE(ur(-1:ni+2))
!      ALLOCATE(um(-1:ni+2))
!      ALLOCATE(ue(-1:ni+2))
!      ALLOCATE(ur0(-1:ni+2))
!      ALLOCATE(um0(-1:ni+2))
!      ALLOCATE(ue0(-1:ni+2))
!      ALLOCATE(pres(-1:ni+2))
!      ALLOCATE(ds(-1:ni+2))
!      ALLOCATE(dt(-1:ni+2))
!      ALLOCATE(dx(-1:ni+2))
!      ALLOCATE(res1(-1:ni+2))
!      ALLOCATE(res2(-1:ni+2))
!      ALLOCATE(res3(-1:ni+2))
!      ALLOCATE(error1(-1:ni+2))
!      ALLOCATE(error2(-1:ni+2))
!      ALLOCATE(error3(-1:ni+2))
!      ALLOCATE(area(-1:ni+2))
!      ALLOCATE(dadx(-1:ni+2))
!
!  COMPUTE MESH GEOMETRIC PROPERTIES
!
      call griddata(ni,x,ds,dx,area,dadx,arearat)
!
!  SET INITIAL FLOW PROPERTIES
!
      call initialise(ni,pmachin,prat,ur,um,ue,pres)
!
!  BEGIN TIMESTEPPING
!
      open(30,file='convergence.plt') 
      write(30,301) 
301   format(' VARIABLES = "Iteration" "Log10(RMS)" ')
!
      resmax=-10.0d0 
!
!----------------------------------------
      do nt=1,maxit
!----------------------------------------
!
!  SET RUNGE-KUTTA DATA
!
        call stagedata(nstages,alp,cflscale,cfl)
!
!  COMPUTE TIMESTEP IN EACH CELL
!
        call timestep(ni,ur,um,pres,dx,CFL,dt)
!
!  SET ZERO SOLUTION AT EACH TIME LEVEL
!
        do i=1,ni
          ur0(i)=ur(i)
          um0(i)=um(i)
          ue0(i)=ue(i)
        enddo
!
!  LOOP OVER R-K STAGES
!
!----------------------------------------
        do ns=1,nstages
!----------------------------------------
          do i=1,ni
            dt(i)=dt(i)*alp(ns)/alp(ns-1)
          enddo
!
!  SET DATA OUTSIDE BOUNDARIES
!  COMPUTE RESIDUAL VECTOR
!  UPDATE SOLUTION USING TIMESTEPPING SCHEME
! 
          call boundaries(ni,pmachin,prat,ur,um,ue,pres)
          call resid(ni,ur,um,ue,dx,ds,area,pk,res1,res2,res3)
          call update(ni,pres,area,dadx,res1,res2,res3,dt, 
     &                ur,um,ue,ur0,um0,ue0)
!
!--------------R-K STAGE LOOP------------
        enddo
!----------------------------------------
!
      rms=0.0d0
      do i=1,ni
cc        print*,i,rms
        rms1=(ur(i)-ur0(i))**2
        rms2=(um(i)-um0(i))**2
        rms3=(ue(i)-ue0(i))**2
        rms=rms+(rms1+rms2+rms3)/(dt(i)*dt(i))
      enddo
        rms=dsqrt(rms/dfloat(ni))
        if(rms.gt.resmax) resmax=rms
        write(6,601) nt,dlog10(rms/resmax)
        write(30,302) nt,dlog10(rms/resmax)
601     format('ITERATION ',i6,' Log10(RMSERROR) ',f6.2)
302     format(i10,f8.3) 
!
      IF((rms/resmax).le.restol) go to 777
!
!------------TIMESTEP LOOP---------------
      enddo
!----------------------------------------
!
777   continue
      open(20,file='solution.plt') 
      write(20,201) 
201   format(' VARIABLES = "X" "rho" "mom" "e" "pres" ')
202   format(5f12.8) 
      do i=0,ni
      write(20,202) x(i),0.5*(ur(i)+ur(i+1)),
     &  0.5*(um(i)+um(i+1)),0.5*(ue(i)+ue(i+1)),0.5*(pres(i)+pres(i+1))
      enddo
      close(20)
      close(30)
!
!      DEALLOCATE(x,ur,um,ue,pres,ds,dt,dx,res1,res2,res3)
!      DEALLOCATE(ur0,um0,ue0,error1,error2,error3,area,dadx)
!
!----------------------------------------------------------------------
!                     A D J O I N T      C O D E 
!----------------------------------------------------------------------
      print*, "    "
      print*, "    "
      print*, "----------Adjoint code----------"




!     debugging statements      
      i = 110
      print*, "     "
      print*, "u = ", ur(i),um(i),ue(i)
      print*, "R = ", res1(i),res2(i),res3(i)




!     primal output of cost function
!----------------------------------------------------------------------
      call cost(J,arearat,ur,um,ue,ni)

!     create seed vectors
      do i=-1,ni+1
         vec0(i)=0.0d0
         vec1(i)=1.0d0
      end do

!     get dJ/da
      call cost_d(j,jd,arearat,1.0d0,ur,vec0,um,vec0,ue,vec0,ni)
      djda = jd

!     get dJ/du
      call cost_d(j,jd,arearat,0.0d0,ur,vec1,um,vec0,ue,vec0,ni)
      djdu(1) = jd
      call cost_d(j,jd,arearat,0.0d0,ur,vec0,um,vec1,ue,vec0,ni)
      djdu(2) = jd
      call cost_d(j,jd,arearat,0.0d0,ur,vec0,um,vec0,ue,vec1,ni)
      djdu(3) = jd







      i = 110;
!     primal output of adapted residual function
!----------------------------------------------------------------------
      call jRES(ni,ur,um,ue,pk,res1,res2,res3,arearat)

!     debugging statements      
      print*, "     "
      i = 110
      print*, "u = ", ur(i),um(i),ue(i)
      print*, "R = ", res1(i),res2(i),res3(i)


!     get dRda
      call JRES_D(ni,ur,vec0,um,vec0,ue,vec0,pk,
     &           res1,res1d,res2,res2d,res3,res3d,arearat,1.0d0)
      dRda(:,1) = res1d
      dRda(:,2) = res2d
      dRda(:,3) = res3d
     

!     get dRdu
      call JRES_D(ni,ur,vec1,um,vec0,ue,vec0,pk,
     &           res1,res1d,res2,res2d,res3,res3d,arearat,0.0d0)
      dRd1(:,1) = res1d
      dRd1(:,2) = res2d
      dRd1(:,3) = res3d
     
      call JRES_D(ni,ur,vec0,um,vec1,ue,vec0,pk,
     &           res1,res1d,res2,res2d,res3,res3d,arearat,0.0d0)
      dRd2(:,1) = res1d
      dRd2(:,2) = res2d
      dRd2(:,3) = res3d
     
      call JRES_D(ni,ur,vec0,um,vec0,ue,vec1,pk,
     &           res1,res1d,res2,res2d,res3,res3d,arearat,0.0d0)
      dRd3(:,1) = res1d
      dRd3(:,2) = res2d
      dRd3(:,3) = res3d
     


!     solve adjoint equation for adjoint vector 
!----------------------------------------------------------------------

      ans = 0.0d0
      do i=1,ni

      adjoint(1,1) = dJdu(1)
      adjoint(1,2) = dJdu(2)
      adjoint(1,3) = dJdu(3)
      
      fluxjac(1,1) = dRd1(i,1)
      fluxjac(2,1) = dRd1(i,2)
      fluxjac(3,1) = dRd1(i,3)
      fluxjac(1,2) = dRd2(i,1)
      fluxjac(2,2) = dRd2(i,2)
      fluxjac(3,2) = dRd2(i,3)
      fluxjac(1,3) = dRd3(i,1)
      fluxjac(2,3) = dRd3(i,2)
      fluxjac(3,3) = dRd3(i,3)

      call DGETRF(3,3,fluxjac,3,IPIV,INFO)
      call DGETRS('N',3,1,fluxjac,3,IPIV,adjoint,3,INFO)

      ans = ans -(adjoint(1,1)*dRda(i,1) + 
     &   adjoint(1,2)*dRda(i,2)  + adjoint(1,3)*dRda(i,3) -djda)

      print*, i, -(adjoint(1,1)*dRda(i,1) + 
     &   adjoint(1,2)*dRda(i,2)  + adjoint(1,3)*dRda(i,3) -djda)

      end do

      print*, ans/ni
      print*, "    "







!      stop
      end subroutine nozzle_main
!
!**********************************************************************
