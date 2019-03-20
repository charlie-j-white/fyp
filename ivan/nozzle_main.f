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
      INTEGER:: i,ni,nt,nstages,maxit,ns,INFO,ii
      INTEGER, DIMENSION(3*ni,3*ni):: IPIV
!      INTEGER, DIMENSION(3*ni+9,3*ni+9):: IPIV
!      INTEGER, DIMENSION(3,3):: IPIV
      DOUBLE PRECISION:: J,jd,arearat,djda,ans
      DOUBLE PRECISION, DIMENSION(-1:ni+1):: 
     &             x,ur,um,ue,pres,ds,dx,res1,res2,res3,
     &             error1,error2,error3,dt,area,dadx,ur0,um0,ue0,
     &             vec0,vec1,
     &             res1d,res2d,res3d
      DOUBLE PRECISION, DIMENSION(0:6):: alp
      DOUBLE PRECISION, DIMENSION(1:3):: djdu
      DOUBLE PRECISION, DIMENSION(3,3):: fluxjac
      DOUBLE PRECISION, DIMENSION(-1:ni+1,3):: dRda,dRd1,dRd2,dRd3
      DOUBLE PRECISION:: CFL,pmachin,prat,pk,cflscale
      DOUBLE PRECISION:: rms,rms1,rms2,rms3,restol,resmax
      DOUBLE PRECISION, DIMENSION(1,1:3*ni+9)::flow,residual,flowd,
     & djdu2, dRda2, dRdu2, adjoint
      DOUBLE PRECISION, DIMENSION(1:3*ni+9,1:3*ni+9):: fluxjac2
      DOUBLE PRECISION, DIMENSION(1:3*ni,1:3*ni+9):: hold1
      DOUBLE PRECISION, DIMENSION(1:3*ni,1:3*ni):: hold2
      DOUBLE PRECISION, DIMENSION(1,1:3*ni):: hold3,hold4,hold5
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

!     concatenate flow and residual vectors
      flow(1,     1:  ni+3) = ur
      flow(1,  ni+4:2*ni+6) = um
      flow(1,2*ni+7:3*ni+9) = ue
      residual(1,     1:  ni+3) = res1
      residual(1,  ni+4:2*ni+6) = res2
      residual(1,2*ni+7:3*ni+9) = res3






!     primal output of cost function
!----------------------------------------------------------------------
      call cost(J,arearat,flow,ni)




!     get dJ/du

      do i=1,3*ni+9
        flowd(1,i) = 0.0d0
      end do

      flowd(1,1) = 1.0d0
      call COST_D(j,jd,arearat, 0.0d0, flow, flowd, ni)
      djdu2(1,1) = jd

      do i=2,3*ni+9
        flowd(1,i-1) = 0.0d0
        flowd(1,i) = 1.0d0
        call COST_D(j,jd,arearat, 0.0d0, flow, flowd, ni)
        djdu2(1,i) = jd
      end do

      flowd(1,3*ni+9) = 0.0d0



!     get dJ/da
      call COST_D(j,jd,arearat, 1.0d0, flow, flowd, ni)
      djda = jd







!     primal output of adapted residual function
!----------------------------------------------------------------------

      call jRES(ni,pk,pmachin,prat,arearat,residual,flow)






!     get dRdu

      do i=1,3*ni+9
        flowd(1,i) = 0.0d0
      end do

      flowd(1,1) = 1.0d0
      call JRES_D(ni, pk,pmachin,prat, 
     & arearat, 0.0d0, residual, dRdu2, flow, flowd)
      fluxjac2(1,:) = dRdu2(1,:)

      do i=2,3*ni+9
        dRdu2 = 0.0d0
        flowd(1,i-1) = 0.0d0
        flowd(1,i) = 1.0d0
        call JRES_D(ni,pk,pmachin,prat,
     &  arearat, 0.0d0, residual, dRdu2, flow, flowd)
        fluxjac2(i,:) = dRdu2(1,:)
      end do

      flowd(1,3*ni+9) = 0.0d0





!     get dRda
      call JRES_D(ni, pk,pmachin,prat,
     &  arearat, 1.0d0, residual, dRda2, flow, flowd)











!     solve adjoint equation for adjoint vector 
!----------------------------------------------------------------------


!!      assign adjoint vector for linear solve
!      adjoint(1,:) = djdu2(1,:)


!     remove halo rows from important vectors

      do i=3,ni+2
      hold1(i-2,:) = fluxjac2(i,:)
      end do

      do i=ni+6,2*ni+5
      hold1(i-5,:) = fluxjac2(i,:)
      end do

      do i=2*ni+9,3*ni+8
      hold1(i-8,:) = fluxjac2(i,:)
      end do

      do i=3,ni+2
      hold2(:,i-2) = fluxjac2(:,i)
      hold3(1,i-2) = djdu2(1,i)
      hold4(1,i-2) = dRda2(1,i)
      hold5(1,i-2) = dRda2(1,i)
      end do

      do i=ni+6,2*ni+5
      hold2(:,i-5) = fluxjac2(:,i)
      hold3(1,i-5) = djdu2(1,i)
      hold4(1,i-5) = dRda2(1,i)
      end do

      do i=2*ni+9,3*ni+8
      hold2(:,i-8) = fluxjac2(:,i)
      hold3(1,i-8) = djdu2(1,i)
      hold4(1,i-8) = dRda2(1,i)
      end do




!     count number of non-zero elements on each row of flux jacobian

      print*, "      "
      print*, "      "
      do ii=1,3*ni

      ans = 0.0d0
      do i=1,3*ni
!        if( fluxjac2(ii,i) /= 0.0d0 ) then
!        if( fluxjac2(i,ii) /= 0.0d0 ) then
        if( hold2(i,ii) /= 0.0d0 ) then
!        if( fluxjac2(ii,i) /= fluxjac2(ii,i) ) then
!        if( fluxjac2(ii,i) == fluxjac2(ii,i) - 1.0d0 ) then
        ans = ans + 1.0d0
        else

        end if
      end do

      print*, "nonzero elements in each flux jacobian column"
      print*, ii, ans
      end do



      print*, "      "
      print*, "Begin factorisation . . ."
!      call DGETRF(3*ni+9,3*ni+9,fluxjac2,3*ni+9,IPIV,INFO)
      call DGETRF(3*ni,3*ni,hold2,3*ni,IPIV,INFO)
      print*, "Begin system solution . . ."
!      call DGETRS('N',3*ni+9,1,fluxjac2,3*ni+9,IPIV,adjoint,3*ni+9,INFO)
      call DGETRS('T',3*ni,1,hold2,3*ni,IPIV,hold3,3*ni,INFO)
      print*, "System solved."
      print*, "      "


      print*, "       i ::      flux_jac dRdu     *      adjoint vector
     &      =        dJdu"
      print*, "      "
      ans = 0.0d0
      do i=1,3*ni
        print*, i, hold2(i,i), hold3(1,i)
!        print*, i, fluxjac2(i,i), adjoint(1,i), djdu2(1,i)
!        ans = ans + adjoint(1,i)*dRda2(1,i)
      end do


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------



!      print*, "      "
!      print*, "      "
!      print*, "       i ::         adjoint         *         dRda"
!      print*, "      "
!      ans = 0.0d0
!      do i=1,3*ni+9
!        if( adjoint(1,i) /= adjoint(1,i) ) then
!          print*, i, adjoint(1,i), dRda2(1,i), "NaN"
!        else if ( adjoint(1,i) - 1.0d0 == adjoint(1,i) ) then
!          print*, i, adjoint(1,i), dRda2(1,i), "inf"
!        else
!          print*, i, adjoint(1,i), dRda2(1,i)
!          ans = ans + adjoint(1,i)*dRda2(1,i)
!        end if
!
!!        print*, i, adjoint(1,i), dRda2(1,i)
!      end do

      ans = 0.0d0
      do i=1,3*ni
      ans = ans + hold3(1,i)*hold4(1,i)
      end do


      print*, "      "
      print*, "dJda = ", djda
      print*, "DJDa = ", djda - ans



      stop







      stop


!!     debugging statements      
!      print*, "     "
!      i = 1
!      print*, "u  = ", ur(i),um(i),ue(i)
!      print*, "uc = ", flow(1,i+2),flow(1,ni+5+i),flow(1,2*ni+8+i)
!      print*, "R  = ", res1(i),res2(i),res3(i)
!      print*, "Rc = ", residual(1,i+2), residual(1,ni+5+i), 
!     & residual(1,2*ni+8+i)

!      stop
      end subroutine nozzle_main
!
!**********************************************************************
