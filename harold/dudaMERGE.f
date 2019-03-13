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
     
      subroutine cost(J,arearat,ur,um,ue,ni)

      implicit none

      integer, intent(in) :: ni
      double precision, intent(in) :: arearat
      double precision, intent(out) :: J
      double precision, dimension(-1:ni+1), intent(in) :: ur,um,ue

      integer :: ii
      double precision, dimension(-1:ni+2) :: x,ds,dx,area,dadx
      double precision, dimension(-1:ni+1) :: pres


      J = 0




!
! for the cost function, integrate normal pressure 
! on both upper and lower nozzle surfaces
!-----------------------------------------------------------------------


      ! create mesh
      ! useful inputs: ni, arearat        useful outputs: dx
      call griddata(ni,x,ds,dx,area,dadx,arearat)

      do ii = 1,ni

             ! calculate pressure
             pres(ii)= (1.4d0-1.0d0)*(ue(ii)-0.5d0*um(ii)*um(ii)/ur(ii))

             ! integrate pressure
             J = J + 2*pres(ii)*dx(ii)*dsqrt(1+dadx(ii)*dadx(ii))
      end do



      end subroutine cost
!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine griddata(ni,x,ds,dx,area,dadx,arearat)
!
      IMPLICIT NONE
!
      INTEGER:: ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: x,ds,dx,area,dadx,xc
!      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: xc
      DOUBLE PRECISION:: arearat,pi
      INTEGER:: i
!
!  CALCULATE CELL LENGTHS AND DISTANCES BETWEEN CELL CENTRES
!
!      ALLOCATE(xc(-1:ni+2))
!
      do i=0,ni
        x(i)=dfloat(i)*1.0d0/dfloat(ni) 
      enddo
!
      do i=1,ni
        xc(i)=0.5d0*(x(i)+x(i-1))
      enddo
        xc(0)=2.0d0*x(0)-xc(1)
        xc(ni+1)=2.0d0*x(ni)-xc(ni)
        x(-1)=2.0d0*x(0)-x(1)
        x(ni+1)=2.0d0*x(ni)-x(ni-1)
      do i=1,ni+1
        ds(i)=xc(i)-xc(i-1)
      enddo
      ds(0)=ds(1)
      ds(ni+2)=ds(ni+1)
      do i=0,ni+1
        dx(i)=x(i)-x(i-1)
      enddo
!
      pi=4.0d0*datan(1.0d0)
!
      do i=0,ni
        if(x(i).lt.0.1d0) then 
          area(i)=arearat
          dadx(i)=0.0d0
        elseif(x(i).gt.0.9d0) then 
          area(i)=arearat
          dadx(i)=0.0d0
        else
       area(i)=1.0d0+(arearat-1.0d0)*((dcos(1.25d0*pi*(x(i)-0.1d0)))**2)
         dadx(i)=-2.5d0*pi*(arearat-1.0d0)*dcos(1.25d0*pi*(x(i)-0.1d0))
     &             *dsin(1.25d0*pi*(x(i)-0.1d0))
        endif
      enddo
!
!      DEALLOCATE(xc)
!
      return 
      end subroutine griddata
!
!**********************************************************************
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
      INTEGER:: i,ni,nt,nstages,maxit,ns
      DOUBLE PRECISION:: J,jd,arearat
      DOUBLE PRECISION, DIMENSION(-1:ni+1):: 
     &             x,ur,um,ue,pres,ds,dx,res1,res2,res3,
     &             error1,error2,error3,dt,area,dadx,ur0,um0,ue0,
     &             vec0,vec1 
      DOUBLE PRECISION, DIMENSION(0:6):: alp
      DOUBLE PRECISION, DIMENSION(1:3):: djdu
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
      ! primal output of cost function
      call cost(J,arearat,ur,um,ue,ni)
      print*, "   J = ", J

      ! create seed vectors
      do i=-1,ni+1
         vec0(i)=0.0d0
         vec1(i)=1.0d0
      end do




      ! get dJ/da
!     call cost_d(j,jd,arearat,arearatd,ur,urd,um,umd,ue,ued,ni)
      call cost_d(j,jd,arearat,1.0d0,ur,vec0,um,vec0,ue,vec0,ni)
      print*, "djda = ", jd

      ! get dJ/du
      call cost_d(j,jd,arearat,0.0d0,ur,vec1,um,vec0,ue,vec0,ni)
      djdu(1) = jd
      call cost_d(j,jd,arearat,0.0d0,ur,vec0,um,vec1,ue,vec0,ni)
      djdu(2) = jd
      call cost_d(j,jd,arearat,0.0d0,ur,vec0,um,vec0,ue,vec1,ni)
      djdu(3) = jd
      print*, "djdu = ", djdu

      ! get du/da
      print*, "duda = "


      ! calculate derivatives relative to the residual
      call jRES(ni,ur,um,ue,pk,res1,res2,res3,arearat)
      print*, "    "
      print*, "dRdu = "
      print*, "dRda = "

      ! form total product to compare with finite difference
      print*, "    "
      print*, "DJDa = "












!      stop
      end subroutine nozzle_main
!
!**********************************************************************
!***********************************************************************
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
      subroutine resid(ni,ur,um,ue,dx,ds,area,pk,res1,res2,res3)
!
      IMPLICIT NONE
!
      INTEGER:: ni,nj
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,ue,dx,ds,area
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: res1,res2,res3 
      DOUBLE PRECISION:: gam,pk,drr,drl,dmr,dml,der,del,sr,sm,se
      DOUBLE PRECISION:: vr,vl,uml,umr,urr,url,uer,uel,pr,cr,pl,cl
      DOUBLE PRECISION:: pmr,pml,smr,sml
!      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: resf1,resf2,resf3
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: resf1,
     &                              resf2,resf3
      INTEGER:: i
!
      res1=0.0d0
      res2=0.0d0
      res3=0.0d0
      gam=1.4d0
!
!      ALLOCATE(resf1(-1:ni+2))
!      ALLOCATE(resf2(-1:ni+2))
!      ALLOCATE(resf3(-1:ni+2))
      resf1=0.0d0
      resf2=0.0d0
      resf3=0.0d0
!
      do i=0,ni+1
        drr=(ur(i+1)-ur(i))/ds(i+1)
        drl=(ur(i)-ur(i-1))/ds(i)
      sr=(2.0d0*drr*drl + 0.0000001d0)/(drr*drr + drl*drl + 0.0000001d0)
        dmr=(um(i+1)-um(i))/ds(i+1)
        dml=(um(i)-um(i-1))/ds(i)
      sm=(2.0d0*dmr*dml + 0.0000001d0)/(dmr*dmr + dml*dml + 0.0000001d0)
        der=(ue(i+1)-ue(i))/ds(i+1)
        del=(ue(i)-ue(i-1))/ds(i)
      se=(2.0d0*der*del + 0.0000001d0)/(der*der + del*del + 0.0000001d0)
        if(sr.lt.0.0d0) sr=0.0d0
        if(sm.lt.0.0d0) sm=0.0d0
        if(se.lt.0.0d0) se=0.0d0
        urr=ur(i)+0.5d0*sr*((1.0d0-pk)*drl + (1.0d0+pk)*drr)*0.5d0*dx(i)
        url=ur(i)-0.5d0*sr*((1.0d0-pk)*drr + (1.0d0+pk)*drl)*0.5d0*dx(i)
        umr=um(i)+0.5d0*sm*((1.0d0-pk)*dml + (1.0d0+pk)*dmr)*0.5d0*dx(i)
        uml=um(i)-0.5d0*sm*((1.0d0-pk)*dmr + (1.0d0+pk)*dml)*0.5d0*dx(i)
        uer=ue(i)+0.5d0*se*((1.0d0-pk)*del + (1.0d0+pk)*der)*0.5d0*dx(i)
        uel=ue(i)-0.5d0*se*((1.0d0-pk)*der + (1.0d0+pk)*del)*0.5d0*dx(i)
        vr=umr/urr
        vl=uml/url
        pr=(gam-1.0d0)*(uer-0.5d0*vr*vr*urr)
        cr=dsqrt(gam*pr/urr)
        pmr=vr/cr
        pl=(gam-1.0d0)*(uel-0.5d0*vl*vl*url)
        cl=dsqrt(gam*pl/url)
        pml=vl/cr
!
        if(pmr.gt.1.0d0) then 
          resf1(i+1)=urr*vr
          resf2(i+1)=pr+urr*vr*vr
          resf3(i+1)=(pr+uer)*vr
        elseif(pmr.lt.-1.0d0) then 
          resf1(i+1)=0.0d0
          resf2(i+1)=0.0d0
          resf3(i+1)=0.0d0
        else
          smr=0.25d0*urr*cr*((pmr+1.0d0)**2)
          resf1(i+1)=smr
          resf2(i+1)=smr*((-vr+2.0d0*cr)/gam + vr) 
          resf3(i+1)=smr*(((gam-1.0d0)*vr+2.0d0*cr)**2)/
     &  (2.0d0*((gam*gam)-1.0d0))
        endif
        if(pml.lt.-1.0d0) then 
          resf1(i)=resf1(i)+url*vl
          resf2(i)=resf2(i)+pl+url*vl*vl
          resf3(i)=resf3(i)+(pl+uel)*vl
        elseif(pml.gt.1.0d0) then 
          resf1(i)=resf1(i)+0.0d0
          resf2(i)=resf2(i)+0.0d0
          resf3(i)=resf3(i)+0.0d0
        else
          sml=-0.25d0*url*cl*((pml-1.0d0)**2)
          resf1(i)=resf1(i)+sml
          resf2(i)=resf2(i)+sml*((-vl-2.0d0*cl)/gam + vl) 
          resf3(i)=resf3(i)+sml*(((gam-1.0d0)*vl-2.0d0*cl)**2)/
     &  (2.0d0*((gam*gam)-1.0d0))
        endif
      enddo
!
      do i=1,ni
        res1(i)=(resf1(i+1)*area(i)-resf1(i)*area(i-1))/dx(i)
        res2(i)=(resf2(i+1)*area(i)-resf2(i)*area(i-1))/dx(i)
        res3(i)=(resf3(i+1)*area(i)-resf3(i)*area(i-1))/dx(i)
      enddo
!
!      DEALLOCATE(resf1,resf2,resf3) 
!
      return
      end subroutine resid
!
!**********************************************************************
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
