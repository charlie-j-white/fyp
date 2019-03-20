!***********************************************************************
!
!                               j R E S 
!      
!  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
!  LAST UPDATE: NOVEMBER 2014.
!***********************************************************************
!
!     subroutine jRES(ni,ur,um,ue,dx,ds,area,pk,res1,res2,res3)
      subroutine jRES(ni,pk,pmachin,prat,arearat,residual,flow)
!
      IMPLICIT NONE
!
      INTEGER:: ni,nj
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: ur,um,ue,dx,ds,area,pres
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: res1,res2,res3 
      DOUBLE PRECISION:: gam,pk,drr,drl,dmr,dml,der,del,sr,sm,se
      DOUBLE PRECISION:: vr,vl,uml,umr,urr,url,uer,uel,pr,cr,pl,cl
      DOUBLE PRECISION:: pmr,pml,smr,sml
      DOUBLE PRECISION:: arearat,pmachin,prat
!      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: resf1,resf2,resf3
      DOUBLE PRECISION, DIMENSION(-1:ni+2):: resf1,
     &                              resf2,resf3,x,dadx
      INTEGER:: i
      DOUBLE PRECISION, DIMENSION(1,1:3*ni+9):: residual,flow
!

      call griddata(ni,x,ds,dx,area,dadx,arearat)

! assign flow solution to variables


      do i=1,ni+3
      ur(i-2) = flow(1,i)
      end do

      do i=1,ni+3
      um(i-2) = flow(1,ni+3+i)
      end do

      do i=1,ni+3
      ue(i-2) = flow(1,2*ni+6+i)
      end do

      do i=-1,ni+2
      pres(i) = 0.4d0*(ue(i)-0.5d0*um(i)*um(i)/ur(i))
      end do

      call boundaries(ni,pmachin,prat,ur,um,ue,pres)



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

!     assign residual values to vector

      do i=1,ni+3
      residual(1,i) = res1(i-2)
      end do

      do i=1,ni+3
      residual(1,ni+3+i) = res2(i-2)
      end do

      do i=1,ni+3
      residual(1,2*ni+6+i) = res3(i-2)
      end do






!
!      DEALLOCATE(resf1,resf2,resf3) 
!
      return
      end subroutine jRES
!
!**********************************************************************
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
