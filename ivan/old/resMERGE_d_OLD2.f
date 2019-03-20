C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.14 (r7260) - 18 Jan 2019 10:11
C
C  Differentiation of jres in forward (tangent) mode:
C   variations   of useful results: residual
C   with respect to varying inputs: flow arearat
C   RW status of diff variables: flow:in residual:out arearat:in
C***********************************************************************
C
C                               j R E S 
C      
C  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
C  LAST UPDATE: NOVEMBER 2014.
C***********************************************************************
C
C     subroutine jRES(ni,ur,um,ue,dx,ds,area,pk,res1,res2,res3)
      SUBROUTINE JRES_D(ni, pk, arearat, arearatd, residual, residuald, 
     +                  flow, flowd)
      IMPLICIT NONE
C
      INTEGER ni, nj
      DOUBLE PRECISION, DIMENSION(-1:ni+2) :: ur, um, ue, dx, ds, area
      DOUBLE PRECISION urd(-1:ni+2), umd(-1:ni+2), ued(-1:ni+2), aread(-
     +                 1:ni+2)
      DOUBLE PRECISION, DIMENSION(-1:ni+2) :: res1, res2, res3
      DOUBLE PRECISION res1d(-1:ni+2), res2d(-1:ni+2), res3d(-1:ni+2)
      DOUBLE PRECISION gam, pk, drr, drl, dmr, dml, der, del, sr, sm, se
      DOUBLE PRECISION drrd, drld, dmrd, dmld, derd, deld, srd, smd, sed
      DOUBLE PRECISION vr, vl, uml, umr, urr, url, uer, uel, pr, cr, pl
     +                 , cl
      DOUBLE PRECISION vrd, vld, umld, umrd, urrd, urld, uerd, ueld, prd
     +                 , crd, pld, cld
      DOUBLE PRECISION pmr, pml, smr, sml
      DOUBLE PRECISION pmrd, pmld, smrd, smld
      DOUBLE PRECISION arearat
      DOUBLE PRECISION arearatd
C      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: resf1,resf2,resf3
      DOUBLE PRECISION, DIMENSION(-1:ni+2) :: resf1, resf2, resf3, x, 
     +                                        dadx
      DOUBLE PRECISION resf1d(-1:ni+2), resf2d(-1:ni+2), resf3d(-1:ni+2)
      INTEGER i
      DOUBLE PRECISION, DIMENSION(1, 3*ni+9) :: residual, flow
      DOUBLE PRECISION residuald(1, 3*ni+9), flowd(1, 3*ni+9)
      INTRINSIC DSQRT
      DOUBLE PRECISION arg1
      DOUBLE PRECISION arg1d
C
C
      CALL GRIDDATA_DD(ni, x, ds, dx, area, aread, dadx, arearat, 
     +                arearatd)
      urd = 0.D0
C
C assign flow solution to variables
C
C
      DO i=1,ni+3
        urd(i-2) = flowd(1, i)
        ur(i-2) = flow(1, i)
      ENDDO
      umd = 0.D0
C
      DO i=1,ni+3
        umd(i-2) = flowd(1, ni+3+i)
        um(i-2) = flow(1, ni+3+i)
      ENDDO
      ued = 0.D0
C
      DO i=1,ni+3
        ued(i-2) = flowd(1, 2*ni+6+i)
        ue(i-2) = flow(1, 2*ni+6+i)
      ENDDO
C
C
C
C
      res1 = 0.0d0
      res2 = 0.0d0
      res3 = 0.0d0
      gam = 1.4d0
C
C      ALLOCATE(resf1(-1:ni+2))
C      ALLOCATE(resf2(-1:ni+2))
C      ALLOCATE(resf3(-1:ni+2))
      resf1 = 0.0d0
      resf2 = 0.0d0
      resf3 = 0.0d0
      resf1d = 0.D0
      resf2d = 0.D0
      resf3d = 0.D0
C
      DO i=0,ni+1
        drrd = (urd(i+1)-urd(i))/ds(i+1)
        drr = (ur(i+1)-ur(i))/ds(i+1)
        drld = (urd(i)-urd(i-1))/ds(i)
        drl = (ur(i)-ur(i-1))/ds(i)
        srd = (2.0d0*(drrd*drl+drr*drld)*(drr*drr+drl*drl+0.0000001d0)-(
     +    2.0d0*drr*drl+0.0000001d0)*(drrd*drr+drr*drrd+drld*drl+drl*
     +    drld))/(drr*drr+drl*drl+0.0000001d0)**2
        sr = (2.0d0*drr*drl+0.0000001d0)/(drr*drr+drl*drl+0.0000001d0)
        dmrd = (umd(i+1)-umd(i))/ds(i+1)
        dmr = (um(i+1)-um(i))/ds(i+1)
        dmld = (umd(i)-umd(i-1))/ds(i)
        dml = (um(i)-um(i-1))/ds(i)
        smd = (2.0d0*(dmrd*dml+dmr*dmld)*(dmr*dmr+dml*dml+0.0000001d0)-(
     +    2.0d0*dmr*dml+0.0000001d0)*(dmrd*dmr+dmr*dmrd+dmld*dml+dml*
     +    dmld))/(dmr*dmr+dml*dml+0.0000001d0)**2
        sm = (2.0d0*dmr*dml+0.0000001d0)/(dmr*dmr+dml*dml+0.0000001d0)
        derd = (ued(i+1)-ued(i))/ds(i+1)
        der = (ue(i+1)-ue(i))/ds(i+1)
        deld = (ued(i)-ued(i-1))/ds(i)
        del = (ue(i)-ue(i-1))/ds(i)
        sed = (2.0d0*(derd*del+der*deld)*(der*der+del*del+0.0000001d0)-(
     +    2.0d0*der*del+0.0000001d0)*(derd*der+der*derd+deld*del+del*
     +    deld))/(der*der+del*del+0.0000001d0)**2
        se = (2.0d0*der*del+0.0000001d0)/(der*der+del*del+0.0000001d0)
        IF (sr .LT. 0.0d0) THEN
          sr = 0.0d0
          srd = 0.D0
        END IF
        IF (sm .LT. 0.0d0) THEN
          sm = 0.0d0
          smd = 0.D0
        END IF
        IF (se .LT. 0.0d0) THEN
          se = 0.0d0
          sed = 0.D0
        END IF
        urrd = urd(i) + 0.5d0**2*dx(i)*(srd*((1.0d0-pk)*drl+(1.0d0+pk)*
     +    drr)+sr*((1.0d0-pk)*drld+(1.0d0+pk)*drrd))
        urr = ur(i) + 0.5d0*sr*((1.0d0-pk)*drl+(1.0d0+pk)*drr)*0.5d0*dx(
     +    i)
        urld = urd(i) - 0.5d0**2*dx(i)*(srd*((1.0d0-pk)*drr+(1.0d0+pk)*
     +    drl)+sr*((1.0d0-pk)*drrd+(1.0d0+pk)*drld))
        url = ur(i) - 0.5d0*sr*((1.0d0-pk)*drr+(1.0d0+pk)*drl)*0.5d0*dx(
     +    i)
        umrd = umd(i) + 0.5d0**2*dx(i)*(smd*((1.0d0-pk)*dml+(1.0d0+pk)*
     +    dmr)+sm*((1.0d0-pk)*dmld+(1.0d0+pk)*dmrd))
        umr = um(i) + 0.5d0*sm*((1.0d0-pk)*dml+(1.0d0+pk)*dmr)*0.5d0*dx(
     +    i)
        umld = umd(i) - 0.5d0**2*dx(i)*(smd*((1.0d0-pk)*dmr+(1.0d0+pk)*
     +    dml)+sm*((1.0d0-pk)*dmrd+(1.0d0+pk)*dmld))
        uml = um(i) - 0.5d0*sm*((1.0d0-pk)*dmr+(1.0d0+pk)*dml)*0.5d0*dx(
     +    i)
        uerd = ued(i) + 0.5d0**2*dx(i)*(sed*((1.0d0-pk)*del+(1.0d0+pk)*
     +    der)+se*((1.0d0-pk)*deld+(1.0d0+pk)*derd))
        uer = ue(i) + 0.5d0*se*((1.0d0-pk)*del+(1.0d0+pk)*der)*0.5d0*dx(
     +    i)
        ueld = ued(i) - 0.5d0**2*dx(i)*(sed*((1.0d0-pk)*der+(1.0d0+pk)*
     +    del)+se*((1.0d0-pk)*derd+(1.0d0+pk)*deld))
        uel = ue(i) - 0.5d0*se*((1.0d0-pk)*der+(1.0d0+pk)*del)*0.5d0*dx(
     +    i)
        vrd = (umrd*urr-umr*urrd)/urr**2
        vr = umr/urr
        vld = (umld*url-uml*urld)/url**2
        vl = uml/url
        prd = (gam-1.0d0)*(uerd-0.5d0*((vrd*vr+vr*vrd)*urr+vr**2*urrd))
        pr = (gam-1.0d0)*(uer-0.5d0*vr*vr*urr)
        arg1d = (gam*prd*urr-gam*pr*urrd)/urr**2
        arg1 = gam*pr/urr
        IF (arg1 .EQ. 0.0) THEN
          crd = 0.D0
        ELSE
          crd = arg1d/(2.D0*DSQRT(arg1))
        END IF
        cr = DSQRT(arg1)
        pmrd = (vrd*cr-vr*crd)/cr**2
        pmr = vr/cr
        pld = (gam-1.0d0)*(ueld-0.5d0*((vld*vl+vl*vld)*url+vl**2*urld))
        pl = (gam-1.0d0)*(uel-0.5d0*vl*vl*url)
        arg1d = (gam*pld*url-gam*pl*urld)/url**2
        arg1 = gam*pl/url
        IF (arg1 .EQ. 0.0) THEN
          cld = 0.D0
        ELSE
          cld = arg1d/(2.D0*DSQRT(arg1))
        END IF
        cl = DSQRT(arg1)
        pmld = (vld*cr-vl*crd)/cr**2
        pml = vl/cr
C
        IF (pmr .GT. 1.0d0) THEN
          resf1d(i+1) = urrd*vr + urr*vrd
          resf1(i+1) = urr*vr
          resf2d(i+1) = prd + (urrd*vr+urr*vrd)*vr + urr*vr*vrd
          resf2(i+1) = pr + urr*vr*vr
          resf3d(i+1) = (prd+uerd)*vr + (pr+uer)*vrd
          resf3(i+1) = (pr+uer)*vr
        ELSE IF (pmr .LT. -1.0d0) THEN
          resf1d(i+1) = 0.D0
          resf1(i+1) = 0.0d0
          resf2d(i+1) = 0.D0
          resf2(i+1) = 0.0d0
          resf3d(i+1) = 0.D0
          resf3(i+1) = 0.0d0
        ELSE
          smrd = 0.25d0*((urrd*cr+urr*crd)*(pmr+1.0d0)**2+urr*cr*2*(pmr+
     +      1.0d0)*pmrd)
          smr = 0.25d0*urr*cr*(pmr+1.0d0)**2
          resf1d(i+1) = smrd
          resf1(i+1) = smr
          resf2d(i+1) = smrd*((-vr+2.0d0*cr)/gam+vr) + smr*((2.0d0*crd-
     +      vrd)/gam+vrd)
          resf2(i+1) = smr*((-vr+2.0d0*cr)/gam+vr)
          resf3d(i+1) = (smrd*((gam-1.0d0)*vr+2.0d0*cr)**2+smr*2*((gam-
     +      1.0d0)*vr+2.0d0*cr)*((gam-1.0d0)*vrd+2.0d0*crd))/(2.0d0*(gam
     +      *gam-1.0d0))
          resf3(i+1) = smr*((gam-1.0d0)*vr+2.0d0*cr)**2/(2.0d0*(gam*gam-
     +      1.0d0))
        END IF
        IF (pml .LT. -1.0d0) THEN
          resf1d(i) = resf1d(i) + urld*vl + url*vld
          resf1(i) = resf1(i) + url*vl
          resf2d(i) = resf2d(i) + pld + (urld*vl+url*vld)*vl + url*vl*
     +      vld
          resf2(i) = resf2(i) + pl + url*vl*vl
          resf3d(i) = resf3d(i) + (pld+ueld)*vl + (pl+uel)*vld
          resf3(i) = resf3(i) + (pl+uel)*vl
        ELSE IF (pml .GT. 1.0d0) THEN
          resf1(i) = resf1(i) + 0.0d0
          resf2(i) = resf2(i) + 0.0d0
          resf3(i) = resf3(i) + 0.0d0
        ELSE
          smld = -(0.25d0*((urld*cl+url*cld)*(pml-1.0d0)**2+url*cl*2*(
     +      pml-1.0d0)*pmld))
          sml = -(0.25d0*url*cl*(pml-1.0d0)**2)
          resf1d(i) = resf1d(i) + smld
          resf1(i) = resf1(i) + sml
          resf2d(i) = resf2d(i) + smld*((-vl-2.0d0*cl)/gam+vl) + sml*((-
     +      vld-2.0d0*cld)/gam+vld)
          resf2(i) = resf2(i) + sml*((-vl-2.0d0*cl)/gam+vl)
          resf3d(i) = resf3d(i) + (smld*((gam-1.0d0)*vl-2.0d0*cl)**2+sml
     +      *2*((gam-1.0d0)*vl-2.0d0*cl)*((gam-1.0d0)*vld-2.0d0*cld))/(
     +      2.0d0*(gam*gam-1.0d0))
          resf3(i) = resf3(i) + sml*((gam-1.0d0)*vl-2.0d0*cl)**2/(2.0d0*
     +      (gam*gam-1.0d0))
        END IF
      ENDDO
      res1d = 0.D0
      res2d = 0.D0
      res3d = 0.D0
C
      DO i=1,ni
        res1d(i) = (resf1d(i+1)*area(i)+resf1(i+1)*aread(i)-resf1d(i)*
     +    area(i-1)-resf1(i)*aread(i-1))/dx(i)
        res1(i) = (resf1(i+1)*area(i)-resf1(i)*area(i-1))/dx(i)
        res2d(i) = (resf2d(i+1)*area(i)+resf2(i+1)*aread(i)-resf2d(i)*
     +    area(i-1)-resf2(i)*aread(i-1))/dx(i)
        res2(i) = (resf2(i+1)*area(i)-resf2(i)*area(i-1))/dx(i)
        res3d(i) = (resf3d(i+1)*area(i)+resf3(i+1)*aread(i)-resf3d(i)*
     +    area(i-1)-resf3(i)*aread(i-1))/dx(i)
        res3(i) = (resf3(i+1)*area(i)-resf3(i)*area(i-1))/dx(i)
      ENDDO
      residuald = 0.D0
C
C     assign residual values to vector
C
      DO i=1,ni+3
        residuald(1, i) = res1d(i-2)
        residual(1, i) = res1(i-2)
      ENDDO
C
      DO i=1,ni+3
        residuald(1, ni+3+i) = res2d(i-2)
        residual(1, ni+3+i) = res2(i-2)
      ENDDO
C
      DO i=1,ni+3
        residuald(1, 2*ni+6+i) = res3d(i-2)
        residual(1, 2*ni+6+i) = res3(i-2)
      ENDDO
C
C
C
C
C
C
C
C      DEALLOCATE(resf1,resf2,resf3) 
C
      RETURN
      END

C  Differentiation of griddata in forward (tangent) mode:
C   variations   of useful results: area
C   with respect to varying inputs: arearat
C
C**********************************************************************
C***********************************************************************
C  CHRISTIAN ALLEN. UNIVERSITY OF BRISTOL.
C  LAST UPDATE: NOVEMBER 2014.
C***********************************************************************
C
      SUBROUTINE GRIDDATA_DD(ni, x, ds, dx, area, aread, dadx, arearat, 
     +                      arearatd)
      IMPLICIT NONE
C
      INTEGER ni
      DOUBLE PRECISION, DIMENSION(-1:ni+2) :: x, ds, dx, area, dadx, xc
      DOUBLE PRECISION aread(-1:ni+2)
C      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: xc
      DOUBLE PRECISION arearat, pi
      DOUBLE PRECISION arearatd
      INTEGER i
      INTRINSIC DFLOAT
      INTRINSIC DATAN
      INTRINSIC DCOS
      INTRINSIC DSIN
      DOUBLE PRECISION arg1
      DOUBLE PRECISION arg2
      INTEGER ii1
C
C  CALCULATE CELL LENGTHS AND DISTANCES BETWEEN CELL CENTRES
C
C      ALLOCATE(xc(-1:ni+2))
C
      DO i=0,ni
        x(i) = DFLOAT(i)*1.0d0/DFLOAT(ni)
      ENDDO
C
      DO i=1,ni
        xc(i) = 0.5d0*(x(i)+x(i-1))
      ENDDO
      xc(0) = 2.0d0*x(0) - xc(1)
      xc(ni+1) = 2.0d0*x(ni) - xc(ni)
      x(-1) = 2.0d0*x(0) - x(1)
      x(ni+1) = 2.0d0*x(ni) - x(ni-1)
      DO i=1,ni+1
        ds(i) = xc(i) - xc(i-1)
      ENDDO
      ds(0) = ds(1)
      ds(ni+2) = ds(ni+1)
      DO i=0,ni+1
        dx(i) = x(i) - x(i-1)
      ENDDO
C
      pi = 4.0d0*DATAN(1.0d0)
      DO ii1=-1,ni+2
        aread(ii1) = 0.D0
      ENDDO
C
      DO i=0,ni
        IF (x(i) .LT. 0.1d0) THEN
          aread(i) = arearatd
          area(i) = arearat
          dadx(i) = 0.0d0
        ELSE IF (x(i) .GT. 0.9d0) THEN
          aread(i) = arearatd
          area(i) = arearat
          dadx(i) = 0.0d0
        ELSE
          arg1 = 1.25d0*pi*(x(i)-0.1d0)
          aread(i) = DCOS(arg1)**2*arearatd
          area(i) = 1.0d0 + (arearat-1.0d0)*DCOS(arg1)**2
          arg1 = 1.25d0*pi*(x(i)-0.1d0)
          arg2 = 1.25d0*pi*(x(i)-0.1d0)
          dadx(i) = -(2.5d0*pi*(arearat-1.0d0)*DCOS(arg1)*DSIN(arg2))
        END IF
      ENDDO
C
C      DEALLOCATE(xc)
C
      RETURN
      END
C
C**********************************************************************

