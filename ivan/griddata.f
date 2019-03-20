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
