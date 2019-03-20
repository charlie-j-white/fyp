     
      subroutine cost(J,arearat,flow,ni)

      implicit none

      integer, intent(in) :: ni
      double precision, intent(in) :: arearat
      double precision, intent(out) :: J
      double precision, dimension(1,1:3*ni+9), intent(in) :: flow

      integer :: ii, i
      double precision, dimension(-1:ni+2) :: x,ds,dx,area,dadx
      double precision, dimension(-1:ni+1) :: pres
      double precision, dimension(-1:ni+1) :: ur,um,ue


      J = 0.0d0

!     assign variables from concatenation

      do i=1,ni+3
      ur(i-2) = flow(1,i)
      end do

      do i=1,ni+3
      um(i-2) = flow(1,ni+3+i)
      end do

      do i=1,ni+3
      ue(i-2) = flow(1,2*ni+6+i)
      end do





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
             J = J +2.0d0*pres(ii)*dx(ii)*dsqrt(1.0d0+dadx(ii)*dadx(ii))
      end do



      end subroutine cost
