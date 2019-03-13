     
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
