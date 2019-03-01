


     
      subroutine jcost(N, J, pres, dx)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: dx
      double precision, intent(out) :: J
      double precision, dimension(0:N+1), intent(in) :: pres

      integer :: ii


      J = 0




!
! for the cost function, integrate vertical pressure (direction normal
! to the gross flow direction) on both upper and lower nozzle surfaces
!-----------------------------------------------------------------------

      do ii = 1,N
             J = J + 2*dx*pres(ii)
      end do



      end subroutine jcost
