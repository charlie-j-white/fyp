      program main
!
      implicit none
!
      integer, external :: colour
!
!
!
      integer :: i,j
      integer :: cerr
      double precision :: R
!
      integer :: nw = 500
      double precision, dimension(500,500) :: fluxjac
!
!
      print*, "Call from the Fortran program"
!
!
      do i = 1,nw
      do j = 1,nw
!
      fluxjac(i,j) = 0.0d0
!
      if (i .EQ. j) then
        fluxjac(i,j) = 1.0d0
      end if
!
      end do
      end do
!
!
      fluxjac(20,60) = 1.0d0
      fluxjac(70,30) = 1.0d0
!
!
!
!
!
      cerr = colour(nw,fluxjac)
!
!
!
      end program main
