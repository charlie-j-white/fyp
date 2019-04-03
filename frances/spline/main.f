      program main
!
      integer :: na = 40
      double precision, dimension(1,40) :: alpha
!
      integer :: nx = 150
      integer :: ny = 40
      integer :: nl = 2
      double precision, dimension(0-2:150+2,0-2:40+2) ::
     &   meshX,meshY
!
      integer :: i
      double precision :: AR,pi,frac
!
!
!
!
!     hard code constants
!
      AR = 2.0d0
      pi = 3. 1415926535897932d0
!
!
!
!     set mesh height shape
!
      do i = 1,na
!
      frac = DBLE(i)/(na+1.0d0)
      alpha(1,i) = 1.0d0/AR + (1.0d0 - 1.0d0/AR)*cos(pi*frac)**2.0d0
!      alpha(1,i)=1.0d0/AR+(1.0d0-1.0d0/AR)*cos(2.0d0*pi*frac)**2.0d0
!
      end do
!
!
      alpha(1,10) = alpha(1,10) + 0.1
!
!
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY)
!
!
!
!
!
!
      stop
!
      end program main
