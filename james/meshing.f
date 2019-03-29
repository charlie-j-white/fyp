!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine meshing(nx,ny,na,nl,alpha,meshX,meshY)
      integer :: nx,ny,na,nl
      integer :: i,j 
      double precision, dimension(1,na) :: alpha
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision :: Lx, Ly, AR, PS, pi
!
!
!
!
      Lx = 1.0d0
      Ly = 1.0d0
      AR = alpha(1,1)
      PS = alpha(1,2)
      pi = 3.141592d0
!
!
      do j = 0-nl,ny+nl
      do i = 0-nl,nx+nl
!
      meshX(i,j) = (i+nl)*Lx/(nx+2.0d0*nl)
      meshY(i,j) = (j+nl)*1.0d0/(ny+2.0d0*nl)*2.0d0 - 1.0d0
      meshY(i,j) = meshY(i,j)*Ly*
     & (1.0d0/AR + (1.0d0-1.0d0/AR)*cos(pi*(i+nl)/(nx+2.0d0*nl))**2.0d0)
!
      end do
!      print*, meshX(:,j)
      end do
!
!      print*, "     "
!
!      do j = 0-nl,ny+nl
!      print*, meshY(:,j)
!      end do
!
!
!
      end subroutine meshing
