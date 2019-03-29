!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine postprocess(nx,ny,nl,flow,u1,u2,u3,u5,meshX,meshY)
      integer :: nx,ny,nl,nh,nb
      integer :: i,j,R
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision :: gam
!      double precision :: 
!
!
!
!
!
!
!
!
!     calculate pressure from flow solution

      gam = 1.4d0

      do i = 0,nx+1
      do j = 1,ny
      pres(i,j) = (gam-1.0d0)*
     & (u5(i,j)-0.5d0*(u2(i,j)**2+u3(i,j)**2)/u1(i,j))
      end do
      end do
!
!
!
!
!
!     write info to file
!
      open(100, file='transect.plt')
      write(100,*) 'VARIABLES = "x" "u1" "u2" "u3" "u5" "pres"'
101   format(6f12.8)
!
!
!      print*, "ny = ", ny
      j = nint((dfloat(ny)+1.0d0)/2.0d0)
!      print*, "transect at", j

      do i = 0,nx
            write(100,101) 
     &                  meshX(i,j),
     &                  0.5*(u1(i,j)+u1(i+1,j)),
     &                  0.5*(u2(i,j)+u2(i+1,j)),
     &                  0.5*(u3(i,j)+u3(i+1,j)),
     &                  0.5*(u5(i,j)+u5(i+1,j)),
     &                  0.5*(pres(i,j)+pres(i+1,j))
      end do
!
!
!
!
!
!
!
!
      end subroutine postprocess
