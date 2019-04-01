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
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
!
!
!
!
!
!     calculate pressure from flow solution
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
      gam = 1.4d0
!
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
      j = nint((dfloat(ny)+1.0d0)/2.0d0)

      do i = 0,nx
            write(100,101) 
     &                  meshX(i,j),
     &                  0.5d0*(u1(i,j)+u1(i+1,j)),
     &                  0.5d0*(u2(i,j)+u2(i+1,j)),
     &                  0.5d0*(u3(i,j)+u3(i+1,j)),
     &                  0.5d0*(u5(i,j)+u5(i+1,j)),
     &                  0.5d0*(pres(i,j)+pres(i+1,j))
      end do
!
!      close(100)
!
!
!
!
!
!
!
!
!
!
!     tecplot formatting
!
      open(300, file='solution.plt')
      write(300,*) 'VARIABLES = "X" "Y" "RHO" "U" "V" "pres"'
      write(300,*) 'Zone I=',ny+1,', J=',nx-1,', F=POINT'
301   format(6f12.8)
!
      print*, "Start file write . . ."
!
      do i = 1,nx-1
      do j = 0,ny
!
      write(300,301)
     &     meshX(i,j),
     &     meshY(i,j),
     &     0.25d0*(u1(i,j)+u1(i+1,j)+u1(i,j+1)+u1(i+1,j+1)),
     &     0.25d0*(u2(i,j)/u1(i,j)+u2(i+1,j)/u1(i+1,j)
     &         +   u2(i,j+1)/u1(i,j+1)+u2(i+1,j+1)/u1(i+1,j+1)),
     &     0.25d0*(u3(i,j)/u1(i,j)+u3(i+1,j)/u1(i+1,j)
     &         +   u3(i,j+1)/u1(i,j+1)+u3(i+1,j+1)/u1(i+1,j+1)),
     &     0.25d0*(pres(i,j)+pres(i+1,j)+pres(i,j+1)+pres(i+1,j+1))
!
      end do
      end do
!
      close(300)
      print*, "Done."
!
!
!
!
!
!
!
      end subroutine postprocess
