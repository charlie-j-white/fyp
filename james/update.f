!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine update(nx,ny,nl,flow,residual,np,params,meshX,meshY,
     & dtmax)
!
!
      integer :: nx,ny,nl,np
      integer :: i,j,imax
!
      double precision :: gam,CFL,vx,vy,dx,dy,dtmax,SoS,vabs,
     & vxm,vym,dxm,dym
!
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,residual
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres,
     &         r1,r2,r3,r5,
     &         dt
!
!
!
!
      gam = params(1,3)
      CFL = params(1,4)
      dtmax = 10000.00000d0
!
!
!
!
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!
!
!
!
!     compute timestep: SoS, local speed and local time
!
      do j = 1,ny
      do i = 1,nx
!
      pres(i,j) = (gam-1.0d0)*
     & (u5(i,j)-0.5d0*(u2(i,j)**2+u3(i,j)**2)/u1(i,j))
      SoS = DSQRT(gam*pres(i,j)/u1(i,j))
      vx = DABS(u2(i,j)/u1(i,j)) +SoS
      vy = DABS(u3(i,j)/u1(i,j)) +SoS
!
      dx = DMAX1(DABS(meshX(i,j) - meshX(i-1,j)),
     &           DABS(meshX(i-1,j) - meshX(i-2,j)))
      dx = DMAX1(DABS(meshX(i+1,j) - meshX(i,j)), dx)
      dy = DMAX1(DABS(meshY(i,j) - meshY(i,j-1)),
     &           DABS(meshY(i,j-1) - meshY(i,j-2)))
      dy = DMAX1(DABS(meshY(i,j+1) - meshY(i,j)), dy)
!
      dt(i,j) = CFL/(vx/dx + vy/dy)
!
      if (dt(i,j) .LT. dtmax) then
      dtmax = dt(i,j)
      dxm = dx
      dym = dy
      vxm = vx
      vym = vy
      imax = i
      end if
!      
      end do
      end do
!
!
!      print*, dtmax, imax, dxm, dym, vxm, vym
!      print*, "  "
!
!
!
!     do the actual update procedure using first order Euler time
!     integration
!
      do j = 1,ny
      do i = 1,nx
!
      u1(i,j) = u1(i,j) - dtmax*r1(i,j)
      u2(i,j) = u2(i,j) - dtmax*r2(i,j)
      u3(i,j) = u3(i,j) - dtmax*r3(i,j)
      u5(i,j) = u5(i,j) - dtmax*r5(i,j)
!
!      u1(i,j) = u1(i,j) - dt(i,j)*r1(i,j)
!      u2(i,j) = u2(i,j) - dt(i,j)*r2(i,j)
!      u3(i,j) = u3(i,j) - dt(i,j)*r3(i,j)
!      u5(i,j) = u5(i,j) - dt(i,j)*r5(i,j)
!
!
      end do
      end do
!
!
!
!
!
!   
!
!
!
      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
!
!
!
!
      end subroutine update
