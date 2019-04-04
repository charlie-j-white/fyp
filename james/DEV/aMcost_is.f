!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine cost_is(nx,ny,nl,na,np,params,flow,alpha,Jcost)
!
!
      integer :: nx,ny,nl,na,np
      integer :: i,j
!
      double precision :: ds, Jcost
!
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,pres
!
!
!
!
!
!
      call split_fwd2(nx,ny,nl,flow,u1,u2,u3,u5)
!
      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
      call meshing2(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!
!
!     cost function - just integrate pressure along walls. assume line
!     of symmetry
!
      do i = 1,nx
!
      ds = (meshX(i,ny)+meshX(i-1,ny))**2.0d0+
     &     (meshY(i,ny)+meshY(i-1,ny))**2.0d0
      ds = DSQRT(ds)
!
      Jcost = Jcost + 2*ds*pres(i,ny)
!      
      end do
!
!
!
!
!
!
!
!
      end subroutine cost_is
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine split_fwd2(nx,ny,nl,flow,u1,u2,u3,u5)
      integer :: nx,ny,nl,nh,nb
      integer :: i,j,R
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)), intent(in) 
     & :: flow
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5
!      double precision :: 
!
!
!                                               1  2  3  4  5
!       1 2 3 4 5 6 7 ... 13 14 15     ==>      6  7  8  9  10
!                                               11 12 13 14 15
!
!
!
!
      nh = 2*nl*(nx+ny)
      nb = nx*ny
!
!
!
!     seed flow vector with test values
!
!      do i = 1,4*nh+4*nb
!      flow(1,i) = i
!      end do
!
!
!
!
!
!
!
!
!     assign first flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 0
      u1(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = nl*nx
      u1(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx
      u1(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx+nl*ny
      u1(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign second flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 1*nh+0
      u2(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 1*nh+nl*nx
      u2(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx
      u2(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx+nl*ny
      u2(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign third flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 2*nh+0
      u3(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 2*nh+nl*nx
      u3(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx
      u3(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx+nl*ny
      u3(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!     assign fifth flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 3*nh+0
      u5(i,1-j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 3*nh+nl*nx
      u5(i,ny+j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx
      u5(1-j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx+nl*ny
      u5(nx+j,i) = flow(1,R+(j-1)*ny+i)
      end do
      end do
!
!
!
!
!     assign body cells
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+0*nb
      u1(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+1*nb
      u2(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+2*nb
      u3(i,j) = flow(1,R+(j-1)*nx+i)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+3*nb
      u5(i,j) = flow(1,R+(j-1)*nx+i)
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
!      print*, "     "
!      print*, "------------------------- fwd --------------------------"
!      print*, "     "
!      print*, flow
!!
!      call debug_cell(nx,ny,nl,u1)
!      call debug_cell(nx,ny,nl,u2)
!      call debug_cell(nx,ny,nl,u3)
!      call debug_cell(nx,ny,nl,u5)
!!
!      print*, "--------------------------------------------------------"
!      print*, "     "
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
      end subroutine split_fwd
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine pressure(nx,ny,nl,pres,u1,u2,u3,u5)
      integer :: nx,ny,nl
      integer :: i,j
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         pres
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl), intent(in) :: 
     &         u1,u2,u3,u5
      double precision :: ga1
!
!
      ga1 = 0.4d0
!
!
!
      do j = 1-nl,ny+nl
      do i = 1-nl,nx+nl
!
      if (u1(i,j).NE.0.0d0) then
      pres(i,j) = ga1*(u5(i,j) - 
     & 0.5d0*(u2(i,j)**2.0d0+u3(i,j)**2.0d0)/u1(i,j))
      end if
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
      end subroutine pressure
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine meshing2(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!    
      integer:: nx,ny,na,nl,np
!    
      integer :: n,i,j,k
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) ::
     &    meshX,meshY
!
!     try increasing all array vectors by one to account for nx inc.
      double precision, dimension(0:na+2) :: xc,yc,a,c,l,z
      double precision, dimension(0:na+1) :: h,mu,b,d
      double precision, dimension(1:na+1) :: be
      double precision, dimension(1:nx+2) :: xt,yt,phi,phi_ig,ps,xi,y0
!
      double precision :: Lx,Ly,s_pos,s_hgt,s_wdt
!
!
!
!
!
!
!     define numbers
!-----------------------------------------------------------------------
!
!     the design variables are the spline control points for the
!     y-coordinate of the mesh. they are initialised in the shape of a
!     squared cosine wave with a height of 1.0
!
!       *                                                               *
!   fixed       *                                               *    fixed
!          1 of na                                           na of na
!
!                       *                               *
!                  2 of na                             . . .
!
!                              *                *
!                       3 of na        *
!                                   . . .
!
!     the start and end points should NOT be provided, they are
!     accounted for by the meshing subroutine
!
!
      n = na + 1
      nx = nx + 1
      Lx = 1.0d0
      Ly = params(1,12)
      s_pos = params(1,9)
      s_hgt = params(1,10)
      s_wdt = params(1,11)
!      j = 0.0d0
!
!
!     initialise control points:
!
      xc(0) = 0.0d0
      xc(n) = 1.0d0
      yc(0) = 1.0d0
      yc(n) = 1.0d0
!
      do i = 1,na
      xc(i) = DBLE(i)/(na+1.0d0)
      yc(i) = alpha(1,i)
      end do
!
!
!
!
!
!
!
!     create spline polynomial
!-----------------------------------------------------------------------
!
!     first arrays
!
      do i = 0,n
      a(i) = yc(i)
      end do
!
      do i = 0,n-1
      h(i) = xc(i+1) - xc(i)
      end do
!
      do i = 1,n-1
      be(i) = 3*(a(i+1)-a(i))/h(i) - 3*(a(i)-a(i-1))/h(i-1)
      end do
!
!
!
!     next stage of arrays
!
      l(0) = 1.0d0
      mu(0) = 0.0d0
      z(0) = 0.0d0
!
      do i = 1,n-1
      l(i) = 2*(xc(i+1)-xc(i-1)) - h(i-1)*mu(i-1)
      mu(i) = h(i)/l(i)
      z(i) = (be(i) - h(i-1)*z(i-1))/l(i)
      end do
!
      l(n) = 1.0d0
      z(n) = 0.0d0
!
!
!     final arrays
!
      c(n) = 0.0d0
!
      do i = 1,n
      j = n-i
      c(j) = z(j) - mu(j)*c(j+1)
      b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2.0d0*c(j))/3.0d0
      d(j) = (c(j+1)-c(j))/(3.0d0*h(j))
      end do
!
!
!!      do i = 0,n-1
!!      print*, a(i), b(i), c(i), d(i)
!!      end do
!
!
!
!
!     form the actual spline using polynomial coefficients
      do i = 1,nx
!
!
!     the actual points that form the spline
      xt(i) = DBLE(i-1)*Lx/(nx-1.0d0)
!
!
!     determine what j value to use; check through control points
      do k = 0,n
      if (xc(k) .GT. xt(i)) then
        j = k-1
        goto 33
      else
        j = 0
      end if
      end do
33    continue
!
!
!     do the actual spline formula
      yt(i) = a(j)+
     &        b(j)*(xt(i)-xc(j))**1.0d0 +
     &        c(j)*(xt(i)-xc(j))**2.0d0 +
     &        d(j)*(xt(i)-xc(j))**3.0d0
!
!
!
!
!
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
!
!
!
!
!     adaptive meshing IF PARAMS SAYS SO
!-----------------------------------------------------------------------
      if (1 .EQ. 1) then
!
!     create potential function and integrate it
!
      phi_ig(1) = 0.0d0
!
      do i = 2,nx
      phi(i) = 1.0d0/yt(i) + s_hgt*3**(-s_wdt*(i-s_pos*nx)**2)
      phi_ig(i) = phi_ig(i-1) + phi(i)
      end do
!
!
!     scale the potential function
!
      do i = 1,nx
      phi_ig(i) = (phi_ig(i) - phi_ig(1))
     &           /(phi_ig(nx) - phi_ig(1))
      end do
!
!
!
!
!
!     seed the CDF and interpolate the values back to the x axis
!
      do i = 1,nx
!
!     calculate seed points to use
      ps(i) = (DBLE(i)-1)/(DBLE(nx)-1)
      y0(i) = 0.0d0
!
!     compare psi_ig values to find points for linear interpolation
      do k = 1,nx
!
!     if a seed point is equal to an integrated point, avoid /0
      if (phi_ig(k) .EQ. ps(i)) then
        xi(i) = xt(k)
        goto 44
!
!     else perform linear interpolation
      else if (phi_ig(k) .GT. ps(i)) then
        xi(i) = xt(k-1) + (ps(i)-phi_ig(k-1))*
     & (xt(k)-xt(k-1))/(phi_ig(k)-phi_ig(k-1))
        goto 44
      end if
!
!     exit checking loop
      end do
44    continue
!
!     continue seed check loop
      end do
!
!
!
!
!
!
!     write to relevant vector for consistency
!
      do i = 1,nx
      xt(i) = xi(i)
      end do
!
!
!
!
!
!
!     have to now re-fit y according to the spline
      do i = 1,nx
!
!     determine what j value to use; check through control points
      do k = 0,n
      if (xc(k) .GT. xt(i)) then
        j = k-1
        goto 55
      else
        j = 0
      end if
      end do
55    continue
!
!     do the actual spline formula
      yt(i) = a(j)+
     &        b(j)*(xt(i)-xc(j))**1.0d0 +
     &        c(j)*(xt(i)-xc(j))**2.0d0 +
     &        d(j)*(xt(i)-xc(j))**3.0d0
!
!
!
      end do
!
!
!
!!      do i = 1,nx
!!      print*, xt(i)
!!      end do 
!
!
!
!     end of the adaptive meshing section ---------
      end if
!
!
!      do i = 1,nx
!      print*, i, xt(i), yt(i)
!      end do
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
!
!
!     create mesh_ using linear spacing; scale to account for halo cells
!-----------------------------------------------------------------------
!
      nx = nx - 1
!
!
!     meshX: body x cells
!
      do i = 0,nx
      do j = 0-nl,ny+nl
      meshX(i,j) = xt(i+1)
      end do
      end do 
!
!
!     halo x cells
!
      do i = 1,nl
      do j = 0-nl,ny+nl
      meshX(0-i,j) = -DBLE(i)*Lx/(nx+1)
      end do
      end do 
!
      do i = 1,nl
      do j = 0-nl,ny+nl
      meshX(nx+i,j) = Lx + DBLE(i)*Lx/(nx+1)
      end do
      end do 
!
!
!
!
!
!     meshY:
!
      do i = 0-nl,nx+nl
      do j = 0-nl,ny+nl
      meshY(i,j) = Ly*(j*1.0d0/DBLE(ny) - 0.5d0)
      end do
      end do 
!
      do i = 0,nx
      do j = 0-nl,ny+nl
      meshY(i,j) = yt(i+1)*Ly*(j*1.0d0/DBLE(ny) - 0.5d0)
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
!
!
!     final closing bits
!-----------------------------------------------------------------------
!
!     write info to file
!
!      open(100, file='mesh.plt')
!      write(100,*) 'VARIABLES = "x" "y"'
!      write(100,*) 'Zone I=',ny+1,', J=',nx+1,', F=POINT'
!101   format(2f14.8)
!!
!      do i = 0,nx
!      do j = 0,ny
!            write(100,101) meshX(i,j), meshY(i,j)
!      end do
!      end do
!!
!      close(100)
!!
!!
!!
!      open(200, file='transect.plt')
!      write(200,*) 'VARIABLES = "x" "y"'
!201   format(2f14.8)
!!
!      do i = 1,nx+1
!            write(200,201) xt(i), yt(i)
!      end do
!!
!      close(200)
!
!
!
!
!
      end subroutine meshing
