!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine aresid(nx,ny,nl,flow,residual,np,params,dt,na,alpha)
!
!
      integer :: nx,ny,nl,np,na
      integer :: i,j
!
      double precision :: xa,xb,xc,xd,ya,yb,yc,yd,area
      double precision :: pi
!
      double precision, dimension(1,na) :: alpha
      double precision, dimension(1,np) :: params
      double precision, dimension(4) ::
     & fa,fb,fc,fd,ga,gb,gc,gd,
     & da,db,dc,dd
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow,residual
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5,
     &         r1,r2,r3,r5,
     &         volume,dt
!
!
!
!
!
!
!
!
!
!     Visualise the cell in three dimensions to see cross product
!
!                                mesh(i-1,j)         mesh(i,j)
!       z|                             .______________.
!        |   /y                       /       B      ^           
!        |  /                        /              /
!        | /                        /C            A/--------->n_A
!        |/______                  /              /
!               x                 /_______D______/
!                         mesh(i-1,j-1)          mesh(i,j-1)
!
!
!     Make sure to form the correct "face triangles" to get the right
!     normal
!      
!                             xA
!              .__________.<____.    
!     y|        \     B    \    ^  
!      |         \          \   |                        xA   0     yA
!      |          \C        A\  |yA     n.dS = A X k  =  yA X 0  = -xA 
!      |_____      \          \ |                         0   1      0
!           x       \____D_____\|
!
!
!     Area of a triangle between vectors P and Q for use in the cell
!     area calculation
!                                           area = 0.5*|PXQ|
!                                                = 0.5*|xP.yQ-xQ.yP|
!        
!
!
!        
!     When caulation the fluxes for each face use the definitions below
!     for example F_A = F(U(i+1/2,j))   
!        
!
!         *-------------*-------------*-------------*
!         |             |             |             |
!         |             |    i,j+1    |             |    A = (i+1/2,j)
!         |             |             |             |
!         *-------------*------B------*-------------*    B = (i,j+1/2)
!         |             |             |             |
!         |    i-1,j    C     i,j     A    i+1,j    |    C = (i-1/2,j)
!         |             |             |             |
!         *-------------*------D------*-------------*    D = (i,j-1/2)
!         |             |             |             |
!         |             |    i,j-1    |             |
!         |             |             |             |
!         *-------------*-------------*-------------*
!
!     A face: i+1/2 , j
!     B face: i     , j+1/2
!     C face: i-1/2 , j
!     D face: i     , j-1/2
!
!
!
!
!
!
!
!
!   
!     apply boundary conditions and calculate pressure field
!
      call meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
!
      call split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_fwd(nx,ny,nl,residual,r1,r2,r3,r5)
!
      call boundaries(nx,ny,nl,np,params,u1,u2,u3,u5,meshX,meshY)
!
!      call pressure(nx,ny,nl,pres,u1,u2,u3,u5)
!
      do i = 1,4
      da(i) = 0.0d0
      db(i) = 0.0d0
      dc(i) = 0.0d0
      dd(i) = 0.0d0
      end do
!
!
!
!     calculate cell area
!
      pi = 3.1415926535897932d0
!
      do j = 1-nl,ny+nl
      do i = 1-nl,nx+nl
!
!
      xa = meshX(i,j) - meshX(i,j-1)
      xb = meshX(i-1,j) - meshX(i,j)
      xc = meshX(i-1,j-1) - meshX(i-1,j)
      xd = meshX(i,j-1) - meshX(i-1,j-1)
      ya = meshY(i,j) - meshY(i,j-1)
      yb = meshY(i-1,j) - meshY(i,j)
      yc = meshY(i-1,j-1) - meshY(i-1,j)
      yd = meshY(i,j-1) - meshY(i-1,j-1)
!
!
!
      volume(i,j) = 0.5d0*DABS(xa*yb-xb*ya) + 0.5d0*DABS(xc*yd-xd*yc)
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
!     Perform finite volume calculation over body cells
!
      do j = 1,ny
      do i = 1,nx
!
!     calculate geometric information
!
      xa = (meshX(i,j) - meshX(i,j-1))
      xb = (meshX(i-1,j) - meshX(i,j))
      xc = (meshX(i-1,j-1) - meshX(i-1,j))
      xd = (meshX(i,j-1) - meshX(i-1,j-1))
      ya = (meshY(i,j) - meshY(i,j-1))
      yb = (meshY(i-1,j) - meshY(i,j))
      yc = (meshY(i-1,j-1) - meshY(i-1,j))
      yd = (meshY(i,j-1) - meshY(i-1,j-1))
!
      area = volume(i,j) 
!
!      print*
!      print*, xa,xb,xc,xd,ya,yb,yc,yd
!      print*, 
!     & 180.0d0*DATAN2(-xa,ya)/pi,
!     & 180.0d0*DATAN2(-xb,yb)/pi,
!     & 180.0d0*DATAN2(-xc,yc)/pi,
!     & 180.0d0*DATAN2(-xd,yd)/pi
!
!
!
!
!
!
!     create F and G vector
!
!     A face: i+1/2 , j
      call fg_vector(fa,ga,
     & 0.5d0*(u1(i+1,j)+u1(i,j)),
     & 0.5d0*(u2(i+1,j)+u2(i,j)),
     & 0.5d0*(u3(i+1,j)+u3(i,j)),
     & 0.5d0*(u5(i+1,j)+u5(i,j)))
!
!     B face: i     , j+1/2
      call fg_vector(fb,gb,
     & 0.5d0*(u1(i,j+1)+u1(i,j)),
     & 0.5d0*(u2(i,j+1)+u2(i,j)),
     & 0.5d0*(u3(i,j+1)+u3(i,j)),
     & 0.5d0*(u5(i,j+1)+u5(i,j)))
!
!     C face: i-1/2 , j
      call fg_vector(fc,gc,
     & 0.5d0*(u1(i-1,j)+u1(i,j)),
     & 0.5d0*(u2(i-1,j)+u2(i,j)),
     & 0.5d0*(u3(i-1,j)+u3(i,j)),
     & 0.5d0*(u5(i-1,j)+u5(i,j)))
!
!     D face: i     , j-1/2
      call fg_vector(fd,gd,
     & 0.5d0*(u1(i,j-1)+u1(i,j)),
     & 0.5d0*(u2(i,j-1)+u2(i,j)),
     & 0.5d0*(u3(i,j-1)+u3(i,j)),
     & 0.5d0*(u5(i,j-1)+u5(i,j)))
!
!
!
!
!
!     do JST stuff
!
      call jst_calcs(dt(i,j),np,params,
     & da,volume(i+1,j),volume(i,j),
     & u1(i+2,j),u1(i+1,j),u1(i,j),u1(i-1,j),
     & u2(i+2,j),u2(i+1,j),u2(i,j),u2(i-1,j),
     & u3(i+2,j),u3(i+1,j),u3(i,j),u3(i-1,j),
     & u5(i+2,j),u5(i+1,j),u5(i,j),u5(i-1,j)
     & )
!
      call jst_calcs(dt(i,j),np,params,
     & db,volume(i,j+1),volume(i,j),
     & u1(i,j+2),u1(i,j+1),u1(i,j),u1(i,j-1),
     & u2(i,j+2),u2(i,j+1),u2(i,j),u2(i,j-1),
     & u3(i,j+2),u3(i,j+1),u3(i,j),u3(i,j-1),
     & u5(i,j+2),u5(i,j+1),u5(i,j),u5(i,j-1)
     & )
!
      call jst_calcs(dt(i,j),np,params,
     & dc,volume(i,j),volume(i-1,j),
     & u1(i+1,j),u1(i,j),u1(i-1,j),u1(i-2,j),
     & u2(i+1,j),u2(i,j),u2(i-1,j),u2(i-2,j),
     & u3(i+1,j),u3(i,j),u3(i-1,j),u3(i-2,j),
     & u5(i+1,j),u5(i,j),u5(i-1,j),u5(i-2,j)
     & )
!
!
      call jst_calcs(dt(i,j),np,params,
     & dd,volume(i,j),volume(i,j-1),
     & u1(i,j+1),u1(i,j),u1(i,j-1),u1(i,j-2),
     & u2(i,j+1),u2(i,j),u2(i,j-1),u2(i,j-2),
     & u3(i,j+1),u3(i,j),u3(i,j-1),u3(i,j-2),
     & u5(i,j+1),u5(i,j),u5(i,j-1),u5(i,j-2)
     & )
!
!
!      print*, dt(i,j)
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
!     calculate residual for specific cell; subtract dissipation terms
!     from original FV section
!
      r1(i,j) = (-ga(1)*xa + fa(1)*ya
     &           -gb(1)*xb + fb(1)*yb
     &           -gc(1)*xc + fc(1)*yc
     &           -gd(1)*xd + fd(1)*yd
     & - (da(1) + db(1) - dc(1) - dd(1))
     & )/area 
!
      r2(i,j) = (-ga(2)*xa + fa(2)*ya
     &           -gb(2)*xb + fb(2)*yb
     &           -gc(2)*xc + fc(2)*yc
     &           -gd(2)*xd + fd(2)*yd
     & - (da(2) + db(2) - dc(2) - dd(2))
     & )/area 
!
      r3(i,j) = (-ga(3)*xa + fa(3)*ya
     &           -gb(3)*xb + fb(3)*yb
     &           -gc(3)*xc + fc(3)*yc
     &           -gd(3)*xd + fd(3)*yd
     & - (da(3) + db(3) - dc(3) - dd(3))
     & )/area 
!
      r5(i,j) = (-ga(4)*xa + fa(4)*ya
     &           -gb(4)*xb + fb(4)*yb
     &           -gc(4)*xc + fc(4)*yc
     &           -gd(4)*xd + fd(4)*yd
     & - (da(4) + db(4) - dc(4) - dd(4))
     & )/area 
!
!
!      print*, r1(i,j), (da(1) + db(1) - dc(1) - dd(1))
!      print*, r2(i,j), (da(2) + db(2) - dc(2) - dd(2))
!      print*, r3(i,j), (da(3) + db(3) - dc(3) - dd(3))
!      print*, r5(i,j), (da(4) + db(4) - dc(4) - dd(4))
!
!
!
!
!
!
!
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
      call split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
      call split_rev(nx,ny,nl,residual,r1,r2,r3,r5)
!
!
!
      end subroutine aresid
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine meshing(nx,ny,na,nl,alpha,meshX,meshY,np,params)
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
      double precision, dimension(0:na+1) :: xc,yc,a,c,l,z
      double precision, dimension(0:na) :: h,mu,b,d
      double precision, dimension(1:na) :: be
      double precision, dimension(1:nx+1) :: xt,yt,phi,phi_ig,ps,xi,y0
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
!      nx = nx + 1
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
      do i = 1,nx+1
!
!
!     the actual points that form the spline
!      xt(i) = DBLE(i-1)*Lx/(nx-1.0d0)
      xt(i) = DBLE(i-1)*Lx/(nx)
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
!     last value of yt seems buggy, add conditional statement
      if ( nx+1 .EQ. i ) then
        j = na
      end if
!
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
!     adaptive meshing IF PARAMS SAYS SO
!-----------------------------------------------------------------------
      if (1 .EQ. 1) then
!
!     create potential function and integrate it
!
      phi_ig(1) = 0.0d0
!
      do i = 2,nx+1
      phi(i) = 1.0d0/yt(i) + s_hgt*3**(-s_wdt*(i-s_pos*(nx+1))**2)
      phi_ig(i) = phi_ig(i-1) + phi(i)
      end do
!
!
!     scale the potential function
!
      do i = 1,nx+1
      phi_ig(i) = (phi_ig(i) - phi_ig(1))
     &           /(phi_ig(nx+1) - phi_ig(1))
      end do
!
!
!
!
!
!
!     seed the CDF and interpolate the values back to the x axis
!
      do i = 1,nx+1
!
!     calculate seed points to use
      ps(i) = (DBLE(i)-1)/DBLE(nx)
      y0(i) = 0.0d0
!
!     compare psi_ig values to find points for linear interpolation
      do k = 1,nx+1
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
      do i = 1,nx+1
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
!     sort out end of yt array
      if ( i .EQ. nx+1) then
        j = na
      end if
!
!
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
      end do
!
!
!
!
!
!
!
!     end of the adaptive meshing section ---------
      end if
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
!     create mesh_ using linear spacing; scale to account for halo cells
!-----------------------------------------------------------------------
!
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
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine split_fwd(nx,ny,nl,flow,u1,u2,u3,u5)
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
      subroutine boundaries(nx,ny,nl,np,params,u1,u2,u3,u5,meshX,meshY)
      integer :: nx,ny,nl,np
      integer :: i,j
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl) :: 
     &         u1,u2,u3,u5
      double precision, dimension(0-nl:nx+nl,0-nl:ny+nl) :: 
     &         meshX,meshY
      double precision, dimension(1,np) :: 
     &         params 
      double precision :: Minf,Prat,gam,ga1,
     &  Pinf,rinf,Pstag,rstag,
     &  Pext,
     &  k,Vw1,Vw2
      double precision, dimension(1:ny) :: 
     &  Pbc,rbc,Mp2,uint,Pext2
!
!
!
!
!
!
!
!     perform calculations
!-----------------------------------------------------------------------
!
!     take the flow parameters from the input file
!
      Minf = params(1,1)
      Prat = params(1,2)
      gam  = params(1,3)
      ga1  = gam - 1.0d0
!
!
!     calculate non-reflecting inflow and outflow conditions
!
      Pinf = 1.0d0/gam
      rinf = 1.0d0
      Pstag=Pinf*(1.0d0+0.5d0*ga1*Minf*Minf)**(gam/ga1)
      rstag=rinf*(1.0d0+0.5d0*ga1*Minf*Minf)**(1.0d0/ga1)
!
      do i = 1,ny
      uint(i) = u2(1,i)/u1(1,i)
      Mp2(i) = (uint(i)**2.0d0)/
     & (ga1*(1.0d0/ga1 + 0.5d0*Minf**2.0d0-0.5d0*uint(i)**2.0d0))
      Pbc(i) = Pstag/(1.0d0+0.5d0*ga1*Mp2(i))**(gam/ga1)
      rbc(i) = rstag/(1.0d0+0.5d0*ga1*Mp2(i))**(1.0d0/ga1)
      Pext2(i) = Pbc(i)*prat
      end do
!
      Pext = Pinf*prat
!
!
!
!
!
!     inflow conditions;
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,ny
      u1(1-j,i) = rbc(i)
      u2(1-j,i) = rbc(i)*uint(i)
      u3(1-j,i) = 0.0d0
!      u3(1-j,i) = u3(1,i)
      u5(1-j,i) = Pbc(i)/ga1 + 0.5d0*rbc(i)*uint(i)**2.0d0
      end do
      end do
!
!
!
!
!
!
!     outflow conditions; transient for all but u5
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,ny
      u1(nx+j,i) = u1(nx,i)
      u2(nx+j,i) = u2(nx,i)
      u3(nx+j,i) = 0.0d0
!      u3(nx+j,i) = u3(nx,i)
      u5(nx+j,i) = Pext/ga1 + 0.5d0*(u2(nx,i)**2.0d0)/u1(nx,i)
      end do
      end do
!
!
!
!
!
!
!
!     solid wall; use transient solutions
!-----------------------------------------------------------------------
!
      do j = 1,nl
      do i = 1,nx
!     bottom wall
      u1(i,1-j) = u1(i,1)
      u5(i,1-j) = u5(i,1)
!     top wall
      u1(i,ny+j) = u1(i,ny)
      u5(i,ny+j) = u5(i,ny)
      end do
      end do
!
!
!
!
!
!
!     solid wall; set normal velocity to zero
!-----------------------------------------------------------------------
!
!
      do j = 1,nl
      do i = 1,nx
!
!
!     bottom wall
!
      Vw1 = meshX(i,0) - meshX(i-1,0)
      Vw2 = meshY(i,0) - meshY(i-1,0)
      k = 2.0d0*(Vw1*u2(i,1)+Vw2*u3(i,1))/(Vw1*Vw1+Vw2*Vw2)
!
      u2(i,1-j) = k*Vw1 - u2(i,1)
      u3(i,1-j) = k*Vw2 - u3(i,1)
!
!
!
!     top wall
!
      Vw1 = meshX(i,ny) - meshX(i-1,ny)
      Vw2 = meshY(i,ny) - meshY(i-1,ny)
      k = 2.0d0*(Vw1*u2(i,ny)+Vw2*u3(i,ny))/(Vw1*Vw1+Vw2*Vw2)
!
      u2(i,ny+j) = k*Vw1 - u2(i,ny)
      u3(i,ny+j) = k*Vw2 - u3(i,ny)
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
!
!
!
!
!
      end subroutine boundaries
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine fg_vector(f,g,u1,u2,u3,u5)
      double precision :: u1,u2,u3,u5,p
      double precision, dimension(4) :: f,g
!
!
      p = 0.4d0*(u5-0.5d0*(u2*u2+u3*u3)/u1)
!
      f(1) = u2
      f(2) = u2*u2/u1 + p
      f(3) = u2*u3/u1
      f(4) = u2*(u5+p)/u1
!
      g(1) = u3
      g(2) = u2*u3/u1
      g(3) = u3*u3/u1 + p
      g(4) = u3*(u5+p)/u1
!
!
!
      end subroutine fg_vector
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine jst_calcs(dt,np,params,da,area1,area0,
     &                    u12,u11,u10,u1i,
     &                    u22,u21,u20,u2i,
     &                    u32,u31,u30,u3i,
     &                    u52,u51,u50,u5i
     & )
!
!
!
!
!
      integer :: np
      double precision, dimension(1,np) :: params
      double precision, dimension(4) :: da
      double precision :: u12,u11,u10,u1i,
     &                    u22,u21,u20,u2i,
     &                    u32,u31,u30,u3i,
     &                    u52,u51,u50,u5i,
     & area1,area0,dt,
     &                     p2, p1, p0, pi,
     & E2,E4,K2,K4,v1,v0,h
!
!
!
!
!
!     define constants and pressure sensors
!
      K2 = params(1,6)
      K4 = params(1,7)
!
      p2 = 0.4d0*(u52-0.5d0*(u22**2.0d0+u32**2.0d0)/u12)
      p1 = 0.4d0*(u51-0.5d0*(u21**2.0d0+u31**2.0d0)/u11)
      p0 = 0.4d0*(u50-0.5d0*(u20**2.0d0+u30**2.0d0)/u10)
      pi = 0.4d0*(u5i-0.5d0*(u2i**2.0d0+u3i**2.0d0)/u1i)
!
      v1 = DABS(p2-2.0d0*p1+p0)/(DABS(p2)+2.0D0*DABS(p1)+DABS(p0))
      v0 = DABS(p1-2.0d0*p0+pi)/(DABS(p1)+2.0D0*DABS(p0)+DABS(pi))
!
      E2 = K2*DMAX1(v1,v0)
      E4 = DMAX1(0.0d0,K4-E2)
!
      h = 0.5d0*(area1+area0)
!
!
!
!
!
!     now calculate the actual dissipation terms
!
      da(1) = (h/dt)*( E2*(u11-u10) - E4*(u12-3.0d0*u11+3.0d0*u10-u1i) )
      da(2) = (h/dt)*( E2*(u21-u20) - E4*(u22-3.0d0*u21+3.0d0*u20-u2i) )
      da(3) = (h/dt)*( E2*(u31-u30) - E4*(u32-3.0d0*u31+3.0d0*u30-u3i) )
      da(4) = (h/dt)*( E2*(u51-u50) - E4*(u52-3.0d0*u51+3.0d0*u50-u5i) )
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
!
!
!
!
!
!
!
!
      end subroutine jst_calcs
!**********************************************************************!
!                          Charlie Anderson                            !
!                  University of Bristol, March 2019                   !
!**********************************************************************!
!    
      subroutine split_rev(nx,ny,nl,flow,u1,u2,u3,u5)
      integer :: nx,ny,nl,nh,nb
      integer :: i,j,R
      double precision, dimension(1,4*(nx+2*nl)*(ny+2*nl)) :: 
     &         flow
      double precision, dimension(1-nl:nx+nl,1-nl:ny+nl), intent(in) :: 
     &         u1,u2,u3,u5
!      double precision :: 
!
!
!
!         1  2  3  4  5
!         6  7  8  9  10        ==>       1 2 3 4 5 6 7 ... 13 14 15    
!         11 12 13 14 15
!
!
      nh = 2*nl*(nx+ny)
      nb = nx*ny
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
!     assign first flow variable (wall_bottom, wall_top, inflow, outflow)
!
      do j = 1,nl
      do i = 1,nx
      R = 0
      flow(1,R+(j-1)*nx+i) = u1(i,1-j) 
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = nl*nx
      flow(1,R+(j-1)*nx+i) = u1(i,ny+j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx
      flow(1,R+(j-1)*ny+i) = u1(1-j,i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nl*nx+nl*ny
      flow(1,R+(j-1)*ny+i) = u1(nx+j,i)
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
      flow(1,R+(j-1)*nx+i) = u2(i,1-j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 1*nh+nl*nx
      flow(1,R+(j-1)*nx+i) = u2(i,ny+j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx
      flow(1,R+(j-1)*ny+i) = u2(1-j,i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 1*nh+2*nl*nx+nl*ny
      flow(1,R+(j-1)*ny+i) = u2(nx+j,i)
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
      flow(1,R+(j-1)*nx+i) = u3(i,1-j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 2*nh+nl*nx
      flow(1,R+(j-1)*nx+i) = u3(i,ny+j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx
      flow(1,R+(j-1)*ny+i) = u3(1-j,i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 2*nh+2*nl*nx+nl*ny
      flow(1,R+(j-1)*ny+i) = u3(nx+j,i)
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
      flow(1,R+(j-1)*nx+i) = u5(i,1-j)
      end do
      end do
!
      do j = 1,nl
      do i = 1,nx
      R = 3*nh+nl*nx
      flow(1,R+(j-1)*nx+i) = u5(i,ny+j) 
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx
      flow(1,R+(j-1)*ny+i) = u5(1-j,i)
      end do
      end do
!
      do j = 1,nl
      do i = 1,ny
      R = 3*nh+2*nl*nx+nl*ny
      flow(1,R+(j-1)*ny+i) = u5(nx+j,i)
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
      flow(1,R+(j-1)*nx+i) = u1(i,j)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+1*nb
      flow(1,R+(j-1)*nx+i) = u2(i,j)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+2*nb
      flow(1,R+(j-1)*nx+i) = u3(i,j)
      end do
      end do
!
      do j = 1,ny
      do i = 1,nx
      R = 4*nh+3*nb
      flow(1,R+(j-1)*nx+i) = u5(i,j)
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
      end subroutine split_rev
