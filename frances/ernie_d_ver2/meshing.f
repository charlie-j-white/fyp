



     
      subroutine meshing(area_rat, N, S, Sx, dx, xpos)

      implicit none

      integer, intent(in) :: N
      double precision, intent(in) :: area_rat
      double precision, intent(out) :: dx
      double precision, dimension(0:N+1), intent(out) :: S, Sx
      double precision, dimension(-1:N+1), intent(out) :: xpos

      integer :: ii
      double precision :: L, x, r, pi







      L = 1.0
      r = area_rat
      pi = 3.1415
      dx = L/N


      do ii = 1,N

            xpos(ii) = ii*L/N
            x = (ii-0.5)*L/N 

            S(ii) = ((r+1)/(2*r))*(1+(r-1)/(r+1)*cos(2*pi*x/L))
            Sx(ii) = pi*((1-r)/(r*L))*sin(2*pi*x/L)

      end do

      S(0) = S(1) 
      S(N+1) = S(N)

      Sx(0) = Sx(1) 
      Sx(N+1) = Sx(N)


      xpos(-1) = -L/N
      xpos(0) = 0
      xpos(N+1) = 1+L/N











      end subroutine meshing
