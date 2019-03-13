!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     
      subroutine update(x1, x2, g, gx, gy)

      implicit none

      double precision :: x1, x2, g, gx, gy


      g = x1*sin(x1*x2) + x2
      gx = sin(x1*x2) + x1*x2*cos(x1*x2)
      gy = x1*x1*cos(x1*x2) + 1
      

      end subroutine update
