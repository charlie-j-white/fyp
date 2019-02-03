
      program main
      
      implicit none

      integer :: Nx, Ny, ii, jj
      double precision :: x, y, f, fx, fy



      Nx = 5
      Ny = 5




! write info to file
!-----------------------------------------------------------------------

      open(100, file='solution.plt')
      write(100,*) 'VARIABLES = "x" "y" "F" "Fx" "Fy"'
101   format(5f12.8)





! start iteration
!-----------------------------------------------------------------------
      do ii = 1,Nx
      do jj = 1,Ny

             x = ii/1.0
             y = jj/1.0
             call update(x, y, f, fx, fy)
            
             write(100,101) x, y,
     &                  f, fx, fy

      end do
      end do



        

      endprogram main
