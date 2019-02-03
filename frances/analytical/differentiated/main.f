
      program main
      
      implicit none

      integer :: Nx, Ny, ii, jj
      double precision :: x, y, f, fx, fy, fxd, fyd



      Nx = 5
      Ny = 5




! write info to file
!-----------------------------------------------------------------------

      open(100, file='solution.plt')
      write(100,*) 'VARIABLES = "x" "y" "F" "Fx" "Fy" "Fx_d" "Fy_d"'
101   format(7f12.8)





! start iteration
!-----------------------------------------------------------------------
      do ii = 1,Nx
      do jj = 1,Ny

             x = ii/1.0
             y = jj/1.0
!
!      first one: get wrt x (dxdt=1, dydt=0)
! second one one: get wrt y (dxdt=0, dydt=1)
!                 analytical cases handled automatically
!
             call update_d(x, 1.0, y, 0.0, f, fxd, fx, fy)
             call update_d(x, 0.0, y, 1.0, f, fyd, fx, fy)
            
             write(100,101) x, y,
     &                  f, fx, fy,
     &                  fxd, fyd

      end do
      end do



        

      endprogram main
