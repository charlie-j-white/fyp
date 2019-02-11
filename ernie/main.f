
      program main
      
      implicit none

      integer :: N, ii, jj, tmax
      double precision :: mach_inf, area_rat, pres_rat, dx, dt, cfl
      double precision, dimension(:), allocatable :: S, Sx, xpos, 
     &          pres, u1, u2, u5 



! define primary flow variables
!-----------------------------------------------------------------------
      
      mach_inf = 0.3
      area_rat = 2.0
      pres_rat = 0.8
      cfl = 0.9
      N = 1000
      tmax = 100000



! allocate sizes for arrays 
!-----------------------------------------------------------------------

      allocate(S(0:N+1))
      allocate(Sx(0:N+1))
      allocate(xpos(-1:N+1))
      allocate(pres(0:N+1))
      allocate(u1(0:N+1))
      allocate(u2(0:N+1))
      allocate(u5(0:N+1))



! create mesh 
!-----------------------------------------------------------------------
!
!             -1   0   1   2  ... N-2 N-1  N  N+1
!
! position         0                       L
! points       x   o   o   o  ...  o   o   o   x 
! cells          H   C   C    ...    C   C   H
!
!                0   1   2    ...   N-1  N  N+1
!

      call meshing(area_rat, N, S, Sx, dx, xpos)


! initialise variables 
!-----------------------------------------------------------------------

      call initialise(N, mach_inf, pres_rat, pres, u1, u2, u5)


! start time stepping
!-----------------------------------------------------------------------
      do ii = 1,tmax
             call boundaries(N, mach_inf, pres_rat, pres, u1, u2, u5)
             call update(N, mach_inf, pres_rat, pres, u1, u2, u5, 
     &                          S, Sx, dt, dx, cfl)
      end do




! correct variables for area
!-----------------------------------------------------------------------
!      do ii = 0,N+1
!             u1(ii) = u1(ii)/S(ii)
!             u2(ii) = u2(ii)/S(ii)
!             u5(ii) = u5(ii)/S(ii)
!      end do





! write info to file
!-----------------------------------------------------------------------
      open(100, file='solution.plt')
  
      write(100,*) 'VARIABLES = "x" "u1" "u2" "u5" "pres"'

101   format(5f12.8)

      do ii = 1,N

            write(100,101) xpos(ii),
     &                  0.5*(u1(ii)+u1(ii+1)),
     &                  0.5*(u2(ii)+u2(ii+1)),
     &                  0.5*(u5(ii)+u5(ii+1)),
     &                  0.5*(pres(ii)+pres(ii+1))
                                

      end do
        



! 
!-----------------------------------------------------------------------

      endprogram main
