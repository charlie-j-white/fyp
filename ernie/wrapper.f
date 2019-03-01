      program wrapper

      integer :: N
      double precision :: area_rat, J, Jp, dhi, Jd

      N = 1000
      area_rat = 2.0d0

      ! COST FUNCTION
      call main(area_rat, J, N)
      print*, J



      ! FINITE DIFFERENCE
      dh = 0.000000001
      area_rat = 2.0d0 + dh
      call main(area_rat, Jp, N)
      print*, (Jp - J)/dh


      ! AUTOMATIC DIFFERENTIATION
      area_rat = 2.0d0
      call main_d(area_rat, 1.0d0, J, Jd, N)
      print*, J
      print*, Jd


      end program wrapper
