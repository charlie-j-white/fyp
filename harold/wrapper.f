      program wrapper

      double precision :: J,j_fdp,j_fdm,arearat_fd,delta

      integer :: nstages,maxit,ni
      double precision :: pmachin,arearat,prat,CFLscale,pk,restol


      pmachin = 0.3
      arearat = 2.0
      prat = 0.8
      CFLscale = 1.0
      pk = 0.3333333333333
      nstages = 5
      maxit = 10000
      ni = 1001
      restol = 0.00000001

! finite differences
      delta = 0.00000000001
      arearat_fd = arearat + delta
      call nozzle_main(pmachin,
     & arearat_fd,prat,CFLscale,pk,nstages,maxit,ni,restol,j_fdp)
      arearat_fd = arearat - delta
      call nozzle_main(pmachin,
     & arearat_fd,prat,CFLscale,pk,nstages,maxit,ni,restol,j_fdm)
      

! adjoint code
      call nozzle_main(
     &  pmachin,arearat,prat,CFLscale,pk,nstages,maxit,ni,restol,J)


! results      
      print*, "   "
      print*, "   "
      print*, "-------Finite difference--------"
      print*, "   J = ", j
      print*, "J_fd = ", j_fdp
      print*, "J_fd = ", j_fdm
      print*, "DJDa = ", (j_fdp-j)/delta
      print*, "DJDa = ", 0.5d0*(j_fdp-j_fdm)/delta
      print*, "   "
      print*, "   "
      print*, "---Automatic Differentiation----"
      print*, "DJDa = "
      print*, "   "
      print*, "   "
      print*, "   "
      print*, "   "
      print*, "---Comparison of DJDa values----"
      print*, "     FD = "
      print*, "     AD = "
      print*, "Adjoint = "
      print*, "   "
      print*, "   "




      end program wrapper
