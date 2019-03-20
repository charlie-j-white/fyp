      program main

      integer :: N, INFO
      integer, dimension (2,2) :: IPIV
      integer, dimension (3,3) :: IPIV2
      double precision, dimension (2,2) :: A
      double precision, dimension (1,2) :: f
      double precision, dimension (3,3) :: B
      double precision, dimension (1,3) :: g

!     A(col, row)
      A(1,1) = 1.0d0
      A(2,1) = 4.0d0
      A(1,2) = 4.0d0
      A(2,2) = 2.0d0
      f(1,1) = 10.0d0
      f(1,2) = 26.0d0

      print*, "       "
      print*, A(:,1)
      print*, A(:,2)
      print*, "       "
      print*, f 

      call DGETRF(2,2,A,2,IPIV,INFO)
      call DGETRS('N', 2, 1, A, 2, IPIV, f, 2, INFO)


      print*, "       "
      print*, f 

      print*, "       "
      print*, INFO 
      print*, "       "
      print*, "       "
      print*, "       "
      print*, "       "

!---------------------------------------------------------
      B(1,1) = 4.0d0
      B(2,1) = 2.0d0
      B(3,1) =-2.0d0
      B(1,2) = 2.0d0
      B(2,2) = 8.0d0
      B(3,2) = 4.0d0
      B(1,3) =30.0d0
      B(2,3) =12.0d0
      B(3,3) =-4.0d0
      g(1,1) =10.0d0
      g(1,2) =32.0d0
      g(1,3) =24.0d0

      print*, "       "
      print*, B(:,1)
      print*, B(:,2)
      print*, B(:,3)
      print*, "       "
      print*, g(:,1)
      print*, g(:,2)
      print*, g(:,3)

      call DGETRF(3,3,B,3,IPIV2,INFO)
      call DGETRS('T', 3, 1, B, 3, IPIV2, g, 3, INFO)


      print*, "       "
      print*, g(:,1)
      print*, g(:,2)
      print*, g(:,3)

      print*, "       "
      print*, INFO 
      print*, "       "
      print*, "       "








      end program main
