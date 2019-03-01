



     
      subroutine qvector(u1, u2, u5, q1, q2, q5, S, Sx)

      implicit none

      double precision, intent(in) :: u1, u2, u5, S, Sx
      double precision, intent(out) :: q1, q2, q5
      
      double precision :: pres, gam

      gam = 1.40


      

! clculate pressure
      pres = (gam-1)*(u5 - 0.5*(u2**2)/u1)/S


! calculate q vector
      q1 = 0.0d0
      q2 = pres*Sx
      q5 = 0.0d0


      end subroutine qvector
