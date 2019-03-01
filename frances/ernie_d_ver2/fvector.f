



     
      subroutine fvector(u1, u2, u5, f1, f2, f5)

      implicit none

      double precision, intent(in) :: u1, u2, u5
      double precision, intent(out) :: f1, f2, f5
      
      double precision :: pres, vel, gam

      gam = 1.40
      ! area terms S inlculded in the u variable, or cancelled


      

! clculate pressure
      pres = (gam-1)*(u5 - 0.5*(u2**2)/u1)

! calculate velocity
      vel = u2/u1


! calculate f vector
      f1 = vel*u1
      f2 = vel*u2 + pres
      f5 = vel*u5 + pres*vel


      end subroutine fvector
