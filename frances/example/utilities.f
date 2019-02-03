      double precision function f(t)
      double precision t
      f = exp(t*t)
      return
      end

      double precision function g(t)
      double precision t
      g = 1
      if (t .ne. 0.d0) then
        g = sin(t) / t
      endif
      return
      end

      subroutine sub3(z,t)
      real z,t,u,v
      integer i
      include 'globals.inc'
c c1
      i = 5
      x(1) = y*z + t
      call sub4(u,x(i),z,v,t)
      t = t + x(1)*z + 3*v
c c2
      end
c c3
c d0
      subroutine sub4(u,y2,z,v,t)
      real u,y2,z
      dimension y2(0:6)
      include 'globals.inc'
c d1
      u = u*y + y2(3)*z
      y = z+v*y
      v = u*y2(5)
c d2
      end

