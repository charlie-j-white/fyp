      program main
      real x1, x2, x3, y1, y2
      x1 = 0.0
      x2 = 2.0
      x3 = 3.0
      call sub1(x1, x2, x3, y1, y2)
      print *, 'results are y1=', y1, ' y2=', y2
      end

      subroutine sub1(first, other, third, out, res)
      real first, other, third, out, res
      include 'globals.inc'
c start mixing all variables :
      out = first*other + 3.0*third
      res = sin(first) + cos(other)
      first = 2.0 * first
      third = 2.0/other
c now overwrite some values...
      if (first.gt.other) then
         first = cos(out)
         other = F(out)
      else
c ...in one or in all cases
         first = exp(out)
         other = G(out)
      endif
      end
      subroutine sub2(x,y,z,N,o)
      real a,x,y,z
      integer n,o,i
      dimension x(2*N),y(2*N),z(2*N)
C a normal loop:
      a = 0.5 * X(20)
      DO 200 i=5+o,N,2
C with an if-then and a goto
         if (Z(i) .gt. 0.0) then
            Z(i) = Z(i-1)-2*Z(i)+Y(i+1)
            X(i) = 3*X(i) - Y(i+1)*Y(i-1)
         else if (Z(i).lt.-10.0) then
            Z(i) = -Z(i)
         else
            Z(i) = cos(X(i))
            goto 100
         endif
         Y(i) = log(Z(i))
 100     Y(i) = Y(i) + 1
 200  continue
      a = X(10) + a
C a while loop:
      i = 1
      DO while(i.lt.N+5)
         i = i+5
         Z(i) = Z(i-1)+Z(i+1)
         X(i) = 3*X(i) - Y(i+1)*Y(i-1)
      enddo
      end
