      program Main
C     This program calculate the Cumulative Distribution Function for a
C     Normal Distribution. Input any real value, the Population Mean and
C     the Population Standard Deviation
      integer n
      real x, mu, sigma, answ, a, upperLimit, s, myFunc
      external myFunc

      parameter (pi=3.14159)
      print *,'Input any real value for x: '
      read(*,*)x
      print *,'Input the Population Mean (mu): '
      read(*,*)mu
      print *,'Input the Population Standard Deviation (sigma):'
      read(*,*)sigma
      answ = 0
      a=0
      upperLimit=(x-mu)/(sqrt(2.0)*sigma)
      n=1
      call trapzd(myFunc, a, upperLimit, s, n)
      
C      call NCDF (x, mu, sigma, answ, func)
      print *,'The NCDF is: ',answ
      read(*,*)
      stop 'Finishing program...'
      end
C     ******************************************************************
      real function func(t)
        real t
        print *,'exponential function func...',t
        func = exp(-(t**2))
      return
      end function
C     ******************************************************************
      real function erf(area)
        real area
        integer j, m
        print *,'PI= ',pi
        erf=2/sqrt(pi)*area
      return
      end
C     ******************************************************************
      subroutine trapzd (func, a, b, s, n)
        integer n
        real a, b, s, func
        external func
        integer it, j
        real del, sum, tnm, x
        if (n.eq.1)then
          s=0.5*(b-a)*(func(a)+func(b))
          print *,'func(a)= ', func(a)
        else
          it=2**(n-2)
          tnm=it
          del=(b-a)/tnm
          x=a+0.5*del
          do j=1, it
            sum=sum+func(x)
            x=x+del
          enddo
          s=0.5*(s+(b-a)*sum/tnm)
        end if
        print *,'trapezoid integration subroutine ...'
        print *,'s= ',s
        return
      end
C     ******************************************************************
      subroutine NCDF (x, mu, sigma, answ, func)
C     Input/Output parameters
      real x, mu, sigma, answ, func
C     Local variables
      external func
      call trapzd(func, a, upperLimit, s, n)
      print *,'The upper limit is: ',upperLimit
      print *,'The area is: ', s
      print *,'The Error is: ', erf(s)
      answ = (1+erf(s))/2
      print *,'The answer is: ',answ
      return
      end
C     ******************************************************************
