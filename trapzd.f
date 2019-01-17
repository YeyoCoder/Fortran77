      program Main
      real area, a, b, s, myFunc
      external myFunc
      integer n
      a=1.0
      b=2.0
      s=0.0
      n=1
C      myFunc=func(b)
      call trapzd(myFunc,a,b,s,n)
      area = s
      print *,'area: ',area
      read (*,*)
      stop
      end
C     ******************************************************************
      real function myFunc(t)
        real t
        myFunc = exp(-(t**2))
C         myFunc = t
      return
      end
C     ******************************************************************
      subroutine trapzd (func, a, b, s, n)
        integer n
        real a, b, s, func
        external func
        integer it, j
        real del, sum, tnm, x
        print *,'Hello from trapzd!= '
        if (n.eq.1)then
          print *,'Hello from N=1'
          s=0.5*(b-a)*(func(a)+func(b))
          print *,'s= ',s
        else
          print *,'Hello from ELSE'
          it=2**(n-2)
          print *,'it=',it
          tnm=it
          print *,'tnm',tnm
          del=(b-a)/tnm
          print *,'del',del
          x=a+0.5*del
          print *,'x',x
          sum=0.
          do j=1, it
            sum=sum+func(x)
            print *,'sum',sum
            x=x+del
          enddo
          s=0.5*(s+(b-a)*sum/tnm)
          print*,'s=',s
        end if
        print *,'trapezoid integration subroutine ...'
        print *,'s= ',s
        read(*,*)
        return
      end
