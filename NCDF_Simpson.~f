      program Main
       integer n
       real a, b, XI, myFunc, err, mu, sigma, x, ncdf
       parameter(pi=3.14159)
       external myFunc
       print *,'Input x: '
       read (*,*)x
       print *,'Input mu: '
       read(*,*)mu
       print *,'Input sigma: '
       read(*,*)sigma
       n=4
       a=0
       b=(x-mu)/sqrt(2.)*sigma
       call SimpsonComp(a,b,n,XI, myFunc)
       print *,'Area = ',XI
       err=(2/sqrt(pi))*XI
       print *,'The error function is= ',err
       ncdf=(1+err)/2
       print *,'The Normal Cumulative Distribution is: ',ncdf
       read(*,*)
       stop 'Finishing program...'
      end
C     ******************************************************************
      subroutine SimpsonComp(a,b,n,XI,func)
        integer n
        real a,b, XI, XI0, X, X1, X2, func
        real h
        integer i
        h=(b-a)/n
        XI1=0
        XI2=0
        XI0=func(a)+func(b)
        do i=1, n-1
          X=a+i*h
          if(mod(i,2).EQ.0) then
            XI2=XI2+func(X)
          else
            XI1=XI1+func(X)
          endif
        enddo
        XI=h/3*(XI0+2*XI2+4*XI1)
      return
      end
C     ******************************************************************
      real function myFunc(x)
        real x
        myFunc = exp(-(x**2))
        return
      end
