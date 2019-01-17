C23456
      program Main
       real mu, sigma, x, answ
       parameter(pi=3.14159)
C       external myFunc
       print *,'Input x: '
       read (*,*)x
       print *,'Input mu: '
       read(*,*)mu
       print *,'Input sigma: '
       read(*,*)sigma
C       print *,'The error function is= ',err(x,mu,sigma)
       answ=0
C       print *,'NCDF already on memory?',ncdf
       answ = what(x, mu, sigma)
       print *,'The Normal Cumulative Distribution is: ',answ

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
C     ******************************************************************
      real function err(x,mu,sigma)
        integer n
        real x,mu,sigma,a,b,XI,myFunc
        external myFunc
        parameter(pi=3.14159)
        n=50
        a=0
        b=(x-mu)/(sqrt(2.)*sigma)
        XI=0
        call SimpsonComp(a,b,n,XI,myFunc)
C        print *,'Area= ', XI
C        print *,'PI= ', pi
        err=2/sqrt(pi)*XI
C        print *,'Error= ',err
      return
      end
C     ******************************************************************
      real function myNCDF(x, mu, sigma)
        real x, mu, sigma, error, tmp
C        print *,'1� ncdf',ncdf
        error = err(x,mu,sigma)
C        print*,'Error in ncdf func: ',err(x,mu,sigma)
        myNCDF=0.
C        print *,'2� ncdf',ncdf
        tmp= (1+ error) /2
        print *,'Check operations ->', tmp
        myNCDF = tmp
C        print *,'3� ncdf',ncdf
C         myNCDF = (1+ err(x,mu,sigma))/2
      return
      end
C     ******************************************************************
      real function what(x,mu,sigma)
        real x, mu, sigma
        what = (1+ err(x,mu,sigma))/2
        return
      end
C     ******************************************************************
C      real function g(x,mu,sigma,p)
C         real x,mu,sigma,p
C         g=ncdf(x,mu,sigma)- p
C      return
C      end
