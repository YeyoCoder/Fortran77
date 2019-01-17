C23456
      PROGRAM Main
       REAL mu, sigma, x, ncdf
       PRINT *,'Input x: '
       READ (*,*)x
       PRINT *,'Input mu: '
       READ(*,*)mu
       PRINT *,'Input sigma: '
       READ(*,*)sigma
       answ=0
       answ = ncdf(x, mu, sigma)
       PRINT *,'The Normal Cumulative Distribution is: ',answ
       READ(*,*)
       STOP 'Finishing program...'
      END
C     ******************************************************************
      SUBROUTINE SimpsonComp(a,b,n,XI,func)
        INTEGER n
        REAL a,b, XI, XI0, X, X1, X2, func
        REAL h
        INTEGER i
        h=(b-a)/n
        XI1=0
        XI2=0
        XI0=func(a)+func(b)
        DO i=1, n-1
          X=a+i*h
          IF(mod(i,2).EQ.0) THEN
            XI2=XI2+func(X)
          ELSE
            XI1=XI1+func(X)
          ENDIF
        ENDDO
        XI=h/3*(XI0+2*XI2+4*XI1)
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION myFunc(x)
        REAL x
        myFunc = exp(-(x**2))
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION err(x,mu,sigma)
        INTEGER n
        REAL x,mu,sigma,a,b,XI,myFunc
        EXTERNAL myFunc
        PARAMETER(pi=3.14159)
        n=50
        a=0
        b=(x-mu)/(sqrt(2.)*sigma)
        XI=0
        CALL SimpsonComp(a,b,n,XI,myFunc)
        err=2/sqrt(pi)*XI
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION ncdf(x,mu,sigma)
        REAL x, mu, sigma
        ncdf = (1+ err(x,mu,sigma))/2
       RETURN
      END
