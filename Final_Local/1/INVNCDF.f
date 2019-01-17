!23456
      PROGRAM Main
        REAL mu, sigma, p, answ, invncdf
        PRINT *,'Input F(x) (probability - p): '
        READ (*,*)p
        PRINT *,'Input mu: '
        READ(*,*)mu
        PRINT *,'Input sigma: '
        READ(*,*)sigma

        PRINT *,'Using Bisection to calculate the inverse of the NCDF wit
C     Continue string in next line with the & character
     &h the p value', p
        answ=invncdf(mu,sigma,p)
        PRINT *,'The value of x should be: ', answ
        READ(*,*)
      STOP 'Finishing program...'
      END
C     ******************************************************************
      REAL FUNCTION invncdf(mu,sigma,p)
        REAL mu,sigma,p,anotherFunc,a,b,accy,rtbis
        EXTERNAL anotherFunc
        PARAMETER(z1=-3, z2=3)
        a=(sigma*z1)+mu
        b=(sigma*z2)+mu
        accy=0.000001
        invncdf=rtbis(anotherFunc,a,b,accy,mu,sigma,p)
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION anotherFunc(x,mu,sigma,p)
         REAL x,mu,sigma,p,ncdf
         anotherFunc=ncdf(x,mu,sigma)- p
C         print*,'anotherFunc=',anotherFunc
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION ncdf(x,mu,sigma)
        REAL x, mu, sigma,err
        ncdf = (1+ err(x,mu,sigma))/2
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
C        print *,'Area= ', XI
C        print *,'PI= ', pi
        err=2/sqrt(pi)*XI
C        print *,'Error= ',err
       RETURN
      END
C     ******************************************************************
      REAL FUNCTION myFunc(x)
        REAL x
        myFunc = exp(-(x**2))
       RETURN
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
      FUNCTION rtbis(func,x1,x2,xacc,mu,sigma,p)
        INTEGER jmax
        REAL rtbis, x1, x2, xacc, func, mu, sigma, p
        EXTERNAL func
        PARAMETER(jmax=40)
        INTEGER j
        REAL dx,f,fmid,xmid
        fmid=func(x2,mu,sigma,p)
C        print *,'fmid=',fmid
        f=func(x1,mu,sigma,p)
C        print *,'f=',f
        IF(f*fmid.GE.0.) STOP 'root must be bracketed in rtbis'
        IF(f.LT.0.)THEN
          rtbis=x1
C          print *,'rtbis=',rtbis
          dx=x2-x1
C          print *,'dx=',dx
        ELSE
          rtbis=x2
C          print *,'rtbis=',rtbis
          dx=x1-x2
C          print *,'dx=',dx
        ENDIF
        DO j=1, jmax
C          print *,'j=',j
          dx=dx*.5
C          print *,'dx=',dx
          xmid=rtbis+dx
C          print *,'xmid=',xmid
          fmid=func(xmid,mu,sigma,p)
C          print *,'fmid=',fmid
        IF(fmid.LE.0.) THEN
          rtbis=xmid
C          print *,'rtbis=',rtbis
        ENDIF
        IF(abs(dx).LT.xacc .OR. fmid.EQ.0.)THEN
C          print *,'2IF -> rtbis=',rtbis
          RETURN
        ENDIF
        ENDDO
        STOP 'too many bisections in rtbis='
        READ(*,*)
      END
