C23456
      program Main
       real mu, sigma, x, fda, invfda
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
       answ = fda(x, mu, sigma)
       print *,'The Normal Cumulative Distribution is: ',answ
       print *,'Using Bisection to calculate the inverse of the NCDF wit
C     Continue string in next line with the & character
     &h the p value', answ
       answ=invfda(mu,sigma,answ)
       print *,'The value of x should be: ', answ
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
C     ******************************************************************
C     NCDF function. For a strange reason if I call this function ncdf
C     the program make bad arithmetics and return an misterious number
      real function fda(x,mu,sigma)
        real x, mu, sigma
        fda = (1+ err(x,mu,sigma))/2
        return
      end
C     ******************************************************************
      real function anotherFunc(x,mu,sigma,p)
         real x,mu,sigma,p
         anotherFunc=fda(x,mu,sigma)- p
C         print*,'anotherFunc=',anotherFunc
      return
      end
C     ******************************************************************
      real function invfda(mu,sigma,p)
        real mu,sigma,p,anotherFunc,a,b,accy,rtbis
        external anotherFunc,rtbis
        parameter(z1=-3, z2=3)
        a=(sigma*z1)+mu
C        print *,'a= ',a
        b=(sigma*z2)+mu
C        print *,'b= ',b
        accy=0.000001
C        print *,'accy=',accy
        invfda=rtbis(anotherFunc,a,b,accy,mu,sigma,p)
C       print *,'inversafda=',invfuncda
       return
      end
      
C     ******************************************************************
C     ******************************************************************
      function rtbis(func,x1,x2,xacc,mu,sigma,p)
        integer jmax
        real rtbis, x1, x2, xacc, func, mu, sigma, p
        external func
        parameter(jmax=40)
        integer j
        real dx,f,fmid,xmid
        fmid=func(x2,mu,sigma,p)
C        print *,'fmid=',fmid
        f=func(x1,mu,sigma,p)
C        print *,'f=',f
        if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'
        if(f.lt.0.)then
          rtbis=x1
C          print *,'rtbis=',rtbis
          dx=x2-x1
C          print *,'dx=',dx
        else
          rtbis=x2
C          print *,'rtbis=',rtbis
          dx=x1-x2
C          print *,'dx=',dx
        endif
        do j=1, jmax
C          print *,'j=',j
          dx=dx*.5
C          print *,'dx=',dx
          xmid=rtbis+dx
C          print *,'xmid=',xmid
          fmid=func(xmid,mu,sigma,p)
C          print *,'fmid=',fmid
        if(fmid.le.0.) then
          rtbis=xmid
C          print *,'rtbis=',rtbis
        endif
        if(abs(dx).lt.xacc .or. fmid.eq.0.)then
C          print *,'2IF -> rtbis=',rtbis
          return
        endif
        enddo
        stop 'too many bisections in rtbis='
        read(*,*)
      end
      
