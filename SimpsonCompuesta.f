      program Main
       integer n
       real a, b, XI
       print *,'Input a: '
       read (*,*)a
       print *,'Input b: '
       read(*,*)b
       print *,'Iput n: '
       read(*,*)n
       call SimpsonComp(a,b,n,XI)
       print *,'Area = ',XI
       read(*,*)
      stop
      end
C     ******************************************************************
      subroutine SimpsonComp(a,b,n,XI)
        integer n
        real a,b, XI, XI0, X, X1, X2
        real h
        integer i
        h=(b-a)/n
        XI1=0
        XI2=0
C        XI0=exp(-(a**2))+exp(-(b**2))
        XI0=exp(a)+exp(b)
        do i=1, n-1
          X=a+i*h
          if(mod(i,2).EQ.0) then
            XI2=XI2+exp(X)
          else
            XI1=XI1+exp(X)
          endif
        enddo
        XI=h/3*(XI0+2*XI2+4*XI1)
      return
      end
C     ******************************************************************
