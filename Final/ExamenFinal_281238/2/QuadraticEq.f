      PROGRAM Main
      REAL a,b,x
      PARAMETER (tol=0.01)
      COMPLEX x1,x2
      x1=(0,0)
      x2=(0,0)
      PRINT *,'Input coeficient a: '
      READ(*,*)a
      PRINT *,'Input coeficient b: '
      READ(*,*)b
      PRINT *,'Input coeficient c: '
      READ(*,*)c
      CALL roots(a,b,c,tol,x1,x2)
      PRINT*,'Root x1 is: ',x1
      PRINT*,'Root x2 is: ',x2
      READ(*,*)
      STOP 'Finishing program'
      END
C     ******************************************************************
      SUBROUTINE roots(a,b,c,tol,x1,x2)
      REAL a,b,c,tol,diff,disct,a1
      COMPLEX x1,x2
      diff=abs(abs(b)-sqrt(disct(a,b,c)))
C      IF(b.EQ.0.AND.disct(a,b,c).LT.0) THEN
      IF(disct(a,b,c).LT.0) THEN
        print*,'Imaginary roots. This program work only for real roots!'
        a1= b/(2*a)
        b1= sqrt(abs(disct(a,b,c)))/2*a
        x1=cmplx(a1,b1)
        x2=cmplx(-1*a1,-1*b1)
      ELSEIF(b.GT.0.AND.diff.LT.tol) THEN
        x1=(-b-(sqrt(disct(a,b,c))))/(2*a)
        x2=c/(a*x1)
        PRINT *,'using the alternative formula to calculate x2...'
      ELSEIF(b.LT.0.AND.diff.LT.tol)THEN
        x1=(-b+(sqrt(disct(a,b,c))))/(2*a)
        x2=c/(a*x1)
        PRINT *,'using the alternative formula to calculate x2...'
      ELSE
        PRINT *,'The difference is less than the tolerance ',tol
        x1=(-b+(sqrt(disct(a,b,c))))/(2*a)
        x2=(-b-(sqrt(disct(a,b,c))))/(2*a)
      ENDIF
      RETURN
      END
C     ******************************************************************
      REAL FUNCTION disct(a,b,c)
        REAL a,b,c
        disct=(b**2)-(4*a*c)
      RETURN
      END
      
