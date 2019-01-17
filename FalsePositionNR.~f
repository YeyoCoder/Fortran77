C     FALSE POSITION
      FUNCTION FalsePosition(func,x1,x2,xacc)
      INTEGER maximumIteration
      REAL root, x1,x2,xacc,func
      EXTERNAL func
C     SET TO THE MAXIMUM ALLOWED NUMBER OF ITERATIONS
      PARAMETER (maximumIteration=30)
C     USING THE FALSE POSITION METHOD, FIND THE ROOT OF A FUNCTION(func)
C     KNOWN TO LIE BETWEEN x1 AND x2. THE ROOT, RETURNED AS root, IS
C     REFINED UNTIL ITS ACCURACY IS +/- xacc
      INTEGER j
      REAL del,dx,f,fh,fl,swap,xh,xl
C     BE SURE THE INTERVAL BRACKET A ROOT
      f1=func(x1)
      fh=func(x2)

      IF(fl*fh.GT.0) PAUSE 'root must be bracketed in FalsePosition'
      IF(fl.LT.0)THEN
        xl=x2
        xh=x1
      ELSE
        xl=x2
        xh=x1
        swap=fl
        fl=fh
        fh=swap
      ENDIF
      dx=xh-hl
C     FALSE POSITION LOOP
      DO j=1,maximumIteration
        root=xl+dx*fl/(fl-fh)
        f=func(root)
C     REPLACE APPROPIATE LIMIT
        IF(f.LT.0)
          del=xl-root
          xl=root
          fl=f
        ELSE
          del=xh-root
          xh=root
          fh=f
        ENDIF
        dx=xh-xl
        IF((del).LT.xacc.OR.f.EQ.0) RETURN
      ENDDO
      PAUSE 'FalsePosition exceed maximum iterations'
      END
      
      
      
      
