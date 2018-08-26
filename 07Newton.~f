C     MêTODO DE NEWTON PARA F(X)=X**2-COS(X)
C     CON APROXIMACION INICIAL X=1
      PROGRAM NEWTONRAPHSON
      
      REAL F, DF, X
      X = 1
      DO I = 1, 10
         WRITE (*,*) X
         X = X - F(X)/DF(X)
      END DO
      READ(*,*)
      STOP
      END
      
      
      REAL FUNCTION F(X)
      F = X**2-COS(X)
      RETURN
      END
      
      REAL FUNCTION DF(X)
      DF = 2X+SIN(X)
      RETURN
      END
