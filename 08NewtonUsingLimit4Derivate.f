C     MêTODO DE NEWTON PARA F(X)=X**2-COS(X)
C     CON APROXIMACION INICIAL X=1
      PROGRAM NEWTONRAPHSON

      REAL F, DF, X
      X = 1
      DO I = 1, 10
         WRITE (*,*) I, X
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
C      REAL H =0.1
C     PARECE QUE NO SE PUEEN DEFINIR LAS VARIABLES DENTRO DE LAS FUNCIONES
      DF = (F(X+.01)-F(X))/.01
      RETURN
      END
