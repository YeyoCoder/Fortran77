C     BISECTION METHOD TO FIND
C     THE ROOT OF F(X) = X**2 -2

C     MêTODO DE LA BISECCI‡N PARA HALLAR
C     LAS RA÷CES DE  F(X) = X**2 -2
      PROGRAM BISECTION
C     VARIABLE DECLARATIONS
      REAL A, B, C, D, FA ,FB, FC
C     INSTRUCTIONS
      A = 1
      B = 2
      C = (A+B)/2
      FA = A**2-2
      FB = B**2-2
      FC = C**2-2
      D = FA*FB
      
      WRITE(*,*)'         i    |a           |b            |c            |d'
C     .LT. = LEST THAN
      IF (D .LT. 0) THEN
        DO I=1, 10
          WRITE(*,*) I, A, B, C, D
          IF (D.GT.0) THEN
            A = C
          ELSE
            B = C
          END IF
          C = (A+B)/2
          FA = A**2-2
          FC = C**2-2
          D = FA*FC
        END DO
      ELSE
          WRITE(*,*) 'No root in the interval'
      END IF
      READ(*,*)
      STOP
      END
