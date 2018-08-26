      PROGRAM INTBISTRING
      INTEGER N
      WRITE(*,*)'Input a Number: '
      READ(*,*) N
      IF(N.LT.0) THEN
        WRITE(*,*)'1'
      ELSE
        WRITE(*,*)'0'
      ENDIF
      STOP
      END
      
