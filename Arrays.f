      PROGRAM Arrays
      REAL a(20), b(0:19), weird(-162:237)
      INTEGER i, sq(10)
      DO 100 i=1, 10
        sq(i)=i**2
        WRITE(*,*)sq(i)
  100 CONTINUE
      READ(*,*)
      STOP
      END
  
