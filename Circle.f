C     THIS PROGRAM READS A REAL NUMBER R AND PRINTS THE AREA OF A CIRCLE
C     WITH RADIUS r
      PROGRAM Circle
      REAL r, area
      WRITE(*,*)'Give radius r: '
      READ(*,*)r
      area=3.14159*r*r
      WRITE(*,*)'Area= ',area
      READ(*,*)
      STOP
      END
      

