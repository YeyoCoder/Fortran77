      PROGRAM LineContinuation
C23456789 (THIS DEMOSTRATES COLUMN POSITION!)
C     THE NEXT STATEMENT GOES OVER TWO PHYSICAL LINES
      REAL r,area
      PARAMETER(pi=3.14159265358979)
      WRITE(*,*)'Input radius r: '
      READ(*,*)r
      area =pi
     +     *r *r
      WRITE(*,*)'Area = ',area
      READ(*,*)
      STOP
      END
      
