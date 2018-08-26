C     THIS PROGRAM SHOW THE EXISTENCE OF INT NUMBERS DIFERENT OF ZERO
      PROGRAM DOBLEZERO
      INTEGER N
      N =2**30
C     UNCOMENT FOR
      N = -N
      WRITE(*,*) '-N= ' N
      N = N + N
      WRITE(*,*)'-2N' N
      N = N + N
      WRITE(*,*)'4N =' N
      READ(*,*)
      STOP
      END
