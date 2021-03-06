!23456
      PROGRAM Main
        INTEGER n,np,i,j
        REAL a(10,10),b(10)
        n=4
        np=10
        a(1,1)= 1
        a(1,2)= 1
        a(1,3)= 0
        a(1,4)= 3
        a(2,1)= 2
        a(2,2)= 1
        a(2,3)= -1
        a(2,4)= 1
        a(3,1)= 3
        a(3,2)= -1
        a(3,3)= -1
        a(3,4)= 2
        a(4,1)= -1
        a(4,2)= 2
        a(4,3)= 3
        a(4,4)= -1
        PRINT *,'+------------   Crout Algorithm    ------------------+'
        PRINT *,'This program uses the Crout LU factorization to solve t
     &he system'
        PRINT *,'Ax=b, with'
        PRINT*,'A = ***************************************************'
        DO i=1,n
          WRITE(*,*) (a(i,j),j=1,n)
        ENDDO
        b(1)=8
        b(2)=7
        b(3)=14
        b(4)=-7
        PRINT *,'b= ***************************************************'
        DO i=1,n
          WRITE(*,*)b(i)
        ENDDO
        CALL ludcmp(a,n,np)
        PRINT*,'LU= ***************************************************'
        DO i=1,n
          WRITE(*,*) (a(i,j),j=1,n)
        ENDDO
        PRINT*,'X= ****************************************************'
        CALL lubksb(a,n,np,b)
        DO i=1,n
          WRITE(*,*) b(i)
        ENDDO
        PRINT *,'******************************************************'
      READ(*,*)
      STOP 'Finishing program'
      END
C     ******************************************************************
      SUBROUTINE ludcmp(a,n,np)
        INTEGER n,np,i,j
        REAL a(np,np),sum
        IF(a(1,1).EQ.0.)STOP 'Factorization Imposible'
        DO j=1,n
          DO i=1,j
            sum=a(i,j)
            DO k=1,i-1
              sum=sum-a(i,k)*a(k,j)
            ENDDO
            a(i,j)=sum
          ENDDO
          DO i=j+1,n
            sum=a(i,j)
            DO k=1,j-1
              sum=sum-a(i,k)*a(k,j)
            ENDDO
            IF(abs(a(j,j)).EQ.0.)THEN
              STOP 'Factorization Impossible'
            ELSE
              a(i,j)=sum/a(j,j)
            ENDIF
         ENDDO
        ENDDO
       RETURN
      END
C     ******************************************************************
      SUBROUTINE lubksb(a,n,np,b)
       INTEGER n,np,i,j
       REAL a(np,np),b(np),sum
       DO i=1,n
         sum=b(i)
         DO j=1,i-1
           sum=sum-a(i,j)*b(j)
         ENDDO
         b(i)=sum
       ENDDO
       DO i=n,1,-1
         sum=b(i)
         DO j=i+1,n
           sum=sum-a(i,j)*b(j)
         ENDDO
         b(i)=sum/a(i,i)
       ENDDO
       RETURN
      END
      
