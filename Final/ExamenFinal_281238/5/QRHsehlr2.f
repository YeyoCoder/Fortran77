      PROGRAM Main
        INTEGER m,n,np,mp,i,j
        REAL a(10,10),b(10)
        m=4
        n=3
        mp=10
        np=10

        a(1,1)=1.
        a(2,1)=7.
        a(3,1)=4.
        a(4,1)=1.

        a(1,2)=2.
        a(2,2)=6.
        a(3,2)=4.
        a(4,2)=0.

        a(1,3)=2.
        a(2,3)=10.
        a(3,3)=6.
        a(4,3)=1.

        b(1)=6.
        b(2)=6.
        b(3)=8.
        b(4)=3.
        PRINT *,'+----------------------------------------------------+'
        DO i=1,m
            WRITE(*,*)(a(i,j),j=1,n)
        ENDDO
        PRINT *,'+----------------------------------------------------+'

        CALL HsehlrQR(a,m,n,mp,np)
        PRINT *,'+----------------------------------------------------+'
        DO i=1,m
            WRITE(*,*)(a(i,j),j=1,n)
        ENDDO
        PRINT *,'+----------------------------------------------------+'
      READ(*,*)
      STOP 'Finishing program...'
      END
C***********************************************************************
      SUBROUTINE HsehlrQR(a,m,n,mp,np)
       REAL a(mp,np)
       INTEGER m,n,mp,np
       REAL q(mp,mp),r(mp,np),x(mp),w(mp),tmp(mp),norm(np),wR(np)
       REAL tauW(mp),mtxpdt(np,mp)
       INTEGER i,j,k,MAXNORMCOL,ii
       REAL sum,normx,s,u1,tau,MAXNORM
C     PIVOTEO ¨Es necesario? TODO
        MAXNORM=0.
        DO i=1,n
          sum=0.
          DO j=1,m
            sum=sum+a(j,i)**2
          ENDDO
          norm(i)=sqrt(sum)
          PRINT *,'Norm=',norm(i),' i=',i
          IF(norm(i).GT.MAXNORM) THEN
            MAXNORM=norm(i)
            MAXNORMCOL=i
          ENDIF
        ENDDO
C       Exchange of rows
        DO i=1,m
          tmp(i)=a(i,1)
          a(i,1)=a(i,MAXNORMCOL)
          a(i,MAXNORMCOL)=tmp(i)
        ENDDO
        
        PRINT *,'The col with max norm is ',MAXNORMCOL
        PRINT*,'Matrix R ---------------------------------------------+'
        DO i=1,m
         DO j=1,n
           r(i,j)=a(i,j)
         ENDDO
         WRITE(*,*)(r(i,j),j=1,n)
        ENDDO
        
        PRINT *,'Matrix Q --------------------------------------------+'
        DO i=1,m
          DO j=1,m
            IF(i.EQ.j)THEN
              q(i,j)=1.
            ELSE
              q(i,j)=0.
            ENDIF
           ENDDO
           WRITE(*,*)(q(i,j),j=1,m)
        ENDDO
        

        DO j=1,n
        PRINT*,j,'----------------------------------------------------+'
          sum=0.
          DO i=j,m
            x(i)=r(i,j)
            sum=sum+x(i)**2
          ENDDO
          normx=sqrt(sum)
          PRINT *,'||x||=',normx
C       SIGN
          s=-sign(1.,r(j,j))
          PRINT *,'s=',s
          u1=r(j,j)-s*normx
          PRINT *,'u1= ',u1
          PRINT *,'w= ------------------------------------------------+'
          DO i=k,m
            w(i)=r(i,j)/u1
            PRINT*,w(i)
          ENDDO
          w(1)=1.
          tau=-s*u1/normx
          PRINT *,'tau=',tau
C         R= & q

C         Multiply  w'*R(j:end,:)
          PRINT *,'wtranp*R(j:end,:)----------------------------------+'
          DO i=j,n
            sum=0.
              DO k=1,m
                sum=sum+(w(k)*r(k,i))
              ENDDO
            wR(i)=sum
            WRITE(*,*) wR(i)
          ENDDO
          PRINT *,'tau*W ---------------------------------------------+'
          DO i=j,m
            tauW(i)=tau*w(i)
            PRINT *,tauW(i)
          ENDDO
          PRINT *,'mtxpdt---------------------------------------------+'
C         Multiply vectors to ger a Matrix
          DO i=j,m
            DO k=1,n
              sum=0.
              sum=sum+(tauW(i)*wR(k))
              mtxpdt(i,k)=sum
            ENDDO
          ENDDO
C         Substrac matrix to get R
          PRINT *,'R -------------------------------------------------+'
          DO i=j,m
            DO k=1,n
              r(i,k)=r(i,k)-(mtxpdt(i,k))
            ENDDO
            WRITE(*,*)(r(i,k),k=1,n)
          ENDDO

          PRINT *,'Q -------------------------------------------------+'
          DO i=1,m
           DO k=j,m
             q(i,k)=q(i,k)-((q(i,k)*w(k))*(tau*w(k)))
           ENDDO
           WRITE(*,*)(q(i,k),k=j,m)
             ENDDO
         PRINT *,'R -FINAL   -----------------------------------------+'



        ENDDO
        DO i=1,m
             WRITE(*,*)(r(i,j),j=1,n)
         ENDDO
       RETURN
      END
C     Function Norm



