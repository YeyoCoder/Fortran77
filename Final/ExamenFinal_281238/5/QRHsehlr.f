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
        
        CALL qrhsehlr(a,m,n,mp,np,b)
        PRINT *,'+----------------------------------------------------+'
        DO i=1,m
            WRITE(*,*)(a(i,j),j=1,n)
        ENDDO
        PRINT *,'+----------------------------------------------------+'
      READ(*,*)
      STOP 'Finishing program...'
      END
C     ******************************************************************
      SUBROUTINE qrhsehlr(a,m,n,mp,np,b)
        REAL a(mp,np),b(mp)
        INTEGER m,n,mp,np
        INTEGER i,j,MAXNORMCOL,l,k
        REAL MAXNORM,sum,h(mp,mp),x(mp),u(mp),e(mp),normx,normu
        REAL p(np,np),norm(np),tmp(np),r(mp,np),q(mp,mp),hh(mp,np)
        REAL s(mp,np)
C       Intercambio de renglon pivot
        MAXNORM=0.
        DO i=1,n
          sum=0.
          DO j=1,m
            sum=sum+a(j,i)**2
            PRINT *,'i=',i,'j=',j,'sum=',sum
          ENDDO
          norm(i)=sqrt(sum)
          PRINT *,'Norm=',norm(i),' i=',i
          IF(norm(i).GT.MAXNORM) THEN
            MAXNORM=norm(i)
            MAXNORMCOL=i
          ENDIF
        ENDDO
        DO i=1,n
          DO j=1,n
            IF(i.EQ.j)THEN
              p(i,j)=1.
            ELSE
              p(i,j)=0.
            ENDIF
          ENDDO
        ENDDO
        DO i=1,n
          WRITE(*,*) (p(i,j),j=1,n)
        ENDDO
        PRINT *,'The col with max norm is ',MAXNORMCOL
C       Exchange of rows

        DO i=1,m
          tmp(i)=a(i,1)
          a(i,1)=a(i,MAXNORMCOL)
          a(i,MAXNORMCOL)=tmp(i)
        ENDDO
        
        DO i=1,n
          tmp(i)=p(i,1)
          p(i,1)=p(i,MAXNORMCOL)
          p(i,MAXNORMCOL)=tmp(i)
        ENDDO
        DO i=1,n
          WRITE(*,*) (p(i,j),j=1,n)
        ENDDO
C       Definimos las matrices R y Q
        PRINT *,'+R --------------------------------------------------+'
        DO i=1,m
          DO j=1,n
            r(i,j)=a(i,j)
          ENDDO
          WRITE(*,*) (r(i,j),j=1,n)
        ENDDO
        PRINT *,'+Q --------------------------------------------------+'
        DO i=1,m
          DO j=1,m
            IF(j.EQ.i)THEN
              q(i,j)=1.
            ELSE
              q(i,j)=0.
            ENDIF
          ENDDO
          WRITE(*,*) (q(i,j),j=1,n)
        ENDDO
C       Definimos sigma
        sigma=0.
        DO i=1,m
          sigma=sigma+a(i,1)**2
        ENDDO
        sigma=-1*sqrt(sigma)
        PRINT *,'sigma=',sigma
C       Algoritmo para Factorizacion QR
        DO k=1,n
          sum=0.
          DO j=k,m
            x(j)=r(j,k)
            sum=sum+x(j)**2
          ENDDO
          normx=sqrt(sum)
          PRINT *,'X= *************************************************'
          DO j=1,m
            PRINT*,x(j)
          ENDDO
          PRINT*,'norm= ',normx
C         creamos e1
          PRINT *,'e',k,'= ********************************************'
          DO j=1,m
            IF(k.EQ.j)THEN
              e(j)=1.
            ELSE
              e(j)=0.
            ENDIF
             PRINT *,e(j)
          ENDDO
          PRINT *,'1.- u',k,'= sign(x1)||x||e+x  **********************'
          DO j=k,m
            u(j)=(sign(1.,x(j))*normx*e(j))+x(j)
            PRINT*,u(j)
          ENDDO
C         Calculamos la norma de U
          sum=0.
          DO j=k,m
            sum=sum+u(j)**2
          ENDDO
          normu=sqrt(sum)
          PRINT *,'||u||=',normu
          
          PRINT *,'2.- u',k,'= ****************************************'
          PRINT *,'u/||u|| ********************************************'
          DO j=k,m
            u(j)=u(j)/normu
            PRINT *,u(j)
          ENDDO
          PRINT *,'uu*',k,'= ******************************************'
          DO i=k,m
            DO j=k,m
              h(i,j)=u(i)*u(j)
            ENDDO
            WRITE(*,*)(h(i,j),j=1,m)
          ENDDO
C       Multiplicarion of H-mxm . R-mxn
          PRINT *,'uu*.R',k,'= ***************************************'
          DO i=k,m
            DO j=k,n
              sum=0.
              DO l=k,m
                sum=sum+(h(i,l)*r(l,j))
              ENDDO
              hh(i,j)=sum
            ENDDO
            WRITE(*,*)(hh(i,j),j=1,n)
          ENDDO
          
          PRINT *,'R',k,'= *******************************************'
          DO i=k,m
            DO j=k,n
              IF(i.GT.j)THEN
                r(i,j)=0.
              ELSE
                r(i,j)=r(i,j)-(2*hh(i,j))
              ENDIF
            ENDDO
           WRITE(*,*)(r(i,j),j=k,n)
          ENDDO
          
          PRINT *,'Q(:,k.:m)uu* --------------------------------------+'
          DO i=1,m
            DO j=k,m
              sum=0.
              DO l=1,m
                sum=sum+(q(i,l)*hh(l,j))
              ENDDO
            q(i,j)=q(i,j)-(2*sum)
            ENDDO
          ENDDO
          PRINT *,'Q',k,'= ********************************************'

        ENDDO
        
        
        PRINT *,'R = **************************************************'
        DO i=1,m
          WRITE(*,*)(r(i,j),j=1,n)
        ENDDO
        
        PRINT *,'Q = **************************************************'
        DO i=1,m
          WRITE(*,*)(q(i,j),j=1,m)
        ENDDO
        PRINT *,'S = **************************************************'
        DO i=1,m
          DO j=1,n
            sum=0.
            DO k=1,m
              sum=sum+(q(i,k)*r(k,j))
            ENDDO
          s(i,j)=sum
          ENDDO
          WRITE(*,*)(s(i,j),j=1,n)
        ENDDO
        
       RETURN
      END

