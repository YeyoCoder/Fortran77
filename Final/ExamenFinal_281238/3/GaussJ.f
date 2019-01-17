!23456
      PROGRAM Main
        INTEGER n,np,m,mp
        REAL a(10,10), b(10,10)
        PRINT*,'Gauss-Jordan'
        PRINT*,'Input the number of ROWS and COLUMNS of the square matri
     &x'
        READ(*,*)n
        PRINT*,'Input the number of SOLUTION Vectors'
        READ(*,*)m
        np=10
        mp=10
        CALL inputMatrix(a,n,np,b,m,mp)
        CALL gaussj(a,n,np,b,m,mp)
        PRINT *,'The solution vector x is: '
        CALL showMatrix(b,n,m)
        READ(*,*)
       STOP 'Finishing program...'
      END
C     ******************************************************************
      SUBROUTINE gaussj(a,n,np,b,m,mp)
        INTEGER m,mp,n,np,NMAX
        REAL a(np,np), b(np,mp)
        PARAMETER(NMAX=50)
        INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
        REAL big,dum,pivinv
        DO j=1,n
          ipiv(j)=0
        ENDDO
        DO i=1,n
          big=0.
          DO j=1,n
            IF(ipiv(j).NE.1)THEN
              DO k=1,n
                IF(ipiv(k).EQ.0)THEN
                  IF(abs(a(j,k)).GE.big) THEN
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                  ENDIF
                ELSE IF(ipiv(k).GT.1)THEN
                  STOP 'singular matrix in gaussj'
                  READ(*,*)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          ipiv(icol)=ipiv(icol)+1
          IF(irow.NE.icol) THEN
            DO l=1,n
              dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
            ENDDO
            DO l=1,m
              dum=b(irow,l)
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
            ENDDO
          ENDIF
          indxr(i)=irow
          indxc(i)=icol
          IF(a(icol,icol).EQ.0.) THEN
            STOP 'singular matrix in gaussj'
            READ(*,*)
          ENDIF
          pivinv=1./a(icol,icol)
          a(icol,icol)=1.
          DO l=1,n
            a(icol,l)=a(icol,l)*pivinv
          ENDDO
          DO l=1,m
            b(icol,l)=b(icol,l)*pivinv
          ENDDO
          DO ll=1,n
            IF(ll.NE.icol)THEN
              dum=a(ll,icol)
              a(ll,icol)=0.
              DO l=1,n
                a(ll,l)=a(ll,l)-a(icol,l)*dum
              ENDDO
              DO l=1,m
                b(ll,l)=b(ll,l)-b(icol,l)*dum
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        DO l=n,1,-1
          IF(indxr(l).NE.indxc(l))THEN
            DO k=1,n
              dum=a(k,indxr(l))
              a(k,indxr(l))=a(k,indxc(l))
              a(k,indxc(l))=dum
            ENDDO
          ENDIF
        ENDDO
       RETURN
      END
C     ******************************************************************
      SUBROUTINE inputMatrix(a,n,np,b,m,mp)
        INTEGER n,m,np,mp
        REAL a(np,np), b(np,mp)
        INTEGER i,j
        PRINT *,'Input the elements of the matrix A'
        DO i=1,n
          DO j=1,n
            WRITE(*,'(1A1,A1,2I1,A3,$)')char(13),'a',i,j,' = '
            READ(*,*) a(i,j)
          ENDDO
        ENDDO
        PRINT*,'*******************************************************'
        DO i=1,n
          WRITE(*,*) (a(i,j),j=1,n)
        ENDDO
        PRINT*,'*******************************************************'
        PRINT *,'Input the solution vector'
        DO i=1,n
          DO j=1,m
            WRITE(*,'(1A1,A1,2I1,A3,$)')char(13),'b',i,j,' = '
            READ(*,*)b(i,j)
          ENDDO
        ENDDO
        CALL showMatrix(b,n,m)
       RETURN
      END
C     ******************************************************************
      SUBROUTINE showMatrix(a,n,m)
        INTEGER n,m
        REAL a(n,m)
        INTEGER i,j
        PRINT*,'*******************************************************'
        DO i=1,n
          WRITE(*,*) (a(i,j),j=1,m)
        ENDDO
        PRINT*,'*******************************************************'
       RETURN
      END
