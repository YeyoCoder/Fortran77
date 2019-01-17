!23456
      PROGRAM Main
        INTEGER m,mp,n,np
        REAL a(10,10),u(10,10),v(10,10),w(10),b(10),x(10)
        REAL wmax,wmin
        INTEGER i,j
        mp=10
        np=10
        n= 3
        m= 5
      
        a(1,1)=1.
        a(1,2)=0.
        a(1,3)=0.
        
        a(2,1)=1.
        a(2,2)=0.25
        a(2,3)=0.0625
        
        a(3,1)=1
        a(3,2)=0.5
        a(3,3)=0.25

        a(4,1)=1
        a(4,2)=0.75
        a(4,3)=0.5625
        
        a(5,1)=1
        a(5,2)=1
        a(5,3)=1
        
        b(1)=1.0000
        b(2)=1.2840
        b(3)=1.6487
        b(4)=2.1170
        b(5)=2.7183
        
        DO i=1,m
          WRITE(*,*)(a(i,j), j=1,n)
        ENDDO
        PRINT *,'******************************************************'
        DO i=1,n
          DO j=1,n
            u(i,j)=a(i,j)
          ENDDO
        ENDDO
        CALL svdcmp(u,m,n,mp,np,w,v)
        wmax=0.
        DO j=1,n
          IF(w(j).GT.wmax)wmax=w(j)
        ENDDO
        wmin=wmax*1.0e-6
        DO j=1,n
          if(w(j).LT.wmin)w(j)=0.
        ENDDO
        CALL svbksb(a,w,v,m,n,mp,np,b,x)
        PRINT *,'Matrix U  ********************************************'
        DO i=1,m
          WRITE(*,*)(u(i,j), j=1,n)
        ENDDO
        PRINT *,'Vector S**********************************************'
        DO i=1,n
          WRITE(*,*)(w(i))
        ENDDO
        PRINT *,'******************************************************'
        PRINT *,'The answer vector X is'
        DO i=1,n
          WRITE(*,*)(x(i))
        ENDDO
        PRINT *,'******************************************************'
      READ(*,*)
      STOP 'Finishing program...'
      END
C     ******************************************************************
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
        INTEGER m,mp,n,np,NMAX
        REAL a(mp,np),v(np,np),w(np)
        PARAMETER(NMAX=500)
C     Uses pythag
        INTEGER i,its,j,jj,k,l,nm
        REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
        g=0.0
        scale=0.0
        anorm=0.0
        DO i=1,n
          l=i+1
          rv1(i)=scale*g
          g=0.0
          s=0.0
          scale=0.0
          IF(i.LE.m)THEN
            DO k=i,m
              scale=scale+abs(a(k,i))
            ENDDO
            IF(scale.NE.0.0)THEN
              DO k=i,m
                a(k,i)=a(k,i)/scale
                s=s+a(k,i)*a(k,i)
              ENDDO
              f=a(i,i)
              g=-sign(sqrt(s),f)
              h=f*g-s
              a(i,i)=f-g
              DO j=l,n
                s=0.0
                DO k=i,m
                  s=s+a(k,i)*a(k,j)
                ENDDO
                f=s/h
                DO k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
                ENDDO
              ENDDO
              DO k=i,m
                a(k,i)=scale*a(k,i)
              ENDDO
            ENDIF
          ENDIF
          w(i)=scale*g
          g=0.0
          s=0.0
          scale=0.0
          IF((i.LE.m).AND.(i.NE.n))THEN
            DO k=l,n
              scale=scale+abs(a(i,k))
            ENDDO
            IF(scale.NE.0.0)THEN
              DO k=l,n
                a(i,k)=a(i,k)/scale
                s=s+a(i,k)*a(i,k)
              ENDDO
              f=a(i,l)
              g=-sign(sqrt(s),f)
              h=f*g-s
              a(i,l)=f-g
              DO k=l,n
                rv1(k)=a(i,k)/h
              ENDDO
              DO j=l,m
                s=0.0
                DO k=1,n
                  a(j,k)=a(j,k)+s*rv1(k)
                ENDDO
              ENDDO
              DO k=l,n
                a(i,k)=scale*a(i,k)
              ENDDO
            ENDIF
          ENDIF
          anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
        ENDDO
        DO i=n,1,-1
          IF(i.LT.n)THEN
            IF(g.NE.0.0)THEN
              DO j=l,n
                v(j,i)=(a(i,j)/a(i,l))/g
              ENDDO
              DO j=l,n
                s=0.0
                DO k=l,n
                  s=s+a(i,k)*v(k,j)
                ENDDO
                DO k=l,n
                  v(k,j)=v(k,j)+s*v(k,i)
                ENDDO
              ENDDO
            ENDIF
            DO j=l,n
              v(i,j)=0.0
              v(j,i)=0.0
            ENDDO
          ENDIF
          v(i,i)=1.0
          g=rv1(i)
          l=i
        ENDDO
        DO i=min(m,n),1,-1
          l=i+1
          g=w(i)
          DO j=l,n
            a(i,j)=0.0
          ENDDO
          IF(g.NE.0.0)THEN
            g=1.0/g
            DO j=l,n
              s=0.0
              DO k=l,m
                s=s+a(k,i)*a(k,j)
              ENDDO
              f=(s/a(i,i))*g
              DO k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              ENDDO
            ENDDO
            DO j=i,m
              a(j,i)=a(j,i)*g
            ENDDO
          ELSE
            DO j=i,m
              a(j,i)=0.0
            ENDDO
          ENDIF
          a(i,i)=a(i,i)+1.0
        ENDDO
        DO k=n,1,-1
          DO its=1,30
            DO l=k,1,-1
              nm=l-1
              IF((abs(rv1(l))+anorm).EQ.anorm) goto 2
              IF((abs(w(nm))+anorm).EQ.anorm) goto 1
            ENDDO
    1       c=0.0
            s=1.0
            DO i=l,k
              f=s*rv1(i)
              IF((abs(f)+anorm).EQ.anorm) goto 2
              g=w(i)
              h=pythag(f,g)
              w(i)=h
              h=1.0/h
              c=(g*h)
              s=-(f*h)
              DO j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
              ENDDO
            ENDDO
    2       z=w(k)
            IF(l.EQ.k)THEN
              IF(z.LT.0.0)THEN
                w(k)=-z
                DO j=1,n
                  v(j,k)=-v(j,k)
                ENDDO
              ENDIF
              GOTO 3
            ENDIF
            IF(its.EQ.30) THEN
              STOP 'No convergence in svdcmp'
              READ(*,*)
            ENDIF
            x=w(l)
            nm=k-1
            y=w(nm)
            g=rv1(nm)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g=pythag(f,1.0)
            f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c=1.0
            s=1.0
            DO j=l,nm
              i=j+1
              g=rv1(i)
              y=w(i)
              h=s*g
              g=c*g
              z=pythag(f,h)
              rv1(j)=z
              c=f/z
              s=h/z
              f=(x*c)+(g*s)
              g=-(x*s)+(g*c)
              h=y*s
              y=y*c
              DO jj=1,n
                x=v(jj,j)
                z=v(jj,i)
                v(jj,j)=(x*c)+(z*s)
                v(jj,i)=-(x*s)+(z*c)
              ENDDO
              z=pythag(f,h)
              w(j)=z
              IF(z.NE.0.0)THEN
                z=1.0/z
                c=f*z
                s=h*z
              ENDIF
              f= (c*g)+(s*y)
              x=-(s*g)+(c*y)
              DO jj=1,m
                y=a(jj,j)
                z=a(jj,i)
                a(jj,j)= (y*c)+(z*s)
                a(jj,i)=-(y*s)+(z*c)
              ENDDO
            ENDDO
            rv1(l)=0.0
            rv1(k)=f
            w(k)=x
          ENDDO
    3     CONTINUE
        ENDDO
       RETURN
      END
C     ******************************************************************
      FUNCTION pythag(a,b)
        REAL a,b,pythag
        REAL absa,absb
        absa=abs(a)
        absb=abs(b)
        IF(absa.GT.absb)THEN
          pythag=absa*sqrt(1.+(absb/absa)**2)
        ELSE
          IF(absb.EQ.0.)THEN
            pythag=0.
          ELSE
            pythag=absb*sqrt(1.+(absa/absb)**2)
          ENDIF
        ENDIF
       RETURN
      END
C     ******************************************************************
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
        INTEGER m,mp,n,np,NMAX
        REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
        PARAMETER(NMAX=500)
        INTEGER i,j,jj
        REAL s,tmp(NMAX)
        DO  j=1,n
          s=0.
          IF(w(j).NE.0.)THEN
            DO i=1,m
              s=s+u(i,j)*b(i)
            ENDDO
            s=s/w(j)
          ENDIF
          tmp(j)=s
        ENDDO
        DO j=1,n
          s=0.
          DO jj=1,n
            s=s+v(j,jj)*tmp(jj)
          ENDDO
          x(j)=s
        ENDDO
       RETURN
      END
        

