      SUBROUTINE scrsho(fx)
      INTEGER ISCR,JSCR
      REAL fx
      EXTERNAL fx
C     NUMBER OF HORIZONTAL AND VERTICAL POSITIONS IN DISPLAY.
C     FOR INTERACTIVE CRT TERMINAL USE.
C     PRODUCE A CRUDE GRAPH OF THE FUNCTION fx OVER THE PROMPTED-FOR
C     INTERVAL  x1,x2.
C     QUERY FOR ANOTHER PLOT UNTIL THE USER SIGNALS SATISFACTION.
      PARAMETER (ISCR=60,JSCR=21)
      INTEGER i, j, jz
      REAL dx, dyj, x, x1, x2, ybig, ysml, y(ISCR)
      CHARACTER*1 scr(ISCR,JSCR),blank,zero,yy,xx,ff
      SAVE blank,zero,yy,xx,ff
      DATA blank,zero,yy,xx,ff/' ','-','1','-','x'/
 1    continue
      write(*,*)' Enter  x1,x2 (= to stop)'
C     QUERY FOR ANOTHER PLOT, QUIT IF x1=x2
      read(*,*)x1,x2
      if(x1.eq.x2)return
      do j=1,JSCR
        scr(1,j)=yy
        scr(ISCR,j)=yy
      enddo
      do i=2,ISCR-1
        scr(i,1)=xx
        scr(i,JSCR)=xx
        do j=2,JSCR-1
          scr(i,j)=blank
        enddo
      enddo
      dx=(x2-x1)/(ISCR-1)
      x=x1
      ybig=0.
      ysml=ybig
      do i=1, ISCR
        y(i)=fx(x)
        if(y(i).lt.ysml) ysml=y(i)
        if(y(i).gt.ybig) ybig=y(i)
        x=x+dx
      enddo
      if(ybig.eq.ysml) ybig=ysml+1.
      dyj=(JSCR-1)/(ybig-ysml)
      jz=1-ysml*dyj
      do i=1, ISCR
        scr(i,jz)=zero
        j=1+(y(i)-ysml)*dyj
        scr(i,j)=ff
      enddo
      write(*,'(1x,1pe10.3,1x,80a1)') ybig,(scr(i,JSCR),i=1,ISCR)
      do j=JSCR-1,2,-1
        write(*,'(12x,80a1)')(scr(i,j),i=1,ISCR)
      enddo
      write(*,'(1x,1pe10.3,1x,80a1)') ysml,(scr(i,1), i=1, ISCR)
      write(*,'(12x,1pe10.3,40x,e10.3)') x1,x2
      goto 1
      END
      
      
