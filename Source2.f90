
	SUBROUTINE NEWCOT2_gfm(H,Y,N1,N2,JG)	
	REAL*8 YR(N2),YI(N2),H,JGR,JGI
	COMPLEX*16 Y(N2),JG
	DO I=N1,N2	
	YR(I)=DREAL(Y(I))
	YI(I)=DIMAG(Y(I))	
	END DO
	JG=DCMPLX(0.D0,0.D0)
	IF(N2-N1.EQ.0)  JG=DCMPLX(0.D0,0.D0)
	IF(N2-N1.EQ.1)  JG=DCMPLX(H*(YR(N2)+YR(N1))/2.D0,&
	 & H*(YI(N2)+YI(N1))/2.D0)
	IF((N2-N1.GE.2)) THEN
	IF(MOD(N2-N1,2).EQ.0) THEN
	  DO I=N1,N2-2,2
	     JG=JG+DCMPLX(H*2.D0/6.D0*(1.D0*YR(I)+4.D0*YR(I+1)+YR(I+2)),&
		 & H*2.D0/6.D0*(1.D0*YI(I)+4.D0*YI(I+1)+YI(I+2)))
	  END DO
	ELSE

      IF(N2-N1.GE.5) THEN
	  DO I=N1,N2-3-2,2
	     JG=JG+DCMPLX(H*2.D0/6.D0*(1.D0*YR(I)+4.D0*YR(I+1)+YR(I+2)),&
		 & H*2.D0/6.D0*(1.D0*YI(I)+4.D0*YI(I+1)+YI(I+2)))
	  END DO
	  ELSE
	  END IF
        JG=JG+DCMPLX(H*3.D0/8.D0*(YR(N2-3)+3.D0*YR(N2-2)+3.D0*YR(N2-1)+YR(N2)),&
		 H*3.D0/8.D0*(YI(N2-3)+3.D0*YI(N2-2)+3.D0*YI(N2-1)+YI(N2)))
	END IF
	END IF	
	RETURN
	END

	SUBROUTINE NEWCOT4_gfm(H,Y,N1,N2,JG)	
	REAL*8 YR(N2),YI(N2),H,JGR,JGI
	COMPLEX*16 Y(N2),JG
	DO I=N1,N2	
	YR(I)=DREAL(Y(I))
	YI(I)=DIMAG(Y(I))	
	END DO
      JG=0.D0
	  DO I=N1,N2-4,4
	     JG=JG+DCMPLX(H*4.D0/90.D0*(7.D0*YR(I)+32.D0*YR(I+1)+12.D0*YR(I+2)+32.D0*YR(I+3)+7.D0*YR(I+4)),&
		 & H*4.D0/90.D0*(7.D0*YI(I)+32.D0*YI(I+1)+12.D0*YI(I+2)+32.D0*YI(I+3)+7.D0*YI(I+4)))
	  END DO
	RETURN
	END



	PARAMETER(NT=13200,NK=1000)

	REAL*8 TIM(NT),AT(NT),T,TTAU,SOMGA,HT,A1,A2,A3,PI,HHG(NK),OMGA(NK),OMGA0 &
	        & ,OMGAL,OMGAR,HOMGA,ASW1,ASW2,FAI,NFAI,ATX(NT),ATY(NT),HHGX(NK),HHGY(NK)
	COMPLEX*16 ATTU(NT),TE(NT),AS(NT),PWT(NK),TK(NK),AS2
	COMPLEX*16 IMAGINARY
	IMAGINARY=DCMPLX(0.D0,1.D0)
	PI=3.1415926535897932384d0
	
    
	OPEN(FILE='r2(t).DAT',UNIT=52)
	OPEN(FILE='OUT2.DAT',UNIT=53)
	OPEN(FILE='HHG.DAT',UNIT=54)
    OPEN(FILE='AS(t).DAT',UNIT=55)
    OPEN(FILE='AS2(t).DAT',UNIT=56)
    OPEN(FILE='ACC-READ-WRITE.DAT',UNIT=57)
     OPEN(FILE='X(t).DAT',UNIT=58)
     OPEN(FILE='HHG-ellipticity.DAT',UNIT=59)
	TTAU=10.D0
	
	NW=16
	do I=1,NT
	READ(52,*)A1,A2,A3
	TIM(I)=A1
	AT(I)=A2*(1.D0-DCOS(DBLE(I-1)*2.D0*PI/DBLE(NT-1)))/2.D0
	ATX(I)=A2*(1.D0-DCOS(DBLE(I-1)*2.D0*PI/DBLE(NT-1)))/2.D0
	ATY(I)=A3*(1.D0-DCOS(DBLE(I-1)*2.D0*PI/DBLE(NT-1)))/2.D0	
!	IF(I.LE.INT(NT/NW)) AT(I)=AT(I)*DEXP(-4.D0*DLOG(2.D0)*DBLE(I-INT(NT/NW))**2.D0/(DBLE(INT(NT/NW))/4.D0)**2.D0)
!	IF(I.GE.(NT-INT(NT/NW)+1)) AT(I)=AT(I)*DEXP(-4.D0*DLOG(2.D0)*DBLE(I-(NT-INT(NT/NW)+1))**2.D0/(DBLE(INT(NT/NW))/4.D0)**2.D0)
	WRITE(57,'(1X,5(1X,F15.8))')TIM(I),A2,ATX(I),A3,ATY(I)
	END DO
!	stop

	OMGA0=0.057D0
	HT=TIM(2)-TIM(1)
	OMGAL=0.0034d0 !-PI/HT
	OMGAR=350.d0*OMGA0 !PI/HT
	HOMGA=(OMGAR-OMGAL)/DBLE(NK-1)
    DO K=1,NK
    OMGA(K)=OMGAL+DBLE(K-1)*HOMGA
    END DO

!    goto 11
! Ð¡²¨±ä»»----------------------------------------------------
    NH=16
	DO K=1,nk,2


	DO I=1,NT,NH
    NIH=INT(3.D0*TTAU*DSQRT(2.D0*4.D0*DLOG(2.D0))/DABS(OMGA(K))/HT)
	NI1=MAX(1,I-NIH)
	NI2=MIN(NT,I+NIH)

	DO J=NI1,NI2
	TE(J)=DSQRT(DABS(OMGA(K)))*AT(J)*(1.D0/DSQRT(TTAU))* &
     & CDEXP(IMAGINARY*OMGA(K)*(TIM(J)-TIM(I))) &
     & *DEXP(-(OMGA(K)*(TIM(J)-TIM(I)))**2.D0/(2.D0*TTAU**2.D0))
 !  write(*,"(I2,5X,D20.12,5X,D20.12)") J,DREAL(TE(J)),DIMAG(TE(J))
	END DO
	CALL NEWCOT2_gfm(HT,TE,NI1,NI2,ATTU(I))

	WRITE(53,*) omga(k)/0.057D0, tim(i)*0.0242d0,(CDABS(ATTU(I)))**2.D0

	END DO
!	IF(MOD(K,20).EQ.0)write(*,*)'wavelet',k
	END DO
	stop
! Ð¡²¨±ä»»----------------------------------------------------
 !   GOTO 13
 
!  HHG--------------------------------------------------------
11   NFAI=0.D0
     DO K=1,NK
     DO I=1,NT
	 TE(I)=AT(I)*CDEXP(IMAGINARY*OMGA(K)*TIM(I))/DSQRT(2.D0*PI)
	 END DO
	 CALL NEWCOT2_gfm(HT,TE,1,NT,PWT(K))
	 HHG(K)=(CDABS(PWT(K))/(TIM(NT)-TIM(1))/(OMGA(K))**2.D0)**2.D0
	 
	 IF(K.GT.1) THEN
	 IF(DREAL(PWT(K-1)).GT.0.D0.AND.DIMAG(PWT(K-1)).GT.0.D0 & 
	    & .AND.DREAL(PWT(K)).LT.0.D0.AND.DIMAG(PWT(K)).GT.0.D0) NFAI=NFAI+1.D0
	 IF(DREAL(PWT(K-1)).LT.0.D0.AND.DIMAG(PWT(K-1)).LT.0.D0 & 
	    & .AND.DREAL(PWT(K)).GT.0.D0.AND.DIMAG(PWT(K)).LT.0.D0) NFAI=NFAI+1.D0
	 IF(DREAL(PWT(K-1)).LT.0.D0.AND.DIMAG(PWT(K-1)).GT.0.D0 & 
	    & .AND.DREAL(PWT(K)).GT.0.D0.AND.DIMAG(PWT(K)).GT.0.D0) NFAI=NFAI-1.D0
	 IF(DREAL(PWT(K-1)).GT.0.D0.AND.DIMAG(PWT(K-1)).LT.0.D0 & 
	    & .AND.DREAL(PWT(K)).LT.0.D0.AND.DIMAG(PWT(K)).LT.0.D0) NFAI=NFAI-1.D0
     ELSE
     END IF
	    	    
	 FAI=NFAI*PI+DATAN(DIMAG(PWT(K))/DREAL(PWT(K)))
     WRITE(54,*) OMGA(K)/OMGA0,DLOG10(HHG(K)),FAI/PI
	if(mod(k,20).eq.0) write(*,*)'hhg',k
   END DO
!  HHG--------------------------------------------------------
   
!   GOTO 12

!  °¢ÃëÂö³å---------------------------------------------------
     ASW1=45.D0*OMGA0
     ASW2=60.D0*OMGA0
 !    DO I=1,NT
 !     DO K=1,NK
 !	   TK(K)=PWT(K)*CDEXP(-IMAGINARY*OMGA(K)*TIM(I))/DSQRT(2.D0*PI) &
 !        &   *DEXP(-4.D0*DLOG(2.D0)*(OMGA(K)-(ASW1+ASW2)/2.D0/27.2D0)**2.D0/((ASW2-ASW1)/27.2D0)**2.D0)
 !	  END DO
!	  CALL NEWCOT2_gfm(HOMGA,TK,1,NK,AS(I))
 !     WRITE(55,*) TIM(I),DREAL(AS(I))**2.D0+DIMAG(AS(I))**2.D0
  !   END DO
     
   DO I=1,NT
    AS(I)=0.D0
    DO K=INT((ASW1-OMGAL)/HOMGA),INT((ASW2-OMGAL)/HOMGA)
    AS(I)=AS(I)+PWT(K)*CDEXP(-IMAGINARY*OMGA(K)*TIM(I))
    TK(K)=PWT(K)*CDEXP(-IMAGINARY*OMGA(K)*TIM(I))/DSQRT(2.D0*PI)
    END DO
    CALL NEWCOT2_gfm(HOMGA,TK,INT((ASW1-OMGAL)/HOMGA),INT((ASW2-OMGAL)/HOMGA),AS2)
    WRITE(56,*) TIM(I),DREAL(AS(I))**2.D0+DIMAG(AS(I))**2.D0,DREAL(AS2)**2.D0+DIMAG(AS2)**2.D0
   END DO
   
   GOTO 13

!  °¢ÃëÂö³å---------------------------------------------------
     
 12    DO K=1,NK,5
 
     DO J=NT/3,NT-NT/3,21
     
     DO I=1,NT
	 TE(I)=ATX(I)*CDEXP(IMAGINARY*OMGA(K)*TIM(I))/DSQRT(2.D0*PI)*DEXP(-4.D0*DLOG(2.D0)*(TIM(I)-TIM(J))**2.D0/13.7789D0**2.D0)
	 END DO
	 CALL NEWCOT2_gfm(HT,TE,1,NT,PWT(K))
	 HHGX(K)=(CDABS(PWT(K))/(TIM(NT)-TIM(1))/(OMGA(K))**2.D0)
	 
     DO I=1,NT
	 TE(I)=ATY(I)*CDEXP(IMAGINARY*OMGA(K)*TIM(I))/DSQRT(2.D0*PI)*DEXP(-4.D0*DLOG(2.D0)*(TIM(I)-TIM(J))**2.D0/13.7789D0**2.D0)
	 END DO
	 CALL NEWCOT2_gfm(HT,TE,1,NT,PWT(K))
	 HHGY(K)=(CDABS(PWT(K))/(TIM(NT)-TIM(1))/(OMGA(K))**2.D0)
	 
	 WRITE(59,*) OMGA(K)/0.057D0,TIM(J),HHGY(K)/DSQRT(HHGY(K)**2.D0+HHGX(K)**2.D0)

	 END DO
	 WRITE(*,*) K,'OK'
     END DO
 13    pause
    stop

 	close(52)


	END



	 
