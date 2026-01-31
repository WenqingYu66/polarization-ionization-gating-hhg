SUBROUTINE NEWCOT2_R_gfm(H,Y,N1,N2,JG) !实数数值积分，Y为被积函数，N1为下限对应的格点，N2为上限对应的格点，JG结果	
	REAL*8 YR(N2),H,JG,Y(N2)
	DO I=N1,N2	
	YR(I)=Y(I)	
	END DO
	JG=0.D0
	IF(N2-N1.EQ.0)  JG=0.D0
	IF(N2-N1.EQ.1)  JG=H*(YR(N2)+YR(N1))/2.D0
	IF((N2-N1.GE.2)) THEN
	IF(MOD(N2-N1,2).EQ.0) THEN
	  DO I=N1,N2-2,2
	     JG=JG+H*2.D0/6.D0*(1.D0*YR(I)+4.D0*YR(I+1)+YR(I+2))
	  END DO
	ELSE
      IF(N2-N1.GE.5) THEN
	  DO I=N1,N2-3-2,2
	     JG=JG+H*2.D0/6.D0*(1.D0*YR(I)+4.D0*YR(I+1)+YR(I+2))
	  END DO
	  ELSE
	  END IF
        JG=JG+H*3.D0/8.D0*(YR(N2-3)+3.D0*YR(N2-2)+3.D0*YR(N2-1)+YR(N2))
	END IF
	END IF	
	RETURN
	END
																		    



	SUBROUTINE NEWCOT2_gfm(H,Y,N1,N2,JG) !复数数值积分，Y为被积函数，N1为下限对应的格点，N2为上限对应的格点，JG结果	
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
	

	SUBROUTINE DS(H,HS,N,DS1)
	COMPLEX*16 HS(N),DS1(N)
	REAL*8 H

	DS1(1)=(HS(2)-HS(1))/H
	DS1(N)=(HS(N)-HS(N-1))/H

	DO I=2,N-1
	DS1(I)=(HS(I+1)-HS(I-1))/2.D0/H
	END DO
	RETURN
	END




	FUNCTION D(PX,PY,PZ,ALFA)
	REAL*8 PX,PY,PZ,ALFA,PI
	COMPLEX*16 D,IMAGINARY
	IMAGINARY=DCMPLX(0.D0,1.D0)
	pi=3.1415926535897932384d0
	D=IMAGINARY*2.D0**(7.D0/2.D0)*(ALFA)**(5.D0/4.D0)/PI*(PX)/(PX**2.D0+PY**2.D0+PZ**2.D0+ALFA)**3.D0
 	RETURN
	END
	
!	FUNCTION D_1(PX,PY,ALFA)
!	REAL*8 PX,PY,ALFA,PI
!	COMPLEX*16 D,IMAGINARY
!	IMAGINARY=DCMPLX(0.D0,1.D0)
!	pi=3.1415926535897932384d0
!	D=IMAGINARY*PX/(PI*ALFA)**(3.D0/4.D0)/ALFA*DEXP(-1.D0*(PX**2.D0+PY**2.D0+ALFA)/2.D0/ALFA)
 !	RETURN
!	END


!	program main
!	implicit none																			   
    PARAMETER (N=13200)
    REAL*8 HT,t0,t1,t(N),pi,UPX,UPY,UPZ,IP,EX0,EY0,EZ0,EX(N),EY(N),EZ(N),OMGAX,OMGAY,OMGAZ,SST,PSTX,PSTY,PSTZ,&
        	& XT(N),F(N),FWHMX,FWHMY,FWHMZ,AX(N),AY(N),AZ(N),A2(N),SMALL,ALFA,DEATA,JG,JG1,RE(N),YT(N),GAMMAX,&
			& GAMMAY,GAMMAZ,JG3,JG5,CN ,epsilon,ax0,ay0,az0,ZT(N),ell,cnx,cny,tao,tim(n),taox,thita,kesai,td,cep,&
			&tplus(n),tminus(N),E0,C0(N),JG0,CNL,GLM,F0,NSTAR,M,WADKx(N),wadky(n),gt(n),taol(n),etot(n),tstar,differ_p,startphase
	COMPLEX*16,EXTERNAL:: D 
	COMPLEX*16 IMAGINARY,CMPL3X(N),CMPL3Y(N),CMPL3Z(N),JG2,JG4,axc(n),ayc(n),ds1(n),JG6
	INTEGER CPUNUMBER,N1,N2
	CPUNUMBER=4
	PI=3.1415926535897932384d0
	IMAGINARY=DCMPLX(0.D0,1.D0)
    
	IP=0.9D0
    ell=1.d0
	Ex0=0.28D0
	EY0=0.28d0
	EZ0=0.d0

	OMGAX=0.057D0
	OMGAY=0.057D0
	OMGAZ=0.057D0


	CN=12.D0

    epsilon=0.1d0

	SMALL=1.0D-16
	ALFA=2.D0*IP


	T0=-2.d0*pi/omgax*CN/2.D0
	T1=2.d0*pi/omgax*CN/2.D0
		HT=(T1-T0)/(N-1.D0)
	  DO I=1,N
	   T(I)=HT*(I-1.D0)+T0
	  END DO



	  OPEN(FILE='E(t).DAT',UNIT=9)
	  OPEN(FILE='A(t).DAT',UNIT=11)
	  OPEN(FILE='r2(t).DAT',UNIT=12)
	  OPEN(FILE='C02.DAT',UNIT=15)

	   

	  N1=INT(2.D0*PI/OMGAX/HT*2.d0)
	  N2=N-N1
      DO I=1,N1
	  F(I)=(T(I)-T(1))/(T(N1)-T(1))
	  END DO
      DO I=N1+1,N2
	  F(I)=1.D0
	  END DO
      DO I=N2+1,N
	  F(I)=(T(N)-T(I))/(T(N)-T(N2))
	  END DO
     

!     FWHMX=2.d0*pi/omgax*6.d0
!	  FWHMY=2.D0*pi/omgay*6.d0
!	  FWHMZ=2.D0*pi/omgaZ*6.d0
! ----------------------------------------------------------------------------------------------------

	
	CNX=12.D0

     TAOX=2.D0*PI/OMGAX*1.5d0                ! FWHMX=1.5T
	
 
	write(*,*)taox*0.0242d0
	
    cep=pi/5.d0
	td=2.D0*PI/OMGAX*1.25d0
    differ_p=td*OMGAX                       !  phase difference
    startphase=PI*0.5d0                     !  start phase

	DO I=1,N
   EX(I)=EX0*dcos(omgax*t(i)+differ_p+startphase)*DEXP(-2.d0*DLOG(2.D0)*(T(I)+td/2.d0)**2.D0/(taox**2.D0))                                     ! O light fast 
!    EX(I)=EX0*dcos(omgax*t(i))*DEXP(-2.d0*DLOG(2.D0)*(T(I))**2.D0/(taox**2.D0)) 
	 EY(I)=EY0*dcos(omgaY*t(i)+startphase)*DEXP(-2.d0*DLOG(2.D0)*(T(I)-td/2.d0)**2.D0/(taoX**2.D0))                                             ! E light slow
	 ETOT(I)=DEXP(-2.d0*DLOG(2.D0)*(T(I)-td/2.d0)**2.D0/(taox**2.D0))/DEXP(-2.d0*DLOG(2.D0)*(T(I)+td/2.d0)**2.D0/(taoX**2.D0))
    
	WRITE(9,'(1X,4(1X,F15.7))') T(I),ex(i),EY(I),ETOT(I)
	   END DO
	 

!	DO I=N/3+1, N/3+N/3
   
!      EX(I)=EX0*DEXP(-2.d0*DLOG(2.D0)*(T(I)+td)**2.D0/(taox**2.D0))*dcos(omgax*t(i))
!	  EY(I)=EY0*DEXP(-2.d0*DLOG(2.D0)*(T(I)-td)**2.D0/(taoX**2.D0))*dSIN(omgaY*t(i))
    
!	WRITE(9,*) T(I),ex(i),EY(I)
!	   END DO

!      DO I=N/3+N/3+1,N
  
 !      EX(I)=EX0*DEXP(-2.d0*DLOG(2.D0)*T(I)**2.D0/(taox**2.D0))*dcos(omgax*t(i))
 !        EX(I)=EX0*DEXP(-2.d0*DLOG(2.D0)*(T(I)+td)**2.D0/(taox**2.D0))*dcos(omgax*t(i))
!	   EY(I)=EY0*DEXP(-2.d0*DLOG(2.D0)*(T(I)-td)**2.D0/(taoX**2.D0))*dSIN(omgaY*t(i))
    
!	WRITE(9,*) T(I),ex(i),EY(I)
!	   END DO

!	   do i=1,n
!	   write(9,*) T(I),ex(i),EY(I)
!	   end do


!	   STOP


 
!---------------------------------------------------------------------------------------------------------


	  DO I=1,N
      RE(I)=-EX(I)
	  END DO
	  DO I=1,N
	  CALL NEWCOT2_R_gfm(HT,RE,1,I,JG)
      AX(I)=JG 
	  END DO
	  
	  DO I=1,N
     RE(I)=-EY(I)
	  END DO
	  DO I=1,N
	  CALL NEWCOT2_R_gfm(HT,RE,1,I,JG)
	  AY(I)=JG
	  END DO
 
	  DO I=1,N
     RE(I)=-EZ(I)
	  END DO
	  DO I=1,N
	  CALL NEWCOT2_R_gfm(HT,RE,1,I,JG)
	  AZ(I)=JG
	  END DO


     DO I=1,N
	 A2(I)=AX(I)**2.D0+AY(I)**2.D0+AZ(I)**2.D0

	 WRITE(11,'(1X,4(1X,F15.7))') T(I),AX(I),AY(I),AZ(i)
	 
!---------------------------------------------------------------------------------------------

	 END DO
	 
	
	  CNL=4.25575D0
      GLM=1.D0
    
      NSTAR=0.74387D0
 !     M=0.D0
     
	     DO I=1,N
  !    wADK(I)=CNL*IP*(2.D0*F0/dabs(EX(I)))**(2.D0*NSTAR-1.D0)*DEXP(-2.D0*F0/(3.D0*dabs(EX(I))))
      wadkx(i)=ip*cnl*glm*(4.d0*ip/(dabs(ex(i)))*dsqrt(2.d0*ip))**(2.d0*dsqrt(13.6d0/24.588d0)-1.d0)*&
      &dexp(-4.d0*ip*dsqrt(2.d0*ip)/(3.d0*(dabs(ex(i)))))
  !    wadk(i)=ip*cnl*(4.d0*ip/(dabs(ex(i)))*dsqrt(2.d0*ip))**(2*dsqrt(13.6d0/24.588d0)-1)&
  !  &*dexp(-4.d0*ip/3.d0/(dabs(ex(i)))*dsqrt(2.d0*ip))
  	  wadky(i)=ip*cnl*glm*(4.d0*ip/(dabs(ey(i)))*dsqrt(2.d0*ip))**(2.d0*dsqrt(13.6d0/24.588d0)-1.d0)*&
      &dexp(-4.d0*ip*dsqrt(2.d0*ip)/(3.d0*(dabs(ey(i)))))
	   call NEWCOT2_r_gfm(ht,wadkx,1,I,jg)
	   call NEWCOT2_r_gfm(ht,wadky,1,I,jg1)
	!	gt(i)=dexp(-JG)
	    gt(I)=dexp(-jg-jg1)
		
		write(15,111) t(I)*0.0242D0,ex(i),gt(I),dexp(-jg)
		

	    111 format(8(1x,e16.10))
!	    WRITE(*,*)JG
	  
       end do
!	  stop

 
	 
   	 XT(1)=0.D0
	 YT(1)=0.D0
	 ZT(1)=0.D0

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(J,I,JG,JG1,JG2,JG3,JG4,PSTX,PSTY,SST,CMPL3X,CMPL3Y)
	 DO J=1,N
      JG=0.D0
      JG3=0.D0
      JG1=0.D0	 
	  JG5=0.D0
	  JG0=0.D0 

	
	  DO I=J-1,1,-1

	  JG=JG+(AX(I)+AX(I+1))/2.D0*HT
	  PSTX=JG/(T(J)-T(I))
	  JG3=JG3+(AY(I)+AY(I+1))/2.D0*HT
	  PSTY=JG3/(T(J)-T(I))
	  JG5=JG5+(AZ(I)+AZ(I+1))/2.D0*HT
	  PSTZ=JG5/(T(J)-T(I))

	  JG1=JG1+(A2(I)+A2(I+1))/2.D0*HT 
!	  jg0=jg0+pstx*(AX(I)+AX(I+1))/2.D0*HT+psty*(Ay(I)+Ay(I+1))/2.D0*HT
	  SST=IP*(T(J)-T(I))-1.D0/2.D0*(PSTX**2.D0+PSTY**2.D0+PSTZ**2.D0)*(T(J)-T(I))+1.D0/2.D0*JG1

	  CMPL3X(I)=(PI/(SMALL+IMAGINARY*(T(J)-T(I))/2.D0))**1.5D0 &
	           & *DCONJG(D(PSTX-AX(J),PSTY-AY(J),PSTZ-AZ(J),ALFA))*CDEXP(-IMAGINARY*SST)* &
	           & (D(PSTX-AX(I),PSTY-AY(I),PSTZ-AZ(I),ALFA)*EX(I)+D(PSTY-AY(I),PSTX-AX(I),PSTZ-AZ(I),ALFA)*EY(I) &
			   & +D(PSTZ-AZ(I),PSTX-AX(I),PSTY-AY(I),ALFA)*EZ(I))*gt(i)
	  CMPL3Y(I)=(PI/(SMALL+IMAGINARY*(T(J)-T(I))/2.D0))**1.5D0 &
	           & *DCONJG(D(PSTY-AY(J),PSTX-AX(J),PSTZ-AZ(J),ALFA))*CDEXP(-IMAGINARY*SST)* &
	           & (D(PSTX-AX(I),PSTY-AY(I),PSTZ-AZ(I),ALFA)*EX(I)+D(PSTY-AY(I),PSTX-AX(I),PSTZ-AZ(I),ALFA)*EY(I) &
			   & +D(PSTZ-AZ(I),PSTX-AX(I),PSTY-AY(I),ALFA)*EZ(I))*gt(i)	   
 			   
	  CMPL3Z(I)=(PI/(SMALL+IMAGINARY*(T(J)-T(I))/2.D0))**1.5D0 &
 	           & *DCONJG(D(PSTZ-AZ(J),PSTX-AX(J),PSTY-AY(J),ALFA))*CDEXP(-IMAGINARY*SST)* &
	           & (D(PSTX-AX(I),PSTY-AY(I),PSTZ-AZ(I),ALFA)*EX(I)+D(PSTY-AY(I),PSTX-AX(I),PSTZ-AZ(I),ALFA)*EY(I) &
			   & +D(PSTZ-AZ(I),PSTX-AX(I),PSTY-AY(I),ALFA)*EZ(I))*gt(i)	   
			   	            
	  END DO
	  
	
     if (j.lt.2.d0*pi/omgax/ht) then
      CALL NEWCOT2_gfm(HT,CMPL3x,1,j-1,JG2)
	  else
      CALL NEWCOT2_gfm(HT,CMPL3x,j-1099,j-1,JG2)
	  end if
	  XT(J)=2.D0*DREAL(IMAGINARY*JG2)
	  
	  
!	  CALL NEWCOT2_gfm(HT,CMPL3Y,MAX(1,INT(J-1-2.D0*PI/OMGAY/HT*13.d0)),MAX(1,INT(j-1-2.D0*PI/OMGAY/HT*12.d0)),JG4)
	  CALL NEWCOT2_gfm(HT,CMPL3Y,1,j-1,JG4)
! 	  CALL NEWCOT2_gfm(HT,CMPL3Y,1,MIN(J-1,INT(2.D0*PI/OMGAY/HT)),JG4)
	  
	  YT(J)=2.D0*DREAL(IMAGINARY*JG4)
	  
	  CALL NEWCOT2_gfm(HT,CMPL3Z,1,j-1,JG6)

	  ZT(J)=2.D0*DREAL(IMAGINARY*JG6)


	  WRITE(12,'(1X,4(1X,F15.7))') T(J),XT(J),YT(J)
      


	  IF(MOD(J,200).EQ.0)WRITE(*,*) J

!      WRITE(*,*) J
	 END DO
!$OMP END PARALLEL DO



	
	OPEN(FILE='R(t).DAT',UNIT=10)
	

	DO J=1,N
	WRITE(10,'(1X,4(1X,F15.7))') T(J),XT(J),YT(J)
	

	END DO
	WRITE(*,*) 'IP=',IP
	WRITE(*,*) 'UP_X=',UPX
	WRITE(*,*) 'OMGAX=',OMGAX							
	WRITE(*,*) 'EX0=',EX0
	WRITE(*,*) 'FWHMX=',FWHMX
	WRITE(*,*) 'CUT OFF_X=',(3.17*UPX+IP)/OMGAX
	WRITE(*,*) 'GAMMAX=',GAMMAX
	WRITE(*,*) 'UP_Y=',UPY
	WRITE(*,*) 'OMGAY=',OMGAY							
	WRITE(*,*) 'EY0=',EY0
	WRITE(*,*) 'FWHMY=',FWHMY
	WRITE(*,*) 'CUT OFF_Y=',(3.17*UPY+IP)/OMGAY
    WRITE(*,*) 'GAMMAY=',GAMMAY
	WRITE(*,*) 'UP_Z=',UPZ
	WRITE(*,*) 'OMGAZ=',OMGAZ							
	WRITE(*,*) 'EZ0=',EZ0
	WRITE(*,*) 'FWHMZ=',FWHMZ
	WRITE(*,*) 'CUT OFF_Z=',(3.17*UPZ+IP)/OMGAZ
    WRITE(*,*) 'GAMMAZ=',GAMMAZ


!	 PAUSE
!	 STOP
	 END

