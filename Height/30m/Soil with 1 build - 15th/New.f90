!     Last change:  JP    5 Sep 2016    8:29 am
      PROGRAM JPN 
!     THE 1D SH-wave response of soil column above the basement using FD method
      implicit none
	  integer, PARAMETER :: IDXMX=5000,JDZMX=5000,FXI=250,FXF=2474,FZI=1281,FZF=2114,INTV=2
      REAL,ALLOCATABLE,DIMENSION (:) :: NPRNT,IRECX,IRECZ,IRECX1,IRECZ1
      real, ALLOCATABLE, DIMENSION (:) :: SU,SU1
      real, dimension (1:5000) :: VV
      real, dimension (1:5000,1:5000) :: TT,VI
      real, dimension (1:90000) ::addx,addy
      real, ALLOCATABLE, DIMENSION (:,:,:) :: riset
	  real, ALLOCATABLE, DIMENSION (:,:):: E,D,GX,GZ,R,RX,RY,AY1,AY2,AY3,AY4,XYY1,XYY2,XYY3,XYY4,ZYY1,ZYY2,ZYY3,ZYY4,slip,slip1,randis,rake,ptr
      real, ALLOCATABLE, DIMENSION (:,:,:) :: V,XY,ZY,GXY1,GXY2,GXY3,GXY4,GZY1,GZY2,GZY3,GZY4
	  real :: F1,T0,ALPHA,TAU,H,HP,HPP,G11,G12,G13,G14,EPSB,DTDH,PI,DM,RM,D1,R1,AY11,AY21,AY31,AY41,D2,R2,AY12,AY22,AY32,AY42,D3,R3,AY13,AY23,AY33,AY43,D4,R4,AY14,AY24,AY34,AY44,D5,R5,AY15,AY25,AY35,AY45,D6,R6,AY16,AY26,AY36,AY46,TMAX,DH,DH1,DT,DX,DX1,F0,TS,GM,PH,AF1,AF2,AF3,AF4,G21,G22,G23,G24,F11,F12,F13,F14,F21,F22,F23,F24,AJ,F,JB,AI,TR,M0,slip_sum,slip_max,slip_avg,num,lam,del,phi
      integer :: HI,HJ,IX,IXM,JZM,JZ,NT,NSRC1,NRECZ,NRECZ1,I,STBF,NREC,J,NT1,NT2,JZP,T,N,KZP,LA,NA,NTPM,JP1,JP2,JM1,JM2,IP,IP2,IM1,IM2,IR,IP1,c1,x,ptr_in,K,counter,ki,kj
      ALLOCATE (riset(FXI:FXF,FZI:FZF,0:700))
      OPEN(1,FILE='randis.dat')
      OPEN(2,FILE='slip.dat')
      OPEN(5,FILE='s.dot')
      OPEN(6,FILE='v.dat')
      OPEN(7,FILE='buildv.dat')
      OPEN(19,FILE='test.dat')
      OPEN(28,FILE='STF.dat')
      OPEN(38,FILE='lav.dat')
      EPSB=.1
      PI=3.1415926535897
      READ(5,*)IX,JZ,NRECZ,NRECZ1
      IXM=IX-1
      READ(5,*)DM,RM,D1,R1,AY11,AY21,AY31,AY41,D2,R2,AY12,AY22,AY32,AY42,D3,R3,AY13,AY23,AY33,AY43,D4,R4,AY14,AY24,AY34,AY44,TMAX,DH,DH1,DT,DX,DX1,F0,TS,GM,PH,AF1,AF2,AF3,AF4
      
	  ALLOCATE (ptr(IDXMX,JDZMX),slip(IDXMX,JDZMX),slip1(IDXMX,JDZMX),randis(IDXMX,JDZMX),rake(IDXMX,JDZMX),E(IDXMX,JDZMX),D(IDXMX,JDZMX),GX(IDXMX,JDZMX),GZ(IDXMX,JDZMX),R(IDXMX,JDZMX),RX(IDXMX,JDZMX),RY(IDXMX,JDZMX),AY1(IDXMX,JDZMX),AY2(IDXMX,JDZMX),AY3(IDXMX,JDZMX),AY4(IDXMX,JDZMX),XYY1(IDXMX,JDZMX),XYY2(IDXMX,JDZMX),XYY3(IDXMX,JDZMX),XYY4(IDXMX,JDZMX),ZYY1(IDXMX,JDZMX),ZYY2(IDXMX,JDZMX),ZYY3(IDXMX,JDZMX),ZYY4(IDXMX,JDZMX))
      ALLOCATE (V(IDXMX,JDZMX,2),XY(IDXMX,JDZMX,2),ZY(IDXMX,JDZMX,2),GXY1(IDXMX,JDZMX,2),GXY2(IDXMX,JDZMX,2),GXY3(IDXMX,JDZMX,2),GXY4(IDXMX,JDZMX,2),GZY1(IDXMX,JDZMX,2),GZY2(IDXMX,JDZMX,2),GZY3(IDXMX,JDZMX,2),GZY4(IDXMX,JDZMX,2))  
	  ALLOCATE (SU(IDXMX),SU1(IDXMX),NPRNT(IDXMX),IRECX(IDXMX),IRECZ(IDXMX),IRECX1(IDXMX),IRECZ1(IDXMX))
	  
      HI=1917
      HJ=1726
      lam=141
      del=87
      phi=140
      
      kj=0
      ki=0
      ptr_in=0
      slip_max=0
      slip_avg=0
      num=0
      slip_sum=0
	  
      do j=1,(((FXF-FXI)/INTV)+1)
        read(2,*) (slip(I,J),I=1,(((FZF-FZI)/INTV)+1))
        read(1,*) (randis(I,J),I=1,(((FZF-FZI)/INTV)+1))
      END DO 
      do i=FXI,FXF,INTV
      do j=FZI,FZF,INTV
        slip(i,j)=slip(j-FZI+1-kj,i-FXI+1-ki)
        randis(i,j)=randis(j-FZI+1-kj,i-FXI+1-ki)
        rake(i,j)=lam+2.5*(randis(j-FZI+1-kj,i-FXI+1-ki))
        riset(i,j,0)=randis(j-FZI+1-kj,i-FXI+1-ki)
        kj=kj+(INTV-1)
      end do
        ki=ki+(INTV-1)
        kj=0
      end do
     
      
       !M0=((112E16*dt)/(dh**3))
        M0=((794E15*dt)/(dh1**3))
      !M0=((16000000000*1000000*dt)/(dh**3)) 
       
      
      JZM=JZ-1
      DTDH=DT/DH
      DT=DH*DTDH
      NT=INT(TMAX/DT)
      STBF=SQRT(RM*DM)*DTDH
      IF (STBF>0.7044) then
      print *,'STABILITY VOILATED'
      STOP
      end if
	  ptr=0
      
	
      NREC=IX
      DO I=1,NREC
      IRECX(I)=I
      IRECZ(I)=NRECZ
      IRECZ1(I)=NRECZ1
	  END DO

      
	  DO  I=1,IX
      DO  J=1,30
         D(I,J)=D1
         R(I,J)=R1
         AY1(I,J)=AY11
         AY2(I,J)=AY21
         AY3(I,J)=AY31
         AY4(I,J)=AY41
      END DO
	  END DO     
	  DO  I=1,IX
      DO  J=31,59
         D(I,J)=D2
         R(I,J)=R2
         AY1(I,J)=AY12
         AY2(I,J)=AY22
         AY3(I,J)=AY32
         AY4(I,J)=AY42
      END DO
	  END DO
      DO  I=1,IX
      DO  J=60,JZ
         D(I,J)=D3
         R(I,J)=R3
         AY1(I,J)=AY13
         AY2(I,J)=AY23
         AY3(I,J)=AY33
         AY4(I,J)=AY43
      END DO
	  END DO
      !DO  K=2800,3525,25
	  DO  I=3150,(3150+10)
      DO  J=21,30
         D(I,J)=D4
         R(I,J)=R4
         AY1(I,J)=AY14
         AY2(I,J)=AY24
         AY3(I,J)=AY34
         AY4(I,J)=AY44
     END DO
     END DO
     !END DO
      
     DO  I=1,5000
     DO  J=1,1670
         VI(I,J)=0.56*SQRT(R(I,J)*D(I,J))
     END DO
     END DO
     !DO  I=1,IX
     !DO  J=253,313
     !    VI(I,J)=(0.56+(0.24/(314-J)))*SQRT(R(I,J)*D(I,J))
     !END DO
     !END DO
     DO  I=1,5000
     DO  J=1671,5000
         VI(I,J)=0.8*SQRT(R(I,J)*D(I,J))
     END DO
     END DO
     
     
     !DO  I=FXI,FXF
     !DO  J=FZI,FZI+222
     !    D(I,J)=0.7*D(I,J)         
     !END DO
     !END DO
     !DO  I=FXI,FXF
     !DO  J=FZI+223,FZF+55
     !    D(I,J)=(0.7+(0.3*0.0015*(J+1-FZI-223)))*D(I,J)         
     !END DO
     !END DO
     !
    

!     Gridding pattern
      DO  I=1,2775
      DO  J=1,1031
          GX(I,J)=DX1
          GZ(I,J)=DH
      END DO
      END DO
      DO  I=2775,IX
      DO  J=1,1031
          GX(I,J)=DX
          GZ(I,J)=DH
      END DO
      END DO
      DO  I=1,2775
      DO  J=1031,JZ
          GX(I,J)=DX1
          GZ(I,J)=DH1
      END DO
      END DO
      DO  I=2775,IX
      DO  J=1031,JZ
          GX(I,J)=DX
          GZ(I,J)=DH1
      END DO
      END DO
      
      DO  J=1,JZ
      DO  I=1,IX
         E(I,J)=SQRT(R(I,J)*D(I,J))*DT
      END DO
      END DO
!     AF IS ATTENUATING FREQUENCY
      G11=2*PI*AF1*DT/(2-2*PI*AF1*DT)
      G12=2*PI*AF2*DT/(2-2*PI*AF2*DT)
      G13=2*PI*AF3*DT/(2-2*PI*AF3*DT)
      G14=2*PI*AF4*DT/(2-2*PI*AF4*DT)

      G21=2/(2-2*PI*AF1*DT)
      G22=2/(2-2*PI*AF2*DT)
      G23=2/(2-2*PI*AF3*DT)
      G24=2/(2-2*PI*AF4*DT)

      DO J=1,JZM
      DO I=1,IXM
        XYY1(I,J)=2*G21*AY1(I,J)*R(I,J)*R(I+1,J)/(R(I,J)+R(I+1,J))
        XYY2(I,J)=2*G22*AY2(I,J)*R(I,J)*R(I+1,J)/(R(I,J)+R(I+1,J))
        XYY3(I,J)=2*G23*AY3(I,J)*R(I,J)*R(I+1,J)/(R(I,J)+R(I+1,J))
        XYY4(I,J)=2*G24*AY4(I,J)*R(I,J)*R(I+1,J)/(R(I,J)+R(I+1,J))

        ZYY1(I,J)=2*G21*AY1(I,J)*R(I,J)*R(I,J+1)/(R(I,J)+R(I,J+1))
        ZYY2(I,J)=2*G22*AY2(I,J)*R(I,J)*R(I,J+1)/(R(I,J)+R(I,J+1))
        ZYY3(I,J)=2*G23*AY3(I,J)*R(I,J)*R(I,J+1)/(R(I,J)+R(I,J+1))
        ZYY4(I,J)=2*G24*AY4(I,J)*R(I,J)*R(I,J+1)/(R(I,J)+R(I,J+1))
      END DO
      END DO
       NT1=1
       NT2=2
   
     CALL time (tt,HI,HJ)
   
     slip1=slip
     slip1=slip1-maxval(slip1)
     slip1=slip1*-1
     do i=FXI,FXF,INTV
     do j=FZI,FZI+277,INTV
         slip1(I,J)=2*slip1(i,j)         
     END DO
     END DO
     do i=FXI,FXF,INTV
     do j=FZI+278,FZI+444,INTV
         slip1(I,J)=(2-(0.005988*(J+1-FZI-278)))*slip1(i,j)         
     END DO
     END DO
     
     
     slip1=slip1/(sum(slip1)/size(slip1))
     print*,sum(slip1)/size(slip1)
     slip1=slip1*0.4
     print*,sum(slip1)/size(slip1)
     print*,maxval(slip1)
      do j=FZI,FZF,INTV
       WRITE(19,615) (slip1(I,J),I=FXI,FXF,INTV)
       615  FORMAT(1113F25.10)
      END DO  
     
      do i=FXI,FXF,INTV
      do j=FZI,FZF,INTV
        !CALL RICKMOD1 (DT,NSRC1,((slip1(i,j)/2)-0.0005),slip1(i,j),riset(i,j,:))
         ! if ((slip1(i,j))>0.048)then
           CALL LAH (DT,NSRC1,slip1(i,j),riset(i,j,:))
          !endif 
      end do
      end do
      riset(:,:,0)=0
   !      do i=FXI,FXF,INTV
   !   do j=FZI,FZF,INTV
   !         WRITE(28,150) (riset(I,J,K),K=0,250)
   !150      FORMAT(251F18.8)
   !   end do
   !   end do
        do i=FXI,FXF,INTV
        do j=FZI,FZF,INTV
            if ((TT(I,J)-(0.1364*randis(i,j)))>0 )then
            TT(I,J)=TT(I,J)-(0.1364*randis(i,j))
            end if
        end do 
        end do
        
        tt=int(tt/dt)
        tt(HI,HJ)=1  
        do j=FZI,FZF,INTV
         WRITE(14,606) (TT(I,J),I=FXI,FXF,INTV)
         606  FORMAT(1113F25.10)
        END DO 
       JZP=JZ
       t=-dt
  DO N=1,NT
      
      counter=0
      do i=FXI,FXF,INTV
      do j=FZI,FZF,INTV
        if(tt(i,j) == N) then
          ptr(i,j)=ptr(i,j)+1
          ptr_in=ptr_in+1
        end if
        if(riset(i,j,ptr(i,j))>0 ) then
         counter=counter+1
         addx(counter)=i
         addy(counter)=j         
        end if
      end do
      end do
      
      do i=1,counter 
       xy(addx(i),addy(i),nt2)=xy(addx(i),addy(i),nt2)-riset(addx(i),addy(i),ptr(addx(i),addy(i)))*M0*(sin(del*pi/180)*cos(rake(addx(i),addy(i))*pi/180)*cos(2*phi*pi/180)+0.5*sin(2*del*pi/180)*sin(rake(addx(i),addy(i))*pi/180)*sin(2*phi*pi/180))*slip(addx(i),addy(i))
       zy(addx(i),addy(i),nt2)=zy(addx(i),addy(i),nt2)-riset(addx(i),addy(i),ptr(addx(i),addy(i)))*M0*(cos(del*pi/180)*cos(rake(addx(i),addy(i))*pi/180)*cos(phi*pi/180)+cos(2*del*pi/180)*sin(rake(addx(i),addy(i))*pi/180)*sin(phi*pi/180))*slip(addx(i),addy(i))
!       WRITE(38,380),N,addx(i),addy(i),riset(addx(i),addy(i),ptr(addx(i),addy(i)))
!380 Format (I,3F25.10)  
       ptr(addx(i),addy(i))=ptr(addx(i),addy(i))+1
      end do
     	
     CALL CNTDF (NT1,NT2)
       t=t+dt
      DO I=1,IX
       LA=IRECX(I)
       NA=IRECZ(I)
       SU(I)=V(LA,NA,NT1)
       SU1(I)=V(LA,IRECZ1(I),NT1)
     END DO
      WRITE(6,604)(SU(i) ,i=2818,3543,25)
 604  FORMAT(30F18.8)

      WRITE(7,608)(SU1(i) ,i=2804,3529,25)
  608  FORMAT(30F18.8)

      NTPM=NT1
      NT1=NT2
      NT2=NTPM
    END DO
    print*,ptr_in

    contains   
    SUBROUTINE CNTDF (NT1,NT2)
       implicit none 
       integer, intent(in) :: NT1,NT2
      
      PI=3.1415926535897

      DO J=3,JZM-1
       JP1=J+1
       JP2=J+2
       JM1=J-1
       JM2=J-2
      DO I=3,IXM-1
       IP1=I+1
       IP2=I+2
       IM1=I-1
       IM2=I-2

      V(I,J,NT1)=V(I,J,NT2)+(DT*D(I,J))*((9/4*(XY(I,J,NT2)-XY(IM1,J,NT2))-1/12*(XY(IP1,J,NT2)-XY(IM2,J,NT2)))/(GX(I-1,J)+GX(I,J))+(9/4*(ZY(I,J,NT2)-ZY(I,JM1,NT2))-1/12*(ZY(I,JP1,NT2)-ZY(I,JM2,NT2)))/(GZ(I,J-1)+GZ(I,J)))

  END DO
  END DO
      CALL BNDZ (NT1,NT2,JZP)
      CALL BNDX (NT1,NT2)
       F11=4*PI*AF1*DT/(2+2*PI*AF1*DT)
       F12=4*PI*AF2*DT/(2+2*PI*AF2*DT)
       F13=4*PI*AF3*DT/(2+2*PI*AF3*DT)
       F14=4*PI*AF4*DT/(2+2*PI*AF4*DT)

       F21=(2-2*PI*AF1*DT)/(2+2*PI*AF1*DT)
       F22=(2-2*PI*AF2*DT)/(2+2*PI*AF2*DT)
       F23=(2-2*PI*AF3*DT)/(2+2*PI*AF3*DT)
       F24=(2-2*PI*AF4*DT)/(2+2*PI*AF4*DT)

      DO J=3,JZM-1
       JP1=J+1
       JP2=J+2
       JM1=J-1
       JM2=J-2
      DO I=3,IXM-1
       IP1=I+1
       IP2=I+2
       IM1=I-1
       IM2=I-2
      GXY1(I,J,NT1)=F11*(9/8*(V(IP1,J,NT1)-V(I,J,NT1))-1/24*(V(IP2,J,NT1)-V(IM1,J,NT1)))/GX(I,J)+F21*GXY1(I,J,NT2)

      GXY2(I,J,NT1)=F12*(9/8*(V(IP1,J,NT1)-V(I,J,NT1))-1/24*(V(IP2,J,NT1)-V(IM1,J,NT1)))/GX(I,J)+F22*GXY2(I,J,NT2)

      GXY3(I,J,NT1)=F13*(9/8*(V(IP1,J,NT1)-V(I,J,NT1))-1/24*(V(IP2,J,NT1)-V(IM1,J,NT1)))/GX(I,J)+F23*GXY3(I,J,NT2)

      GXY4(I,J,NT1)=F14*(9/8*(V(IP1,J,NT1)-V(I,J,NT1))-1/24*(V(IP2,J,NT1)-V(IM1,J,NT1)))/GX(I,J)+F24*GXY4(I,J,NT2)


      GZY1(I,J,NT1)=F11*(9/8*(V(I,JP1,NT1)-V(I,J,NT1))-1/24*(V(I,JP2,NT1)-V(I,JM1,NT1)))/GZ(I,J)+F21*GZY1(I,J,NT2)

      GZY2(I,J,NT1)=F12*(9/8*(V(I,JP1,NT1)-V(I,J,NT1))-1/24*(V(I,JP2,NT1)-V(I,JM1,NT1)))/GZ(I,J)+F22*GZY2(I,J,NT2)

      GZY3(I,J,NT1)=F13*(9/8*(V(I,JP1,NT1)-V(I,J,NT1))-1/24*(V(I,JP2,NT1)-V(I,JM1,NT1)))/GZ(I,J)+F23*GZY3(I,J,NT2)

      GZY4(I,J,NT1)=F14*(9/8*(V(I,JP1,NT1)-V(I,J,NT1))-1/24*(V(I,JP2,NT1)-V(I,JM1,NT1)))/GZ(I,J)+F24*GZY4(I,J,NT2)
  END DO
  END DO
    DO  J=1,JZ-1
    DO  I=1,IX-1
      RX(I,J)=(2*R(I,J)*R(I+1,J)/(R(I,J)+R(I+1,J)))*(1+G11*AY1(I,J)+G12*AY2(I,J)+G13*AY3(I,J)+G14*AY4(I,J))

      RY(I,J)=(2*R(I,J)*R(I,J+1)/(R(I,J)+R(I,J+1)))*(1+G11*AY1(I,J)+G12*AY2(I,J)+G13*AY3(I,J)+G14*AY4(I,J))
   END DO
   END DO
  DO J=3,JZM-1
      JP1=J+1
      JP2=J+2
      JM1=J-1
      JM2=J-2
  DO I=3,IXM-1
      IP1=I+1
      IP2=I+2
      IM1=I-1
      IM2=I-2
      XY(I,J,NT1)=XY(I,J,NT2)+DT*((RX(I,J)*(9/8*(V(IP1,J,NT1)-V(I,J,NT1))-1/24*(V(IP2,J,NT1)-V(IM1,J,NT1))))/GX(I,J)-(XYY1(I,J)*GXY1(I,J,NT1)+XYY2(I,J)*GXY2(I,J,NT1)+XYY3(I,J)*GXY3(I,J,NT1)+XYY4(I,J)*GXY4(I,J,NT1)))

      ZY(I,J,NT1)=ZY(I,J,NT2)+DT*((RY(I,J)*(9/8*(V(I,JP1,NT1)-V(I,J,NT1))-1/24*(V(I,JP2,NT1)-V(I,JM1,NT1))))/GZ(I,J)-(ZYY1(I,J)*GZY1(I,J,NT1)+ZYY2(I,J)*GZY2(I,J,NT1)+ZYY3(I,J)*GZY3(I,J,NT1)+ZYY4(I,J)*GZY4(I,J,NT1)))
   END DO
   END DO
      CALL BNDZ1 (NT1,NT2,JZP)
      CALL BNDX1 (NT1,NT2)
END SUBROUTINE CNTDF
    SUBROUTINE BNDZ (NT1,NT2,JZP)
      implicit none
      integer, intent(in) :: NT1,NT2,JZP
           
      DO I=2,IXM
         V(I,2,NT1)=V(I,3,NT1)+(V(I,4,NT1)-V(I,3,NT1))/10000
      END DO

      DO  I=2,IXM
           V(I,JZM,NT1)=V(I,JZM,NT2)-(E(I,JZM)/GZ(I,JZM))*(V(I,JZM,NT2)-V(I,JZM-1,NT2))
      END DO

      DO I=2,IXM
         V(I,JZ,NT1)=V(I,JZ,NT2)-(E(I,JZ)/GZ(I,JZM))*(V(I,JZ,NT2)-V(I,JZM,NT2))
      END DO

      DO J=2,199
       AJ=J
       F=((200.-AJ)/100)**3
       JB=JZ+1-J
      DO I=2,IXM
!      V(I,J,NT1)=V(I,J,NT1)-(EPSB*F)*(V(I,J,NT1)
!     $-V(I,J,NT2)-(E(I,J)/GZ(I,J))*
!     $(V(I,J+1,NT2)-V(I,J,NT2)))

         V(I,JB,NT1)=V(I,JB,NT1)-(EPSB*F)*(V(I,JB,NT1)-V(I,JB,NT2)+(E(I,JB)/GZ(I,JB))*(V(I,JB,NT2)-V(I,JB-1,NT2)))
      END DO
      END DO
  
      END SUBROUTINE BNDZ
	  
    SUBROUTINE BNDX (NT1,NT2)
      implicit none
      integer, intent(in) :: NT1,NT2
 
      

      DO J=2,JZ
         V(2,J,NT1)=V(2,J,NT2)+(E(2,J)/GX(2,J))*(V(3,J,NT2)-V(2,J,NT2))
      END DO

      DO J=2,JZ
         V(1,J,NT1)=V(1,J,NT2)+(E(1,J)/GX(1,J))*(V(2,J,NT2)-V(1,J,NT2))
      END DO

      DO J=1,JZ
         V(IXM,J,NT1)=V(IXM,J,NT2)-(E(IXM,J)/GX(IXM,J))*(V(IXM,J,NT2)-V(IXM-1,J,NT2))
      END DO

      DO J=1,JZ
         V(IX,J,NT1)=V(IX,J,NT2)-(E(IX,J)/GX(IX,J))*(V(IX,J,NT2)-V(IXM,J,NT2))
      END DO
      DO I=2,199
         AI=I
         F=((200.-AI)/100)**3
         IR=IX+1-I
      DO J=2,JZM
      V(I,J,NT1)=V(I,J,NT1)-(F*EPSB)*(V(I,J,NT1)-V(I,J,NT2)-(E(I,J)/GX(I,J))*(V(I+1,J,NT2)-V(I,J,NT2)))

      V(IR,J,NT1)=V(IR,J,NT1)-(F*EPSB)*(V(IR,J,NT1)-V(IR,J,NT2)+(E(IR,J)/GX(IR,J))*(V(IR,J,NT2)-V(IR-1,J,NT2)))
     END DO
     END DO
  
    END SUBROUTINE BNDX
    SUBROUTINE BNDZ1(NT1,NT2,JZP)
      implicit none
      integer, intent(in) :: NT1,NT2,JZP
 
      
      DO I=2,IXM
         ZY(I,1,NT1)=-ZY(I,4,NT1)/10000
         ZY(I,2,NT1)=-ZY(I,3,NT1)/10000
      END DO 

      DO I=2,IXM
       XY(I,JZM,NT1)=XY(I,JZM,NT2)-(E(I,JZM)/GZ(I,JZM))*(XY(I,JZM,NT2)-XY(I,JZM-1,NT2))

       ZY(I,JZM,NT1)=ZY(I,JZM,NT2)-(E(I,JZM)/GZ(I,JZM))*(ZY(I,JZM,NT2)-ZY(I,JZM-1,NT2))
      END DO

      DO  I=2,IXM
       XY(I,JZ,NT1)=XY(I,JZ,NT2)-(E(I,JZM)/GZ(I,JZ))*(XY(I,JZ,NT2)-XY(I,JZM,NT2))

       ZY(I,JZ,NT1)=ZY(I,JZ,NT2)-(E(I,JZM)/GZ(I,JZ))*(ZY(I,JZ,NT2)-ZY(I,JZM,NT2))
      END DO

      DO J=2,199
         AJ=J
         F=((200.-AJ)/100)**3
         JB=JZ+1-J
      DO I=2,IXM
!      xy(I,J,NT1)=xy(I,J,NT1)-(EPSB*F)*(xy(I,J,NT1)
!     $-xy(I,J,NT2)-(E(I,J)/GZ(I,J))*
!     $(xy(I,J+1,NT2)-xy(I,J,NT2)))

!      zy(I,J,NT1)=zy(I,J,NT1)-(EPSB*F)*(zy(I,J,NT1)
!     $-zy(I,J,NT2)-(E(I,J)/GZ(I,J))*
!     $(zy(I,J+1,NT2)-zy(I,J,NT2)))

        xy(I,JB,NT1)=xy(I,JB,NT1)-(EPSB*F)*(xy(I,JB,NT1)-xy(I,JB,NT2)+(E(I,JB)/GZ(I,JB))*(xy(I,JB,NT2)-xy(I,JB-1,NT2)))

        zy(I,JB,NT1)=zy(I,JB,NT1)-(EPSB*F)*(zy(I,JB,NT1)-zy(I,JB,NT2)+(E(I,JB)/GZ(I,JB))*(zy(I,JB,NT2)-zy(I,JB-1,NT2)))
      END DO
      END DO
      
      END SUBROUTINE BNDZ1
    SUBROUTINE BNDX1 (NT1,NT2)
     implicit none
     integer, intent(in) :: NT1,NT2
  
      
      DO J=1,JZ
      XY(I,2,NT1)=XY(I,2,NT2)+(E(I,2)/GX(I,2))*(XY(I,3,NT2)-XY(I,2,NT2))

      ZY(I,2,NT1)=ZY(I,2,NT2)+(E(I,2)/GZ(I,2))*(ZY(I,3,NT2)-ZY(I,2,NT2))
      END DO

      DO J=1,JZ
      XY(I,1,NT1)=XY(I,1,NT2)+(E(I,1)/GX(I,1))*(XY(I,2,NT2)-XY(I,1,NT2))

      ZY(I,1,NT1)=ZY(I,1,NT2)+(E(I,1)/GZ(I,1))*(ZY(I,2,NT2)-ZY(I,1,NT2))
      END DO

      DO J=1,JZ
      XY(IXM,J,NT1)=XY(IXM,J,NT2)-(E(IXM,J)/GX(IXM,J))*(XY(IXM,J,NT2)-XY(IXM-1,J,NT2))

      ZY(IXM,J,NT1)=ZY(IXM,J,NT2)-(E(IXM,J)/GX(IXM,J))*(ZY(IXM,J,NT2)-ZY(IXM-1,J,NT2))
      END DO

      DO J=1,JZ
      XY(IX,J,NT1)=XY(IX,J,NT2)-(E(IX,J)/GX(IXM,J))*(XY(IX,J,NT2)-XY(IXM,J,NT2))

      ZY(IX,J,NT1)=ZY(IX,J,NT2)-(E(IX,J)/GX(IXM,J))*(ZY(IX,J,NT2)-ZY(IXM,J,NT2))
      END DO

      DO I=2,199
       AI=I
       F=((200.-AI)/100)**3
       IR=IX+1-I
      DO J=2,JZM
      xy(I,J,NT1)=xy(I,J,NT1)-(F*EPSB)*(xy(I,J,NT1)-xy(I,J,NT2)-(E(I,J)/GX(I,J))*(xy(I+1,J,NT2)-xy(I,J,NT2)))
      zy(I,J,NT1)=zy(I,J,NT1)-(F*EPSB)*(zy(I,J,NT1)-zy(I,J,NT2)-(E(I,J)/GX(I,J))*(zy(I+1,J,NT2)-zy(I,J,NT2)))

      xy(IR,J,NT1)=xy(IR,J,NT1)-(F*EPSB)*(xy(IR,J,NT1)-xy(IR,J,NT2)+(E(IR,J)/GX(IR,J))*(xy(IR,J,NT2)-xy(IR-1,J,NT2)))
      zy(IR,J,NT1)=zy(IR,J,NT1)-(F*EPSB)*(zy(IR,J,NT1)-zy(IR,J,NT2)+(E(IR,J)/GX(IR,J))*(zy(IR,J,NT2)-zy(IR-1,J,NT2)))
      END DO
      END DO
      END SUBROUTINE BNDX1
 
 
 SUBROUTINE RICKMOD1 (DT,LS,TS,TR,F)
   implicit none
   INTEGER, INTENT(OUT) :: LS
   REAL, INTENT (IN) :: DT,TS,TR
   real, dimension (0:500), INTENT (OUT) :: F
   real :: t,  S, K, C1, C2, C3, C4, C5, C6, pi,nt,i,ut,tot
   
   print*,ts,tr
    ! Body of Console16
     
    !Ts=0.2
    ut=1
    pi=3.1415926
    !Tr=0.5
    
    LS=(Tr+(2*Ts))/dt
     K=(2/(pi*Tr*Ts**2))*ut
     
     OPEN(8,FILE='STF.dat')
     do i=0,LS
        t=i*dt   
        C1=((0.5*t+0.25*Tr)*((t*(Tr-t))**0.5))+((t*Tr-Tr**2)*(asin((t/Tr)**0.5)))-((0.75*Tr**2)*(atan(((Tr-t)/t)**0.5)))
        C2=(0.375*pi*Tr**2)
        C3=(Ts-t-0.5*Tr)*(((t-Ts)*(Tr-t+Ts))**0.5)+(Tr*(2*Tr-2*t+2*Ts))*(asin(((t-Ts)/Tr)**0.5))+((1.5)*(Tr**2)*atan(((Tr-t+Ts)/(t-Ts))**0.5))
        C4=((-Ts+(0.5*t)+(0.25*Tr))*(((t-2*Ts)*(Tr-t+2*Ts))**0.5))+(Tr*(-Tr+t-2*Ts)*asin(((t-2*Ts)/Tr)**0.5))-(0.75*(Tr**2)*atan(((Tr-t+2*Ts)/(t-2*Ts))**0.5))
        C5=(pi/2)*(Tr)*(t-Tr)
        C6=(pi/2)*(Tr)*(2*Ts-t+Tr)
       
       
            if (Tr>2*Ts) then
              
                if (t<0) then
                     S=0
                else if (0<=t.AND.t<Ts) then
                    S=K*(C1+C2)
                else if (Ts<=t.AND.t<(2*Ts)) then
                    S=K*(C1-C2+C3)
                else if ((2*Ts)<=t.AND.t<Tr) then
                    S=K*(C1+C3+C4)
                else if (Tr<=t.AND.t<(Tr+Ts)) then
                    S=K*(C5+C3+C4)
                else if ((Tr+Ts)<=t.AND.t<(Tr+2*Ts)) then
                    S=K*(C4+C6)
                else if ((Tr+2*Ts)<=t) then
                    S=0
                end if
            else if (Ts<Tr.AND.Tr<(2*Ts)) then
                if (t<0) then
                    S=0
                else if (0<=t.AND.t<Ts) then
                    S=K*(C1+C2)
                else if (Ts<=t.AND.t<(Tr)) then
                    S=K*(C1-C2+C3)
                else if (Tr<=t.AND.t<(2*Ts)) then
                    S=K*(C5+C3-C2)
                else if ((2*Ts)<=t.AND.(t<Tr+Ts)) then
                    S=K*(C5+C3+C4)
                else if ((Tr+Ts)<=t.AND.t<(Tr+(2*Ts))) then
                    S=K*(C4+C6)
                else if ((Tr+2*Ts)<=t) then
                    S=0
                end if
            end if
            
   !          WRITE(5,150) t,s
   !150       FORMAT(2F18.8)
            tot=tot+(s*dt)
            f(i)=s
        end do
		do i=0,LS
            f(i)=f(i)
            WRITE(8,150) F(I)
   150      FORMAT(F18.8)
		end do
       !print *,tot
      END SUBROUTINE RICKMOD1
	 
	!  SUBROUTINE RICKMOD (DT,F0,LS,TS,F)
 !     implicit none
 !     INTEGER, INTENT(OUT) :: LS
 !     REAL, INTENT (IN) :: DT,F0,TS
 !     real, dimension (1:5000), INTENT (OUT) :: F
 !     REAL :: T,G,GP,GPP,FMAX  
 !     F1=.87*F0
 !     T0=1./F1
 !     LS=.001+(T0/DT)
 !     T0=LS*DT
 !     LS=2*LS+1
 !     T0=TS
 !     FMAX=0
 !     ALPHA=(PI*F0)**2
 !     DO I=1,LS
 !     T=DT*(I-1)
 !     TAU=T/T0
 !
 !     G=(1.0-(TAU-1)**2)**3
 !     GP=-6.0*(TAU-1)*(1.0-(TAU-1)**2)**2/T0
 !     GPP=6.0*(1.0-(TAU-1)**2)*(5.0*(TAU-1)**2-1.0)/(T0*T0)
 !
 !     H=-EXP(-ALPHA*(T-T0)**2)
 !     HP=-2*ALPHA*(T-T0)*H
 !     HPP=-2*ALPHA*(1-2*ALPHA*(T-T0)**2)*H
 !     F(I)=GPP*H+2.0*GP*HP+G*HPP
 !     IF(ABS(F(I))>FMAX) then 
 !     FMAX=ABS(F(I))
 !     END IF
 !End do
 !     FMAX=1.0E+10/FMAX
 !     DO I=1,LS
 !     F(I)= F(I)*FMAX
 !WRITE(9,608)F(I)
 !    608  FORMAT(F15.2)
 !    
 !     end do
        
      !END SUBROUTINE RICKMOD
   SUBROUTINE RICKMOD (DT,F0,LS,TS,GM,PH,F)
   implicit none
   INTEGER, INTENT(OUT) :: LS
   REAL, INTENT (IN) :: DT,F0,TS,GM,PH
   real, dimension (1:5000), INTENT (OUT) :: F
   REAL :: T,G,GP,GPP,TOT   
     !! DATA PI/3.1415927/
      LS=INT(2*TS/DT)      !lS is number of iteration for gabor
      !WRITE(6,*)DT,F0,LS,TS,GM,PH
      DO I=0,LS
         T=DT*(I-1)
         G=(2*PI*F0*(T-TS)/GM)**2
         GP=EXP(-G)
         GPP=COS(2*PI*F0*(T-TS)+PH)
         F(I)=GP*GPP
         TOT=TOT+F(I)
!      WRITE(6,101) F(I)
!  101 FORMAT (F10.6)
      END DO
       OPEN(8,FILE='STF.dat')
        DO I=0,LS
         F(I)= F(I)/TOT
         WRITE(8,150) F(I)
   150      FORMAT(F18.8)
      END DO
      
 END SUBROUTINE RICKMOD 
SUBROUTINE LAH (DT,LS,TR,F)
   implicit none
   INTEGER, INTENT(OUT) :: LS
   REAL, INTENT (IN) :: DT,TR
   real, dimension (0:700), INTENT (OUT) :: F
   real :: t,  S, K, CN, pi,nt,i,ut,tot,t1,t2
   
    !tr=0.35
    t1=0.13*tr
    t2=tr-t1
    pi=3.1415926
    LS=(tr)/dt
   
     
     do i=0,LS
        t=i*dt   
        CN=pi/((1.4*pi*t1)+(1.2*t1)+(0.3*pi*t2))
       
            if (0<=t.AND.t<t1) then
              S=CN*(0.7-0.7*cos(pi*t/t1)+0.6*sin(0.5*pi*t/t1))
            
            else if (t1<=t.AND.T<(2*T1)) then
              S=CN*(1.0-0.7*cos(pi*t/t1)+0.3*cos(pi*(t-t1)/t2))
            
            else if (2*t1<=t.AND.T<tr) then
              S=CN*(0.3+0.3*cos(pi*(t-t1)/t2))
            end if
            
   !          WRITE(5,150) t,s
   !150       FORMAT(2F18.8)
            tot=tot+(s*dt)
            f(i)=s
        end do
		   !print *,tot
      !END SUBROUTINE RICKMOD1
 END SUBROUTINE LAH
 SUBROUTINE time(t,ci,cj)
    implicit none
    integer :: ri,i,j,h
    integer,INTENT (IN) :: CI,CJ
    real, dimension (1:5000,1:5000), INTENT (OUT) :: t
    
    h=18
    
    OPEN(14,FILE='atime.dat')
    
    t(CI,CJ)=0    
   
    
    
    t(ci+1,cj)=h/vi(ci+1,cj)
    t(ci-1,cj)=h/vi(ci-1,cj)
    t(ci,cj+1)=h/vi(ci,cj+1)
    t(ci,cj-1)=h/vi(ci,cj-1)
    
    t(ci+1,cj-1)=(2*(h/vi(ci+1,cj-1))**2-(t(ci+1,cj)-t(ci,cj-1))**2)**0.5   
    t(ci-1,cj-1)=(2*(h/vi(ci-1,cj-1))**2-(t(ci-1,cj)-t(ci,cj-1)))**0.5   
    t(ci+1,cj+1)=(2*(h/vi(ci+1,cj+1))**2-(t(ci+1,cj)-t(ci,cj+1)))**0.5   
    t(ci-1,cj+1)=(2*(h/vi(ci-1,cj+1))**2-(t(ci-1,cj)-t(ci,cj+1)))**0.5   
    
    do ri=3,1920
     
	   if ((cj-ri+1)>0) then 
           t((ci),(cj-ri+1))=t((ci),(cj-ri+2))+((h/vi((ci),(cj-ri+1)))**2-0.25*(t(ci+1,cj-ri+2)-t(ci-1,cj-ri+2))**2)**0.5 
       endif
     
	     t((ci),(cj+ri-1))=t((ci),(cj+ri-2))+((h/vi((ci),(cj+ri-1)))**2-0.25*(t(ci+1,cj+ri-2)-t(ci-1,cj+ri-2))**2)**0.5
	 
         t((ci+ri-1),(cj))=t((ci+ri-2),(cj))+((h/vi((ci+ri-1),(cj)))**2-0.25*(t(ci+ri-2,cj+1)-t(ci+ri-2,cj-1))**2)**0.5
	    
         t((ci-ri+1),(cj))=t((ci-ri+2),(cj))+((h/vi((ci-ri+1),(cj)))**2-0.25*(t(ci-ri+2,cj+1)-t(ci-ri+2,cj-1))**2)**0.5
	    
      if ((cj-ri+1)>0) then 
       if (((h/vi((ci),(cj-ri+1)))**2-0.25*(t(ci+1,cj-ri+2)-t(ci-1,cj-ri+2))**2) < 0) then
	       t((ci),(cj-ri+1))=t((ci),(cj-ri+2))+((h/vi((ci),(cj-ri+1))))   
       end if
      endif
      
	  if (((h/vi((ci),(cj+ri-1)))**2-0.25*(t(ci+1,cj+ri-2)-t(ci-1,cj+ri-2))**2) < 0) then
	       t((ci),(cj+ri-1))=t((ci),(cj+ri-2))+((h/vi((ci),(cj+ri-1))))
	  end if
	  if (((h/vi((ci+ri-1),(cj)))**2-0.25*(t(ci+ri-2,cj+1)-t(ci+ri-2,cj-1))**2) < 0) then
	      t((ci+ri-1),(cj))=t((ci+ri-2),(cj))+((h/vi((ci+ri-1),(cj))))
      end if
	  if (((h/vi((ci-ri+1),(cj)))**2-0.25*(t(ci-ri+2,cj+1)-t(ci-ri+2,cj-1))**2) < 0) then
	     t((ci-ri+1),(cj))=t((ci-ri+2),(cj))+((h/vi((ci-ri+1),(cj))))
	  end if
      do i=(ci+1),(ci+(ri-2))
         if ((cj-ri+1)>0) then 
          t(i,(cj-ri+1))=t(i-1,cj-ri+2)+(2*(h/vi(i,(cj-ri+1)))**2-(t(i-1,cj-ri+1)-t(i,cj-ri+2))**2)**0.5   
         end if 
		   t(i,(cj+ri-1))=t(i-1,cj+ri-2)+(2*(h/vi(i,(cj+ri-1)))**2-(t(i-1,cj+ri-1)-t(i,cj+ri-2))**2)**0.5
          if ((cj-ri+1)>0) then  
		  if ((2*(h/vi(i,(cj-ri+1)))**2-(t(i-1,cj-ri+1)-t(i,cj-ri+2))**2) < 0) then
		    t(i,(cj-ri+1))=t(i-1,cj-ri+2)+((h/vi(i,(cj-ri+1))))
          end if
          end if
		  if ((2*(h/vi(i,(cj+ri-1)))**2-(t(i-1,cj+ri-1)-t(i,cj+ri-2))**2) < 0) then 
		    t(i,(cj+ri-1))=t(i-1,cj+ri-2)+((h/vi(i,(cj+ri-1)))) 
		  end if
      end do
      do i=(ci-1),(ci-(ri-2)),-1
          if ((cj-ri+1)>0) then 
           t(i,(cj-ri+1))=t(i+1,cj-ri+2)+(2*(h/vi(i,(cj-ri+1)))**2-(t(i+1,cj-ri+1)-t(i,cj-ri+2))**2)**0.5   
		   end if
           t(i,(cj+ri-1))=t(i+1,cj+ri-2)+(2*(h/vi(i,(cj+ri-1)))**2-(t(i+1,cj+ri-1)-t(i,cj+ri-2))**2)**0.5
           if ((cj-ri+1)>0) then 
		    if ((2*(h/vi(i,(cj-ri+1)))**2-(t(i+1,cj-ri+1)-t(i,cj-ri+2))**2) < 0) then
		    t(i,(cj-ri+1))=t(i+1,cj-ri+2)+((h/vi(i,(cj-ri+1))))
    end if
    end if
		  if ((2*(h/vi(i,(cj+ri-1)))**2-(t(i+1,cj+ri-1)-t(i,cj+ri-2))**2) < 0) then 
		   t(i,(cj+ri-1))=t(i+1,cj+ri-2)+((h/vi(i,(cj+ri-1)))) 
		  end if
      end do
      
	  do j=(cj+1),(cj+(ri-2))
           t((ci-ri+1),j)=t(ci+ri-2,j-1)+(2*(h/vi((ci-ri+1),j))**2-(t(ci+ri-1,j-1)-t(ci+ri-2,j))**2)**0.5   
		   t((ci+ri-1),j)=t(ci-ri+2,j-1)+(2*(h/vi((ci+ri-1),j))**2-(t(ci-ri+1,j-1)-t(ci-ri+2,j))**2)**0.5   
           if ((2*(h/vi((ci-ri+1),j))**2-(t(ci+ri-1,j-1)-t(ci+ri-2,j))**2) < 0) then
		    t((ci-ri+1),j)=t(ci+ri-2,j-1)+((h/vi((ci-ri+1),j)))
		  end if
		  if ((2*(h/vi((ci+ri-1),j))**2-(t(ci-ri+1,j-1)-t(ci-ri+2,j))**2) < 0) then 
		   t((ci+ri-1),j)=t(ci-ri+2,j-1)+((h/vi((ci+ri-1),j))) 
		  end if
      end do
      do j=(cj-1),(cj-(ri-2)),-1
          if(j>0) then
           t((ci-ri+1),j)=t(ci+ri-2,j+1)+(2*(h/vi((ci-ri+1),j))**2-(t(ci+ri-1,j+1)-t(ci+ri-2,j))**2)**0.5   
	       t((ci+ri-1),j)=t(ci-ri+2,j+1)+(2*(h/vi((ci+ri-1),j))**2-(t(ci-ri+1,j+1)-t(ci-ri+2,j))**2)**0.5   
           if ((2*(h/vi((ci-ri+1),j))**2-(t(ci+ri-1,j+1)-t(ci+ri-2,j))**2) < 0) then
		    t((ci-ri+1),j)=t(ci+ri-2,j+1)+((h/vi((ci-ri+1),j)))
		  end if
		  if ((2*(h/vi((ci+ri-1),j))**2-(t(ci-ri+1,j+1)-t(ci-ri+2,j))**2) < 0) then 
		   t((ci+ri-1),j)=t(ci-ri+2,j+1)+((h/vi((ci+ri-1),j))) 
          end if
          end if
      end do
       
        if ((cj-ri+1)>0) then 
		t((ci+ri-1),(cj-ri+1))=t((ci+ri-2),(cj-ri+2))+(2*(h/vi((ci+ri-1),(cj-ri+1)))**2-(t(ci+ri-2,cj-ri+1)-t(ci+ri-1,cj-ri+2))**2)**0.5   
        end if
        t((ci+ri-1),(cj+ri-1))=t((ci+ri-2),(cj+ri-2))+(2*(h/vi((ci+ri-1),(cj+ri-1)))**2-(t(ci+ri-2,cj+ri-1)-t(ci+ri-1,cj+ri-2))**2)**0.5
		t((ci-ri+1),(cj+ri-1))=t((ci-ri+2),(cj+ri-2))+(2*(h/vi((ci-ri+1),(cj+ri-1)))**2-(t(ci-ri+2,cj+ri-1)-t(ci-ri+1,cj+ri-2))**2)**0.5
		if ((cj-ri+1)>0) then 
        t((ci-ri+1),(cj-ri+1))=t((ci-ri+2),(cj-ri+2))+(2*(h/vi((ci-ri+1),(cj-ri+1)))**2-(t(ci-ri+2,cj-ri+1)-t(ci-ri+1,cj-ri+2))**2)**0.5   
        end if
        if ((cj-ri+1)>0) then 
        if ((2*(h/vi((ci+ri-1),(cj-ri+1)))**2-(t(ci+ri-2,cj-ri+1)-t(ci+ri-1,cj-ri+2))**2) < 0) then 
		   t((ci+ri-1),(cj-ri+1))=t((ci+ri-2),(cj-ri+2))+((h/vi((ci+ri-1),(cj-ri+1)))) 
        end if
        end if
        if ((2*(h/vi((ci+ri-1),(cj+ri-1)))**2-(t(ci+ri-2,cj+ri-1)-t(ci+ri-1,cj+ri-2))**2) < 0) then 
		   t((ci+ri-1),(cj+ri-1))=t((ci+ri-2),(cj+ri-2))+((h/vi((ci+ri-1),(cj+ri-1)))) 
        end if
        if ((2*(h/vi((ci-ri+1),(cj+ri-1)))**2-(t(ci-ri+2,cj+ri-1)-t(ci-ri+1,cj+ri-2))**2) < 0) then 
		   t((ci-ri+1),(cj+ri-1))=t((ci-ri+2),(cj+ri-2))+((h/vi((ci-ri+1),(cj+ri-1)))) 
        end if
        if ((cj-ri+1)>0) then 
        if ((2*(h/vi((ci-ri+1),(cj-ri+1)))**2-(t(ci-ri+2,cj-ri+1)-t(ci-ri+1,cj-ri+2))**2) < 0) then 
		   t((ci-ri+1),(cj-ri+1))=t((ci-ri+2),(cj-ri+2))+((h/vi((ci-ri+1),(cj-ri+1)))) 
        end if
        end if
	end do 
    
    
    END SUBROUTINE time
	  END PROGRAM JPN 