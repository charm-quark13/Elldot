      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER NP,NP2,I,J,L,LN,CORR,ITER,PC,PC_COUNT,LWORK,INFO
      INCLUDE 'dim.inc'   
      PARAMETER (LWORK = 100, NP2 = 2*NP)
      INTEGER IPIV(NP2),KLI

      DOUBLE PRECISION T,TOL,EOLD,U0,U1,MIX,EX,CRIT,
     &                 EH,EVXC,EXC,ETOT,TX(NP),TY(NP),TZ(NP),TT(3),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXD(NP),
     &                 BXCZ(NP),VHXC0(NP),BXCX0(NP),BXCY0(NP),BXCZ0(NP),
     &                 VHXC1(NP),BXCX1(NP),BXCY1(NP),BXCZ1(NP),
     &                 V(NP),BX(NP),BY(NP),BZ(NP),MTOT(3),
     &                 VT(NP),BTX(NP),BTY(NP),BTZ(NP),MM(NP),
     &                 N(NP),MX(NP),MY(NP),MZ(NP),E(NP2),RWORK(100)
      DOUBLE PRECISION TIME,DT,OMEGA

      DOUBLE COMPLEX M(NP2,NP2),WORK(LWORK),ZERO,ONE,IONE,
     &               A(NP2,NP2),R(NP2),PSI1(NP2),PSI2(NP2),
     &               PSI1P(NP2),PSI2P(NP2)
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      T = 0.5D0
      WRITE(*,*)'U0?'
      READ(*,*)U0
      U1 = U0/2.D0

C      WRITE(*,*)'Number of PC steps?'
C      READ(*,*)PC
      PC=2

      WRITE(*,*)'Slater (0) or KLI (1)?'
      READ(*,*)KLI

      WRITE(*,*)'XC torque? (1 yes, 2 no)'
      READ(*,*)L
      IF (L.EQ.1) THEN
         WRITE(*,*)'Torque correction? (0=no, 1=yes)'
         READ(*,*)CORR
      ENDIF

      LN = 20
      IF (KLI.EQ.1) LN = 50
      IF (CORR.EQ.1) LN = LN+10
      IF (L.EQ.2) LN = LN+20

      OMEGA = 0.5D0
      DT = 0.01D0
C**---------------------------------------------------------------------
C**   First, calculate the ground state
C**---------------------------------------------------------------------
      MIX = 0.20D0
      TOL = 1.D-10
      EOLD = 0.D0
      ITER = 0
      
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO    

      DO I=1,NP
         VHXC(I) = 0.D0
         BXCX(I) = 0.D0
         BXCY(I) = 0.D0
         BXCZ(I) = 0.D0
      ENDDO    

1     CONTINUE
      ITER = ITER + 1
      WRITE(*,*)
      WRITE(*,*)'************',ITER,'*************'
      WRITE(*,*)
      IF (ITER.GT.100000) STOP

      DO I=1,NP
         VT(I) = V(I) + VHXC(I)
         BTX(I) = BX(I) + BXCX(I)
         BTY(I) = BY(I) + BXCY(I)
         BTZ(I) = BZ(I) + BXCZ(I)
      ENDDO   

      CALL MATRIX(M,T,VT,BTX,BTY,BTZ)
      
      CALL ZHEEV( 'V', 'U', NP2, M, NP2, E, WORK, LWORK, RWORK, INFO )

      DO I=1,NP2
         PSI1(I) = M(I,1)
         PSI2(I) = M(I,2)
         WRITE(*,*)E(I)
      ENDDO
C**----------------------------------------------------------------------
C**   calculate the densities and update the potentials
C**----------------------------------------------------------------------
      CALL DENCALC(PSI1,PSI2,N,MX,MY,MZ,MM,MTOT)

      DO I=1,NP
         VHXC0(I) = VHXC(I)
         BXCX0(I) = BXCX(I)
         BXCY0(I) = BXCY(I)
         BXCZ0(I) = BXCZ(I)
      ENDDO                

      CALL XCPOT_SLATER(KLI,PSI1,PSI2,U0,U1,VHXC,BXCX,BXCY,BXCZ)

      IF (L.EQ.2) THEN
      DO I=1,NP
         BXCX1(I) = BXCX(I)
         BXCY1(I) = BXCY(I)
         BXCZ1(I) = BXCZ(I)
      ENDDO                

C**   project out transverse BXC

      DO I=1,NP
         BXCX(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MX(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCY(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MY(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCZ(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MZ(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
      ENDDO    

      ENDIF

      IF (CORR.EQ.1) THEN
         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
         CALL BCORR(MX,MY,MZ,BXCX,BXCY,BXCZ,TT)
      ENDIF

      DO 21 I=1,NP
         VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXC0(I)
         BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCX0(I)
         BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCY0(I)
         BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZ0(I)
21    CONTINUE

      CRIT = 0.D0
      DO I=1,NP
         CRIT = CRIT + DSQRT(MX(I)**2 + MY(I)**2 + MZ(I)**2)         
      ENDDO

      IF (DABS((CRIT- EOLD)/CRIT).GT.TOL) THEN
          EOLD = CRIT
          GOTO 1
      ENDIF    

      CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)

      TIME = 0.D0
      WRITE(LN,999)TIME,N(1),N(2),N(3),N(4)
      WRITE(LN+1,999)TIME,MX(1),MX(2),MX(3),MX(4)
      WRITE(LN+2,999)TIME,MY(1),MY(2),MY(3),MY(4)
      WRITE(LN+3,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      WRITE(LN+4,999)TIME,MM(1),MM(2),MM(3),MM(4)
      WRITE(LN+5,999)TIME,MTOT(1),MTOT(2),MTOT(3),
     &             DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
      IF (L.EQ.1) THEN
      WRITE(LN+6,999)TIME,TX(1),TX(2),TX(3),TX(4)
      WRITE(LN+7,999)TIME,TY(1),TY(2),TY(3),TY(4)
      WRITE(LN+8,999)TIME,TZ(1),TZ(2),TZ(3),TZ(4)
      WRITE(LN+9,997)TIME,TT(1),TT(2),TT(3)
      ENDIF
C**---------------------------------------------------------------------
C**   Now start the time propagation with a sudden switch
C**   (read in the new potential and magnetic field)
C**---------------------------------------------------------------------
      READ(1,*)
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO    

100   CONTINUE
      TIME = TIME + DT
C      WRITE(*,*)TIME

C      DO 16 I=1,NP
C         BXD(I) = BX(I) + 0.01D0*DSIN(OMEGA*(TIME-DT/2.D0))
C16    CONTINUE

      DO I=1,NP
         VHXC0(I) = VHXC(I)
         BXCX0(I) = BXCX(I)
         BXCY0(I) = BXCY(I)
         BXCZ0(I) = BXCZ(I)
      ENDDO    
      IF (CORR.EQ.1) THEN
         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX0,BXCY0,BXCZ0)
         CALL BCORR(MX,MY,MZ,BXCX0,BXCY0,BXCZ0,TT)
      ENDIF

      PC_COUNT=-1
120   CONTINUE
      PC_COUNT = PC_COUNT+1

      DO I=1,NP
            VT(I) = V(I) + VHXC(I)
            BTX(I) = BX(I) + BXCX(I)
C            BTX(I) = BXD(I) + BXCX(I)
            BTY(I) = BY(I) + BXCY(I)
            BTZ(I) = BZ(I) + BXCZ(I)
      ENDDO    

      CALL MATRIX(M,T,VT,BTX,BTY,BTZ)

      DO I=1,NP2
      DO J=1,NP2
         A(I,J) = -0.5D0*IONE*DT*M(I,J)
         IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
      ENDDO    

      DO I=1,NP2
         R(I) = ZERO
         DO J=1,NP2
            R(I) = R(I) + A(I,J)*PSI1(J)
         ENDDO
      ENDDO    

      DO I=1,NP2
      DO J=1,NP2
         A(I,J) = 0.5D0*IONE*DT*M(I,J)
         IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
      ENDDO    

      CALL ZGESV( NP2, 1, A, NP2, IPIV, R, NP2, INFO )
      DO I=1,NP2
         PSI1P(I) = R(I)
      ENDDO

      DO I=1,NP2
      DO J=1,NP2
         A(I,J) = -0.5D0*IONE*DT*M(I,J)
         IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
      ENDDO    

      DO I=1,NP2
         R(I) = ZERO
         DO J=1,NP2
            R(I) = R(I) + A(I,J)*PSI2(J)
         ENDDO
      ENDDO    

      DO I=1,NP2
      DO J=1,NP2
         A(I,J) = 0.5D0*IONE*DT*M(I,J)
         IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
      ENDDO    

      CALL ZGESV( NP2, 1, A, NP2, IPIV, R, NP2, INFO )
      DO I=1,NP2
         PSI2P(I) = R(I)
      ENDDO

      CALL XCPOT_SLATER(KLI,PSI1P,PSI2P,U0,U1,VHXC1,BXCX1,BXCY1,BXCZ1)
      IF (CORR.EQ.1) THEN
         CALL DENCALC(PSI1P,PSI2P,N,MX,MY,MZ,MM,MTOT)
         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX1,BXCY1,BXCZ1)
         CALL BCORR(MX,MY,MZ,BXCX1,BXCY1,BXCZ1,TT)
      ENDIF

      IF (PC_COUNT.LT.PC) THEN

         DO I=1,NP
            VHXC(I) = 0.5D0*(VHXC0(I) + VHXC1(I))
            BXCX(I) = 0.5D0*(BXCX0(I) + BXCX1(I))
            BXCY(I) = 0.5D0*(BXCY0(I) + BXCY1(I))
            BXCZ(I) = 0.5D0*(BXCZ0(I) + BXCZ1(I))
         ENDDO    
         GOTO 120

      ELSE

         DO I=1,NP2
            PSI1(I) = PSI1P(I)
            PSI2(I) = PSI2P(I)
         ENDDO    
         DO I=1,NP
            VHXC(I) = VHXC1(I)
            BXCX(I) = BXCX1(I)
            BXCY(I) = BXCY1(I)
            BXCZ(I) = BXCZ1(I)
         ENDDO    

      ENDIF

      CALL DENCALC(PSI1,PSI2,N,MX,MY,MZ,MM,MTOT)
      If (L.EQ.2) THEN
C**   project out transverse BXC
      DO I=1,NP
         BXCX(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MX(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCY(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MY(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCZ(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MZ(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
      ENDDO    

C      DO 90 I=1,NP
C         BB = DSQRT(BXCX1(I)**2 + BXCY1(I)**2 + BXCZ1(I)**2)
C         BXCX(I) = BB*MX(I)/MM(I)
C         BXCY(I) = BB*MY(I)/MM(I)
C         BXCZ(I) = BB*MZ(I)/MM(I)
C90    CONTINUE
      ENDIF

      CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
      WRITE(LN,999)TIME,N(1),N(2),N(3),N(4)
      WRITE(LN+1,999)TIME,MX(1),MX(2),MX(3),MX(4)
      WRITE(LN+2,999)TIME,MY(1),MY(2),MY(3),MY(4)
      WRITE(LN+3,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      WRITE(LN+4,999)TIME,MM(1),MM(2),MM(3),MM(4)
      WRITE(LN+5,999)TIME,MTOT(1),MTOT(2),MTOT(3),
     &             DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
      IF (L.EQ.1) THEN
      WRITE(LN+6,999)TIME,TX(1),TX(2),TX(3),TX(4)
      WRITE(LN+7,999)TIME,TY(1),TY(2),TY(3),TY(4)
      WRITE(LN+8,999)TIME,TZ(1),TZ(2),TZ(3),TZ(4)
      WRITE(LN+9,997)TIME,TT(1),TT(2),TT(3)
      ENDIF
      IF (TIME.LT.199.9999D0) GOTO 100

997   FORMAT(F8.3,3F10.5)
998   FORMAT(F8.3,5F10.5)
999   FORMAT(F8.3,4F10.5)
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
      IMPLICIT NONE
 
      INTEGER I,NP
      INCLUDE 'dim.inc'   

      DOUBLE PRECISION TT(3),TX(NP),TY(NP),TZ(NP),MX(NP),MY(NP),
     &                 MZ(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

      TT(1)=0.D0
      TT(2)=0.D0
      TT(3)=0.D0
      DO I=1,NP
         TX(I) = MY(I)*BXCZ(I) - MZ(I)*BXCY(I)
         TY(I) = MZ(I)*BXCX(I) - MX(I)*BXCZ(I)
         TZ(I) = MX(I)*BXCY(I) - MY(I)*BXCX(I)
         TT(1) = TT(1) + TX(I)
         TT(2) = TT(2) + TY(I)
         TT(3) = TT(3) + TZ(I)
      ENDDO    

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE MATRIX(M,T,V,BX,BY,BZ)
      IMPLICIT NONE

      INTEGER I,J,NP
      INCLUDE 'dim.inc'   
      DOUBLE PRECISION T,V(NP),BX(NP),BY(NP),BZ(NP)
      DOUBLE COMPLEX M(2*NP,2*NP),ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      DO I=1,2*NP
      DO J=1,2*NP
         M(I,J) = ZERO
      ENDDO
      ENDDO

      DO I=1,NP
         M(I,I) = (V(I) + BZ(I))*ONE
         M(I,I+NP) = BX(I)*ONE - BY(I)*IONE
         M(I+NP,I+NP) = (V(I) - BZ(I))*ONE
         M(I+NP,I) = BX(I)*ONE + BY(I)*IONE
      ENDDO    

      DO I=1,NP-1
         M(I,I+1) = -T*ONE
         M(I+NP,I+NP+1) = -T*ONE
         M(I+1,I) = -T*ONE
         M(I+NP+1,I+NP) = -T*ONE
      ENDDO    

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE BCORR (MX,MY,MZ,BX,BY,BZ,TT)
      IMPLICIT NONE

      INTEGER I,J,II,JJ,NP,NM,INFO,COP
      INCLUDE 'dim.inc'
      PARAMETER (NM = 3*NP-3)
      INTEGER IPIV(NM),LWORK
      PARAMETER (LWORK=100)
      DOUBLE PRECISION MX(NP),MY(NP),MZ(NP),BX(NP),BY(NP),BZ(NP),
     &                 MAT(NM,NM),RVEC(NM),WORK(LWORK),
     &                 BBX(NP),BBY(NP),BBZ(NP)
      DOUBLE PRECISION AX(NP),AY(NP),AZ(NP),AD,MSX,MSY,MSZ,M2,TT(3),
     &                 PX(NP),PY(NP),PZ(NP),QX(NP),QY(NP),QZ(NP)

      COP=0
C**--------------------------------------------------------------------***
C**   Do we have a coplanar situation?
C**   COP=1,2,3: torque along x,y,z
C**--------------------------------------------------------------------***
      MSX = 0.D0
      MSY = 0.D0
      MSZ = 0.D0
      M2 = 0.D0
      DO I=1,NP
         MSX = MSX + DABS(MX(I))
         MSY = MSY + DABS(MY(I))
         MSZ = MSZ + DABS(MZ(I))
         M2 = M2 + MX(I)**2 + MY(I)**2 + MZ(I)**2
      ENDDO    

      IF (MSX.LT.1.D-10) COP=1
      IF (MSY.LT.1.D-10) COP=2
      IF (MSZ.LT.1.D-10) COP=3

      IF (COP.EQ.0) THEN
C**--------------------------------------------------------------------***
C**   Not coplanar, so we have to solve a system of equations
C**--------------------------------------------------------------------***

      DO I=1,NP
         AD = MY(NP)*MZ(1) - MZ(NP)*MY(1)
         AX(I) = -(MY(NP)*MZ(I) - MZ(NP)*MY(I))/AD
         AY(I) = -(MZ(NP)*MX(I) - MX(NP)*MZ(I))/AD
         AZ(I) = -(MX(NP)*MY(I) - MY(NP)*MX(I))/AD
      ENDDO    

      DO I=1,NP
         PX(I) = MY(1)*AX(I) + MY(I)
         PY(I) = MY(1)*AY(I) - MX(I)
         PZ(I) = MY(1)*AZ(I)

         QX(I) = MZ(1)*AX(I) + MZ(I)
         QY(I) = MZ(1)*AY(I)
         QZ(I) = MZ(1)*AZ(I) - MX(I)
      ENDDO    

      DO I=1,NP-1
         II = 3*(I-1) + 1

         RVEC(II)   = -MX(NP)**2*BY(I) - MX(NP)**2*AY(I)*BX(1)
     &              - (MY(1)*AY(I)-MX(I))*MX(NP)*BY(NP)
     &              - MZ(1)*AY(I)*MX(NP)*BZ(NP)

         RVEC(II+1) = -MX(NP)**2*BZ(I) - MX(NP)**2*AZ(I)*BX(1)
     &              - MY(1)*AZ(I)*MX(NP)*BY(NP)
     &              - (MZ(1)*AZ(I)-MX(I))*MX(NP)*BZ(NP)

         RVEC(II+2) = -MX(NP)**2*BX(I+1) - MX(NP)**2*AX(I+1)*BX(1)
     &              - (MY(1)*AX(I+1)+MY(I+1))*MX(NP)*BY(NP)
     &              - (MZ(1)*AX(I+1)+MZ(I+1))*MX(NP)*BZ(NP)
      ENDDO    

      DO I=1,NP-1
      DO J=1,NP-1
         II = 3*(I-1) + 1
         JJ = 3*(J-1) + 1

         MAT(II,JJ) = MX(NP)**2*AY(I)*AY(J)
     &              + (MY(1)*AY(I)-MX(I))*PY(J)
     &              + MZ(1)*AY(I)*QY(J)

         MAT(II,JJ+1) = MX(NP)**2*AY(I)*AZ(J)
     &              + (MY(1)*AY(I)-MX(I))*PZ(J)
     &              + MZ(1)*AY(I)*QZ(J)

         MAT(II,JJ+2) = MX(NP)**2*AY(I)*AX(J+1)
     &              + (MY(1)*AY(I)-MX(I))*PX(J+1)
     &              + MZ(1)*AY(I)*QX(J+1)

         MAT(II+1,JJ) = MX(NP)**2*AZ(I)*AY(J)
     &              + MY(1)*AZ(I)*PY(J)
     &              + (MZ(1)*AZ(I)-MX(I))*QY(J)

         MAT(II+1,JJ+1) = MX(NP)**2*AZ(I)*AZ(J)
     &              + MY(1)*AZ(I)*PZ(J)
     &              + (MZ(1)*AZ(I)-MX(I))*QZ(J)

         MAT(II+1,JJ+2) = MX(NP)**2*AZ(I)*AX(J+1)
     &              + MY(1)*AZ(I)*PX(J+1)
     &              + (MZ(1)*AZ(I)-MX(I))*QX(J+1)

         MAT(II+2,JJ) = MX(NP)**2*AX(I+1)*AY(J)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PY(J)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QY(J)

         MAT(II+2,JJ+1) = MX(NP)**2*AX(I+1)*AZ(J)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PZ(J)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QZ(J)

         MAT(II+2,JJ+2) = MX(NP)**2*AX(I+1)*AX(J+1)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PX(J+1)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QX(J+1)

      ENDDO    
      ENDDO    

      DO I=1,NM
         MAT(I,I) = MAT(I,I) + MX(NP)**2
      ENDDO    

      CALL DSYSV('L',NM,1,MAT,NM,IPIV,RVEC,NM,WORK,LWORK,INFO )

      DO I=1,NP-1
         II = 3*(I-1)+1
         BBY(I) = RVEC(II)
         BBZ(I) = RVEC(II+1)
         BBX(I+1) = RVEC(II+2)
      ENDDO    

      BBX(1) = 0.D0
      DO I=1,NP-1
         BBX(1) = BBX(1) + AY(I)*BBY(I) + AZ(I)*BBZ(I)
     &                   + AX(I+1)*BBX(I+1)
      ENDDO    

      BBY(NP) = MY(NP)*BBX(NP)/MX(NP)
      BBZ(NP) = MZ(NP)*BBX(NP)/MX(NP)
      DO I=1,NP-1
         BBY(NP) = BBY(NP) + MY(I)*BBX(I)/MX(NP)
     &                     - MX(I)*BBY(I)/MX(NP)
         BBZ(NP) = BBZ(NP) + MZ(I)*BBX(I)/MX(NP)
     &                     - MX(I)*BBZ(I)/MX(NP)
      ENDDO    

      DO I=1,NP
         BX(I) = - BBX(I)
         BY(I) = - BBY(I)
         BZ(I) = - BBZ(I)
      ENDDO    

      ELSE
C**--------------------------------------------------------------------***
C**   coplanar, so we have an explicit solution
C**--------------------------------------------------------------------***
      DO I=1,NP
         BX(I) = BX(I) - (TT(2)*MZ(I) - TT(3)*MY(I))/M2
         BY(I) = BY(I) - (TT(3)*MX(I) - TT(1)*MZ(I))/M2
         BZ(I) = BZ(I) - (TT(1)*MY(I) - TT(2)*MX(I))/M2
      ENDDO    

      ENDIF

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE DENCALC(PSI1,PSI2,N,MX,MY,MZ,MM,MTOT)
      IMPLICIT NONE
      
      INTEGER I,NP
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX PSI1(2*NP),PSI2(2*NP),
     &               NUU(NP),NUD(NP),NDU(NP),NDD(NP)
      DOUBLE PRECISION N(NP),MX(NP),MY(NP),MZ(NP),MM(NP),MTOT(3)
      
      DO I=1,NP
         NUU(I) = CDABS(PSI1(I))**2 + CDABS(PSI2(I))**2
         NUD(I) = PSI1(I)*DCONJG(PSI1(I+NP))+PSI2(I)*DCONJG(PSI2(I+NP))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(PSI1(I+NP))**2 + CDABS(PSI2(I+NP))**2
      ENDDO    
      
      DO I=1,NP
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
         MM(I) = DSQRT(MX(I)**2+MY(I)**2+MZ(I)**2)
      ENDDO       

      MTOT(1) = 0.D0
      MTOT(2) = 0.D0
      MTOT(3) = 0.D0
      DO I=1,NP
         MTOT(1) = MTOT(1) + MX(I)
         MTOT(2) = MTOT(2) + MY(I)
         MTOT(3) = MTOT(3) + MZ(I)
      ENDDO    

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE XCPOT_SLATER(KLI,PSI1,PSI2,U0,U1,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE
      
      INTEGER NP,I,J,K,L,R,KLI
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP),
     &               BUU(NP),BUD(NP),BDU(NP),BDD(NP),IONE,
     &               GAMMA(2,2,NP,NP),PSI1(2*NP),PSI2(2*NP),
     &               PHI(2,2,NP)
      PARAMETER (IONE=(0.D0,1.D0))
      DOUBLE PRECISION U0,U1,N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP)

      DOUBLE COMPLEX VXCMAT_SLATER(2,2,NP),A(2),B(2),RMAT(2,2,2,NP),
     &               VRMAT(2,4,NP),VNMR(2,4,NP),DUM,NM(4,4,NP),
     &               NMR(2,2,2,NP),VXCMAT(2,2,NP)

      DO I=1,NP
         NUU(I) = CDABS(PSI1(I))**2 + CDABS(PSI2(I))**2
         NUD(I) = PSI1(I)*DCONJG(PSI1(I+NP))+PSI2(I)*DCONJG(PSI2(I+NP))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(PSI1(I+NP))**2 + CDABS(PSI2(I+NP))**2
      ENDDO    

      DO I=1,NP
         PHI(1,1,I) = PSI1(I)
         PHI(1,2,I) = PSI1(I+NP)
         PHI(2,1,I) = PSI2(I)
         PHI(2,2,I) = PSI2(I+NP)
      ENDDO    

      DO I=1,2
      DO J=1,2
      DO K=1,NP
      DO L=1,NP
         GAMMA(I,J,K,L) = PHI(1,I,K)*DCONJG(PHI(1,J,L))
     &                  + PHI(2,I,K)*DCONJG(PHI(2,J,L))
      ENDDO    
      ENDDO    
      ENDDO    
      ENDDO    

      
      DO I=1,NP
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
      ENDDO    
      
      VH(1) = U0*N(1) + U1*N(2)
      DO I=2,NP-1
         VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
      ENDDO    
      VH(NP) = U0*N(NP) + U1*N(NP-1)

      DO I=1,NP

      BUU(I) = -2.D0*U0*(NUU(I)*NUU(I) + NUD(I)*NDU(I))
      BDU(I) = -2.D0*U0*(NDU(I)*NUU(I) + NDD(I)*NDU(I))
      BUD(I) = -2.D0*U0*(NUU(I)*NUD(I) + NUD(I)*NDD(I))
      BDD(I) = -2.D0*U0*(NDU(I)*NUD(I) + NDD(I)*NDD(I))

      IF (I.LT.NP) THEN
      BUU(I) = BUU(I) -2.D0*U1*(GAMMA(1,1,I,I+1)*GAMMA(1,1,I+1,I)
     &                         +GAMMA(1,2,I,I+1)*GAMMA(2,1,I+1,I))
      BDU(I) = BDU(I) -2.D0*U1*(GAMMA(2,1,I,I+1)*GAMMA(1,1,I+1,I)
     &                         +GAMMA(2,2,I,I+1)*GAMMA(2,1,I+1,I))
      BUD(I) = BUD(I) -2.D0*U1*(GAMMA(1,1,I,I+1)*GAMMA(1,2,I+1,I)
     &                         +GAMMA(1,2,I,I+1)*GAMMA(2,2,I+1,I))
      BDD(I) = BDD(I) -2.D0*U1*(GAMMA(2,1,I,I+1)*GAMMA(1,2,I+1,I)
     &                         +GAMMA(2,2,I,I+1)*GAMMA(2,2,I+1,I))
      ENDIF

      IF (I.GT.1) THEN
      BUU(I) = BUU(I) -2.D0*U1*(GAMMA(1,1,I,I-1)*GAMMA(1,1,I-1,I)
     &                         +GAMMA(1,2,I,I-1)*GAMMA(2,1,I-1,I))
      BDU(I) = BDU(I) -2.D0*U1*(GAMMA(2,1,I,I-1)*GAMMA(1,1,I-1,I)
     &                         +GAMMA(2,2,I,I-1)*GAMMA(2,1,I-1,I))
      BUD(I) = BUD(I) -2.D0*U1*(GAMMA(1,1,I,I-1)*GAMMA(1,2,I-1,I)
     &                         +GAMMA(1,2,I,I-1)*GAMMA(2,2,I-1,I))
      BDD(I) = BDD(I) -2.D0*U1*(GAMMA(2,1,I,I-1)*GAMMA(1,2,I-1,I)
     &                         +GAMMA(2,2,I,I-1)*GAMMA(2,2,I-1,I))
      ENDIF

      ENDDO    

      DO R=1,NP
         DEN(R) = 2.D0*N(R)*(NUU(R)*NDD(R)-NUD(R)*NDU(R))

         NM(1,1,R) = (N(R)*NDD(R) - NUD(R)*NDU(R))/DEN(R)
         NM(1,2,R) = -NDD(R)*NUD(R)/DEN(R)
         NM(1,3,R) = -NDD(R)*NDU(R)/DEN(R)
         NM(1,4,R) = NUD(R)*NDU(R)/DEN(R)

         NM(2,1,R) = -NDD(R)*NDU(R)/DEN(R)
         NM(2,2,R) = (2.D0*NUU(R)*NDD(R) - NUD(R)*NDU(R))/DEN(R)
         NM(2,3,R) = NDU(R)**2/DEN(R)
         NM(2,4,R) = -NUU(R)*NDU(R)/DEN(R)

         NM(3,1,R) = -NDD(R)*NUD(R)/DEN(R)
         NM(3,2,R) = NUD(R)**2/DEN(R)
         NM(3,3,R) = (2.D0*NUU(R)*NDD(R) - NUD(R)*NDU(R))/DEN(R)
         NM(3,4,R) = -NUU(R)*NUD(R)/DEN(R)

         NM(4,1,R) = NUD(R)*NDU(R)/DEN(R)
         NM(4,2,R) = -NUU(R)*NUD(R)/DEN(R)
         NM(4,3,R) = -NUU(R)*NDU(R)/DEN(R)
         NM(4,4,R) = (N(R)*NUU(R) - NUD(R)*NDU(R))/DEN(R)

         VUU(R) = NM(1,1,R)*BUU(R) + NM(1,2,R)*BDU(R)
     &          + NM(1,3,R)*BUD(R) + NM(1,4,R)*BDD(R)
         VDU(R) = NM(2,1,R)*BUU(R) + NM(2,2,R)*BDU(R)
     &          + NM(2,3,R)*BUD(R) + NM(2,4,R)*BDD(R)
         VUD(R) = NM(3,1,R)*BUU(R) + NM(3,2,R)*BDU(R)
     &          + NM(3,3,R)*BUD(R) + NM(3,4,R)*BDD(R)
         VDD(R) = NM(4,1,R)*BUU(R) + NM(4,2,R)*BDU(R)
     &          + NM(4,3,R)*BUD(R) + NM(4,4,R)*BDD(R)

         VXC (R) = DREAL(VUU(R) + VDD(R))/2.D0
         BXCX(R) = DREAL(VDU(R) + VUD(R))/2.D0
         BXCY(R) = DREAL(-IONE*VDU(R) + IONE*VUD(R))/2.D0
         BXCZ(R) = DREAL(VUU(R) - VDD(R))/2.D0
      ENDDO    
C---------------------------------------------------------------------
C**   If KLI=1, construct the KLI XC potential.
C---------------------------------------------------------------------
      IF (KLI.EQ.1) THEN

      DO R=1,NP
         VXCMAT_SLATER(1,1,R) = VUU(R)
         VXCMAT_SLATER(1,2,R) = VUD(R)
         VXCMAT_SLATER(2,1,R) = VDU(R)
         VXCMAT_SLATER(2,2,R) = VDD(R)
         DO I=1,2
         DO J=1,2
            RMAT(1,I,J,R) = PHI(1,I,R)*DCONJG(PHI(1,J,R))
            RMAT(2,I,J,R) = PHI(2,I,R)*DCONJG(PHI(2,J,R))
         ENDDO
         ENDDO
      ENDDO    

      DO I=1,2
      DO R=1,NP
         VRMAT(I,1,R) = RMAT(I,1,1,R)
         VRMAT(I,2,R) = RMAT(I,2,1,R)
         VRMAT(I,3,R) = RMAT(I,1,2,R)
         VRMAT(I,4,R) = RMAT(I,2,2,R)
      ENDDO    
      ENDDO    

      DO I=1,2
      DO R=1,NP
      DO J=1,4
         DUM = (0.D0,0.D0)
         DO  K=1,4
            DUM = DUM + NM(J,K,R)*VRMAT(I,K,R)
         ENDDO
         VNMR(I,J,R) = DUM
      ENDDO    
      ENDDO    
      ENDDO    

      DO I=1,2
      DO R=1,NP
         NMR(I,1,1,R) = VNMR(I,1,R)
         NMR(I,2,1,R) = VNMR(I,2,R)
         NMR(I,1,2,R) = VNMR(I,3,R)
         NMR(I,2,2,R) = VNMR(I,4,R)
      ENDDO    
      ENDDO    

      CALL BCALC(PHI,GAMMA,U0,U1,B)
      CALL ACALC(RMAT,NMR,VXCMAT_SLATER,A,B)

      DO I=1,2
      DO J=1,2
      DO R=1,NP
         VXCMAT(I,J,R) = VXCMAT_SLATER(I,J,R)
     &                 + (A(1)-B(1))*NMR(1,I,J,R)
     &                 + (A(2)-B(2))*NMR(2,I,J,R)
      ENDDO   
      ENDDO   
      ENDDO   

      DO R=1,NP
         VXC(R) = DREAL((VXCMAT(1,1,R) + VXCMAT(2,2,R))/2.D0)
         BXCX(R) = DREAL((VXCMAT(1,2,R) + VXCMAT(2,1,R))/2.D0)
         BXCY(R) = DREAL(IONE*(VXCMAT(1,2,R) - VXCMAT(2,1,R))/2.D0)
         BXCZ(R) = DREAL((VXCMAT(1,1,R) - VXCMAT(2,2,R))/2.D0)
      ENDDO    

      ENDIF

      DO I=1,NP
         VHXC(I) = VH(I) + VXC(I)
      ENDDO

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE ACALC(RMAT,NMR,VXCMAT_SLATER,A,B)
      IMPLICIT NONE

      INTEGER NP,I,J,K,L,R
      INCLUDE 'dim.inc'  
      DOUBLE COMPLEX A(2),B(2),RMAT(2,2,2,NP),NMR(2,2,2,NP),M12,
     &               VXCMAT_SLATER(2,2,NP),Q(2),P(2,2),DUM,ZERO,ONE

      ZERO = (0.D0,0.D0)
      ONE = (1.D0,0.D0)

      DO I=1,2
      DO J=1,2
         DUM = ZERO
         DO K=1,2
         DO L=1,2
         DO R=1,NP
            DUM = DUM + RMAT(J,L,K,R)*NMR(I,K,L,R)
         ENDDO    
         ENDDO    
         ENDDO    
         P(J,I) = -2.D0*DUM
      ENDDO   
      ENDDO   

      DO J=1,2
         DUM = ZERO
         DO K=1,2
         DO L=1,2
         DO R=1,NP
            DUM = DUM + RMAT(J,L,K,R)*(VXCMAT_SLATER(K,L,R)
     &                - B(1)*NMR(1,K,L,R) - B(2)*NMR(2,K,L,R))
         ENDDO   
         ENDDO   
         ENDDO   
         Q(J) = 2.D0*DUM
      ENDDO    

      M12 = NMR(1,1,1,1)/NMR(2,1,1,1)

      A(1) = (Q(2) - (M12*B(1) + B(2))*(ONE + P(2,2)))
     &       /(P(2,1) - M12*(ONE+P(2,2)))

      A(2) = (B(1) - A(1))*M12 + B(2)

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE BCALC(PHI,GAMMA,U0,U1,B)
      IMPLICIT NONE

      INTEGER NP,I,K,ALPHA,BETA
      INCLUDE 'dim.inc'
      DOUBLE PRECISION U0,U1
      DOUBLE COMPLEX B(2),GAMMA(2,2,NP,NP),PHI(2,2,NP),DUM

      DO I=1,2
         DUM = (0.D0,0.D0)
            DO ALPHA=1,2
            DO BETA=1,2
               DO K=1,NP
                  DUM = DUM + U0*DCONJG(PHI(I,ALPHA,K))
     &                          *GAMMA(ALPHA,BETA,K,K)*PHI(I,BETA,K)
               ENDDO
               DO K=1,NP-1
                  DUM = DUM + U1*DCONJG(PHI(I,ALPHA,K))
     &                        *GAMMA(ALPHA,BETA,K,K+1)*PHI(I,BETA,K+1)
     &                      + U1*DCONJG(PHI(I,ALPHA,K+1))
     &                        *GAMMA(ALPHA,BETA,K+1,K)*PHI(I,BETA,K)
               ENDDO   
            ENDDO   
            ENDDO
         B(I) = -DUM - DCONJG(DUM)
      ENDDO    

      RETURN
      END
