      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER NP,I,ITER,REST,XC,CORR,LWORK,INFO
      INCLUDE 'dim.inc'   
      PARAMETER (LWORK = 100)

      DOUBLE PRECISION T,TOL,EOLD,U0,U1,MIX,EX,CRIT,
     &                 EH,EVXC,EXC,ETOT,TX(NP),TY(NP),TZ(NP),TT(3),
     &                 VHXC(NP),VXC(NP),BXCX(NP),BXCY(NP),
     &                 BXCZ(NP),VHXCO(NP),BXCXO(NP),BXCYO(NP),BXCZO(NP),
     &                 V(NP),BX(NP),BY(NP),BZ(NP),
     &                 VT(NP),BTX(NP),BTY(NP),BTZ(NP),
     &                 N(NP),MX(NP),MY(NP),MZ(NP),E(2*NP),RWORK(100)

      DOUBLE COMPLEX M(2*NP,2*NP),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               WORK(LWORK)

      COMMON /EVALS/ E,PHI,GAMMA

      WRITE(*,*)'RESTART? (0=no, 1=yes)'
      READ(*,*)REST
      WRITE(*,*)'Torque correction? (0=no, 1=yes)'
      READ(*,*)CORR

      WRITE(*,*)'XC=1: Slater'
      WRITE(*,*)'XC=2: OEP'
      WRITE(*,*)'XC=3: KLI'
      READ(*,*)XC

      T = 0.5D0
      WRITE(*,*)'U0?'
      READ(*,*)U0
      U1 = U0/2.D0
      
      MIX = 0.10D0
      TOL = 1.D-10
      EOLD = 0.D0
      ITER = 0
      
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO   

      IF (REST.EQ.0) THEN
         DO I=1,NP
            VHXC(I) = 0.D0
            BXCX(I) = 0.D0
            BXCY(I) = 0.D0
            BXCZ(I) = 0.D0
         ENDDO    
      ELSEIF (REST.EQ.1) THEN
         READ(2)VHXC,BXCX,BXCY,BXCZ
         REWIND 2
      ENDIF

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
      
      CALL ZHEEV( 'V', 'U', 2*NP, M, 2*NP, E, WORK, LWORK, RWORK, INFO )
     
      DO I=1,2*NP
         WRITE(*,*)E(I)
      ENDDO
      WRITE(99,*)ITER,real(E(1))
C**----------------------------------------------------------------------
C**   calculate the densities and update the potentials
C**----------------------------------------------------------------------
      CALL GCALC(M,PHI,GAMMA)
      CALL DENCALC(M,N,MX,MY,MZ)

      DO I=1,NP
         VHXCO(I) = VHXC(I)
         BXCXO(I) = BXCX(I)
         BXCYO(I) = BXCY(I)
         BXCZO(I) = BXCZ(I)
      ENDDO                

      IF (XC.EQ.1.OR.XC.EQ.3) THEN
         CALL XCPOT_SLATER(XC,U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      ELSEIF (XC.EQ.2) THEN
         CALL XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      ENDIF

C**   symmetrization, if needed
C      VHXC(3)=VHXC(2)
C      VHXC(4)=VHXC(1)
C      BXCZ(1) = BXCX(4)
C      BXCZ(2) = BXCX(3)
C      BXCZ(3) = BXCX(2)
C      BXCZ(4) = BXCX(1)

      CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)

C**   do the xc torque correction

      IF (CORR.EQ.1) THEN
         WRITE(*,*)'total torque:'
         WRITE(*,*)TT
         WRITE(*,*)
         CALL BCORR(MX,MY,MZ,BXCX,BXCY,BXCZ,TT)
         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
         WRITE(*,*)'total torque:'
         WRITE(*,*)TT
         WRITE(*,*)
      ENDIF

      DO I=1,NP
         VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXCO(I)
         BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCXO(I)
         BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCYO(I)
         BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZO(I)
      ENDDO    

      CRIT = 0.D0
      DO I=1,NP
         CRIT = CRIT + DSQRT(MX(I)**2 + MY(I)**2 + MZ(I)**2)         
      ENDDO

      IF (DABS((CRIT- EOLD)/CRIT).GT.TOL) THEN
          EOLD = CRIT
          GOTO 1
      ENDIF    

      WRITE(2)VHXC,BXCX,BXCY,BXCZ

      WRITE(*,*)'       N             MX               MY            MZ'
      
      DO I=1,NP
      WRITE(*,*)real(N(I)),real(MX(I)),real(MY(I)),real(MZ(I))
      ENDDO
      
      WRITE(*,*)
      WRITE(*,*)'magnetization magnitude:'
      DO I=1,NP
      WRITE(*,*)sqrt(real(MX(I)**2 + MY(I)**2 + MZ(I)**2)) 
      ENDDO
     
      WRITE(*,*)
      WRITE(*,*)'      VXC         BXCX           BXCY             BXCZ'
      
      DO I=1,NP
      WRITE(*,*)real(VXC(I)),real(BXCX(I)),real(BXCY(I)),real(BXCZ(I))
      ENDDO

      WRITE(*,*)
      WRITE(*,*)'=====================GS Energy====================='

      EH = 0.D0
      DO I=1,NP
         EH = EH - 0.5D0*U0*N(I)**2
      ENDDO   
      DO I=1,NP-1
         EH = EH - U1*N(I)*N(I+1)
      ENDDO   

      EVXC = 0.D0
      DO I=1,NP
         EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &         -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)
      ENDDO

      CALL EX_CALC(U0,U1,EX)

      EXC = EX
      
      ETOT = E(1) + E(2) + EH + EVXC + EXC

      WRITE(*,*)' EKS = ',E(1) + E(2)
      WRITE(*,*)'  EH = ',EH   
      WRITE(*,*)'EVXC = ',EVXC 
      WRITE(*,*)' EXC = ',EXC
      WRITE(*,*)
      WRITE(*,*)'ETOT = ',real(ETOT)

      WRITE(*,*)
      WRITE(*,*)'(xc torque)_x    (xc torque)_y    (xc torque)_z'
      DO I=1,NP
         WRITE(*,*)real(TX(I)),' ',real(TY(I)),' ',real(TZ(I))
      ENDDO
      WRITE(*,*)'-----------------------------------------'
      WRITE(*,*)real(TT(1)),' ',real(TT(2)),' ',real(TT(3))
      WRITE(*,*)
   
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
      SUBROUTINE DENCALC(M,N,MX,MY,MZ)
      IMPLICIT NONE
      
      INTEGER I,NP
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX M(2*NP,2*NP),NUU(NP),NUD(NP),NDU(NP),NDD(NP)
      DOUBLE PRECISION N(NP),MX(NP),MY(NP),MZ(NP)
      
      DO I=1,NP
         NUU(I) = CDABS(M(I,1))**2 + CDABS(M(I,2))**2
         NUD(I) = M(I,1)*DCONJG(M(I+NP,1)) + M(I,2)*DCONJG(M(I+NP,2))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(M(I+NP,1))**2 + CDABS(M(I+NP,2))**2
      ENDDO   
      
      DO I=1,NP
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
      ENDDO   

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE GCALC(M,PHI,GAMMA)
      IMPLICIT NONE

      INTEGER I,J,K,L,NP
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX M(2*NP,2*NP),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)

C**   Define the orbitals phi(m,sigma,x):
C**   m = 1...2*NP is the orbital index
C**   sigma = 1,2 (up, down) is the spin index
C**   x = 1,...,NP (lattice points) is the spatial coordinate

      DO I=1,2*NP 
      DO J=1,NP
         PHI(I,1,J) = M(J,I)
         PHI(I,2,J) = M(J+NP,I)
      ENDDO    
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

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE XCPOT_SLATER(XC,U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE
      
      INTEGER NP,I,J,K,R,XC
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP),
     &               BUU(NP),BUD(NP),BDU(NP),BDD(NP),IONE,
     &               GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)
      PARAMETER (IONE=(0.D0,1.D0))
      DOUBLE PRECISION U0,U1,N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP)

      DOUBLE COMPLEX VXCMAT_SLATER(2,2,NP),A(2),B(2),RMAT(2,2,2,NP),
     &               VRMAT(2,4,NP),VNMR(2,4,NP),DUM,NM(4,4,NP),
     &               NMR(2,2,2,NP),VXCMAT(2,2,NP)

      COMMON /EVALS/ E,PHI,GAMMA

      DO I=1,NP
         NUU(I) = CDABS(PHI(1,1,I))**2 + CDABS(PHI(2,1,I))**2
         NUD(I) = PHI(1,1,I)*DCONJG(PHI(1,2,I))
     &          + PHI(2,1,I)*DCONJG(PHI(2,2,I))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(PHI(1,2,I))**2 + CDABS(PHI(2,2,I))**2
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
C**   If XC=3, construct the KLI XC potential.
C---------------------------------------------------------------------
      IF (XC.EQ.3) THEN

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
         DO K=1,4
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

      CALL BCALC(U0,U1,B)
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

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE BCALC(U0,U1,B)
      IMPLICIT NONE

      INTEGER NP,I,K,ALPHA,BETA
      INCLUDE 'dim.inc'
      DOUBLE PRECISION U0,U1,E(2*NP)
      DOUBLE COMPLEX B(2),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),DUM

      COMMON /EVALS/ E,PHI,GAMMA

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

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,ND,I,J,K,MU,NU,ALPHA,BETA,R,RP,INFO,LWORK
      INCLUDE 'dim.inc'
      PARAMETER (ND=4*NP,LWORK=100)
      DOUBLE PRECISION SING(ND)
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),DUM,
     &               LM(ND,ND),BMAT(2,2*NP),RVEC(ND),VXCMAT(2,2,NP),
     &               X(ND),X1(ND),UMAT(ND,ND),VTMAT(ND,ND),WORK(LWORK)
      DOUBLE PRECISION U0,U1,E(2*NP),N(NP),VH(NP),VXC(NP),RWORK(100),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

      COMMON /EVALS/ E,PHI,GAMMA

      DO I=1,NP
         N(I) = DREAL(GAMMA(1,1,I,I) + GAMMA(2,2,I,I))
      ENDDO

      VH(1) = U0*N(1) + U1*N(2)
      DO I=2,NP-1
         VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
      ENDDO   
      VH(NP) = U0*N(NP) + U1*N(NP-1)

C**   calculate the (ND x ND) OEP matrix

      DO MU=1,2
      DO NU=1,2
      DO ALPHA=1,2
      DO BETA=1,2
      DO R=1,NP
      DO RP=1,NP
         DUM = (0.D0,0.D0)
         DO I=1,2
         DO J=1,2*NP
            IF (I.NE.J) THEN
               DUM = DUM + PHI(I,BETA,RP)*DCONJG(PHI(J,ALPHA,RP))
     &               *DCONJG(PHI(I,MU,R))*PHI(J,NU,R)/(E(J)-E(I))
     &                   + DCONJG(PHI(I,ALPHA,RP))*PHI(J,BETA,RP)
     &               *PHI(I,NU,R)*DCONJG(PHI(J,MU,R))/(E(J)-E(I))
            ENDIF
         ENDDO    
         ENDDO    
         LM(MU+(NU-1)*2+(R-1)*4,ALPHA+(BETA-1)*2+(RP-1)*4) = DUM
      ENDDO    
      ENDDO    
      ENDDO    
      ENDDO    
      ENDDO    
      ENDDO    

C**   calculate the OEP right-hand side

      DO I=1,2
      DO J=1,2*NP
         DUM = (0.D0,0.D0)
         IF (I.NE.J) THEN
            DO ALPHA=1,2
            DO BETA=1,2
               DO K=1,NP
                  DUM = DUM + (U0/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K))
     &                          *GAMMA(ALPHA,BETA,K,K)*PHI(J,BETA,K)
               ENDDO   
               DO K=1,NP-1
                  DUM = DUM + (U1/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K))
     &                        *GAMMA(ALPHA,BETA,K,K+1)*PHI(J,BETA,K+1)
     &                      + (U1/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K+1))
     &                        *GAMMA(ALPHA,BETA,K+1,K)*PHI(J,BETA,K)
               ENDDO   
            ENDDO    
            ENDDO    
         ENDIF
         BMAT(I,J) = DUM
      ENDDO    
      ENDDO    

      DO MU=1,2
      DO NU=1,2
      DO R=1,NP
         DUM = (0.D0,0.D0)
         DO I=1,2
         DO J=1,2*NP
            DUM = DUM + DCONJG(BMAT(I,J)*PHI(I,MU,R))*PHI(J,NU,R)
     &                + BMAT(I,J)*DCONJG(PHI(J,MU,R))*PHI(I,NU,R)
         ENDDO    
         ENDDO    
         RVEC(MU+(NU-1)*2+(R-1)*4) = DUM
      ENDDO    
      ENDDO    
      ENDDO    

C**   now do the singular value decomposition

      CALL ZGESVD( 'A', 'A', ND, ND, LM, ND, SING, UMAT, ND, VTMAT, ND,
     &                   WORK, LWORK, RWORK, INFO )

      DO I=1,ND
         DUM = (0.D0,0.D0)
         DO J=1,ND
            DUM = DUM + DCONJG(UMAT(J,I))*RVEC(J)
         ENDDO   
         X1(I) = DUM
      ENDDO    

      WRITE(*,*)
C      WRITE(*,*)'Lowest OEP singular values:',real(SING(ND)),
C     &          real(SING(ND-1)),real(SING(ND-2)),real(SING(ND-3))
       WRITE(*,*)SING
      WRITE(*,*)

C**   The singular values are ordered from largest to smallest. 
C**   By default, we drop the last singular value SING(ND).
C**   But for strong correlations (large U), it may happen that
C**   other singular values need to be dropped as well. This needs
C**   to be checked from case to case, and, if necessary, the code
C**   below needs to be modified.

      DO I=1,ND
         IF (I.LT.ND) THEN
C         IF (I.LT.ND-3) THEN
            SING(I) = 1.D0/SING(I)
         ELSE
            SING(I) = 0.D0
         ENDIF
         X1(I) = SING(I)*X1(I)
      ENDDO    

      DO I=1,ND
         DUM = (0.D0,0.D0)
         DO J=1,ND
            DUM = DUM + DCONJG(VTMAT(J,I))*X1(J)
         ENDDO   
         X(I) = DUM
      ENDDO    

      DO MU=1,2
      DO NU=1,2
      DO R=1,NP
         VXCMAT(MU,NU,R) = X(MU+(NU-1)*2+(R-1)*4)
      ENDDO   
      ENDDO   
      ENDDO   

      DO R=1,NP
         VXC(R) = DREAL(VXCMAT(1,1,R) + VXCMAT(2,2,R))/2.D0
         VHXC(R) = VH(R) + VXC(R)
         BXCX(R) = DREAL(VXCMAT(1,2,R) + VXCMAT(2,1,R))/2.D0
         BXCY(R) = DREAL((0.D0,1.D0)*(VXCMAT(1,2,R) - VXCMAT(2,1,R)))
     &             /2.D0
         BXCZ(R) = DREAL(VXCMAT(1,1,R) - VXCMAT(2,2,R))/2.D0
      ENDDO   

      RETURN
      END 
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE EX_CALC(U0,U1,EX)
      IMPLICIT NONE
      
      INTEGER NP,I,TAU,SIGMA
      INCLUDE 'dim.inc'   
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)
      DOUBLE PRECISION E(2*NP),U0,U1,EX                        

      COMMON /EVALS/ E,PHI,GAMMA

      EX = 0.D0
      DO I=1,NP
      DO TAU=1,2
      DO SIGMA=1,2
         EX = EX -0.5D0*U0*GAMMA(SIGMA,TAU,I,I)*GAMMA(TAU,SIGMA,I,I)
      ENDDO   
      ENDDO   
      ENDDO   

      DO I=1,NP-1
      DO TAU=1,2
      DO SIGMA=1,2
         EX = EX - U1*GAMMA(SIGMA,TAU,I,I+1)*GAMMA(TAU,SIGMA,I+1,I)
      ENDDO   
      ENDDO   
      ENDDO   

      RETURN
      END
