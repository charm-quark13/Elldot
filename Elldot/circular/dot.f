      PROGRAM MAIN
      IMPLICIT NONE

      INTEGER I,J,K,L,ITER,IERR,INFO

      INTEGER POINTER,QNUDIM,QNDDIM,NUDIM,NDDIM

      DOUBLE PRECISION TOL,R0,NORM,HA0,BB,MUBOHR,G,AA,OMEGA,
     &                 EOLD,ETOT,EVH,EXC,EVXC,EAXC,EPS,THR,FODA

      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'

      PARAMETER (QNUDIM=MAX(NLU,1),QNDDIM=MAX(NLD,1),
     &           NUDIM =MAX(NU,1), NDDIM=MAX(ND,1))

      INTEGER QNU(QNUDIM),QND(QNDDIM),OCCU(NUDIM),OCCD(NDDIM)
      DOUBLE PRECISION R(NR),RR(NR),MAINDIAG(NR),
     &                 OFFDIAG(NR-1),WORK(2*NR),X(NR,NR),
     &                 UX(NR,NUDIM),DX(NR,NDDIM),
     &                 UEIG(NUDIM,2),DEIG(NDDIM,2),
     &                 UEPS(NUDIM),DEPS(NDDIM),VORT(NR),
     &                 VEXT(NR),VH(NR),VXCU(NR),VXCD(NR),AXC(NR),
     &                 VOLDU(NR),VOLDD(NR),VU(NR),VD(NR),JP(NR),
     &                 RHO(NR),RHOUP(NR),RHODOWN(NR),ZETA(NR),KK(NR)

C***********************************************************************
C********************** Begin of the main program **********************
C***********************************************************************
      TOL=RESID_0
      HA0=27.2116D3*M0/E0**2

      BB = B * 4.254341D-6 * E0**1.5D0 / M0**2
      G = -0.44D0
      MUBOHR = M0 * DSQRT(E0) /2.D0
      OMEGA = DSQRT( (OMEGA0/HA0)**2 + 0.25D0*BB**2*E0 )
C**------------------------------------------------------------------***
C**   read in the quantum numbers                                    ***
C**------------------------------------------------------------------***
      OPEN(1,FILE='quantumnumbers',STATUS='OLD',FORM='FORMATTED')
      DO 1 I=1,NLU
1        READ(1,899)QNU(I)
      DO 2 I=1,NLD
2        READ(1,899)QND(I)
C**------------------------------------------------------------------***
C**   read in the occupation numbers                                 ***
C**------------------------------------------------------------------***
      OPEN(2,FILE='occupation',STATUS='OLD',FORM='FORMATTED')
      DO 5 I=1,NU
5        READ(2,899)OCCU(I)
      DO 6 I=1,ND
6        READ(2,899)OCCD(I)
C**------------------------------------------------------------------***
C**   set up the logarithmic grid                                    ***
C**------------------------------------------------------------------***
      DO 10 I=1,NR  
10       R(I)= (DBLE(I)-N0)*DR   
         R0 = R(1)-DR

      DO 11 I=1,NR
11       RR(I) = DEXP(R(I))

C      DO 12 I=1,NR
C12       WRITE(10,*)I,R(I),RR(I)
C
C      DO 13 I=1,NR 
C         WRITE(13,*)RR(I),real(FODA(0,0,RR(I)))
C13       RHO(I) = RR(I)**2 * FODA(0,0,RR(I))**2
C
C      CALL NORMCALC(RHO,NORM)
C      WRITE(*,*)NORM
C**------------------------------------------------------------------***
C**   set up the bare potential of the dot                           ***
C**------------------------------------------------------------------***
      CALL DOTPOT(RR,VEXT)
C**------------------------------------------------------------------***
C**   Set the initial guesses for the eigenvectors :                 ***
C**------------------------------------------------------------------***
      IF (REST.EQ.0) CALL START(RR,UX,DX,QNU,QND)
      IF (REST.EQ.1) THEN
         READ(10)UX,DX
         REWIND 10
      ENDIF
C**------------------------------------------------------------------***
C**                                                                  ***
C**   here begins the Kohn-Sham self-consistency loop:               ***
C**                                                                  ***
C**------------------------------------------------------------------***
      EOLD=0.D0
      ITER=0
15    CONTINUE
      ITER=ITER+1
      WRITE(*,907)ITER,TOL
C**------------------------------------------------------------------***
C**   Now calculate the Kohn-Sham potential.                         ***
C**------------------------------------------------------------------***
      CALL DENCALC(RHOUP,RHODOWN,RHO,ZETA,JP,VORT,UX,DX,QNU,QND,
     &             OCCU,OCCD,RR)          
      
      IF (XC.EQ.0) THEN

         CALL POTENTIAL(R,RR,RHO,ZETA,VH,VXCU,VXCD,EXC)

         DO 222 I=1,NR
222         AXC(I) = 0.D0

      ELSEIF (XC.EQ.1) THEN

         CALL SCALPOT(R,RR,RHO,ZETA,VORT,VH,VXCU,VXCD,EXC)
         CALL VECPOT(RR,RHO,ZETA,VORT,AXC)

         DO 25 I=1,NR
            AA = AXC(I)*JP(I)/RHO(I)
            VXCU(I) = VXCU(I) - AA
25          VXCD(I) = VXCD(I) - AA

      ENDIF

      DO 23 I=1,NR
         KK(I) = -0.5D0*RR(I)*RHO(I)*VH(I)
23    CONTINUE
      CALL DAVINT(RR,KK,NR,0.D0,RR(NR),EVH,IERR)

      DO 24 I=1,NR
24       KK(I) = -RR(I)*(RHOUP(I)*VXCU(I)+RHODOWN(I)*VXCD(I))
      CALL DAVINT(RR,KK,NR,0.D0,RR(NR),EVXC,IERR)

      EAXC=0.D0
      IF (XC.EQ.1) THEN
         DO 26 I=1,NR
26          KK(I) = -RR(I)*JP(I)*AXC(I)
         CALL DAVINT(RR,KK,NR,0.D0,RR(NR),EAXC,IERR)
      ENDIF

      ETOT=2.D0*PI*(EVH+EVXC+EAXC+EXC)

      DO 800 I=1,NR
C         VU(I) = VEXT(I) + VH(I) + VXCU(I) + 0.5D0*G*MUBOHR*BB
C800      VD(I) = VEXT(I) + VH(I) + VXCD(I) - 0.5D0*G*MUBOHR*BB
         VU(I) = VEXT(I) + 0.5D0*G*MUBOHR*BB
800      VD(I) = VEXT(I) - 0.5D0*G*MUBOHR*BB

      IF (ITER.EQ.1.AND.REST.EQ.0) THEN
         DO 801 I=1,NR
            VOLDU(I)=VU(I)
801         VOLDD(I)=VD(I)
      ELSEIF (ITER.EQ.1.AND.REST.EQ.1) THEN
         READ(11)VOLDU,VOLDD
         REWIND 11
      ENDIF

      DO 803 I=1,NR
         VU(I) = MIX*VU(I) + (1.D0-MIX)*VOLDU(I)
         VD(I) = MIX*VD(I) + (1.D0-MIX)*VOLDD(I)
         VOLDU(I)=VU(I)
803      VOLDD(I)=VD(I)
C**------------------------------------------------------------------***
C**   Set up the diagonal matrix elements and call the eigenvalue solver
C**------------------------------------------------------------------***
      POINTER=1
      DO 1001 J=1,NLU

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 1009 I=1,NR-1
1009        OFFDIAG(I) = -(0.5D0/DR**2)*DEXP(-R(I)-R(I+1))

         DO 1010 I=1,NR
1010        MAINDIAG(I) = (1.D0/DR**2+L**2*0.5D0)*DEXP(-2.D0*R(I))
     &                  + VU(I) - 0.5D0*L*DSQRT(E0)*BB - L*AXC(I)/RR(I)

         MAINDIAG(1) = MAINDIAG(1) - (0.5D0/DR**2)*DEXP(-2.D0*R(1))
     &                  * DEXP(-ABS(L)*(R(1)-R0))

         CALL DSTEV( 'V', NR, MAINDIAG, OFFDIAG, X, NR, WORK, INFO )

         DO 1012 I=1,QNU(J)
            DO 1011 K=1,NR
1011           UX(K,POINTER+I-1) = X(K,I)
1012        UEPS(POINTER+I-1) = MAINDIAG(I) 

         POINTER = POINTER + QNU(J)

1001  CONTINUE

      POINTER=1
      DO 1002 J=1,NLD

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 1019 I=1,NR-1
1019        OFFDIAG(I) = -(0.5D0/DR**2)*DEXP(-R(I)-R(I+1))

         DO 1020 I=1,NR
1020        MAINDIAG(I) = (1.D0/DR**2+L**2*0.5D0)*DEXP(-2.D0*R(I))
     &                  + VD(I) - 0.5D0*L*DSQRT(E0)*BB - L*AXC(I)/RR(I)

         MAINDIAG(1) = MAINDIAG(1) - (0.5D0/DR**2)*DEXP(-2.D0*R(1))
     &                  * DEXP(-ABS(L)*(R(1)-R0))

         CALL DSTEV( 'V', NR, MAINDIAG, OFFDIAG, X, NR, WORK, INFO )

         DO 1022 I=1,QND(J)
            DO 1021 K=1,NR
1021           DX(K,POINTER+I-1) = X(K,I)
1022        DEPS(POINTER+I-1) = MAINDIAG(I) 

         POINTER = POINTER + QND(J)

1002  CONTINUE

      CALL NORMALIZE(RR,UX,DX,QNU,QND)
C**------------------------------------------------------------------***
C**   calculate the total energy:                                    ***
C**------------------------------------------------------------------***
      DO 240 I=1,NU
240      ETOT = ETOT + UEPS(I)*OCCU(I)
      DO 241 I=1,ND
241      ETOT = ETOT + DEPS(I)*OCCD(I)
C**------------------------------------------------------------------***
C**   output of the results                                          ***
C**------------------------------------------------------------------***
      
      IF (NU.GT.0) WRITE(3,*)'spin up:'

      POINTER=0
      DO 250 J=1,NLU

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 251 I=1,QNU(J)
251       WRITE(3,904)OCCU(POINTER+I),QNU(J)-I,L,UEPS(POINTER+I)*HA0

         POINTER = POINTER + QNU(J)
250   CONTINUE

      IF (ND.GT.0) WRITE(3,*)'spin down:'

      POINTER=0
      DO 260 J=1,NLD

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 261 I=1,QND(J)
261       WRITE(3,904)OCCD(POINTER+I),QND(J)-I,L,DEPS(POINTER+I)*HA0

         POINTER = POINTER + QND(J)
260   CONTINUE
  
      REWIND 3
C**-----------------------------------------------------------------***
C**   Note: fort.3 contains the quantum numbers OCC, N, L           ***
C**   and the energy eigenvalues in exactly the same order in       ***
C**   which the eigenvectors are saved in UX and DX                 ***
C**-----------------------------------------------------------------***
      WRITE(*,*)' '
      WRITE(*,*)'    TOTAL ENERGY (meV):   ',ETOT*HA0
      WRITE(*,*)' '
C**-----------------------------------------------------------------***
C**   Jump back to the entry point for the self-consistency loop.   ***
C**   The residual is divided by the factor REF as soon as          ***
C**   the relative change of the total energy in the present        ***
C**   iteration step is less that THR, where THR=10*TOL.            ***
C**                                                                 ***
C**   The self-consistency iteration stops when TOL is lower        ***
C**   than the CRITERION.                                           ***
C**-----------------------------------------------------------------***
      EPS=DABS(1.D0-EOLD/ETOT)
      EOLD=ETOT
      THR=1.D1*TOL        
      IF (EPS.GT.THR) THEN
         WRITE(*,910)EPS,THR
      ELSE
         WRITE(*,911)EPS,THR
         TOL=TOL/REF
      ENDIF
      IF (TOL.GT.CRITERION) GOTO 15
C**-----------------------------------------------------------------***
C**   End of the self-consistency loop. Store results in files.     ***
C**-----------------------------------------------------------------***
      WRITE(*,*)
      WRITE(10)UX,DX
      WRITE(11)VOLDU,VOLDD
      WRITE(12)RHOUP,RHODOWN
      DO 2251 I=1,NR
         WRITE(25,*)real(RR(I)),real(VEXT(I))
         WRITE(26,*)real(RR(I)),real(VU(I)),real(VD(I)),real(AXC(I))
         WRITE(27,*)real(RR(I)),real(JP(I)),real(VORT(I))
2251     WRITE(30,*)real(RR(I)),real(RHOUP(I)),real(RHODOWN(I))

899   FORMAT(I1)
900   FORMAT(/,'convergence after ',I3,' iterations')
901   FORMAT('NOCC =',I2,'    DELTA =',E10.3)
902   FORMAT(E20.12)
903   FORMAT(4X,'n =',I2,'  l =',I2,'  e = ',F10.6,'     occ = ',I2)
904   FORMAT(3I5,F14.8)
907   FORMAT(/4X,
     &'*****************************************************',/,4X,
     &'***************** ITERATION NO.',I3,' ******************',/,4X,
     &'*****************************************************',/,4X,
     &'************** RESIDUAL = ',D12.5,' ***************',/,4X,
     &'*****************************************************',/)
910   FORMAT(4X,'EPS = ',D12.5,' > THR = ',D12.5)
911   FORMAT(4X,'EPS = ',D12.5,' < THR = ',D12.5)
1000  CONTINUE
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE START(RR,UX,DX,QNU,QND)
      IMPLICIT NONE

      INTEGER I,J,K,L
      INTEGER POINTER,QNUDIM,QNDDIM,NUDIM,NDDIM
      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'
      PARAMETER (QNUDIM=MAX(NLU,1),QNDDIM=MAX(NLD,1),
     &           NUDIM =MAX(NU,1), NDDIM=MAX(ND,1))
      INTEGER QNU(QNUDIM),QND(QNDDIM)
      DOUBLE PRECISION FODA,UX(NR,NUDIM),DX(NR,NDDIM),RR(NR)

      POINTER=0
      DO 1 J=1,NLU

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 10 I=1,QNU(J)
         DO 10 K=1,NR
10          UX(K,POINTER+I) = FODA(L,I,RR(K))*RR(K)

         POINTER = POINTER + QNU(J)
1     CONTINUE

      POINTER=0
      DO 2 J=1,NLD

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 20 I=1,QND(J)
         DO 20 K=1,NR
20          DX(K,POINTER+I) = FODA(L,I,RR(K))*RR(K)

         POINTER = POINTER + QND(J)
2     CONTINUE

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE DENCALC(RHOUP,RHODOWN,RHO,ZETA,JP,VORT,UX,DX,
     &                   QNU,QND,OCCU,OCCD,RR)          
      IMPLICIT NONE

      INTEGER I,J,K,L
      INTEGER POINTER,QNUDIM,QNDDIM,NUDIM,NDDIM
      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'
      PARAMETER (QNUDIM=MAX(NLU,1),QNDDIM=MAX(NLD,1),
     &           NUDIM =MAX(NU,1), NDDIM=MAX(ND,1))
      INTEGER QNU(QNUDIM),QND(QNDDIM),OCCU(NUDIM),OCCD(NDDIM)
      DOUBLE PRECISION UX(NR,NUDIM),DX(NR,NDDIM),RHOUP(NR),VORT(NR),
     &              DER(NR),RHODOWN(NR),RHO(NR),ZETA(NR),JP(NR),RR(NR)

      DO 10 I=1,NR
10       RHOUP(I)=0.D0

      DO 11 I=1,NR
11       RHODOWN(I)=0.D0

      DO 12 I=1,NR
12       JP(I) = 0.D0

      POINTER=0
      DO 2 J=1,NLU

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 20 I=1,QNU(J)
         DO 20 K=1,NR
20          RHOUP(K) = RHOUP(K) + 
     &                 OCCU(POINTER+I)*UX(K,POINTER+I)**2/RR(K)**2

         DO 25 I=1,QNU(J)
         DO 25 K=1,NR
25          JP(K) = JP(K) - 
     &              OCCU(POINTER+I)*L*UX(K,POINTER+I)**2/RR(K)**3

         POINTER = POINTER + QNU(J)
2     CONTINUE

      POINTER=0
      DO 3 J=1,NLD

         IF (J.EQ.1) L=0
         IF (J.GT.1.AND.MOD(J,2).EQ.0) L=J/2
         IF (J.GT.1.AND.MOD(J,2).EQ.1) L=-J/2

         DO 30 I=1,QND(J)
         DO 30 K=1,NR
30          RHODOWN(K) = RHODOWN(K) + 
     &                   OCCD(POINTER+I)*DX(K,POINTER+I)**2/RR(K)**2

         DO 35 I=1,QND(J)
         DO 35 K=1,NR
35          JP(K) = JP(K) - 
     &              OCCD(POINTER+I)*L*DX(K,POINTER+I)**2/RR(K)**3

         POINTER = POINTER + QND(J)
3     CONTINUE

      DO 50 I=1,NR
         RHO(I) = RHOUP(I) + RHODOWN(I)
50       ZETA(I)= (RHOUP(I) - RHODOWN(I))/RHO(I)

      DO 60 I=1,NR
60       VORT(I) = RR(I)*JP(I)/RHO(I)

      CALL DERIVATIVE(RR,VORT,DER)

      DO 70 I=1,NR
70       VORT(I) = DER(I)/(RR(I)*DSQRT(E0))

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE NORMALIZE(RR,UX,DX,QNU,QND)             
      IMPLICIT NONE

      INTEGER I,J,K,IERR
      DOUBLE PRECISION NORM
      INTEGER POINTER,QNUDIM,QNDDIM,NUDIM,NDDIM
      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'
      PARAMETER (QNUDIM=MAX(NLU,1),QNDDIM=MAX(NLD,1),
     &           NUDIM =MAX(NU,1), NDDIM=MAX(ND,1))
      INTEGER QNU(QNUDIM),QND(QNDDIM)
      DOUBLE PRECISION UX(NR,NUDIM),DX(NR,NDDIM),RR(NR),KK(NR)

      POINTER=0
      DO 1 J=1,NLU

         DO 10 I=1,QNU(J) 
         DO 11 K=1,NR
11          KK(K) = UX(K,POINTER+I)**2/RR(K)
         CALL DAVINT(RR,KK,NR,0.D0,RR(NR),NORM,IERR)
         DO 12 K=1,NR
12          UX(K,POINTER+I) = UX(K,POINTER+I)/DSQRT(NORM*2.D0*PI)
10       CONTINUE

         POINTER = POINTER + QNU(J)
1     CONTINUE


      POINTER=0
      DO 2 J=1,NLD

         DO 20 I=1,QND(J) 
         DO 21 K=1,NR
21          KK(K) = DX(K,POINTER+I)**2/RR(K)
         CALL DAVINT(RR,KK,NR,0.D0,RR(NR),NORM,IERR)
         DO 22 K=1,NR
22          DX(K,POINTER+I) = DX(K,POINTER+I)/DSQRT(NORM*2.D0*PI)
20       CONTINUE

         POINTER = POINTER + QND(J)
2     CONTINUE
      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE DOTPOT(R,VEXT)
      IMPLICIT NONE

      INTEGER I,PAR
      DOUBLE PRECISION RDISK,HA0,OM,BB,NPLUS,V0
      INCLUDE 'grid.inc'
      INCLUDE 'ellpar.inc'

      DOUBLE PRECISION R(NR),VEXT(NR)
      DOUBLE PRECISION M1,FIT1,FIT2
    
      BB = B * 4.254341D-6 * E0**1.5D0 / M0**2
      HA0=27.2116D3*M0/E0**2
      OM=OMEGA0/HA0

      PAR=1

      IF (PAR.EQ.1) THEN 


      DO 10 I=1,NR
10       VEXT(I) = 0.5D0 * OM**2 * R(I)**2   

      ELSEIF (PAR.EQ.2) THEN

      RDISK = 1.D2 
C      RDISK = 1.D3 

C      NPLUS = 4.D11
C**-----------------------------------------------------------------***
C**
C**   Here, RDISK is the radius of the jellium disk (in Angstrom) and
C**   NPLUS is the jellium density (in cm**-2)
C**
C**-----------------------------------------------------------------***
      RDISK = RDISK / A0
C      NPLUS = NPLUS * (E0/M0)**2  / 3.571068D16
      NPLUS = RDISK*OM**2/PI

      DO 20 I=1,NR
   
         IF (R(I).LT.RDISK) THEN
            M1 = 1.D0 - R(I)**2/RDISK**2
            FIT2 = 1.D0 + M1*(CC1 + M1*(CC2 + M1*(CC3 + M1*CC4))) -
     &             DLOG(M1)*M1*(DD1 + M1*(DD2 + M1*(DD3 + M1*DD4)))
            VEXT(I) = 4.D0*NPLUS*RDISK*FIT2
         ELSE
            M1 = 1.D0 - RDISK**2/R(I)**2
            FIT1 = AA0 + M1*(AA1 + M1*(AA2 + M1*(AA3 + M1*AA4))) - 
     &          DLOG(M1)*(BB0 + M1*(BB1 + M1*(BB2 + M1*(BB3 + M1*BB4))))
            FIT2 = 1.D0 + M1*(CC1 + M1*(CC2 + M1*(CC3 + M1*CC4))) -
     &             DLOG(M1)*M1*(DD1 + M1*(DD2 + M1*(DD3 + M1*DD4)))
            VEXT(I) = 4.D0*NPLUS*R(I)
     &                *(FIT2 - (1.D0-RDISK**2/R(I)**2)*FIT1)
         ENDIF
20    CONTINUE

      V0 = 2.*PI*NPLUS*RDISK

      DO 30 I=1,NR
30       VEXT(I) = -VEXT(I) + V0

      ENDIF

      DO 100 I=1,NR
100      VEXT(I) = VEXT(I) + 0.125D0 * E0 * BB**2 * R(I)**2      

C      DO 200 I=1,NR
C200      WRITE(28,*)real(R(I)),real(VEXT(I))

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE NORMCALC(PSI,NORM)
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION NORM,CC
      INCLUDE 'grid.inc'

      DOUBLE PRECISION PSI(NR),CO(NR)

      CC = DR*5.D0/288.D0

      CO(1)  = CC*19.D0
      DO 15 I=0,(NR-6)/5
         CO(5*I+2) = CC*75.D0
         CO(5*I+3) = CC*50.D0
         CO(5*I+4) = CC*50.D0
         CO(5*I+5) = CC*75.D0
15       CO(5*I+6) = CC*38.D0
      CO(NR) = CC*19.D0

      NORM = 0.D0
      DO 20 I=1,NR
20       NORM = NORM + CO(I)*PSI(I)

      NORM = NORM * 2.D0*PI

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE NORMCALC1(PSI,NORM)
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION NORM,CC
      INCLUDE 'grid.inc'

      DOUBLE PRECISION PSI(NR),CO(NR)

      CC = DR*5.D0/288.D0

      CO(1)  = CC*19.D0
      DO 15 I=0,(NR-6)/5
         CO(5*I+2) = CC*75.D0
         CO(5*I+3) = CC*50.D0
         CO(5*I+4) = CC*50.D0
         CO(5*I+5) = CC*75.D0
15       CO(5*I+6) = CC*38.D0
      CO(NR) = CC*19.D0

      NORM = 0.D0
      DO 20 I=1,NR
20       NORM = NORM + CO(I)*PSI(I)

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION FACTORIAL(N)
      IMPLICIT NONE

      INTEGER N,I

      IF (N.LE.1) THEN
         FACTORIAL = 1.D0
      ELSE
         FACTORIAL = 1.D0
         DO 10 I=1,N
10          FACTORIAL = FACTORIAL*I
      ENDIF

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION LAGUERRE(L,N,X)
      IMPLICIT NONE
C**------------------------------------------------------------------***
C**                                                       L
C**   this subroutine calculates the Laguerre polynomial L (x)
C**                                                       n       
C**------------------------------------------------------------------***

      INTEGER L,M,N
      DOUBLE PRECISION X,FACTORIAL

      IF (N.EQ.0) THEN
         LAGUERRE = 1.D0
      ELSEIF (N.EQ.1) THEN
         LAGUERRE = DBLE(L) + 1.D0 - X
      ELSE
         LAGUERRE = (FACTORIAL(N+L)/FACTORIAL(N))/FACTORIAL(L)
         DO 10 M=1,N
            LAGUERRE = LAGUERRE + (-X)**M
     &       *(FACTORIAL(N+L)/FACTORIAL(N-M))/
     &        (FACTORIAL(L+M)*FACTORIAL(M))
10       CONTINUE
      ENDIF

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION FODA(L,N,X)
      IMPLICIT NONE
C**------------------------------------------------------------------***
C**                                             
C**   this subroutine gives the Fock-Darwin states in a magnetic field,
C**   i.e. the eigenstates in a 2D oscillator potential. The normalization
C**   factor is included.                                   
C**                                                         
C**------------------------------------------------------------------***

      INTEGER L,LABS,N
      DOUBLE PRECISION HA0,OM,BB
      INCLUDE 'grid.inc'
      DOUBLE PRECISION X,LAGUERRE,NORM,FACTORIAL

      HA0=27.2116D3*M0/E0**2
      BB = B * 4.254341D-6 * E0**1.5D0 / M0**2
      OM = DSQRT( (OMEGA0/HA0)**2 + 0.25D0*E0*BB**2 )

      IF (L.LT.0) LABS = -L
      IF (L.GE.0) LABS = L

      FODA = X**LABS * DEXP(-OM*X**2/2.D0) * LAGUERRE(LABS,N,OM*X**2)

      NORM = PI*FACTORIAL(LABS+N)/(FACTORIAL(N)*OM**(LABS+1))

      FODA = FODA/DSQRT(NORM)

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE COULOMB(R,RR,RHO,VH)
      IMPLICIT NONE

      INTEGER I,J
      DOUBLE PRECISION X1,X2
      INCLUDE 'grid.inc'
      INCLUDE 'ellpar.inc'

      DOUBLE PRECISION R(NR),RR(NR),RHO(NR),VH(NR)
      DOUBLE PRECISION M1,FIT1,FIT2

      FIT1(M1) = AA0 + M1*(AA1 + M1*(AA2 + M1*(AA3 + M1*AA4))) -
     &        DLOG(M1)*(BB0 + M1*(BB1 + M1*(BB2 + M1*(BB3 + M1*BB4))))
      FIT2(M1) = 1.D0 + M1*(CC1 + M1*(CC2 + M1*(CC3 + M1*CC4))) -
     &             DLOG(M1)*M1*(DD1 + M1*(DD2 + M1*(DD3 + M1*DD4)))
C**-----------------------------------------------------------------***
      M1 = 1.D0 - RR(1)**2/RR(2)**2
      X1 = FIT2(M1)*RR(2)/RR(1)
      VH(1) = X1*(RHO(1)+RHO(2))/2.D0

      DO 5 J=2,NR-1
            M1 = 1.D0 - RR(1)**2/RR(J+1)**2
            X1 = FIT2(M1)*RR(J+1)/RR(1)
            M1 = 1.D0 - RR(1)**2/RR(J)**2
            X2 = FIT2(M1)*RR(J)/RR(1)
            VH(1) = VH(1) + (X1-X2)*(RHO(J)+RHO(J+1))/2.D0
5        CONTINUE

      DO 10 I=2,NR-1

         M1 = 1.D0 - RR(1)**2/RR(I)**2
         X1 = FIT2(M1) - M1*FIT1(M1)
         VH(I) = RHO(1)*X1

         DO 20 J=1,I-2
            M1 = 1.D0 - RR(J+1)**2/RR(I)**2
            X1 = FIT2(M1) - M1*FIT1(M1)
            M1 = 1.D0 - RR(J)**2/RR(I)**2
            X2 = FIT2(M1) - M1*FIT1(M1)
            VH(I) = VH(I) + (X1-X2)*(RHO(J)+RHO(J+1))/2.D0
20       CONTINUE

         DO 30 J=I+1,NR-1
            M1 = 1.D0 - RR(I)**2/RR(J+1)**2
            X1 = FIT2(M1)*RR(J+1)/RR(I)
            M1 = 1.D0 - RR(I)**2/RR(J)**2
            X2 = FIT2(M1)*RR(J)/RR(I)
            VH(I) = VH(I) + (X1-X2)*(RHO(J)+RHO(J+1))/2.D0
30       CONTINUE

         M1 = 1.D0 - RR(I-1)**2/RR(I)**2
         X1 = -FIT2(M1) + M1*FIT1(M1)
         VH(I) = VH(I) + X1*(RHO(I-1)+RHO(I))/2.D0

         M1 = 1.D0 - RR(I)**2/RR(I+1)**2
         X1 = FIT2(M1)*RR(I+1)/RR(I)
         VH(I) = VH(I) + X1*(RHO(I+1)+RHO(I))/2.D0

10    CONTINUE
   
      M1 = 1.D0 - RR(1)**2/RR(NR)**2
      X1 = FIT2(M1) - M1*FIT1(M1)
      VH(NR) = RHO(NR)*X1

      DO 40 J=1,NR-2
            M1 = 1.D0 - RR(J+1)**2/RR(NR)**2
            X1 = FIT2(M1) - M1*FIT1(M1)
            M1 = 1.D0 - RR(J)**2/RR(NR)**2
            X2 = FIT2(M1) - M1*FIT1(M1)
            VH(NR) = VH(NR) + (X1-X2)*(RHO(J)+RHO(J+1))/2.D0
40    CONTINUE

      DO 50 I=1,NR
50       VH(I) = 4.D0*RR(I)*VH(I)

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE COULOMB1(R,RR,RHO,VH)
      IMPLICIT NONE

      INTEGER I,IER,NEVAL,LENIW,LENW,LAST
      PARAMETER (LENIW=100,LENW=1000)
      DOUBLE PRECISION RPRIME
      INCLUDE 'grid.inc'
      INCLUDE 'ellpar.inc'

      INTEGER IWORK(LENIW)
      DOUBLE PRECISION R(NR),RR(NR),RHO(NR),VH(NR),Y2(NR),
     &                 WORK(LENW),X(NR),Y(NR)
      DOUBLE PRECISION F,POINTS(3),RESULT,ABSERR
      EXTERNAL F

      COMMON /RRR/ RPRIME
      COMMON /SPLI/ X,Y,Y2

      CALL SPLINE(RR,RHO,NR,1.D30,1.D30,Y2)

      DO 5 I=1,NR
          X(I)=RR(I)
5         Y(I)=RHO(I)

      DO 10 I=1,NR

         RPRIME = RR(I)

         POINTS(1) = RPRIME

         CALL DQAGP(F,0.D0,RR(NR)+DR,3,POINTS,0.D0,1.D-6,RESULT,ABSERR,
     &   NEVAL,IER,LENIW,LENW,LAST,IWORK,WORK)

         VH(I) = 4.D0*RESULT

10    CONTINUE

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      DOUBLE PRECISION FUNCTION F(RR)
      IMPLICIT NONE

      INCLUDE 'grid.inc'
      INCLUDE 'ellpar.inc'

      DOUBLE PRECISION M1,FIT1,DEN,RR,RPRIME
      DOUBLE PRECISION R(NR),RHO(NR),Y2(NR)

      COMMON /RRR/ RPRIME
      COMMON /SPLI/ R,RHO,Y2

      CALL SPLINT(R,RHO,Y2,NR,RR,DEN)

      M1 = 1.D0 - 4.D0*RR*RPRIME/(RR+RPRIME)**2
      FIT1 = AA0 + M1*(AA1 + M1*(AA2 + M1*(AA3 + M1*AA4))) -
     &      DLOG(M1)*(BB0 + M1*(BB1 + M1*(BB2 + M1*(BB3 + M1*BB4))))
      F = RR*DEN*FIT1/(RR+RPRIME)

      RETURN
      END
C****************************************************************
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION x,y,xa(n),ya(n),y2a(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if (xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.D0) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.D0
      return
      END
C**********************************************************************
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT NONE
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER(NMAX=5000)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      if (yp1.gt..99D30) then
         y2(1)=0.D0
         u(1)=0.D0
      else
         y2(1)=-0.5D0
         u(1)=(3.D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.D0
         y2(i)=(sig-1.D0)/p
         u(i)=(6.D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99D30) then
         qn=0.D0
         un=0.D0
      else
         qn=0.5D0
         un=(3.D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
      do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE POTENTIAL(R,RR,RHO,ZETA,VH,VXCU,VXCD,EXC)
      IMPLICIT NONE

      INTEGER I,IERR
      DOUBLE PRECISION RS,AA,EXC
      INCLUDE 'grid.inc'

      DOUBLE PRECISION A00,A01,A02,A03,A10,A11,A12,A13,X,
     &                 POL,PADE0,PADE1,DPADE0,DPADE1
      PARAMETER (A00 = -0.3568D0, A10 = -0.0515D0,
     &          A01 = 1.1300D0,  A11 = 340.5813D0,
     &          A02 = 0.9052D0,  A12 = 75.2293D0,
     &          A03 = 0.4165D0,  A13 = 37.0170D0)

      DOUBLE PRECISION R(NR),RR(NR),RHO(NR),ZETA(NR),KK(NR),
     &                 VXCU(NR),VXCD(NR),VH(NR)

      PADE0(X) = 0.5D0*A00*(1.D0 + A01*X)/
     &            (1.D0 + X*(A01 + X*(A02  + X*A03)))
      PADE1(X) = 0.5D0*A10*(1.D0 + A11*X)/
     &            (1.D0 + X*(A11 + X*(A12  + X*A13)))
      POL(X) = ( (1.D0+X)**1.5D0 + (1.D0-X)**1.5D0 - 2.D0 )
     &         /(2.D0*(DSQRT(2.D0) - 1.D0))
      DPADE0(X) = 0.5D0*A00*(
     &            A01/(1.D0 + X*(A01 + X*(A02  + X*A03))) - 
     &            (1.D0 + A01*X)*(A01 + 2.D0*A02*X + 3.D0*A03*X**2)/
     &            (1.D0 + X*(A01 + X*(A02  + X*A03)))**2 )
      DPADE1(X) = 0.5D0*A10*(
     &            A11/(1.D0 + X*(A11 + X*(A12  + X*A13))) - 
     &            (1.D0 + A11*X)*(A11 + 2.D0*A12*X + 3.D0*A13*X**2)/
     &            (1.D0 + X*(A11 + X*(A12  + X*A13)))**2 )

      CALL COULOMB1(R,RR,RHO,VH)

      DO 10 I=1,NR
         RS = 1.D0/DSQRT(PI*RHO(I))
         AA = -DSQRT(2.D0)/(PI*RS)
         VXCU(I) = AA*((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0) 
     &          + AA*(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )*
     &               (1.D0 - ZETA(I))
         VXCD(I) = AA*((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0) 
     &          - AA*(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )*
     &               (1.D0 + ZETA(I))
         AA = -DSQRT(8.D0)/(3.D0*PI*RS)
         KK(I) = AA*RR(I)*RHO(I)*
     &                 ((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0)
10    CONTINUE

      DO 20 I=1,NR
         X = 1.D0/(PI*RHO(I))**0.25D0
         VXCU(I)=VXCU(I) + PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) 
     &          - RHO(I)*PI*0.25D0*X**5 *
     &          ( DPADE0(X) + POL(ZETA(I))*(DPADE1(X)-DPADE0(X)) )
     &          + (1.D0 - ZETA(I))* ( 0.75D0 / (DSQRT(2.D0)-1.D0) )
     &           *(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )
     &           *(PADE1(X)-PADE0(X))
         VXCD(I)=VXCD(I) + PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) 
     &          - RHO(I)*PI*0.25D0*X**5 *
     &          ( DPADE0(X) + POL(ZETA(I))*(DPADE1(X)-DPADE0(X)) )
     &          - (1.D0 + ZETA(I))* ( 0.75D0 / (DSQRT(2.D0)-1.D0) )
     &           *(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )
     &           *(PADE1(X)-PADE0(X))
         KK(I) = KK(I) + RR(I)*RHO(I)*
     &           ( PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) )
20    CONTINUE

      CALL DAVINT(RR,KK,NR,0.D0,RR(NR),EXC,IERR)

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE SCALPOT(R,RR,RHO,ZETA,VORT,VH,VXCU,VXCD,EXC)
      IMPLICIT NONE

      INTEGER I,IERR,WARN
      DOUBLE PRECISION RS,AA,EXC
      INCLUDE 'grid.inc'

      DOUBLE PRECISION A00,A01,A02,A03,A10,A11,A12,A13,X,
     &                 POL,PADE0,PADE1,DPADE0,DPADE1
      PARAMETER (A00 = -0.3568D0, A10 = -0.0515D0,
     &          A01 = 1.1300D0,  A11 = 340.5813D0,
     &          A02 = 0.9052D0,  A12 = 75.2293D0,
     &          A03 = 0.4165D0,  A13 = 37.0170D0)

      DOUBLE PRECISION R(NR),RR(NR),RHO(NR),ZETA(NR),KK(NR),NU(NR),
     &                 VXCU(NR),VXCD(NR),VH(NR),VORT(NR)

      PADE0(X) = 0.5D0*A00*(1.D0 + A01*X)/
     &            (1.D0 + X*(A01 + X*(A02  + X*A03)))
      PADE1(X) = 0.5D0*A10*(1.D0 + A11*X)/
     &            (1.D0 + X*(A11 + X*(A12  + X*A13)))
      POL(X) = ( (1.D0+X)**1.5D0 + (1.D0-X)**1.5D0 - 2.D0 )
     &         /(2.D0*(DSQRT(2.D0) - 1.D0))
      DPADE0(X) = 0.5D0*A00*(
     &            A01/(1.D0 + X*(A01 + X*(A02  + X*A03))) - 
     &            (1.D0 + A01*X)*(A01 + 2.D0*A02*X + 3.D0*A03*X**2)/
     &            (1.D0 + X*(A01 + X*(A02  + X*A03)))**2 )
      DPADE1(X) = 0.5D0*A10*(
     &            A11/(1.D0 + X*(A11 + X*(A12  + X*A13))) - 
     &            (1.D0 + A11*X)*(A11 + 2.D0*A12*X + 3.D0*A13*X**2)/
     &            (1.D0 + X*(A11 + X*(A12  + X*A13)))**2 )

      CALL COULOMB1(R,RR,RHO,VH)

      WARN=0
      DO 5 I=1,NR
         NU(I) = -2.D0*PI*RHO(I)/(VORT(I)*DSQRT(E0))
         IF (NU(I).LT.0.D0) THEN
            WARN=1
            NU(I) = -NU(I)
         ENDIF
5     CONTINUE
      IF (WARN.EQ.1) WRITE(*,*)'warning: negative filling factor!'

      DO 10 I=1,NR
         RS = 1.D0/DSQRT(PI*RHO(I))
         AA = -DSQRT(2.D0)/(PI*RS)
         VXCU(I) = AA*((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0) 
     &          + AA*(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )*
     &               (1.D0 - ZETA(I))
         VXCD(I) = AA*((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0) 
     &          - AA*(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )*
     &               (1.D0 + ZETA(I))
         AA = -DSQRT(8.D0)/(3.D0*PI*RS)
         KK(I) = AA*RR(I)*RHO(I)*
     &                 ((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0)
10    CONTINUE

      DO 20 I=1,NR
         X = 1.D0/(PI*RHO(I))**0.25D0
         VXCU(I)=VXCU(I) + PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) 
     &          - RHO(I)*PI*0.25D0*X**5 *
     &          ( DPADE0(X) + POL(ZETA(I))*(DPADE1(X)-DPADE0(X)) )
     &          + (1.D0 - ZETA(I))* ( 0.75D0 / (DSQRT(2.D0)-1.D0) )
     &           *(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )
     &           *(PADE1(X)-PADE0(X))
         VXCD(I)=VXCD(I) + PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) 
     &          - RHO(I)*PI*0.25D0*X**5 *
     &          ( DPADE0(X) + POL(ZETA(I))*(DPADE1(X)-DPADE0(X)) )
     &          - (1.D0 + ZETA(I))* ( 0.75D0 / (DSQRT(2.D0)-1.D0) )
     &           *(DSQRT(1.D0+ZETA(I)) - DSQRT(1.D0-ZETA(I)) )
     &           *(PADE1(X)-PADE0(X))
         KK(I) = KK(I) + RR(I)*RHO(I)*
     &           ( PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X)) )
20    CONTINUE

      DO 30 I=1,NR
         AA = - 0.782133D0*DSQRT(2.D0*PI*RHO(I)) * (1.5D0 - 
     &             0.47264D0*NU(I)**0.74D0 + 0.0384D0*NU(I)**1.7D0 ) 
     &             /(1.D0 + NU(I)**4)
     &           + 4.D0*NU(I)**4 * ( 0.782133D0*DSQRT(2.D0*PI*RHO(I))
     &           *(1.D0 - 0.211D0*NU(I)**0.74D0 + 0.012D0*NU(I)**1.7D0 )
     &            + KK(I)/(RR(I)*RHO(I)) )/(1.D0 + NU(I)**4)**2
         VXCU(I) = VXCU(I)*NU(I)**4/(1.D0 + NU(I)**4) + AA
         VXCD(I) = VXCD(I)*NU(I)**4/(1.D0 + NU(I)**4) + AA
30    CONTINUE

      DO 100 I=1,NR
100      KK(I) = NU(I)**4*KK(I)/(1.D0 + NU(I)**4)
     &         - 0.782133D0*RR(I)*RHO(I)*DSQRT(2.D0*PI*RHO(I))
     &          *(1.D0 - 0.211D0*NU(I)**0.74D0 + 0.012D0*NU(I)**1.7D0 )
     &          /(1.D0 + NU(I)**4)

      CALL DAVINT(RR,KK,NR,0.D0,RR(NR),EXC,IERR)

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE VECPOT(RR,RHO,ZETA,VORT,AXC)
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION RS,AA
      INCLUDE 'grid.inc'

      DOUBLE PRECISION A00,A01,A02,A03,A10,A11,A12,A13,X,
     &                 POL,PADE0,PADE1
      PARAMETER (A00 = -0.3568D0, A10 = -0.0515D0,
     &          A01 = 1.1300D0,  A11 = 340.5813D0,
     &          A02 = 0.9052D0,  A12 = 75.2293D0,
     &          A03 = 0.4165D0,  A13 = 37.0170D0)

      DOUBLE PRECISION RR(NR),RHO(NR),ZETA(NR),NU(NR),
     &                 AXC(NR),VORT(NR),ETC(NR),DER(NR)

      PADE0(X) = 0.5D0*A00*(1.D0 + A01*X)/
     &            (1.D0 + X*(A01 + X*(A02  + X*A03)))
      PADE1(X) = 0.5D0*A10*(1.D0 + A11*X)/
     &            (1.D0 + X*(A11 + X*(A12  + X*A13)))
      POL(X) = ( (1.D0+X)**1.5D0 + (1.D0-X)**1.5D0 - 2.D0 )
     &         /(2.D0*(DSQRT(2.D0) - 1.D0))

      DO 5 I=1,NR
5        NU(I) = DABS(2.D0*PI*RHO(I)/(VORT(I)*DSQRT(E0)))

      DO 10 I=1,NR 
         RS = 1.D0/DSQRT(PI*RHO(I))
         AA = -DSQRT(8.D0)/(3.D0*PI*RS)
         X = DSQRT(RS)
         ETC(I) = AA*((1.D0+ZETA(I))**1.5D0 + (1.D0-ZETA(I))**1.5D0) 
     &         + PADE0(X) + POL(ZETA(I))*(PADE1(X)-PADE0(X))
10    CONTINUE

      DO 50 I=1,NR
         AXC(I) = 0.782133D0*(DSQRT(2.D0*PI*RHO(I))/VORT(I))*
     &             (-0.15614D0*NU(I)**0.74D0 + 0.0204D0*NU(I)**1.7D0)
     &            /(1.D0 + NU(I)**4)
     &          + 4.D0*(NU(I)**4/VORT(I)) 
     &            *( -0.782133D0*DSQRT(2.D0*PI*RHO(I))*
     &             (1.D0 - 0.211D0*NU(I)**0.74D0 + 0.012D0*NU(I)**1.7D0) 
     &             - ETC(I) ) /(1.D0 + NU(I)**4)**2
50    CONTINUE

      DO 60 I=1,NR
60       AXC(I) = RHO(I)*AXC(I)

      CALL DERIVATIVE(RR,AXC,DER)

      DO 100 I=1,NR
100      AXC(I) = -DER(I)/(RHO(I)*DSQRT(E0))

      RETURN
      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE  DERIVATIVE(RR,F,DF)
      IMPLICIT NONE
  
      INTEGER I
      INCLUDE 'grid.inc'

      DOUBLE PRECISION RR(NR),F(NR),DF(NR)

      DF(1) = (-25.D0*F(1) + 48.D0*F(2) - 36.D0*F(3) + 16.D0*F(4)
     &         - 3.D0*F(5) )/(12.D0*DR*RR(1))

      DF(2) = (- 3.D0*F(1) - 10.D0*F(2) + 18.D0*F(3) - 6.D0*F(4)
     &         + F(5) )/(12.D0*DR*RR(2))

      DO 10 I=3,NR-2
10       DF(I) = (F(I-2) - 8.D0*(F(I-1)-F(I+1)) - F(I+2) )
     &          /(12.D0*DR*RR(I))

      DF(NR-1) = (-F(NR-4) + 6.D0*F(NR-3) - 18.D0*F(NR-2)
     &          + 10.D0*F(NR-1) + 3.D0*F(NR) )/(12.D0*DR*RR(NR-1))

      DF(NR) = (3.D0*F(NR-4) - 16.D0*F(NR-3) + 36.D0*F(NR-2) 
     &         - 48.D0*F(NR-1) + 25.D0*F(NR) )/(12.D0*DR*RR(NR))

      RETURN
      END
C***********************************************************************
C***********************************************************************
C*********************** Begin of library subroutines ******************
C*************** (numerical integration, eigenvalue solver) ************
C***********************************************************************
C***********************************************************************
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module DAVINT from package SLATEC.
* Retrieved from CAMSUN on Thu Nov 20 10:58:05 1997.
* ======================================================================
*DECK DAVINT
      SUBROUTINE DAVINT (X, Y, N, XLO, XUP, ANS, IERR)
C***BEGIN PROLOGUE  DAVINT
C***PURPOSE  Integrate a function tabulated at arbitrarily spaced
C            abscissas using overlapping parabolas.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1B2
C***TYPE      DOUBLE PRECISION (AVINT-S, DAVINT-D)
C***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         DAVINT integrates a function tabulated at arbitrarily spaced
C         abscissas.  The limits of integration need not coincide
C         with the tabulated abscissas.
C
C         A method of overlapping parabolas fitted to the data is used
C         provided that there are at least 3 abscissas between the
C         limits of integration.  DAVINT also handles two special cases.
C         If the limits of integration are equal, DAVINT returns a
C         result of zero regardless of the number of tabulated values.
C         If there are only two function values, DAVINT uses the
C         trapezoid rule.
C
C     Description of Parameters
C         The user must dimension all arrays appearing in the call list
C              X(N), Y(N)
C
C         Input--
C      X    - DOUBLE PRECISION array of abscissas, which must be in
C             increasing order.
C      Y    - DOUBLE PRECISION array of function values. i.e.,
C                Y(I)=FUNC(X(I))
C      N    - The integer number of function values supplied.
C                N .GE. 2 unless XLO = XUP.
C      XLO  - DOUBLE PRECISION lower limit of integration
C      XUP  - DOUBLE PRECISION upper limit of integration.  Must have
C              XLO.LE.XUP
C
C         Output--
C      ANS  - Double Precision computed approximate value of integral
C      IERR - A status code
C           --Normal Code
C                =1 Means the requested integration was performed.
C           --Abnormal Codes
C                =2 Means XUP was less than XLO.
C                =3 Means the number of X(I) between XLO and XUP
C                   (inclusive) was less than 3 and neither of the two
C                   special cases described in the abstract occurred.
C                   No integration was performed.
C                =4 Means the restriction X(I+1).GT.X(I) was violated.
C                =5 Means the number N of function values was .lt. 2.
C                   ANS is set to zero if IERR=2,3,4,or 5.
C
C    DAVINT is documented completely in SC-M-69-335
C    Original program from *Numerical Integration* by Davis & Rabinowitz
C    Adaptation and modifications by Rondall E Jones.
C
C***REFERENCES  R. E. Jones, Approximate integrator of functions
C                 tabulated at arbitrarily spaced abscissas,
C                 Report SC-M-69-335, Sandia Laboratories, 1969.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   690901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAVINT
C
      INTEGER I, IERR, INLFT, INRT, ISTART, ISTOP, N
      DOUBLE PRECISION A, ANS, B, C, CA, CB, CC, FL, FR, R3, RP5,
     1     SLOPE, SUM, SYL, SYL2, SYL3, SYU, SYU2, SYU3, TERM1, TERM2,
     2     TERM3, X, X1, X12, X13, X2, X23, X3, XLO, XUP, Y
      DIMENSION X(*),Y(*)
C     BEGIN BLOCK PERMITTING ...EXITS TO 190
C        BEGIN BLOCK PERMITTING ...EXITS TO 180
C***FIRST EXECUTABLE STATEMENT  DAVINT
            IERR = 1
            ANS = 0.0D0
            IF (XLO .GT. XUP) GO TO 160
               IF (XLO .EQ. XUP) GO TO 150
                  IF (N .GE. 2) GO TO 10
                     IERR = 5
                     CALL XERMSG ('SLATEC', 'DAVINT',
     +                  'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.',
     +                  4, 1)
C     ...............EXIT
                     GO TO 190
   10             CONTINUE
                  DO 20 I = 2, N
C        ............EXIT
                     IF (X(I) .LE. X(I-1)) GO TO 180
C                 ...EXIT
                     IF (X(I) .GT. XUP) GO TO 30
   20             CONTINUE
   30             CONTINUE
                  IF (N .GE. 3) GO TO 40
C
C                    SPECIAL N=2 CASE
                     SLOPE = (Y(2) - Y(1))/(X(2) - X(1))
                     FL = Y(1) + SLOPE*(XLO - X(1))
                     FR = Y(2) + SLOPE*(XUP - X(2))
                     ANS = 0.5D0*(FL + FR)*(XUP - XLO)
C     ...............EXIT
                     GO TO 190
   40             CONTINUE
                  IF (X(N-2) .GE. XLO) GO TO 50
                     IERR = 3
                     CALL XERMSG ('SLATEC', 'DAVINT',
     +                  'THERE WERE LESS THAN THREE FUNCTION VALUES ' //
     +                  'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
C     ...............EXIT
                     GO TO 190
   50             CONTINUE
                  IF (X(3) .LE. XUP) GO TO 60
                     IERR = 3
                     CALL XERMSG ('SLATEC', 'DAVINT',
     +                  'THERE WERE LESS THAN THREE FUNCTION VALUES ' //
     +                  'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
C     ...............EXIT
                     GO TO 190
   60             CONTINUE
                  I = 1
   70             IF (X(I) .GE. XLO) GO TO 80
                     I = I + 1
                  GO TO 70
   80             CONTINUE
                  INLFT = I
                  I = N
   90             IF (X(I) .LE. XUP) GO TO 100
                     I = I - 1
                  GO TO 90
  100             CONTINUE
                  INRT = I
                  IF ((INRT - INLFT) .GE. 2) GO TO 110
                     IERR = 3
                     CALL XERMSG ('SLATEC', 'DAVINT',
     +                  'THERE WERE LESS THAN THREE FUNCTION VALUES ' //
     +                  'BETWEEN THE LIMITS OF INTEGRATION.', 4, 1)
C     ...............EXIT
                     GO TO 190
  110             CONTINUE
                  ISTART = INLFT
                  IF (INLFT .EQ. 1) ISTART = 2
                  ISTOP = INRT
                  IF (INRT .EQ. N) ISTOP = N - 1
C
                  R3 = 3.0D0
                  RP5 = 0.5D0
                  SUM = 0.0D0
                  SYL = XLO
                  SYL2 = SYL*SYL
                  SYL3 = SYL2*SYL
C
                  DO 140 I = ISTART, ISTOP
                     X1 = X(I-1)
                     X2 = X(I)
                     X3 = X(I+1)
                     X12 = X1 - X2
                     X13 = X1 - X3
                     X23 = X2 - X3
                     TERM1 = Y(I-1)/(X12*X13)
                     TERM2 = -Y(I)/(X12*X23)
                     TERM3 = Y(I+1)/(X13*X23)
                     A = TERM1 + TERM2 + TERM3
                     B = -(X2 + X3)*TERM1 - (X1 + X3)*TERM2
     1                   - (X1 + X2)*TERM3
                     C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
                     IF (I .GT. ISTART) GO TO 120
                        CA = A
                        CB = B
                        CC = C
                     GO TO 130
  120                CONTINUE
                        CA = 0.5D0*(A + CA)
                        CB = 0.5D0*(B + CB)
                        CC = 0.5D0*(C + CC)
  130                CONTINUE
                     SYU = X2
                     SYU2 = SYU*SYU
                     SYU3 = SYU2*SYU
                     SUM = SUM + CA*(SYU3 - SYL3)/R3
     1                     + CB*RP5*(SYU2 - SYL2) + CC*(SYU - SYL)
                     CA = A
                     CB = B
                     CC = C
                     SYL = SYU
                     SYL2 = SYU2
                     SYL3 = SYU3
  140             CONTINUE
                  SYU = XUP
                  ANS = SUM + CA*(SYU**3 - SYL3)/R3
     1                  + CB*RP5*(SYU**2 - SYL2) + CC*(SYU - SYL)
  150          CONTINUE
            GO TO 170
  160       CONTINUE
               IERR = 2
               CALL XERMSG ('SLATEC', 'DAVINT',
     +            'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER ' //
     +            'THAN THE LOWER LIMIT.', 4, 1)
  170       CONTINUE
C     ......EXIT
            GO TO 190
  180    CONTINUE
         IERR = 4
         CALL XERMSG ('SLATEC', 'DAVINT',
     +      'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE ' //
     +      'X(I-1) .LT. X(I) FOR ALL I.', 4, 1)
  190 CONTINUE
      RETURN
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END
*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END
*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END
*DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
C***BEGIN PROLOGUE  XGETUA
C***PURPOSE  Return unit number(s) to which error messages are being
C            sent.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XGETUA-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  J4SAVE
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***PURPOSE  Symbolic dump (should be locally written).
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (FDUMP-A)
C***KEYWORDS  ERROR, XERMSG
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
*DECK J4SAVE
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
C***BEGIN PROLOGUE  J4SAVE
C***SUBSIDIARY
C***PURPOSE  Save or recall global variables needed by error
C            handling routines.
C***LIBRARY   SLATEC (XERROR)
C***TYPE      INTEGER (J4SAVE-I)
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                = 6 Refers to the 2nd unit for error messages
C                = 7 Refers to the 3rd unit for error messages
C                = 8 Refers to the 4th unit for error messages
C                = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C***SEE ALSO  XERMSG
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900205  Minor modifications to prologue.  (WRB)
C   900402  Added TYPE section.  (WRB)
C   910411  Added KEYWORDS section.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END
*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module DQAGP from package CMLIB.
* Retrieved from CAMSUN on Mon Jul 20 17:11:05 1998.
* ======================================================================
      SUBROUTINE DQAGP(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,RESULT,ABSERR,
     1   NEVAL,IER,LENIW,LENW,LAST,IWORK,WORK)
C***BEGIN PROLOGUE  DQAGP
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A2A1
C***KEYWORDS  AUTOMATIC INTEGRATOR,EXTRAPOLATION,GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE,SINGULARITIES AT USER SPECIFIED POINTS
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            break points of the integration interval, where local
C            difficulties of the integrand may occur (e.g.
C            SINGULARITIES, DISCONTINUITIES), are provided by the user.
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            NPTS2  - Integer
C                     Number equal to two more than the number of
C                     user-supplied break points within the integration
C                     range, NPTS.GE.2.
C                     If NPTS2.LT.2, The routine will end with IER = 6.
C
C            POINTS - Double precision
C                     Vector of dimension NPTS2, the first (NPTS2-2)
C                     elements of which are the user provided break
C                     points. If these points do not constitute an
C                     ascending sequence there will be an automatic
C                     sorting.
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. it is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. one can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (i.e. SINGULARITY,
C                             DISCONTINUITY within the interval), it
C                             should be supplied to the routine as an
C                             element of the vector points. If necessary
C                             an appropriate special-purpose integrator
C                             must be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             roundoff error is detected in the
C                             extrapolation table.
C                             It is presumed that the requested
C                             tolerance cannot be achieved, and that
C                             the returned RESULT is the best which
C                             can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. it must be noted that
C                             divergence can occur with any other value
C                             of IER.GT.0.
C                         = 6 The input is invalid because
C                             NPTS2.LT.2 or
C                             break points are specified outside
C                             the integration range or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero. Exept when LENIW or LENW or NPTS2 is
C                             invalid, IWORK(1), IWORK(LIMIT+1),
C                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1)
C                             are set to zero.
C                             WORK(1) is set to A anD WORK(LIMIT+1)
C                             to B (where LIMIT = (LENIW-NPTS2)/2).
C
C         DIMENSIONING PARAMETERS
C            LENIW - Integer
C                    Dimensioning parameter for IWORK
C                    LENIW determines LIMIT = (LENIW-NPTS2)/2,
C                    which is the maximum number of subintervals in the
C                    partition of the given integration interval (A,B),
C                    LENIW.GE.(3*NPTS2-2).
C                    If LENIW.LT.(3*NPTS2-2), the routine will end with
C                    IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for WORK
C                    LENW must be at least LENIW*2-NPTS2.
C                    If LENW.LT.LENIW*2-NPTS2, the routine will end
C                    with IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least LENIW. on return,
C                    the first K elements of which contain
C                    pointers to the error estimates over the
C                    subintervals, such that WORK(LIMIT*3+IWORK(1)),...,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
C                    K = LIMIT+1-LAST otherwise
C                    IWORK(LIMIT+1), ...,IWORK(LIMIT+LAST) Contain the
C                     subdivision levels of the subintervals, i.e.
C                     if (AA,BB) is a subinterval of (P1,P2)
C                     where P1 as well as P2 is a user-provided
C                     break point or integration LIMIT, then (AA,BB) has
C                     level L if ABS(BB-AA) = ABS(P2-P1)*2**(-L),
C                    IWORK(LIMIT*2+1), ..., IWORK(LIMIT*2+NPTS2) have
C                     no significance for the user,
C                    note that LIMIT = (LENIW-NPTS2)/2.
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
C                     the right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
C                     contain the corresponding error estimates,
C                    WORK(LIMIT*4+1), ..., WORK(LIMIT*4+NPTS2)
C                     contain the integration limits and the
C                     break points sorted in an ascending sequence.
C                    note that LIMIT = (LENIW-NPTS2)/2.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGPE,XERROR
C***END PROLOGUE  DQAGP
C
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,POINTS,RESULT,WORK
      INTEGER IER,IWORK,LAST,LENIW,LENW,LIMIT,LVL,L1,L2,L3,L4,NEVAL,
     1  NPTS2
C
      DIMENSION IWORK(LENIW),POINTS(NPTS2),WORK(LENW)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGP
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LENIW.LT.(3*NPTS2-2).OR.LENW.LT.(LENIW*2-NPTS2).OR.NPTS2.LT.2)
     1  GO TO 10
C
C         PREPARE CALL FOR DQAGPE.
C
      LIMIT = (LENIW-NPTS2)/2
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
      L4 = LIMIT+L3
C
      CALL DQAGPE(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1  NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),WORK(L4),
     2  IWORK(1),IWORK(L1),IWORK(L2),LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
      LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF(IER.NE.0) CALL XERROR( 'ABNORMAL RETURN FROM DQAGP',26,IER,LVL)
      RETURN
      END
      SUBROUTINE DQAGPE(F,A,B,NPTS2,POINTS,EPSABS,EPSREL,LIMIT,RESULT,
     1   ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,PTS,IORD,LEVEL,NDIN,
     2   LAST)
C***BEGIN PROLOGUE  DQAGPE
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  840425   (YYMMDD)
C***CATEGORY NO.  H2A2A1
C***KEYWORDS  AUTOMATIC INTEGRATOR,EXTRAPOLATION,GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE.,SINGULARITIES AT USER SPECIFIED POINTS
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = Integral of F over (A,B), hopefully
C            satisfying following claim for accuracy ABS(I-RESULT).LE.
C            MAX(EPSABS,EPSREL*ABS(I)). Break points of the integration
C            interval, where local difficulties of the integrand may
C            occur(e.g. singularities,discontinuities),provided by user.
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Double precision version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Double precision
C                     Lower limit of integration
C
C            B      - Double precision
C                     Upper limit of integration
C
C            NPTS2  - Integer
C                     Number equal to two more than the number of
C                     user-supplied break points within the integration
C                     range, NPTS2.GE.2.
C                     If NPTS2.LT.2, the routine will end with IER = 6.
C
C            POINTS - Double precision
C                     Vector of dimension NPTS2, the first (NPTS2-2)
C                     elements of which are the user provided break
C                     POINTS. If these POINTS do not constitute an
C                     ascending sequence there will be an automatic
C                     sorting.
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.NPTS2
C                     If LIMIT.LT.NPTS2, the routine will end with
C                     IER = 6.
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine.
C                             The estimates for integral and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (i.e. SINGULARITY,
C                             DISCONTINUITY within the interval), it
C                             should be supplied to the routine as an
C                             element of the vector points. If necessary
C                             an appropriate special-purpose integrator
C                             must be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             At some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table. It is presumed that
C                             the requested tolerance cannot be
C                             achieved, and that the returned result is
C                             the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.GT.0.
C                         = 6 The input is invalid because
C                             NPTS2.LT.2 or
C                             Break points are specified outside
C                             the integration range or
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             or LIMIT.LT.NPTS2.
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             and ELIST(1) are set to zero. ALIST(1) and
C                             BLIST(1) are set to A and B respectively.
C
C            ALIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left end points
C                     of the subintervals in the partition of the given
C                     integration range (A,B)
C
C            BLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right end points
C                     of the subintervals in the partition of the given
C                     integration range (A,B)
C
C            RLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            PTS    - Double precision
C                     Vector of dimension at least NPTS2, containing the
C                     integration limits and the break points of the
C                     interval in ascending sequence.
C
C            LEVEL  - Integer
C                     Vector of dimension at least LIMIT, containing the
C                     subdivision levels of the subinterval, i.e. if
C                     (AA,BB) is a subinterval of (P1,P2) where P1 as
C                     well as P2 is a user-provided break point or
C                     integration limit, then (AA,BB) has level L if
C                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
C
C            NDIN   - Integer
C                     Vector of dimension at least NPTS2, after first
C                     integration over the intervals (PTS(I)),PTS(I+1),
C                     I = 0,1, ..., NPTS2-2, the error estimates over
C                     some of the intervals may have been increased
C                     artificially, in order to put their subdivision
C                     forward. If this happens for the subinterval
C                     numbered K, NDIN(K) is put to 1, otherwise
C                     NDIN(K) = 0.
C
C            IORD   - Integer
C                     Vector of dimension at least LIMIT, the first K
C                     elements of which are pointers to the
C                     error estimates over the subintervals,
C                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
C                     form a decreasing sequence, with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise
C
C            LAST   - Integer
C                     Number of subintervals actually produced in the
C                     subdivisions process
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DQELG,DQK21,DQPSRT
C***END PROLOGUE  DQAGPE
      DOUBLE PRECISION A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,B,BLIST,B1,B2,CORREC,DABS,DEFABS,DEFAB1,DEFAB2,DMAX1,DMIN1,
     2  DRES,D1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND,
     3  ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW,POINTS,PTS,
     4  RESA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SIGN,TEMP,UFLOW
      INTEGER I,ID,IER,IERRO,IND1,IND2,IORD,IP1,IROFF1,IROFF2,IROFF3,J,
     1  JLOW,JUPBND,K,KSGN,KTMIN,LAST,LEVCUR,LEVEL,LEVMAX,LIMIT,MAXERR,
     2  NDIN,NEVAL,NINT,NINTP1,NPTS,NPTS2,NRES,NRMAX,NUMRL2
      LOGICAL EXTRAP,NOEXT
C
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     1  LEVEL(LIMIT),NDIN(NPTS2),POINTS(NPTS2),PTS(NPTS2),RES3LA(3),
     2  RLIST(LIMIT),RLIST2(52)
C
      EXTERNAL F
C
C            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
C            (LIMEXP+2) AT LEAST).
C
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
C                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
C                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
C           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
C                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
C                       BEEN OBTAINED, IT IS PUT IN RLIST2(NUMRL2) AFTER
C                       NUMRL2 HAS BEEN INCREASED BY ONE.
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
C                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
C                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
C                       TRY TO DECREASE THE VALUE OF ERLARG.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
C                       NO LONGER ALLOWED (TRUE-VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGPE
      EPMACH = D1MACH(4)
C
C            TEST ON VALIDITY OF PARAMETERS
C            -----------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      LEVEL(1) = 0
      NPTS = NPTS2-2
      IF(NPTS2.LT.2.OR.LIMIT.LE.NPTS.OR.(EPSABS.LE.0.0D+00.AND.
     1  EPSREL.LT.DMAX1(0.5D+02*EPMACH,0.5D-28))) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
C            ASCENDING SEQUENCE.
C
      SIGN = 1.0D+00
      IF(A.GT.B) SIGN = -1.0D+00
      PTS(1) = DMIN1(A,B)
      IF(NPTS.EQ.0) GO TO 15
      DO 10 I = 1,NPTS
        PTS(I+1) = POINTS(I)
   10 CONTINUE
   15 PTS(NPTS+2) = DMAX1(A,B)
      NINT = NPTS+1
      A1 = PTS(1)
      IF(NPTS.EQ.0) GO TO 40
      NINTP1 = NINT+1
      DO 20 I = 1,NINT
        IP1 = I+1
        DO 20 J = IP1,NINTP1
          IF(PTS(I).LE.PTS(J)) GO TO 20
          TEMP = PTS(I)
          PTS(I) = PTS(J)
          PTS(J) = TEMP
   20 CONTINUE
      IF(PTS(1).NE.DMIN1(A,B).OR.PTS(NINTP1).NE.DMAX1(A,B)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.
C            ------------------------------------------------
C
   40 RESABS = 0.0D+00
      DO 50 I = 1,NINT
        B1 = PTS(I+1)
        CALL DQK21(F,A1,B1,AREA1,ERROR1,DEFABS,RESA)
        ABSERR = ABSERR+ERROR1
        RESULT = RESULT+AREA1
        NDIN(I) = 0
        IF(ERROR1.EQ.RESA.AND.ERROR1.NE.0.0D+00) NDIN(I) = 1
        RESABS = RESABS+DEFABS
        LEVEL(I) = 0
        ELIST(I) = ERROR1
        ALIST(I) = A1
        BLIST(I) = B1
        RLIST(I) = AREA1
        IORD(I) = I
        A1 = B1
   50 CONTINUE
      ERRSUM = 0.0D+00
      DO 55 I = 1,NINT
        IF(NDIN(I).EQ.1) ELIST(I) = ABSERR
        ERRSUM = ERRSUM+ELIST(I)
   55 CONTINUE
C
C           TEST ON ACCURACY.
C
      LAST = NINT
      NEVAL = 21*NINT
      DRES = DABS(RESULT)
      ERRBND = DMAX1(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.0.1D+03*EPMACH*RESABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(NINT.EQ.1) GO TO 80
      DO 70 I = 1,NPTS
        JLOW = I+1
        IND1 = IORD(I)
        DO 60 J = JLOW,NINT
          IND2 = IORD(J)
          IF(ELIST(IND1).GT.ELIST(IND2)) GO TO 60
          IND1 = IND2
          K = J
   60   CONTINUE
        IF(IND1.EQ.IORD(I)) GO TO 70
        IORD(K) = IORD(I)
        IORD(I) = IND1
   70 CONTINUE
      IF(LIMIT.LT.NPTS2) IER = 1
   80 IF(IER.NE.0.OR.ABSERR.LE.ERRBND) GO TO 210
C
C           INITIALIZATION
C           --------------
C
      RLIST2(1) = RESULT
      MAXERR = IORD(1)
      ERRMAX = ELIST(MAXERR)
      AREA = RESULT
      NRMAX = 1
      NRES = 0
      NUMRL2 = 1
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      ERLARG = ERRSUM
      ERTEST = ERRBND
      LEVMAX = 1
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      IERRO = 0
      UFLOW = D1MACH(1)
      OFLOW = D1MACH(2)
      ABSERR = OFLOW
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*RESABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 160 LAST = NPTS2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
C           ESTIMATE.
C
        LEVCUR = LEVEL(MAXERR)+1
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK21(F,A1,B1,AREA1,ERROR1,RESA,DEFAB1)
        CALL DQK21(F,A2,B2,AREA2,ERROR2,RESA,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+42
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 95
        IF(DABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*DABS(AREA12)
     1  .OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 90
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   90   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   95   LEVEL(MAXERR) = LEVCUR
        LEVEL(LAST) = LEVCUR
        RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE
C
        IF(DMAX1(DABS(A1),DABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*
     1  (DABS(A2)+0.1D+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 100
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 110
  100   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
  110   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(ERRSUM.LE.ERRBND) GO TO 190
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0) GO TO 170
        IF(NOEXT) GO TO 160
        ERLARG = ERLARG-ERLAST
        IF(LEVCUR+1.LE.LEVMAX) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 120
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
        EXTRAP = .TRUE.
        NRMAX = 2
  120   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 140
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER
C           THE LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 130 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
C ***JUMP OUT OF DO-LOOP
          IF(LEVEL(MAXERR)+1.LE.LEVMAX) GO TO 160
          NRMAX = NRMAX+1
  130   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
  140   NUMRL2 = NUMRL2+1
        RLIST2(NUMRL2) = AREA
        IF(NUMRL2.LE.2) GO TO 155
        CALL DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
        KTMIN = KTMIN+1
        IF(KTMIN.GT.5.AND.ABSERR.LT.0.1D-02*ERRSUM) IER = 5
        IF(ABSEPS.GE.ABSERR) GO TO 150
        KTMIN = 0
        ABSERR = ABSEPS
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = DMAX1(EPSABS,EPSREL*DABS(RESEPS))
C ***JUMP OUT OF DO-LOOP
        IF(ABSERR.LT.ERTEST) GO TO 170
C
C           PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
  150   IF(NUMRL2.EQ.1) NOEXT = .TRUE.
        IF(IER.GE.5) GO TO 170
  155   MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        LEVMAX = LEVMAX+1
        ERLARG = ERRSUM
  160 CONTINUE
C
C           SET THE FINAL RESULT.
C           ---------------------
C
C
  170 IF(ABSERR.EQ.OFLOW) GO TO 190
      IF((IER+IERRO).EQ.0) GO TO 180
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00)GO TO 175
      IF(ABSERR.GT.ERRSUM)GO TO 190
      IF(AREA.EQ.0.0D+00) GO TO 210
      GO TO 180
  175 IF(ABSERR/DABS(RESULT).GT.ERRSUM/DABS(AREA))GO TO 190
C
C           TEST ON DIVERGENCE.
C
  180 IF(KSGN.EQ.(-1).AND.DMAX1(DABS(RESULT),DABS(AREA)).LE.
     1  RESABS*0.1D-01) GO TO 210
      IF(0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03.OR.
     1  ERRSUM.GT.DABS(AREA)) IER = 6
      GO TO 210
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  190 RESULT = 0.0D+00
      DO 200 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  200 CONTINUE
      ABSERR = ERRSUM
  210 IF(IER.GT.2) IER = IER-1
      RESULT = RESULT*SIGN
  999 RETURN
      END
      SUBROUTINE DQELG(N,EPSTAB,RESULT,ABSERR,RES3LA,NRES)
C***BEGIN PROLOGUE  DQELG
C***REFER TO  DQAGIE,DQAGOE,DQAGPE,DQAGSE
C***ROUTINES CALLED  D1MACH
C***REVISION DATE  830518   (YYMMDD)
C***KEYWORDS  CONVERGENCE ACCELERATION,EPSILON ALGORITHM,EXTRAPOLATION
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine determines the limit of a given sequence of
C            approximations, by means of the Epsilon algorithm of
C            P.Wynn. An estimate of the absolute error is also given.
C            The condensed Epsilon table is computed. Only those
C            elements needed for the computation of the next diagonal
C            are preserved.
C***DESCRIPTION
C
C           Epsilon algorithm
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C              N      - Integer
C                       EPSTAB(N) contains the new element in the
C                       first column of the epsilon table.
C
C              EPSTAB - Double precision
C                       Vector of dimension 52 containing the elements
C                       of the two lower diagonals of the triangular
C                       epsilon table. The elements are numbered
C                       starting at the right-hand corner of the
C                       triangle.
C
C              RESULT - Double precision
C                       Resulting approximation to the integral
C
C              ABSERR - Double precision
C                       Estimate of the absolute error computed from
C                       RESULT and the 3 previous results
C
C              RES3LA - Double precision
C                       Vector of dimension 3 containing the last 3
C                       results
C
C              NRES   - Integer
C                       Number of calls to the routine
C                       (should be zero at first call)
C***END PROLOGUE  DQELG
C
      DOUBLE PRECISION ABSERR,DABS,DELTA1,DELTA2,DELTA3,DMAX1,D1MACH,
     1  EPMACH,EPSINF,EPSTAB,ERROR,ERR1,ERR2,ERR3,E0,E1,E1ABS,E2,E3,
     2  OFLOW,RES,RESULT,RES3LA,SS,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
      DIMENSION EPSTAB(52),RES3LA(3)
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
C           E1       ELEMENT IN THE EPSILON TABLE IS BASED
C           E2
C           E3                 E0
C                        E3    E1    NEW
C                              E2
C           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
C                    DIAGONAL
C           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
C           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
C                    OF ERROR
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
C           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
C           DIAGONAL OF THE EPSILON TABLE IS DELETED.
C
C***FIRST EXECUTABLE STATEMENT  DQELG
      EPMACH = D1MACH(4)
      OFLOW = D1MACH(2)
      NRES = NRES+1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF(N.LT.3) GO TO 100
      LIMEXP = 50
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 40 I = 1,NEWELM
        K2 = K1-1
        K3 = K1-2
        RES = EPSTAB(K1+2)
        E0 = EPSTAB(K3)
        E1 = EPSTAB(K2)
        E2 = RES
        E1ABS = DABS(E1)
        DELTA2 = E2-E1
        ERR2 = DABS(DELTA2)
        TOL2 = DMAX1(DABS(E2),E1ABS)*EPMACH
        DELTA3 = E1-E0
        ERR3 = DABS(DELTA3)
        TOL3 = DMAX1(E1ABS,DABS(E0))*EPMACH
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
C
C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
C           ACCURACY, CONVERGENCE IS ASSUMED.
C           RESULT = E2
C           ABSERR = ABS(E1-E0)+ABS(E2-E1)
C
        RESULT = RES
        ABSERR = ERR2+ERR3
C ***JUMP OUT OF DO-LOOP
        GO TO 100
   10   E3 = EPSTAB(K1)
        EPSTAB(K1) = E1
        DELTA1 = E1-E3
        ERR1 = DABS(DELTA1)
        TOL1 = DMAX1(E1ABS,DABS(E3))*EPMACH
C
C           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
C           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
C
        IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
        SS = 0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
        EPSINF = DABS(SS*E1)
C
C           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
C           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
C           OF N.
C
        IF(EPSINF.GT.0.1D-03) GO TO 30
   20   N = I+I-1
C ***JUMP OUT OF DO-LOOP
        GO TO 50
C
C           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
C           THE VALUE OF RESULT.
C
   30   RES = E1+0.1D+01/SS
        EPSTAB(K1) = RES
        K1 = K1-2
        ERROR = ERR2+DABS(RES-E2)+ERR3
        IF(ERROR.GT.ABSERR) GO TO 40
        ABSERR = ERROR
        RESULT = RES
   40 CONTINUE
C
C           SHIFT THE TABLE.
C
   50 IF(N.EQ.LIMEXP) N = 2*(LIMEXP/2)-1
      IB = 1
      IF((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM+1
      DO 60 I=1,IE
        IB2 = IB+2
        EPSTAB(IB) = EPSTAB(IB2)
        IB = IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      INDX = NUM-N+1
      DO 70 I = 1,N
        EPSTAB(I)= EPSTAB(INDX)
        INDX = INDX+1
   70 CONTINUE
   80 IF(NRES.GE.4) GO TO 90
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 100
C
C           COMPUTE ERROR ESTIMATE
C
   90 ABSERR = DABS(RESULT-RES3LA(3))+DABS(RESULT-RES3LA(2))
     1  +DABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  100 ABSERR = DMAX1(ABSERR,0.5D+01*EPMACH*DABS(RESULT))
      RETURN
      END
      SUBROUTINE DQK21(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C***BEGIN PROLOGUE  DQK21
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A1A2
C***KEYWORDS  21-POINT GAUSS-KRONROD RULES
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the driver program.
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 21-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 10-POINT GAUSS RULE (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DQK21
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.0666713443 0868813759 3568809893 332 D0 /
      DATA WG  (  2) / 0.1494513491 5058059314 5776339657 697 D0 /
      DATA WG  (  3) / 0.2190863625 1598204399 5534934228 163 D0 /
      DATA WG  (  4) / 0.2692667193 0999635509 1226921569 469 D0 /
      DATA WG  (  5) / 0.2955242247 1475287017 3892994651 338 D0 /
C
      DATA XGK (  1) / 0.9956571630 2580808073 5527280689 003 D0 /
      DATA XGK (  2) / 0.9739065285 1717172007 7964012084 452 D0 /
      DATA XGK (  3) / 0.9301574913 5570822600 1207180059 508 D0 /
      DATA XGK (  4) / 0.8650633666 8898451073 2096688423 493 D0 /
      DATA XGK (  5) / 0.7808177265 8641689706 3717578345 042 D0 /
      DATA XGK (  6) / 0.6794095682 9902440623 4327365114 874 D0 /
      DATA XGK (  7) / 0.5627571346 6860468333 9000099272 694 D0 /
      DATA XGK (  8) / 0.4333953941 2924719079 9265943165 784 D0 /
      DATA XGK (  9) / 0.2943928627 0146019813 1126603103 866 D0 /
      DATA XGK ( 10) / 0.1488743389 8163121088 4826001129 720 D0 /
      DATA XGK ( 11) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0116946388 6737187427 8064396062 192 D0 /
      DATA WGK (  2) / 0.0325581623 0796472747 8818972459 390 D0 /
      DATA WGK (  3) / 0.0547558965 7435199603 1381300244 580 D0 /
      DATA WGK (  4) / 0.0750396748 1091995276 7043140916 190 D0 /
      DATA WGK (  5) / 0.0931254545 8369760553 5065465083 366 D0 /
      DATA WGK (  6) / 0.1093871588 0229764189 9210590325 805 D0 /
      DATA WGK (  7) / 0.1234919762 6206585107 7958109831 074 D0 /
      DATA WGK (  8) / 0.1347092173 1147332592 8054001771 707 D0 /
      DATA WGK (  9) / 0.1427759385 7706008079 7094273138 717 D0 /
      DATA WGK ( 10) / 0.1477391049 0133849137 4841515972 068 D0 /
      DATA WGK ( 11) / 0.1494455540 0291690566 4936468389 821 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK21
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
      FC = F(CENTR)
      RESK = WGK(11)*FC
      RESABS = DABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(11)*DABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*DMIN1(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = DMAX1
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
      SUBROUTINE DQPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C***BEGIN PROLOGUE  DQPSRT
C***REFER TO  DQAGE,DQAGIE,DQAGPE,DQAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION DATE  810101   (YYMMDD)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  This routine maintains the descending ordering in the
C            list of the local error estimated resulting from the
C            interval subdivision process. At each call two error
C            estimates are inserted using the sequential search
C            method, top-down for the largest error estimate and
C            bottom-up for the smallest error estimate.
C***DESCRIPTION
C
C           Ordering routine
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - Integer
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - Integer
C                       Number of error estimates currently in the list
C
C              MAXERR - Integer
C                       Maxerr points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - Double precision
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - Double precision
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - Integer
C                       Vector of dimension LAST, the first K elements
C                       of which contain pointers to the error
C                       estimates, such that
C                       ELIST(IORD(1)),...,  ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - Integer
C                       MAXERR = IORD(NRMAX)
C***END PROLOGUE  DQPSRT
C
      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(LAST),IORD(LAST)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
C           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
C           SUBDIVISIONS STILL ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
C*********************************************************************
C*********************************************************************
C*********************************************************************
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
*  Purpose
*  =======
*
*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of DLAE2
*
      END
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of DLAEV2
*
      END
      DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANST  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric tridiagonal matrix A.
*
*  Description
*  ===========
*
*  DLANST returns the value
*
*     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANST as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
*          set to zero.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of A.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) sub-diagonal or super-diagonal elements of A.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $              ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+
     $                 ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL DLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL DLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      DLANST = ANORM
      RETURN
*
*     End of DLANST
*
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of DLASET
*
      END
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASR applies a sequence of plane rotations to a real matrix A,
*  from either the left or the right.
*  
*  When SIDE = 'L', the transformation takes the form
*  
*     A := P*A
*  
*  and when SIDE = 'R', the transformation takes the form
*  
*     A := A*P**T
*  
*  where P is an orthogonal matrix consisting of a sequence of z plane
*  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
*  and P**T is the transpose of P.
*  
*  When DIRECT = 'F' (Forward sequence), then
*  
*     P = P(z-1) * ... * P(2) * P(1)
*  
*  and when DIRECT = 'B' (Backward sequence), then
*  
*     P = P(1) * P(2) * ... * P(z-1)
*  
*  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
*  
*     R(k) = (  c(k)  s(k) )
*          = ( -s(k)  c(k) ).
*  
*  When PIVOT = 'V' (Variable pivot), the rotation is performed
*  for the plane (k,k+1), i.e., P(k) has the form
*  
*     P(k) = (  1                                            )
*            (       ...                                     )
*            (              1                                )
*            (                   c(k)  s(k)                  )
*            (                  -s(k)  c(k)                  )
*            (                                1              )
*            (                                     ...       )
*            (                                            1  )
*  
*  where R(k) appears as a rank-2 modification to the identity matrix in
*  rows and columns k and k+1.
*  
*  When PIVOT = 'T' (Top pivot), the rotation is performed for the
*  plane (1,k+1), so P(k) has the form
*  
*     P(k) = (  c(k)                    s(k)                 )
*            (         1                                     )
*            (              ...                              )
*            (                     1                         )
*            ( -s(k)                    c(k)                 )
*            (                                 1             )
*            (                                      ...      )
*            (                                             1 )
*  
*  where R(k) appears in rows and columns 1 and k+1.
*  
*  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
*  performed for the plane (k,z), giving P(k) the form
*  
*     P(k) = ( 1                                             )
*            (      ...                                      )
*            (             1                                 )
*            (                  c(k)                    s(k) )
*            (                         1                     )
*            (                              ...              )
*            (                                     1         )
*            (                 -s(k)                    c(k) )
*  
*  where R(k) appears in rows and columns k and z.  The rotations are
*  performed without ever forming P(k) explicitly.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P**T
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
*          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C       (input) DOUBLE PRECISION array, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          The cosines c(k) of the plane rotations.
*
*  S       (input) DOUBLE PRECISION array, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          The sines s(k) of the plane rotations.  The 2-by-2 plane
*          rotation part of the matrix P(k), R(k), has the form
*          R(k) = (  c(k)  s(k) )
*                 ( -s(k)  c(k) ).
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The M-by-N matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P**T if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLASR
*
      END
      SUBROUTINE DLASRT( ID, N, D, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * )
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band symmetric matrix can also be found
*  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
*  tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  symmetric matrix.  On entry, Z must contain the
*                  orthogonal matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is orthogonally similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR,
     $                   DLASRT, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      GO TO 190
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
*
  190 CONTINUE
      RETURN
*
*     End of DSTEQR
*
      END
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
     $                   NMAXIT
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
     $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
     $                   SIGMA, SSFMAX, SSFMIN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $       1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 150
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
     $         GO TO 120
  110    CONTINUE
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use DLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 150
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 150
*
      END IF
*
*     Undo scaling if necessary
*
  150 CONTINUE
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.LT.NMAXIT )
     $   GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      GO TO 180
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
*
  180 CONTINUE
      RETURN
*
*     End of DSTERF
*
      END
      SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric tridiagonal matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix A, stored in elements 1 to N-1 of E.
*          On exit, the contents of E are destroyed.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with D(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If JOBZ = 'N', WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of E did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTZ
      INTEGER            IMAX, ISCALE
      DOUBLE PRECISION   BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM,
     $                   TNRM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTEQR, DSTERF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -6
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      TNRM = DLANST( 'M', N, D, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL DSCAL( N, SIGMA, D, 1 )
         CALL DSCAL( N-1, SIGMA, E( 1 ), 1 )
      END IF
*
*     For eigenvalues only, call DSTERF.  For eigenvalues and
*     eigenvectors, call DSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, D, E, INFO )
      ELSE
         CALL DSTEQR( 'I', N, D, E, Z, LDZ, WORK, INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, D, 1 )
      END IF
*
      RETURN
*
*     End of DSTEV
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) INTEGER 
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
 
* BLAS REQUIRED BY LAPACK ROUTINE:  dstev
* -----------------------------------------------------------
* Note: Link to BLAS optimized for your system, if available.
* -----------------------------------------------------------
 
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C    I1MACH( 7) = A, THE BASE.
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C               WHERE  EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, THE BASE.
C  SINGLE-PRECISION
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C  DOUBLE-PRECISION
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
      INTEGER IMACH(16), OUTPUT, SC, SMALL(2)
      SAVE IMACH, SC
      REAL RMACH
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
      INTEGER I3, J, K, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
C  INCLUDING AUTO-DOUBLE COMPILERS.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SC/987/
C
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532
     *           .AND. SMALL(2) .EQ. -448790528)
     *       .OR.     (SMALL(2) .EQ. 1117925532
     *           .AND. SMALL(1) .EQ. -448790528)) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SC = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
            ELSE
               DO 10 I3 = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(I3)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
*              *** CRAY T3E ***
               IMACH( 1) = 5
               IMACH( 2) = 6
               IMACH( 3) = 0
               IMACH( 4) = 0
               IMACH( 5) = 64
               IMACH( 6) = 8
               IMACH( 7) = 2
               IMACH( 8) = 63
               CALL I1MCR1(IMACH(9), K, 32767, 16777215, 16777215)
               IMACH(10) = 2
               IMACH(11) = 53
               IMACH(12) = -1021
               IMACH(13) = 1024
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               GO TO 35
 20            CALL I1MCR1(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 46
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 45
               CALL I1MCR1(IMACH(9), K, 0, 4194303, 16777215)
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8188
               IMACH(13) = 8189
               IMACH(14) = 94
               IMACH(15) = -8141
               IMACH(16) = 8189
               GO TO 35
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
 35      SC = 987
         END IF
 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/
     * ' statements appropriate for your machine and setting'/
     * ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/
     * ' appropriate for your machine.')
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 40
      I1MACH = IMACH(I)
      RETURN
 40   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
      STOP
* /* C source for I1MACH -- remove the * in column 1 */
* /* Note that some values may need changing. */
*#include <stdio.h>
*#include <float.h>
*#include <limits.h>
*#include <math.h>
*
*long i1mach_(long *i)
*{
*	switch(*i){
*	  case 1:  return 5;	/* standard input */
*	  case 2:  return 6;	/* standard output */
*	  case 3:  return 7;	/* standard punch */
*	  case 4:  return 0;	/* standard error */
*	  case 5:  return 32;	/* bits per integer */
*	  case 6:  return sizeof(int);
*	  case 7:  return 2;	/* base for integers */
*	  case 8:  return 31;	/* digits of integer base */
*	  case 9:  return LONG_MAX;
*	  case 10: return FLT_RADIX;
*	  case 11: return FLT_MANT_DIG;
*	  case 12: return FLT_MIN_EXP;
*	  case 13: return FLT_MAX_EXP;
*	  case 14: return DBL_MANT_DIG;
*	  case 15: return DBL_MIN_EXP;
*	  case 16: return DBL_MAX_EXP;
*	  }
*	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
*	exit(1);return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCR1(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
