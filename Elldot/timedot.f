      PROGRAM MAIN

      IMPLICIT NONE

      INCLUDE 'add.inc'
      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'

      DOUBLE PRECISION TOL,R0,HA0,BB,MUBOHR,G,OMEGA
      INTEGER I,J,K,N,STEP,MUP,MDN,NUDIM,NDDIM,INFO
      
      PARAMETER (NUDIM =MAX(NU,1), NDDIM=MAX(ND,1))
      PARAMETER (MUP=NR*NLU, MDN=NR*NLD)

      DOUBLE COMPLEX HAMUP(MUP,MUP),HAMDN(MDN,MDN),
     &               AUP(MUP,MUP),ADN(MDN,MDN),
     &               ONE,IONE,ZERO,
     &               OFFDIAGU(MUP),OFFDIAGD(MDN),
     &               MAINDIAGU(MUP),MAINDIAGD(MDN),
     &               MATU(MUP,NU),MATD(MDN,ND),
     &               IDENU(MUP,MUP),IDEND(MDN,MDN),
     &               IDIAGU(MUP), IDIAGD(MDN),
     &               UX(MUP,NUDIM),UX_TI(MUP,NUDIM),UXT(MUP,NUDIM),
     &               DX(MDN,NDDIM),DX_TI(MDN,NDDIM),DXT(MDN,NDDIM),
     &               NDX(MDN,MDN),NUX(MUP,MUP),
     &               UXTC(MUP,NUDIM),DXTC(MDN,NDDIM)

      DOUBLE PRECISION T(NT),R(NR),RR(NR),POTU(MUP),POTD(MDN)
      DOUBLE PRECISION INVD,C
      INTEGER IPIVU(MUP),IPIVD(MDN),FILES

      PARAMETER(IONE=(0.D0,1.D0),ONE=(1.D0,0.D0),ZERO=(0.D0,0.D0))
      PARAMETER(FILES = 6)
      
      OPEN(1,FILE='fort.41',STATUS='OLD',FORM='FORMATTED')
      OPEN(2,FILE='fort.42',STATUS='OLD',FORM='FORMATTED')
      OPEN(3,FILE='fort.43',STATUS='OLD',FORM='FORMATTED')

      OPEN(4,FILE='fort.51',STATUS='OLD',FORM='FORMATTED')
      OPEN(5,FILE='fort.52',STATUS='OLD',FORM='FORMATTED')
      OPEN(6,FILE='fort.53',STATUS='OLD',FORM='FORMATTED')


      DO 12 I=1,NT   
         T(I) = 0.D0
12    CONTINUE
      
      INVD = 1.D0/DELTA
C*****===================================================================*****
C**** BEGINNING OF THE TIME ITERATION LOOP                                ****
C*****===================================================================*****

      DO 13 STEP=1,NT
         T(STEP)= T(STEP)+0.1D-3*INVD
         write(*,*) step
C*****===================================================================*****
C     SETTING UP PARABOLIC WELL [ f(r)=0.5*C*r**2 ]
C        C0 READ IN (DEFINES INITIAL SHAPE OF WELL) 
C        2ND TERM IN C DEFINES THE OSCILLATION OF WELL IN TIME
C*****===================================================================*****

C             C = C0+B*1*DSIN(OMEGAC*T(STEP))
              C = C0

C*****===================================================================*****
C     SETTING THE VALUE OF R ACCORDING TO LOG GRID
C*****===================================================================*****
      
             DO 10 I=1,NR
                R(I) = (DBLE(I)-N0)*DR
10           CONTINUE
      
             DO 11 I=1,NR
                RR(I) = DEXP(R(I))
11           CONTINUE
C*****--------------------------------------------------------------------*****
C     CALCULATION OF THE POTENTIAL WELL ACROSS ALL POINTS
C*****--------------------------------------------------------------------*****
             DO 20 I=1,NR
                POTD(I) = 0.5D0*C*RR(I)**2
                POTU(I) = 0.5D0*C*RR(I)**2
20           CONTINUE
                 

C*****====================================================================*****
C***  GENERATING IDENTITY MATRIX                                            ***
C*****====================================================================*****
             DO 30 I=1,MUP
                IDIAGU(I)=ONE
30           CONTINUE

             DO 31 I=1,MUP
             DO 31 K=1,MUP
                IDENU(I,K)=ZERO
31           CONTINUE

             DO 32 I=1,MUP
                IDENU(I,I) = IDIAGU(I)
32           CONTINUE
C*****====================================================================*****
C     SETTING UP HAMILTONIAN                                                 **
C*****====================================================================*****   
             DO 40 I=1,MUP
             DO 40 K=1,MUP
                HAMUP(I,K) = ZERO
40           CONTINUE

             DO 41 J=1,NLU
             DO 41 I=1,NR
                MAINDIAGU((J-1)*NR+I)=ZERO
                OFFDIAGU((J-1)*NR+I)=ZERO
41           CONTINUE

             DO 42 J=1,NLU
             DO 42 I=1,NR
                MAINDIAGU((J-1)*NR+I)= INVD + POTU(I)
                OFFDIAGU((J-1)*NR+I) =(-0.5D0)*INVD
42           CONTINUE
C*****====================================================================*****
C***  COMPILING THE HAMILTONIAN ACCORDING TO CRANK-NICHOLSON                ***
C***     (1+iH*t*0.5)PSI(n+1) = (1-iH*t*0.5)PSI(n)                          ***
C*****====================================================================*****

             DO 43 I=1,MUP
                HAMUP(I,I) = MAINDIAGU(I)*IONE*0.5D0*T(STEP)
43           CONTINUE
      
             DO 44 I=1,MUP-1
                HAMUP(I,I+1) = OFFDIAGU(I)*IONE*0.5D0*T(STEP)
44           CONTINUE

             AUP=IDENU-HAMUP

C*****====================================================================*****
C***   READING IN WAVEFUNCTIONS AND SETTING WAVEFUNCTION TO PREVIOUS        *** 
C***     OUTPUT FROM LAST ITERATION                                         ***
C*****====================================================================*****

             IF (STEP.EQ.1) THEN
                DO 45 K=1,MUP
                   READ(1,*) UX_TI(K,1)
                   READ(2,*) UX_TI(K,2)
                   READ(3,*) UX_TI(K,3) 
45              CONTINUE

                DO 46 I=1,MUP
                DO 46 K=1,NUDIM
                   UX(I,K) = UX_TI(I,K)
46              CONTINUE

             ELSE
                DO 47 I=1,MUP
                DO 47 K=1,NUDIM
                   UX(I,K) = UXT(I,K)
47              CONTINUE
             ENDIF  

             DO 50 I=1,MUP
             DO 50 K=1,NUDIM
                MATU(I,K) = ZERO
50           CONTINUE

             MATU=MATMUL(AUP,UX)             
C*****====================================================================*****
C***   SETTING UP THE LHS OF THE CRANK-NICHOLSON ALGORITHM                  ***
C*****====================================================================*****
             AUP = IDENU+HAMUP             

C*****====================================================================*****
C***   SOLVING THE LINEAR EQUATIONS AND DUMPING RESULTS TO FORT FILES       ***
C*****====================================================================*****
             CALL ZGESV(MUP,NUDIM,AUP,MUP,IPIVU,UXT,MUP,INFO)

             DO 60 I=1,MUP
             DO 60 K=1,NUDIM
                WRITE(100+K,*) UXT(I,K)
60           CONTINUE


*******************************************************************************
*******************************************************************************
C           BEGINNING OF THE SPIN DOWN TIME EVOLUTION                         *
C                                                                             *  
****                                                                       ****
***************     **********************************     ********************
*******************************************************************************
***************     **********************************     ********************


C*****====================================================================*****
C***  GENERATING IDENTITY MATRIX                                            ***
C*****====================================================================*****
             DO 130 I=1,MDN
                IDIAGD(I)=ONE
130           CONTINUE

             DO 131 I=1,MDN
             DO 131 K=1,MDN
                IDEND(I,K)=ZERO
131           CONTINUE

             DO 132 I=1,MDN
                IDEND(I,I) = IDIAGD(I)
132           CONTINUE

C*****====================================================================*****
C     SETTING UP HAMILTONIAN                                                 **
C*****====================================================================*****   
             DO 140 I=1,MDN
             DO 140 K=1,MDN
                HAMDN(I,K) = ZERO
140           CONTINUE

             DO 141 J=1,NLU
             DO 141 I=1,NR
                MAINDIAGD((J-1)*NR+I)=ZERO
                OFFDIAGD((J-1)*NR+I)=ZERO
141           CONTINUE

             DO 142 J=1,NLD
             DO 142 I=1,NR
                MAINDIAGD((J-1)*NR+I)= INVD + POTD(I)
                OFFDIAGD((J-1)*NR+I) =(-0.5D0)*INVD
142           CONTINUE

C*****====================================================================*****
C***  COMPILING THE HAMILTONIAN ACCORDING TO CRANK-NICHOLSON                ***
C***     (1+iH*t*0.5)PSI(n+1) = (1-iH*t*0.5)PSI(n)                          ***
C*****====================================================================*****
             DO 143 I=1,MDN
                HAMDN(I,I) = MAINDIAGD(I)*IONE*0.5D0*T(STEP)
143           CONTINUE

             DO 144 I=1,MDN-1
                HAMDN(I,I+1) = OFFDIAGD(I)*IONE*0.5D0*T(STEP)
144           CONTINUE

             ADN=IDEND-HAMDN

C*****====================================================================*****
C***   READING IN WAVEFUNCTIONS AND SETTING WAVEFUNCTION TO PREVIOUS        *** 
C***     OUTPUT FROM LAST ITERATION                                         ***
C*****====================================================================*****

             IF (STEP.EQ.1) THEN
                 DO 145 K=1,MDN
                   READ(4,*) DX_TI(K,1)
                   READ(5,*) DX_TI(K,2)
                   READ(6,*) DX_TI(K,3)  
145              CONTINUE

                DO 146 I=1,MDN
                DO 146 K=1,NDDIM
                   DX(I,K) = DX_TI(I,K)
146              CONTINUE

             ELSE
                DO 147 I=1,MDN
                DO 147 K=1,NDDIM
                   DX(I,K) = DXT(I,K)
147             CONTINUE
             ENDIF

             DO 150 I=1,MDN
             DO 150 K=1,NDDIM
                MATD(I,K) = ZERO
150           CONTINUE

             MATD=MATMUL(ADN,DX)

C*****====================================================================*****
C***   SETTING UP THE LHS OF THE CRANK-NICHOLSON ALGORITHM                  ***
C*****====================================================================*****
             ADN = IDEND+HAMDN

C*****====================================================================*****
C***   SOLVING THE LINEAR EQUATIONS AND DUMPING RESULTS TO FORT FILES       ***
C*****====================================================================*****
             CALL ZGESV(MDN,NDDIM,ADN,MDN,IPIVD,DXT,MDN,INFO)

             DO 160 I=1,MDN
             DO 160 K=1,NDDIM
                WRITE(200+K,*) DXT(I,K)
160           CONTINUE

C             DO 170 I=1,MDN
C             DO 170 J=1,MUP
C             DO 170 K=1,NDDIM
C             DO 170 N=1,NUDIM
C                DXTC(I,K)=CONJG(DXT(I,K))
C                UXTC(J,N)=CONJG(UXT(J,N))
C170           CONTINUE

C             NUX=MATMUL(UXT,TRANSPOSE(UXTC))
C             NDX=MATMUL(DXT,TRANSPOSE(DXTC))

             WRITE(*,*) 'iteration completed'
   
13    CONTINUE

C      DO 999 I=1,FILES
C         CLOSE(I)
C999   CONTINUE
      
      END

      
