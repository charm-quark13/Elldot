      PROGRAM MAIN

      IMPLICIT NONE

      INCLUDE 'add.inc'
      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'

      DOUBLE PRECISION TOL,R0,HA0,BB,MUBOHR,G,OMEGA
      INTEGER I,J,K,STEP,MUP,MDN,NUDIM,NDDIM,INFO
      
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
     &               DX(MDN,NDDIM),DX_TI(MDN,NDDIM),DXT(MDN,NDDIM)

      DOUBLE PRECISION T(NT),R(NR),RR(NR),POTU(MUP),POTD(MDN)
      DOUBLE PRECISION INVD,C
      INTEGER IPIVU(MUP),IPIVD(MDN)

      PARAMETER(IONE=(0.D0,1.D0),ONE=(1.D0,0.D0),ZERO=(0.D0,0.D0))

C      OPEN(1,FILE='fort.71',STATUS='OLD',FORM='FORMATTED')
C      OPEN(2,FILE='fort.72',STATUS='OLD',FORM='FORMATTED')
C      OPEN(3,FILE='fort.73',STATUS='OLD',FORM='FORMATTED')

      DO 12 I=1,NT   
         T(I) = 0.D0
12    CONTINUE
      
      INVD = 1.D0/DELTA
C*****===================================================================*****
C**** BEGINNING OF THE TIME ITERATION LOOP                                ****
C*****===================================================================*****

      DO 13 STEP=1,NT
         T(STEP)= T(STEP)+0.1D0*INVD

C*****===================================================================*****
C     SETTING UP PARABOLIC WELL [ f(r)=0.5*C*r**2 ]
C        C0 READ IN (DEFINES INITIAL SHAPE OF WELL) 
C        2ND TERM IN C DEFINES THE OSCILLATION OF WELL IN TIME
C*****===================================================================*****

             C = C0+B*1*DSIN(OMEGAC*T(STEP))

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
                 write(101,*)RR(I),POTD(I)
20           CONTINUE
                 
                 stop
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
                DO 45 I=1,NUDIM
                   READ(I) UX_TI(MUP,I)
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
             CALL CGESV( MUP, NUDIM, AUP, MUP, IPIVU, UXT, MUP, INFO)

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
                DO 145 I=1,NDDIM
                   READ(I) DX_TI(MDN,I)
145              CONTINUE

                DO 146 I=1,MDN
                DO 146 K=1,NDDIM
                   DX(I,K) = DX_TI(I,K)
146              CONTINUE

                DO 147 I=1,MDN
                DO 147 K=1,NDDIM
                   UX(I,K) = UXT(I,K)
147              CONTINUE
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
             CALL CGESV( MDN, NDDIM, ADN, MDN, IPIVD, DXT, MDN, INFO)

             DO 160 I=1,MDN
             DO 160 K=1,NDDIM
                WRITE(200+K,*) DXT(I,K)
160           CONTINUE


13    CONTINUE

      END



******====================================================================*****
******====================================================================*****
******====================================================================*****

      SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  CGESV computes the solution to a complex system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           CGETRF, CGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL CGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL CGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     End of CGESV
*
      END 
