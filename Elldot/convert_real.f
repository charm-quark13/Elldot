      Program CONVERSION

      INTEGER NUDIM,NDDIM,I,J,K,M,FILES

      INCLUDE 'grid.inc'
      INCLUDE 'mode.inc'

      PARAMETER(NUDIM=MAX(NU,1),NDDIM=MAX(ND,1))
      PARAMETER(M=NR*NLU)

      DOUBLE PRECISION UXOLD(M,NUDIM),DXOLD(M,NDDIM)

      DOUBLE COMPLEX UXCOM(M,NUDIM),DXCOM(M,NDDIM)
      DOUBLE COMPLEX ONE,IONE,ZERO

      PARAMETER (ONE=(1.D0,0.D0),IONE=(0.D0,1.D0),
     &           ZERO=(0.D0,0.D0))

C*****************************************************************
C***    FILES DENOTES THE NUMBER OF FILES BEING READ IN AND    ***
C***       OUTPUTTED (SO THE PROGRAM CAN CLOSE THE FILES)      ***
c***                                                           ***
C*****************************************************************
      PARAMETER( FILES = 6 )

      OPEN(1,FILE='fort.71',STATUS='OLD',FORM='FORMATTED')
      OPEN(2,FILE='fort.72',STATUS='OLD',FORM='FORMATTED')
      OPEN(3,FILE='fort.73',STATUS='OLD',FORM='FORMATTED')

      OPEN(4,FILE='fort.81',STATUS='OLD',FORM='FORMATTED')
      OPEN(5,FILE='fort.82',STATUS='OLD',FORM='FORMATTED')
      OPEN(6,FILE='fort.83',STATUS='OLD',FORM='FORMATTED')

      DO 10 I=1,M
      DO 10 K=1,NUDIM

         READ(K,*) UXOLD(I,K)

10    CONTINUE

      DO 20 I=1,M
      DO 20 K=1,NDDIM

         READ(3+K,*) DXOLD(I,K)

20    CONTINUE
    
      DO 30 I=1,M
      DO 30 K=1,NUDIM

         UXCOM(I,K) = UXOLD(I,K)*ONE

30    CONTINUE

      DO 40 I=1,M
      DO 40 K=1,NDDIM

         DXCOM(I,K) = DXOLD(I,K)*ONE

40    CONTINUE

      DO 50 I=1,M
      DO 50 K=1,NUDIM

         WRITE(40+K,*) UXCOM(I,K)

50    CONTINUE

      DO 60 I=1,M
      DO 60 K=1,NDDIM

         WRITE(50+K,*) DXCOM(I,K) 

60    CONTINUE

      DO 70 I=1,FILES 
        CLOSE(I)
70    CONTINUE

      END
