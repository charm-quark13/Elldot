
      INTEGER I
      DOUBLE PRECISION NORM,FODA

      INCLUDE 'grid.inc'

      DOUBLE PRECISION R(NR),RR(NR),RHO(NR)

C**------------------------------------------------------------------***
C**   set up the logarithmic grid                                    ***
C**------------------------------------------------------------------***
      DO 10 I=1,NR  
10       R(I)= (DBLE(I)-N0)*DR   

      DO 11 I=1,NR
11       RR(I) = DEXP(R(I))

      DO 12 I=1,NR
12       WRITE(10,*)I,R(I),RR(I)

      DO 13 I=1,NR 
         WRITE(13,*)RR(I),real(FODA(0,1,RR(I)))
13       RHO(I) = RR(I)**2 * FODA(0,1,RR(I))**2

      CALL NORMCALC(RHO,NORM)
      WRITE(*,*)NORM

      END
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE NORMCALC(RHO,NORM)
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION NORM,CC
      INCLUDE 'grid.inc'

      DOUBLE PRECISION RHO(NR),CO(NR)

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
20       NORM = NORM + CO(I)*RHO(I)

      NORM = NORM * 2.D0*PI

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
      End
