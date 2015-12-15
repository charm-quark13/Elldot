      PROGRAM DOTPOT
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

      PAR=2

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

      WRITE(*,*) VEXT, V0

      RETURN
      END

