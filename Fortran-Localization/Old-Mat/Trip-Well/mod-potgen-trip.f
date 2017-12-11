      Module POTGEN3

      IMPLICIT NONE

      contains


      SUBROUTINE POTENTIAL(V,GRID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: GRID

      DOUBLE PRECISION, INTENT(OUT) :: V(GRID)

      INTEGER :: I

      DO I=1,grid
         V(I) = 0.D0
         IF (I.GE.10.AND.I.LE.40) THEN
            V(I) = -6.D0
         ELSEIF (I.GE.60.AND.I.LE.75) THEN
            V(I) = 2.0D0
         ELSEIF (I.GE.80.AND.I.LE.95) THEN
            V(I) = -2.0D0
         ELSEIF (I.GE.125.AND.I.LE.180) THEN
            V(I) = -3.0D0
         ELSE
            V(I) = 0.D0
         ENDIF
      end do

      END SUBROUTINE

      END MODULE

