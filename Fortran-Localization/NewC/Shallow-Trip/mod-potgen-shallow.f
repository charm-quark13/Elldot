      Module POTGEN3

      IMPLICIT NONE

      contains


      SUBROUTINE POTENTIAL(V,GRID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: GRID

      DOUBLE PRECISION, INTENT(OUT) :: V(GRID)

      INTEGER :: I

      DO I=1,GRID
         V(I) = 0.D0
         IF (I.GE.10.AND.I.LE.60) THEN
            V(I) = -1.5D0
!         ELSEIF (I.GE.60.AND.I.LE.75) THEN
!            V(I) = 0.0D0
         ELSEIF (I.GE.80.AND.I.LE.140) THEN
            V(I) = -4.0D0
         ELSEIF (I.GE.155.AND.I.LE.195) THEN
            V(I) = -1.0D0
         ELSE
            V(I) = 0.D0
         ENDIF

      end do

      END SUBROUTINE

      END MODULE

