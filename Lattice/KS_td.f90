PROGRAM TWOSPIN
IMPLICIT NONE

INTEGER :: I,INFO,MODE,ITER, xc_slater, corr, t_it, &
            pc_count, L
!integer :: J,k, vec, krev, xc_guess
integer, parameter :: lwork=100, np=3, pc=2

real(8) :: C,CP,T,TP,V(3),BX(3),BY(3),BZ(3), &
      VT(3),BXT(3),BYT(3),BZT(3),VXC(3), &
      vhxc(3),BXCX(3),BXCY(3),BXCZ(3), &
      vhxc0(3),BXCX0(3),BXCY0(3),BXCZ0(3), &
      vhxc1(3),BXCX1(3),BXCY1(3),BXCZ1(3), &
      N(3),MX(3),MY(3),MZ(3),E(6),RWORK(100), &
      dt, time, &
      TT(3),TX(3),TY(3),TZ(3)
!real(8) :: mmat(500, 6, 3), mback(500, 6, 3)
real(8) :: U0,MIX,TOL,EOLD,EH,ETOT,EVXC,EXC,EX,EC
!real(8) :: CRIT

character(100) ::  wmat, wsite
! character(100) :: vec_files(3),

complex(8) :: M(6,6),WORK(LWORK),GAM(2,2,3,3),PHI(6,6), m_temp(6,6)
complex(8), parameter :: one = (1.d0, 0.d0), ione = (0.d0, 1.d0), zero = (0.d0, 0.d0)
COMMON /EVALS/ E,PHI,GAM

write(wmat,'(a)') '(5f15.9)'

dt = .01d0
time = 0.d0
U0 = 1.d0
L = 1

!      T = 0.5D0
!      WRITE(*,*)'C?'
!      READ(*,*)C
!      WRITE(*,*)'Chain or triangle? (1,2)'
!      READ(*,*)MODE
!      IF (MODE.EQ.1) THEN
!         TP = 0.D0
!         CP = 0.D0
!      ELSEIF (MODE.EQ.2) THEN
!         TP = T
!         CP = C
!      ELSE
!         WRITE(*,*)'Invalid choice'
!         STOP
!      ENDIF

! do i = 1, 3
!   v(i) = 0.d0
!   bx(i) = 0.d0
!   by(i) = 0.d0
!   bz(i) = 0.d0
! end do

xc_slater = 1
corr = 0
mode = 2

MIX = 0.1D0
TOL = 1.D-8

open(1, file='planar_field.txt')


DO 5 I=1,3
      READ(1,*)V(I),BX(I),BY(I),BZ(I)
      N(I) = 0.D0
      VHXC(I) = 0.D0
      BXCX(I) = 0.D0
      BXCY(I) = 0.D0
      BXCZ(I) = 0.D0
5     CONTINUE

close(1)

open(1, file='constants.txt')

read(1,*) c, t

U0 = 1.d0

TP = T

cp = c

E(1) = 1.d0
EOLD = 0.D0
ITER = 0

do while (iter.le.100001.and.DABS((E(1) - EOLD)/E(1)).GT.TOL)
      ITER = ITER + 1
      IF (ITER.GT.100000) then
            write(*,*) 'max iterations reached.'
            call exit(-1)
      end if

      DO I=1,3
            VT(I) = V(I) + VHXC(I)
            BXT(I) = BX(I) + BXCX(I)
            BYT(I) = BY(I) + BXCY(I)
            BZT(I) = BZ(I) + BXCZ(I)
      end do

      CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

      CALL ZHEEV( 'V', 'U', 6, M, 6, E, WORK, LWORK, RWORK, INFO )

!      DO 10 I=1,6
!10       WRITE(*,*)E(I)
!**----------------------------------------------------------------------
!**   calculate the densities and update the potentials
!**----------------------------------------------------------------------
      CALL GCALC(M,PHI,GAM)
      CALL DENCALC(M,N,MX,MY,MZ)

      DO I=1,3
            VHXC0(I) = VHXC(I)
            BXCX0(I) = BXCX(I)
            BXCY0(I) = BXCY(I)
            BXCZ0(I) = BXCZ(I)
      end do

!      N(3) = N(1)
!      MZ(3) = MZ(1)

      if (xc_slater.eq.1) then
            CALL XCPOT_SLATER(U0,0.d0,VXC,VHXC,BXCX,BXCY,BXCZ)
      else
            CALL XCPOT_BALDA( &
                  U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ,EC)
      end if

      DO I=1,3
            VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXC0(I)
            BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCX0(I)
            BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCY0(I)
            BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZ0(I)
      end do
            
end do

EH = 0.D0
DO I=1,3
      EH = EH - 0.5D0*U0*N(I)**2
end do

EVXC = 0.D0
DO I=1,3
      EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I) &
            -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)
end do

CALL EX_CALC(U0,EX)

EXC = EX + EC
ETOT = E(1) + E(2) + E(3) + EH + EVXC + EXC

write(*,*) m

do i = 1, 3
      write(wsite, '(a, i1, a)') 'site-', i, '_KStddensity.txt'
      open(20+i, file=wsite)
end do

do i = 1, 3
      write(20+i, wmat) 0.d0, n(i), mx(i), my(i), mz(i)
end do

open(100, file='planar_field.txt')

do i=1, np
      read(100, *) v(i), bx(i), by(i), bz(i)
end do

close(100)
!------------------------------------------------------------------------
!**   Begin time-dependent iterations.
!------------------------------------------------------------------------
m_temp = m
! ^ set updating density matrix to the converged groud state.

do t_it = 1, 1000
      TIME = TIME + DT

!      DO 16 I=1,NP
!         BXD(I) = BX(I) + 0.01D0*DSIN(OMEGA*(TIME-DT/2.D0))
!16    CONTINUE

      DO I=1,NP
            vhxc0(I) = vhxc(I)
            BXCX0(I) = BXCX(I)
            BXCY0(I) = BXCY(I)
            BXCZ0(I) = BXCZ(I)
      ENDDO

      ! IF (CORR.EQ.1) THEN
      !       CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX0,BXCY0,BXCZ0)
      !       CALL BCORR(MX,MY,MZ,BXCX0,BXCY0,BXCZ0,TT)
      ! ENDIF

      PC_COUNT=0
      do while (pc_count.lt.pc)

            DO I=1,NP
                  VT(I) = V(I) + vhxc(I)
                  BXT(I) = BX(I) + BXCX(I)
                  BYT(I) = BY(I) + BXCY(I)
                  BZT(I) = BZ(I) + BXCZ(I)
            ENDDO

            CALL MATRIX(M,C,CP,T,TP,V,BX,BY,BZ)

            call time_step(m, m_temp, c, t, u0, vxc, vhxc, bxcx1, bxcy1, bxcz1, dt)

            ! IF (CORR.EQ.1) THEN
            !       CALL DENCALC(M_temp,N,MX,MY,MZ)
            !       CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX1,BXCY1,BXCZ1)
            !       CALL BCORR(MX,MY,MZ,BXCX1,BXCY1,BXCZ1,TT)
            ! ENDIF
      
            DO I=1,NP
                  vhxc(I) = 0.5D0*(vhxc0(I) + vhxc1(I))
                  BXCX(I) = 0.5D0*(BXCX0(I) + BXCX1(I))
                  BXCY(I) = 0.5D0*(BXCY0(I) + BXCY1(I))
                  BXCZ(I) = 0.5D0*(BXCZ0(I) + BXCZ1(I))
            ENDDO

            PC_COUNT = PC_COUNT+1

      end do

      DO I=1,NP
            vhxc(I) = vhxc1(I)
            BXCX(I) = BXCX1(I)
            BXCY(I) = BXCY1(I)
            BXCZ(I) = BXCZ1(I)
      ENDDO

      CALL DENCALC(M_temp,N,MX,MY,MZ)
!       If (L.EQ.2) THEN

! !**   project out transverse BXC

!             DO I=1,NP
!                   BXCX(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I)) &
!                         *MX(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
!                   BXCY(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I)) &
!                         *MY(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
!                   BXCZ(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I)) &
!                         *MZ(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
!             ENDDO    

!       !      DO I=1,NP
!       !         BB = DSQRT(BXCX1(I)**2 + BXCY1(I)**2 + BXCZ1(I)**2)
!       !         BXCX(I) = BB*MX(I)/MM(I)
!       !         BXCY(I) = BB*MY(I)/MM(I)
!       !         BXCZ(I) = BB*MZ(I)/MM(I)
!       !      end do
!       ENDIF

      ! CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
      ! WRITE(LN,999)TIME,N(1),N(2),N(3),N(4)
      ! WRITE(LN+1,999)TIME,MX(1),MX(2),MX(3),MX(4)
      ! WRITE(LN+2,999)TIME,MY(1),MY(2),MY(3),MY(4)
      ! WRITE(LN+3,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      ! WRITE(LN+4,999)TIME,MM(1),MM(2),MM(3),MM(4)
      ! WRITE(LN+5,999)TIME,MTOT(1),MTOT(2),MTOT(3), &
      !       DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
      ! IF (L.EQ.1) THEN
      !       WRITE(LN+6,999)TIME,TX(1),TX(2),TX(3),TX(4)
      !       WRITE(LN+7,999)TIME,TY(1),TY(2),TY(3),TY(4)
      !       WRITE(LN+8,999)TIME,TZ(1),TZ(2),TZ(3),TZ(4)
      !       WRITE(LN+9,997)TIME,TT(1),TT(2),TT(3)
      ! ENDIF

      do i = 1, 3
            write(20+i, wmat) time, dble(n(i)), dble(mx(i)), dble(my(i)), dble(mz(i))
            ! write(20+i,*) '##############################'
      end do

      do i = 1, 3
            close(20+i)
      end do

end do

END

!************************************************************************
!************************************************************************
!************************************************************************

subroutine time_step(m, m_time, c, t, u0, vxc, vhxc, bxcx1, bxcy1, bxcz1, dt)    
implicit none

reaL(8) :: vxc(3), vhxc(3), BXCX1(3), BXCY1(3), BXCZ1(3)
real(8) :: cp, tp, c, t, dt, u0
complex(8) :: m(6,6), m_time(6,6), a(6,6), r(6)
complex(8), parameter :: one = (1.d0, 0.d0), ione = (0.d0, 1.d0), zero = (0.d0, 0.d0)
integer :: i, j, ipiv(6), info
complex(8) :: gam(2,2,3,3),PHI(6,6)
real(8) :: E(6)

COMMON / evals / e, phi, gam

cp = c
tp = t

DO I=1,6
      DO J=1,6
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

DO I=1,6
      R(I) = ZERO
      DO J=1,6
            R(I) = R(I) + A(I,J)*m_time(J, 1)
      ENDDO
ENDDO    

DO I=1,6
      DO J=1,6
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

DO I=1,6
      m_time(I, 1) = R(I)
ENDDO

DO I=1,6
      DO J=1,6
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

DO I=1,6
      R(I) = ZERO
      DO J=1,6
            R(I) = R(I) + A(I,J)*m_time(J, 2)
      ENDDO
ENDDO    

DO I=1,6
      DO J=1,6
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

DO I=1,6
      m_time(I, 2) = R(I)
ENDDO

DO I=1,6
      DO J=1,6
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

DO I=1,6
      R(I) = ZERO
      DO J=1,6
            R(I) = R(I) + A(I,J)*m_time(J, 3)
      ENDDO
ENDDO    

DO I=1,6
      DO J=1,6
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
      ENDDO    
ENDDO    

CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

DO I=1,6
      m_time(I, 3) = R(I)
ENDDO

CALL GCALC(M_time,PHI,gam)
CALL XCPOT_SLATER(U0,0.d0,vxc,vhxc,BXCX1,BXCY1,BXCZ1)

end subroutine

!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE MATRIX(M,C,CP,T,TP,V,BX,BY,BZ)
IMPLICIT NONE

INTEGER :: I,J
real(8) :: C,CP,T,TP,V(3),BX(3),BY(3),BZ(3)
complex(8) :: M(6,6),ZERO,ONE,IONE
!complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0)

DO I=1,6
      DO J=1,6
            M(I,J) = ZERO
      end do
end do

DO I=1,3
      M(I,I) = (V(I) + BZ(I))*ONE
      M(I,I+3) = BX(I)*ONE - BY(I)*IONE
      M(I+3,I+3) = (V(I) - BZ(I))*ONE
      M(I+3,I) = BX(I)*ONE + BY(I)*IONE
end do

M(1,2) = -T*ONE + C*IONE
M(2,3) = -T*ONE + C*IONE
M(1,3) = -TP*ONE - CP*IONE
M(4,5) = -T*ONE - C*IONE
M(5,6) = -T*ONE - C*IONE
M(4,6) = -TP*ONE + CP*IONE

DO I=1,2
      DO J=I+1,3
            M(J,I) = DCONJG(M(I,J))
      end do
end do

END subroutine 
!************************************************************************
!************************************************************************
!************************************************************************
SUBROUTINE GCALC(M,PHI,gam)
IMPLICIT NONE

INTEGER :: I,J,K,L,NP
complex(8) :: M(2*NP,2*NP),gam(2,2,NP,NP),PHI(2*NP,2,NP)

!**   Define the orbitals phi(m,sigma,x):
!**   m = 1...2*NP is the orbital index
!**   sigma = 1,2 (up, down) is the spin index
!**   x = 1,...,NP (lattice points) is the spatial coordinate

DO I=1,2*NP
      DO J=1,NP
            PHI(I,1,J) = M(J,I)
            PHI(I,2,J) = M(J+NP,I)
      end do
end do

DO I=1,2
      DO J=1,2
            DO K=1,NP
                  DO L=1,NP
                        gam(I,J,K,L) = PHI(1,I,K)*DCONJG(PHI(1,J,L)) & 
                                    + PHI(2,I,K)*DCONJG(PHI(2,J,L)) & 
                                    + PHI(3,I,K)*DCONJG(PHI(3,J,L))
                  end do
            end do
      end do
end do

END subroutine
!************************************************************************
!************************************************************************
!************************************************************************
SUBROUTINE DENCALC(M,N,MX,MY,MZ)
IMPLICIT NONE

INTEGER :: I
complex(8) :: M(6,6),NUU(3),NUD(3),NDU(3),NDD(3)
real(8) :: N(3),MX(3),MY(3),MZ(3)

DO I=1,3
      NUU(I) = CDABS(M(I,1))**2+CDABS(M(I,2))**2+CDABS(M(I,3))**2
      NUD(I) = M(I,1)*DCONJG(M(I+3,1)) + M(I,2)*DCONJG(M(I+3,2)) &
          + M(I,3)*DCONJG(M(I+3,3))
      NDU(I) = DCONJG(NUD(I))
      NDD(I) = CDABS(M(I+3,1))**2 + CDABS(M(I+3,2))**2 &
          + CDABS(M(I+3,3))**2
end do

DO I=1,3
      N(I) = DREAL(NUU(I) + NDD(I))
      MX(I) = DREAL(NUD(I) + NDU(I))
      MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
      MZ(I) = DREAL(NUU(I) - NDD(I))
end do

END subroutine
!**********************************************************************
!**********************************************************************
!**********************************************************************
SUBROUTINE XCPOT_BALDA(U0,T,VXC,VHXC,BXCX,BXCY,BXCZ, &
                       NN,MX,MY,MZ,EC)
IMPLICIT NONE

INTEGER NP,I
PARAMETER (NP=3)
real(8) :: U0,UU,PI,EDR,T,EC,A,B,C,BETU,S, &
                 ALPHA,BETA,gam,BETA_DM,BETA_DN,ALPHA_DM, &
                 ALPHA_DN,gam_DM,gam_DN,GARG,EHOM
real(8) :: N(NP),NN(NP),MX(NP),MY(NP),MZ(NP),M(NP), &
                 ECD(NP),VC,BC,BCX,BCY,BCZ
real(8) :: VXC(NP),VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

PI = 3.141592653589793D0
EDR = 1.D0/3.D0

UU = U0/T

A = 0.7504D0
B = 0.1479D0
C = 0.5574D0
BETU = (2.D0 + A*UU + B*UU**2)/(1.D0 + C*UU + B*UU**2)

DO I=1,NP
      N(I) = NN(I)
      M(I) = DSQRT(MX(I)**2+MY(I)**2+MZ(I)**2)

!**      This is the exchange potential and magnetic field:
      VXC(I) = -U0*N(I)/2.D0
      BXCX(I) = -U0*MX(I)/2.D0
      BXCY(I) = -U0*MY(I)/2.D0
      BXCZ(I) = -U0*MZ(I)/2.D0

      S = 1.D0
      IF (N(I).GE.1.D0) THEN
      N(I) = 2.D0-N(I)
      S = -1.D0
      ENDIF

      ALPHA = ( (N(I)**2-M(I)**2)/N(I)**1.875D0)**(UU**EDR)
      ALPHA_DN = UU**EDR*(N(I)**0.125D0 - M(I)**2/N(I)**1.875D0) &
                  **(UU**EDR - 1.D0) &
                  * (1.D0 + 15D0*M(I)**2/N(I)**2)/(8.D0*N(I)**0.875D0)
      ALPHA_DM = -UU**EDR*(N(I)**0.125D0 - M(I)**2/N(I)**1.875D0) &
            **(UU**EDR - 1.D0) * 2.D0*M(I)/N(I)**1.875D0

      BETA = BETU**ALPHA
      BETA_DN = ALPHA*BETU**(ALPHA-1.D0) * ALPHA_DN
      BETA_DM = ALPHA*BETU**(ALPHA-1.D0) * ALPHA_DM

      GARG = DSQRT(UU)/(1.D0 - (M(I)/N(I))**1.5D0)

      IF (GARG.LT.46.D0) THEN
            gam = 2.D0*DEXP(GARG)
      ELSE
            gam = 1.D20
      ENDIF
      gam_DN = -gam*DSQRT(UU)/(1.D0-(M(I)/N(I))**1.5D0)**2 &
            * (1.5D0*M(I)/N(I)**2)*DSQRT(M(I)/N(I))
      gam_DM =  gam*DSQRT(UU)/(1.D0-(M(I)/N(I))**1.5D0)**2 &
            * (1.5D0/N(I))*DSQRT(M(I)/N(I))

      EHOM = -(2.D0*T/PI)*BETA*DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/gam)

      ECD(I) = EHOM + U0*(M(I)**2-N(I)**2)/4.D0 + (4.D0*T/PI) &
                * DSIN(PI*N(I)/2.D0)*DCOS(PI*M(I)/2.D0)

      VC = -(2.D0/PI)*BETA_DN &
               *DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/gam) &
             - 2.D0*(1.D0 - N(I)*BETA_DN/BETA) &
               *DCOS(PI*N(I)/BETA)*DCOS(PI*M(I)/gam) &
             - 2.D0*BETA*(M(I)*gam_DN/gam**2) &
               *DSIN(PI*N(I)/BETA)*DSIN(PI*M(I)/gam)

      BC = -(2.D0/PI)*BETA_DM &
         *DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/gam) &
       + 2.D0*(N(I)*BETA_DM/BETA) &
         *DCOS(PI*N(I)/BETA)*DCOS(PI*M(I)/gam) &
       + 2.D0*BETA*(1.D0/gam - M(I)*gam_DM/gam**2) &
         *DSIN(PI*N(I)/BETA)*DSIN(PI*M(I)/gam)

      VC = T*VC - U0*N(I)/2.D0 &
        + 2.D0*T*DCOS(PI*N(I)/2.D0)*DCOS(PI*M(I)/2.D0)

!         VXC(I) = VXC(I) + S*VC
      VHXC(I) = VXC(I) + U0*NN(I)

      IF (M(I).GT.1.D-15) THEN
            BCX = T*BC*MX(I)/M(I) + U0*MX(I)/2.D0 &
                        - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MX(I)/M(I)

            BCY = T*BC*MY(I)/M(I) + U0*MY(I)/2.D0 &
                        - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MY(I)/M(I)

            BCZ = T*BC*MZ(I)/M(I) + U0*MZ(I)/2.D0 &
                        - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MZ(I)/M(I)

!         BXCX(I) = BXCX(I) + BCX
!         BXCY(I) = BXCY(I) + BCY
!         BXCZ(I) = BXCZ(I) + BCZ
      ENDIF

end do

EC = 0.D0
!      DO I=1,NP
!           EC = EC + ECD(I)
!      end do

END subroutine

!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE EX_CALC(U0,EX)
IMPLICIT NONE

INTEGER NP,I,TAU,SIGMA
PARAMETER (NP=3)
complex(8) :: gam(2,2,NP,NP),PHI(2*NP,2,NP)
real(8) :: E(2*NP),U0,EX

COMMON /EVALS/ E,PHI,gam

EX = 0.D0
DO I=1,NP
      DO TAU=1,2
            DO SIGMA=1,2
                  EX = EX - 0.5D0*U0*gam(SIGMA,TAU,I,I)*gam(TAU,SIGMA,I,I)
            end do
      end do
end do

END subroutine

!************************************************************************
!************************************************************************
!************************************************************************

SUBROUTINE XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
IMPLICIT NONE

INTEGER :: I
integer, PARAMETER :: NP = 3
complex(8) :: NUU(NP),NUD(NP),NDU(NP),NDD(NP), &
                  VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP), &
                  BUU(NP),BUD(NP),BDU(NP),BDD(NP),MAT(4,4), &
                  gam(2,2,NP,NP),PHI(2*NP,2,NP)
complex(8), PARAMETER :: IONE=(0.D0,1.D0)
real(8) :: U0,U1,N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP), &
                  VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP)

COMMON /EVALS/ E,PHI,gam

DO I=1,NP
      NUU(I) = CDABS(PHI(1,1,I))**2 + CDABS(PHI(2,1,I))**2
      NUD(I) = PHI(1,1,I)*DCONJG(PHI(1,2,I)) &
                + PHI(2,1,I)*DCONJG(PHI(2,2,I))
      NDU(I) = DCONJG(NUD(I))
      NDD(I) = CDABS(PHI(1,2,I))**2 + CDABS(PHI(2,2,I))**2
end do

DO I=1,NP
      N(I) = DREAL(NUU(I) + NDD(I))
      MX(I) = DREAL(NUD(I) + NDU(I))
      MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
      MZ(I) = DREAL(NUU(I) - NDD(I))
end do

VH(1) = U0*N(1) + U1*N(2)
DO I=2,NP-1
      VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
end do
VH(NP) = U0*N(NP) + U1*N(NP-1)

DO I=1,NP

      BUU(I) = -2.D0*U0*(NUU(I)*NUU(I) + NUD(I)*NDU(I))
      BDU(I) = -2.D0*U0*(NDU(I)*NUU(I) + NDD(I)*NDU(I))
      BUD(I) = -2.D0*U0*(NUU(I)*NUD(I) + NUD(I)*NDD(I))
      BDD(I) = -2.D0*U0*(NDU(I)*NUD(I) + NDD(I)*NDD(I))

      IF (I.LT.NP) THEN
            BUU(I) = BUU(I) -2.D0*U1*(gam(1,1,I,I+1)*gam(1,1,I+1,I) &
                                     +gam(1,2,I,I+1)*gam(2,1,I+1,I))
            BDU(I) = BDU(I) -2.D0*U1*(gam(2,1,I,I+1)*gam(1,1,I+1,I) &
                                     +gam(2,2,I,I+1)*gam(2,1,I+1,I))
            BUD(I) = BUD(I) -2.D0*U1*(gam(1,1,I,I+1)*gam(1,2,I+1,I) &
                                     +gam(1,2,I,I+1)*gam(2,2,I+1,I))
            BDD(I) = BDD(I) -2.D0*U1*(gam(2,1,I,I+1)*gam(1,2,I+1,I) &
                                     +gam(2,2,I,I+1)*gam(2,2,I+1,I))
      ENDIF

      IF (I.GT.1) THEN
            BUU(I) = BUU(I) -2.D0*U1*(gam(1,1,I,I-1)*gam(1,1,I-1,I) &
                                     +gam(1,2,I,I-1)*gam(2,1,I-1,I))
            BDU(I) = BDU(I) -2.D0*U1*(gam(2,1,I,I-1)*gam(1,1,I-1,I) &
                                     +gam(2,2,I,I-1)*gam(2,1,I-1,I))
            BUD(I) = BUD(I) -2.D0*U1*(gam(1,1,I,I-1)*gam(1,2,I-1,I) &
                                     +gam(1,2,I,I-1)*gam(2,2,I-1,I))
            BDD(I) = BDD(I) -2.D0*U1*(gam(2,1,I,I-1)*gam(1,2,I-1,I) &
                                     +gam(2,2,I,I-1)*gam(2,2,I-1,I))
      ENDIF

end do

DO I=1,NP
      DEN(I) = 2.D0*N(I)*(NUU(I)*NDD(I)-NUD(I)*NDU(I))

      MAT(1,1) = N(I)*NDD(I) - NUD(I)*NDU(I)
      MAT(1,2) = -NDD(I)*NUD(I)
      MAT(1,3) = -NDD(I)*NDU(I)
      MAT(1,4) = NUD(I)*NDU(I)

      MAT(2,1) = -NDD(I)*NDU(I)
      MAT(2,2) = 2.D0*NUU(I)*NDD(I) - NUD(I)*NDU(I)
      MAT(2,3) = NDU(I)**2
      MAT(2,4) = -NUU(I)*NDU(I)

      MAT(3,1) = -NDD(I)*NUD(I)
      MAT(3,2) = NUD(I)**2
      MAT(3,3) = 2.D0*NUU(I)*NDD(I) - NUD(I)*NDU(I)
      MAT(3,4) = -NUU(I)*NUD(I)

      MAT(4,1) = NUD(I)*NDU(I)
      MAT(4,2) = -NUU(I)*NUD(I)
      MAT(4,3) = -NUU(I)*NDU(I)
      MAT(4,4) = N(I)*NUU(I) - NUD(I)*NDU(I)

      VUU(I) = ( MAT(1,1)*BUU(I) + MAT(1,2)*BDU(I) &
            + MAT(1,3)*BUD(I) + MAT(1,4)*BDD(I) )/DEN(I)
      VDU(I) = ( MAT(2,1)*BUU(I) + MAT(2,2)*BDU(I) &
            + MAT(2,3)*BUD(I) + MAT(2,4)*BDD(I) )/DEN(I)
      VUD(I) = ( MAT(3,1)*BUU(I) + MAT(3,2)*BDU(I) &
            + MAT(3,3)*BUD(I) + MAT(3,4)*BDD(I) )/DEN(I)
      VDD(I) = ( MAT(4,1)*BUU(I) + MAT(4,2)*BDU(I) &
            + MAT(4,3)*BUD(I) + MAT(4,4)*BDD(I) )/DEN(I)

      VXC (I) = DREAL(VUU(I) + VDD(I))/2.D0
      BXCX(I) = DREAL(VDU(I) + VUD(I))/2.D0
      BXCY(I) = DREAL(-IONE*VDU(I) + IONE*VUD(I))/2.D0
      BXCZ(I) = DREAL(VUU(I) - VDD(I))/2.D0

      VHXC(I) = VH(I) + VXC(I)
end do

END subroutine
