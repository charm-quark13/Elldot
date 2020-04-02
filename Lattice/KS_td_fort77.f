      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER :: I,j,k,INFO,MODE,ITER, xc_slater, corr, t_it, rest,
     &            pc_count, L
      !integer :: J,k, vec, krev, xc_guess
      integer, parameter :: lwork=100, np=3, pc=2

      real(8) :: C,CP,T,TP,V(3),BX(3),BY(3),BZ(3),
     &           VT(3),BXT(3),BYT(3),BZT(3),VXC(3),
     &           VHXC(3),BXCX(3),BXCY(3),BXCZ(3),
     &           vxc0(3),VHXC0(3),BXCX0(3),BXCY0(3),BXCZ0(3),
     &           N(3),MX(3),MY(3),MZ(3),E(6),RWORK(100), 
     &           dt, time, vss(3),bxss(3),byss(3),bzss(3),
     &           TT(3),TX(3),TY(3),TZ(3), flag_vec(3),
     &           vxc1(3),vhxc1(3),BXCX1(3),BXCY1(3),BXCZ1(3)
      real(8) :: U0,MIX,TOL,EOLD,EH,ETOT,EVXC,EXC,EX,EC,CRIT

      character(100) ::  wmat, wsite

      complex(8) :: M(6,6),WORK(LWORK),GAMMA(2,2,3,3),PHI(2*NP,2,NP), 
     &              m_temp(6,6), m_temp1(6,6)
      complex(8), parameter :: one = (1.d0, 0.d0), 
     &                         ione = (0.d0, 1.d0), 
     &                         zero = (0.d0, 0.d0)
 
      COMMON /EVALS/ E,PHI,GAMMA

      write(wmat,'(a)') '(5f15.9)'

      dt = .01d0
      time = 0.d0
      U0 = 1.d0
      L = 1
      
      phi = zero

      xc_slater = 1
      corr = 0
      mode = 2

      MIX = 0.1D0
      TOL = 1.D-9

      open(1, file='planar_field.txt')
      DO I=1,3
            READ(1,*)V(I),BX(I),BY(I),BZ(I)
            N(I) = 0.D0
            VHXC(I) = 0.D0
            BXCX(I) = 0.D0
            BXCY(I) = 0.D0
            BXCZ(I) = 0.D0
      end do
      close(1)

      open(1, file='constants.txt')
      read(1,*) c, t
      close(1)

      U0 = 1.d0

      TP = T

      cp = c

      ! E(1) = 1.d0
      EOLD = 0.D0
      ITER = 0

      itloop : do iter = 1, 200001 
            ! ITER = ITER + 1
            IF (ITER.GT.200000) then
                  write(*,*) 'max iterations reached.'
                  call exit(-1)
            end if

            DO I=1,3
                  VT(I) = V(I) + VHXC(I) * 0.d0
                  BXT(I) = BX(I) + BXCX(I) * 0.d0
                  BYT(I) = BY(I) + BXCY(I) * 0.d0
                  BZT(I) = BZ(I) + BXCZ(I) * 0.d0
            end do

            CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

            CALL ZHEEV( 'V', 'U', 6, M, 6, E, WORK, LWORK, RWORK, INFO )

      !      DO 10 I=1,6
      !10       WRITE(*,*)E(I)
      !**----------------------------------------------------------------------
      !**   calculate the densities and update the potentials
      !**----------------------------------------------------------------------
            CALL GCALC(M,PHI,gamma)
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
                  CALL XCPOT_BALDA(
     &                  U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ,EC)
            end if

            DO I=1,3
                  VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXC0(I)
                  BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCX0(I)
                  BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCY0(I)
                  BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZ0(I)
            end do
            
            if (DABS((E(1) - EOLD)/E(1)).LT.TOL) then
                  write(*,*) DABS((E(1) - EOLD)/E(1))
                  write(*,*) 'convergence met,', 
     &                              'exiting convergence loop...'
                  exit itloop
            end if
            eold = e(1)
      end do itloop

      EH = 0.D0
      DO I=1,3
            EH = EH - 0.5D0*U0*N(I)**2
      end do

      EVXC = 0.D0
      DO I=1,3
            EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &             -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)
      end do

      CALL EX_CALC(U0,EX)

      EXC = EX + EC
      ETOT = E(1) + E(2) + E(3) + EH + EVXC + EXC

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
      ! ^ set density matrix to the converged groud state.


      timeprop : do t_it = 1, 1000
            TIME = TIME + DT

C      !      DO 16 I=1,NP
C      !         BXD(I) = BX(I) + 0.01D0*DSIN(OMEGA*(TIME-DT/2.D0))
C      !16    CONTINUE

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

            DO I=1,NP
                  bxcx = 0.d0
                  bxcy = 0.d0
                  bxcz = 0.d0
                  vhxc = 0.d0
            ENDDO


            pc_step : do pc_count = 1, 3
! C                  write(*,*) pc_count
                  DO I=1,NP
                        VT(I) = V(I) + vhxc(I)
                        BXT(I) = BX(I) + BXCX(I)
                        BYT(I) = BY(I) + BXCY(I)
                        BZT(I) = BZ(I) + BXCZ(I)
                  ENDDO

                  CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

                  call time_step(m, m_temp, m_temp1, dt)

                  call time_slater(m_temp1,u0,vhxc1,bxcx1,bxcy1,bxcz1)

                  ! write(*,*) m_temp1

!                   if (xc_slater.eq.1) then
!                         CALL XCPOT_SLATER(U0,0.d0,vxc1,vhxc1,BXCX1,BXCY1
!      &                                          ,BXCZ1)
!                   else
!                         CALL XCPOT_BALDA(
!      &                        U0,T,VXC1,VHXC1,BXCX1,BXCY1,BXCZ1
!      &                        ,N,MX,MY,MZ,EC
!      &                  )
!                   end if
                  ! IF (CORR.EQ.1) THEN
                  !       CALL DENCALC(M_temp,N,MX,MY,MZ)
                  !       CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX1,BXCY1,BXCZ1)
                  !       CALL BCORR(MX,MY,MZ,BXCX1,BXCY1,BXCZ1,TT)
                  ! ENDIF
                  if (pc_count.lt.3) then
                        DO I=1,NP
                              vhxc(I) = 0.5D0*(vhxc0(I) + vhxc1(I))
                              BXCX(I) = 0.5D0*(BXCX0(I) + BXCX1(I))
                              BXCY(I) = 0.5D0*(BXCY0(I) + BXCY1(I))
                              BXCZ(I) = 0.5D0*(BXCZ0(I) + BXCZ1(I))
                        ENDDO
                  end if
            end do pc_step

            ! m_temp = m_temp1

            DO I=1,NP
                  vhxc(I) = vhxc1(I)
                  BXCX(I) = BXCX1(I)
                  BXCY(I) = BXCY1(I)
                  BXCZ(I) = BXCZ1(I)
            ENDDO

            CALL DENCALC(M_temp1,N,MX,MY,MZ)
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

            do i = 1, 3
                  write(20+i, wmat) time, dble(n(i)), dble(mx(i)),
     &                                  dble(my(i)), dble(mz(i))
            end do


            m_temp = m_temp1
      end do timeprop

      do i = 1, 3
            close(20+i)
      end do


      END

!************************************************************************
!************************************************************************
!************************************************************************

      subroutine time_step(H, m_time, m_time1, dt)    
      implicit none

      real(8) :: cp, tp, c, t, dt, u0
      complex(8) :: H(6,6), m_time(6,6), m_time1(6,6), a(6,6), r(6)
      complex(8), parameter :: 
     &      one = (1.d0, 0.d0), ione = (0.d0, 1.d0), zero = (0.d0, 0.d0)
      integer :: i, j, ipiv(6), info

      m_time1 = zero

      ! rhs state 1

      DO I=1,6
            DO J=1,6
                  A(I,J) = -0.5D0*IONE*DT*H(I,J) 
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
                  A(I,J) = 0.5D0*IONE*DT*H(I,J) 
                  IF (I.EQ.J) A(I,J) = ONE + A(I,J)
            ENDDO    
      ENDDO    

      CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

      DO I=1,6
            m_time1(I, 1) = R(I)
      ENDDO

      DO I=1,6
            DO J=1,6
                  A(I,J) = -0.5D0*IONE*DT*H(I,J) 
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
                  A(I,J) = 0.5D0*IONE*DT*H(I,J) 
                  IF (I.EQ.J) A(I,J) = ONE + A(I,J)
            ENDDO    
      ENDDO    

      CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

      DO I=1,6
            m_time1(I, 2) = R(I)
      ENDDO

      DO I=1,6
            DO J=1,6
                  A(I,J) = -0.5D0*IONE*DT*H(I,J) 
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
                  A(I,J) = 0.5D0*IONE*DT*H(I,J) 
                  IF (I.EQ.J) A(I,J) = ONE + A(I,J)
            ENDDO    
      ENDDO    

      CALL ZGESV( 6, 1, A, 6, IPIV, R, 6, INFO )

      DO I=1,6
            m_time1(I, 3) = R(I)
      ENDDO

      end subroutine

C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE MATRIX(M,C,CP,T,TP,V,BX,BY,BZ)
      IMPLICIT NONE

      INTEGER I,J
      DOUBLE PRECISION C,CP,T,TP,V(3),BX(3),BY(3),BZ(3)
      DOUBLE COMPLEX M(6,6),ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      DO 1 I=1,6
      DO 1 J=1,6
1        M(I,J) = ZERO

      DO 10 I=1,3
         M(I,I) = (V(I) + BZ(I))*ONE
         M(I,I+3) = BX(I)*ONE - BY(I)*IONE
         M(I+3,I+3) = (V(I) - BZ(I))*ONE
         M(I+3,I) = BX(I)*ONE + BY(I)*IONE
10    CONTINUE

      M(1,2) = -T*ONE + C*IONE
      M(2,3) = -T*ONE + C*IONE
      M(1,3) = -TP*ONE - CP*IONE
      M(4,5) = -T*ONE - C*IONE
      M(5,6) = -T*ONE - C*IONE
      M(4,6) = -TP*ONE + CP*IONE

      DO 2 I=1,2
      DO 2 J=I+1,3
         M(J,I) = DCONJG(M(I,J))
2     CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE GCALC(M,PHI,GAMMA)
      IMPLICIT NONE

      INTEGER I,J,K,L,NP
      PARAMETER (NP=3)
      DOUBLE COMPLEX M(2*NP,2*NP),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)

C**   Define the orbitals phi(m,sigma,x):
C**   m = 1...2*NP is the orbital index
C**   sigma = 1,2 (up, down) is the spin index
C**   x = 1,...,NP (lattice points) is the spatial coordinate

      DO 9 I=1,2*NP
      DO 9 J=1,NP
         PHI(I,1,J) = M(J,I)
         PHI(I,2,J) = M(J+NP,I)
9     CONTINUE

      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,NP
      DO 10 L=1,NP
         GAMMA(I,J,K,L) = PHI(1,I,K)*DCONJG(PHI(1,J,L))
     &                  + PHI(2,I,K)*DCONJG(PHI(2,J,L))
     &                  + PHI(3,I,K)*DCONJG(PHI(3,J,L))
10    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE DENCALC(M,N,MX,MY,MZ)
      IMPLICIT NONE

      INTEGER I
      DOUBLE COMPLEX M(6,6),NUU(3),NUD(3),NDU(3),NDD(3)
      DOUBLE PRECISION N(3),MX(3),MY(3),MZ(3)

      DO 5 I=1,3
         NUU(I) = CDABS(M(I,1))**2+CDABS(M(I,2))**2+CDABS(M(I,3))**2
         NUD(I) = M(I,1)*DCONJG(M(I+3,1)) + M(I,2)*DCONJG(M(I+3,2))
     &          + M(I,3)*DCONJG(M(I+3,3))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(M(I+3,1))**2 + CDABS(M(I+3,2))**2
     &          + CDABS(M(I+3,3))**2
5     CONTINUE

      DO 10 I=1,3
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
10    CONTINUE

      RETURN
      END
C**********************************************************************
C**********************************************************************
C**********************************************************************
      SUBROUTINE XCPOT_BALDA(U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,
     &                       NN,MX,MY,MZ,EC)
      IMPLICIT NONE

      INTEGER NP,I
      PARAMETER (NP=3)
      DOUBLE PRECISION U0,UU,PI,EDR,T,EC,A,B,C,BETU,S,
     &                 ALPHA,BETA,GAMMA,BETA_DM,BETA_DN,ALPHA_DM,
     &                 ALPHA_DN,GAMMA_DM,GAMMA_DN,GARG,EHOM
      DOUBLE PRECISION N(NP),NN(NP),MX(NP),MY(NP),MZ(NP),M(NP),
     &                 ECD(NP),VC,BC,BCX,BCY,BCZ
      DOUBLE PRECISION VXC(NP),VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

      PI = 3.141592653589793D0
      EDR = 1.D0/3.D0

      UU = U0/T

      A = 0.7504D0
      B = 0.1479D0
      C = 0.5574D0
      BETU = (2.D0 + A*UU + B*UU**2)/(1.D0 + C*UU + B*UU**2)

      DO 1 I=1,NP
         N(I) = NN(I)
         M(I) = DSQRT(MX(I)**2+MY(I)**2+MZ(I)**2)

C**      This is the exchange potential and magnetic field:
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
         ALPHA_DN = UU**EDR*(N(I)**0.125D0 - M(I)**2/N(I)**1.875D0)
     &            **(UU**EDR - 1.D0)
     &            * (1.D0 + 15D0*M(I)**2/N(I)**2)/(8.D0*N(I)**0.875D0)
         ALPHA_DM = -UU**EDR*(N(I)**0.125D0 - M(I)**2/N(I)**1.875D0)
     &            **(UU**EDR - 1.D0) * 2.D0*M(I)/N(I)**1.875D0

         BETA = BETU**ALPHA
         BETA_DN = ALPHA*BETU**(ALPHA-1.D0) * ALPHA_DN
         BETA_DM = ALPHA*BETU**(ALPHA-1.D0) * ALPHA_DM

         GARG = DSQRT(UU)/(1.D0 - (M(I)/N(I))**1.5D0)

         IF (GARG.LT.46.D0) THEN
            GAMMA = 2.D0*DEXP(GARG)
         ELSE
            GAMMA = 1.D20
         ENDIF
         GAMMA_DN = -GAMMA*DSQRT(UU)/(1.D0-(M(I)/N(I))**1.5D0)**2
     &            * (1.5D0*M(I)/N(I)**2)*DSQRT(M(I)/N(I))
         GAMMA_DM =  GAMMA*DSQRT(UU)/(1.D0-(M(I)/N(I))**1.5D0)**2
     &            * (1.5D0/N(I))*DSQRT(M(I)/N(I))

         EHOM = -(2.D0*T/PI)*BETA*DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/GAMMA)

         ECD(I) = EHOM + U0*(M(I)**2-N(I)**2)/4.D0 + (4.D0*T/PI)
     &          * DSIN(PI*N(I)/2.D0)*DCOS(PI*M(I)/2.D0)

         VC = -(2.D0/PI)*BETA_DN
     &         *DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/GAMMA)
     &       - 2.D0*(1.D0 - N(I)*BETA_DN/BETA)
     &         *DCOS(PI*N(I)/BETA)*DCOS(PI*M(I)/GAMMA)
     &       - 2.D0*BETA*(M(I)*GAMMA_DN/GAMMA**2)
     &         *DSIN(PI*N(I)/BETA)*DSIN(PI*M(I)/GAMMA)

         BC = -(2.D0/PI)*BETA_DM
     &         *DSIN(PI*N(I)/BETA)*DCOS(PI*M(I)/GAMMA)
     &       + 2.D0*(N(I)*BETA_DM/BETA)
     &         *DCOS(PI*N(I)/BETA)*DCOS(PI*M(I)/GAMMA)
     &       + 2.D0*BETA*(1.D0/GAMMA - M(I)*GAMMA_DM/GAMMA**2)
     &         *DSIN(PI*N(I)/BETA)*DSIN(PI*M(I)/GAMMA)

         VC = T*VC - U0*N(I)/2.D0
     &        + 2.D0*T*DCOS(PI*N(I)/2.D0)*DCOS(PI*M(I)/2.D0)

C         VXC(I) = VXC(I) + S*VC
         VHXC(I) = VXC(I) + U0*NN(I)

         IF (M(I).GT.1.D-15) THEN
         BCX = T*BC*MX(I)/M(I) + U0*MX(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MX(I)/M(I)

         BCY = T*BC*MY(I)/M(I) + U0*MY(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MY(I)/M(I)

         BCZ = T*BC*MZ(I)/M(I) + U0*MZ(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MZ(I)/M(I)

C         BXCX(I) = BXCX(I) + BCX
C         BXCY(I) = BXCY(I) + BCY
C         BXCZ(I) = BXCZ(I) + BCZ
         ENDIF

1     CONTINUE

      EC = 0.D0
C      DO 55 I=1,NP
C         EC = EC + ECD(I)
C55    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE EX_CALC(U0,EX)
      IMPLICIT NONE

      INTEGER NP,I,TAU,SIGMA
      PARAMETER (NP=3)
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)
      DOUBLE PRECISION E(2*NP),U0,EX

      COMMON /EVALS/ E,PHI,GAMMA

      EX = 0.D0
      DO 10 I=1,NP
      DO 10 TAU=1,2
      DO 10 SIGMA=1,2
         EX = EX -0.5D0*U0*GAMMA(SIGMA,TAU,I,I)*GAMMA(TAU,SIGMA,I,I)
10    CONTINUE

      RETURN
      END
C************************************************************************

      SUBROUTINE XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,I
      PARAMETER (NP = 3)
      DOUBLE COMPLEX NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP),
     &               BUU(NP),BUD(NP),BDU(NP),BDD(NP),MAT(4,4),IONE,
     &               GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)
      PARAMETER (IONE=(0.D0,1.D0))
      DOUBLE PRECISION U0,U1,N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP)

      COMMON /EVALS/ E,PHI,GAMMA

      DO 1 I=1,NP
         NUU(I) = CDABS(PHI(1,1,I))**2 + CDABS(PHI(2,1,I))**2
         NUD(I) = PHI(1,1,I)*DCONJG(PHI(1,2,I))
     &          + PHI(2,1,I)*DCONJG(PHI(2,2,I))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(PHI(1,2,I))**2 + CDABS(PHI(2,2,I))**2
1     CONTINUE

      DO 2 I=1,NP
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
2     CONTINUE

      VH(1) = U0*N(1) + U1*N(2)
      DO 3 I=2,NP-1
         VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
3     CONTINUE
      VH(NP) = U0*N(NP) + U1*N(NP-1)

      DO 5 I=1,NP

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

5     CONTINUE

      DO 10 I=1,NP
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

         VUU(I) = ( MAT(1,1)*BUU(I) + MAT(1,2)*BDU(I)
     &            + MAT(1,3)*BUD(I) + MAT(1,4)*BDD(I) )/DEN(I)
         VDU(I) = ( MAT(2,1)*BUU(I) + MAT(2,2)*BDU(I)
     &            + MAT(2,3)*BUD(I) + MAT(2,4)*BDD(I) )/DEN(I)
         VUD(I) = ( MAT(3,1)*BUU(I) + MAT(3,2)*BDU(I)
     &            + MAT(3,3)*BUD(I) + MAT(3,4)*BDD(I) )/DEN(I)
         VDD(I) = ( MAT(4,1)*BUU(I) + MAT(4,2)*BDU(I)
     &            + MAT(4,3)*BUD(I) + MAT(4,4)*BDD(I) )/DEN(I)

         VXC (I) = DREAL(VUU(I) + VDD(I))/2.D0
         BXCX(I) = DREAL(VDU(I) + VUD(I))/2.D0
         BXCY(I) = DREAL(-IONE*VDU(I) + IONE*VUD(I))/2.D0
         BXCZ(I) = DREAL(VUU(I) - VDD(I))/2.D0

         VHXC(I) = VH(I) + VXC(I)
10    CONTINUE

      RETURN
      END

C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE time_slater(M,U0,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE
      
      INTEGER I,J,K,L,R,occ
      integer, parameter :: np = 3
      DOUBLE COMPLEX NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP),
     &               BUU(NP),BUD(NP),BDU(NP),BDD(NP),IONE,
     &               G(2,2,NP,NP),M(2*NP, 2*NP),zero,
     &               P(2*np,2,NP), psi1(2*np), psi2(2*np), psi3(2*np)
      PARAMETER (IONE=(0.D0,1.D0), zero=(0.d0,0.d0))
      DOUBLE PRECISION U0,U1,N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP)

      DOUBLE COMPLEX VXCMAT_SLATER(2,2,NP),A(2),B(2),RMAT(2,2,2,NP),
     &               VRMAT(2,4,NP),VNMR(2,4,NP),DUM,NM(4,4,NP),
     &               NMR(2,2,2,NP),VXCMAT(2,2,NP)

      U1 = 0.d0
      do i = 1, 2*np
            psi1(i) = m(i, 1)
            psi2(i) = m(i, 2)
            psi3(i) = m(i, 3)
      end do

      DO I=1,NP
         p(1,1,I) = PSI1(I)
         p(1,2,I) = PSI1(I+NP)
         p(2,1,I) = PSI2(I)
         p(2,2,I) = PSI2(I+NP)
         p(3,1,I) = PSI3(I)
         p(3,2,I) = PSI3(I+NP)
      ENDDO    

      DO I=1,2
      DO J=1,2
      DO K=1,NP
      DO L=1,NP
         g(I,J,K,L) = p(1,I,K)*DCONJG(p(1,J,L))
     &                  + p(2,I,K)*DCONJG(p(2,J,L))
     &                  + p(3,I,K)*DCONJG(p(3,J,L))
      ENDDO    
      ENDDO    
      ENDDO    
      ENDDO    

      nuu = zero
      ndd = zero
      nud = zero
      ndu = zero

      DO 1 I=1,NP
            NUU(I) = CDABS(PHI(1,1,I))**2 + CDABS(PHI(2,1,I))**2
            NUD(I) = PHI(1,1,I)*DCONJG(PHI(1,2,I))
     &               + PHI(2,1,I)*DCONJG(PHI(2,2,I))
            NDU(I) = DCONJG(NUD(I))
            NDD(I) = CDABS(PHI(1,2,I))**2 + CDABS(PHI(2,2,I))**2
1     CONTINUE
   

      ! do i = 1, np
      !       DO occ=1,3
      !             NUU(I) = NUU(I) + CDABS(p(occ,1,i))**2
      !             NUD(I) = NUD(I) + p(occ,1,I)*DCONJG(p(occ,2,I))
      !             NDU(I) = NDU(I) + DCONJG(NUD(I))
      !             NDD(I) = ndd(i) + cdabs(p(occ,2,i))**2
      !       end do
      ! end do
   
         DO 2 I=1,NP
            N(I) = DREAL(NUU(I) + NDD(I))
            MX(I) = DREAL(NUD(I) + NDU(I))
            MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
            MZ(I) = DREAL(NUU(I) - NDD(I))
2     CONTINUE
   
         VH(1) = U0*N(1) + U1*N(2)
         DO 3 I=2,NP-1
            VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
3     CONTINUE
         VH(NP) = U0*N(NP) + U1*N(NP-1)
   
         DO 5 I=1,NP
   
         BUU(I) = -2.D0*U0*(NUU(I)*NUU(I) + NUD(I)*NDU(I))
         BDU(I) = -2.D0*U0*(NDU(I)*NUU(I) + NDD(I)*NDU(I))
         BUD(I) = -2.D0*U0*(NUU(I)*NUD(I) + NUD(I)*NDD(I))
         BDD(I) = -2.D0*U0*(NDU(I)*NUD(I) + NDD(I)*NDD(I))
   
         IF (I.LT.NP) THEN
         BUU(I) = BUU(I) -2.D0*U1*(g(1,1,I,I+1)*g(1,1,I+1,I)
     &                         +g(1,2,I,I+1)*g(2,1,I+1,I))
         BDU(I) = BDU(I) -2.D0*U1*(g(2,1,I,I+1)*g(1,1,I+1,I)
     &                         +g(2,2,I,I+1)*g(2,1,I+1,I))
         BUD(I) = BUD(I) -2.D0*U1*(g(1,1,I,I+1)*g(1,2,I+1,I)
     &                         +g(1,2,I,I+1)*g(2,2,I+1,I))
         BDD(I) = BDD(I) -2.D0*U1*(g(2,1,I,I+1)*g(1,2,I+1,I)
     &                         +g(2,2,I,I+1)*g(2,2,I+1,I))
         ENDIF
   
         IF (I.GT.1) THEN
         BUU(I) = BUU(I) -2.D0*U1*(g(1,1,I,I-1)*g(1,1,I-1,I)
     &                         +g(1,2,I,I-1)*g(2,1,I-1,I))
         BDU(I) = BDU(I) -2.D0*U1*(g(2,1,I,I-1)*g(1,1,I-1,I)
     &                         +g(2,2,I,I-1)*g(2,1,I-1,I))
         BUD(I) = BUD(I) -2.D0*U1*(g(1,1,I,I-1)*g(1,2,I-1,I)
     &                         +g(1,2,I,I-1)*g(2,2,I-1,I))
         BDD(I) = BDD(I) -2.D0*U1*(g(2,1,I,I-1)*g(1,2,I-1,I)
     &                         +g(2,2,I,I-1)*g(2,2,I-1,I))
         ENDIF
   
5     CONTINUE
   
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
      
      DO I=1,NP
            VHXC(I) = VH(I) + VXC(I)
      ENDDO
   

      RETURN
      END
