      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER I,J,k,LWORK,INFO,MODE,ITER, vec, krev, xc_guess, xc_slater
      PARAMETER (LWORK = 100)

      DOUBLE PRECISION C,CP,T,TP,V(3),BX(3),BY(3),BZ(3),
     &                 VT(3),BXT(3),BYT(3),BZT(3),VXC(3),
     &                 VHXC(3),BXCX(3),BXCY(3),BXCZ(3),
     &                 VHXCO(3),BXCXO(3),BXCYO(3),BXCZO(3),
     &                 N(3),MX(3),MY(3),MZ(3),E(6),RWORK(100),
     &                 mmat(500, 6, 3), mback(500, 6, 3)
      DOUBLE PRECISION U0,MIX,TOL,EOLD,CRIT,EH,ETOT,EVXC,EXC,EX,EC

      character(100) :: vec_files(3)

      DOUBLE COMPLEX M(6,6),WORK(LWORK),GAMMA(2,2,3,3),PHI(6,6)
      COMMON /EVALS/ E,PHI,GAMMA

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

      xc_slater = 0

      mode = 2

      do vec = 1, 3
        write(vec_files(vec), '(a,i1,a)')
     &              'Site', vec, '-discretePoints_Pdiff_Slater.txt'

        open(1000+vec, file=vec_files(vec))
      end do

      MIX = 0.1D0
      TOL = 1.D-8

      open(111, file='discrete_E-Pdiff_Slater.txt')

      open(1, file='0field_PDiff.txt')

      DO 5 I=1,3
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
         N(I) = 0.D0
         VHXC(I) = 0.D0
         BXCX(I) = 0.D0
         BXCY(I) = 0.D0
         BXCZ(I) = 0.D0
5     CONTINUE

      close(1)

      U0 = 1.d0

        T = 0.6d0

        TP = T

          c = 1.15d0

          cp = c

          EOLD = 0.D0
          ITER = 0

1       CONTINUE
          ITER = ITER + 1
          IF (ITER.GT.100000) then
            write(*,*) 'max iterations in loop #:'
            call exit(-1)
          end if

          DO 3 I=1,3
                VT(I) = V(I) + VHXC(I)
                BXT(I) = BX(I) + BXCX(I)
                BYT(I) = BY(I) + BXCY(I)
                BZT(I) = BZ(I) + BXCZ(I)
3         CONTINUE

          CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

          CALL ZHEEV( 'V', 'U', 6, M, 6, E, WORK, LWORK, RWORK, INFO )

C      DO 10 I=1,6
C10       WRITE(*,*)E(I)
C**----------------------------------------------------------------------
C**   calculate the densities and update the potentials
C**----------------------------------------------------------------------
          CALL GCALC(M,PHI,GAMMA)
          CALL DENCALC(M,N,MX,MY,MZ)

          DO 20 I=1,3
             VHXCO(I) = VHXC(I)
             BXCXO(I) = BXCX(I)
             BXCYO(I) = BXCY(I)
             BXCZO(I) = BXCZ(I)
20        CONTINUE

          if (xc_slater.eq.1) then
            CALL XCPOT_SLATER(U0,0.d0,VXC,VHXC,BXCX,BXCY,BXCZ)
          else
            CALL XCPOT_BALDA(
     &             U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ,EC)
          end if

          DO 21 I=1,3
             VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXCO(I)
             BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCXO(I)
             BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCYO(I)
             BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZO(I)
21      CONTINUE

          IF (DABS((E(1) - EOLD)/E(1)).GT.TOL) THEN
              EOLD = E(1)
              GOTO 1
          ENDIF

          EH = 0.D0
          DO 85 I=1,3
             EH = EH - 0.5D0*U0*N(I)**2
85        CONTINUE

          EVXC = 0.D0
          DO 90 I=1,3
90          EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &           -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)


          CALL EX_CALC(U0,EX)

          EXC = EX + EC
          ETOT = E(1) + E(2) + E(3) + EH + EVXC + EXC

          call time_prop()


          do vec = 1, 3
            mmat(k, 1, vec) = t/u0
            mmat(k, 2, vec) = c/u0
            mmat(k, 3, vec) = real(etot)
            mmat(k, 3+vec, 1) = mx(vec)
            mmat(k, 3+vec, 2) = my(vec)
            mmat(k, 3+vec, 3) = mz(vec)
          end do

      do k = 1, 4
        write(111, *) mmat(k, 3, 1)
      end do

      close(111)

      do vec=1, 3
        close(2000+vec)
      end do

      END
C************************************************************************
C************************************************************************
C************************************************************************
C**---------------------------------------------------------------------
C**   Now start the time propagation with a sudden switch
C**   (read in the new potential and magnetic field)
C**---------------------------------------------------------------------

      Subroutine time_prop(C, T, psi1, psi2)
      IMPLICIT NONE


      OPEN(1, file='td_field.txt')
      DO I=1,NP
         READ(1,*)V(I),BX(I),BY(I),BZ(I)
      ENDDO    
      CLOSE(1)

100   CONTINUE
      TIME = TIME + DT
C      WRITE(*,*)TIME

C      DO 16 I=1,NP
C         BXD(I) = BX(I) + 0.01D0*DSIN(OMEGA*(TIME-DT/2.D0))
C16    CONTINUE

      DO I=1,NP
         VHXC0(I) = VHXC(I)
         BXCX0(I) = BXCX(I)
         BXCY0(I) = BXCY(I)
         BXCZ0(I) = BXCZ(I)
      ENDDO    

!      IF (CORR.EQ.1) THEN
!         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX0,BXCY0,BXCZ0)
!         CALL BCORR(MX,MY,MZ,BXCX0,BXCY0,BXCZ0,TT)
!      ENDIF

      PC_COUNT=-1
120   CONTINUE
      PC_COUNT = PC_COUNT+1

      DO I=1,NP
            VT(I) = V(I) + VHXC(I)
            BTX(I) = BX(I) + BXCX(I)
C            BTX(I) = BXD(I) + BXCX(I)
            BTY(I) = BY(I) + BXCY(I)
            BTZ(I) = BZ(I) + BXCZ(I)
      ENDDO    

      CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

      DO I=1,NP2
        DO J=1,NP2
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
        ENDDO    
      ENDDO    

      DO I=1,NP2
         R(I) = ZERO
         DO J=1,NP2
            R(I) = R(I) + A(I,J)*PSI1(J)
         ENDDO
      ENDDO    

      DO I=1,NP2
        DO J=1,NP2
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
        ENDDO    
      ENDDO    

      CALL ZGESV( NP2, 1, A, NP2, IPIV, R, NP2, INFO )
      DO I=1,NP2
         PSI1P(I) = R(I)
      ENDDO

      DO I=1,NP2
        DO J=1,NP2
            A(I,J) = -0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
        ENDDO    
      ENDDO    

      DO I=1,NP2
         R(I) = ZERO
         DO J=1,NP2
            R(I) = R(I) + A(I,J)*PSI2(J)
         ENDDO
      ENDDO    

      DO I=1,NP2
        DO J=1,NP2
            A(I,J) = 0.5D0*IONE*DT*M(I,J)
            IF (I.EQ.J) A(I,J) = ONE + A(I,J)
        ENDDO    
      ENDDO    

      CALL ZGESV( NP2, 1, A, NP2, IPIV, R, NP2, INFO )
      DO I=1,NP2
         PSI2P(I) = R(I)
      ENDDO

      CALL XCPOT_SLATER(KLI,PSI1P,PSI2P,U0,U1,VHXC1,BXCX1,BXCY1,BXCZ1)

!      IF (CORR.EQ.1) THEN
!         CALL DENCALC(PSI1P,PSI2P,N,MX,MY,MZ,MM,MTOT)
!         CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX1,BXCY1,BXCZ1)
!         CALL BCORR(MX,MY,MZ,BXCX1,BXCY1,BXCZ1,TT)
!      ENDIF

      IF (PC_COUNT.LT.PC) THEN

         DO I=1,NP
            VHXC(I) = 0.5D0*(VHXC0(I) + VHXC1(I))
            BXCX(I) = 0.5D0*(BXCX0(I) + BXCX1(I))
            BXCY(I) = 0.5D0*(BXCY0(I) + BXCY1(I))
            BXCZ(I) = 0.5D0*(BXCZ0(I) + BXCZ1(I))
         ENDDO    
         GOTO 120

      ELSE

         DO I=1,NP2
            PSI1(I) = PSI1P(I)
            PSI2(I) = PSI2P(I)
         ENDDO    
         DO I=1,NP
            VHXC(I) = VHXC1(I)
            BXCX(I) = BXCX1(I)
            BXCY(I) = BXCY1(I)
            BXCZ(I) = BXCZ1(I)
         ENDDO    

      ENDIF

      CALL DENCALC(PSI1,PSI2,N,MX,MY,MZ,MM,MTOT)
      If (L.EQ.2) THEN
C**   project out transverse BXC
      DO I=1,NP
         BXCX(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MX(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCY(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MY(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
         BXCZ(I) = (BXCX1(I)*MX(I) + BXCY1(I)*MY(I) + BXCZ1(I)*MZ(I))
     &             *MZ(I)/(MX(I)**2+MY(I)**2+MZ(I)**2)
      ENDDO    

C      DO 90 I=1,NP
C         BB = DSQRT(BXCX1(I)**2 + BXCY1(I)**2 + BXCZ1(I)**2)
C         BXCX(I) = BB*MX(I)/MM(I)
C         BXCY(I) = BB*MY(I)/MM(I)
C         BXCZ(I) = BB*MZ(I)/MM(I)
C90    CONTINUE
      ENDIF

      CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
      WRITE(LN,999)TIME,N(1),N(2),N(3),N(4)
      WRITE(LN+1,999)TIME,MX(1),MX(2),MX(3),MX(4)
      WRITE(LN+2,999)TIME,MY(1),MY(2),MY(3),MY(4)
      WRITE(LN+3,999)TIME,MZ(1),MZ(2),MZ(3),MZ(4)
      WRITE(LN+4,999)TIME,MM(1),MM(2),MM(3),MM(4)
      WRITE(LN+5,999)TIME,MTOT(1),MTOT(2),MTOT(3),
     &             DSQRT(MTOT(1)**2+MTOT(2)**2+MTOT(3)**2)
      IF (L.EQ.1) THEN
      WRITE(LN+6,999)TIME,TX(1),TX(2),TX(3),TX(4)
      WRITE(LN+7,999)TIME,TY(1),TY(2),TY(3),TY(4)
      WRITE(LN+8,999)TIME,TZ(1),TZ(2),TZ(3),TZ(4)
      WRITE(LN+9,997)TIME,TT(1),TT(2),TT(3)
      ENDIF
      IF (TIME.LT.199.9999D0) GOTO 100




      END SUBROUTINE

C***********************************************************************
C***********************************************************************
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
