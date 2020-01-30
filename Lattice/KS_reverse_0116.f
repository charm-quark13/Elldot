      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER I,J,k,LWORK,INFO,MODE,ITER, vec, krev
      PARAMETER (LWORK = 100)

      DOUBLE PRECISION C,CP,T,TP,V(3),BX(3),BY(3),BZ(3),
     &                 VT(3),BXT(3),BYT(3),BZT(3),VXC(3),
     &                 VHXC(3),BXCX(3),BXCY(3),BXCZ(3),
     &                 VHXCO(3),BXCXO(3),BXCYO(3),BXCZO(3),
     &                 N(3),MX(3),MY(3),MZ(3),E(6),RWORK(100),
     &                 mmat(500, 6, 3), mback(500, 6, 3)
      DOUBLE PRECISION U0,MIX,TOL,EOLD,CRIT,EH,ETOT,EVXC,EXC,EX,EC

      character(30) :: vec_files(3)

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

      do i = 1, 3
        v(i) = 0.d0
        bx(i) = 0.d0
        by(i) = 0.d0
        bz(i) = 0.d0
      end do

      mode = 2

      do vec = 1, 3
        write(vec_files(vec), '(a,i1,a)') 'Site', vec,
     &                                    '-KSRev-VectorField.txt'
        open(1000+vec, file=vec_files(vec))

        write(vec_files(vec), '(a,i1,a)') 'Site', vec,
     &                                    '-KSForward-Vfield.txt'
        open(2000+vec, file=vec_files(vec))

        write(vec_files(vec), '(a,i1,a)') 'Site', vec,
     &                                    '-KSBackward-Vfield.txt'
        open(3000+vec, file=vec_files(vec))

      end do

      MIX = 0.1D0
      TOL = 1.D-8

      DO 5 I=1,3
         !READ(1,*)V(I),BX(I),BY(I),BZ(I)
         N(I) = 0.D0
         VHXC(I) = 0.D0
         BXCX(I) = 0.D0
         BXCY(I) = 0.D0
         BXCZ(I) = 0.D0
5     CONTINUE


      U0 = 1.d0


      DO 101 j=1,100

        T = dble(j)*.01d0

        TP = T

        do k = 1, 500
          if (k.lt.251) then
            c = (250 - k) * 0.01d0
          else
            c = (k - 250) * .01d0
          end if

          cp = c

          EOLD = 0.D0
          ITER = 0

1       CONTINUE
          ITER = ITER + 1
          IF (ITER.GT.100000) call exit(-1)

          DO 3 I=1,3
                VT(I) = V(I) + VHXC(I)
                BXT(I) = BX(I) + BXCX(I)
                BYT(I) = BY(I) + BXCY(I)
C                BZT(I) = BZ(I)
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

C      N(3) = N(1)
C      MZ(3) = MZ(1)

          CALL XCPOT_BALDA(
     &             U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ,EC)

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

!          WRITE(10,*)t/U0, C/U0,real(N(1)),real(N(2)),real(N(3))
!          WRITE(11,*)t/U0, C/U0,real(MX(1)),real(MX(2)),real(MX(3))
!          WRITE(12,*)t/U0, C/U0,real(MY(1)),real(MY(2)),real(MY(3))
!          WRITE(13,*)t/U0, C/U0,real(MZ(1)),real(MZ(2)),real(MZ(3))

          !WRITE(11,*)t/U0, C/U0,real(MX(1)),real(My(1)),real(Mz(1))
          !WRITE(12,*)t/U0, C/U0,real(MX(2)),real(My(2)),real(Mz(2))
          !WRITE(13,*)t/U0, C/U0,real(MX(3)),real(My(3)),real(Mz(3))

!          if (dabs(mz(1)).lt.1d-7.and.
!     &        dabs(mz(2)).lt.1d-7.and.
!     &        dabs(mz(3)).lt.1d-7) then
!              phase = 1

!            write(100+vec,*) t/u0, c/u0, real(mx(vec)),
!     &                       real(my(vec)), real(mz(vec))

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

!          WRITE(15,*)real(U0),real(ETOT)
!
!          WRITE(16,*)real(U0),real(E(1)),real(E(2)),real(E(3)),
!     &                         real(E(4))

          do vec = 1, 3
            mmat(k, 1, vec) = t/u0
            mmat(k, 2, vec) = c/u0
            mmat(k, 3, vec) = real(etot)
            mmat(k, 3+vec, 1) = mx(vec)
            mmat(k, 3+vec, 2) = my(vec)
            mmat(k, 3+vec, 3) = mz(vec)
          end do

          if (k.eq.500) then
            do iter=1,500
              do vec=1,3
                write(2000+vec,*)
     &            mmat(iter, 1, 1),
     &            mmat(iter, 2, 1),
     &            mmat(iter, 3+vec, 1),
     &            mmat(iter, 3+vec, 2),
     &            mmat(iter, 3+vec, 3)
              end do
            end do
          end if

!!************************************************************!!
!!    Beginning of the reverse iteration loop                 !!
!!************************************************************!!

          if (k.eq.500) then
            do krev = 500, 1, -1
              if (krev.lt.251) then
                c = (250 - krev) * 0.01d0
              else
                c = (krev - 250) * .01d0
              end if

              cp = c

              EOLD = 0.D0
              ITER = 0

10            CONTINUE
              ITER = ITER + 1
              IF (ITER.GT.100000) call exit(-1)

              DO 30 I=1,3
                    VT(I) = V(I) + VHXC(I)
                    BXT(I) = BX(I) + BXCX(I)
                    BYT(I) = BY(I) + BXCY(I)
C                   BZT(I) = BZ(I)
                    BZT(I) = BZ(I) + BXCZ(I)
30             CONTINUE

              CALL MATRIX(M,C,CP,T,TP,VT,BXT,BYT,BZT)

              CALL ZHEEV( 'V', 'U', 6, M, 6, E,
     &                    WORK, LWORK, RWORK, INFO )

C      DO 10 I=1,6
C10       WRITE(*,*)E(I)
C**----------------------------------------------------------------------
C**   calculate the densities and update the potentials
C**----------------------------------------------------------------------
              CALL GCALC(M,PHI,GAMMA)
              CALL DENCALC(M,N,MX,MY,MZ)

              DO 200 I=1,3
                 VHXCO(I) = VHXC(I)
                 BXCXO(I) = BXCX(I)
                 BXCYO(I) = BXCY(I)
                 BXCZO(I) = BXCZ(I)
200           CONTINUE

C      N(3) = N(1)
C      MZ(3) = MZ(1)

              CALL XCPOT_BALDA(
     &              U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ,EC)

              DO 210 I=1,3
                 VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXCO(I)
                 BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCXO(I)
                 BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCYO(I)
                 BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZO(I)
210            CONTINUE

              IF (DABS((E(1) - EOLD)/E(1)).GT.TOL) THEN
                  EOLD = E(1)
                  GOTO 10
              ENDIF

              !WRITE(11,*)t/U0, C/U0,real(MX(1)),real(My(1)),real(Mz(1))
              !WRITE(12,*)t/U0, C/U0,real(MX(2)),real(My(2)),real(Mz(2))
              !WRITE(13,*)t/U0, C/U0,real(MX(3)),real(My(3)),real(Mz(3))

    !          if (dabs(mz(1)).lt.1d-7.and.
    !     &        dabs(mz(2)).lt.1d-7.and.
    !     &        dabs(mz(3)).lt.1d-7) then
    !              phase = 1

              EH = 0.D0
              DO 850 I=1,3
                 EH = EH - 0.5D0*U0*N(I)**2
850            CONTINUE

              EVXC = 0.D0
              DO 900 I=1,3
900               EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &           -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)

              CALL EX_CALC(U0,EX)

              EXC = EX + EC
              ETOT = E(1) + E(2) + E(3) + EH + EVXC + EXC

!              WRITE(15,*)real(U0),real(ETOT)
!
!              WRITE(16,*)real(U0),real(E(1)),real(E(2)),real(E(3)),
!     &                         real(E(4))

              do vec = 1, 3
                mback(krev, 1, vec) = t/u0
                mback(krev, 2, vec) = c/u0
                mback(krev, 3, vec) = real(etot)
                mback(krev, 3+vec, 1) = mx(vec)
                mback(krev, 3+vec, 2) = my(vec)
                mback(krev, 3+vec, 3) = mz(vec)
              end do


              if (dabs(etot).lt.dabs(mmat(krev, 3, 1))) then
                do vec = 1, 3
                  mmat(k, 1, vec) = t/u0
                  mmat(k, 2, vec) = c/u0
                  mmat(k, 3, vec) = real(etot)
                  mmat(k, 3+vec, 1) = real(mx(vec))
                  mmat(k, 3+vec, 2) = real(my(vec))
                  mmat(k, 3+vec, 3) = real(mz(vec))
                end do
              end if

            end do

            do iter=1,500
              do vec=1,3
                write(3000+vec,*)
     &            mback(iter, 1, 1),
     &            mback(iter, 2, 1),
     &            mback(iter, 3+vec, 1),
     &            mback(iter, 3+vec, 2),
     &            mback(iter, 3+vec, 3)
              end do
            end do

          end if

        end do

        do k=1,500
          do vec=1,3
            write(1000+vec,*)
     &        mmat(k, 1, 1),
     &        mmat(k, 2, 1),
     &        mmat(k, 3+vec, 1),
     &        mmat(k, 3+vec, 2),
     &        mmat(k, 3+vec, 3)
          end do
        end do

101   CONTINUE

      do vec=1, 3
        close(100+vec)
      end do

      END
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
