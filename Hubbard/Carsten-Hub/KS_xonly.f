      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER NP,I,ITER,REST,CORR,LWORK,INFO,it,xc,bmag,
     &        wmflg,wtflg,weflg
      INCLUDE 'dim.inc'
      PARAMETER (LWORK = 100)

      DOUBLE PRECISION T,TOL,EOLD,U1,U0,MIX,EX,CRIT,
     &                 EH,EVXC,EXC,ETOT,TX(NP),TY(NP),TZ(NP),TT(3),
     &                 VHXC(NP),VXC(NP),BXCX(NP),BXCY(NP),
     &                 BXCZ(NP),VHXCO(NP),BXCXO(NP),BXCYO(NP),BXCZO(NP),
     &                 V(NP),BX(NP),BY(NP),BZ(NP),
     &                 VT(NP),BTX(NP),BTY(NP),BTZ(NP),
     &                 N(NP),MX(NP),MY(NP),MZ(NP),E(2*NP),RWORK(100),
     &                 tplot(np+1)

      DOUBLE COMPLEX M(2*NP,2*NP),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               WORK(LWORK)

      character(30) :: mwn, mwx, mwy, mwz, tprint, twx, twy, twz

      COMMON /EVALS/ E,PHI,GAMMA

      write(tprint,'(a, i3, a)') '(', np + 1, 'e16.6)'

      xc = 2
      corr = 0

      wmflg = 1
      wtflg = 1
      weflg = 1

      do bmag = 1, 11

        v = 0.d0
        Bx = 0.d0
        by = 0.d0
        bz = 0.d0

        v(1) = 1.d0
        do i = 2, 3
          v(i) = -1.d0
        end do
        v(4) = 1.d0

        Bx(1) = dble(bmag/10.d0)
        Bx(2) = dble(bmag/10.d0)
!        Bx(3) = .01d0
!        Bx(4) = -.25d0

!        Bx = 0.d0
        By(2) = -dble(bmag/10.d0)
        By(3) = -dble(bmag/10.d0)

!        Bz(1) = -2.50d-1
!        Bz(2) = 2.50d-2
        Bz(3) = dble(bmag/10.d0)
        Bz(4) = dble(bmag/10.d0)

        if (bmag.eq.11) then

          Bx(1) = 10.d-3
          Bx(2) = 10.d-3

          By(2) = -10.d-3
          By(3) = -10.d-3

          Bz(3) = 10.d-3
          Bz(4) = 10.d-3

        end if

        if (wmflg.eq.1) then
          if (bmag.lt.10) then
            write(mwn,'(a,i1,a)') '4pt-B', bmag, '-OEP-n.txt'
            write(mwx,'(a,i1,a)') '4pt-B', bmag, '-OEP-mx.txt'
            write(mwy,'(a,i1,a)') '4pt-B', bmag, '-OEP-my.txt'
            write(mwz,'(a,i1,a)') '4pt-B', bmag, '-OEP-mz.txt'
          elseif (bmag.eq.10) then
            write(mwn,'(a,i2,a)') '4pt-B', bmag, '-OEP-n.txt'
            write(mwx,'(a,i2,a)') '4pt-B', bmag, '-OEP-mx.txt'
            write(mwy,'(a,i2,a)') '4pt-B', bmag, '-OEP-my.txt'
            write(mwz,'(a,i2,a)') '4pt-B', bmag, '-OEP-mz.txt'
          elseif (bmag.eq.11) then
            write(mwn,'(a)') '4pt-B01-OEP-n.txt'
            write(mwx,'(a)') '4pt-B01-OEP-mx.txt'
            write(mwy,'(a)') '4pt-B01-OEP-my.txt'
            write(mwz,'(a)') '4pt-B01-OEP-mz.txt'
          end if

          open(100,file = mwn)
          open(101,file = mwx)
          open(102,file = mwy)
          open(103,file = mwz)
        end if

        if (wtflg.eq.1) then
          if (bmag.lt.10) then
            write(twx,'(a,i1,a)') '4pt-B', bmag, '-OEP-tx.txt'
            write(twy,'(a,i1,a)') '4pt-B', bmag, '-OEP-ty.txt'
            write(twz,'(a,i1,a)') '4pt-B', bmag, '-OEP-tz.txt'
          elseif (bmag.eq.10) then
            write(twx,'(a,i2,a)') '4pt-B', bmag, '-OEP-tx.txt'
            write(twy,'(a,i2,a)') '4pt-B', bmag, '-OEP-ty.txt'
            write(twz,'(a,i2,a)') '4pt-B', bmag, '-OEP-tz.txt'
          elseif (bmag.eq.11) then
            write(twx,'(a)') '4pt-B01-OEP-tx.txt'
            write(twy,'(a)') '4pt-B01-OEP-ty.txt'
            write(twz,'(a)') '4pt-B01-OEP-tz.txt'
          end if

          open(201,file = twx)
          open(202,file = twy)
          open(203,file = twz)
        end if

        if (weflg.eq.1) then
          if (bmag.lt.10) then
            write(twx,'(a,i1,a)') '4pt-B', bmag, '-OEP-E.txt'
          elseif (bmag.eq.10) then
            write(twx,'(a,i2,a)') '4pt-B', bmag, '-OEP-E.txt'
          elseif (bmag.eq.11) then
            write(twx,'(a)') '4pt-B01-OEP-E.txt'
          end if

          open(300,file = twx)
        end if

        do it = 1, 100

        if (it.ge.2) then
          rest = 1
        else
          rest = 0
        end if

        T = 0.5D0
        U0 = dble(it)/10.d0
        U1 = U0/2.D0

        MIX = 0.1D0
        TOL = 1.D-10
        EOLD = 0.D0
        ITER = 0

        IF (REST.EQ.0) THEN
           DO 2 I=1,NP
              VHXC(I) = 0.D0
              BXCX(I) = 0.D0
              BXCY(I) = 0.D0
              BXCZ(I) = 0.D0
2        CONTINUE
        ELSEIF (REST.EQ.1) THEN
           READ(2)VHXC,BXCX,BXCY,BXCZ
           REWIND 2
        ENDIF

1     CONTINUE
        ITER = ITER + 1

        IF (ITER.GT.500000) STOP

        DO 3 I=1,NP
            VT(I) = V(I) + VHXC(I)
            BTX(I) = BX(I) + BXCX(I)
            BTY(I) = BY(I) + BXCY(I)
            BTZ(I) = BZ(I) + BXCZ(I)
3     CONTINUE

        CALL MATRIX(M,T,VT,BTX,BTY,BTZ)

        CALL ZHEEV( 'V', 'U', 2*NP, M, 2*NP, E, WORK, LWORK,
     &               RWORK, INFO )

!        DO 10 I=1,2*NP
!10       WRITE(*,*)E(I)
C**----------------------------------------------------------------------
C**   calculate the densities and update the potentials
C**----------------------------------------------------------------------
        CALL GCALC(M,PHI,GAMMA)
        CALL DENCALC(M,N,MX,MY,MZ)

        DO 20 I=1,NP
           VHXCO(I) = VHXC(I)
           BXCXO(I) = BXCX(I)
           BXCYO(I) = BXCY(I)
           BXCZO(I) = BXCZ(I)
20    CONTINUE

      IF (XC.EQ.1) THEN
         CALL XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      ELSEIF (XC.EQ.2) THEN
         CALL XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      ENDIF

C**   symmetrization, if needed
        VHXC(3)=VHXC(2)
        VHXC(4)=VHXC(1)
        BXCZ(1) = BXCX(4)
        BXCZ(2) = BXCX(3)
        BXCZ(3) = BXCX(2)
        BXCZ(4) = BXCX(1)

        IF (corr.eq.1) then
          CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
        end if

C**   do the xc torque correction

!        IF (CORR.EQ.1) THEN
!           WRITE(*,*)'total torque:'
!           WRITE(*,*)TT
!           WRITE(*,*)
!           CALL BCORR(MX,MY,MZ,BXCX,BXCY,BXCZ,TT)
!           CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
!           WRITE(*,*)'total torque:'
!           WRITE(*,*)TT
!           WRITE(*,*)
!        ENDIF

        DO 21 I=1,NP
           VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXCO(I)
           BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCXO(I)
           BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCYO(I)
           BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZO(I)
21    CONTINUE

        CRIT = 0.D0
        DO 25 I=1,NP
!25       CRIT = CRIT + N(I) + DABS(MX(I)) + DABS(MY(I)) + DABS(MZ(I))
25       CRIT = CRIT + dsqrt(dble(MX(I)**2 + MY(I)**2 + MZ(I)**2))
        IF (DABS((CRIT- EOLD)/CRIT).GT.TOL) THEN
            EOLD = CRIT
            GOTO 1
        ENDIF

        WRITE(*,*)
        WRITE(*,*)'************',ITER,'*************'
        WRITE(*,*)

        open(2,form='unformatted')

        WRITE(2)VHXC,BXCX,BXCY,BXCZ

        close(2)


        !write(*,*)'       N             MX               MY         MZ'

!        DO 80 I=1,NP
!80    write(*,*)real(N(I)),real(MX(I)),real(MY(I)),real(MZ(I))

        !write(*,*)
        !write(*,*)'magnetization magnitude:'
!        DO 81 I=1,NP
!81    write(*,*)sqrt(real(MX(I)**2 + MY(I)**2 + MZ(I)**2))

        !write(*,*)
        !write(*,*)'      VXC         BXCX           BXCY           BXCZ'

!        DO 82 I=1,NP
!82    write(*,*)real(VXC(I)),real(BXCX(I)),real(BXCY(I)),real(BXCZ(I))

        !write(*,*)
        !write(*,*)'=====================GS Energy====================='

        EH = 0.D0
        DO 85 I=1,NP
           EH = EH - 0.5D0*U0*N(I)**2
85    CONTINUE
        DO 86 I=1,NP-1
           EH = EH - U1*N(I)*N(I+1)
86    CONTINUE

        EVXC = 0.D0
        DO 90 I=1,NP
90       EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &         -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)

        CALL EX_CALC(U0,U1,EX)

        EXC = EX

        ETOT = E(1) + E(2) + EH + EVXC + EXC

        !write(*,*)' EKS = ',E(1) + E(2)
        !write(*,*)'  EH = ',EH
        !write(*,*)'EVXC = ',EVXC
        !write(*,*)' EXC = ',EXC
        !write(*,*)
        !write(*,*)'ETOT = ',real(ETOT)

        !write(*,*)
        !write(*,*)'(xc torque)_x    (xc torque)_y    (xc torque)_z'
!        DO 95 I=1,NP
!95       write(*,*)real(TX(I)),' ',real(TY(I)),' ',real(TZ(I))
        !write(*,*)'-----------------------------------------'
        !write(*,*)real(TT(1)),' ',real(TT(2)),' ',real(TT(3))
        !write(*,*)

        tplot(1) = u0

        if (wmflg.eq.1) then
          do i=1, np
            tplot(i+1) = n(i)
          end do
          write(100,tprint) tplot

          do i=1, np
            tplot(i+1) = mx(i)
          end do
          write(101,tprint) tplot

          do i=1, np
            tplot(i+1) = my(i)
          end do
          write(102,tprint) tplot

          do i=1, np
            tplot(i+1) = mz(i)
          end do
          write(103,tprint) tplot
        end if

        if (wtflg.eq.1) then
          do i=1, np
            tplot(i+1) = tx(i)
          end do
          write(201,tprint) tplot

          do i=1, np
            tplot(i+1) = ty(i)
          end do
          write(202,tprint) tplot

          do i=1, np
            tplot(i+1) = tz(i)
          end do
          write(203,tprint) tplot
        end if

        if (weflg.eq.1) then
          do i=1, np
            tplot(i+1) = E(i)
!            tplot(i+1) = Etot
          end do
          write(300,tprint) tplot
        end if

        write(*,*) '************************************'
        write(*,*) it

        end do

        if (wmflg.eq.1) then
          close(100)
          close(101)
          close(102)
          close(103)
        end if

        if (wtflg.eq.1) then
          close(201)
          close(202)
          close(203)
        end if

      end do

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE BCORR (MX,MY,MZ,BX,BY,BZ,TT)
      IMPLICIT NONE

      INTEGER I,J,II,JJ,NP,NM,INFO,COP
      INCLUDE 'dim.inc'
      PARAMETER (NM = 3*NP-3)
      INTEGER IPIV(NM),LWORK
      PARAMETER (LWORK=100)
      DOUBLE PRECISION MX(NP),MY(NP),MZ(NP),BX(NP),BY(NP),BZ(NP),
     &                 MAT(NM,NM),RVEC(NM),WORK(LWORK),
     &                 BBX(NP),BBY(NP),BBZ(NP)
      DOUBLE PRECISION AX(NP),AY(NP),AZ(NP),AD,MSX,MSY,MSZ,M2,TT(3),
     &                 PX(NP),PY(NP),PZ(NP),QX(NP),QY(NP),QZ(NP)

      COP=0
C**--------------------------------------------------------------------***
C**   Do we have a coplanar situation?
C**   COP=1,2,3: torque along x,y,z
C**--------------------------------------------------------------------***
      MSX = 0.D0
      MSY = 0.D0
      MSZ = 0.D0
      M2 = 0.D0
      DO 1 I=1,NP
         MSX = MSX + DABS(MX(I))
         MSY = MSY + DABS(MY(I))
         MSZ = MSZ + DABS(MZ(I))
         M2 = M2 + MX(I)**2 + MY(I)**2 + MZ(I)**2
1     CONTINUE

      IF (MSX.LT.1.D-10) COP=1
      IF (MSY.LT.1.D-10) COP=2
      IF (MSZ.LT.1.D-10) COP=3

      IF (COP.EQ.0) THEN
C**--------------------------------------------------------------------***
C**   Not coplanar, so we have to solve a system of equations
C**--------------------------------------------------------------------***

      DO 10 I=1,NP
         AD = MY(NP)*MZ(1) - MZ(NP)*MY(1)
         AX(I) = -(MY(NP)*MZ(I) - MZ(NP)*MY(I))/AD
         AY(I) = -(MZ(NP)*MX(I) - MX(NP)*MZ(I))/AD
         AZ(I) = -(MX(NP)*MY(I) - MY(NP)*MX(I))/AD
10    CONTINUE

      DO 11 I=1,NP
         PX(I) = MY(1)*AX(I) + MY(I)
         PY(I) = MY(1)*AY(I) - MX(I)
         PZ(I) = MY(1)*AZ(I)

         QX(I) = MZ(1)*AX(I) + MZ(I)
         QY(I) = MZ(1)*AY(I)
         QZ(I) = MZ(1)*AZ(I) - MX(I)
11    CONTINUE

      DO 20 I=1,NP-1
         II = 3*(I-1) + 1

         RVEC(II)   = -MX(NP)**2*BY(I) - MX(NP)**2*AY(I)*BX(1)
     &              - (MY(1)*AY(I)-MX(I))*MX(NP)*BY(NP)
     &              - MZ(1)*AY(I)*MX(NP)*BZ(NP)

         RVEC(II+1) = -MX(NP)**2*BZ(I) - MX(NP)**2*AZ(I)*BX(1)
     &              - MY(1)*AZ(I)*MX(NP)*BY(NP)
     &              - (MZ(1)*AZ(I)-MX(I))*MX(NP)*BZ(NP)

         RVEC(II+2) = -MX(NP)**2*BX(I+1) - MX(NP)**2*AX(I+1)*BX(1)
     &              - (MY(1)*AX(I+1)+MY(I+1))*MX(NP)*BY(NP)
     &              - (MZ(1)*AX(I+1)+MZ(I+1))*MX(NP)*BZ(NP)

20    CONTINUE

      DO 30 I=1,NP-1
      DO 30 J=1,NP-1
         II = 3*(I-1) + 1
         JJ = 3*(J-1) + 1

         MAT(II,JJ) = MX(NP)**2*AY(I)*AY(J)
     &              + (MY(1)*AY(I)-MX(I))*PY(J)
     &              + MZ(1)*AY(I)*QY(J)

         MAT(II,JJ+1) = MX(NP)**2*AY(I)*AZ(J)
     &              + (MY(1)*AY(I)-MX(I))*PZ(J)
     &              + MZ(1)*AY(I)*QZ(J)

         MAT(II,JJ+2) = MX(NP)**2*AY(I)*AX(J+1)
     &              + (MY(1)*AY(I)-MX(I))*PX(J+1)
     &              + MZ(1)*AY(I)*QX(J+1)

         MAT(II+1,JJ) = MX(NP)**2*AZ(I)*AY(J)
     &              + MY(1)*AZ(I)*PY(J)
     &              + (MZ(1)*AZ(I)-MX(I))*QY(J)

         MAT(II+1,JJ+1) = MX(NP)**2*AZ(I)*AZ(J)
     &              + MY(1)*AZ(I)*PZ(J)
     &              + (MZ(1)*AZ(I)-MX(I))*QZ(J)

         MAT(II+1,JJ+2) = MX(NP)**2*AZ(I)*AX(J+1)
     &              + MY(1)*AZ(I)*PX(J+1)
     &              + (MZ(1)*AZ(I)-MX(I))*QX(J+1)

         MAT(II+2,JJ) = MX(NP)**2*AX(I+1)*AY(J)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PY(J)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QY(J)

         MAT(II+2,JJ+1) = MX(NP)**2*AX(I+1)*AZ(J)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PZ(J)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QZ(J)

         MAT(II+2,JJ+2) = MX(NP)**2*AX(I+1)*AX(J+1)
     &              + (MY(1)*AX(I+1)+MY(I+1))*PX(J+1)
     &              + (MZ(1)*AX(I+1)+MZ(I+1))*QX(J+1)

30    CONTINUE

      DO 35 I=1,NM
         MAT(I,I) = MAT(I,I) + MX(NP)**2
35    CONTINUE

      CALL DSYSV('L',NM,1,MAT,NM,IPIV,RVEC,NM,WORK,LWORK,INFO )

      DO 40 I=1,NP-1
         II = 3*(I-1)+1
         BBY(I) = RVEC(II)
         BBZ(I) = RVEC(II+1)
         BBX(I+1) = RVEC(II+2)
40    CONTINUE

      BBX(1) = 0.D0
      DO 45 I=1,NP-1
         BBX(1) = BBX(1) + AY(I)*BBY(I) + AZ(I)*BBZ(I)
     &                   + AX(I+1)*BBX(I+1)
45    CONTINUE

      BBY(NP) = MY(NP)*BBX(NP)/MX(NP)
      BBZ(NP) = MZ(NP)*BBX(NP)/MX(NP)
      DO 50 I=1,NP-1
         BBY(NP) = BBY(NP) + MY(I)*BBX(I)/MX(NP)
     &                     - MX(I)*BBY(I)/MX(NP)
         BBZ(NP) = BBZ(NP) + MZ(I)*BBX(I)/MX(NP)
     &                     - MX(I)*BBZ(I)/MX(NP)
50    CONTINUE

      DO 100 I=1,NP
         BX(I) = - BBX(I)
         BY(I) = - BBY(I)
         BZ(I) = - BBZ(I)
100   CONTINUE

      ELSE
C**--------------------------------------------------------------------***
C**   coplanar, so we have an explicit solution
C**--------------------------------------------------------------------***
      DO 200 I=1,NP
         BX(I) = BX(I) - (TT(2)*MZ(I) - TT(3)*MY(I))/M2
         BY(I) = BY(I) - (TT(3)*MX(I) - TT(1)*MZ(I))/M2
         BZ(I) = BZ(I) - (TT(1)*MY(I) - TT(2)*MX(I))/M2
200   CONTINUE

      ENDIF

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER I,NP
      INCLUDE 'dim.inc'

      DOUBLE PRECISION TT(3),TX(NP),TY(NP),TZ(NP),MX(NP),MY(NP),
     &                 MZ(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

      TT(1)=0.D0
      TT(2)=0.D0
      TT(3)=0.D0
      DO 10 I=1,NP
         TX(I) = MY(I)*BXCZ(I) - MZ(I)*BXCY(I)
         TY(I) = MZ(I)*BXCX(I) - MX(I)*BXCZ(I)
         TZ(I) = MX(I)*BXCY(I) - MY(I)*BXCX(I)
         TT(1) = TT(1) + TX(I)
         TT(2) = TT(2) + TY(I)
         TT(3) = TT(3) + TZ(I)
10    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE MATRIX(M,T,V,BX,BY,BZ)
      IMPLICIT NONE

      INTEGER I,J,NP
      INCLUDE 'dim.inc'
      DOUBLE PRECISION T,V(NP),BX(NP),BY(NP),BZ(NP)
      DOUBLE COMPLEX M(2*NP,2*NP),ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      DO 1 I=1,2*NP
      DO 1 J=1,2*NP
1        M(I,J) = ZERO

      DO 10 I=1,NP
         M(I,I) = (V(I) + BZ(I))*ONE
         M(I,I+NP) = BX(I)*ONE - BY(I)*IONE
         M(I+NP,I+NP) = (V(I) - BZ(I))*ONE
         M(I+NP,I) = BX(I)*ONE + BY(I)*IONE
10    CONTINUE

      DO 20 I=1,NP-1
         M(I,I+1) = -T*ONE
         M(I+NP,I+NP+1) = -T*ONE
         M(I+1,I) = -T*ONE
         M(I+NP+1,I+NP) = -T*ONE
20    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE DENCALC(M,N,MX,MY,MZ)
      IMPLICIT NONE

      INTEGER I,NP
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX M(2*NP,2*NP),NUU(NP),NUD(NP),NDU(NP),NDD(NP)
      DOUBLE PRECISION N(NP),MX(NP),MY(NP),MZ(NP)

      DO 5 I=1,NP
         NUU(I) = CDABS(M(I,1))**2 + CDABS(M(I,2))**2
         NUD(I) = M(I,1)*DCONJG(M(I+NP,1)) + M(I,2)*DCONJG(M(I+NP,2))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(M(I+NP,1))**2 + CDABS(M(I+NP,2))**2
5     CONTINUE

      DO 10 I=1,NP
         N(I) = DREAL(NUU(I) + NDD(I))
         MX(I) = DREAL(NUD(I) + NDU(I))
         MY(I) = DREAL((0.D0,1.D0)*(NUD(I) - NDU(I)))
         MZ(I) = DREAL(NUU(I) - NDD(I))
10    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE GCALC(M,PHI,GAMMA)
      IMPLICIT NONE

      INTEGER I,J,K,L,NP
      INCLUDE 'dim.inc'
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
10    CONTINUE

      RETURN
      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,I
      INCLUDE 'dim.inc'
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

      SUBROUTINE XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,ND,I,J,K,MU,NU,ALPHA,BETA,R,RP,INFO,LWORK
      INCLUDE 'dim.inc'
      PARAMETER (ND=4*NP,LWORK=100)
      DOUBLE PRECISION SING(ND)
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),DUM,
     &               LM(ND,ND),BMAT(2,2*NP),RVEC(ND),VXCMAT(2,2,NP),
     &               X(ND),X1(ND),UMAT(ND,ND),VTMAT(ND,ND),WORK(LWORK)
      DOUBLE PRECISION U0,U1,E(2*NP),N(NP),VH(NP),VXC(NP),RWORK(100),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP)

      COMMON /EVALS/ E,PHI,GAMMA

      DO 1 I=1,NP
1        N(I) = DREAL(GAMMA(1,1,I,I) + GAMMA(2,2,I,I))

      VH(1) = U0*N(1) + U1*N(2)
      DO 2 I=2,NP-1
         VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
2     CONTINUE
      VH(NP) = U0*N(NP) + U1*N(NP-1)

C**   calculate the (ND x ND) OEP matrix

      DO 4 MU=1,2
      DO 4 NU=1,2
      DO 4 ALPHA=1,2
      DO 4 BETA=1,2
      DO 4 R=1,NP
      DO 4 RP=1,NP
         DUM = (0.D0,0.D0)
         DO 5 I=1,2
         DO 5 J=1,2*NP
            IF (I.NE.J) THEN
               DUM = DUM + PHI(I,BETA,RP)*DCONJG(PHI(J,ALPHA,RP))
     &               *DCONJG(PHI(I,MU,R))*PHI(J,NU,R)/(E(J)-E(I))
     &                   + DCONJG(PHI(I,ALPHA,RP))*PHI(J,BETA,RP)
     &               *PHI(I,NU,R)*DCONJG(PHI(J,MU,R))/(E(J)-E(I))
            ENDIF
5        CONTINUE
         LM(MU+(NU-1)*2+(R-1)*4,ALPHA+(BETA-1)*2+(RP-1)*4) = DUM
4     CONTINUE

C**   calculate the OEP right-hand side

      DO 10 I=1,2
      DO 10 J=1,2*NP
         DUM = (0.D0,0.D0)
         IF (I.NE.J) THEN
            DO 11 ALPHA=1,2
            DO 11 BETA=1,2
               DO 12 K=1,NP
                  DUM = DUM + (U0/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K))
     &                          *GAMMA(ALPHA,BETA,K,K)*PHI(J,BETA,K)
12             CONTINUE
               DO 13 K=1,NP-1
                  DUM = DUM + (U1/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K))
     &                        *GAMMA(ALPHA,BETA,K,K+1)*PHI(J,BETA,K+1)
     &                      + (U1/(E(I)-E(J)))*DCONJG(PHI(I,ALPHA,K+1))
     &                        *GAMMA(ALPHA,BETA,K+1,K)*PHI(J,BETA,K)
13             CONTINUE
11          CONTINUE
         ENDIF
         BMAT(I,J) = DUM
10    CONTINUE

      DO 15 MU=1,2
      DO 15 NU=1,2
      DO 15 R=1,NP
         DUM = (0.D0,0.D0)
         DO 16 I=1,2
         DO 16 J=1,2*NP
            DUM = DUM + DCONJG(BMAT(I,J)*PHI(I,MU,R))*PHI(J,NU,R)
     &                + BMAT(I,J)*DCONJG(PHI(J,MU,R))*PHI(I,NU,R)
16       CONTINUE
         RVEC(MU+(NU-1)*2+(R-1)*4) = DUM
15    CONTINUE

C**   now do the singular value decomposition

      CALL ZGESVD( 'A', 'A', ND, ND, LM, ND, SING, UMAT, ND, VTMAT, ND,
     &                   WORK, LWORK, RWORK, INFO )

      DO 70 I=1,ND
         DUM = (0.D0,0.D0)
         DO 71 J=1,ND
            DUM = DUM + DCONJG(UMAT(J,I))*RVEC(J)
71       CONTINUE
         X1(I) = DUM
70    CONTINUE

!      WRITE(*,*)
!      WRITE(*,*)'Lowest OEP singular values:',real(SING(ND)),
!     &          real(SING(ND-1)),real(SING(ND-2))
!      WRITE(*,*)

C**   The singular values are ordered from largest to smallest.
C**   By default, we drop the last singular value SING(ND).
C**   But for strong correlations (large U), it may happen that
C**   other singular values need to be dropped as well. This needs
C**   to be checked from case to case, and, if necessary, the code
C**   below needs to be modified.

      DO 72 I=1,ND
!         IF (I.LT.ND) THEN
         IF (SING(I).lt.1d-11) then
!         IF (I.LT.ND-1) THEN
            SING(I) = 0.D0
         ELSE
            SING(I) = 1.D0/SING(I)
         ENDIF
         X1(I) = SING(I)*X1(I)
72    CONTINUE

      DO 75 I=1,ND
         DUM = (0.D0,0.D0)
         DO 76 J=1,ND
            DUM = DUM + DCONJG(VTMAT(J,I))*X1(J)
76       CONTINUE
         X(I) = DUM
75    CONTINUE

      DO 80 MU=1,2
      DO 80 NU=1,2
      DO 80 R=1,NP
         VXCMAT(MU,NU,R) = X(MU+(NU-1)*2+(R-1)*4)
80    CONTINUE

      DO 90 R=1,NP
         VXC(R) = DREAL(VXCMAT(1,1,R) + VXCMAT(2,2,R))/2.D0
         VHXC(R) = VH(R) + VXC(R)
         BXCX(R) = DREAL(VXCMAT(1,2,R) + VXCMAT(2,1,R))/2.D0
         BXCY(R) = DREAL((0.D0,1.D0)*(VXCMAT(1,2,R) - VXCMAT(2,1,R)))
     &             /2.D0
         BXCZ(R) = DREAL(VXCMAT(1,1,R) - VXCMAT(2,2,R))/2.D0
90    CONTINUE

      RETURN
      END

C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE EX_CALC(U0,U1,EX)
      IMPLICIT NONE

      INTEGER NP,I,TAU,SIGMA
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)
      DOUBLE PRECISION E(2*NP),U0,U1,EX

      COMMON /EVALS/ E,PHI,GAMMA

      EX = 0.D0
      DO 10 I=1,NP
      DO 10 TAU=1,2
      DO 10 SIGMA=1,2
         EX = EX -0.5D0*U0*GAMMA(SIGMA,TAU,I,I)*GAMMA(TAU,SIGMA,I,I)
10    CONTINUE

      DO 20 I=1,NP-1
      DO 20 TAU=1,2
      DO 20 SIGMA=1,2
         EX = EX - U1*GAMMA(SIGMA,TAU,I,I+1)*GAMMA(TAU,SIGMA,I+1,I)
20    CONTINUE

      RETURN
      END
