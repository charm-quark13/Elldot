      PROGRAM TWOSPIN
      IMPLICIT NONE

      INTEGER NP,I,ITER,REST,CORR,LWORK,INFO,it,xc,bmag,
     &        wmflg,wtflg,weflg,nd,revit,j,k,mp2,qrev,bxcflg
      INCLUDE 'dim.inc'
      PARAMETER (LWORK = 100, ND = 4*NP)

      DOUBLE PRECISION T,TOL,EOLD,U1,U0,MIX,EX,CRIT,
     &                 EH,EVXC,EXC,ETOT,TX(NP),TY(NP),TZ(NP),TT(3),
     &                 VHXC(NP),VXC(NP),BXCX(NP),BXCY(NP),
     &                 BXCZ(NP),VHXCO(NP),BXCXO(NP),BXCYO(NP),BXCZO(NP),
     &                 V(NP),BX(NP),BY(NP),BZ(NP),
     &                 VT(NP),BTX(NP),BTY(NP),BTZ(NP),
     &                 N(NP),MX(NP),MY(NP),MZ(NP),E(2*NP),RWORK(100),
     &                 tplot(np+1),sing(nd),eback, etplot(np+2),
     &                 emplot(np+2),dl,ec,excl,lambda, explot(np+2),
     &                 mdotb(np+2)

      DOUBLE PRECISION etmats(101,np+2,3), emmats(101,np+2,4),
     &                 exmats(101,np+2,4), mbmats(101,np+2)

      DOUBLE COMPLEX M(2*NP,2*NP),GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               WORK(LWORK)

      character(30) :: mwn, mwx, mwy, mwz, tprint, twx, twy, twz,
     &                 etprint

      COMMON /EVALS/ E,PHI,GAMMA

      write(tprint,'(a, i3, a)') '(', np + 1, 'e16.6)'
      write(etprint,'(a, i2, a)') '(', np + 2, 'e16.6)'

      xc = 1
      corr = 0
      rest = 0
      mp2 = 0

      if (xc.lt.4) then
        qrev = 1
      else
        qrev = 0
      end if

      wmflg = 1
      wtflg = 1
      weflg = 1
      bxcflg = 1

      do bmag = 1, 1

        etmats = 0.d0
        emmats = 0.d0
        exmats = 0.d0
        mbmats = 0.d0

        v = 0.d0
        Bx = 0.d0
        by = 0.d0
        bz = 0.d0

        v(1) = 1.d0
        do i = 2, 3
          v(i) = -1.d0
        end do
        v(4) = 1d0


        Bx(1) = 2d-1
        Bx(2) = dble(bmag/10.d0)

        Bz(3) = 3d-1
        Bz(4) = -2d-1


        if (bmag.eq.11) then

          Bx(1) = 1.d-2
          Bx(2) = 1.d-2
!          Bx(4) = -1.d-1
!          Bx(1) = 1.d-3
!          Bx(2) = 1.d-3

!          By(1) = 1.d-2
          By(2) = -1.d-2
          By(3) = -1.d-2
!          By(2) = -1.d-3
!          By(3) = -1.d-3

          Bz(3) = 1.d-2
          Bz(4) = 1.d-2
!          Bz(3) = 1.d-3
!          Bz(4) = 1.d-3

        end if


        do it = 0, 100

        if (it.ge.1) then
!        if (it.eq.) then
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
           DO 200 I=1,NP
              VHXC(I) = 0.D0
              BXCX(I) = 0.D0
              BXCY(I) = 0.D0
              BXCZ(I) = 0.D0
200        CONTINUE
        ELSEIF (REST.EQ.1) THEN
           READ(2)VHXC,BXCX,BXCY,BXCZ
           REWIND 2
        ENDIF

1     CONTINUE
        ITER = ITER + 1
        IF (ITER.GT.1000000.and.rest.eq.1) then
          write(*,*) 'E not converging. Starting from initial condition'
          iter = 1
          rest = 0
          VHXC(I) = 0.D0
          BXCX(I) = 0.D0
          BXCY(I) = 0.D0
          BXCZ(I) = 0.D0
        elseif (iter.gt.600000.and.rest.eq.0) then
          write(*,*) 'System unconverged. Exiting...'
          call exit(-1)
        end if

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

        DO 2 I=1,NP
           VHXCO(I) = VHXC(I)
           BXCXO(I) = BXCX(I)
           BXCYO(I) = BXCY(I)
           BXCZO(I) = BXCZ(I)
2     CONTINUE

      IF (XC.EQ.1) THEN
         CALL XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      ELSEIF (XC.EQ.2) THEN
         CALL XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,sing)
       ELSEIF (XC.EQ.3) THEN
         CALL XCPOT_BALDA(U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ
     &                    ,EC)
       ELSEIF (XC.GE.4) THEN
         CALL XCPOT_STLS(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,EXC,1.D0,
     &                   ITER,XC)
       ENDIF

       IF (MP2.EQ.1) THEN
          CALL SLATER_INTEGRALS(U0,U1)
          CALL XCPOT_SLATER_MP2(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
       ENDIF

C**   blong calculates the longitudinal component of bxc and returns it
C**   as the new value for each bxc vector

        call blong(n, mx, my, mz, bxcx, bxcy, bxcz)


        CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)

C**   do the xc torque correction

        IF (CORR.EQ.1) THEN
!           WRITE(*,*)'total torque:'
!           WRITE(*,*)TT
!           WRITE(*,*)
           CALL BCORR(MX,MY,MZ,BXCX,BXCY,BXCZ,TT)
           CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
!           write(*,*) 'TC'
!           WRITE(*,*)'total torque:'
!           WRITE(*,*)TT
!           WRITE(*,*)
        ENDIF

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

!        open(201, file='crit_values', status='replace')
!        write(201,*) crit, eold
!        close(201)

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

        IF (XC.LE.3) CALL EX_CALC(U0,U1,EX)

        IF (MP2.EQ.1) THEN
            CALL SLATER_INTEGRALS(U0,U1)
            CALL EC_MP2_CALC(EC,VXC,BXCX,BXCY,BXCZ)
        ENDIF

        IF (XC.EQ.1.OR.XC.EQ.2) THEN
           EXC = EX
           IF (MP2.EQ.1) EXC = EXC + EC
        ELSEIF (XC.EQ.3) THEN
           EXC = EX + EC
        ENDIF

        IF (XC.GE.4) THEN
          EXC = 0.D0
          DL = 0.02D0
          DO 50 J=50,1,-1
C         DO 50 J=1,50
              LAMBDA = (J*1.D0-0.5D0)*DL
              CALL XCPOT_STLS(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,EXCL,
     &                      LAMBDA,0,XC)
              EXC = EXC + DL*EXCL
50        CONTINUE
        ENDIF

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
!95      write(*,*)real(TX(I)),' ',real(TY(I)),' ',real(TZ(I))
        !write(*,*)'-----------------------------------------'
        !write(*,*)real(TT(1)),' ',real(TT(2)),' ',real(TT(3))
        !write(*,*)

        tplot(1) = u0

        mbmats(it+1,1) = u0
        do i = 1, np
          mbmats(it+1,2+i) = mx(i)*bxcx(i)
     &                   + my(i)*bxcy(i)
     &                   + mz(i)*bxcz(i)
        end do
        mbmats(it+1,2) = eback

        do i = 1, 3
          etmats(it+1,1,i) = u0
          etmats(it+1,2,i) = etot

          emmats(it+1,1,i) = u0
          emmats(it+1,2,i) = etot

          exmats(it+1,1,i) = u0
          exmats(it+1,2,i) = etot
        end do

        emmats(it+1,1,4) = u0
        emmats(it+1,2,4) = etot

        exmats(it+1,1,4) = u0
        exmats(it+1,2,4) = etot

        do i=1, np
          etmats(it+1,2+i,1) = tx(i)
          etmats(it+1,2+i,2) = ty(i)
          etmats(it+1,2+i,3) = tz(i)

          emmats(it+1,2+i,1) = n(i)
          emmats(it+1,2+i,2) = mx(i)
          emmats(it+1,2+i,3) = my(i)
          emmats(it+1,2+i,4) = mz(i)

          exmats(it+1,2+i,1) = bxcx(i)
          exmats(it+1,2+i,2) = bxcy(i)
          exmats(it+1,2+i,3) = bxcz(i)
          exmats(it+1,2+i,4) = vhxc(i)
        end do

        ! write(*,etprint) etmats(:,:,1)

        write(*,*) '************************************'
        write(*,*) it

**************************************

        if (it.eq.100.and.qrev.eq.1) then
          do revit = 100, 0, -1

!            if (revit.le.99) then
!    !        if (it.eq.) then
              rest = 1
!            else
!              rest = 0
!            end if

            T = 0.5D0
            U0 = dble(revit)/10.d0
            U1 = U0/2.D0

            MIX = 0.1D0
            TOL = 1.D-10
            EOLD = 0.D0
            ITER = 0

            IF (REST.EQ.0) THEN
               DO I=1,NP
                  VHXC(I) = 0.D0
                  BXCX(I) = 0.D0
                  BXCY(I) = 0.D0
                  BXCZ(I) = 0.D0
               end do
            ELSEIF (REST.EQ.1) THEN
               READ(2)VHXC,BXCX,BXCY,BXCZ
               REWIND 2
            ENDIF

10          CONTINUE
            ITER = ITER + 1
            IF (ITER.GT.1000000.and.rest.eq.1) then
              write(*,*)
     &'E not converging. Starting from initial condition'
              iter = 1
              rest = 0
              VHXC(I) = 0.D0
              BXCX(I) = 0.D0
              BXCY(I) = 0.D0
              BXCZ(I) = 0.D0
            elseif (iter.gt.600000.and.rest.eq.0) then
              write(*,*) 'System unconverged. Exiting...'
              call exit(-1)
            end if

            DO I=1,NP
                VT(I) = V(I) + VHXC(I)
                BTX(I) = BX(I) + BXCX(I)
                BTY(I) = BY(I) + BXCY(I)
                BTZ(I) = BZ(I) + BXCZ(I)
            end do

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

            DO I=1,NP
               VHXCO(I) = VHXC(I)
               BXCXO(I) = BXCX(I)
               BXCYO(I) = BXCY(I)
               BXCZO(I) = BXCZ(I)
            end do

          IF (XC.EQ.1) THEN
             CALL XCPOT_SLATER(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
          ELSEIF (XC.EQ.2) THEN
             CALL XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,sing)
          ELSEIF (XC.EQ.3) THEN
              CALL XCPOT_BALDA(U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,N,MX,MY,MZ
     &                         ,EC)
          ELSEIF (XC.GE.4) THEN
              CALL XCPOT_STLS(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,EXC,1.D0,
     &                   ITER,XC)
          ENDIF

          IF (MP2.EQ.1) THEN
             CALL SLATER_INTEGRALS(U0,U1)
             CALL XCPOT_SLATER_MP2(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
          ENDIF

C**   blong calculates the longitudinal component of bxc and returns it
C**   as the new value for each bxc vector

          call blong(n, mx, my, mz, bxcx, bxcy, bxcz)


            CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)

C**   do the xc torque correction

            IF (CORR.EQ.1) THEN
    !           WRITE(*,*)'total torque:'
    !           WRITE(*,*)TT
    !           WRITE(*,*)
               CALL BCORR(MX,MY,MZ,BXCX,BXCY,BXCZ,TT)
               CALL TCALC (TT,TX,TY,TZ,MX,MY,MZ,BXCX,BXCY,BXCZ)
    !           write(*,*) 'TC'
    !           WRITE(*,*)'total torque:'
    !           WRITE(*,*)TT
    !           WRITE(*,*)
            ENDIF

            DO I=1,NP
               VHXC(I) = MIX*VHXC(I) + (1.D0-MIX)*VHXCO(I)
               BXCX(I) = MIX*BXCX(I) + (1.D0-MIX)*BXCXO(I)
               BXCY(I) = MIX*BXCY(I) + (1.D0-MIX)*BXCYO(I)
               BXCZ(I) = MIX*BXCZ(I) + (1.D0-MIX)*BXCZO(I)
            end do

            CRIT = 0.D0
            DO I=1,NP
    !25       CRIT = CRIT + N(I) + DABS(MX(I)) + DABS(MY(I)) + DABS(MZ(I))
              CRIT = CRIT + dsqrt(dble(MX(I)**2 + MY(I)**2 + MZ(I)**2))
            end do
    !        open(201, file='crit_values', status='replace')
    !        write(201,*) crit, eold
    !        close(201)

            IF (DABS((CRIT- EOLD)/CRIT).GT.TOL) THEN
                EOLD = CRIT
                GOTO 10
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
            DO I=1,NP
               EH = EH - 0.5D0*U0*N(I)**2
            end do

            DO I=1,NP-1
               EH = EH - U1*N(I)*N(I+1)
            end do

            EVXC = 0.D0
            DO I=1,NP
              EVXC = EVXC-N(I)*VXC(I)-MX(I)*BXCX(I)
     &                     -MY(I)*BXCY(I)-MZ(I)*BXCZ(I)
            end do

            IF (XC.LE.3) CALL EX_CALC(U0,U1,EX)

            IF (MP2.EQ.1) THEN
                CALL SLATER_INTEGRALS(U0,U1)
                CALL EC_MP2_CALC(EC,VXC,BXCX,BXCY,BXCZ)
            ENDIF

            IF (XC.EQ.1.OR.XC.EQ.2) THEN
               EXC = EX
               IF (MP2.EQ.1) EXC = EXC + EC
            ELSEIF (XC.EQ.3) THEN
               EXC = EX + EC
            ENDIF

            IF (XC.GE.4) THEN
               EXC = 0.D0
               DL = 0.02D0
               DO 500 J=50,1,-1
C             DO 50 J=1,50
                  LAMBDA = (J*1.D0-0.5D0)*DL
                  CALL XCPOT_STLS(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,EXCL,
     &                      LAMBDA,0,XC)
                  EXC = EXC + DL*EXCL
500          CONTINUE
            ENDIF

            eback = E(1) + E(2) + EH + EVXC + EXC

            !write(*,*)' EKS = ',E(1) + E(2)
            !write(*,*)'  EH = ',EH
            !write(*,*)'EVXC = ',EVXC
            !write(*,*)' EXC = ',EXC
            !write(*,*)
            !write(*,*)'ETOT = ',real(ETOT)

            !write(*,*)
            !write(*,*)'(xc torque)_x    (xc torque)_y    (xc torque)_z'
    !        DO 95 I=1,NP
    !95      write(*,*)real(TX(I)),' ',real(TY(I)),' ',real(TZ(I))
            !write(*,*)'-----------------------------------------'
            !write(*,*)real(TT(1)),' ',real(TT(2)),' ',real(TT(3))
            !write(*,*)

            tplot(1) = u0

            if (eback.lt.etot) then
              do i = 1, 3
                etmats(revit+1,2,i) = eback

                emmats(revit+1,2,i) = eback

                exmats(revit+1,2,i) = eback
              end do

              emmats(revit+1,2,4) = eback

              exmats(revit+1,2,4) = eback

              do i = 1, np
                mbmats(revit+1,2+i) = mx(i)*bxcx(i)
     &                   + my(i)*bxcy(i)
     &                   + mz(i)*bxcz(i)
              end do
              mbmats(revit+1,2) = eback

              do i = 1, np
                etmats(revit+1,2+i,1) = tx(i)
                etmats(revit+1,2+i,2) = ty(i)
                etmats(revit+1,2+i,3) = tz(i)

                emmats(revit+1,2+i,1) = n(i)
                emmats(revit+1,2+i,2) = mx(i)
                emmats(revit+1,2+i,3) = my(i)
                emmats(revit+1,2+i,4) = mz(i)

                exmats(revit+1,2+i,1) = bxcx(i)
                exmats(revit+1,2+i,2) = bxcy(i)
                exmats(revit+1,2+i,3) = bxcz(i)
                exmats(revit+1,2+i,4) = vhxc(i)
              end do
            end if

            write(*,*) '************************************'
            write(*,*) revit

          end do

        end if

**************************************

        end do

        if (wmflg.eq.1) then
          if (bmag.lt.10) then
            write(mwn,'(a,i1,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-n.txt'
            write(mwx,'(a,i1,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-mx.txt'
            write(mwy,'(a,i1,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-my.txt'
            write(mwz,'(a,i1,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-mz.txt'
          elseif (bmag.eq.10) then
            write(mwn,'(a,i2,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-n.txt'
            write(mwx,'(a,i2,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-mx.txt'
            write(mwy,'(a,i2,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-my.txt'
            write(mwz,'(a,i2,a)') '4pt-B', bmag,
     &                                '-Slat-AsymLong-mz.txt'
          elseif (bmag.eq.11) then
            write(mwn,'(a)') '4pt-B01-Slat-AsymLong-n.txt'
            write(mwx,'(a)') '4pt-B01-Slat-AsymLong-mx.txt'
            write(mwy,'(a)') '4pt-B01-Slat-AsymLong-my.txt'
            write(mwz,'(a)') '4pt-B01-Slat-AsymLong-mz.txt'
          end if

          open(100,file = mwn)
          open(101,file = mwx)
          open(102,file = mwy)
          open(103,file = mwz)
        end if


        if (wtflg.eq.1) then
          if (bmag.lt.10) then
            write(twx,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-tx.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-tx.txt'
            write(twy,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-ty.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-ty.txt'
            write(twz,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-tz.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-tz.txt'
          elseif (bmag.eq.10) then
            write(twx,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-tx.txt'
     &         '-Slat-AsymLong-tx.txt'
            write(twy,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-ty.txt'
     &         '-Slat-AsymLong-ty.txt'
            write(twz,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-tz.txt'
     &         '-Slat-AsymLong-tz.txt'
          elseif (bmag.eq.11) then
            write(twx,'(a)') '4pt-B01-Slat-AsymLong-tx.txt'
            write(twy,'(a)') '4pt-B01-Slat-AsymLong-ty.txt'
            write(twz,'(a)') '4pt-B01-Slat-AsymLong-tz.txt'
          end if

          open(201,file = twx)
          open(202,file = twy)
          open(203,file = twz)
        end if

        if (bxcflg.eq.1) then
          if (bmag.lt.10) then
            write(mwn,'(a)') '4pt-B1-Slat-AsymLong-vhxc.txt'
            write(twx,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-tx.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-bxcx.txt'
            write(twy,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-ty.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-bxcy.txt'
            write(twz,'(a,i1,a)')
!     &        '4pt-B', bmag, '-Slat-AsymLong-tz.txt'
     &         '4pt-B', bmag, '-Slat-AsymLong-bxcz.txt'
          elseif (bmag.eq.10) then
            write(twx,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-tx.txt'
     &         '-Slat-AsymLong-bxcx.txt'
            write(twy,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-ty.txt'
     &         '-Slat-AsymLong-bxcy.txt'
            write(twz,'(a,i2,a)') '4pt-B', bmag,
!     &        '-Slat-AsymLong-tz.txt'
     &         '-Slat-AsymLong-bxcz.txt'
          elseif (bmag.eq.11) then
            write(twx,'(a)') '4pt-B01-Slat-AsymLong-bxcx.txt'
            write(twy,'(a)') '4pt-B01-Slat-AsymLong-bxcy.txt'
            write(twz,'(a)') '4pt-B01-Slat-AsymLong-bxcz.txt'
          end if

          open(301,file = twx)
          open(302,file = twy)
          open(303,file = twz)
          open(304,file = mwn)

          open(400,file = '4pt-Slat-AsymLong-mdotb.txt')
        end if

        do i = 1, 4
          do j = 1,101
            do k = 1, 2+np
              if (i.lt.4) then
                if (i.eq.1) then
                  mdotb(k) = mbmats(j,k)
                endif
                etplot(k) = etmats(j,k,i)
              end if
              explot(k) = exmats(j,k,i)
              emplot(k) = emmats(j,k,i)
            end do
          if (i.lt.4) THEN
            if (i.eq.1) then
              write(399+i,etprint) mdotb
            end if
            write(200+i,etprint) etplot
          end if
          write(300+i,etprint) explot
          write(99+i,etprint) emplot
          end do
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

        if (bxcflg.eq.1) then
          close(301)
          close(302)
          close(303)
          close(304)
          close(400)
        end if

      end do

      END

C************************************************************************

      subroutine blong(n,mx,my,mz,bxcx,bxcy,bxcz)
      implicit none

      integer :: np
      INCLUDE 'dim.inc'

      real(8),intent(in) :: n(np), mx(np), my(np), mz(np)

      real(8),intent(inout) :: bxcx(np), bxcy(np), bxcz(np)

      real(8) :: blx(np), bly(np), blz(np)

      integer :: I,x,y,z

      real(8) :: mdotb(np), mm(np)


      do i = 1,np
        mm(i) = dsqrt(mx(i)**2 + my(i)**2 + mz(i)**2)
      end do

      do i = 1, np
        mdotb(i) = mx(i)*bxcx(i) + my(i)*BXCY(I) + mz(i)*bxcz(i)
      end do

      do i = 1, np
        blx(i) = mdotb(i)*mx(i)/mm(i)**2
        bly(i) = mdotb(i)*my(i)/mm(i)**2
        blz(i) = mdotb(i)*mz(i)/mm(i)**2
      end do

      do i=1, np
        bxcx(i) = blx(i)
        bxcy(i) = bly(i)
        bxcz(i) = blz(i)
      end do

      end subroutine

C************************************************************************

      SUBROUTINE XCPOT_STLS(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,EXC,
     &                      LAMBDA,ITER,XC)
      IMPLICIT NONE

      INTEGER I,J,K,L,R,RP,NP,ITER_STLS,ITER,XC
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX NN(2,2,NP),NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               NM(4,4,NP),DEN(NP),NMAT(16,16,NP,NP),
     &               S(2,2,2,2,NP,NP),D(2,2,2,2,NP,NP),SS,
     &               DVEC(16,NP,NP),FXCVEC(16,NP,NP),
     &               FXC(2,2,2,2,NP,NP),VXCMAT(2,2,NP),
     &               FHXC(2,2,2,2,NP,NP),FHXCO(2,2,2,2,NP,NP),
     &               GAMMA(2,2,NP,NP),PHI(2*NP,2,NP)

      DOUBLE PRECISION U0,U1,N(NP),VH(NP),VXC(NP),MIX,TOL,CRIT,LAMBDA,
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),HUB(NP,NP),
     &                 E(2*NP),EXC,FFO,FFN

      SAVE FHXC,FHXCO,FFO

      COMMON /EVALS/ E,PHI,GAMMA

C**   calculate the inverse density matrix:

      DO 1 I=1,NP
         NUU(I) = CDABS(PHI(1,1,I))**2 + CDABS(PHI(2,1,I))**2
         NUD(I) = PHI(1,1,I)*DCONJG(PHI(1,2,I))
     &          + PHI(2,1,I)*DCONJG(PHI(2,2,I))
         NDU(I) = DCONJG(NUD(I))
         NDD(I) = CDABS(PHI(1,2,I))**2 + CDABS(PHI(2,2,I))**2
         NN(1,1,I) = NUU(I)
         NN(1,2,I) = NUD(I)
         NN(2,1,I) = NDU(I)
         NN(2,2,I) = NDD(I)
1     CONTINUE

      DO 2 R=1,NP
         N(R) = DREAL(NUU(R) + NDD(R))
         DEN(R) = 2.D0*N(R)*( NUU(R)*NDD(R) - NUD(R)*NDU(R) )

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
2     CONTINUE

      DO 5 R=1,NP
      DO 5 RP=1,NP
         DO 6 I=1,4
         DO 6 J=1,4
         DO 6 K=1,4
         DO 6 L=1,4
            NMAT(4*(I-1)+K,4*(J-1)+L,R,RP)=NM(I,J,R)*DCONJG(NM(K,L,RP))
6        CONTINUE
         IF (R.EQ.RP) THEN
            HUB(R,RP) = U0
         ELSEIF (IABS(R-RP).EQ.1) THEN
            HUB(R,RP) = U1
         ELSE
            HUB(R,RP) = 0.D0
         ENDIF
5     CONTINUE

      IF (ITER.EQ.1) THEN
         DO 12 I=1,2
         DO 12 J=1,2
         DO 12 K=1,2
         DO 12 L=1,2
         DO 12 R=1,NP
         DO 12 RP=1,NP
            FHXC(I,J,K,L,R,RP) = (0.D0,0.D0)
            FHXCO(I,J,K,L,R,RP) = (0.D0,0.D0)
12       CONTINUE

         FFO = 0.D0
      ENDIF

C**   Start the STLS self-consistent iteration. The first step
C**   is the Slater exchange kernel

      ITER_STLS = 0
      MIX = 0.20D0
      IF (U0.GT.6.D0) MIX=0.1D0
      TOL = 1.D-8

900   CONTINUE
      ITER_STLS = ITER_STLS+1
      IF (ITER_STLS.EQ.100000) STOP

C**   Define the Hubbard interaction (including U0 and U1),
C**   calculate the static structure factor, and get the d-vector:

C      CALL SCALC_NONINT(S)
      CALL SCALC(FHXC,S,LAMBDA)

      DO 10 R=1,NP
      DO 10 RP=1,NP
         DO 11 I=1,2
         DO 11 J=1,2
         DO 11 K=1,2
         DO 11 L=1,2
            IF ((I.EQ.K).AND.(R.EQ.RP)) THEN
               SS = NN(L,J,R)
            ELSE
               SS = (0.D0,0.D0)
            ENDIF

            D(I,J,K,L,R,RP) = 4.D0*HUB(R,RP)*(S(I,J,K,L,R,RP) - SS)
11       CONTINUE
10    CONTINUE

      DO 15 R=1,NP
      DO 15 RP=1,NP
         DVEC(1,R,RP) = D(1,1,1,1,R,RP)
         DVEC(2,R,RP) = D(1,1,2,1,R,RP)
         DVEC(3,R,RP) = D(1,1,1,2,R,RP)
         DVEC(4,R,RP) = D(1,1,2,2,R,RP)
         DVEC(5,R,RP) = D(2,1,1,1,R,RP)
         DVEC(6,R,RP) = D(2,1,2,1,R,RP)
         DVEC(7,R,RP) = D(2,1,1,2,R,RP)
         DVEC(8,R,RP) = D(2,1,2,2,R,RP)
         DVEC(9,R,RP) = D(1,2,1,1,R,RP)
         DVEC(10,R,RP) = D(1,2,2,1,R,RP)
         DVEC(11,R,RP) = D(1,2,1,2,R,RP)
         DVEC(12,R,RP) = D(1,2,2,2,R,RP)
         DVEC(13,R,RP) = D(2,2,1,1,R,RP)
         DVEC(14,R,RP) = D(2,2,2,1,R,RP)
         DVEC(15,R,RP) = D(2,2,1,2,R,RP)
         DVEC(16,R,RP) = D(2,2,2,2,R,RP)
15    CONTINUE

C**   calculate fxc from the static structure factor

      DO 20 R=1,NP
      DO 20 RP=1,NP
      DO 20 I=1,16
         SS = (0.D0,0.D0)
         DO 19 J=1,16
            SS = SS + NMAT(I,J,R,RP)*DVEC(J,R,RP)
19       CONTINUE
         FXCVEC(I,R,RP) = SS
20    CONTINUE

      DO 21 R=1,NP
      DO 21 RP=1,NP
         FXC(1,1,1,1,R,RP) = FXCVEC(1,R,RP)
         FXC(1,1,2,1,R,RP) = FXCVEC(2,R,RP)
         FXC(1,1,1,2,R,RP) = FXCVEC(3,R,RP)
         FXC(1,1,2,2,R,RP) = FXCVEC(4,R,RP)
         FXC(2,1,1,1,R,RP) = FXCVEC(5,R,RP)
         FXC(2,1,2,1,R,RP) = FXCVEC(6,R,RP)
         FXC(2,1,1,2,R,RP) = FXCVEC(7,R,RP)
         FXC(2,1,2,2,R,RP) = FXCVEC(8,R,RP)
         FXC(1,2,1,1,R,RP) = FXCVEC(9,R,RP)
         FXC(1,2,2,1,R,RP) = FXCVEC(10,R,RP)
         FXC(1,2,1,2,R,RP) = FXCVEC(11,R,RP)
         FXC(1,2,2,2,R,RP) = FXCVEC(12,R,RP)
         FXC(2,2,1,1,R,RP) = FXCVEC(13,R,RP)
         FXC(2,2,2,1,R,RP) = FXCVEC(14,R,RP)
         FXC(2,2,1,2,R,RP) = FXCVEC(15,R,RP)
         FXC(2,2,2,2,R,RP) = FXCVEC(16,R,RP)
21    CONTINUE

      DO 22 I=1,2
      DO 22 J=1,2
      DO 22 K=1,2
      DO 22 L=1,2
      DO 22 R=1,NP
      DO 22 RP=1,NP
         FHXC(I,J,K,L,R,RP) = FXC(I,J,K,L,R,RP)
         IF ((I.EQ.J).AND.(K.EQ.L)) THEN
            FHXC(I,J,K,L,R,RP) = FHXC(I,J,K,L,R,RP) + HUB(R,RP)
         ENDIF

         IF (XC.EQ.4.OR.XC.EQ.5) THEN
         FHXC(I,J,K,L,R,RP) = MIX*FHXC(I,J,K,L,R,RP)
     &                     + (1.D0-MIX)*FHXCO(I,J,K,L,R,RP)
         ENDIF

         FHXCO(I,J,K,L,R,RP) = FHXC(I,J,K,L,R,RP)
22    CONTINUE

      FFN = 0.D0
      DO 25 I=1,16
      DO 25 R=1,NP
      DO 25 RP=1,NP
25       FFN = FFN + CDABS(FXCVEC(I,R,RP))
      CRIT = DABS((FFO-FFN)/FFN)
C      WRITE(*,*)CRIT
C      WRITE(30,*)CRIT
      FFO = FFN
      IF (XC.EQ.4.OR.XC.EQ.5) THEN
         IF (CRIT.GT.TOL) GOTO 900
         IF (ITER.GT.1) THEN
            WRITE(*,*)'STLS iterations: ',ITER_STLS
            WRITE(*,*)
         ENDIF
      ENDIF
C---------------------------------------------------------------------
C**   The STLS self-consistent iteration is completed. Now use the
C**   STLS xc kernel to construct the XC potential and magnetic field.
C---------------------------------------------------------------------
      IF (ITER.GE.1) THEN

         IF (XC.EQ.4.OR.XC.EQ.6) THEN
            DO 30 I=1,2
            DO 30 J=1,2
            DO 30 R=1,NP
               SS = (0.D0,0.D0)
               DO 31 K=1,2
               DO 31 L=1,2
               DO 31 RP=1,NP
                  SS = SS + FXC(I,J,K,L,R,RP)*NN(K,L,RP)
31             CONTINUE
               VXCMAT(I,J,R) = SS
30          CONTINUE
         ELSEIF (XC.EQ.5) THEN
            CALL STLS_OEP(U0,U1,FHXC,VXCMAT)
         ENDIF

         VH(1) = U0*N(1) + U1*N(2)
         DO 33 I=2,NP-1
            VH(I) = U0*N(I) + U1*N(I-1) + U1*N(I+1)
33       CONTINUE
         VH(NP) = U0*N(NP) + U1*N(NP-1)

         DO 40 R=1,NP
            VXC(R) = (VXCMAT(1,1,R) + VXCMAT(2,2,R))/2.D0
            VHXC(R) = VH(R) + VXC(R)
            BXCX(R) = (VXCMAT(1,2,R) + VXCMAT(2,1,R))/2.D0
            BXCY(R) = (0.D0,1.D0)*(VXCMAT(1,2,R) - VXCMAT(2,1,R))/2.D0
            BXCZ(R) = (VXCMAT(1,1,R) - VXCMAT(2,2,R))/2.D0
40       CONTINUE

C**   calculate EXC (for a given Lambda)
      ELSEIF (ITER.EQ.0) THEN

         DO 92 R=1,NP
            NM(1,1,R) = 2.D0*NUU(R)
            NM(1,2,R) = NUD(R)
            NM(1,3,R) = NDU(R)
            NM(1,4,R) = (0.D0,0.D0)

            NM(2,1,R) = NDU(R)
            NM(2,2,R) = NUU(R)+NDD(R)
            NM(2,3,R) = (0.D0,0.D0)
            NM(2,4,R) = NDU(R)

            NM(3,1,R) = NUD(R)
            NM(3,2,R) = (0.D0,0.D0)
            NM(3,3,R) = NUU(R)+NDD(R)
            NM(3,4,R) = NUD(R)

            NM(4,1,R) = (0.D0,0.D0)
            NM(4,2,R) = NUD(R)
            NM(4,3,R) = NDU(R)
            NM(4,4,R) = 2.D0*NDD(R)
92       CONTINUE

         DO 95 R=1,NP
         DO 95 RP=1,NP
            DO 96 I=1,4
            DO 96 J=1,4
            DO 96 K=1,4
            DO 96 L=1,4
               NMAT(4*(I-1)+K,4*(J-1)+L,R,RP)
     &              =NM(I,J,R)*DCONJG(NM(K,L,RP))
96          CONTINUE
95       CONTINUE

         DO 100 R=1,NP
         DO 100 RP=1,NP
         DO 100 I=1,16
            SS = (0.D0,0.D0)
            DO 101 J=1,16
               SS = SS + NMAT(I,J,R,RP)*FXCVEC(J,R,RP)
101         CONTINUE
         DVEC(I,R,RP) = SS
100      CONTINUE

         SS = (0.D0,0.D0)
         DO 200 R=1,NP
         DO 200 RP=1,NP
            SS = SS + DVEC(1,R,RP) + DVEC(4,R,RP)
     &              + DVEC(13,R,RP) + DVEC(16,R,RP)
200      CONTINUE

         EXC = 0.125D0*DREAL(SS)

      ENDIF

      RETURN
      END

C************************************************************************

      SUBROUTINE STLS_OEP(U0,U1,FHXC,VXCMAT)
      IMPLICIT NONE

      INTEGER NP,ND,I,J,MU,NU,ALPHA,BETA,R,RP,LWORK,INFO
      INCLUDE 'dim.inc'
      PARAMETER (ND=4*NP,LWORK=100)
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),DUM,
     &               LM(ND,ND),RVEC(ND),VXCMAT(2,2,NP),
     &               X(ND),X1(ND),UMAT(ND,ND),VTMAT(ND,ND),WORK(LWORK),
     &               FHXC(2,2,2,2,NP,NP)
      DOUBLE PRECISION U0,U1,E(2*NP),RWORK(100),SING(ND)

      COMMON /EVALS/ E,PHI,GAMMA

      DO 10 MU=1,2
      DO 10 NU=1,2
      DO 10 ALPHA=1,2
      DO 10 BETA=1,2
      DO 10 R=1,NP
      DO 10 RP=1,NP
         DUM = (0.D0,0.D0)
         DO 11 I=1,2
         DO 11 J=1,2*NP
            IF (I.NE.J) THEN
               DUM = DUM + (PHI(I,BETA,RP)*DCONJG(PHI(J,ALPHA,RP))
     &                    *DCONJG(PHI(I,MU,R))*PHI(J,NU,R)
     &                   + DCONJG(PHI(I,ALPHA,RP))*PHI(J,BETA,RP)
     &               *PHI(I,NU,R)*DCONJG(PHI(J,MU,R)))/(E(J)-E(I))
            ENDIF
11       CONTINUE
         LM(MU+(NU-1)*2+(R-1)*4,ALPHA+(BETA-1)*2+(RP-1)*4) = DUM
10    CONTINUE

      CALL BCALC(U0,U1,FHXC,RVEC)

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

C      WRITE(*,*)
C      WRITE(*,*)SING(ND),SING(ND-1),SING(ND-2),SING(ND-3)
C      WRITE(*,*)

      DO 72 I=1,ND
C         IF (DABS(SING(I)).GT.1.D-10) THEN
         IF (I.LT.ND-1) THEN
            SING(I) = 1.D0/SING(I)
         ELSE
            SING(I) = 0.D0
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

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE BCALC(U0,U1,FHXC,RVEC)
      IMPLICIT NONE

      INTEGER NP,ND,NEX,NEXH,M,R,RP,SI,TAU,NU,MU,I,J,A,IP,AP
      INCLUDE 'dim.inc'
      PARAMETER (ND=4*NP,NEX=4*NP-4,NEXH=2*NP-2)
      DOUBLE COMPLEX X(NEX,NEX),Y(NEX,NEX),DM1,DUM,IDR,
     &               DEXC(2,2,NP),B(2,2*NP),ZERO
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               FHXC(2,2,2,2,NP,NP),RVEC(ND)
      DOUBLE PRECISION E(2*NP),U0,U1
      PARAMETER (ZERO=(0.D0,0.D0))

      COMMON /EVALS/ E,PHI,GAMMA

      CALL CASIDA(FHXC,X,Y,1.D0)

      DO 10 I=1,2
      DO 10 SI=1,2
      DO 10 R=1,NP

         DM1 = ZERO

         DO 11 M=1,NEX
         DO 11 A=1,NEXH
            DUM = ZERO
            DO 12 RP=1,NP
               IDR = ZERO
               DO 13 IP=1,2
               DO 13 AP=1,NEXH
               DO 13 TAU=1,2
                  IDR = IDR + DCONJG(PHI(IP,TAU,RP))*PHI(AP+2,TAU,RP)
     &                        *DCONJG(X(AP+(IP-1)*NEXH,M))
     &                       + PHI(IP,TAU,RP)*DCONJG(PHI(AP+2,TAU,RP))
     &                        *DCONJG(Y(AP+(IP-1)*NEXH,M))
13             CONTINUE
               IF (RP.EQ.R) THEN
                  DUM = DUM + U0*IDR
               ELSEIF ( (RP.EQ.R-1).OR.(RP.EQ.R+1) ) THEN
                  DUM = DUM + U1*IDR
               ENDIF
12          CONTINUE
            DM1 = DM1 + DCONJG(PHI(A+2,SI,R))*X(A+(I-1)*NEXH,M)*DUM
11       CONTINUE

         DO 15 M=1,NEX
         DO 15 A=1,NEXH
            DUM = ZERO
            DO 16 RP=1,NP
               IDR = ZERO
               DO 17 IP=1,2
               DO 17 AP=1,NEXH
               DO 17 TAU=1,2
                  IDR = IDR + PHI(IP,TAU,RP)*DCONJG(PHI(AP+2,TAU,RP))
     &                       *X(AP+(IP-1)*NEXH,M)
     &                      + DCONJG(PHI(IP,TAU,RP))*PHI(AP+2,TAU,RP)
     &                       *Y(AP+(IP-1)*NEXH,M)
17             CONTINUE
               IF (RP.EQ.R) THEN
                  DUM = DUM + U0*IDR
               ELSEIF ( (RP.EQ.R-1).OR.(RP.EQ.R+1) ) THEN
                  DUM = DUM + U1*IDR
               ENDIF
16          CONTINUE
            DM1 = DM1
     &          + DCONJG(PHI(A+2,SI,R))*DCONJG(Y(A+(I-1)*NEXH,M))*DUM
15       CONTINUE

         DEXC(I,SI,R) = DM1 - U0*DCONJG(PHI(I,SI,R))
10    CONTINUE

      DO 30 I=1,2
      DO 30 J=1,2*NP
         IF (I.NE.J) THEN
            DUM = ZERO
            DO 31 R=1,NP
            DO 31 SI=1,2
               DUM = DUM + DEXC(I,SI,R)*PHI(J,SI,R)
31          CONTINUE
            B(I,J) = DUM/(E(J)-E(I))
         ELSE
            B(I,J) = ZERO
         ENDIF
30    CONTINUE

      DO 50 MU=1,2
      DO 50 NU=1,2
      DO 50 R=1,NP
         DUM = ZERO
         DO 51 I=1,2
         DO 51 J=1,2*NP
            DUM = DUM + DCONJG(B(I,J))*DCONJG(PHI(I,MU,R))*PHI(J,NU,R)
     &                + B(I,J)*DCONJG(PHI(J,MU,R))*PHI(I,NU,R)
51       CONTINUE
         RVEC(MU+(NU-1)*2+(R-1)*4) = DUM
50    CONTINUE

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE CASIDA(FHXC,X,Y,LAMBDA)
      IMPLICIT NONE

      INTEGER NP,NEX,NEXH,INFO,LWORK,I,J,K,R,RP,IP,A,AP,
     &        ALPHA,ALPHAP,SIGMA,SIGMAP
      PARAMETER (LWORK=100)
      INCLUDE 'dim.inc'
      PARAMETER (NEX=4*NP-4, NEXH=2*NP-2)
      DOUBLE COMPLEX PHI(2*NP,2,NP),GAMMA(2,2,NP,NP)
      DOUBLE COMPLEX OMEGA(NEX),K11(NEX,NEX),K12(NEX,NEX),K21(NEX,NEX),
     &               K22(NEX,NEX),X(NEX,NEX),Y(NEX,NEX),
     &               MAT(2*NEX,2*NEX),W(2*NEX),VL(2*NEX,2*NEX),
     &               VR(2*NEX,2*NEX),FHXC(2,2,2,2,NP,NP),WORK(LWORK),
     &               SWP,SW(2*NEX),D11,D12,D21,D22
      DOUBLE PRECISION E(2*NP),RWORK(4*NEX),NN,LAMBDA,SHIFT
      DOUBLE COMPLEX ZERO,ONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0))

      COMMON /EVALS/ E,PHI,GAMMA
      COMMON /STLS_SHIFT/ SHIFT

C      SHIFT = 0.D0

      DO 1 I=1,2
      DO 1 A=3,2*NP
         OMEGA(A-2 + (I-1)*NEXH) = (E(A)-E(I)+SHIFT)*ONE
1     CONTINUE

      DO 9 I=1,NEX
      DO 9 J=1,NEX
         K11(I,J) = ZERO
         K12(I,J) = ZERO
         K21(I,J) = ZERO
9        K22(I,J) = ZERO

C**   K11 = K_{ia,i'a'}
C**   K12 = K_{ia,a'i'}
C**   K21 = K_{ai,i'a'}
C**   K22 = K_{ai,a'i'}

      DO 10 I=1,2
      DO 10 A=1,NEXH
      DO 10 IP=1,2
      DO 10 AP=1,NEXH
         D11 = ZERO
         D12 = ZERO
         D21 = ZERO
         D22 = ZERO
         DO 11 ALPHA=1,2
         DO 11 ALPHAP=1,2
         DO 11 SIGMA=1,2
         DO 11 SIGMAP=1,2
         DO 11 R=1,NP
         DO 11 RP=1,NP
            D11 = D11 + DCONJG(PHI(I,ALPHA,R))*PHI(A+2,ALPHAP,R)
     &          * FHXC(ALPHA,ALPHAP,SIGMA,SIGMAP,R,RP)*LAMBDA
     &          * PHI(IP,SIGMA,RP)*DCONJG(PHI(AP+2,SIGMAP,RP))

            D12 = D12 + DCONJG(PHI(I,ALPHA,R))*PHI(A+2,ALPHAP,R)
     &          * FHXC(ALPHA,ALPHAP,SIGMA,SIGMAP,R,RP)*LAMBDA
     &          * PHI(AP+2,SIGMA,RP)*DCONJG(PHI(IP,SIGMAP,RP))

            D21 = D21 + DCONJG(PHI(A+2,ALPHA,R))*PHI(I,ALPHAP,R)
     &          * FHXC(ALPHA,ALPHAP,SIGMA,SIGMAP,R,RP)*LAMBDA
     &          * PHI(IP,SIGMA,RP)*DCONJG(PHI(AP+2,SIGMAP,RP))

            D22 = D22 + DCONJG(PHI(A+2,ALPHA,R))*PHI(I,ALPHAP,R)
     &          * FHXC(ALPHA,ALPHAP,SIGMA,SIGMAP,R,RP)*LAMBDA
     &          * PHI(AP+2,SIGMA,RP)*DCONJG(PHI(IP,SIGMAP,RP))
11       CONTINUE
         K11(NEXH*(I-1)+A,NEXH*(IP-1)+AP) = D11
         K12(NEXH*(I-1)+A,NEXH*(IP-1)+AP) = D12
         K21(NEXH*(I-1)+A,NEXH*(IP-1)+AP) = D21
         K22(NEXH*(I-1)+A,NEXH*(IP-1)+AP) = D22
10    CONTINUE

      DO 20 I=1,NEX
      DO 20 J=1,NEX
         MAT(I,J) = -K11(I,J)
         MAT(I,J+NEX) = -K12(I,J)
         MAT(I+NEX,J) = K21(I,J)
         MAT(I+NEX,J+NEX) = K22(I,J)
20    CONTINUE

      DO 21 I=1,NEX
         MAT(I,I) = MAT(I,I) - OMEGA(I)
         MAT(I+NEX,I+NEX) = MAT(I+NEX,I+NEX) + OMEGA(I)
21    CONTINUE

      CALL ZGEEV('N','V',2*NEX,MAT,2*NEX,W,VL,2*NEX,VR,2*NEX,
     &           WORK,LWORK,RWORK,INFO)

      DO 25 I=1,2*NEX
         IF (DABS(DIMAG(W(I))).GT.1.D-8) THEN
            WRITE(*,*)W
            WRITE(*,*)'xxxxxxxxxx   Casida collapses   xxxxxxxxxx'
            STOP
         ENDIF
25    CONTINUE

C**   now we need to sort the eigenvalues and eigenvectors

49    I=0
50    I=I+1
      IF (DREAL(W(I)).LE.DREAL(W(I+1))) THEN
         IF (I.LT.2*NEX-1) THEN
            GOTO 50
         ELSE
            GOTO 60
         ENDIF
      ELSE
         SWP = W(I+1)
         W(I+1) = W(I)
         W(I) = SWP

         DO 100 J=1,2*NEX
            SW(J) = VR(J,I+1)
            VR(J,I+1) = VR(J,I)
            VR(J,I) = SW(J)
100      CONTINUE

      ENDIF
      GOTO 49
60    CONTINUE

C**  Now we need to make sure the eigenvectors are properly normalized:
C**  -X_m X_m + Y_m Y_m = +/- 1 for positive and negative eigenvalues,
C**  respectively.

      DO 80 J=1,NEX
         NN = 0.D0
         DO 81 I=1,NEX
81          NN = NN + CDABS(VR(I,J))**2 - CDABS(VR(I+NEX,J))**2
         DO 82 I=1,2*NEX
82          VR(I,J) = VR(I,J)/DSQRT(DABS(NN))
80    CONTINUE

      DO 85 J=NEX+1,2*NEX
         NN = 0.D0
         DO 86 I=1,NEX
86          NN = NN - CDABS(VR(I,J))**2 + CDABS(VR(I+NEX,J))**2
         DO 87 I=1,2*NEX
87          VR(I,J) = VR(I,J)/DSQRT(DABS(NN))
85    CONTINUE

C**   We are interested in the positive eigenvalues only, and the
C**   associated eigenvectors, which we call X and Y

      DO 200 K=1,NEX
      DO 200 J=1,NEX
         X(J,K) = VR(J,K+NEX)
         Y(J,K) = VR(J+NEX,K+NEX)
200   CONTINUE

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE SCALC(FHXC,S,LAMBDA)
      IMPLICIT NONE

      INTEGER NP,MM,R,RP,SI,SIP,GA,GAP,I,A,IP,AP,NEX,NEXH
      INCLUDE 'dim.inc'
      PARAMETER (NEX=4*NP-4, NEXH=2*NP-2)
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               S(2,2,2,2,NP,NP),X(NEX,NEX),Y(NEX,NEX),SS,
     &               FHXC(2,2,2,2,NP,NP)
      DOUBLE PRECISION E(2*NP),LAMBDA

      COMMON /EVALS/ E,PHI,GAMMA

      CALL CASIDA(FHXC,X,Y,LAMBDA)

C      DO 5 I=1,NEX
C      DO 5 A=1,NEX
C         X(I,A) = (0.D0,0.D0)
C         IF (I.EQ.A) THEN
C            Y(I,A) = (1.D0,0.D0)
C         ELSE
C            Y(I,A) = (0.D0,0.D0)
C         ENDIF
C5     CONTINUE

      DO 10 SI=1,2
      DO 10 SIP=1,2
      DO 10 GA=1,2
      DO 10 GAP=1,2
      DO 10 R=1,NP
      DO 10 RP=1,NP
          SS = (0.D0,0.D0)
          DO 11 MM=1,NEX
             DO 12 I=1,2
             DO 12 A=1,NEXH
             DO 12 IP=1,2
             DO 12 AP=1,NEXH
                SS = SS + (PHI(I,SI,R)*DCONJG(PHI(A+2,SIP,R))
     &                    *X(A+(I-1)*NEXH,MM)
     &                    +DCONJG(PHI(I,SIP,R))*PHI(A+2,SI,R)
     &                    *Y(A+(I-1)*NEXH,MM))
     &                  * (DCONJG(PHI(IP,GA,RP))*PHI(AP+2,GAP,RP)
     &                    *DCONJG(X(AP+(IP-1)*NEXH,MM))
     &                    +PHI(IP,GAP,RP)*DCONJG(PHI(AP+2,GA,RP))
     &                    *DCONJG(Y(AP+(IP-1)*NEXH,MM)))
12           CONTINUE
11        CONTINUE
          S(SI,SIP,GA,GAP,R,RP) = SS
10    CONTINUE

      END
C************************************************************************
C************************************************************************
C************************************************************************
      SUBROUTINE SCALC_NONINT(S)
      IMPLICIT NONE

      INTEGER I,J,K,NP,SIG,SIGP,TAU,TAUP
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),
     &               S(2,2,2,2,NP,NP),SS
      DOUBLE PRECISION E(2*NP)

      COMMON /EVALS/ E,PHI,GAMMA

C**   This subroutine calculates the noninteracting structure factor.
C**   We need this as the first step of the STLS iteration.

      DO 20 I=1,NP
      DO 20 J=1,NP
      DO 20 SIG=1,2
      DO 20 SIGP=1,2
      DO 20 TAU=1,2
      DO 20 TAUP=1,2
         SS = (0.D0,0.D0)
         DO 30 K=3,2*NP
            SS = SS + PHI(K,SIG,I)*DCONJG(PHI(K,TAU,J))
30       CONTINUE
         SS = SS*(DCONJG(PHI(1,SIGP,I))*PHI(1,TAUP,J)
     &          + DCONJG(PHI(2,SIGP,I))*PHI(2,TAUP,J))
         S(SIG,SIGP,TAU,TAUP,I,J)=SS
20    CONTINUE

      END

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

      COP=2
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
10      CONTINUE

        DO 11 I=1,NP
           PX(I) = MY(1)*AX(I) + MY(I)
           PY(I) = MY(1)*AY(I) - MX(I)
           PZ(I) = MY(1)*AZ(I)

           QX(I) = MZ(1)*AX(I) + MZ(I)
           QY(I) = MZ(1)*AY(I)
           QZ(I) = MZ(1)*AZ(I) - MX(I)
11      CONTINUE

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

20      CONTINUE

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
     &                + (MY(1)*AX(I+1)+MY(I+1))*PZ(J)
     &                + (MZ(1)*AX(I+1)+MZ(I+1))*QZ(J)

           MAT(II+2,JJ+2) = MX(NP)**2*AX(I+1)*AX(J+1)
     &                + (MY(1)*AX(I+1)+MY(I+1))*PX(J+1)
     &                + (MZ(1)*AX(I+1)+MZ(I+1))*QX(J+1)

30      CONTINUE

        DO 35 I=1,NM
           MAT(I,I) = MAT(I,I) + MX(NP)**2
35      CONTINUE

        CALL DSYSV('L',NM,1,MAT,NM,IPIV,RVEC,NM,WORK,LWORK,INFO )

        DO 40 I=1,NP-1
           II = 3*(I-1)+1
           BBY(I) = RVEC(II)
           BBZ(I) = RVEC(II+1)
           BBX(I+1) = RVEC(II+2)
40      CONTINUE

        BBX(1) = 0.D0
        DO 45 I=1,NP-1
           BBX(1) = BBX(1) + AY(I)*BBY(I) + AZ(I)*BBZ(I)
     &                     + AX(I+1)*BBX(I+1)
45      CONTINUE

        BBY(NP) = MY(NP)*BBX(NP)/MX(NP)
        BBZ(NP) = MZ(NP)*BBX(NP)/MX(NP)
        DO 50 I=1,NP-1
           BBY(NP) = BBY(NP) + MY(I)*BBX(I)/MX(NP)
     &                       - MX(I)*BBY(I)/MX(NP)
           BBZ(NP) = BBZ(NP) + MZ(I)*BBX(I)/MX(NP)
     &                       - MX(I)*BBZ(I)/MX(NP)
50      CONTINUE

        DO 100 I=1,NP
           BX(I) = - BBX(I)
           BY(I) = - BBY(I)
           BZ(I) = - BBZ(I)
100     CONTINUE

      ELSE
C**--------------------------------------------------------------------***
C**   coplanar, so we have an explicit solution
C**--------------------------------------------------------------------***
        DO 200 I=1,NP
           BX(I) = BX(I) - (TT(2)*MZ(I) - TT(3)*MY(I))/M2
           BY(I) = BY(I) - (TT(3)*MX(I) - TT(1)*MZ(I))/M2
           BZ(I) = BZ(I) - (TT(1)*MY(I) - TT(2)*MX(I))/M2
200     CONTINUE

      ENDIF

      RETURN
      END

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

      SUBROUTINE XCPOT_OEP(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ,sing)
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

!      open(101, file='svd_values', status='replace')
!      write(101,*) sing
!      close(101)

C**   The singular values are ordered from largest to smallest.
C**   By default, we drop the last singular value SING(ND).
C**   But for strong correlations (large U), it may happen that
C**   other singular values need to be dropped as well. This needs
C**   to be checked from case to case, and, if necessary, the code
C**   below needs to be modified.

      DO 72 I=1,ND
!         IF (I.LT.ND) THEN
         IF (dabs(SING(I)).gt.1d-14) then
!         IF (I.LT.ND-1) THEN
           SING(I) = 1.D0/SING(I)
         ELSE
           SING(I) = 0.D0
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

      SUBROUTINE XCPOT_SLATER_MP2(U0,U1,VXC,VHXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,I,J,K,R,MU,NU
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX NUU(NP),NUD(NP),NDU(NP),NDD(NP),
     &               VUU(NP),VUD(NP),VDU(NP),VDD(NP),DEN(NP),
     &               BUU(NP),BUD(NP),BDU(NP),BDD(NP),MAT(4,4),IONE,
     &               ZERO,GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),DR,
     &               SL(2,2,2*NP,2*NP),RHS(2,2,NP),
     &               DEDPS(2*NP,2,NP)
      PARAMETER (ZERO=(0.D0,0.D0),IONE=(0.D0,1.D0))
      DOUBLE PRECISION N(NP),MX(NP),MY(NP),MZ(NP),VH(NP),VXC(NP),
     &                 VHXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP),E(2*NP),
     &                 U0,U1

      COMMON /EVALS/ E,PHI,GAMMA
      COMMON /SLINT/ SL,DEDPS

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
         MY(I) = DREAL(IONE*(NUD(I) - NDU(I)))
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

C**   Define the MP2 right-hand side
      DO 50 R=1,NP
      DO 50 NU=1,2
      DO 50 MU=1,2
         DR=ZERO
         DO 60 I=1,2*NP
            DR = DR + DEDPS(I,NU,R)*DCONJG(PHI(I,MU,R))
     &              + DCONJG(DEDPS(I,MU,R))*PHI(I,NU,R)
60       CONTINUE
         RHS(NU,MU,R)=DR
50    CONTINUE

C**   Now construct the Slater potential
      DO 65 R=1,NP
         BUU(R) = BUU(R) + RHS(1,1,R)
         BDU(R) = BDU(R) + RHS(2,1,R)
         BUD(R) = BUD(R) + RHS(1,2,R)
         BDD(R) = BDD(R) + RHS(2,2,R)
65    CONTINUE

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

C**********************************************************************

      SUBROUTINE XCPOT_BALDA(U0,T,VXC,VHXC,BXCX,BXCY,BXCZ,
     &                       NN,MX,MY,MZ,EC)
      IMPLICIT NONE

      INTEGER NP,I
      INCLUDE 'dim.inc'
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

         VXC(I) = VXC(I) + S*VC
         VHXC(I) = VXC(I) + U0*NN(I)

         IF (M(I).GT.1.D-15) THEN
         BCX = T*BC*MX(I)/M(I) + U0*MX(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MX(I)/M(I)

         BCY = T*BC*MY(I)/M(I) + U0*MY(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MY(I)/M(I)

         BCZ = T*BC*MZ(I)/M(I) + U0*MZ(I)/2.D0
     &       - 2.D0*T*DSIN(PI*N(I)/2.D0)*DSIN(PI*M(I)/2.D0)*MZ(I)/M(I)

         BXCX(I) = BXCX(I) + BCX
         BXCY(I) = BXCY(I) + BCY
         BXCZ(I) = BXCZ(I) + BCZ
         ENDIF

1     CONTINUE

      EC = 0.D0
      DO 55 I=1,NP
         EC = EC + ECD(I)
55    CONTINUE

      RETURN
      END

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

C************************************************************************

      SUBROUTINE EC_MP2_CALC(EC,VXC,BXCX,BXCY,BXCZ)
      IMPLICIT NONE

      INTEGER NP,I,J,K,L,TAU,SIGMA,R,RP
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),SL(2,2,2*NP,2*NP),
     &               SD,ECC,DEDPS(2*NP,2,NP),V(2,2,NP),VIL(2,2*NP),SLM
      DOUBLE PRECISION E(2*NP),EC,VXC(NP),BXCX(NP),BXCY(NP),BXCZ(NP)
      DOUBLE COMPLEX ZERO,ONE,IONE
      PARAMETER (ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),IONE=(0.D0,1.D0))

      COMMON /EVALS/ E,PHI,GAMMA
      COMMON /SLINT/ SL,DEDPS

C**   Define the orbitals phi(m,sigma,x):
C**   m = 1...2*NP is the orbital index
C**   sigma = 1,2 (up, down) is the spin index
C**   x = 1,...,NP (lattice points) is the spatial coordinate

C**   First, calculate the MP2 contribution

      ECC = ZERO
      DO 40 I=1,2
      DO 40 J=1,2
      DO 40 K=3,2*NP
      DO 40 L=3,2*NP
         ECC = ECC + SL(I,J,K,L)*DCONJG(SL(I,J,K,L)-SL(J,I,K,L))
     &               /(E(I)+E(J)-E(K)-E(L))
40    CONTINUE

      EC = 0.5D0*DREAL(ECC)
      WRITE(*,*)'        MP2:',EC

C**   Now calculate the DeltaHF contribution

      DO 50 R=1,NP
         V(1,1,R) = (VXC(R) + BXCZ(R))*ONE
         V(1,2,R) = BXCX(R)*ONE - BXCY(R)*IONE
         V(2,1,R) = BXCX(R)*ONE + BXCY(R)*IONE
         V(2,2,R) = (VXC(R) - BXCZ(R))*ONE
50    CONTINUE

      DO 60 I=1,2
      DO 60 L=1,2*NP
         ECC = ZERO
         DO 61 R=1,NP
            ECC = ECC + DCONJG(PHI(I,1,R))*V(1,1,R)*PHI(L,1,R)
     &                + DCONJG(PHI(I,1,R))*V(1,2,R)*PHI(L,2,R)
     &                + DCONJG(PHI(I,2,R))*V(2,1,R)*PHI(L,1,R)
     &                + DCONJG(PHI(I,2,R))*V(2,2,R)*PHI(L,2,R)
61       CONTINUE
         VIL(I,L) = ECC
60    CONTINUE

      ECC = ZERO
      DO 70 I=1,2
      DO 70 L=3,2*NP
         ECC = ECC + CDABS(VIL(I,L) + SL(I,1,1,L) + SL(I,2,2,L))**2
     &               /(E(I)-E(L))
70    CONTINUE

      EC = EC + DREAL(ECC)
      WRITE(*,*)'MP2+DeltaHF:',EC

      RETURN
      END

C************************************************************************

      SUBROUTINE SLATER_INTEGRALS(U0,U1)
      IMPLICIT NONE

      INTEGER NP,H,I,J,K,L,TAU,SIGMA,R,RP
      INCLUDE 'dim.inc'
      DOUBLE COMPLEX GAMMA(2,2,NP,NP),PHI(2*NP,2,NP),SL(2,2,2*NP,2*NP),
     &               SD,HI(2,2*NP,NP),DEDPS(2*NP,2,NP)
      DOUBLE PRECISION E(2*NP),U0,U1

      COMMON /EVALS/ E,PHI,GAMMA
      COMMON /SLINT/ SL,DEDPS

C**   Define the orbitals phi(m,sigma,x):
C**   m = 1...2*NP is the orbital index
C**   sigma = 1,2 (up, down) is the spin index
C**   x = 1,...,NP (lattice points) is the spatial coordinate

!     m is orbital index (running from ground state to highest calculated level)
!     with x representing the spatial coordinate and sigma representing
!     spin_up (1) or spin_down (2).
!     ex: the ground-state, spin_up electron on the second lattice site:
!     phi(1,1,2)

C**   define Slater integrals

      DO 10 I=1,2
      DO 10 J=1,2
      DO 10 K=1,2*NP
      DO 10 L=1,2*NP
         SD = (0.D0,0.D0)
         DO 20 SIGMA=1,2
         DO 20 TAU=1,2
            DO 30 R=1,NP
            DO 30 RP=1,NP
               IF (R.EQ.RP) THEN
                  SD = SD + U0*DCONJG(PHI(I,SIGMA,R))*PHI(K,SIGMA,R)
     &                        *DCONJG(PHI(J,TAU,R))*PHI(L,TAU,R)
               ELSEIF (R.EQ.RP+1.OR.R.EQ.RP-1) THEN
                  SD = SD + U1*DCONJG(PHI(I,SIGMA,R))*PHI(K,SIGMA,R)
     &                        *DCONJG(PHI(J,TAU,RP))*PHI(L,TAU,RP)
               ENDIF
30          CONTINUE
20       CONTINUE
         SL(I,J,K,L) = SD
10    CONTINUE

C**   define H-integrals

      DO 50 I=1,2
      DO 50 K=1,2*NP
      DO 50 R=1,NP
         SD = (0.D0,0.D0)
         IF (K.GE.3) THEN
         DO 51 SIGMA=1,2
         DO 52 RP=1,NP
            IF (R.EQ.RP) THEN
               SD = SD + U0*DCONJG(PHI(I,SIGMA,R))*PHI(K,SIGMA,R)
            ELSEIF (R.EQ.RP+1.OR.R.EQ.RP-1) THEN
               SD = SD + U1*DCONJG(PHI(I,SIGMA,RP))*PHI(K,SIGMA,RP)
            ENDIF
52       CONTINUE
51       CONTINUE
         ENDIF
         HI(I,K,R)=SD
50    CONTINUE

C**   calculate the functional derivatives, dE/dPsi*

!     H within the occupied orbitals

      DO 100 H=1,2
      DO 100 SIGMA=1,2
      DO 100 R=1,NP
         SD = (0.D0,0.D0)
         DO 101 J=1,2
         DO 101 K=3,2*NP
         DO 101 L=3,2*NP
            SD = SD+(PHI(K,SIGMA,R)*HI(J,L,R)-PHI(L,SIGMA,R)*HI(J,K,R))
     &               *DCONJG(SL(H,J,K,L)-SL(J,H,K,L))
     &               /(E(H)+E(J)-E(K)-E(L))
101      CONTINUE
         DEDPS(H,SIGMA,R) = 0.5D0*SD
100   CONTINUE


!     H in the unoccupied orbitals.

      DO 150 H=3,2*NP
      DO 150 SIGMA=1,2
      DO 150 R=1,NP
         SD = (0.D0,0.D0)
         DO 151 I=1,2
         DO 151 J=1,2
         DO 151 K=3,2*NP
            SD = SD  +(PHI(I,SIGMA,R)*DCONJG(HI(J,K,R))
     &                -PHI(J,SIGMA,R)*DCONJG(HI(I,K,R)))
     &               *(SL(I,J,H,K)-SL(I,J,K,H))
     &               /(E(I)+E(J)-E(H)-E(K))
151      CONTINUE
         DEDPS(H,SIGMA,R) = 0.5D0*SD
150   CONTINUE

      RETURN
      END
