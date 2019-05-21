      program interHam
      implicit none

      integer, parameter :: spin=2, sites=4
      integer, parameter :: dim = spin*sites, intd=2*sites**2-sites
      integer,parameter :: lwork=600
      real(8), parameter :: t=0.5d0
      complex(8), parameter :: zero=(0.d0,0.d0), ione=(0.d0,1.d0),
     &                         one=(1.d0,0.d0)

      integer :: Ntrp,Nsng,i,j,k,l,ii,jj,info,iter
      real(8) :: rwork(100),cn(intd), eplot(intd+1)
      real(8) :: v(dim*2),U0,U1,bx(sites),by(sites),bz(sites)
      complex(8) :: ham(intd,intd)
      complex(8) :: work(lwork)
      !      complex(8) :: dum,vl(6,6),vr(6,6),cn(6)
      complex(8) :: Bp(sites),Bm(sites)

      character(25) :: ewrite

      write(ewrite, '(a, i2, a)') '(', intd + 1, 'e16.6)'

!      open(101, file='SymCP-Exact-E.txt')
      open(101, file='AsymCP-Exact-E.txt')


      do iter=0,100

        u0 = dble(iter)/10.d0
        u1 = u0/2.d0

        bx = 0.d0
        by = 0.d0
        bz = 0.d0

        v(1) = 1.d0
        do i = 2, 3
          v(i) = -1.d0
        end do
        v(sites) = 1.d0

***************** Symmetric Case *******************

!        Bx(1) = (1.d0/10.d0)
!        Bx(2) = (1.d0/10.d0)

!        Bz(3) = (1.d0/10.d0)
!        Bz(4) = (1.d0/10.d0)

****************************************************

***************** Asymmetric Case ******************
        Bx(1) = 2d-1
        Bx(2) = 1.d0/10.d0

        Bz(3) = 3d-1
        Bz(4) = -2d-1

****************************************************

        do i=1,sites
          v(sites+i) = Bx(i)
          v(sites*2+i) = By(i)
          v(sites*3+i) = Bz(i)
        end do

        ham = ZERO

        DO I=1,sites
          BP(I) = (v(sites+i)*ONE + v(sites*2+I)*IONE)/DSQRT(2.D0)
          BM(I) = (v(sites+i)*ONE - v(sites*2+I)*IONE)/DSQRT(2.D0)
        end do

        NSng = sites*(sites+1)/2
        NTrp = 3*sites*(sites-1)/2

C**------------------------------------------------------------
C**   Singlet block
C**------------------------------------------------------------

C**   First the diagonal elements

        II = -sites
        DO I=1,sites
           II = II + sites+2-I
           ham(II,II) = 2.D0*V(I) + U0
        end do

        II = -sites+1
        DO I=1,sites-1
           II = II + sites+2-I
           ham(II,II) = V(I) + V(I+1) + U1
        end do

        IF (sites.GT.2) THEN
          JJ=-sites+2
          DO J=1,sites-2
            JJ = JJ + sites+2-J
            DO I=1,sites-J-1
              ham(JJ+I-1,JJ+I-1) = V(J) + V(J+I+1)
            end do
          end do
        ENDIF

C**   Now the off-diagonal elements

        II = -sites
        DO I=1,sites-1
          II = II + sites+2-I
          ham(II,II+1) = -DSQRT(2.D0)*T
          ham(II+1,II) = -DSQRT(2.D0)*T
        end do

        II = -sites+1
        DO I=1,sites-1
          II = II + sites+2-I
          ham(II,II+sites-I) = -DSQRT(2.D0)*T
          ham(II+sites-I,II) = -DSQRT(2.D0)*T
        end do

        IF (sites.GT.2) THEN
          JJ=-sites+2
          DO J=1,sites-2
            JJ = JJ + sites+2-J
            DO I=1,sites-J-1
              ham(JJ+I-2,JJ+I-1) = -T
              ham(JJ+I-1,JJ+I-2) = -T

              ham(JJ+I-1,JJ+I-1+(sites-J)) = -T
              ham(JJ+I-1+(sites-J),JJ+I-1) = -T
            end do
          end do
        ENDIF

C**------------------------------------------------------------
C**   Singlet-Triplet and Triplet-Singlet blocks
C**------------------------------------------------------------

        II = -sites+1
        JJ = NSng+1
        DO J=1,sites-1
           II = II + sites-J+2
           DO I=1,sites-J

              ham(II+I-1,JJ) = -BP(J) + BP(J+I)
              ham(II+I-1,JJ+1) = v(sites*3+j) - v(sites*3+J+I)
              ham(II+I-1,JJ+2) = BM(J) - BM(J+I)

              ham(JJ,II+I-1) = -BM(J) + BM(J+I)
              ham(JJ+1,II+I-1) = v(sites*3+j) - v(sites*3+J+I)
              ham(JJ+2,II+I-1) = BP(J) - BP(J+I)

              JJ=JJ+3
           end do
        end do
C**------------------------------------------------------------
C**   Triplet block
C**------------------------------------------------------------

        II = Nsng-1
        DO I=1,sites-1
           DO J=I+1,sites
              II = II + 3
              ham(II-1,II-1) = V(I) + V(J)
     &                          + v(sites*3+I) + v(sites*3+J)
              ham(II,II) =  V(I) + V(J)
              ham(II+1,II+1) = V(I) + V(J)
     &                          - v(sites*3+I) - v(sites*3+J)

              IF ((J-I).EQ.1) THEN
                 ham(II-1,II-1) = ham(II-1,II-1) + U1
                 ham(II,II) = ham(II,II) + U1
                 ham(II+1,II+1) = ham(II+1,II+1) + U1
              ENDIF

              ham(II-1,II) = BM(I) + BM(J)
              ham(II,II-1) = BP(I) + BP(J)
              ham(II,II+1) = BM(I) + BM(J)
              ham(II+1,II) = BP(I) + BP(J)
           end do
        end do

        II = NSng-2
        DO I=1,sites-1
           DO J=I+1,sites
              II = II+3

              JJ = NSng-2
              DO K=1,sites-1
                 DO L=K+1,sites
                    JJ = JJ+3
                    IF ((I.EQ.K).AND.(IABS(J-L).EQ.1).OR.
     &                (J.EQ.L).AND.(IABS(I-K).EQ.1)) THEN
                       ham(II,JJ) = -T
                       ham(II+1,JJ+1) = -T
                       ham(II+2,JJ+2) = -T
                    ENDIF
                 end do
              end do

           end do
        end do

        call ZHEEV('v','l', intd, ham, intd, cn, work, lwork, rwork,
     &             info)

        eplot = 0
        eplot(1) = u0
        do i = 1, intd
          eplot(i+1) = cn(i)
        end do

        write(101,ewrite) eplot

      end do

      close(101)

      end
