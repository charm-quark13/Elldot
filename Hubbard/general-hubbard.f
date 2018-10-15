      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE generalCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it

      real(8) :: ftol,fret,x,v(dim*2),dens(dim*2)
      real(8) :: Bx(sites), By(sites), Bz(sites)
      real(8) :: txc1,txc2
      real(8) :: vstart(dim*2),vhxc(2),bxc(intd)
      complex(8) :: htest(intd,intd)
      real(8) :: tau1(3),tau2(3)

      write(matrix,'(a, i3, a)') '(', dim, 'f14.8)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

      call Pauli(sig)

      MIX = 0.10D0
      TOL = 1.D-12

      write(*,*)'Is this a restart of a previous calculation?'
      write(*,*)'0 == No ... 1 == Yes'
      read(*,*) rest

      write(*,*)'Should a torque correction be made?'
      write(*,*)'0 == No ... 1 == Yes'
      read(*,*) corr

      write(*,*)'XC=1: Slater'
      write(*,*)'XC=2: OEP'
      write(*,*)'XC=3: BALDA (U1=0)'
      read(*,*)XC

      v = 0.d0
      Bx = 0.d0
      by = 0.d0
      bz = 0.d0

      Bx(1) = -2.5d0
      Bx(2) = 1.0d0
!      Bx = 0.d0
!      By(1) = .0d0
      By(2) = .75d0
      Bz(1) = .0d0
      Bz(2) = .0d0

***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
      do i=1,sites
        v(1) = 2.5d0
        if (i.ne.1) then
          v(i) = -2.5d0
        end if
      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

      vstart = v

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      do it=10,10

!        do i=1,sites
!          v(1) = 1.d0
!          if (i.ne.1) then
!            v(i) = -1.d0
!          end if
!        end do

!        do i=1,sites
!          v(sites+i) = Bx(i)
!          v(sites*2+i) = By(i)
!          v(sites*3+i) = Bz(i)
!        end do

        U0 = it/10.d0
        U1 = U0/dsqrt(2.d0)
        ntarget = 0.d0

        call interHam(v,U0,U1,htest)
        call intdens(ntarget,htest)

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
!        number = 0
        call frprmn(v,dim*2,ftol,iter,fret)

        write(*,*) '*********************************'
        write(*,*) iter, ftol, fret

***   Final check to ensure ntarget is equivalent to the density found by our
***   conjugate gradient method.
        call hbuild(v,hmat)
        call densvec(dens,hmat)

        write(*,vector) dens
        write(*,*) '^^^^^ dens ^^^^^'
        write(*,vector) ntarget
        write(*,*) '^^^^^ n_target ^^^^^'

        vhxc = 0.d0
        do i=1,2
          vhxc(i) = v(i) - vstart(i)
        end do

        bxc = 0.d0
        do i=1,6
          k=i+2
          bxc(i) = v(k) - vstart(k)
        end do

        tau1(1) = dens(5)*bxc(5) - dens(7)*bxc(3)
        tau1(2) = -(dens(3)*bxc(5)-dens(7)*bxc(1))
        tau1(3) = dens(3)*bxc(3) - dens(5)*bxc(1)

        tau2(1) = dens(6)*bxc(6) - dens(8)*bxc(4)
        tau2(2) = -(dens(4)*bxc(6)-dens(8)*bxc(2))
        tau2(3) = dens(4)*bxc(4) - dens(6)*bxc(2)

        txc1 = 0.d0
        txc2 = 0.d0
        do i=1,3
          txc1 = txc1 + tau1(i)*tau1(i)
          txc2 = txc2 + tau2(i)*tau2(i)
        end do

        txc1 = dsqrt(txc1)
        txc2 = dsqrt(txc2)

        open(100,file='xc-torque-04.txt',access='append')

        if (it.eq.1) then
          write(100,*) 'U','xc torque 1', 'xc torque 2'
        end if
        write(100,*) u0,txc1,-txc2
        write(*,*) '*** end of loop ***'
      end do

      end
