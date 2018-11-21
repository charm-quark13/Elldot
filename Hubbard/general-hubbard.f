      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE generalCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it,steps

      real(8) :: ftol,fret,fretlo,x
      real(8) :: Bx(sites), By(sites), Bz(sites),mags(sites+1)
      real(8) :: U0,U1
      real(8) :: vhxc(sites),bxc(intd),vlo(dim*2),
     &           vstart(dim*2),v(dim*2),dens(dim*2),vprev(dim*2)
      complex(8) :: htest(intd,intd)
      real(8) :: tau(3),tx(sites),ty(sites),tz(sites),tplot(sites+1),
     &           bt(sites),btx(sites),bty(sites),btz(sites),
     &           bmet, dx

      write(matrix,'(a, i3, a)') '(', dim, 'f14.8)'
      write(tprint,'(a, i3, a)') '(', sites + 1, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

      call Pauli(sig)

      dx = 1d-5
      steps = 5

      v = 0.d0
      Bx = 0.d0
      by = 0.d0
      bz = 0.d0

      v(1) = 1.d0
      do i = 2, 3
        v(i) = -1.d0
      end do
      v(4) = 1.d0

      Bx(1) = .4d0
      Bx(2) = .4d0
!      Bx(3) = .01d0
!      Bx(4) = -.25d0

!      Bx = 0.d0
!      By(2) = 1.d-4

!      Bz(1) = -2.50d-1
!      Bz(2) = 2.50d-2
      Bz(3) = .4d0
      Bz(4) = .4d0


***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************
!      do i=1,sites
!        v(1) = 2.5d0
!        if (i.ne.1) then
!          v(i) = -2.5d0
!        end if
!      end do

      do i=1,sites
        v(sites+i) = Bx(i)
        v(sites*2+i) = By(i)
        v(sites*3+i) = Bz(i)
      end do

      vstart = v
      vprev = 0.d0
***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      do it=1,100

        fretlo = 1.d0
        vlo = 0.d0

        v = vstart

!        write(*,*)
        write(*,*) '************************************'
        write(*,*) it

        U0 = dble(it)/10.d0
!        U0 = 0.d0
        U1 = U0/2.d0
!        U1 = 0.d0
        ntarget = 0.d0

        call interHam(v,U0,U1,htest)
        call intdens(ntarget,htest)

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
!        number = 0
        if (it.gt.1) v=vprev
        write(*,vector) vprev

        call frprmn(v,dim*2,ftol,iter,fret)
!        flag = -10*(dim*2)-1
!        flag = dim*2
!        if (iter.eq.1) then
!          v = vprev
!          call frprmn(v,dim*2,ftol,iter,fret)
!        end if

        write(*,*) '*********************************'
        write(*,*) iter, ftol, fret

!        fretlo = 1.d0

        fretlo = fret
        vlo = v

        if (fret.gt.1d-14) then

!          call StepMin(steps, dx, v, ftol, iter, fret, fretlo)

          v = vstart
          write(*,*) '2nd minimization attempt'
          write(*,vector) v
          call frprmn(v,dim*2,ftol,iter,fret)

          if (fretlo.gt.fret) then
            fretlo = fret
            vlo = v
          end if

        end if

        if (fret.gt.1d-6) then
          call StepMin(steps, dx, v, ftol, iter, fret, fretlo)
          vlo = v
        end if

        vprev = v


***   Final check to ensure ntarget is equivalent to the density found by our
***   conjugate gradient method.
        call hbuild(v,hmat)
        call densvec(dens,hmat)

        vhxc = 0.d0
        do i=1,sites
          vhxc(i) = v(i) - vstart(i)
        end do

        call tcalc(dens,v,vstart,tau,tx,ty,tz)

        tplot(1) = u0
        do i=1, sites
          tplot(i+1) = ty(i)
        end do

        open(200,file='4pt-YTorq-4BxBz-wStep.txt')
        write(200,tprint) tplot
!        write(*,tprint) tplot
!          write(*,*) 'torque x            torque y             torque z'
!          write(*,*)

!        write(*,*) tau(1), tau(2), tau(3)
!        write(*,*) 'total tx            total ty             total tz'
!        write(*,*)

        call btrans(dens,v,vstart,bt,btx,bty,btz)

        call metric(bmet,btx,bty,btz,v,vstart)

        open(100,file='4pt-Symm-MzB4-wS.txt')

        mags(1) = u0
        do i=1, sites
          mags(i+1) = dens(sites*3+i)
        end do
        write(100,tprint) mags

        open(300,file='4pt-Symm-MxB4-wS.txt')

        do i=1, sites
          mags(i+1) = dens(sites+i)
        end do
        write(300,tprint) mags
!        write(100,*) u0,bmet,fret!,iter
        !write(*,*) u0,bmet

        write(*,*) '*** end of loop ***'
        write(*,*)
!        write(*,*)
      end do

      end
