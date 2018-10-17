      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE generalCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it

      real(8) :: ftol,fret,x
      real(8) :: Bx(sites), By(sites), Bz(sites)
      real(8) :: U0,U1
      real(8) :: vhxc(sites),bxc(intd),
     &           vstart(dim*2),v(dim*2),dens(dim*2),vprev(dim*2)
      complex(8) :: htest(intd,intd)
      real(8) :: tau(3),tx(sites),ty(sites),tz(sites),
     &           bt(sites),btx(sites),bty(sites),btz(sites),
     &           bmet

      write(matrix,'(a, i3, a)') '(', dim, 'f14.8)'
      write(matprint,'(a, i3, a)') '(', sites, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

      call Pauli(sig)

      v = 0.d0
      Bx = 0.d0
      by = 0.d0
      bz = 0.d0

      v(1) = 1.d0
      v(4) = -1.d0

      Bx(1) = -0.5d0
      Bx(2) = 1.0d0
      Bx(3) = .5d0
      Bx(4) = -1.d0

!      Bx = 0.d0
!      By(2) = .75d0
      Bz(1) = -1.0d0
      Bz(2) = 1.0d0
      Bz(3) = -1.d0
      Bz(4) = 1.d0

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

***************************************************************************
***   Solving the initial Schrodinger equation for our system via hbuild.
***************************************************************************
      do it=18,18

        v = vstart

        write(*,*) it

        U0 = dble(it)/10.d0
        U1 = U0/dsqrt(2.d0)
        ntarget = 0.d0

        call interHam(v,U0,U1,htest)
        call intdens(ntarget,htest)

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
!        number = 0
        !if (it.gt.1) v=vprev

        call frprmn(v,dim*2,ftol,iter,fret)

        if (iter.eq.1) then
          v = vstart
          call frprmn(v,dim*2,ftol,iter,fret)
        end if

        vprev = v

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
        do i=1,sites
          vhxc(i) = v(i) - vstart(i)
        end do

        call tcalc(dens,v,vstart,tau,tx,ty,tz)

        call btrans(dens,v,vstart,bt,btx,bty,btz)

        call metric(bmet,btx,bty,btz,v,vstart)

        open(100,file='metric-4pt.txt')

        !write(100,*) u0,bmet!,fret,iter
        write(*,*) u0,bmet

        write(*,*) '*** end of loop ***'
      end do

      end
