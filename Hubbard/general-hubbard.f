      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE generalCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it,steps,bmag

      real(8) :: ftol,fret,fretlo,x
      real(8) :: Bx(sites), By(sites), Bz(sites),mags(sites+1)
      real(8) :: U0,U1
      real(8) :: vhxc(sites),bxc(intd),vlo(dim*2),vtest(dim*2),
     &           vstart(dim*2),v(dim*2),dens(dim*2),vprev(dim*2),
     &           longd(dim*2)
      complex(8) :: htest(intd,intd),hlong(dim,dim)
      real(8) :: tau(3),tx(sites),ty(sites),tz(sites),tplot(sites+1),
     &           bt(sites),btx(sites),bty(sites),btz(sites),
     &           bl(sites),blx(sites),bly(sites),blz(sites),
     &           bmet, dx, ap(sites), bta(sites), gm(sites)

      character(20) :: mwx, mwy, mwz, mwn, met

      write(matrix,'(a, i3, a)') '(', dim, 'f14.8)'
      write(tprint,'(a, i3, a)') '(', sites + 1, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

      call Pauli(sig)

      dx = 1d-5
      steps = 5


      do bmag = 1, 9

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
  !      Bx(3) = .01d0
  !      Bx(4) = -.25d0

  !      Bx = 0.d0
  !      By(2) = dx

  !      Bz(1) = -2.50d-1
  !      Bz(2) = 2.50d-2
        Bz(3) = dble(bmag/10.d0)
        Bz(4) = dble(bmag/10.d0)

        write(mwn,'(a,i1,a)') '4pt-B', bmag, '-BRot-n.txt'
        write(mwx,'(a,i1,a)') '4pt-B', bmag, '-BRot-mx.txt'
        write(mwy,'(a,i1,a)') '4pt-B', bmag, '-BRot-my.txt'
        write(mwz,'(a,i1,a)') '4pt-B', bmag, '-BRot-mz.txt'
!        write(met,'(a,i1,a)') '4pt-B', bmag, '-Metric.txt'
!        write(*,*) mwn, mwx, mwy, mwz
***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************

!        open(200,file = met)
        open(100,file = mwn)
        open(101,file = mwx)
        open(102,file = mwy)
        open(103,file = mwz)

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

            write(*,*) '********************************'
            write(*,*) iter, ftol, fret

            if (fret.lt.fretlo) then
              fretlo = fret
              vlo = v
            end if

          end if

          if (fret.gt.1d-12) then
            call StepMin(steps, dx, v, ftol, iter, fret, fretlo)
            vlo = v
            write(*,*) '*******************************'
            write(*,*) iter, ftol, fret
          end if

          vprev = vlo


***   Final check to ensure ntarget is equivalent to the density found by our
***   conjugate gradient method.
!          call hbuild(v,hmat)
          call hbuild(vlo,hmat)
          call densvec(dens,hmat)

          vhxc = 0.d0
          do i=1,sites
            vhxc(i) = vlo(i) - vstart(i)
          end do

!          call tcalc(dens,v,vstart,tau,tx,ty,tz)
          call tcalc(dens,vlo,vstart,tau,tx,ty,tz)

          tplot(1) = u0

  !        open(200,file='4pt-YTorq-4BxBz-wStep.txt')
  !        write(200,tprint) tplot
  !        write(*,tprint) tplot
  !          write(*,*) 'torque x            torque y             torque z'
  !          write(*,*)

  !        write(*,*) tau(1), tau(2), tau(3)
  !        write(*,*) 'total tx            total ty             total tz'
  !        write(*,*)

!          call btrans(dens,v,vstart,bt,btx,bty,btz)
!          call btrans(dens,vlo,vstart,bt,btx,bty,btz)

!          call metric(bmet,btx,bty,btz,v,vstart)
          call metric(bmet,btx,bty,btz,vlo,vstart)

!          call blong(dens,v,vstart,bl,blx,bly,blz)
!          call blong(dens,vlo,vstart,bl,blx,bly,blz)

          call bpar(dens,vlo,vstart,bl,blx,bly,blz)

          vtest = vlo
          do i=1,sites
            vtest(sites + i) = vtest(sites + i) + blx(i)
            vtest(sites*2 + i) = vtest(2*sites + i) + bly(i)
            vtest(sites*3 + i) = vtest(3*sites + i) + blz(i)
          end do

          call hbuild(vtest,hlong)
          call densvec(longd,hlong)

          call CalcAngle(ap, bta, gm, dens, longd)

          do i=1, sites
!            tplot(i+1) = dens(i)
            tplot(i+1) = longd(i)
          end do
          write(100,tprint) tplot

          do i=1, sites
!            tplot(i+1) = dens(sites + i)
            tplot(i+1) = longd(sites + i)
          end do
          write(101,tprint) tplot

          do i=1, sites
!            tplot(i+1) = dens(sites*2 + i)
            tplot(i+1) = longd(sites*2 + i)
          end do
          write(102,tprint) tplot

          do i=1, sites
!            tplot(i+1) = dens(sites*3 + i)
            tplot(i+1) = longd(sites*3 + i)
          end do
          write(103,tprint) tplot

          mags(1) = u0
          do i=1, sites
            mags(i+1) = dens(sites*3+i)
          end do
  !        write(100,tprint) mags

  !        open(300,file='4pt-Symm-MxB4-wS.txt')

          do i=1, sites
            mags(i+1) = dens(sites+i)
          end do

!          write(200,*) u0,bmet,fret!,iter
          !write(*,*) u0,bmet

          write(*,*) '*** end of loop ***'
          write(*,*)
  !        write(*,*)
        end do

        close(100)
        close(101)
        close(102)
        close(103)
!        close(200)

      end do

      end
