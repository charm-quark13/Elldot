      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE generalCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it,steps,bmag,
     &           wxcflg, wdflg, wmflg, weflg, wextflg

      real(8) :: ftol,fret,fretlo,x
      real(8) :: Bx(sites), By(sites), Bz(sites),mags(sites+1)
      real(8) :: U0,U1
      real(8) :: vhxc(sites),bxc(intd),vlo(dim*2),vtest(dim*2),
     &           vstart(dim*2),v(dim*2),dens(dim*2),vprev(dim*2),
     &           longd(dim*2)
      complex(8) :: htest(intd,intd),hlong(dim,dim)
      real(8) :: tau(3),tx(sites),ty(sites),tz(sites),tplot(sites+2),
     &           bt(sites),btx(sites),bty(sites),btz(sites),
     &           bl(sites),blx(sites),bly(sites),blz(sites),
     &           texx(sites), texy(sites), texz(sites),
     &           bmet, dx, theta(sites), titer(sites+1)

      character(30) :: mwx, mwy, mwz, mwn, met, twx, twy, twz,
     &                 ttotwx, ttotwy, ttotwz, ttotp

      write(matrix,'(a, i3, a)') '(', dim, 'f14.8)'
      write(tprint,'(a, i3, a)') '(', sites + 1, 'e16.6)'
      write(ttotp,'(a, i3, a)') '(', sites + 2, 'e16.6)'
      write(vector,'(a, i3, a)') '(', 1, 'f16.10)'
      write(intprint,'(a, i3, a)') '(', 6, 'e16.6)'
      write(dmat,'(a,i3,a)') '(', dim*2,'f14.10)'

      call Pauli(sig)

      dx = 1d-5
      steps = 5

      wextflg = 1
      weflg = 1
      wxcflg = 1
      wdflg = 0
      wmflg = 0

      if (wmflg.eq.1) then

        write(mwn,'(a)') '4pt-NCUnsym-Map-n.txt'
        write(mwx,'(a)') '4pt-NCUnsym-Map-x.txt'
        write(mwy,'(a)') '4pt-NCUnsym-Map-y.txt'
        write(mwz,'(a)') '4pt-NCUnsym-Map-z.txt'

        write(ttotwx,'(a)') '4pt-NCUnsym-TorqT-x.txt'
        write(ttotwy,'(a)') '4pt-NCUnsym-TorqT-y.txt'
        write(ttotwz,'(a)') '4pt-NCUnsym-TorqT-z.txt'

        open(100,file = mwn)
        open(101,file = mwx)
        open(102,file = mwy)
        open(103,file = mwz)

        open(201,file = ttotwx)
        open(202,file = ttotwy)
        open(203,file = ttotwz)

      end if

      do bmag = 1, 1!1
!      do bmag = 2, 10

        v = 0.d0
        Bx = 0.d0
        by = 0.d0
        bz = 0.d0

        v(1) = 1.d0
        do i = 2, 3
          v(i) = -1.d0
        end do
        v(4) = 1.d0

        Bx(1) = 2d-1
!        Bx(1) = dble(bmag/10.d0)
        Bx(2) = dble(bmag/10.d0)
!        Bx(4) = -dble(bmag/5.d0)

!        By(1) = dble(bmag/10.d0)
!        By(3) = -dble(bmag/10.d0)
!        By(2) = -dble(bmag/10.d1)
!        By(3) = -dble(bmag/10.d1)

!        Bz(3) = dble(bmag/10.d0)
!        Bz(4) = dble(bmag/10.d0)
        Bz(3) = 3d-1
        Bz(4) = -2d-1
!        Bz(3) = dble(bmag/10.d1)
!        Bz(4) = dble(bmag/10.d1)



        if (bmag.eq.11) then

          Bx(1) = 1.d-2
          Bx(2) = 1.d-2
          Bx(4) = -1.d-1
!          Bx(1) = 1.d-3
!          Bx(2) = 1.d-3

          By(1) = 1.d-2
          By(3) = -1.d-2
!          By(2) = -1.d-3
!          By(3) = -1.d-3

          Bz(3) = 1.d-2
          Bz(4) = 1.d-2
!          Bz(3) = 1.d-3
!          Bz(4) = 1.d-3

        end if

!        write(met,'(a,i1,a)') '4pt-B', bmag, '-CNC-Metric.txt'
!        write(*,*) mwn, mwx, mwy, mwz
***************************************************************************
***   Setting the initial potentials which will generate our target density.
***************************************************************************

!        open(200,file = met)

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
        do it=0,100

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

          if (fret.gt.1d-6) then
            call StepMin(steps, dx, v, ftol, iter, fret, fretlo)
            vlo = v
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
          call tortc(dens, vstart, texx, texy, texz)

          tplot(1) = u0
          tplot(2) = bx(1)
          write(*,*) tplot(2)

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
          call blong(dens,vlo,vstart,bl,blx,bly,blz)

!          call bpar(dens,vlo,vstart,bl,blx,bly,blz)

          vtest = vlo
          do i=1,sites
            vtest(sites + i) = vtest(sites + i) + blx(i)
            vtest(sites*2 + i) = vtest(sites*2 + i) + bly(i)
            vtest(sites*3 + i) = vtest(sites*3 + i) + blz(i)
          end do

          call hbuild(vtest,hlong)
          call densvec(longd,hlong)

!          call CalcAngle(theta, dens, longd)

!          write(mwx,'(a,i1,a)') '4pt-B', bmag, '-CNC-theta.txt'

          if (wdflg.eq.1) then

            if (bmag.lt.10) then
              write(mwn,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-n.txt'
              write(mwx,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-mx.txt'
              write(mwy,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-my.txt'
              write(mwz,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-mz.txt'
            elseif (bmag.eq.10) then
              write(mwn,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-n.txt'
              write(mwx,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-mx.txt'
              write(mwy,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-my.txt'
              write(mwz,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-mz.txt'
            elseif (bmag.eq.11) then
              write(mwn,'(a)') '4pt-B01-XC-Exact-CP-n.txt'
              write(mwx,'(a)') '4pt-B01-XC-Exact-CP-mx.txt'
              write(mwy,'(a)') '4pt-B01-XC-Exact-CP-my.txt'
              write(mwz,'(a)') '4pt-B01-XC-Exact-CP-mz.txt'
            end if

            open(30,file = mwn)
            open(31,file = mwx)
            open(32,file = mwy)
            open(33,file = mwz)

          end if

          if (wxcflg.eq.1) then

            if (bmag.lt.10) then
              write(twx,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-tx.txt'
              write(twy,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-ty.txt'
              write(twz,'(a,i1,a)') '4pt-B', bmag, '-XC-Exact-CP-tz.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-tx.txt'
              write(twy,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-ty.txt'
              write(twz,'(a,i2,a)') '4pt-B', bmag, '-XC-Exact-CP-tz.txt'
            else
              write(twx,'(a)') '4pt-B01-XC-Exact-CP-tx.txt'
              write(twy,'(a)') '4pt-B01-XC-Exact-CP-ty.txt'
              write(twz,'(a)') '4pt-B01-XC-Exact-CP-tz.txt'
            end if

            open(301,file = twx)
            open(302,file = twy)
            open(303,file = twz)

          end if

!          open(201,file = mwx)

          titer(1) = u0

          if (wdflg.eq.1) then
            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(30,tprint) titer
          end if

          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(i) - longd(i)
!              tplot(i+1) = longd(i)
            end do

            write(100,ttotp) tplot
          endif


          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(sites + i) - longd(sites + i)
!              tplot(i+1) = longd(sites + i)
            end do

            write(101,ttotp) tplot

            do i=1, sites
              tplot(i+2) = tx(i)
            end do
            write(201,ttotp) tplot
          endif

          if (wdflg.eq.1) then
            do i=1, sites
              titer(i+1) = longd(sites + i)
            end do
            write(31,tprint) titer
          end if

          if (wxcflg.eq.1) then
            do i=1, sites
              titer(i+1) = tx(i)
            end do
            write(301,tprint) titer
          end if

!          do i=1, sites
!            tplot(i+2) = theta(i)
!          end do
!          write(201,tprint) tplot


          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(sites*2 + i) - longd(sites*2 + i)
!              tplot(i+1) = longd(sites*2 + i)
            end do

            write(102,tprint) tplot

            do i=1, sites
              tplot(i+2) = ty(i)
            end do
            write(202,ttotp) tplot

          end if

          if (wdflg.eq.1) then
            do i=1, sites
              titer(i+1) = longd(sites*2 + i)
            end do
            write(32,tprint) titer
          end if

          if (wxcflg.eq.1) then
            do i=1, sites
              titer(i+1) = ty(i)
            end do
            write(302,tprint) titer
          end if

          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(sites*3 + i) - longd(sites*3 + i)
!              tplot(i+1) = longd(sites*3 + i)
            end do

            write(103,tprint) tplot

            do i=1, sites
              tplot(i+2) = tz(i)
            end do
            write(203,ttotp) tplot
          end if

          if (wdflg.eq.1) then
            do i=1, sites
              titer(i+1) = longd(sites*3 + i)
            end do
            write(33,tprint) titer
          end if

          if (wxcflg.eq.1) then
            do i=1, sites
              titer(i+1) = tz(i)
            end do
            write(303,tprint) titer
          end if

          if (wextflg.eq.1) then

            if (bmag.lt.10) then
              write(twx,'(a,i1,a)') '4pt-B', bmag, '-Ext-CP-tx.txt'
              write(twy,'(a,i1,a)') '4pt-B', bmag, '-Ext-CP-ty.txt'
              write(twz,'(a,i1,a)') '4pt-B', bmag, '-Ext-CP-tz.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') '4pt-B', bmag, '-Ext-CP-tx.txt'
              write(twy,'(a,i2,a)') '4pt-B', bmag, '-Ext-CP-ty.txt'
              write(twz,'(a,i2,a)') '4pt-B', bmag, '-Ext-CP-tz.txt'
            elseif (bmag.eq.11) then
              write(twx,'(a)') '4pt-B01-Ext-CP-tx.txt'
              write(twy,'(a)') '4pt-B01-Ext-CP-ty.txt'
              write(twz,'(a)') '4pt-B01-Ext-CP-tz.txt'
            end if

            open(401,file = twx)
            open(402,file = twy)
            open(403,file = twz)

            do i=1, sites
              titer(i+1) = texx(i)
            end do
            write(401,tprint) titer

            do i=1, sites
              titer(i+1) = texy(i)
            end do
            write(402,tprint) titer

            do i=1, sites
              titer(i+1) = texz(i)
            end do
            write(403,tprint) titer

          end if

          if (weflg.eq.1) then

            if (bmag.lt.10) then
              write(twx,'(a,i1,a)') '4pt-B', bmag, '-Exact-CP-E.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') '4pt-B', bmag, '-Exact-CP-E.txt'
            elseif (bmag.eq.11) then
              write(twx,'(a)') '4pt-B01-Exact-CP-E.txt'
            end if

            open(11,file = twx)

            do i=1, sites
              titer(i+1) = en(i)
            end do
            write(11,tprint) titer

          end if

          write(*,*) '*** end of loop ***'
          write(*,*)
!        write(*,*)
        end do

        if (wxcflg.eq.1) then
          close(301)
          close(302)
          close(303)
        end if

        if (wdflg.eq.1) then
          close(31)
          close(32)
          close(33)
        end if

        if (wextflg.eq.1) then
          close(401)
          close(402)
          close(403)
        end if

        if (weflg.eq.1) then
          close(11)
        end if

      end do

      if (wmflg.eq.1) then
        close(100)
        close(101)
        close(102)
        close(103)

        close(201)
        close(202)
        close(203)
      endif

      end
