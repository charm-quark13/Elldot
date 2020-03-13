      Program Inversion
C  Here I use a module to allow for the passing of global variables as well
C  as portability of the code to future programs.
      USE coplanarCG
!      USE interactingHD

      Implicit none

      integer :: i,j,k,iter,it,steps,bmag,
     &           wxcflg, wdflg, wmflg, weflg, wextflg, txcflg, wbflg,
     &           wlflg

      real(8) :: ftol,fret,fretlo,x
      real(8) :: Bx(sites), By(sites), Bz(sites),mags(sites+1)
      real(8) :: U0,U1, dU, delta
      real(8) :: vhxc(dim*2),bxc(intd),vlo(dim*2),vtest(dim*2),
     &           vstart(dim*2),v(dim*2),dens(dim*2),vprev(dim*2),
     &           longd(dim*2)
      complex(8) :: htest(intd,intd),hlong(dim,dim), dumh(intd,intd)
      real(8) :: tau(3),tx(sites),ty(sites),tz(sites),tplot(sites+2),
     &           bt(sites),btx(sites),bty(sites),btz(sites),
     &           bl(sites),blx(sites),bly(sites),blz(sites),
     &           texx(sites), texy(sites), texz(sites), cn(intd),
     &           bmet, dx, theta(sites), titer(sites+1), mdotb(sites+1)

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

      dU = 0.1d0

      txcflg = 0
      wextflg = 1
      weflg = 1
      wxcflg = 1
      wdflg = 1
      wmflg = 0
      wbflg = 1
      wlflg = 0

      if (wmflg.eq.1) then

        write(mwn,'(a)') 'Test-Asym-Map-n.txt'
        write(mwx,'(a)') 'Test-Asym-Map-x.txt'
        write(mwy,'(a)') 'Test-Asym-Map-y.txt'
        write(mwz,'(a)') 'Test-Asym-Map-z.txt'

        write(ttotwx,'(a)') 'Test-Asym-TorqT-x.txt'
        write(ttotwy,'(a)') 'Test-Asym-TorqT-y.txt'
        write(ttotwz,'(a)') 'Test-Asym-TorqT-z.txt'

        open(100,file = mwn)
        open(101,file = mwx)
        open(102,file = mwy)
        open(103,file = mwz)

        open(201,file = ttotwx)
        open(202,file = ttotwy)
        open(203,file = ttotwz)

      end if

      do bmag = 1, 1
!      do bmag = 2, 10

        v = 0.d0
        Bx = 0.d0
        by = 0.d0
        bz = 0.d0

        v(1) = 1.d0
        do i = 2, 3
          v(i) = -1.d0
        end do
        v(sites) = 1.d0


        Bx(1) = 2d-1
        Bx(2) = dble(bmag/10.d0)

        Bz(3) = 3d-1
        Bz(4) = -2d-1


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

!        write(met,'(a,i1,a)') 'Test-B', bmag, '-CNC-Metric.txt'
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
        U0 = 0.d0
        do it=0, 100

!          if (it.lt.80.and.it.eq.0) then
!            U0 = U0
!          elseif (it.lt.80.and.it.ne.0) then
!            U0 = U0 + dU
!          else
!            U0 = U0 + dU/2.d0
!          endif

          fretlo = 1.d0
          vlo = 0.d0

          v = vstart

          U0 = dble(it)/10.d0

!        write(*,*)
          write(*,*) '************************************'
          write(*,*) U0

!          U0 = dble(it)/10.d-1
!        U0 = 0.d0
          U1 = U0/2.d0
!        U1 = 0.d0
          ntarget = 0.d0

          call interHam(v,U0,U1,htest,cn)
          call intdens(ntarget,htest)

***************************************************************************
***   Calling Numerical Recipes' conjugate gradient method optimiztion
***   subroutine.
***************************************************************************
!        number = 0
          if (it.gt.0) v=vprev
!          write(*,vector) vprev

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

          if (fret.gt.1d-10) then
            v=vprev
            call StepMin(steps, dx, v, ftol, iter, fret, fretlo)
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

!          write(*,vector) vlo
!          write(*,*) '******** ^^^^^^^ v ^^^^^^^ *************'
!          write(*,vector) dens
!          write(*,*) '******** ^^^^^^^ m ^^^^^^^ *************'
!          write(*,vector) vlo - vstart
!          write(*,*) '******** ^^^^^^^ vs ^^^^^^^ *************'


          vhxc = 0.d0
          do i=1,sites*4
            vhxc(i) = vlo(i) - vstart(i)
          end do

          if (wbflg.eq.1) then
            if (bmag.lt.10) then
              write(mwn,'(a,i1,a)') 'Test-B', bmag, '-Asym-vxc.txt'
              write(mwx,'(a,i1,a)') 'Test-B', bmag, '-Asym-bxcx.txt'
              write(mwy,'(a,i1,a)') 'Test-B', bmag, '-Asym-bxcy.txt'
              write(mwz,'(a,i1,a)') 'Test-B', bmag, '-Asym-bxcz.txt'
            elseif (bmag.eq.10) then
              write(mwx,'(a,i2,a)') 'Test-B', bmag, '-Asym-U20-bxcx.txt'
              write(mwy,'(a,i2,a)') 'Test-B', bmag, '-Asym-U20-bxcy.txt'
              write(mwz,'(a,i2,a)') 'Test-B', bmag, '-Asym-U20-bxcz.txt'
            elseif (bmag.eq.11) then
              write(mwx,'(a)') 'Test-B01-Asym-U20-bxcx.txt'
              write(mwy,'(a)') 'Test-B01-Asym-U20-bxcy.txt'
              write(mwz,'(a)') 'Test-B01-Asym-U20-bxcz.txt'
            end if

            open(110,file = mwn)
            open(111,file = mwx)
            open(112,file = mwy)
            open(113,file = mwz)

            titer(1) = U0

            do i = 1, sites
              titer(i+1) = vhxc(i)
            end do
            write(110,tprint) titer

            do i = 1, sites
              titer(i+1) = vhxc(sites + i)
            end do
            write(111,tprint) titer

            do i = 1, sites
              titer(i+1) = vhxc(sites*2 + i)
            end do
            write(112,tprint) titer

            do i = 1, sites
              titer(i+1) = vhxc(sites*3 + i)
            end do
            write(113,tprint) titer

            mdotb(1) = u0
            open(99, file = 'Test-B1-Asym-mdotb.txt')
            do i = 1, sites
              mdotb(i+1) = dens(sites + i)*vhxc(sites + i)
     &                   + dens(sites*2 + i)*vhxc(sites*2 + i)
     &                   + dens(sites*3 + i)*vhxc(sites*3 + i)
            end do

            write(99,tprint) mdotb

          end if

!          write(*,vector) vhxc
!          write(*,*) '*************** vhxc ^^^^^^^^^^^^ *********'
!          call tcalc(dens,v,vstart,tau,tx,ty,tz)
          call tcalc(dens,vlo,vstart,tau,tx,ty,tz)
          call tortc(dens, vstart, texx, texy, texz)

          tplot(1) = u0
          tplot(2) = bx(1)
!          write(*,*) tplot(2)

!        open(200,file='Test-YTorq-4BxBz-wStep.txt')
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
!          call metric(bmet,btx,bty,btz,vlo,vstart)

!          call blong(dens,v,vstart,bl,blx,bly,blz)
          call blong(dens,vlo,vstart,bl,blx,bly,blz)

!          call bpar(dens,vlo,vstart,bl,blx,bly,blz)

!          vtest = vlo
          vtest = 0.d0
          do i=1,sites
            vtest(i) = vstart(i) + vhxc(i)
            vtest(sites + i) = vstart(sites + i) + blx(i)
            vtest(sites*2 + i) = vstart(sites*2 + i) + bly(i)
            vtest(sites*3 + i) = vstart(sites*3 + i) + blz(i)
          end do

          if (wlflg.eq.1) then

            call hbuild(vtest,hlong)
            call densvec(longd,hlong)

            if (bmag.lt.10) then
            write(mwn,'(a,i1,a)') 'Test-B', bmag, '-SymCP-nlong-t.txt'
            write(mwx,'(a,i1,a)') 'Test-B', bmag, '-SymCP-mxlong-t.txt'
            write(mwy,'(a,i1,a)') 'Test-B', bmag, '-SymCP-mylong-t.txt'
            write(mwz,'(a,i1,a)') 'Test-B', bmag, '-SymCP-mzlong-t.txt'
            elseif (bmag.eq.10) then
            write(mwn,'(a,i2,a)') 'Test-B', bmag, '-SymCP-nlong-t.txt'
            write(mwx,'(a,i2,a)') 'Test-B', bmag, '-SymCP-mxlong-t.txt'
            write(mwy,'(a,i2,a)') 'Test-B', bmag, '-SymCP-mylong-t.txt'
            write(mwz,'(a,i2,a)') 'Test-B', bmag, '-SymCP-mzlong-t.txt'
            elseif (bmag.eq.11) then
              write(mwn,'(a)') 'Test-B01-SymCP-n.txt'
              write(mwx,'(a)') 'Test-B01-SymCP-mx.txt'
              write(mwy,'(a)') 'Test-B01-SymCP-my.txt'
              write(mwz,'(a)') 'Test-B01-SymCP-mz.txt'
            end if

            open(500,file = mwn)
            open(501,file = mwx)
            open(502,file = mwy)
            open(503,file = mwz)

            open(504,file = 'Test-SymCP-Elong.txt')

            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(500, tprint) titer

            do i=1, sites
              titer(i+1) = longd(sites+i)
            end do
            write(501, tprint) titer

            do i=1, sites
              titer(i+1) = longd(sites*2+i)
            end do
            write(502, tprint) titer

            do i=1, sites
              titer(i+1) = longd(sites*3+i)
            end do
            write(503, tprint) titer

            call interHam(vtest,0.d0,U1,dumh,cn)
            write(504, *) u0, cn(1)

          end if


!          call CalcAngle(theta, dens, longd)

!          write(mwx,'(a,i1,a)') 'Test-B', bmag, '-CNC-theta.txt'

          if (wdflg.eq.1) then

            if (bmag.lt.10) then
              write(mwn,'(a,i1,a)') 'Test-B', bmag, '-Asym-n.txt'
              write(mwx,'(a,i1,a)') 'Test-B', bmag, '-Asym-mx.txt'
              write(mwy,'(a,i1,a)') 'Test-B', bmag, '-Asym-my.txt'
              write(mwz,'(a,i1,a)') 'Test-B', bmag, '-Asym-mz.txt'
            elseif (bmag.eq.10) then
              write(mwn,'(a,i2,a)') 'Test-B', bmag, '-Asym-n.txt'
              write(mwx,'(a,i2,a)') 'Test-B', bmag, '-Asym-mx.txt'
              write(mwy,'(a,i2,a)') 'Test-B', bmag, '-Asym-my.txt'
              write(mwz,'(a,i2,a)') 'Test-B', bmag, '-Asym-mz.txt'
            elseif (bmag.eq.11) then
              write(mwn,'(a)') 'Test-B01-Asym-n.txt'
              write(mwx,'(a)') 'Test-B01-Asym-mx.txt'
              write(mwy,'(a)') 'Test-B01-Asym-my.txt'
              write(mwz,'(a)') 'Test-B01-Asym-mz.txt'
            end if

            open(30,file = mwn)
            open(31,file = mwx)
            open(32,file = mwy)
            open(33,file = mwz)

          end if

          if (wxcflg.eq.1) then

            if (bmag.lt.10) then
              write(twx,'(a,i1,a)') 'Test-B', bmag,'-XC-Asym-tx.txt'
              write(twy,'(a,i1,a)') 'Test-B', bmag,'-XC-Asym-ty.txt'
              write(twz,'(a,i1,a)') 'Test-B', bmag,'-XC-Asym-tz.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') 'Test-B', bmag,'-XC-Asym-tx.txt'
              write(twy,'(a,i2,a)') 'Test-B', bmag,'-XC-Asym-ty.txt'
              write(twz,'(a,i2,a)') 'Test-B', bmag,'-XC-Asym-tz.txt'
            else
              write(twx,'(a)') 'Test-B01-XC-Asym-tx.txt'
              write(twy,'(a)') 'Test-B01-XC-Asym-ty.txt'
              write(twz,'(a)') 'Test-B01-XC-Asym-tz.txt'
            end if

            open(301,file = twx)
            open(302,file = twy)
            open(303,file = twz)

          end if

!          open(201,file = mwx)

          titer(1) = u0

          if (wdflg.eq.1) then
            do i=1, sites
              titer(i+1) = dens(i)
            end do
            write(30,tprint) titer
          end if

          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(i)
!              tplot(i+1) = longd(i)
            end do

            write(100,ttotp) tplot
          endif


          if (wmflg.eq.1) then
            do i=1, sites
              tplot(i+2) = dens(sites + i)
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
              titer(i+1) = dens(sites + i)
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
              tplot(i+2) = dens(sites*2 + i)
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
              titer(i+1) = dens(sites*2 + i)
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
              tplot(i+2) = dens(sites*3 + i)
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
              titer(i+1) = dens(sites*3 + i)
            end do
            write(33,tprint) titer
          end if

          if (wxcflg.eq.1) then
            do i=1, sites
              titer(i+1) = tz(i)
            end do
            write(303,tprint) titer
          end if

          if (txcflg.eq.1) then

            write(mwn,'(a,i1,a)') 'Test-B1-Asym-0txc-n.txt'
            write(mwx,'(a,i1,a)') 'Test-B1-Asym-0txc-mx.txt'
            write(mwy,'(a,i1,a)') 'Test-B1-Asym-0txc-my.txt'
            write(mwz,'(a,i1,a)') 'Test-B1-Asym-0txc-mz.txt'

            open(900, file=mwn)
            open(901, file=mwx)
            open(902, file=mwy)
            open(903, file=mwz)

            titer(1) = u0

            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(900,tprint) titer

            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(901,tprint) titer

            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(902,tprint) titer

            do i=1, sites
              titer(i+1) = longd(i)
            end do
            write(903,tprint) titer

          end if


          if (wextflg.eq.1) then

            if (bmag.lt.10) then
              write(twx,'(a,i1,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-tx.txt'
              write(twy,'(a,i1,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-ty.txt'
              write(twz,'(a,i1,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-tz.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-tx.txt'
              write(twy,'(a,i2,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-ty.txt'
              write(twz,'(a,i2,a)') 'Test-B', bmag,
     &                                      '-Ext-Asym-tz.txt'
            elseif (bmag.eq.11) then
              write(twx,'(a)') 'Test-B01-Ext-Asym-tx.txt'
              write(twy,'(a)') 'Test-B01-Ext-Asym-ty.txt'
              write(twz,'(a)') 'Test-B01-Ext-Asym-tz.txt'
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
              write(twx,'(a,i1,a)') 'Test-B', bmag,
     &                                      '-Exact-Asym-E.txt'
            elseif (bmag.eq.10) then
              write(twx,'(a,i2,a)') 'Test-B', bmag,
     &                                      '-Exact-Asym-E.txt'
            elseif (bmag.eq.11) then
              write(twx,'(a)') 'Test-B01-Exact-Asym-E.txt'
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

        if (txcflg.eq.1) then
          close(900)
          close(901)
          close(902)
          close(903)
        end if

        if (wbflg.eq.1) then
          close(110)
          close(111)
          close(112)
          close(113)
          close(99)
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
