      Module ComplexCGMethod

      IMPLICIT NONE

      integer,parameter :: spin=2, sites = 2, occ=1
      integer, parameter :: dim=spin*sites!,lwork=136
      real(kind=8), parameter :: t=0.5d0
      !real(kind=8), parameter :: dab=2.d0,dbc=-3.d0,dcd=1.5d0

      complex(8) :: hmat(dim,dim), ntarget(dim*2), En(dim)
      complex(8) :: tBx(sites), tBy(sites), tBz(sites)
      complex (8), parameter :: Zero = (0.d0,0.d0), One = (1.d0,0.d0)
      complex (8), parameter :: IOne = (0.d0,1.d0)


      character (len=30) :: matrix

      contains

C  Beginning of the function being minimized through the CG method.
C  In this case, func = \int (n{v(x)} - ntarget{v(x)})**2 dx
      !      function func(diffv) result(integral)
      function func(v) result(integral)
      implicit none
      integer :: i,j
      complex(8) :: v(dim*2), dens(dim*2), integral
      complex(8) :: phi(spin,sites)
      !      write(*,matrix) v

C  Recalculating the eigenvectors through the Schrodiner Eqn.
      call hbuild(v,hmat)

      !      write(*,matrix) v

      dens = zero
      phi = zero

      call wfmatrix(phi,hmat)

      call densvec(dens,phi)

C  Calculating the difference integral between the target and current densities.
      integral = 0.d0
      do i =1,dim*2
        integral = integral + (dens(i)-ntarget(i))**2
      end do

      end function
*************************************************************************
C  Derivative subroutine of func.
C  Derivative is of the form:
C  dF(n{v(x)}) = 2 * \int (n{v(x)} - ntarget{v(x)})*(dn/dv) dx
      !subroutine dfunc(diffv,diffn)
      subroutine dfunc(v,dSdV)

      implicit none
      integer :: a,i,j,k,alpha,counter
      complex(8) :: v(dim*2),x,dens(dim*2),dSdV(dim*2),wf(spin,sites)
      complex(8) :: test(spin,spin), dn(dim*2)
      complex(8) :: dpsi(spin), dpsis(spin), dumpsi(spin)

      wf=zero
      call wfmatrix(wf,hmat)
      call densvec(dens,wf)

      alpha = occ

      dn = zero

      counter = 0

      do j=1,dim*2
        x = 0.d0
        do a=1,sites
          do i=1,spin
            dumpsi(i) = hmat((a-1)*spin+i,occ)
          end do
          call derivative(v,1,a,j,dpsi)
          dpsis = conjg(dpsi)
          do i=1,spin
            do k=1,spin
              test(i,k) = dpsi(i)*conjg(dumpsi(k))
     &                      + dumpsi(i)*dpsis(k)
            end do
          end do

          do i=1,spin
            do k=1,spin
              counter = counter + 1
              dn(counter) = test(k,i)
            end do
          end do

        end do

        do i=1,dim*2
          x = x + cmplx(2.,kind=8)*(dens(i)-ntarget(i))*dn(i)
        end do
        dSdV(j) = x
      end do

      write(*,matrix) dSdV

      end subroutine

***************************************************************************
      subroutine derivative(v,alpha,i,j,dphi)
      implicit none

      integer :: j,k,l,m,xu,xd,beta
      integer, intent(in) :: i,alpha
      complex(8) :: phi(spin)
      complex(8), intent(in) :: v(dim*2)
      complex(8) :: vmat(spin,spin), num, denom
      complex(8), intent(out) :: dphi(spin)

!      vmat(1,1) = 2.d0*one*(v(i) + v(dim*2-sites+i))
!      vmat(2,1) = 2.d0*one*(v(sites+i) + ione*v(dim+i))
!      vmat(1,2) = 2.d0*one*(v(sites+i) - ione*v(dim+i))
!      vmat(2,2) = 2.d0*one*(v(i) - v(dim*2)-sites+i)


***************************************************************************
***   xu and xd point to the proper eigenvector value inside hmat according
***   to the site being calculated. xu = spin up at position x, xd = spin
***   down at position x.
***************************************************************************
      xu = (i-1)*spin + 1
      xd = (i-1)*spin + 2

      phi = zero

      do beta = 1, dim
        if (alpha.eq.beta) then
          dphi = zero
        end if
***************************************************************************
***   Carrying out perturbation theory multiplication, using potential
***   matrix and spinors gives:
***   ( {[V0 + V3]Phi_aup + [V1-iV2]Phi*_bup + ... } * Phi_j spinor
***
***   Taking the derivative with respect to each potential, gives a set of
***   equations for potential derivative. (i.e. - Va = V(1) so j<=sites,
***   which means it would use the first equation. Likewise,
***   Bxb=V(sites + 2), putting it in the second equation.)
***************************************************************************

        if (j.le.sites) then
          num = cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xu,beta))
     &       + cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

        elseif ((j.gt.sites).and.(j.le.sites*2)) then
          num = cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xu,beta))
     &      + cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

        elseif ((j.gt.sites*2).and.(j.le.sites*3)) then
          num = cmplx(-2.,kind(8))*ione*hmat(xd,alpha)
     &                                    *conjg(hmat(xu,beta))
     &      + cmplx(2.,kind(8))*ione*hmat(xu,alpha)*conjg(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

        else
          num = cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xu,beta))
     &      - cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)
        end if


        phi(1) = phi(1) + dphi(1)
        phi(2) = phi(2) + dphi(2)

      end do

      do beta=1,spin
        dphi(beta) = phi(beta)
      end do

!      elseif (mod(j,4).eq.1) then
!        num = cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xu,beta))
!     &      + cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xd,beta))
!        denom = En(beta) - En(alpha)
!
!        dphi(1) = num/denom*hmat(xu,beta)
!        dphi(2) = num/denom*hmat(xd,beta)
!
!      elseif (mod(j,4).eq.2) then
!        num = cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xu,beta))
!     &      + cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xd,beta))
!        denom = En(beta) - En(alpha)
!
!        dphi(1) = num/denom*hmat(xu,beta)
!        dphi(2) = num/denom*hmat(xd,beta)
!
!      elseif (mod(j,4).eq.3) then
!        num = cmplx(-2.,kind(8))*ione*hmat(xd,alpha)
!     &                                       *conjg(hmat(xu,beta))
!     &      + cmplx(2.,kind(8))*ione*hmat(xu,alpha)*conjg(hmat(xd,beta))
!        denom = En(beta) - En(alpha)
!
!        dphi(1) = num/denom*hmat(xu,beta)
!        dphi(2) = num/denom*hmat(xd,beta)
!
!      elseif (mod(j,4).eq.0) then
!        num = cmplx(2.,kind(8))*hmat(xu,alpha)*conjg(hmat(xu,beta))
!     &      - cmplx(2.,kind(8))*hmat(xd,alpha)*conjg(hmat(xd,beta))
!        denom = En(beta) - En(alpha)
!
!        dphi(1) = num/denom*hmat(xu,beta)
!        dphi(2) = num/denom*hmat(xd,beta)
!      end if

      end subroutine

****************************************************************************
***   Density vector calculation subroutine
****************************************************************************
***   Calculates the spin density matrix in vectorized form. Each site is
***   listed as a 4-vector stacked on top of one another.
***   (i.e. - n(1) = |a_up|^2, n(3) = a_up x a_down*, n(5) = |b_up|^2)
****************************************************************************

      subroutine densvec(ntarget,phi)
      implicit none

      integer :: i,j,k,counter
      complex(8) :: phi(spin,sites), ntarget(dim*2)

      counter = 0
      do i=1,sites
        do k=1,spin
          do j=1,spin
            counter = counter + 1
            ntarget(counter) = dble(phi(j,i)*conjg(phi(k,i)))
          end do
        end do
      end do

      end subroutine densvec

****************************************************************************
***   Wavefunction matrix population subroutine.
****************************************************************************
***   Matrix values are wave amplitudes populated according to
***   spin and occupation site in row and column, respectively.
***   (i.e.- phi_up(a) = m_11, phi_dn(c) = m_23)
****************************************************************************
      subroutine wfmatrix(phi,h0)
      implicit none

      integer :: counter, i, j
      complex(8) :: phi(spin,sites), h0(dim,dim)

      counter = 0
      do i=1,spin
        do j=1,sites
          counter = counter + 1
          phi(i,j) = h0(counter,1)
        end do
      end do

      end subroutine wfmatrix
****************************************************************************

      subroutine hbuild(v,mat)
      implicit none

      integer, parameter :: lwork=136
      integer :: i,j,info

      real(8) :: rwork(2*dim)

      complex(8) :: v(dim), x, mat(dim,dim), wn(dim)
      complex(8) :: Bx(sites), By(sites), Bz(sites)
      complex(8) :: work(lwork),vr(dim),vl(dim),test(sites,sites)


C  Setting up noninteracting Hubbard-like system Hamiltonian with t as the
C  hopping constant between lattice sites.
      do i = 1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = v(i)*one + Bz(i)*one
          elseif (abs(i-j).eq.1) Then
            test(i,j) = -t*one
          else
            test(i,j) = zero
          endif
        end do
      end do

      mat = zero

      do i=1,sites
        do j=1,sites
          mat(i,j) = test(i,j)
        end do
      end do

      do i = 1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = v(i) - Bz(i)
          elseif (abs(i-j).eq.1) Then
            test(i,j) = -t
          else
            test(i,j) = zero
          endif
        end do
      end do

      do i=1,sites
        do j=1,sites
          mat(dim/2+i,dim/2+j) = test(i,j)
        end do
      end do

      do i=1,sites
        mat(i,dim/2+i) = Bx(i)-Ione*By(i)
        mat(dim/2+i,i) = Bx(i)+Ione*By(i)
      end do

C  Eigenvalue solver for a complex, non-symmetric matrix.

      call CGEEV('n','v',dim, mat, dim, En, VL, dim, VR, dim,
     &                             WORK, LWORK, RWORK, INFO )

      write(*,*) work(1)

C  Normalization step to normalize each orbital to 1.
!      do j = 1,dim
!        x = 0.d0
!          do i = 1, dim
!            x = x + mat(i,j)**2
!          end do
!          do i = 1, dim
!            mat(i,j) = mat(i,j)/dsqrt(x)
!          end do
!      end do

      end subroutine hbuild

****************************************************************************
****  Numerical Recipes subroutines and functions                       ****
****************************************************************************
        SUBROUTINE frprmn(p,n,ftol,iter,fret)
        INTEGER iter,n,NMAX,ITMAX
        REAL (8) :: fret,ftol,p(n),EPS
        PARAMETER (NMAX=50,ITMAX=200,EPS=1.d-10)
cU    USES dfunc,func,linmin
        INTEGER its,j
        REAL (8) :: dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)

        fp=func(p)
        call dfunc(p,xi)

        do 11 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)
          xi(j)=h(j)
11    continue
        do 14 its=1,ITMAX
          iter=its
          call linmin(p,xi,n,fret)
          if(2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
          fp=func(p)
          call dfunc(p,xi)
          gg=0.d0
          dgg=0.d0
          do 12 j=1,n
            gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
            dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
          if(gg.eq.0.d0)return
          gam=dgg/gg
          do 13 j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
            xi(j)=h(j)
13      continue
14    continue
        pause 'frprmn maximum iterations exceeded'
        return
        END
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
****************************************************************************
      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL (8) :: x, f1dim
      PARAMETER (NMAX=50)
C     USES func
      INTEGER j,ncom
      REAL (8) :: pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      f1dim=func(xt)
      return
      END function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
*****************************************************************************
        SUBROUTINE linmin(p,xi,n,fret)
        INTEGER n,NMAX
        REAL (8) :: fret,p(n),xi(n),TOL
        PARAMETER (NMAX=50,TOL=1.d-4)
CU    USES brent,f1dim,mnbrak
        INTEGER j,ncom
        REAL (8) :: ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
        COMMON /f1com/ pcom,xicom,ncom
        !real (8), EXTERNAL :: f1dim
        ncom=n
        do 11 j=1,n
          pcom(j)=p(j)
          xicom(j)=xi(j)
11    continue
        ax=0.d0
        xx=1.d0
        call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
        fret=brent(ax,xx,bx,f1dim,TOL,xmin)
        do 12 j=1,n
          xi(j)=xmin*xi(j)
          p(j)=p(j)+xi(j)
12    continue
        return
        END subroutine
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
****************************************************************************

        SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
        REAL (8) :: ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
        EXTERNAL func
        PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
        REAL (8) :: dum,fu,q,r,u,ulim
        fa=func(ax)
        fb=func(bx)
        if(fb.gt.fa)then
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
        endif
        cx=bx+GOLD*(bx-ax)
        fc=func(cx)
1     if(fb.ge.fc)then
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
          ulim=bx+GLIMIT*(cx-bx)
          if((bx-u)*(u-cx).gt.0.d0)then
            fu=func(u)
            if(fu.lt.fc)then
              ax=bx
              fa=fb
              bx=u
              fb=fu
              return
            else if(fu.gt.fb)then
              cx=u
              fc=fu
              return
            endif
            u=cx+GOLD*(cx-bx)
            fu=func(u)
          else if((cx-u)*(u-ulim).gt.0.d0)then
            fu=func(u)
            if(fu.lt.fc)then
              bx=cx
              cx=u
              u=cx+GOLD*(cx-bx)
              fb=fc
              fc=fu
              fu=func(u)
            endif
          else if((u-ulim)*(ulim-cx).ge.0.d0)then
            u=ulim
            fu=func(u)
          else
            u=cx+GOLD*(cx-bx)
            fu=func(u)
          endif
          ax=bx
          bx=cx
          cx=u
          fa=fb
          fb=fc
          fc=fu
          goto 1
        endif
        return
        END subroutine
****************************************************************************

        FUNCTION brent(ax,bx,cx,f,tol,xmin)
        INTEGER ITMAX
        REAL (8) :: brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
        EXTERNAL f
        PARAMETER (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
        INTEGER iter
        REAL (8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.
        fx=f(x)
        fv=fx
        fw=fx
        do 11 iter=1,ITMAX
          xm=0.5d0*(a+b)
          tol1=tol*dabs(x)+ZEPS
          tol2=2.d0*tol1
          if(dabs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
          if(dabs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.d0*(q-r)
            if(q.gt.0.d0) p=-p
            q=dabs(q)
            etemp=e
            e=d
      if(dabs(p).ge.dabs(.5d0*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))
     *goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
          endif
1       if(x.ge.xm) then
            e=a-x
          else
            e=b-x
          endif
          d=CGOLD*e
2       if(abs(d).ge.tol1) then
            u=x+d
          else
            u=x+sign(tol1,d)
          endif
          fu=f(u)
          if(fu.le.fx) then
            if(u.ge.x) then
              a=x
            else
              b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          else
            if(u.lt.x) then
              a=u
            else
              b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
              v=w
              fv=fw
              w=u
              fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
              v=u
              fv=fu
            endif
          endif
11    continue
        pause 'brent exceed maximum iterations'
3     xmin=x
        brent=fx
        return
      END function
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.
      End Module
