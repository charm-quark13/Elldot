      Module realCGMethod

      IMPLICIT NONE

      integer,parameter :: spin=2, sites = 2, occ=1
      integer, parameter :: dim=spin*sites!,lwork=136
      real(kind=8), parameter :: t=0.5d0
      !real(kind=8), parameter :: dab=2.d0,dbc=-3.d0,dcd=1.5d0

      real(8) :: hmat(dim,dim), ntarget(dim*2), En(dim)
      real(8) :: tBx(sites), tBy(sites), tBz(sites)
      real(8) :: sig(spin*4,spin*4)
!      real (8), parameter :: Zero = (0.d0,0.d0), One = (1.d0,0.d0)
!      real (8), parameter :: IOne = (0.d0,1.d0)


      character (len=30) :: matrix, matprint, vector

      contains

C  Beginning of the function being minimized through the CG method.
C  In this case, func = \int (n{v(x)} - ntarget{v(x)})**2 dx
!      function func(diffv) result(integral)
      function func(v) result(integral)
      implicit none
      integer :: i,j
      real(8) :: v(dim*2), dens(dim*2), integral
      real(8) :: phi(spin,sites)
!      write(*,matrix) v

      do i=3,6
        v(i) = 0.d0
      end do

      write(*,vector) v

C  Recalculating the eigenvectors through the Schrodinger Eqn.
      call hbuild(v,hmat)

      !write(*,*) 'En             '
      !write(*,vector) En
      !write(*,*) ' '
      !write(*,matrix) transpose(hmat)

      dens = 0.d0
      phi = 0.d0

      call wfmatrix(phi,hmat)

      call densvec(dens,phi)
      !write(*,*) '^^^^^^^^^^^^^^^^^^^func nvec'


***   Using normalized density to our advantage to reduce the dimensionality of
***   the problem. This allows us to pin Bz(b) to zero. We can do the same with
***   the density, leading us to just the differene, or n(b) = 1 - n(a).
      dens(2) = 1-dens(1)
      dens(sites*4) = 0.d0
      do i=sites*2+1,sites*3
        dens(i) = 0.d0
      end do

C  Calculating the difference integral between the target and current densities.
      integral = 0.d0
      do i =1,dim*2
        integral = integral + (dens(i)-ntarget(i))**2
      end do

!      write(*,*) integral

      end function
*************************************************************************
C  Derivative subroutine of func.
C  Derivative is of the form:
C  dF(n{v(x)}) = 2 * \int (n{v(x)} - ntarget{v(x)})*(dn/dv) dx
      !subroutine dfunc(diffv,diffn)
      subroutine dfunc(v,dSdV)

      implicit none
      integer :: a,i,j,k,alpha,counter
      real(8) :: v(dim*2),x,dens(dim*2),dSdV(dim*2),wf(spin,sites)
      real(8) :: dn(dim*2), dpmat(spin,sites*2)!, wfmat(spin,sites)
      real(8) :: dpsi(spin), dpsis(spin), dumpsi(spin)

      wf=0.d0
      call wfmatrix(wf,hmat)
      call densvec(dens,wf)

      alpha = occ

      dn = 0.d0
      dpmat = 0.d0

      do j=1,dim*2
        x = 0.d0
        counter = 0
        do a=1,sites
          do i=1,spin
            dumpsi(i) = hmat((i-1)*spin+a,occ)
          end do
          call derivative(v,1,a,j,dpsi)
!          dpsis = conjg(dpsi)


          do i=1,spin
            dpmat(i,(a-1)*2 + 1) = dumpsi(i)
            dpmat(i,(a-1)*2 + 2) = dpsi(i)
          end do

          !call collect(a,dumpsi,dpsi,dpmat)

!          do i=1,spin
!            do k=1,spin
!              dnmat(i,k) = dpsi(i)*dumpsi(k)
!     &                      + dumpsi(i)*dpsi(k)
!            end do
!          end do

!          write(*,matprint) test

!          do i=1,spin
!            do k=1,spin
!              counter = counter + 1
!              dn(counter) = dnmat(k,i)
!            end do

        end do

        write(*,*) '******', j, '*********'
        !write(*,*) '******', a, '*********'
        write(*,matrix) transpose(dpmat)


        call dnvec(dn,dpmat)

        !write(*,vector) dens
        !write(*,*)'***************'

        do i=1,dim*2
          x = x + (dens(i)-ntarget(i))*dn(i)
        end do

        write(1,*) '**********',j,'*********'
        do i = 1, dim*2
          write(1,*) dens(i) - ntarget(i)
        end do
        write(1,*) '************************'
        do i = 1, dim*2
          write(1,*) dn(i)
        end do
        write(1,*) '************************'

!        if ((j.ge.3).and.(j.le.6)) then
!          dSdV(j) = 0.d0
!        else
          dSdV(j) = 2.d0*x
!        endif

      end do

      write(1,vector) dSdV
      write(1,*) '************************'
      end subroutine

***************************************************************************
      subroutine derivative(v,alpha,i,j,dphi)
      implicit none

      integer :: j,k,l,m,xu,xd,beta
      integer, intent(in) :: i,alpha
      real(8) :: phi(spin)
      real(8), intent(in) :: v(dim*2)
      real(8) :: vmat(spin,spin), num, denom
      real(8), intent(out) :: dphi(spin)

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
      xd = xu + spin

      phi = 0.d0
      dphi = 0.d0

      do beta = 1, dim
        if (alpha.eq.beta) then
          dphi = 0.d0

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

        elseif (j.le.sites) then
          num = 2.d0*hmat(xu,alpha)*(hmat(xu,beta))
     &       + 2.d0*hmat(xd,alpha)*(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = (num/denom)*hmat(xu,beta)
          dphi(2) = (num/denom)*hmat(xd,beta)

          !write(*,*) num, denom

        elseif ((j.gt.sites).and.(j.le.sites*2)) then
          num = 2.d0*hmat(xd,alpha)*(hmat(xu,beta))
     &      + 2.d0*hmat(xu,alpha)*(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

          !write(*,*) num, denom

        elseif ((j.gt.sites*2).and.(j.le.sites*3)) then
          num = 2.d0*hmat(xd,alpha)*(hmat(xu,beta))
     &      + 2.d0*hmat(xu,alpha)*(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

          !write(*,*) num, denom

        else
          num = 2.d0*hmat(xu,alpha)*(hmat(xu,beta))
     &      - 2.d0*hmat(xd,alpha)*(hmat(xd,beta))
          denom = En(beta) - En(alpha)

          dphi(1) = num/denom*hmat(xu,beta)
          dphi(2) = num/denom*hmat(xd,beta)

          !write(*,*) num, denom

        end if


        phi(1) = phi(1) + dphi(1)
        phi(2) = phi(2) + dphi(2)

      end do

      dphi = phi


      end subroutine

***************************************************************************

      subroutine dnvec(n0,phi)
      implicit none

      integer :: i,j,k,m
      real(8) :: x
      real(8) :: phi(spin,sites*2), n0(dim*2), test(spin,spin)
      real(8) :: sigma(spin,spin), dum(spin,spin)

      do i=1,sites
        do j=1,spin
          do k=1,spin
            test(j,k) = phi(j,(i-1)*2+2)*phi(k,(i-1)*2+1)
     &                      + phi(k,(i-1)*2+1)*phi(j,(i-1)*2+2)
          end do
        end do


**************************************************************************
***   Beginning of calculation of n, mx, my, mz. The index m runs from 1-
***   to-4 for the density and then each magnetization direction.
**************************************************************************

        do m=0,3
          do j=1,spin
            do k=1,spin
              sigma(j,k) = sig(m*spin+j,m*spin+k)
            end do
          end do

          dum = matmul(sigma,test)

          x = 0.d0
          do j=1,spin
            do k=1,spin
              if (j.eq.k) then
                x = x + dum(j,k)
              endif
            end do
          end do

          n0(i+(m)*2) = x

        end do

      end do

!      write(*,vector) n0

      end subroutine dnvec

****************************************************************************
***   Density vector calculation subroutine
****************************************************************************
***   Calculates the spin density matrix in vectorized form. Each site is
***   listed as a 4-vector stacked on top of one another.
***   (i.e. - n(1) = |a_up|^2, n(3) = a_up x a_down*, n(5) = |b_up|^2)
****************************************************************************

      subroutine densvec(n0,phi)
      implicit none

      integer :: i,j,k,m
      real(8) :: x
      real(8) :: phi(spin,sites), n0(dim*2), test(spin,spin)
      real(8) :: sigma(spin,spin), dum(spin,spin)

      do i=1,sites
        do j=1,spin
          do k=1,spin
            test(j,k) = phi(j,i)*phi(k,i)
          end do
        end do


**************************************************************************
***   Beginning of calculation of n, mx, my, mz. The index m runs from 1-
***   to-4 for the density and then each magnetization direction.
**************************************************************************

        do m=0,3
          do j=1,spin
            do k=1,spin
              sigma(j,k) = sig(m*spin+j,m*spin+k)
            end do
          end do

          dum = matmul(sigma,test)

          x = 0.d0
          do j=1,spin
            do k=1,spin
              if (j.eq.k) then
                x = x + dum(j,k)
              endif
            end do
          end do

          n0(i+m*2) = x

        end do

      end do

      !write(*,*) '*********'
      !write(*,vector)n0
      !write(*,*) '^^^^^^^^^^^ntarget'
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
      real(8) :: phi(spin,sites), h0(dim,dim)

      counter = 0
      do i=1,spin
        do j=1,sites
          counter = counter + 1
          phi(i,j) = h0(counter,1)
        end do
      end do

      end subroutine wfmatrix

****************************************************************************

      subroutine collect(x,phi,dphi,wfmat)
      implicit none

      integer :: j,x
      real(8) :: phi(spin),dphi(spin),wfmat(spin,sites*spin)

      do j=1,spin
        wfmat(j,(x-1)*2 + 1) = phi(j)
        wfmat(j,(x-1)*2 + 2) = dphi(j)
      end do

      end subroutine collect

****************************************************************************

      subroutine hbuild(v,mat)
      implicit none

      integer, parameter :: lwork=136
      integer :: i,j,jj,info

      real(8) :: v(dim*2), mat(dim,dim), wn(dim), dum, x
      real(8) :: work(lwork),vr(dim,dim),vl(dim,dim),test(sites,sites)

C  Setting up noninteracting Hubbard-like system Hamiltonian with t as the
C  hopping constant between lattice sites.
      do i=1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = v(i) + v((dim-1)*2 + i)
          elseif (abs(i-j).eq.1) Then
            test(i,j) = -t
          else
            test(i,j) = 0.d0
          endif
        end do
      end do

      mat = 0.d0

      do i=1,sites
        do j=1,sites
          mat(i,j) = test(i,j)
        end do
      end do

      do i=1,sites
        do j=1,sites
          if (i.eq.j) Then
            test(i,j) = v(i) - v((dim-1)*2 + i)
          elseif (abs(i-j).eq.1) Then
            test(i,j) = -t
          else
            test(i,j) = 0.d0
          endif
        end do
      end do

      do i=1,sites
        do j=1,sites
          mat(dim/2+i,dim/2+j) = test(i,j)
        end do
      end do

      do i=1,sites
        mat(i,dim/2+i) = v(sites+i)
        mat(dim/2+i,i) = v(sites+i)
      end do

C  Eigenvalue solver for a complex, non-symmetric matrix.

!      call CGEEV('n','v',dim, mat, dim, En, VL, dim, VR, dim,
!     &                             WORK, LWORK, RWORK, INFO )

!      write(*,matrix) mat
!      write(*,*)'***********************'

      write(*,*) '****************************'
      write(*,matrix) mat

      call dgeev('n','v',dim,mat,dim,En,wn,VL,dim,VR,dim
     &                                   ,WORK,LWORK,INFO)

*********************************************************************
***     Sort the eigenvalues in ascending order
*********************************************************************
      DO 51 I=1,dim
        DO 52 J=i+1,dim
            IF (En(I).GE.En(J)) THEN
               DUM = En(I)
               En(I) = En(J)
               En(J) = DUM

               DUM = Wn(I)
               Wn(I) = Wn(J)
               Wn(J) = DUM

               DO 53 JJ=1,dim
                  DUM = VR(JJ,I)
                  VR(JJ,I) = VR(JJ,J)
                  VR(JJ,J) = DUM

                  DUM = VL(JJ,I)
                  VL(JJ,I) = VL(JJ,J)
                  VL(JJ,J) = DUM
53             CONTINUE

            ENDIF
52       CONTINUE
51    CONTINUE

      mat = vr

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

***************************************************************************

      subroutine Pauli(sigma)
      implicit none

      integer :: i,j,k,counter
      real(8), intent(out) :: sigma(spin*4,spin*4)

      sigma = 0.d0
      do i=1,spin
        do j=1,spin
          if (i.eq.j) then
            sigma(i,j) = 1.d0
            sigma(i+3*spin,j+3*spin) = (-1.d0)**(i-1)
          else
            sigma(i+spin,j+spin) = 1.d0
            !sigma(i+2*spin,j+2*spin) = 0.d0
          end if
        end do
      end do

      end subroutine Pauli

****************************************************************************


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
        fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
        do 12 j=1,n
          xi(j)=xmin*xi(j)
          p(j)=p(j)+xi(j)
12    continue
        return
        END subroutine
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

*************************************************************************
      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      REAL(8) dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=100,ZEPS=1.0d-10)
      INTEGER iter
      REAL(8) a,b,d,d1,d2,du,dv,dw,dx,e
      REAL(8) fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2       if(abs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.


*************************************************************************
      FUNCTION df1dim(x)
      INTEGER NMAX
      REAL(8) df1dim,x
      PARAMETER (NMAX=50)
CU    USES dfunc
      INTEGER j,ncom
      REAL(8) df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(xt,df)
      df1dim=0.d0
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
      END
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
