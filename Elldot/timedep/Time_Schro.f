      program Psi0_t0
      
      IMPLICIT NONE
      
      double precision E, V(100), DInv, N, c
      double precision psi0(2), a, b, total
      integer i
      DOUBLE COMPLEX PSI(2) 
      
      do 10 i=1,100
          v(i)=51-i
      
      DInv = 0.05D0
c     Integrand for Binomial to find E_gs:      
      C = v(i)**2 + 4.D0*DInv**2
      
c     Ground state energy:      
      E= 0.5D0*(V(i) - dsqrt(c))

c     defining normalization constant:      
      N= 1.D0/dsqrt(((Dinv**2+(V(i)-E)**2)/DInv**2))
      
      Psi0(1)=1*N

c     Term after N is wavefunction:      
      Psi0(2)=((V(i)-E)/Dinv)*N
      
      Total= psi0(1)**2 + psi0(2)**2
      
c      Write(*,*) i, PSI0(1), PSI0(2)
      
       write(*,*) i, 'Discriminant:',C, 'Normalization Constant:',
     & N, total
      
10    continue
      
      PSI(1) = PSI0(1)*(1.D0,0.D0)
      PSI(2) = PSI0(2)*(1.D0,0.D0)

      end