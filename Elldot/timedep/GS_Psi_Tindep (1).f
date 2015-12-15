      program Psi0_t0
      
      IMPLICIT NONE
      
      double precision E, V, DInv, N, c
      double precision psi0(2), a, b, total
      DOUBLE COMPLEX PSI(2) 
      parameter(v=.1D0)
      
      DInv = 0.05D0
c     Integrand for Binomial to find E_gs:      
      C = v**2 + 4.D0*DInv**2
      
c     Ground state energy:      
      E= 0.5D0*(V - dsqrt(c))

c     defining normalization constant:
      N= 1.D0/dsqrt(((Dinv**2+(V-E)**2)/DInv**2))
      
      Psi0(1)=1*N

c     Term after N is wavefunction:      
      Psi0(2)=((V-E)/Dinv)*N
      
      Total= psi0(1)**2 + psi0(2)**2
      
      Write(1,*) PSI0(1), PSI0(2)
      
      close(1)
      
      write(*,*) total, psi0(1)**2, psi0(2)**2
      
      PSI(1) = PSI0(1)*(1.D0,0.D0)
      PSI(2) = PSI0(2)*(1.D0,0.D0)
      
      WRITE(1) PSI

      end
