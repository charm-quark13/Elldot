      program Psi0_t0
      
      IMPLICIT NONE
      
      double precision E1, V2, DInv, N, SQC
      double precision psi0(2), a, b, total
      parameter(v2=-15)
      
      open(unit=11,file='Psi0_t0.txt',action='write',status='replace')
      
      DInv = 0.1D0
      
      SQC = dsqrt(v2**2 + DInv**2)
      
      a = (v2+SQC)
      b = (v2/dinv**2)
      
      N= 1/(2.*(1+b*a))
      
      Psi0(1)=(1*dsqrt(N))
      
      Psi0(2)=(1./Dinv)*(V2 + sqc)*dsqrt(N)
      
      Total= psi0(1)**2 + psi0(2)**2
      
      Write(11,*) PSI0(1), PSI0(2)
      
      close(11)
      
      write(*,*) total, psi0(1)**2, psi0(2)**2
      
      end