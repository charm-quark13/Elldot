      Program Time_Schro
      
      IMPLICIT NONE
      
      double precision E1, E2, j, diff
      double precision DELTA, TOTAL
      double precision Prob, N, Norm1, Norm2
      double precision a, b, c
      integer V2
            
      Write(*,*) 'Give value for Delta:'
      read(*,*) DELTA
      
      Open(unit=20,file="Schro1.txt",action="write",status="replace")
      
      v2 = 100
      
          j = 1.D0/delta

c     Simplification to make E's easier to read and to follow:      

          c = v2**2 + j**2

          E1= 0.5D0*(2*j - V2 + dsqrt(c))
          E2= 0.5D0*(2*j - V2 - dsqrt(c))
      
           
c     Simplificaiton to make N easier to read and follow:
      
          a = (V2 + dsqrt(c))
          b = (V2/j**2)

c     normalization constant       
      
          N= 1/(2.*(1+b*a))      

c     Value of n(1) for Psi_1 
          norm1= (1.*dsqrt(N))**2

c     Value of n(2) for Psi_1 
          norm2= ((1./j)*(V2 + dsqrt(c))*dsqrt(N))**2
      
      
            
      Write(*,*) V2, norm1, norm2
      
      diff = E1 - E2
      write(*,*) Diff          
            
      End