      program Schro_Lattice
      
      double precision E1, E2, V2, DELTA, M(2,2), STEP1, STEP2
      
      Write(*,*) 'Give value for Delta:'
      read(*,*) DELTA
      
      Write(*,*) 'Give potential at Lattice Point 2' 
      read(*,*) V2
      
      step1 = 1./delta
      step2 = -1./(2.*delta)
      
      M(1,1) = step1
      M(1,2) = step2
      M(2,1) = step2
      M(2,2) = step1+V2
      
      write(*,*) 'The Hamiltonian Matrix is:', M
      
      E1 = ((step1 + V2) + sqrt(v2**2 + 4.*(step2)**2))/2.
      E2 = ((step1 + V2) - sqrt(v2**2 + 4.*(step2)**2))/2.
      
      write(*,*) 'Calculated Energy:', E1, E2
      
      end    
      