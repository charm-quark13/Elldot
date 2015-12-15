      Program TD_SchrodT
      
      IMPLICIT NONE
      
      Double complex i, H(2,2), Psi(2), D, C,  R(2), A(4)
      
      Double precision T, DeltaT, Ide(2,2), V2
      
      Integer j, k, l, m
      
      parameter(i=(0,1))
      
c     opening read in files      
      open(unit=99,file='Psi0_t0.txt',action='read',status='old')
c     psi gs t=0 values (V=15 & delta=10)      
      psi(1) = 0.99999444461345821
      psi(2) = 3.3332777793650408E-003
      open(unit=98,file='Psi_dt.txt',action='write',status='replace')
      open(unit=97,file='Psi_dt_im.txt',action='write',status='replace')
      D=(0.1D0,0)
      V2=-15.0D0
      
      DeltaT = 0.01D0
      T = 0.0D0
      
c     iteration for identity matrix      
      Do 10 j=1,2
          Do 11, k=1,2
              If(j .eq. k) Ide(j,k)=1
              if(j .ne. k) Ide(j,k)=0
11        continue
10    continue
                
      do 20 l=1,1000
          T=T+DeltaT
      C = (i*T)/2.0D0
      write(*,*)
c     values for hamiltonian matrix
      H(1,1) = 0
      H(1,2) = -0.5*D*C
      H(2,1) = -0.5*D*C
      H(2,2) = (D-V2)*C
c      write(*,*) H 

      A(1)=Ide(1,1)+H(1,1)
      A(2)=Ide(1,2)+H(1,2)
      A(3)=Ide(2,1)+H(2,1)
      A(4)=Ide(2,2)+H(2,2)

c     iteration to solve for the right side of equaiton      
      Do 30 m=1,2
      R(m)=(Ide(m,1)-H(m,1))*Psi(1)+(Ide(m,2)-H(m,2))*psi(2)
30    continue
      
c     iteration to solve for the individual psi function
      
      psi(2)=(A(1)*R(2)-A(3)*R(1))/(A(2)*A(3)+A(4)*A(1))
      
      psi(1)=(R(1)-A(2)*psi(2))/A(1)
      
c     printing to separate txt files, complex and real 
      
      Write(98,*) Dreal(psi(1)), dreal(psi(2))
      write(97,*) dimag(psi(1)), dimag(psi(2))
      
20    continue
      
      close(97)
      close(98)
      
      Write(*,*) 'file complete'
      
      end