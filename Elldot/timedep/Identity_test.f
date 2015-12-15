      program identity
      
      integer j, k
      
      real Ide(2,2)
      
      Do 10 j=1,2
          Do 11, k=1,2
              If(j .eq. k) Ide(j,k)=1
              if(j .ne. k) Ide(j,k)=0
11    continue
10    continue
      
      write(*,*)'2,2',Ide(2,2)
      write(*,*)'2,1',Ide(2,1)
      write(*,*)'1,1',Ide(1,1)
      write(*,*)'1,2',Ide(1,2)
      
      end
          