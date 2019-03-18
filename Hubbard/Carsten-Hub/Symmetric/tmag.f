      Program RMSD
      implicit none

      integer, parameter :: sites=4
      integer :: i, j, k, d
      real(8) :: tx(sites + 1), ty(sites + 1),
     &            tz(sites + 1), tt(sites + 1)
      character(150) :: bvector, FileList
      character(2) :: dir(4)
      character(150) :: path

      write(bvector, '(a, i3, a)') '(', sites + 1, 'e16.6)'

!      dir(1) = 'mx'
!      dir(2) = 'my'
!      dir(3) = 'mz'
!      dir(4) = 'n'

      dir(1) = 'tx'
      dir(2) = 'ty'
      dir(3) = 'tz'

      path =
     &'/home/ep3/Documents/Fortran/Physics/Hubbard/4pt-Lattice-NonCol/Te
     &mp/'

      do i=1,9

        write(FileList,'(a,a,i1,a)') 
     &                    trim(path),'4pt-B',i,'-ExactXC-tx.txt'
        FileList=trim(filelist)
        open(100, file=FileList)
  
        write(FileList,'(a,a,i1,a)')
     &                    trim(path),'4pt-B',i,'-ExactXC-ty.txt'
        FileList=trim(filelist)
        open(101, file=filelist)
  
        write(FileList,'(a,a,i1,a)')
     &                    trim(path),'4pt-B',i,'-ExactXC-tz.txt'
        FileList=trim(filelist)
        open(102, file=filelist)
  
        write(FileList,'(a,i1,a)') 'Exact-Torq-SymT',i,'.txt'
        FileList=trim(filelist)
        open(103, file=filelist)
  
        do j=1, 100
          read(100, bvector) tx
          read(101, bvector) ty
          read(102, bvector) tz
  
          tt = 0.d0
          tt(1) = dble(j)/10.d0
          do k=2, sites + 1
            tt(k) = dsqrt(tx(k)**2 + ty(k)**2 + tz(k)**2)
          end do
  
          write(103, bvector) tt
        end do

      close(100)
      close(101)
      close(102)
      close(103)

      end do  


      end program
