      Program RMSD
      implicit none

      integer, parameter :: sites=4
      integer :: i, j, k, d
      real(8) :: Mex(sites + 1), MAp(sites + 1),
     &            MDif(sites + 1)
      character(30) :: bvector, FileList
      character(2) :: dir(4) 

      write(bvector, '(a, i3, a)') '(', sites + 1, 'e16.6)'

      dir(1) = 'mx'
      dir(2) = 'my'
      dir(3) = 'mz'
      dir(4) = 'n'

      do d = 1, 4

        do i = 1, 9

          write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-BTot-',
     &                                           trim(dir(d)), '.txt'
          open(100, file = '/home/ep3/Documents/Fortran/Physics/Hubbard/
     &4pt-Lattice-NonCol/NonCol-MFiles/'FileList)
  
          write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-Slater-',
     &                                           trim(dir(d)), '.txt'
          open(101, file = FileList)
  
          write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-SvsExComp-',
     &                                           trim(dir(d)), '.txt'
          open(102, file = FileList)
  
          do j=1, 100
            read(100, bvector) MEx
  
            read(101, bvector) MAp
  
            BDiff = 0.d0
            BDiff(1) = dble(j)/10.d0
            do k=2, sites + 1
              MDif(k) = MEx(k) - MAp(k)
            end do
  
            write(102, bvector) BDiff
          end do
  
          close(100)
          close(101)
          close(102)

        end do

      end do

      end program
