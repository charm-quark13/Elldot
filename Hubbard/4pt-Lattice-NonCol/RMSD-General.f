      Program RMSD
      implicit none

      integer, parameter :: sites=4
      integer :: i, j, k, d
      real(8) :: BLong(sites + 1), BTot(sites + 1),
     &            BDiff(sites + 1)
      character(30) :: bvector, FileList
      character(2) :: dir(4) 

      write(bvector, '(a, i3, a)') '(', sites + 1, 'e16.6)'

!      dir(1) = 'mx'
!      dir(2) = 'my'
!      dir(3) = 'mz'
!      dir(4) = 'n'

      dir(1) = 'tx'
      dir(2) = 'ty'
      dir(3) = 'tz' 

      do d = 1, 3

        do i = 1, 11

          if (i.lt.10) then
            write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-Exact-',
     &                                           trim(dir(d)), '.txt'
            open(100, file = FileList)
  
            write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-Slater-',
     &                                           trim(dir(d)), '.txt'
            open(101, file = FileList)
  
            write(FileList,'(a,i1,a,a,a)') '4pt-B', i, '-TDiff-',
     &                                           trim(dir(d)), '.txt'
            open(102, file = FileList)
          else
            write(FileList,'(a,i2,a,a,a)') '4pt-B', i, '-Exact-',
     &                                           trim(dir(d)), '.txt'
            open(100, file = FileList)

            write(FileList,'(a,i2,a,a,a)') '4pt-B', i, '-Slater-',
     &                                           trim(dir(d)), '.txt'
            open(101, file = FileList)

            write(FileList,'(a,i2,a,a,a)') '4pt-B', i, '-TDiff-',
     &                                           trim(dir(d)), '.txt'
            open(102, file = FileList)
          end if            


          do j=1, 100
            read(100, bvector) BTot
  
            read(101, bvector) BLong
  
            BDiff = 0.d0
            BDiff(1) = dble(j)/10.d0
            do k=2, sites + 1
              BDiff(k) = BTot(k) - BLong(k)
!              BDiff(k) = BDiff(k)/Btot(k) 
            end do
  
            write(102, bvector) BDiff
          end do
  
          close(100)
          close(101)
          close(102)

        end do

      end do

      end program
