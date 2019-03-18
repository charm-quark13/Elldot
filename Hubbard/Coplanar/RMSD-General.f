      Program RMSD
      implicit none

      integer, parameter :: sites=4
      integer :: i, j, k, d
      real(8) :: tx(sites + 1), ty(sites + 1),
     &            tz(sites + 1), tt(sites + 1)
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

      write(FileList,'(a)') '4pt-B1-Ext-CP-tx.txt'
      open(100, file=FileList)

      write(FileList,'(a)') '4pt-B1-Ext-CP-ty.txt'
      open(101, file=filelist)

      write(FileList,'(a)') '4pt-B1-Ext-CP-tz.txt'
      open(102, file=filelist)

      write(FileList,'(a)') 'Ext-torq-total.txt'
      open(103, file=filelist)

      do j=0, 100
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

      end program
