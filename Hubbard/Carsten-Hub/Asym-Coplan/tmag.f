      Program RMSD
      implicit none

      integer, parameter :: sites=4
      integer :: i, j, k, d
      real(8) :: tx(sites + 1), ty(sites + 1),
     &            tz(sites + 1), tt(sites + 1)
      real(8) :: txa(sites + 2), tya(sites + 2), tza(sites + 2),
     &           tta(sites + 1)
      character(150) :: bvector, FileList, avector
      character(2) :: dir(4)
      character(150) :: path

      write(bvector, '(a, i3, a)') '(', sites + 1, 'e16.6)'
      write(avector, '(a, i3, a)') '(', sites + 2, 'e16.6)'

!      dir(1) = 'mx'
!      dir(2) = 'my'
!      dir(3) = 'mz'
!      dir(4) = 'n'

      dir(1) = 'tx'
      dir(2) = 'ty'
      dir(3) = 'tz'

!      path =
!     &'/home/ep3/Documents/Fortran/Physics/Hubbard/4pt-Lattice-NonCol/Te
!     &mp/'

      do i=1,1

        write(FileList,'(a,i1,a)')
     &                    '4pt-B',i,'-XC-Asym-tx.txt'
        FileList=trim(filelist)
        open(100, file=FileList)

        write(FileList,'(a,i1,a)')
     &                    '4pt-B',i,'-XC-Asym-ty.txt'
        FileList=trim(filelist)
        open(101, file=filelist)

        write(FileList,'(a,i1,a)')
     &                    '4pt-B',i,'-XC-Asym-tz.txt'
        FileList=trim(filelist)
        open(102, file=filelist)

        write(FileList,'(a,i1,a)') 'Asym-XC-TMag.txt'
        FileList=trim(filelist)
        open(103, file=filelist)

        do j=1, 201
          read(100, bvector) tx
          read(101, bvector) ty
          read(102, bvector) tz

!          read(100, avector) txa
!          read(101, avector) tya
!          read(102, avector) tza

          tt = 0.d0
          tt(1) = dble(j-1)/10.d0
          do k=2, sites + 1
            tt(k) = dsqrt(tx(k)**2 + ty(k)**2 + tz(k)**2)
!            tt(k) = dsqrt(txa(k+1)**2 + tya(k+1)**2 + tza(k+1)**2)
          end do

          write(103, bvector) tt
        end do

      close(100)
      close(101)
      close(102)
      close(103)

      end do


      end program
