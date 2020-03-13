      Program mdiff

      Implicit none

      integer, parameter :: sites = 4
      integer :: i, j
      real(8) :: diffv(sites + 1), mlong(sites + 1)
     &           , mex(sites + 1)

      character(30) :: mlread, mxread

      write(mlread, '(a,i3,a)') '(', sites + 1, 'e16.6)'
!      write(mxread, '(a,i3,a)') '(', sites + 2, 'f14.8)'

      open(1, file='4pt-B1-SymCP-mxlong-t.txt')
      open(11, file='4pt-B1-SymCP-mylong-t.txt')
      open(111, file='4pt-B1-SymCP-mzlong-t.txt')
      open(1111, file='4pt-B1-SymCP-nlong-t.txt')

      open(2, file='4pt-B1-Asym-mx.txt')
      open(22, file='4pt-B1-Asym-my.txt')
      open(222, file='4pt-B1-Asym-mz.txt')
      open(2222, file='4pt-B1-Asym-n.txt')

      open(3, file='asym-mx-diff-t.txt')
      open(33, file='asym-my-diff-t.txt')
      open(333, file='asym-mz-diff-t.txt')
      open(3333, file='asym-n-diff-t.txt')

      do i = 1, 101
        read(1, mlread) mlong
        read(2, mlread) mex

        diffv(1) = mlong(1)

        do j=2,sites + 1
          diffv(j) = mex(j) - mlong(j)
        end do
        write(3, mlread) diffv

        read(11, mlread) mlong
        read(22, mlread) mex

        diffv(1) = mlong(1)

        do j=2,sites + 1
          diffv(j) = mex(j) - mlong(j)
        end do
        write(33, mlread) diffv

        read(111, mlread) mlong
        read(222, mlread) mex

        diffv(1) = mlong(1)

        do j=2,sites + 1
          diffv(j) = mex(j) - mlong(j)
        end do
        write(333, mlread) diffv

        read(1111, mlread) mlong
        read(2222, mlread) mex

        diffv(1) = mlong(1)

        do j=2,sites + 1
          diffv(j) = mex(j) - mlong(j)
        end do
        write(3333, mlread) diffv

      end do

      end
