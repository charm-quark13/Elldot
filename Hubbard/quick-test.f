      Program WhileTest
      implicit none

      integer :: i, j
      real(8) :: v(4), vprev(4), vstart(4), y

      y = 2.d0
      do i = 1, 4
        vstart(i) = dble(i)
      end do

      write(*,*) 'V Start'
      write(*,*) vstart
      vprev = vstart
      v = vstart
      do i = 1, 4 - 1
        do j = i+1, 4
          v = vprev
          v(i) = v(i) + 1
          v(j) = v(j) + 1
          write(*,*) 'V'
          write(*,*) v
        end do
        vprev(i) = vstart(i) + 1
        write(*,*) 'V Previous'
        write(*,*) vprev
      end do
      end
