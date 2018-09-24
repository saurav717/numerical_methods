program mainProgram
  implicit none

  integer :: n

  real, dimension(:,:), allocatable :: a
  real, dimension(:), allocatable :: x

  integer :: rowCount

  n = 5
  allocate(a(n,n+1))
  allocate(x(n))

  a(1,1) = 1
  a(1,2) = -1
  a(1,3) = 2
  a(1,4) = -3
  a(1,5) = 4
  a(1,6) = -35.4

  a(2,1) = 2
  a(2,2) = 3
  a(2,3) = -1
  a(2,4) = 5
  a(2,5) = -2
  a(2,6) = 32.4

  a(3,1) = -1
  a(3,2) = 3
  a(3,3) = 2
  a(3,4) = -5
  a(3,5) = 1
  a(3,6) = -17.9

  a(4,1) = 1
  a(4,2) = 2
  a(4,3) = 1
  a(4,4) = 2
  a(4,5) = 3
  a(4,6) = -13.9

  a(5,1) = -4
  a(5,2) = -6
  a(5,3) = -2
  a(5,4) = 8
  a(5,5) = -3
  a(5,6) = 4.9

  call gaussElimination(a, n, x)

  write(*,*) "solution"

  do rowCount = 1,n
    write(*,10) x(rowCount)
  end do

  10 format(f5.2)

end program mainProgram
