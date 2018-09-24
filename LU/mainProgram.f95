program mainProgram
  implicit none

  integer :: n
  real, dimension(:,:), allocatable :: a
  real, dimension(:), allocatable :: x

  n = 5   !!

  allocate(a(n, n+1))
  allocate(x(n))

  a(1,1) = 15
  a(1,2) = -1
  a(1,3) = 2
  a(1,4) = -3
  a(1,5) = 4
  a(1,6) = 8

  a(2,1) = 2
  a(2,2) = 23
  a(2,3) = -1
  a(2,4) = 5
  a(2,5) = -2
  a(2,6) = 82.4

  a(3,1) = -1
  a(3,2) = 3
  a(3,3) = 92
  a(3,4) = -5
  a(3,5) = 1
  a(3,6) = -764.9

  a(4,1) = 1
  a(4,2) = 2
  a(4,3) = 1
  a(4,4) = 27
  a(4,5) = 3
  a(4,6) = -8.9

  a(5,1) = -4
  a(5,2) = -6
  a(5,3) = -2
  a(5,4) = 8
  a(5,5) = 41
  a(5,6) = -201.9



  call luDecomp(a, n, x)
  write(*,*) "Solution by LU Decomposition:"

  call printMatrix(x, n, 1)
end program mainProgram
