subroutine backSubstitution(u, b, n, x)
  implicit none

  integer, intent(in) :: n

  real, dimension(n, n), intent(in) :: u
  real, dimension(n), intent(in) :: b
  real, dimension(n), intent(out) :: x

  integer :: rowCount
  real :: summation
  integer :: termCount

  x(n) = b(n) / u(n,n)
  do rowCount = (n - 1), 1, -1
    summation = 0
    do termCount = (rowCount + 1), n
      summation = summation + (u(rowCount, termCount) * x(termCount))
    end do
    x(rowCount) = (b(rowCount) - summation) / u(rowCount, rowCount)
  end do

end subroutine backSubstitution
