subroutine forwardSubstitution(l, b, n, x)
  implicit none
  integer, intent(in) :: n

  real, dimension(n, n), intent(in) :: l
  real, dimension(n), intent(in) :: b
  real, dimension(n), intent(out) :: x

  integer :: rowCount
  real :: summation
  integer :: termCount

  x(1) = b(1) / l(1,1)

  do rowCount = 2, n
    summation = 0
    do termCount = 1, (rowCount - 1)
      summation = summation + (l(rowCount, termCount) * x(termCount))
    end do
    x(rowCount) = (b(rowCount) - summation) / l(rowCount, rowCount)
  end do

end subroutine forwardSubstitution
