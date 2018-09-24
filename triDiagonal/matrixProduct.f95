subroutine matrixProduct(a, b, rowsA, columnsA, rowsB, columnsB, ab)
  implicit none
  integer, intent(in) :: rowsA, columnsA, rowsB, columnsB

  real, dimension(rowsA, columnsA), intent(in) :: a
  real, dimension(rowsB, columnsB), intent(in) :: b
  real, dimension(rowsA, columnsB), intent(out) :: ab

  integer :: rowCount, columnCount
  real :: summation
  integer :: termCount

  if(columnsA /= rowsB) then
    stop  "Error -- Matrixes cannot be multiplied."
  end if

  do rowCount = 1, rowsA
    do columnCount = 1, columnsB
      summation = 0
      do termCount = 1, columnsA
        summation = summation + ((a(rowCount, termCount)) * (b(termCount, columnCount)))
      end do
      ab(rowCount, columnCount) = summation
    end do
  end do

end subroutine matrixProduct
