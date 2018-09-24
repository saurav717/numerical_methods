logical function diagonallyDominant(matrix, rows, columns)
  implicit none

  integer, intent(in) :: rows, columns
  real, dimension(rows,columns), intent(in) :: matrix

  integer :: rowCount, columnCount

  real :: diagonalElement
  real :: sumOfNonDiagonalElements

  diagonallyDominant = .true.

  do rowCount = 1,rows
    diagonalElement = abs(matrix(rowCount,rowCount))
    sumOfNonDiagonalElements = 0
    do columnCount = 1, columns
      if(columnCount /= rowCount) then
        sumOfNonDiagonalElements = sumOfNonDiagonalElements + abs(matrix(rowCount, columnCount))
      end if
    end do

    if(sumOfNonDiagonalElements > diagonalElement) then
      diagonallyDominant = .false.
      return
    end if
  end do
end function diagonallyDominant
