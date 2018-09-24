subroutine triDiagonal(aIn, n, x)
  implicit none

  real, dimension(n,n+1), intent(in) :: aIn
  integer, intent(in) :: n
  real, dimension(n), intent(inout) :: x

  real, dimension(n,n+1) :: a
  integer :: stepCount, rowCount, columnCount
  integer :: pivotRow, pivotColumn
  real :: factor, pivot

  a = aIn


  do stepCount = 1, (n-1)

    pivotRow = stepCount
    pivotColumn = stepCount
    pivot = a(pivotRow, pivotColumn)

    do rowCount = (pivotRow + 2), n
      factor =  a(rowCount, pivotColumn)/pivot
      do columnCount = pivotColumn, (n+1)
        a(rowCount,columnCount) = a(rowCount, columnCount) - factor*a(pivotRow,columnCount)
      end do
    end do
    call printMatrix(a, n ,(n+1))
  end do

  do stepCount = 1, (n-1)

    pivotRow = stepCount
    pivotColumn = stepCount
    pivot = a(pivotRow, pivotColumn)

    do rowCount = (pivotRow + 2), n
      factor =  a(pivotColumn, rowCount)/pivot
      do columnCount = pivotColumn, (n+1)
        a(columnCount,rowCount) = a(columnCount, rowCount) - factor*a(columnCount,pivotRow)
      end do
    end do
    call printMatrix(a, n ,(n+1))
  end do



end subroutine triDiagonal
