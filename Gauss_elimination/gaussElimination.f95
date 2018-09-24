subroutine gaussElimination(aIn, n, x)
  implicit none

  integer, intent(in) :: n
  real, dimension(n,n+1), intent(in) :: aIn
  real, dimension(n), intent(out) :: x

  real, dimension(n,n+1) :: a

  integer :: stepCount, rowCount, columnCount

  integer :: pivotRow, pivotColumn
  real :: pivot
  real :: factor

  write(*,*)
  write(*,*) "system of linear algebraic equations"

  a = aIn
  write(*,*) "the augmented matrix"
  call printMatrix(a,n,(n+1))

  write(*,*) "part -1 : forward elimination"

  do stepCount = 1, (n-1)

    pivotRow = stepCount
    pivotColumn = stepCount
    pivot = a(pivotRow, pivotColumn)

    do rowCount = (pivotRow + 1), n
      factor =  a(rowCount, pivotColumn)/pivot
      do columnCount = pivotColumn, (n+1)
        a(rowCount,columnCount) = a(rowCount, columnCount) - factor*a(pivotRow,columnCount)
      end do
    end do
    call printMatrix(a, n ,(n+1))
  end do

  write(*,*) "hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee"

  write(*,50) "evaluate x(", n,")"
  x(n) = a(n,n+1)/a(n,n)
  write(*,40) "x(n) = ", x(n)
  write(*,*)

  do rowCount = (n - 1), 1, -1
  !!  write(*,50) "evaluate x(",rowCount,")"
    factor = 0
    do stepCount = (rowCount + 1), n
      factor = factor + a(rowCount, stepCount) * x(stepCount)
    end do
    x(rowCount) = (1/a(rowCount,rowCount)) * (a(rowCount, n+1) - factor)
    write(*,60) "x(",rowCount,") = ", x(rowCount)
    write(*,*)
  end do

  40 format(a6, f5.2)
  50 format(a11, i1, a1)
  60 format(a3, i1, a4, f5.2)

end subroutine gaussElimination
