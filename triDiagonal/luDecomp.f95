subroutine luDecomp(aIn, n, x)
  implicit none
  integer, intent(in) :: n
  real, dimension(n, n+1), intent(inout) :: aIn
  real, dimension(n), intent(out) :: x
  real, dimension(n, n) :: a
  real, dimension(n) :: b


  real, dimension(n, n) :: l
  real, dimension(n, n) :: u
  real, dimension(n, n) :: productOfLU
  real, dimension(n) :: y
  real, dimension(n,n+1) :: aa

  integer :: rowCount, stepCount, columnCount, pivotRow, pivotColumn
  real :: factor, pivot
  real :: summation
  integer :: termCount


  write(*,*)
  write(*,*) "System of Linear Algebraic Equations"
  write(*,*) "Method: LU Decomposition "
  write(*,*)


    aa = aIn

    do stepCount = 1, (n-1)

      pivotRow = stepCount
      pivotColumn = stepCount
      pivot = aa(pivotRow, pivotColumn)

      do rowCount = (pivotRow + 2), n
        factor =  aa(rowCount, pivotColumn)/pivot
        do columnCount = pivotColumn, (n+1)
          aa(rowCount,columnCount) = aa(rowCount, columnCount) - factor*aa(pivotRow,columnCount)
        end do
      end do
      call printMatrix(aa, n ,(n+1))
    end do

    do stepCount = 1, (n-1)

      pivotRow = stepCount
      pivotColumn = stepCount
      pivot = aa(pivotRow, pivotColumn)

      do rowCount = (pivotRow + 2), n
        factor =  aa(pivotColumn, rowCount)/pivot
        do columnCount = pivotColumn, (n+1)
          aa(columnCount,rowCount) = aa(columnCount, rowCount) - factor*aa(columnCount,pivotRow)
        end do
      end do
      call printMatrix(aa, n ,(n+1))
    end do

    aIn = aa


  do rowCount = 1, n
      do columnCount = 1, n
        a(rowCount, columnCount) = aIn(rowCount, columnCount)
      end do
      b(rowCount) = aIn(rowCount, (n + 1))
  end do

  l = 0
  u = 0

  write(*,*)  "The Coefficient Matrix "
  call printMatrix(a, n, n)

  write(*,*)  "The RHS"
  call printMatrix(b, n, 1)

  write(*,*)  "The Lower Triangular Matrix "
  call printMatrix(l, n, n)

  write(*,*)  "The Upper Triangular Matrix"
  call printMatrix(u, n, n)

  write(*,*)  "Part 1 of 3: Finding L and U "
  write(*,*)  "Part 1A: Set u(i,i) = 1 (diagonal elements of u)"

  do rowCount = 1, n
    u(rowCount, rowCount) = 1.0
  end do

  write(*,*)  "The Lower Triangular Matrix "
  call printMatrix(l, n, n)

  write(*,*)  "The Upper Triangular Matrix "
  call printMatrix(u, n, n)

  write(*,*)  " Setring l(i,1) = a(i,1) "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do rowCount = 1, n
    l(rowCount, 1) = a(rowCount, 1)
  end do

  write(*,*) "The Lower Triangular Matrix"
  call printMatrix(l, n, n)

  write(*,*) "The Upper Triangular Matrix"
  call printMatrix(u, n, n)

  write(*,*)  " Setting u(1,j) = a(1,j) / l(1,1) "

  do columnCount = 2, n
    u(1, columnCount) = a(1, columnCount) / l(1, 1)
  end do

  write(*,*) "The Lower Triangular Matrix"
  call printMatrix(l, n, n)

  write(*,*) "The Upper Triangular Matrix"
  call printMatrix(u, n, n)

  write(*,*)  " Evaluate the remaining u(i,j) and l(i,j) "
  do  rowCount = 2, n
    do columnCount = 2, rowCount
      summation = 0
      do termCount = 1, (columnCount - 1)
        summation = summation + (l(rowCount, termCount) * u(termCount,columnCount))
      end do
      l(rowCount, columnCount) = a(rowCount, columnCount) - summation
    end do

    do columnCount = (rowCount + 1), n
        summation = 0
        do termCount = 1, (rowCount - 1)
          summation = summation + (l(rowCount, termCount) * u(termCount, columnCount))
        end do
        u(rowCount, columnCount) = (a(rowCount, columnCount) - summation) / l(rowCount, rowCount)
    end do
  end do

  write(*,*) "The Lower Triangular Matrix (Final)"
  call printMatrix(l, n, n)

  write(*,*) "The Upper Triangular Matrix (Final)"
  call printMatrix(u, n, n)

  call matrixProduct(l, u, n, n, n, n, productOfLU)

  write(*,*) "The Coefficient Matrix"
  call printMatrix(a, n, n)

  write(*,*) "The Product of L and U"
  call printMatrix(productOfLU, n, n)

  write(*,*) " Using forward substitution to solve ly = b"
  call forwardSubstitution(l, b, n, y)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)  "Intermediate Column Vector y"
  call printMatrix(y, n, 1)

  write(*,*)  " Using backward substitution to solve ux = y "
  call backSubstitution(u, y, n, x)

  write(*,*)  "Solution x "
  call printMatrix(x, n, 1)

end subroutine luDecomp
