subroutine gaussSiedel(aIn, n, xGuess, tolerance, x)
  implicit none

  logical, external :: diagonallyDominant

  integer, intent(in) :: n
  real, dimension(n, (n+1)), intent(in) :: aIn
  real, dimension(n), intent(in) :: xGuess
  real, intent(in) :: tolerance

  real, dimension(n), intent(out) :: x

  real,dimension(n,(n+1)) :: a
  real, dimension(n,n) :: coefficientMatrix
  real, dimension(n) :: xPrevious
  real, dimension(n) :: error

  real :: factor
  logical :: solutionFound = .false.

  integer :: iteration = 0
  integer :: rowCount, columnCount

  write(*,*)
  write(*,*) " system of linear algebraic equations"
  write(*,*) "method - jacobi iterative"
  write(*,*)

  a = aIn
  write(*,*) "augmented matrix"
  call printMatrix(a,n,(n+1))

  do rowCount = 1,n
    do columnCount = 1,n
      coefficientMatrix(rowCount,columnCount) =  a(rowCount, columnCount)
    end do
  end do

  if(diagonallyDominant(coefficientMatrix,n,n) .eqv. .false.) then
    stop "Error (Jacobi method  ): The matrix is not diagonally dominant"
  end if

  x = xGuess
  xPrevious = x

  do while(solutionFound .eqv. .false.)
    iteration = iteration + 1

    do rowCount = 1,n
      factor = 0
      do columnCount = 1,n
        if(columnCount /= rowCount) then
          factor = factor + a(rowCount, columnCount) * x(columnCount)
        end if
      end do
      x(rowCount) = (1/a(rowCount,rowCount)) * (a(rowCount, (n+1)) -factor)
    end do

    do rowCount = 1,n
        error(rowCount) = abs((x(rowCount)  -  xPrevious(rowCount))/(x(rowCount)))
    end do

    write(*,10) "iteration #", iteration
    write(*,*) "estimated solution : "
    call printMatrix(x,n,1)

    write(*,*) "Error :"
    call printMatrix(error, n, 1)

    if(maxval(error) < tolerance) then
        solutionFound = .true.
    end if

    xPrevious = x
  end do

  10 format(a11, i2)

end subroutine gaussSiedel
