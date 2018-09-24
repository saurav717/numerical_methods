subroutine printMatrix(matrix, rows, columns)
  implicit none

  integer, intent(in) :: rows, columns
  real, dimension(rows,columns), intent(in) :: matrix

  integer :: rowCount, columnCount

  do rowCount = 1, rows
    do columnCount = 1, columns
      write(*,10,advance = 'no') matrix(rowcount, columnCount)
    end do
    write(*,*)
  end do
  write(*,*)

  10 format(f7.2)
end subroutine printMatrix
