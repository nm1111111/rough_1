SUBROUTINE SlidingWindow(x, n, shift, m, output, fun)
IMPLICIT NONE

INTEGER      :: n, i, m 
REAL(8)      :: x(n), output(n), shift(m)
INTEGER      :: ishift(m), startPosition(n), endPosition(n)
CHARACTER(3) :: fun

do i = 1, n
  
  endPosition(i) = i

end do

if(m == 1) then

  ishift(1) = idnint(shift(1))
  ishift(1) = abs(ishift(1))

  startPosition = endPosition - ishift(1)
  
else 

  ishift = idnint(shift)
  ishift = abs(ishift)

  startPosition = endPosition - ishift

end if

where(startPosition < 1)

  startPosition = 1

end where

if(fun .eq. 'DLY') then

  output = x(startPosition)

end if

if(fun .eq. 'SUM' .or. fun .eq. 'SMA') then

  do i = 1, n

    output(i) = sum(x(startPosition(i):endPosition(i)))

  end do

end if

if(fun .eq. 'SMA') then

  output = output / (endPosition - startPosition + 1)

end if

END SUBROUTINE SlidingWindow


SUBROUTINE Reservoir(x, n, a, output)
IMPLICIT NONE

INTEGER :: n
REAL(8) :: x(n), output(n), a(n)
REAL(8) :: storage, sumOfStorageAndInflow
INTEGER :: i

a = abs(a)
a = 1_8 / a
storage = x(1)

do i = 1, n

  if(IsNaN(storage) .or. storage < 0) then

    storage = 0
  
  endif

  sumOfStorageAndInflow = storage + x(i)

  output(i) = a(i) * sumOfStorageAndInflow
  
  storage = sumOfStorageAndInflow - output(i)

end do

END SUBROUTINE Reservoir
