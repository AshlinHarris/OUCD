MODULE standard

!---------------------------------------------------------------------------------------
! standard module for double precision
!---------------------------------------------------------------------------------------

	INTEGER, parameter:: dp=kind(0.d0)

!---------------------------------------------------------------------------------------
! convert a coordinate to a node index
!---------------------------------------------------------------------------------------

	CONTAINS

	SUBROUTINE COORD_TO_INDEX(i,n,x,y)

		IMPLICIT NONE

		INTEGER, INTENT(OUT) :: i ! index
		INTEGER, INTENT(IN) :: n ! columns in the array
		INTEGER, INTENT(IN) :: x,y ! coordinate

		i = (x-1)*n+y

	END SUBROUTINE COORD_TO_INDEX

!---------------------------------------------------------------------------------------
! Convert time in seconds to an array of hours, minutes, seconds, and milliseconds
!---------------------------------------------------------------------------------------

SUBROUTINE FORMAT_TIME(val,res)

	! USE :: STANDARD

	IMPLICIT NONE

	REAL(dp), INTENT(IN) :: val
	INTEGER, INTENT(OUT) :: res(4)
	INTEGER :: hours, minutes, seconds, milliseconds

	hours = floor(val / 3600)
	minutes = mod(floor(val / 60), 60)
	seconds = mod(floor(val),60)
	milliseconds = floor((val-floor(val))*1000.0_dp)

	res = (/ hours, minutes, seconds, milliseconds /)

END SUBROUTINE FORMAT_TIME

END MODULE standard

