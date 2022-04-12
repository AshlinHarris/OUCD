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

END MODULE standard

