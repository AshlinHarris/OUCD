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

! --------------------------------------------------------------------
! given a tile, return a section of the infinite tiling
! --------------------------------------------------------------------

SUBROUTINE GENERATE_TILING(a,m,n,x,b)

	! To check all the adjacencies of a tile for unique complete
	! dominance, points outside the tile must be checked.

	! This function creates an m+2x by n+2x array and puts an
	! m by n tile in the center.

	! The remaining sides and corners of the larger array are
	! filled with slices of the tile so that a small area of the
	! tiling is simulated.

	! I think setting x=2 is sufficient to reliably simulate the
	! behavior of the infinite tiling. This gives a layer of nodes
	! to dominate the interior nodes an an additional layer of
	! interior nodes to check for matches

	! Without a border, the exterior nodes will not be assigned
	! the correct adjacency list.

	! At x=1, spurious solutions may be produced such as:
	! /----------\
	! |**********|
	! |          |
	! |          |
	! \----------/

	IMPLICIT NONE

	LOGICAL, INTENT(IN) :: a(:,:)

	INTEGER, INTENT(IN) :: m ! rows in the tile
	INTEGER, INTENT(IN) :: n ! columns in the tile
	INTEGER, INTENT(IN) :: x ! width of border;
	                         ! should not exceed m or n

	LOGICAL, INTENT(OUT) :: b(:,:) ! tiling

	b(    1:x    ,     1:x    ) = a(m-x+1:m, n-x+1:n)
	b(    1:x    ,   x+1:n+x  ) = a(m-x+1:m,      : )
	b(    1:x    , n+x+1:2*x+n) = a(m-x+1:m,     1:x)

	b(  x+1:x+m  ,     1:x    ) = a(     : , n-x+1:n)
	b(  x+1:x+m  ,   x+1:n+x  ) = a(     : ,      : )
	b(  x+1:x+m  , n+x+1:2*x+n) = a(     : ,     1:x)

	b(x+m+1:2*x+m,     1:x    ) = a(    1:x, n-x+1:n)
	b(x+m+1:2*x+m,   x+1:n+x  ) = a(    1:x,      : )
	b(x+m+1:2*x+m, n+x+1:2*x+n) = a(    1:x,     1:x)

END SUBROUTINE GENERATE_TILING


END MODULE standard

