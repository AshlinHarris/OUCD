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

!---------------------------------------------------------------------------------------
! create the initial set in the chain of dominating node sets
!---------------------------------------------------------------------------------------

SUBROUTINE GET_FIRST_SET(d,n)

	IMPLICIT NONE
	
	INTEGER, INTENT(INOUT) :: d(:)
	INTEGER, INTENT(IN)	:: n

	INTEGER :: i

	DO i=1,n
		d(i) = i
	END DO

END SUBROUTINE GET_FIRST_SET

!---------------------------------------------------------------------------------------
! check if a set is the final in the chain
!---------------------------------------------------------------------------------------

SUBROUTINE GET_NEXT_SET(d,size,m,n,updated)
	
	IMPLICIT NONE

	INTEGER, INTENT(INOUT) :: d(:)
	INTEGER, INTENT(IN) :: size,m,n
	LOGICAL, INTENT(OUT) :: updated

	INTEGER :: i,j
	LOGICAL :: diagonal

	updated = .FALSE.

	DO j=1,n-1
		IF (d(j) .ne. d(j)+n+1) THEN
			diagonal = .FALSE.
		END IF
	END DO

	! do not update if the set is diagonal
	! any subsequent sets are transposes of previous sets
	IF (diagonal) GO TO 35

	i=size
	
	IF( d(1) .ge. m*n-size+1) GO TO 35

	! go to the next solution in the chain
	DO while (.not. updated)
		IF (d(i) .lt. m*n-size+i) THEN
			d(i) = d(i)+1
			DO j=i+1,size
				d(j)=d(i)+j-i
			END DO
			updated = .TRUE.
		ELSE
			i=i-1
		END IF
		IF (i .lt. 1) THEN
			GO TO 35 ! RETURN without updating
		END IF
	END DO

35	RETURN

END SUBROUTINE GET_NEXT_SET

!---------------------------------------------------------------------------------------
! convert a node index to a coordinate
!---------------------------------------------------------------------------------------

SUBROUTINE INDEX_TO_COORD(i,n,x,y)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: i ! index
	INTEGER, INTENT(IN) :: n ! columns in the array
	INTEGER, INTENT(OUT) :: x,y ! coordinate

	y = mod(i+n-1,n) + 1
	x = (i-y)/n + 1

END SUBROUTINE INDEX_TO_COORD

!---------------------------------------------------------------------------------------
! determine if 2 integer arrays of size n are identical
!---------------------------------------------------------------------------------------

LOGICAL FUNCTION INT_COMPARE(a,b,n) result(matching)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: a(:),b(:)
	INTEGER, INTENT(IN) :: n
	INTEGER :: i

	matching = .TRUE.
	DO i=1,n
		IF (a(i) .ne. b(i)) THEN
			matching = .FALSE.
			GO TO 37
		END IF
37	END DO

END FUNCTION INT_COMPARE

!---------------------------------------------------------------------------------------
! Sort an INTEGER array in nondecreasing order
!---------------------------------------------------------------------------------------

SUBROUTINE INT_SORT(x, size)

	IMPLICIT NONE

	INTEGER, INTENT(INOUT) :: x(:)
	INTEGER, INTENT(IN)	:: size

	INTEGER :: i,j

	DO i = 1, size-1
		j = MIN_INDEX(x, i, size)
		call INT_SWAP(x(i), x(j))
	END DO

END SUBROUTINE  INT_SORT

!---------------------------------------------------------------------------------------
! Swap two integers
!---------------------------------------------------------------------------------------

SUBROUTINE INT_SWAP(a, b)

	IMPLICIT NONE

	INTEGER, INTENT(INOUT) :: a, b

	INTEGER :: x

	x = a
	a = b
	b = x

END SUBROUTINE INT_SWAP

!---------------------------------------------------------------------------------------
! find the index of the minimum value in the integer array x(i:j)
!---------------------------------------------------------------------------------------

PURE INTEGER FUNCTION MIN_INDEX(x, i, j)

	IMPLICIT  NONE

	INTEGER, INTENT(IN) :: x(:)
	INTEGER, INTENT(IN) :: i,j

	INTEGER :: k,current,location

	current  = x(i)
	location = i

	DO k = i+1, j
		IF (x(k) < current) THEN
			current  = x(k)
			location = k
		END IF
	END DO

	MIN_INDEX = location

END FUNCTION  MIN_INDEX


END MODULE standard

