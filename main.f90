
! ======================================================================================
! main program
! ======================================================================================

PROGRAM MAIN

	USE :: standard

	IMPLICIT NONE

		LOGICAL, ALLOCATABLE :: a(:,:)
		INTEGER, ALLOCATABLE :: adj(:,:,:)
		INTEGER :: i,j,x,h,m,n, ALLOCATED,k
		LOGICAL :: is_complete

		CHARACTER(LEN=20) FMT

		k=8
		m=4
		n=4
		ALLOCATE(a(m,n))

		a=.FALSE.

!		DO j=1,3
!				DO i=1,3
!						WRITE(*,'(L)',advance='no') a(i,j)
!				END DO
!		END DO

		!WRITE(*,*) a

		WRITE(FMT,*) m*n
		OPEN(UNIT=100,FILE='file')

		DO i=2**(m*n-1),2**(m*n)-1
				!write(100,*) i
				WRITE(100,"(B" // ADJUSTL(FMT) // ")") i
				!write(100,'(B9)') i
		END DO

		CLOSE(100)

		OPEN(UNIT=101,FILE='file',ACTION='read')
		DO h =2**(m*n-1),2**(m*n)-1
				DO j=1,n
						DO i=1,m
								READ (101,'(I1)',advance='no') x
								a(i,j) = (x == 1)
						END DO
				END DO
				READ (101,*)
				CALL CHECK_COMPLETENESS(a,m,n,is_complete)
				IF (is_complete) THEN
			allocate(adj(k,m-2,n-2),stat=ALLOCATED)
						CALL GENERATE_ADJACENCIES(m,n,k,adj,a)
						IF(CHECK_UNIQUENESS(adj,m,n,k)) THEN
								IF(1.0*COUNT(a)/m/n < 0.25) THEN
										! what cutoffvalues should be used?
										WRITE(*,*) a, ", ", 1.0*COUNT(a)/m/n, ", ", 1.0*COUNT(adj/=0)/COUNT(a)
										CALL PRINT_TILE(a,m,n)
								END IF
						END IF
				END IF
		END DO

		CLOSE(101)

! add to each tile

contains

LOGICAL FUNCTION CHECK_UNIQUENESS(a,m,n,k) result(unique)

! --------------------------------------------------------------------------------------
! check if every node's dominance is unique
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: a(:,:,:)
	INTEGER, INTENT(IN) :: m,n,k
	INTEGER :: i,j,i1,i2

	INTEGER :: MAX_ROW, MAX_COLUMN

	! Adjacencies are not accurate for exterior nodes.
	! Only interior nodes of the augmented tile are considered.
	! These are the only nodes stored in the adjacency lists.
	! So, the size of the adjacency matrix is (k,m+2*x-2,n+2*x-2)
	
	MAX_ROW	= m-2
	MAX_COLUMN = n-2

	! for each node ...
	DO j=1,MAX_COLUMN
	DO i=1,MAX_ROW

		! ... consider all subsequent nodes
		! the loop is broken into 2 parts to avoid the case i = j
		! first the nodes in the same column
		i2=j
		DO i1=i+1,MAX_ROW
			IF ( INT_COMPARE(a(:,i,j), a(:,i1,i2),k) ) THEN
				unique = .FALSE.
				GO TO 39 ! return immediately
			END IF
		END DO

			   ! second, all others
		DO i2=j+1,MAX_COLUMN
		DO i1=1,MAX_ROW

			IF ( INT_COMPARE(a(:,i,j), a(:,i1,i2),k) ) THEN
				unique = .FALSE.
				GO TO 39 ! return immediately
			END IF

		END DO
		END DO

	END DO
	END DO

	! since every adjacency was different, the dominance must be unique
	unique = .TRUE.

39	RETURN

END FUNCTION CHECK_UNIQUENESS

SUBROUTINE CHECK_COMPLETENESS(s,m,n,completeDominance)

! --------------------------------------------------------------------
! determine if a tiling exhibits complete dominance
! by examining interior points of an m by n tiling
! --------------------------------------------------------------------

	IMPLICIT NONE

	LOGICAL, INTENT(IN) :: s(:,:) ! augmented tile
	INTEGER, INTENT(IN) :: m,n

	LOGICAL, INTENT(OUT) :: completeDominance

	INTEGER :: i,j,i1,j2 ! counters
	LOGICAL :: dominated

	INTEGER :: MAX_ROW, MAX_COLUMN

	! loop through each interior point
	DO j=2,n-1
	DO i=2,m-1

		dominated = .false.

		! loop through each adjacent point
		DO j2=j-1,j+1
		DO i1=i-1,i+1

			! check for dominating nodes
			! node cannot dominate self
			IF (s(i1,j2) .and. .not. ( i1 .eq. i .and. j2 .eq. j)  ) THEN
				dominated = .true.
			END IF

		END DO
		END DO

		! cancel the search if dominance is not complete
		IF (.not. dominated) THEN
			completeDominance = .FALSE.
			GO TO 370
		END IF

	END DO
	END DO

	! each point is dominated
	completeDominance = .TRUE.

370	RETURN

END SUBROUTINE CHECK_COMPLETENESS

SUBROUTINE GENERATE_ADJACENCIES(m,n,k,a,d)

! --------------------------------------------------------------------
! return the adjacency arrays of an m by n tile with
! n_d dominating nodes and connectivity k
! --------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: m ! rows
	INTEGER, INTENT(IN) :: n ! columns
	INTEGER, INTENT(IN) :: k ! connectivity

	! the adjacencies array uses coordinates from the augmented grid
	INTEGER, INTENT(INOUT) :: a(:,:,:) ! adjacencies
	LOGICAL, INTENT(IN) :: d(:,:) ! augmented tile

	INTEGER :: i,j,i2,j2! counters
	INTEGER :: MAX_ROW,MAX_COLUMN

	!INTEGER, allocatable :: adjacent_nodes(:)
	INTEGER :: index

	MAX_ROW = m-2
	MAX_COLUMN = n-2

	! initialize adjacency arrays
	a=0

!$OMP PARALLEL SHARED(d,a,MAX_COLUMN,MAX_ROW,n,k) PRIVATE(i,j,i2,j2,index)
!$OMP DO

	! for each interior node of the augmented graph ...
	DO j=1,MAX_COLUMN
	DO i=1,MAX_ROW
		
		! ... record the indices of adjacent dominant nodes.
		! Array indices for the space are off by one.

		! nodes to the left
		j2=j
		DO i2=i,i+2
			IF (d(i2,j2)) THEN
				call COORD_TO_INDEX(index,n,i2,j2)
				a(i2-i+1,i,j) = index
			END IF
		END DO

		! node above
		j2=j+1
		i2=i
		IF (d(i2,j2)) THEN
			call COORD_TO_INDEX(index,n,i2,j2)
			a(4,i,j) = index
		END IF

		! node below
		j2=j+1
		i2=i+2
		IF (d(i2,j2)) THEN
			call COORD_TO_INDEX(index,n,i2,j2)
			a(5,i,j) = index
		END IF

		! nodes to the right
		j2=j+2
		DO i2=i,i+2
			IF (d(i2,j2)) THEN
				call COORD_TO_INDEX(index,n,i2,j2)
				a(i2-i+6,i,j) = index
			END IF
		END DO

		! INT_SORT the adjacencies to enable comparison
		call INT_SORT(a(:,i,j),k)
		
	END DO
	END DO

!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE GENERATE_ADJACENCIES

END PROGRAM MAIN

