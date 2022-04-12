
! ======================================================================================
! main program
! ======================================================================================

PROGRAM MAIN

! --------------------------------------------------------------------------------------
! calls the SEARCH_FOR_SOLUTION SUBROUTINE, which RETURNs the optimal solution, for each size
! of tile i x j such that MIN_SIZE <= j <= i <= MAX_SIZE and for any number of
! dominating nodes from i*j/8 (fewer nodes cannot completely dominate the space)
! and up. A solution is guaranteed with i*j/2 dominating nodes, so no upper
! bound is implemented. 
! --------------------------------------------------------------------------------------

	USE :: standard

	IMPLICIT NONE

	INTEGER, allocatable :: solution(:) ! array of dominating indices
	INTEGER, parameter :: MIN_SIZE = 3  ! minimum number of rows or columns
	INTEGER, parameter :: MAX_SIZE = 7 ! maximum number of rows or columns
	INTEGER :: MIN_D ! fewest dominating nodes to put on a tile
		LOGICAL, allocatable :: matrix(:,:)
	INTEGER :: i,j,k
	LOGICAL :: solved
	INTEGER :: ALLOCATED, DEALLOCATED

	REAL(dp) :: start, finish, time_elapsed
	REAL(dp) :: PROGRAM_start, program_finish, total_time_elapsed
	INTEGER :: time(4)
	INTEGER :: count_first, count_last, count_start, count_finish, count_rate, count_max

	write(*,*)

	CALL SYSTEM_CLOCK(count_first, count_rate, count_max)
	call cpu_time(PROGRAM_start)
	DO j=MIN_SIZE,MAX_SIZE
		DO i=j,MAX_SIZE

			! It takes at least 1/8 of all nodes to dominate the space.
			! It takes at least 7/64 additional nodes to uniquely dominate the space
			MIN_D = floor(i*j*0.234375_dp) 

			! increment number of nodes until a solution can be found
			k = MIN_D

			CALL SYSTEM_CLOCK(count_start, count_rate, count_max)
			call cpu_time(start)
			solved = .FALSE.
			DO while (.not. solved) ! solution is guaranteed
				k = k +1
				! optimal unique complete dominance
				
				call SEARCH_FOR_SOLUTION(solution,i,j,k,solved)
			END DO
			call cpu_time(finish)
			CALL SYSTEM_CLOCK(count_finish, count_rate, count_max)

			! print results

67			format("OPTIMAL SOLUTION FOR "i2," x ",i2," TILING REQUIRES ",i2, " NODES:")
			write(*,67) i,j,k

			call PRINT_NODES(solution,k)

			allocate(matrix(i,j),stat=ALLOCATED)
			IF (ALLOCATED /= 0) STOP "*** Not enough memory ***"
			call LIST_TO_MATRIX(solution,k,matrix,j)
			call PRINT_TILE(matrix,i,j)
			deallocate(matrix,stat=DEALLOCATED)
			IF (DEALLOCATED /= 0) STOP "*** Could not deallocate ***"
68			format('PERCENT SIZE OF DOMINATING SET: ',i2,' / ',i3,' = ',f4.1,'%')
			write(*,68) k,i*j,100.0_dp*k/i/j

			time_elapsed=finish-start
			call FORMAT_TIME(time_elapsed,time)
69			format(a10,2(i2.2,':'),i2.2,'.',i3.3)
	
			write(*,69) 'CPU TIME:', time
			call FORMAT_TIME((count_finish-count_start)/1000.0_dp,time)
			write(*,69) 'REAL TIME:', time
			write(*,*)

		END DO
	END DO
	call cpu_time(PROGRAM_finish)
	CALL SYSTEM_CLOCK(count_last, count_rate, count_max)

	total_time_elapsed = PROGRAM_finish - program_start
71	format('FINISHED ALL TILINGS FROM ',i2,'x',i2,' to ',i2,'x',i2)
	write(*,71) MIN_SIZE, MIN_SIZE, MAX_SIZE, MAX_SIZE
	call FORMAT_TIME(total_time_elapsed,time)
	write(*,69) 'TOTAL CPU TIME: ', time
	call FORMAT_TIME((count_last-count_first)/1000.0_dp,time)
	write(*,69) 'TOTAL REAL TIME: ', time
	write(*,*)

contains

LOGICAL FUNCTION CHECK_UNIQUENESS(a,m,n,x,k) result(unique)

! --------------------------------------------------------------------------------------
! check if every node's dominance is unique
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: a(:,:,:)
	INTEGER, INTENT(IN) :: m,n,x,k
	INTEGER :: i,j,i1,i2

	INTEGER :: MAX_ROW, MAX_COLUMN

	! Adjacencies are not accurate for exterior nodes.
	! Only interior nodes of the augmented tile are considered.
	! These are the only nodes stored in the adjacency lists.
	! So, the size of the adjacency matrix is (k,m+2*x-2,n+2*x-2)
	
	MAX_ROW	= m+2*x-2
	MAX_COLUMN = n+2*x-2

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

SUBROUTINE CHECK_COMPLETENESS(s,m,n,x,completeDominance)

! --------------------------------------------------------------------
! determine if a tiling exhibits complete dominance
! by examining interior points of an m+2*x-2 by n+*2x-2 tiling
! --------------------------------------------------------------------

	IMPLICIT NONE

	LOGICAL, INTENT(IN) :: s(:,:) ! augmented tile
	INTEGER, INTENT(IN) :: m,n,x

	LOGICAL, INTENT(OUT) :: completeDominance

	INTEGER :: i,j,i1,j2 ! counters
	LOGICAL :: dominated

	INTEGER :: MAX_ROW, MAX_COLUMN

	MAX_ROW	= m+2*x-1
	MAX_COLUMN = n+2*x-1

	! loop through each interior point
	DO j=2,MAX_COLUMN
	DO i=2,MAX_ROW

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

SUBROUTINE GENERATE_ADJACENCIES(m,n,k,a,x,d)

! --------------------------------------------------------------------
! return the adjacency arrays of an m+x by n+x tile with
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

	INTEGER, INTENT(IN) :: x
	!INTEGER, allocatable :: adjacent_nodes(:)
	INTEGER :: index

	MAX_ROW = m+2*x-2
	MAX_COLUMN = n+2*x-2

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

SUBROUTINE LIST_TO_MATRIX(d,size,g,n)

! --------------------------------------------------------------------------------------
! convert a list of nodes to a boolean matrix
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER :: d(:) ! array of dominating nodes
	INTEGER, INTENT(IN) :: size  ! elements of array d
	INTEGER, INTENT(IN) :: n	 ! columns in tile
	LOGICAL :: g(:,:) ! boolean matrix of nodes

	INTEGER :: i,j,k
	INTEGER :: index

	g = .FALSE.

	! mark each node as TRUE
	DO k=1,size
		index = d(k)
		IF (index .ge. 1) THEN
			call INDEX_TO_COORD(index,n,i,j)
			g(i,j) = .TRUE.
		END IF
	END DO

END SUBROUTINE LIST_TO_MATRIX

SUBROUTINE PRINT_NODES(d,n)

! --------------------------------------------------------------------------------------
! print a set of node indices to standard output
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: d(:)
	INTEGER, INTENT(IN) :: n

	INTEGER :: i

	write(*,'(a)',advance="no") '{'

	DO i=1,n-1
		write(*,'(i3,", ")',advance="no") d(i)
	END DO

	write(*,'(i3,a)') d(n),'}'

END SUBROUTINE PRINT_NODES

SUBROUTINE PRINT_TILE( tile, m, n )

! --------------------------------------------------------------------------------------
! print a tile as a character grid
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: m,n
	LOGICAL :: tile(:,:)
	INTEGER :: i,j

	! print a line on top
	write (*,"(a)",advance="no") '/'
	DO i=1,m
			write (*,"(a)",advance="no") '-'
	END DO
	write (*,"(a)") '\'
		
	! print the tile with lines on the left and right
	DO j=1,n

		write (*,"(a)",advance="no") '|'

		DO i=1,m

			IF(tile(i,j)) THEN
				write (*,"(a)",advance="no") '@'
			ELSE
				write (*,"(a)",advance="no") ' '
			END IF
				
		END DO
		write (*,"(a)") '|'

	END DO

	! print a line at the base
	write (*,"(a)",advance="no") '\'
	DO i=1,m
			write (*,"(a)",advance="no") '-'
	END DO
	write (*,"(a)") '/'

END SUBROUTINE PRINT_TILE

SUBROUTINE SEARCH_FOR_SOLUTION(d_nodes,m,n,n_d,solved)

! --------------------------------------------------------------------------------------
! find a tiling that dominates all nodes uniquely and completely, if it exists
! for a given size of tile and number of dominating nodes
! --------------------------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: m ! rows in tile
	INTEGER, INTENT(IN) :: n ! columns in tile
	INTEGER, INTENT(IN) :: n_d

	INTEGER, allocatable, INTENT(OUT) :: d_nodes(:)
	LOGICAL, INTENT(OUT)			  :: solved

	LOGICAL, allocatable :: g(:,:)	 ! tile
	LOGICAL, allocatable :: aug(:,:)   ! augmented tile
	INTEGER, allocatable :: adj(:,:,:) ! adjacencies

	INTEGER, parameter :: k = 8 ! connectivity of graph
	INTEGER, parameter :: x = 2 ! border width

	INTEGER :: ALLOCATED	 ! status of allocation
	INTEGER :: DEALLOCATED   ! status of deallocation
	LOGICAL :: continue_flag ! can a solution still be found?
	LOGICAL :: isDominated   ! does a set exhibit complete dominance?

	allocate(d_nodes(n_d),stat=ALLOCATED)
	IF (ALLOCATED /= 0) STOP "*** Not enough memory ***"

	call GET_FIRST_SET(d_nodes,n_d)

	continue_flag = .TRUE.
	solved = .FALSE.

	DO while (continue_flag)

		allocate(g(m,n),stat=ALLOCATED)
		IF (ALLOCATED /= 0) STOP "*** Not enough memory ***"

		call LIST_TO_MATRIX(d_nodes,n_d,g,n)
	
		! create m+2x by n+2x graph containing g and a border
		allocate(aug(m+2*x,n+2*x),stat=ALLOCATED)
		IF (ALLOCATED /= 0) STOP "*** Not enough memory ***"

		call GENERATE_TILING(g,m,n,x,aug)

		!check for completeness
		call CHECK_COMPLETENESS(aug,m,n,x,isDominated)

		IF (isDominated) THEN

			allocate(adj(k,m+2*x-2,n+2*x-2),stat=ALLOCATED)
			IF (ALLOCATED /= 0) STOP "*** Not enough memory ***"

			! check for CHECK_UNIQUENESS
			call GENERATE_ADJACENCIES(m,n,k,adj,x,aug)

			IF (CHECK_UNIQUENESS(adj,m,n,x,k)) THEN
				solved = .TRUE.
			ELSE
				solved = .FALSE.
			END IF

			deallocate(adj,stat=DEALLOCATED)
			IF (DEALLOCATED /= 0) STOP "*** Could not deallocate ***"

		ELSE
			solved = .FALSE.
		END IF

		IF(solved) THEN
			continue_flag = .FALSE.
		ELSE
			! try to update d_nodes
			call GET_NEXT_SET(d_nodes,n_d,m,n,continue_flag)
		END IF

		deallocate(g,stat=DEALLOCATED)
		IF (DEALLOCATED /= 0) STOP "*** Could not deallocate ***"
		deallocate(aug,stat=DEALLOCATED)
		IF (DEALLOCATED /= 0) STOP "*** Could not deallocate ***"

	END DO

END SUBROUTINE SEARCH_FOR_SOLUTION

END PROGRAM MAIN

