
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

END PROGRAM MAIN

