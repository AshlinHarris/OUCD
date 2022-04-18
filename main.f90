
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

INTEGER :: ZERO = 0 ! Dummy  argument for generalized subroutines

CHARACTER(LEN=20) FMT

k=8
m=4
n=4
ALLOCATE(a(m,n))

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
	CALL CHECK_COMPLETENESS(a,m,n,ZERO,is_complete)
	IF (is_complete) THEN
		allocate(adj(k,m-2,n-2),stat=ALLOCATED)
		CALL GENERATE_ADJACENCIES(m,n,k,adj,ZERO,a)
		IF(CHECK_UNIQUENESS(adj,m,n,ZERO,k)) THEN
			IF(1.0*COUNT(a)/m/n < 0.25) THEN
				! What cut-off values should be used?
				WRITE(*,*) a, ", ", 1.0*COUNT(a)/m/n, ", ", 1.0*COUNT(adj/=0)/COUNT(a)
				CALL PRINT_TILE(a,m,n)
			END IF
		END IF
	END IF
END DO

CLOSE(101)

END PROGRAM MAIN

