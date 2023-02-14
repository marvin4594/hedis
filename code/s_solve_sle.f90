SUBROUTINE solve_sle
!solves the sle H du = - A

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

CHARACTER*1 :: FACT, TRANS, EQUED

INTEGER :: i

REAL (KIND = pcs) :: RCOND

REAL (KIND = pcs), DIMENSION(k*N) :: R, C, x, right

REAL (KIND = pcs), DIMENSION(1) :: FERR, BERR

INTEGER, DIMENSION(k*N) :: IPIV

INTEGER, DIMENSION(k*N) :: IWORK

REAL (KIND = pcs), DIMENSION(4*k*N) :: WORK

REAL (KIND = pcs), DIMENSION(1:k*N,1:k*N) :: AF


FACT  = 'E'
TRANS = 'N'
EQUED = 'N'

!----------------------------------------------------------------------------------------------------



DO i = 1, k
  right((i-1)*N+1:i*N) = -A(i,1:N)
END DO


CALL DGESVX( FACT, TRANS, k*N, 1, H, k*N, AF, k*N, IPIV, EQUED, R, C, right, k*N, &
             x, k*N, RCOND, FERR, BERR, WORK, IWORK, error_code )

IF (error_code /= 0) error = .TRUE.


DO i = 1, k
  d_u(i,1:N) = x((i-1)*N+1:i*N)
END DO



!----------------------------------------------------------------------------------------------------

END SUBROUTINE solve_sle