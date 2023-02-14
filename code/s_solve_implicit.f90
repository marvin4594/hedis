SUBROUTINE solve_implicit

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: iter

!----------------------------------------------------------------------------------------------------



d_u = u - u_old
u_old = u
u = u_old + d_u
D_old = D
d_u = ZERO


iter = 0


DO

  iter = iter + 1

  CALL calc_variables

  CALL henyey

  CALL solve_sle

  u = u + d_u

  IF ((MAXVAL(ABS(d_u/u)) < accuracy) .OR. (iter == iter_max)) EXIT

END DO


CALL calc_variables



!----------------------------------------------------------------------------------------------------

END SUBROUTINE solve_implicit