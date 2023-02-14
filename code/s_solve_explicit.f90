SUBROUTINE solve_explicit

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

REAL (KIND = pcs), DIMENSION(1:k,1:N) :: k1, k2, k3, k4

!----------------------------------------------------------------------------------------------------



u_old = u
D_old = D


!-----Runge-Kutta 4th order----

k1 = - d_time*D
u(1:k,1:N)   = u_old(1:k,1:N) + HALF*k1


CALL calc_variables

k2 = -d_time*D
u(1:k,1:N)   = u_old(1:k,1:N) + HALF*k2


CALL calc_variables

k3 = -d_time*D
u(1:k,1:N)   = u_old(1:k,1:N) + k3


CALL calc_variables

k4 = -d_time*D
u(1:k,1:N)   = u_old(1:k,1:N) + (ONE/SIX)*( k1 + TWO*k2 + TWO*k3 + k4 )



CALL calc_variables



!----------------------------------------------------------------------------------------------------

END SUBROUTINE solve_explicit