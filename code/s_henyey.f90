SUBROUTINE henyey
!computes the Henyey-Matrix H

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i, j, l, m, r, q

REAL (KIND = pcs) :: delta


!----------------------------------------------------------------------------------------------------



DO j = 1, k*N
  DO i = 1, k
    H((i-1)*N+1:i*N,j) = A(i,1:N)
  END DO
END DO


DO q = 1, k*N
  j = (q-1)/N + 1
  i = q - (j-1)*N
 
  delta = epsil*ABS(u(j,i))
  IF (abs(u(j,i)) < 1.D-100) delta = epsil
  u(j,i) = u(j,i) + delta
  
    CALL calc_variables

  DO r = 1, k*N
    l = (r-1)/N + 1
    m = r - (l-1)*N
    H(r,q) = (A(l,m)-H(r,q))/delta
  END DO
  
  u(j,i) = u(j,i) - delta
  
END DO


CALL calc_variables



!----------------------------------------------------------------------------------------------------

END SUBROUTINE henyey
