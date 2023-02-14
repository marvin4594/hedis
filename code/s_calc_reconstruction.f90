SUBROUTINE calc_reconstruction

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i

REAL (KIND = pcs), EXTERNAL :: minmod

!----------------------------------------------------------------------------------------------------



SELECT CASE (mode_average)

  CASE (1)     !geometric mean

    u_tilde(1,:) = SQRT( u(1,0:N)*u(1,1:N+1) )
    u_tilde(2,:) = SQRT( MAX(u(2,0:N)*u(2,1:N+1),ZERO) ) * SIGN(ONE,u(2,1:N+1))
    u_tilde(3,:) = SQRT( u(3,0:N)*u(3,1:N+1) )



  CASE (2)     !arithmetic mean

    u_tilde = HALF*( u(:,0:N)*u(:,1:N+1) )

END SELECT





SELECT CASE (mode_flux)

  CASE (1)

    u_tilde_p = u_tilde
    u_tilde_m = u_tilde



  CASE (2)

    DO i = 1, N
      u_dev(:,i) = (/ MINMOD((u(1,i+1)-u(1,i))/d_radius_tilde(i),(u(1,i)-u(1,i-1))/d_radius_tilde(i-1)), &          !u_dev = 0 at i=0 and i=N+1
                      MINMOD((u(2,i+1)-u(2,i))/d_radius_tilde(i),(u(2,i)-u(2,i-1))/d_radius_tilde(i-1)), &
                      MINMOD((u(3,i+1)-u(3,i))/d_radius_tilde(i),(u(3,i)-u(3,i-1))/d_radius_tilde(i-1))  /)
    END DO
 
    DO i = 0, N
      u_tilde_p(:,i) = u(:,i+1) + u_dev(:,i+1) *(radius_tilde(i)-radius(i+1))
      u_tilde_m(:,i) = u(:,i)   + u_dev(:,i)   *(radius_tilde(i)-radius(i))
    END DO

END SELECT



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_reconstruction