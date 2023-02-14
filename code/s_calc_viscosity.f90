SUBROUTINE calc_viscosity

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



SELECT CASE (mode_viscosity)

  CASE (1)     !beta-viscosity, bulk-viscosity=zero

    eta(0:N+1) = beta*u(3,0:N+1)/radius(0:N+1)
    eta_tilde(0:N) = beta*u_tilde(3,0:N)/radius_tilde(0:N)



  CASE (2)     !Pringle-viscosity, bulk-viscosity=zero

    eta(0:N+1) = beta*u(1,0:N+1)
    eta_tilde(0:N) = beta*u_tilde(1,0:N)

END SELECT



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_viscosity