SUBROUTINE calc_timestep

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

REAL (KIND = pcs) :: eps, t_1, t_2, t_3

!----------------------------------------------------------------------------------------------------



t_1 = MINVAL(d_radius(0:N+1)/( ABS(u(2,0:N+1)/u(1,0:N+1)) + cs ))                                   !cfl-condition for advection
t_2 = MINVAL(SQRT(ABS(d_radius(0:N+1)/( - c_grav*mass_central/(radius(0:N+1)**2) + a_g(0:N+1) ))))  !cfl-condition for grav. acceleration
t_3 = MINVAL(ABS( (d_radius(0:N+1)**2)*u(1,0:N+1)/(beta*u(3,0:N+1)) ))                              !cfl-condition for beta-viscosity

d_time_cfl = c_cfl*MIN(t_1,t_2,t_3)



eps = MAXVAL(ABS(ONE - u_old(1,0:N+1)/u(1,0:N+1)))

d_time_itw = MAX(d_time_cfl, (c_itw/eps)*d_time)



SELECT CASE (mode_time_integr)
  CASE (1)
    d_time = d_time_cfl
  CASE (2)
    d_time = d_time_itw
END SELECT



IF (d_time < d_time_min) error = .TRUE.


IF ((store_by_time) .AND. (( time+d_time ) > index_store_results*time_store)) THEN               !passed time to store?
  d_time = index_store_results*time_store - time                                                 !cut timestep
  store_flag = .TRUE.                                                                            !store results
END IF


IF (SW_store_general .AND. (.NOT. general_stored_yet) .AND. ((time+d_time) > time_store_general)) &     !cut timestep for storing general results
   d_time = time_store_general - time



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_timestep
