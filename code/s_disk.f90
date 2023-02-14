SUBROUTINE disk

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



CALL alloc_arrays

CALL init_parameters
CALL init_grid
CALL init_gravitation
CALL initial_conditions


CALL calc_variables

IF (SW_store_timedev) OPEN ( UNIT = 69, FILE = "results/"//TRIM(filename)//"_timedev.dat" )
IF (SW_store_distrib) OPEN ( UNIT = 64, FILE = "results/"//TRIM(filename)//"_distrib.dat" )

CALL store_results



DO


  CALL calc_timestep

  timestep = timestep + 1
  time = time + d_time


  SELECT CASE (mode_time_integr)
    CASE (1)
      CALL solve_explicit
    CASE (2)
      CALL solve_implicit
  END SELECT
 

  CALL calc_masses

  CALL store_results

  IF ((time > time_end) .OR. error) EXIT


END DO



IF (error) Write(*,*) "ERROR, calculation aborted. Error-code: ", error_code              !error-code=0: timestep too small, else: solving sle failed


CALL dealloc_arrays


IF (SW_store_distrib) CLOSE ( UNIT = 64 )
IF (SW_store_timedev) CLOSE ( UNIT = 69 )



!----------------------------------------------------------------------------------------------------

END SUBROUTINE disk
