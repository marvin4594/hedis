SUBROUTINE store_results

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i

!----------------------------------------------------------------------------------------------------



IF (.NOT. store_by_time) THEN

  index_store_results = index_store_results + 1
  IF (index_store_results == steps_store) store_flag = .TRUE.

END IF



IF (store_flag) THEN

  store_flag = .FALSE.

  IF (      store_by_time) index_store_results = index_store_results + 1
  IF (.NOT. store_by_time) index_store_results = 0


  IF (SW_store_distrib) THEN

    DO i = 1, N
      Write(64,*) radius(i), u(1:3,i), u(2,i)/(cs*u(1,i)), - c_grav*mass_bh/(radius(i)**2), + a_g(i)
    END DO
    Write(64,*) ""
    Write(64,*) ""

  END IF


  IF (SW_store_timedev) Write(69,*) timestep, time, d_time, d_time_cfl, &
                                    mass_disk, mass_central, mass_bh, mass_lost, Mdot_disk, Mdot_edd, Mdot_bh

  Write(*,'( i9, 8es14.3, i4 )') timestep, time, d_time, d_time_cfl, mass_disk, mass_central, mass_bh, &
                                 ONE-(mass_lost+mass_disk)/mass_disk_init

END IF



IF ((time >= time_store_general) .AND. (SW_store_general) .AND. (.NOT. general_stored_yet)) THEN

  general_stored_yet = .TRUE.
  Write(71,*) disk_outer_radius, mass_disk_init, mass_bh, mass_central, mass_lost_edd

END IF



!----------------------------------------------------------------------------------------------------

END SUBROUTINE store_results
