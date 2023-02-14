PROGRAM main

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: j, l

!----------------------------------------------------------------------------------------------------
!parameterstudy: M_BH after 10^9 years in dependence of s_out and M_disk_init using the Pringle-model


!----------physical parameters-----------

 disk_inner_radius  =  1.D-3
 disk_outer_radius  =  300.0d0
 log_mass_central   =  4.0d0
 log_mass_disk      =  10.0d0
 beta               =  0.01d0
 cs                 =  9.28D-7
 accr_eff           =  0.1d0
 time_end           =  1.D+9

 SW_ag              =  .TRUE.                               !allow growth of central object
 SW_cem             =  .TRUE.                               !consider escaped mass for gravitational acceleration

 mode_gravitation   =  2                                    !1: Kepler, 2: monopole (low order), 3: monopole (high order), 4: selfgravitation
 k                  =  1                                    !specify disk-model, 1: simple model (one equation), 3: extended model (three equations)
 mode_viscosity     =  1                                    !1: beta-viscosity, 2: Pringle-viscosity
 mode_initial_cond  =  1

!----------numerical parameters----------

 N                  =  50                                   !number of grid points
 part               =  0.5d0                                !part of grid points between disk_inner_radius and disk_outer_radius
 c_cfl              =  0.1d0
 c_itw              =  0.01d0
 d_time_min         =  1.D-10
 mu                 =  0.8d0                                !Crank-Nicholson-Parameter

 mode_flux          =  1                                    !1: averages, 2: upwind-scheme
 mode_average       =  1                                    !1: geometric, 2: arithmetic
 mode_boundarys     =  4                                    !1: non-reflecting, 2: reflecting, 3: no-gradients, 4: simple model, 5: pringle-disk, 6: set by initial conditions
 mode_time_integr   =  2                                    !1: explicit, 2: implicit

 accuracy           =  1.D-4
 epsil              =  1.D-4                                !for calculating H=dA/du in subroutine "henyey"
 iter_max           =  10                                   !maximum iterations for solving nonlinear system of equations
 num_visc           =  1.0d0

!----------store options-----------------

 store_by_time      =  .FALSE.                              !store after specified time or after specified number of timesteps
 time_store         =  1.D+6                                !store after this time
 steps_store        =  100                                  !store after this number of timesteps
 time_store_general =  1.D+9                                !store general results after this time

 SW_store_distrib   =  .FALSE.
 SW_store_timedev   =  .FALSE.
 SW_store_general   =  .TRUE.

 filename           =  "parameterstudy_pri"

!----------------------------------------


!----------------------------------------------------------------------------------------------------





IF (SW_store_general) OPEN ( UNIT = 71, FILE = "results/"//TRIM(filename)//"_general.dat" )

  log_mass_disk = 8.0d0
  DO l = 1, 5

    disk_outer_radius = 200.0d0
    DO j = 1, 25

      CALL disk

    disk_outer_radius = disk_outer_radius + 100.0d0
    END DO

    Write(71,*) ""
    Write(71,*) ""

  log_mass_disk = log_mass_disk + 0.5d0
  END DO


IF (SW_store_general) CLOSE ( UNIT = 71 )





!----------------------------------------------------------------------------------------------------

END PROGRAM main










SUBROUTINE initial_conditions

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

REAL (KIND = pcs) :: x

REAL (KIND = pcs), DIMENSION(1:N) :: sigma

!----------------------------------------------------------------------------------------------------



SELECT CASE (mode_initial_cond)



CASE(1)

   sigma(1:N_disk) = SQRT(ONE - ((radius(1:N_disk)-disk_inner_radius)/(disk_outer_radius-disk_inner_radius))**2)
   x = mass_disk*(ONE-1.D-3)/SUM(sigma(1:N_disk)*ring_area(1:N_disk))
   sigma(1:N_disk) = x*sigma(1:N_disk)

   sigma(N_disk+1:N) = ONE/(radius(N_disk+1:N)**2)
   x = mass_disk*1.D-3/SUM(sigma(N_disk+1:N)*ring_area(N_disk+1:N))                                 !distribution for outer region
   sigma(N_disk+1:N) = x*sigma(N_disk+1:N)                                                          !put 1 promille of total mass into outer region


   u(1,1:N) = sigma(1:N)*radius(1:N)
   a_g = MATMUL(G,u(1,0:N+1)/radius(0:N+1))
   u(3,1:N) = u(1,1:N)*SQRT(c_grav*mass_central*radius(1:N) - a_g(1:N)*(radius(1:N)**3))
   u(2,1:N) = ZERO




CASE (2)

   sigma(1:N_disk) = mass_disk*(ONE-1.D-3)/(c_pi*(disk_outer_radius**2 - disk_inner_radius**2))     !constant distribution for the disk

   sigma(N_disk+1:N) = ONE/(radius(N_disk+1:N)**2)
   x = mass_disk*1.D-3/SUM(sigma(N_disk+1:N)*ring_area(N_disk+1:N))                                 !distribution for outer region
   sigma(N_disk+1:N) = x*sigma(N_disk+1:N)                                                          !put 1 promille of total mass into outer region


   u(1,1:N) = sigma(1:N)*radius(1:N)
   a_g = MATMUL(G,u(1,0:N+1)/radius(0:N+1))
   u(3,1:N) = u(1,1:N)*SQRT(c_grav*mass_central*radius(1:N) - a_g(1:N)*(radius(1:N)**3))
   u(2,1:N) = ZERO



END SELECT



u_old = u


!----------------------------------------------------------------------------------------------------

END SUBROUTINE initial_conditions