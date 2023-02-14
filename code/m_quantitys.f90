MODULE quantitys

USE constants

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



CHARACTER (LEN=128) :: filename

INTEGER :: mode_average, mode_boundarys, mode_flux, mode_gravitation, &
           mode_initial_cond, mode_time_integr, mode_viscosity
INTEGER :: N, N_disk, k
INTEGER :: index_store_results, steps_store
INTEGER :: error_code
INTEGER :: iter_max
INTEGER :: timestep

LOGICAL :: error
LOGICAL :: general_stored_yet, store_by_time, store_flag
LOGICAL :: SW_store_distrib, SW_store_general, SW_store_timedev
LOGICAL :: SW_ag, SW_cem

REAL (KIND = pcs) :: accuracy, epsil, mu, num_visc, part
REAL (KIND = pcs) :: accr_eff, beta, cs
REAL (KIND = pcs) :: c_cfl, c_itw, d_time, d_time_cfl, d_time_itw, d_time_min, time, time_end, time_store, time_store_general
REAL (KIND = pcs) :: disk_inner_radius, disk_outer_radius
REAL (KIND = pcs) :: mass_bh, mass_central, mass_disk, mass_disk_init, mass_lost, mass_lost_edd
REAL (KIND = pcs) :: Mdot_i, Mdot_o, Mdot_i_old, Mdot_o_old, Mdot_bh, Mdot_disk, Mdot_edd
REAL (KIND = pcs) :: log_mass_central, log_mass_disk
REAL (KIND = pcs) :: s_0, sigma_0


REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: a_g

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: eta
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: eta_tilde
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: zeta
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: zeta_tilde

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: d_radius
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: d_radius_tilde
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: radius
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: radius_tilde
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:) :: ring_area

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: d_u
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u_dev
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u_old

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u_tilde
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u_tilde_m
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: u_tilde_p

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: A
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: D
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: D_old

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: F_num
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: V
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: W
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: S

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: G
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) :: H

REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:)   ::         f
REAL (KIND = pcs), ALLOCATABLE, DIMENSION(:,:) ::         d2              



!----------------------------------------------------------------------------------------------------

END MODULE quantitys
