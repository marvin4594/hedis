SUBROUTINE init_parameters

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



Write(*,*) ""
Write(*,*) "timestep     time/yr       d_t/yr        d_t(cfl)/yr   M_disk/M_sun  M_z/M_sun     M_BH/M_sun    Mass conserv."
Write(*,*) "--------------------------------------------------------------------------------------------------------------"



N_disk = part*N

index_store_results = 0
store_flag = .TRUE.
general_stored_yet = .FALSE.

error = .FALSE.

timestep = 0

error_code = 0

d_time     = ZERO
d_time_cfl = ZERO
d_time_itw = ZERO
time       = ZERO


mass_disk    = TEN**log_mass_disk
mass_central = TEN**log_mass_central

mass_disk_init = mass_disk

mass_bh = mass_central

mass_lost     = ZERO
mass_lost_edd = ZERO

Mdot_i = ZERO
Mdot_o = ZERO

eta        = ZERO
eta_tilde  = ZERO
zeta       = ZERO
zeta_tilde = ZERO


A = ZERO

D     = ZERO
D_old = ZERO

S  = ZERO
V  = ZERO
W  = ZERO

d_u     = ZERO
u       = ZERO
u_dev   = ZERO
u_old   = ZERO
u_tilde = ZERO



!----------------------------------------------------------------------------------------------------

END SUBROUTINE init_parameters
