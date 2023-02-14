SUBROUTINE calc_masses

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



Mdot_i_old = Mdot_i
Mdot_o_old = Mdot_o

Mdot_i = TWO*c_pi*u(2,1)
Mdot_o = TWO*c_pi*u(2,N+1)


Mdot_disk = - HALF*(Mdot_i+Mdot_i_old)
Mdot_edd = mass_bh/(accr_eff*tau_s)
Mdot_bh = MAX(MIN(Mdot_disk, Mdot_edd),ZERO)

mass_bh = mass_bh + d_time*Mdot_bh


IF (SW_ag) THEN
  IF (SW_cem) THEN
    mass_central = mass_central + d_time*Mdot_disk
  ELSE
    mass_central = mass_bh
  END IF
END IF


mass_disk = TWO*c_pi*SUM(u(1,1:N)*d_radius(1:N))


mass_lost = mass_lost + d_time*HALF*( Mdot_o - Mdot_i + Mdot_o_old - Mdot_i_old )     !total mass loss

mass_lost_edd = mass_lost_edd + d_time*MAX(Mdot_disk-Mdot_edd,ZERO)                   !mass loss due to Eddington-limit



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_masses
