SUBROUTINE init_grid

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i

!----------------------------------------------------------------------------------------------------



DO i = -2, N+1
  radius_tilde(i) = disk_inner_radius * TEN**( (i) * LOG10(disk_outer_radius/disk_inner_radius) / N_disk )
END DO


radius(-1:N+1) = SQRT(radius_tilde(-1:N+1)*radius_tilde(-2:N))

d_radius(0:N+1) = radius_tilde(0:N+1) - radius_tilde(-1:N)

d_radius_tilde(-1:N) = radius(0:N+1) - radius(-1:N)

ring_area(0:N+1) = TWO*c_pi*radius(0:N+1)*d_radius(0:N+1)



DO i = 0, N
  d2(1,i) =     ONE/(d_radius(i)*d_radius_tilde(i))
  d2(2,i) = - ( ONE/(d_radius(i)*d_radius_tilde(i)) + ONE/(d_radius(i)*d_radius_tilde(i-1)) )
  d2(3,i) =     ONE/(d_radius(i)*d_radius_tilde(i-1))
END DO


f(1:N) = HALF*(d_radius_tilde(0:N-1)**2)*radius_tilde(0:N-1)



!----------------------------------------------------------------------------------------------------

END SUBROUTINE init_grid
