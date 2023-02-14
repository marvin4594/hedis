SUBROUTINE init_gravitation

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i, j

REAL (KIND = pcs), EXTERNAL :: H2F1

!----------------------------------------------------------------------------------------------------



G = ZERO     !Kepler-gravitation


SELECT CASE (mode_gravitation)

  CASE (2)     !monopole approximation, low order

    DO i = 0, N + 1
      DO j = 0, i
        G(i,j) = - c_grav*ring_area(j)/(radius(i)**2)
      END DO
    END DO



  CASE (3)     !monopole approximation, high order

    DO i = 1, N

      DO j = 1, i
       G(i,j-1) = TWO*c_pi*(  - (ONE/THREE)*(radius(j)**3-radius(j-1)**3)/(radius(j)-radius(j-1)) &
                     + (ONE/FOUR)*(radius(j)**2-radius(j-1)**2) &
                     + (ONE/FOUR)*(radius(j)+radius(j-1))*(radius(j)**2-radius(j-1)**2)/(radius(j)-radius(j-1)) )
      END DO
      DO j = 1, i
       G(i,j)   = G(i,j) + TWO*c_pi*(  + (ONE/THREE)*(radius(j)**3-radius(j-1)**3)/(radius(j)-radius(j-1)) &
                   + (ONE/FOUR)*(radius(j)**2-radius(j-1)**2) &
                   - (ONE/FOUR)*(radius(j)+radius(j-1))*(radius(j)**2-radius(j-1)**2)/(radius(j)-radius(j-1))  )
      END DO

      G(i,0:N+1) = -c_grav*G(i,0:N+1)/(radius(i)**2)

    END DO



  CASE (4)     !full selfgravitation

    DO i = 0, N + 1
      DO j = 0, N + 1

        G(i,j) = - TWO*c_pi*c_grav* &
                   (  (radius(i)*(radius_tilde(j)**2)/(TWO*((radius(i)+radius_tilde(j))**3))) * &
                           H2F1( FOUR*radius(i)*radius_tilde(j)/((radius(i)+radius_tilde(j))**2)) &
                    - (radius(i)*(radius_tilde(j-1)**2)/(TWO*((radius(i)+radius_tilde(j-1))**3))) * &
                           H2F1( FOUR*radius(i)*radius_tilde(j-1)/((radius(i)+radius_tilde(j-1))**2))  )

      END DO
    END DO

END SELECT



!----------------------------------------------------------------------------------------------------

END SUBROUTINE init_gravitation