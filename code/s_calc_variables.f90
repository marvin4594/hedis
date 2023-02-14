SUBROUTINE calc_variables

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i

REAL (KIND = pcs) :: a_p, a_m

REAL (KIND = pcs), DIMENSION(-1:N+1) :: angmom
REAL (KIND = pcs), DIMENSION(0:N)    :: d2_angmom

!----------------------------------------------------------------------------------------------------



CALL calc_boundarys

a_g = MATMUL(G,u(1,0:N+1)/radius(0:N+1))



SELECT CASE (k)

  CASE (3)

    CALL calc_reconstruction

    CALL calc_viscosity


    DO i = 0, N

      a_p = MAX( u_tilde_p(2,i)/u_tilde_p(1,i)+cs, u_tilde_m(2,i)/u_tilde_m(1,i)+cs, ZERO )
      a_m = MIN( u_tilde_p(2,i)/u_tilde_p(1,i)-cs, u_tilde_m(2,i)/u_tilde_m(1,i)-cs, ZERO )
 
      F_num(:,i) = ( a_p*(/ u_tilde_m(2,i), u_tilde_m(2,i)**2/u_tilde_m(1,i) + (cs**2)*u_tilde_m(1,i), &
                             u_tilde_m(2,i)*u_tilde_m(3,i)/u_tilde_m(1,i) /) &
                    - a_m*(/ u_tilde_p(2,i), u_tilde_p(2,i)**2/u_tilde_p(1,i) + (cs**2)*u_tilde_p(1,i), &
                             u_tilde_p(2,i)*u_tilde_p(3,i)/u_tilde_p(1,i) /) &
                    - a_p*a_m*num_visc*(u_tilde_m(1:3,i) - u_tilde_p(1:3,i)) )/(a_p - a_m)
 
    END DO


    DO i = 0, N

      V(2,i) =  + (zeta_tilde(i) + FOUR_THIRD*eta_tilde(i))*radius_tilde(i) &
                   *(u_tilde(1,i)*(u(2,i+1)-u(2,i)) - u_tilde(2,i)*(u(1,i+1)-u(1,i)))/((u_tilde(1,i)**2)*d_radius_tilde(i)) &
                + (zeta_tilde(i) - TWO_THIRD*eta_tilde(i))*u_tilde(2,i)/u_tilde(1,i) 

      V(3,i) =  eta_tilde(i)*( radius_tilde(i)*(u_tilde(1,i)*(u(3,i+1)-u(3,i)) - u_tilde(3,i)*(u(1,i+1)-u(1,i))) &
                                       /(d_radius_tilde(i)*(u_tilde(1,i)**2)) - TWO*u_tilde(3,i)/u_tilde(1,i) )
 
    END DO


    DO i = 1, N

      W(2,i) =  - (zeta(i) + FOUR_THIRD*eta(i))*u(2,i)/(u(1,i)*radius(i)) &
                - (zeta(i) - TWO_THIRD*eta(i)) &
                   *(u(1,i)*(u_tilde(2,i) - u_tilde(2,i-1)) - u(2,i)*(u_tilde(1,i) - u_tilde(1,i-1)))/(d_radius(i)*(u(1,i)**2))

      S(2,i) =  u(1,i)*( + (u(3,i)/u(1,i))**2/(radius(i)**3) + (cs**2)/radius(i) - c_grav*mass_central/(radius(i)**2) + a_g(i) )

      D(:,i) = + ((F_num(:,i) - F_num(:,i-1)) - (V(:,i) - V(:,i-1))) /d_radius(i) - W(:,i) - S(:,i)

    END DO


    A(:,1:N) = ((u(:,1:N)-u_old(:,1:N))  +  d_time*( mu*D + (ONE-mu)*D_old ))



  CASE (1)

    angmom(0:N+1) = SQRT(c_grav*mass_central*radius(0:N+1) - a_g(0:N+1)*(radius(0:N+1)**3))
    angmom(-1)    = SQRT(c_grav*mass_central*radius(-1))

    u(3,0:N+1) = angmom(0:N+1)*u(1,0:N+1)

    CALL calc_viscosity

    d2_angmom(0:N) = d2(1,0:N)*angmom(1:N+1) + d2(2,0:N)*angmom(0:N) + d2(3,0:N)*angmom(-1:N-1)      !second derivative of angmom


    DO i = 1, N

      u(2,i) = (   (eta(i)*radius(i) - eta(i-1)*radius(i-1))*(angmom(i) - angmom(i-1)) &
                  - TWO*d_radius_tilde(i-1)*(eta(i)*angmom(i) - eta(i-1)*angmom(i-1)) &
                  + f(i)*( eta(i)*d2_angmom(i) + eta(i-1)*d2_angmom(i-1) ) ) &
                /((angmom(i)-angmom(i-1))*d_radius_tilde(i-1))

    END DO


    D(1,1:N) = ( mu*(u(2,2:N+1)-u(2,1:N)) + (ONE-mu)*(u_old(2,2:N+1)-u_old(2,1:N)) )/d_radius(1:N)

    A(1,1:N) =   (u(1,1:N)-u_old(1,1:N)) + d_time*D(1,1:N)

END SELECT



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_variables