SUBROUTINE calc_boundarys

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: i

REAL (KIND = pcs) :: v_in, v_out, xi, tau

REAL (KIND = pcs), DIMENSION (1:3,0:N+1) :: charak          !characteristic variables

!----------------------------------------------------------------------------------------------------



SELECT CASE (mode_boundarys)

  CASE(1)     !non-reflecting boundary conditions

    charak = ZERO

    i = 2
    charak(:,i) = (/ cs*LOG((u(1,i)*radius(i-1))/(u(1,i-1)*radius(i))) - (u(2,i)/u(1,i)-u(2,i-1)/u(1,i-1)), &
                     ( (u(3,i)/(u(1,i)*radius(i))) - (u(3,i-1)/(u(1,i-1)*radius(i-1))) ), &
                     cs*LOG((u(1,i)*radius(i-1))/(u(1,i-1)*radius(i))) + (u(2,i)/u(1,i)-u(2,i-1)/u(1,i-1)) /)

    i = N
    charak(:,i) = (/ cs*LOG((u(1,i)*radius(i-1))/(u(1,i-1)*radius(i))) - (u(2,i)/u(1,i)-u(2,i-1)/u(1,i-1)), &
                     ( (u(3,i)/(u(1,i)*radius(i))) - (u(3,i-1)/(u(1,i-1)*radius(i-1))) ), &
                     cs*LOG((u(1,i)*radius(i-1))/(u(1,i-1)*radius(i))) + (u(2,i)/u(1,i)-u(2,i-1)/u(1,i-1)) /)


    v_in = u(2,2)/u(1,2) + (u(2,2)/u(1,2) - u(2,1)/u(1,1))*(radius_tilde(0)-radius(1))/d_radius_tilde(1)          !velocity at inner radius

    IF ( ( v_in - cs ) <= ZERO) charak(1,1) = charak(1,2)
    IF (   v_in        <= ZERO) charak(2,1) = charak(2,2)
    IF ( ( v_in + cs ) <= ZERO) charak(3,1) = charak(3,2)


    v_out = u(2,N)/u(1,N) + (u(2,N)/u(1,N) - u(2,N-1)/u(1,N-1))*(radius_tilde(N)-radius(N))/d_radius_tilde(N-1)   !velocity at outer radius

    IF ( ( v_out - cs ) >= ZERO ) charak(1,N+1) = charak(1,N)
    IF (   v_out        >= ZERO ) charak(2,N+1) = charak(2,N)
    IF ( ( v_out + cs ) >= ZERO ) charak(3,N+1) = charak(3,N)


    u(1,0)   = ( radius(0)/radius(1)) * u(1,1) * EXP ( - (HALF/cs) * (charak(1,1) + charak(3,1)) )
    u(2,0)   = ( u(2,1)/u(1,1)             - HALF * (charak(3,1) - charak(1,1)) ) * u(1,0)
    u(3,0)   = ( u(3,1)/(u(1,1)*radius(1)) - charak(2,1) ) * radius(0) * u(1,0)

    u(1,N+1) = (radius(N+1)/radius(N)) * u(1,N) * EXP ( + (HALF/cs) * (charak(3,N+1) + charak(1,N+1)) )
    u(2,N+1) = ( u(2,N)/u(1,N)               + HALF * (charak(3,N+1) - charak(1,N+1))) * u(1,N+1)
    u(3,N+1) = ( (u(3,N)/(radius(N)*u(1,N))) + charak(2,N+1) )*radius(N+1)*u(1,N+1)



  CASE(2)     !reflecting boundary conditions

    u(1,0) =   (radius(0)/radius(1)) * u(1,1)
    u(2,0) = - (u(2,1)/u(1,1)) * u(1,0)
    u(3,0) =   (radius(0)/radius(1)) * (u(1,0)/u(1,1)) * u(3,1)

    u(1,N+1) =   (radius(N+1)/radius(N)) * u(1,N)
    u(2,N+1) = - (u(2,N)/u(1,N)) * u(1,N+1)
    u(3,N+1) =   (radius(N+1)/radius(N)) * (u(1,N+1)/u(1,N)) * u(3,N)



  CASE(3)     !no-gradients boundary conditions

    u(1,0) = (radius(0)/radius(1)) * u(1,1)
    u(2,0) = (u(2,1)/u(1,1)) * u(1,0)
    u(3,0) = (radius(0)/radius(1)) * (u(1,0)/u(1,1)) * u(3,1)

    u(1,N+1) = (radius(N+1)/radius(N)) * u(N,1)
    u(2,N+1) = (u(2,N)/u(1,N)) * u(1,N+1)
    u(3,N+1) = (radius(N+1)/radius(N)) * (u(1,N+1)/u(1,N)) * u(3,N)



  CASE (4)     !boundarys for the pringle-model

    u(1,0)   = 1.D-10 * u(1,1)
    u(1,N+1) = 1.D-10 * u(1,N)

    u(2,0)   = 1.D-10 * u(2,1)
    u(2,N+1) = 1.D-10 * u(2,N)



  CASE (5)     !boundarys for the pringle-disk

    tau     = (THREE*beta/(FOUR*s_0))*time

    u(1,0)   = radius(0)  *(sigma_0/(TWO*(SQRT(radius(0)/s_0)**3)*SQRT(c_pi*tau))) &
                 * (EXP(-(SQRT(radius(0)/s_0)-1)**2/(FOUR*tau)) - EXP(-(SQRT(radius(0)/s_0)+1)**2/(FOUR*tau)))
    u(1,N+1) = radius(N+1)  *(sigma_0/(TWO*(SQRT(radius(N+1)/s_0)**3)*SQRT(c_pi*tau))) &
                 * (EXP(-(SQRT(radius(N+1)/s_0)-1)**2/(FOUR*tau)) - EXP(-(SQRT(radius(N+1)/s_0)+1)**2/(FOUR*tau)))

    xi = SQRT(radius(0)/s_0)
    u(2,0)   = - ((THREE*beta*s_0*sigma_0)/(EIGHT*SQRT(c_pi*(tau**3)))) &
                  * (  (xi+ONE)*EXP(- ((xi+ONE)**2/(FOUR*tau))) - (xi-ONE)*EXP(- ((xi-ONE)**2/(FOUR*tau))))

    xi = SQRT(radius(N+1)/s_0)
    u(2,N+1) = - ((THREE*beta*s_0*sigma_0)/(EIGHT*SQRT(c_pi*(tau**3)))) &
                  * (  (xi+ONE)*EXP(- ((xi+ONE)**2/(FOUR*tau))) - (xi-ONE)*EXP(- ((xi-ONE)**2/(FOUR*tau))))

END SELECT



!----------------------------------------------------------------------------------------------------

END SUBROUTINE calc_boundarys
