FUNCTION H2F1(z)
!computes the hypergeometrical function 2F1

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------

INTEGER :: n

REAL (KIND = KIND(0.0d0)), INTENT (IN)     :: z
REAL (KIND = KIND(0.0d0))                  :: a, b, c, R
REAL (KIND = KIND(0.0d0))                  :: H2F1

!----------------------------------------------------------------------------------------------------



n = 0
a = 0.5d0
b = -1.0d0
c = 0.69314718d0-0.25d0*LOG(1.0d0-z)
H2F1 = 0.0d0


IF (z > 0.50d0) THEN

  DO

    H2F1 = H2F1 + (a**2)*(b+c)*((1.0d0-z)**n)
    n = n + 1
    a = a*(1.0d0 + 1.0d0/(2.0d0*n))
    b = b + 1.0d0/(2.0d0*n) - 1.0d0/(2.0d0*n+1.0d0)
    R = ABS(((a**2)*b*((1.0d0-z)**n))/(1.0d0-((1.0d0+1.0d0/(2.0d0*(n+1.0d0)))**2)*(1.0d0-z)))*(ABS(b)+ABS(c))

    IF (R == 0.0d0) EXIT

  END DO

  H2F1 = (128.0d0/3.141592654d0)*H2F1

ELSE

  DO

    H2F1 = H2F1 + (a**2)*((n+1.0d0)/(n+2.0d0))*z**n
    R = a*(8.0d0*(n+2.0d0)/(n+3.0d0))*(z**(n+1.0d0)/(1.0d0-z))
    n = n + 1
    a = a*((n+0.5d0)/(n+1.0d0))

    IF (R == 0.0d0) EXIT

  END DO

  H2F1 = 8.0d0*H2F1

END IF



!----------------------------------------------------------------------------------------------------

END FUNCTION H2F1
