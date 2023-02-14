DOUBLE PRECISION FUNCTION minmod(a,b)

USE constants

REAL (KIND = KIND(0.0d0)) :: a,b

!----------------------------------------



IF (ABS(a) < ABS(b)) THEN
  minmod = a
ELSE
  minmod = b
END IF

IF (a*b < 0.d0) minmod = 0.d0



!----------------------------------------

END FUNCTION minmod 
