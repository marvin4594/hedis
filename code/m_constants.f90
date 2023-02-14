MODULE constants

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



INTEGER,           PARAMETER :: pcs = KIND(0.0d0)

REAL (KIND = pcs), PARAMETER :: c_pi   = 3.141592654d0

REAL (KIND = pcs), PARAMETER :: c_grav = 4.493568D-15          ![pc^3 / m_sun yr^2]
REAL (KIND = pcs), PARAMETER :: tau_s  = 4.50882D+8            ![yrs]

REAL (KIND = pcs), PARAMETER :: ZERO   = 0.0d0
REAL (KIND = pcs), PARAMETER :: HALF   = 0.5d0
REAL (KIND = pcs), PARAMETER :: ONE    = 1.0d0
REAL (KIND = pcs), PARAMETER :: TWO    = 2.0d0
REAL (KIND = pcs), PARAMETER :: THREE  = 3.0d0
REAL (KIND = pcs), PARAMETER :: FOUR   = 4.0d0
REAL (KIND = pcs), PARAMETER :: SIX    = 6.0d0
REAL (KIND = pcs), PARAMETER :: EIGHT  = 8.0d0
REAL (KIND = pcs), PARAMETER :: TEN    = 10.0d0

REAL (KIND = pcs), PARAMETER :: TWO_THIRD  = 2.0d0/3.0d0
REAL (KIND = pcs), PARAMETER :: FOUR_THIRD = 4.0d0/3.0d0



!----------------------------------------------------------------------------------------------------

END MODULE constants
