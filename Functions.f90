!===================================================================================================
! In this file are functions and ancillary subroutines
!
! Functions
! ----------
! FUNCTION STDEV
! FUNCTION GAUSS
! FUNCTION RMS
! FUNCTION INT_FH (varl, varu, k, Dims)
! FUNCTION INT_HF (varl, varu, k, Dims)
! SUBROUTINE InterpolationTest

! R. Petrie, 3.0: 30-7-13

!===================================================================================================
 REAL(ZREAL8) FUNCTION STDEV (field)
!********************************************
!* Function to calculate standard deviation *
!* of a 2d field that has dimension         *
!* 1:nlongs, 1:nlevs                      *
!********************************************

USE DefConsTypes, ONLY  : &
  ZREAL8,                 &
  nlongs,                 &
  nlevs

IMPLICIT NONE

!Declare parameters
REAL(ZREAL8), INTENT(IN)          :: field(1:nlongs, 1:nlevs)
REAL(ZREAL8)                      :: ave, std, recip

recip = 1. / REAL(nlongs * nlevs)
ave   = recip * SUM(field(1:nlongs,1:nlevs))
std   = recip * SUM((field(1:nlongs,1:nlevs) - ave) * (field(1:nlongs,1:nlevs) - ave))

STDEV = sqrt(std)

END FUNCTION STDEV
!===================================================================================================



!===================================================================================================
REAL(ZREAL8) FUNCTION GAUSS (std)
!********************************************
!* Function to calculate a gaussian         *
!* distributed random variable using the    *
!* Box-Mueller algorithm where the          *
!* distribution has a variance of std^2     *
!********************************************

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  pi

IMPLICIT NONE

!Declare parameters
REAL(ZREAL8), INTENT(IN)         :: std
REAL                             :: RAND
REAL(ZREAL8)                     :: u1, u2

u1 = RAND(0)
u2 = RAND(0)

gauss = sqrt(-2*log(u1))*cos(2*pi*u2)*std;

END FUNCTION GAUSS
!===================================================================================================



!===================================================================================================
REAL(ZREAL8) FUNCTION RMS (state)

!********************************************
!* Function to calculate rms of a 2d-state  *
!* of dims (1:nlongs, 1:nlevs)              *
!********************************************

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  nlongs,                &
  nlevs


IMPLICIT NONE

! Declare parameters
!------------------
REAL(ZREAL8), INTENT(IN)    :: state(1:nlongs, 1:nlevs)

! Declare variables
!-------------------
REAL(ZREAL8)                :: eps, recip

eps   = 10E-19
recip = 1. / REAL(nlongs*nlevs)

! Calculate rmse
!-----------------
RMS = SQRT(recip * SUM(state(1:nlongs,1:nlevs) * state(1:nlongs,1:nlevs)))

IF (RMS .LT. eps) PRINT*, 'RMS SMALL'

END FUNCTION RMS
!==================================================================================================




!==================================================================================================
REAL(ZREAL8) FUNCTION INT_FH (varl, varu, z, Dims)
! Function interpolates variable from full to half levs
! at a given grid point
! varl at (i,k-1)
! varu at (i,k)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  Dimensions_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(IN) :: varl
REAL(ZREAL8),          INTENT(IN) :: varu
INTEGER,               INTENT(IN) :: z
TYPE(Dimensions_type), INTENT(IN) :: dims

! Calculate variable interpolated to full levs
INT_FH = Dims % a1(z) * varu + Dims % b1(z) * varl

END FUNCTION INT_FH
!!==================================================================================================



!!==================================================================================================
REAL(ZREAL8) FUNCTION INT_HF (varl, varu, z, Dims)
! Function interpolates variable from half to full levs
! at a given grid point
! varl at (i,k)
! varu at (i,k+1)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  Dimensions_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(IN) :: varl
REAL(ZREAL8),          INTENT(IN) :: varu
INTEGER,               INTENT(IN) :: z
TYPE(Dimensions_type), INTENT(IN) :: Dims

!Calculate variable interpolated to half levs
INT_HF = Dims % a2(z) * varu + Dims % b2(z) * varl

END FUNCTION INT_HF
!===================================================================================================



!==========================================================================
SUBROUTINE InterpolationTest(state, Dims)

! Subroutine to test interpolation functions

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  model_vars_type,       &
  Dimensions_type,       &
  nlongs,                &
  nlevs

IMPLICIT NONE

!Declare parameters
!------------------
TYPE(model_vars_type), INTENT(IN)   :: state
TYPE(Dimensions_type), INTENT(IN)   :: Dims

!Declare variables
!-----------------
INTEGER                 :: x,z
REAL(ZREAL8)            :: variable_test(1:nlongs, 1:nlevs)
REAL(ZREAL8)            :: INT_HF, INT_FH

PRINT*, ' File directory not known'
PRINT*, 'STOP'
STOP

! change file output locations

!Output original data
OPEN(51, file ='/home/wx019276/Modelling/Matlab/Data/varorig.dat')
DO z = 1, nlevs
  WRITE (51, *) state % w(1:nlongs,z)
ENDDO
CLOSE(51)

! Perform full to half lev interpolation
DO z = 1, nlevs
  DO x = 1, nlongs
   variable_test(x,z) =  INT_FH(state % w(x,z-1), state % w(x,z),z, Dims)
  ENDDO
ENDDO

! Perform half to full lev interpolation
DO z = 2, nlevs-1
  DO x = 1, nlongs
    variable_test(x,z) =  INT_HF(state % w(x,z), state % w(x,z+1),z, Dims)
  ENDDO
ENDDO

!Output interpolated data
OPEN (52, file ='/home/wx019276/Modelling/Matlab/Data/varnew.dat')
DO z = 1, nlevs-1
  WRITE (52, *) variable_test(1:nlongs, z)
ENDDO
CLOSE(52)

END SUBROUTINE InterpolationTest
!==========================================================================
