! This file contains
! SUBROUTINE Convection

! R. Petrie, 3.0: 30-7-13
!!==================================================================================================
SUBROUTINE Convection (state, x, z, xscale, zscale)

! Convection forcing at x, z
! with lengthscales xscale, zscale

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  model_vars_type,        &
  nlongs,                 &
  nlevs


IMPLICIT NONE

!Declare global parameters
!------------------
TYPE(model_vars_type), INTENT(INOUT) :: state
INTEGER,               INTENT(IN)    :: x, z, xscale, zscale
REAL(ZREAL8)                         :: temp1, temp2, zscale2, xscale2
REAL(ZREAL8)                         :: tempvalue
INTEGER                              :: i, k

zscale2 = zscale * zscale
xscale2 = xscale * xscale

DO i = 1, nlongs
  DO k = 1, nlevs
    temp1 = (x - i) * (x - i)
    temp2 = (z - k) * (z - k)
    tempvalue = 0.00001 * EXP (-(temp1/xscale2) - (temp2/zscale2)) * ABS(COS(REAL(4*i)))
    state % b(i,k) = state % b(i,k) + tempvalue
  END DO
END DO

END SUBROUTINE Convection
!!==================================================================================================
