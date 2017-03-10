SUBROUTINE BoundaryMod (param)

! Subroutine to modify the boundary such that they are appropriate for periodic conditions
USE DefConsTypes, ONLY :    &
    ZREAL8,                 &
    nlongs, nlevs,          &
    BoundSpread

IMPLICIT NONE

! Declare parameters
! ------------------
REAL(ZREAL8), INTENT(INOUT) :: param(1:nlongs, 1:nlevs)

! Declare other variables
! -----------------------
REAL(ZREAL8) :: typical_diff, discont
REAL(ZREAL8) :: sign_param(1:nlongs, 1:nlevs)
REAL(ZREAL8) :: temp, total, diff, recipBS2
INTEGER      :: x, z

recipBS2 = 1. / (BoundSpread * BoundSpread)
DO z = 1, nlevs
  total = 0.
  DO x = 1, nlongs-1
    total = total + ABS(param(x+1,z) - param(x,z))
  END DO
  typical_diff = total / REAL(nlongs-1)
  discont      = param(1,z) - param(nlongs,z)
  diff         = discont - typical_diff
  DO x = 1, nlongs/2 - 1
    temp       = -1. * REAL(x) * REAL(x) * recipBS2
    param(x,z) = param(x,z) - EXP(temp) * diff / 2.
  END DO
  DO x = nlongs/2, nlongs
    temp       = -1. * REAL(nlongs-x) * REAL(nlongs-x) * recipBS2
    param(x,z) = param(x,z) + EXP(temp) * diff / 2.
  END DO
END DO

END SUBROUTINE BoundaryMod
