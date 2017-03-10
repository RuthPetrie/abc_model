SUBROUTINE Boundaries (state, set_u, set_v, set_w, set_r, set_b, set_rho, set_beff, set_tracer)

!************************************
!* Subroutine to apply boundary     *
!* conditions to model variables    *
!*                                  *
!* Input flags to choose which      *
!* variables to apply boundary      *
!* conditions to                    *
!*                                  *
!* Set no boundaries to do them all *
!************************************

USE DefConsTypes, ONLY :     &
    model_vars_type,         &
    nlongs, nlevs

IMPLICIT NONE

TYPE(model_vars_type), INTENT(INOUT) :: state
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_u
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_v
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_w
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_r
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_b
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_rho
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_beff
LOGICAL, OPTIONAL,     INTENT(IN)    :: set_tracer

!Local variables
!---------------
LOGICAL                              :: flag_u, flag_v, flag_w, flag_r
LOGICAL                              :: flag_b, flag_rho, flag_beff, flag_tracer
LOGICAL                              :: some_flags

IF (PRESENT(set_u)) THEN
  flag_u = set_u
ELSE
  flag_u = .FALSE.
END IF
IF (PRESENT(set_v)) THEN
  flag_v = set_v
ELSE
  flag_v = .FALSE.
END IF
IF (PRESENT(set_w)) THEN
  flag_w = set_w
ELSE
  flag_w = .FALSE.
END IF
IF (PRESENT(set_r)) THEN
  flag_r = set_r
ELSE
  flag_r = .FALSE.
END IF
IF (PRESENT(set_b)) THEN
  flag_b = set_b
ELSE
  flag_b = .FALSE.
END IF
IF (PRESENT(set_rho)) THEN
  flag_rho = set_rho
ELSE
  flag_rho = .FALSE.
END IF
IF (PRESENT(set_beff)) THEN
  flag_beff = set_beff
ELSE
  flag_beff = .FALSE.
END IF
IF (PRESENT(set_tracer)) THEN
  flag_tracer = set_tracer
ELSE
  flag_tracer = .FALSE.
END IF


some_flags = PRESENT(set_u) .OR. PRESENT(set_v) .OR. PRESENT(set_w) .OR. PRESENT(set_r) .OR. &
             PRESENT(set_b) .OR. PRESENT(set_rho) .OR. PRESENT(set_beff) .OR. PRESENT(set_tracer)

IF (.NOT.some_flags) THEN
  ! No flags set.  This is shorthand for all flags
  flag_u      = .TRUE.
  flag_v      = .TRUE.
  flag_w      = .TRUE.
  flag_r      = .TRUE.
  flag_b      = .TRUE.
  flag_rho    = .TRUE.
  flag_beff   = .TRUE.
  flag_tracer = .TRUE.
END IF

!PRINT*, 'boundaries called'

IF (flag_u) THEN
  ! Horizontal boundaries
  state % u(0, 0:nlevs+1)        = state % u(nlongs, 0:nlevs+1)
  state % u(nlongs+1, 0:nlevs+1) = state % u(1, 0:nlevs+1)
  ! Vertical boundaries
  state % u(0:nlongs+1,0)        = - 1.0 * state % u(:,1)
  state % u(0:nlongs+1,nlevs+1)  = state % u(0:nlongs+1, nlevs)
ENDIF

IF (flag_v) THEN
  ! Horizontal boundaries
  state % v(0, 0:nlevs+1)        = state % v(nlongs, 0:nlevs+1)
  state % v(nlongs+1, 0:nlevs+1) = state % v(1, 0:nlevs+1)
  ! Vertical boundaries
  state % v(0:nlongs+1,0)        = - 1.0 * state % v(:,1)
  state % v(0:nlongs+1,nlevs+1)  = state % v(0:nlongs+1, nlevs)
ENDIF

IF (flag_w) THEN
  ! Horizontal boundaries
  state % w(0, 0:nlevs+1)        = state % w(nlongs, 0:nlevs+1)
  state % w(nlongs+1, 0:nlevs+1) = state % w(1, 0:nlevs+1)
  ! Vertical boundaries
  state % w(0:nlongs+1,0)        = 0.0
  state % w(0:nlongs+1,nlevs)    = 0.0
  state % w(0:nlongs+1,nlevs+1)  = 0.0
ENDIF

IF (flag_r) THEN
  ! Horizontal boundaries
  state % r(0, 0:nlevs+1)        = state % r(nlongs, 0:nlevs+1)
  state % r(nlongs+1, 0:nlevs+1) = state % r(1, 0:nlevs+1)
  ! Vertical boundaries
  state % r(0:nlongs+1,0)        = state % r(:,1)
  state % r(0:nlongs+1,nlevs+1)  = state % r(0:nlongs+1,nlevs)
ENDIF

IF (flag_b) THEN
  ! Horizontal boundaries
  state % b(0, 0:nlevs+1)        = state % b(nlongs, 0:nlevs+1)
  state % b(nlongs+1, 0:nlevs+1) = state % b(1, 0:nlevs+1)
  ! Vertical boundaries
  state % b(0:nlongs+1,0)        = 0.0
  state % b(0:nlongs+1,nlevs)    = 0.0
  state % b(0:nlongs+1,nlevs+1)  = 0.0
ENDIF

IF (flag_rho) THEN
  ! Horizontal boundaries
  state % rho(0, 0:nlevs+1)        = state % rho(nlongs, 0:nlevs+1)
  state % rho(nlongs+1, 0:nlevs+1) = state % rho(1, 0:nlevs+1)
  ! Vertical boundaries
  state % rho(0:nlongs+1,0)        = state % rho(:,1)
  state % rho(0:nlongs+1,nlevs+1)  = state % rho(0:nlongs+1,nlevs)
ENDIF

IF (flag_beff) THEN
  ! Horizontal boundaries
  state % b_ef(0, 0:nlevs+1)        = state % b_ef(nlongs, 0:nlevs+1)
  state % b_ef(nlongs+1, 0:nlevs+1) = state % b_ef(1, 0:nlevs+1)
  ! Vertical boundaries
  state % b_ef(0:nlongs+1,0)        = 0.0
  state % b_ef(0:nlongs+1,nlevs)    = 0.0
  state % b_ef(0:nlongs+1,nlevs+1)  = 0.0
ENDIF

IF (flag_tracer) THEN
  ! Horizontal boundaries
  state % tracer(0, 0:nlevs+1)        = state % tracer(nlongs, 0:nlevs+1)
  state % tracer(nlongs+1, 0:nlevs+1) = state % tracer(1, 0:nlevs+1)
  state % tracer(0:nlongs+1, 0)       = state % tracer(0:nlongs+1, 1)
  state % tracer(0:nlongs+1, nlevs+1) = state % tracer(0:nlongs+1, nlevs)
ENDIF

END SUBROUTINE Boundaries
!===================================================================================================
