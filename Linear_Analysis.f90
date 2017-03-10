SUBROUTINE Linear_Analysis (out_dir, exp_id)

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  wp,                     &
  nlongs,                 &
  nlevs,                  &
  Lx,                     &
  pi,                     &
  H,                      &
  f,                      &
  A,                      &
  B,                      &
  C

USE nag_sym_lin_sys

USE nag_sym_eig,  ONLY :  &
  nag_sym_eig_all

IMPLICIT NONE

!*****************************************************
!*   Program to calculate the eigenvalues and        *
!*   speeds fo the waves in the model and to         *
!*   detect the sensitivity of these to model        *
!*   parameters                                      *
!*                                                   *
!*   R. Petrie, R. Bannister                         *
!*   Version 3                                       *
!*   17-05-14                                        *
!*   Now using time units with dimensions            *
!*   (not scaled by f) 08-12-16                      *
!*****************************************************

! DECLARE PARAMETERS
!-------------------
CHARACTER(*),   INTENT(IN)  :: out_dir, exp_id


! DECLARE VARIABLES
!------------------
REAL(wp)                  :: inv_L, inv_H, sqrtBC, kfact, mfact
REAL(wp)                  :: mall(0:nlevs-1), kall(0:nlongs-1)
REAL(wp)                  :: L(1:5, 1:5)
REAL(wp)                  :: lambda(1:5)
REAL(wp)                  :: grav_frequency (0:nlongs-1, 0:nlevs-1)
REAL(wp)                  :: acou_frequency (0:nlongs-1, 0:nlevs-1)
REAL(wp)                  :: hori_grav_speed(0:nlongs-1, 0:nlevs-1)
REAL(wp)                  :: hori_acou_speed(0:nlongs-1, 0:nlevs-1)
REAL(wp)                  :: vert_grav_speed(0:nlongs-1, 0:nlevs-1)
REAL(wp)                  :: vert_acou_speed(0:nlongs-1, 0:nlevs-1)
CHARACTER                 :: uplo*1
INTEGER                   :: k, m, i, j


PRINT *, 'Inside linear analysis routine.'
PRINT *, 'Values for A, B, C'

PRINT *, 'A = ', A
PRINT *, 'B = ', B
PRINT *, 'C = ', C
PRINT *, 'f = ', f

inv_L   = 1.0 / Lx
inv_H   = 1.0 / H
sqrtBC  = SQRT(B * C)
kfact   = Lx / (2. * pi)
mfact   = H / (2. * pi)
uplo    = 'u'


! Allocate horizontal and vertical wavenumbers
!---------------------------------------------
DO k = 0, nlongs-1
  kall(k) = 2. * REAL(k) * pi
END DO

DO m = 0, nlevs-1
  mall(m) = 2. * REAL(m) * pi
END DO



! Elements of matrix are computed below
DO m = 0, nlevs-1

  DO k = 0, nlongs-1

    L(1,1) = 0.0
    L(1,2) = f
    L(1,3) = 0.0
    L(1,4) = -1. * kall(k) * sqrtBC * inv_L
    L(1,5) = 0.0

    L(2,1) = L(1,2)
    L(2,2) = 0.0
    L(2,3) = 0.0
    L(2,4) = 0.0
    L(2,5) = 0.0

    L(3,1) = 0.0
    L(3,2) = 0.0
    L(3,3) = 0.0
    L(3,4) = mall(m) * sqrtBC * inv_H
    L(3,5) = A

    L(4,1) = L(1,4)
    L(4,2) = 0.0
    L(4,3) = L(3,4)
    L(4,4) = 0.0
    L(4,5) = 0.0

    L(5,1) = 0.0
    L(5,2) = 0.0
    L(5,3) = L(3,5)
    L(5,4) = 0.0
    L(5,5) = 0.0

    IF ( ((m == 0) .AND. (k == 0)) .OR. ((m == 3) .AND. (k == 7)) ) THEN
      PRINT *, 'Example matrix'
      DO i = 1, 5
        PRINT '(5F8.3)', (L(i,j), j=1,5)
      END DO
    END IF

    CALL nag_sym_eig_all (uplo, L(1:5,1:5), lambda(1:5))
    acou_frequency(k,m) = lambda(5)  ! nag routine lists frequencies in order - - 0 + +
    grav_frequency(k,m) = lambda(4)  ! nag routine lists frequencies in order - - 0 + +

  END DO
END DO


! Calculate wave speeds
!----------------------

DO m = 0, nlevs-1
  DO k = 1, nlongs-1
    hori_grav_speed(k,m) = (grav_frequency(k,m) - grav_frequency(k-1,m))  * kfact
    hori_acou_speed(k,m) = (acou_frequency(k,m) - acou_frequency(k-1,m))  * kfact
  ENDDO
  hori_grav_speed(0,m)   = 0.0
  hori_acou_speed(0,m)   = 0.0
ENDDO

DO k = 0, nlongs-1
  DO m = 1, nlevs-1
    vert_grav_speed(k,m) = (grav_frequency(k,m)  - grav_frequency(k,m-1))  * mfact
    vert_acou_speed(k,m) = (acou_frequency(k,m)  - acou_frequency(k,m-1))  * mfact
  ENDDO
  vert_grav_speed(k,0)   = 0.0
  vert_acou_speed(k,0)   = 0.0

ENDDO

! Dump data
!------------


OPEN (51, file = out_dir // '/grav_frequency_'  // exp_id // '.dat')
OPEN (52, file = out_dir // '/acou_frequency_'  // exp_id // '.dat')
OPEN (53, file = out_dir // '/hori_grav_speed_' // exp_id // '.dat')
OPEN (54, file = out_dir // '/hori_acou_speed_' // exp_id // '.dat')
OPEN (55, file = out_dir // '/vert_grav_speed_' // exp_id // '.dat')
OPEN (56, file = out_dir // '/vert_acou_speed_' // exp_id // '.dat')


DO k = 0, nlongs-1
  WRITE (51, *) grav_frequency(k, 0:nlevs-1)
  WRITE (52, *) acou_frequency(k, 0:nlevs-1)
ENDDO

DO m = 0, nlevs-1
  WRITE (53, *) hori_grav_speed(0:nlongs-1, m)
  WRITE (54, *) hori_acou_speed(0:nlongs-1, m)
ENDDO

DO k = 0, nlongs-1
  WRITE (55, *) vert_grav_speed(k, 0:nlevs-1)
  WRITE (56, *) vert_acou_speed(k, 0:nlevs-1)
ENDDO


CLOSE(51)
CLOSE(52)
CLOSE(53)
CLOSE(54)
CLOSE(55)
CLOSE(56)

END SUBROUTINE Linear_Analysis
