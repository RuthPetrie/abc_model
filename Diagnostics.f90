! This file contains

! SUBROUTINE Effective_buoyancy
! SUBROUTINE Energy
! SUBROUTINE Calc_hydro
! SUBROUTINE Calc_geost
! SUBROUTINE Lscale_from_fft
! SUBROUTINE Magnitude_rms

! R. Petrie, 3.0 30-7-2013

!===================================================================================
 SUBROUTINE Effective_buoyancy (state)

!**************************************************
!* Subroutine to calculate the effective buoyancy *
!**************************************************

USE DefConsTypes, ONLY :   &
    model_vars_type,       &
    nlongs, nlevs,         &
    A, dz

IMPLICIT NONE

TYPE(model_vars_type),  INTENT(INOUT)   :: state

! Declare local variables
!-------------------
INTEGER                     :: i, k

! Routine written assuming constant dz
! modification required if using Charney-Philips

! NOT NORMALIZED

DO k = 1, nlevs
  DO i = 1, nlongs

   state % b_ef(i,k) = ( ( state % b(i,k+1) - state % b(i,k-1) ) / dz ) + A*A

  END DO
END DO

END SUBROUTINE Effective_buoyancy
!===================================================================================


!===================================================================================
SUBROUTINE Energy (state)

! To calculate components of energy

USE DefConsTypes, ONLY :  &
    model_vars_type,      &
    ZREAL8,               &
    nlongs, nlevs,        &
    A, B, C

IMPLICIT NONE

TYPE(model_vars_type),  INTENT(INOUT)   :: state
INTEGER                                 :: x, z
REAL(ZREAL8)                            :: u, v, w, rho, bp, r
REAL(ZREAL8)                            :: E_k, E_b, E_e


! Initialise
E_k = 0.0
E_e = 0.0
E_b = 0.0

! Calculate energy
DO z = 1, nlevs
  DO x = 1, nlongs
    r    = state % r(x,z)
    u    = state % u(x,z)
    v    = state % v(x,z)
    w    = state % w(x,z)
    bp   = state % b(x,z)
    rho  = state % rho(x,z)

    E_k  = E_k + rho * (u*u + v*v + w*w) / 2.0
    E_e  = E_e + C * r*r / (2. * B)
    E_b  = E_b + rho * bp*bp / (2. * A*A)
  END DO
END DO

state % Kinetic_Energy = E_k
state % Buoyant_Energy = E_b
state % Elastic_Energy = E_e
state % Total_Energy   = E_k + E_e + E_b

END SUBROUTINE Energy



!===================================================================================
SUBROUTINE Calc_geost (state, Dims)

! Calculate geostrophic imbalance in real space
USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  model_vars_type,        &
  Dimensions_type,        &
  nlongs, nlevs,          &
  f, C, recipdx

IMPLICIT NONE

! Declare parameters
!-------------------
TYPE(model_vars_type), INTENT(INOUT)  :: state
TYPE(Dimensions_type), INTENT(IN)     :: Dims

! Declare local variables
! -----------------------
INTEGER                               :: x, z
REAL(ZREAL8)                          :: rms1, rms2, norm
REAL(ZREAL8)                          :: term1(1:nlongs,1:nlevs), term2(1:nlongs,1:nlevs)

! Functions
! ---------
REAL(ZREAL8)                          :: RMS

! Calculate each term in the geostrophic imbalance equation
DO z = 1, nlevs
  DO x = 1, nlongs
    term1(x,z) = C * (state % r(x,z) - state % r(x-1,z)) * recipdx
    term2(x,z) = -1. * f * (state % v(x,z) + state % v(x-1,z)) / 2.
  END DO
END DO

! Find the RMS of each
rms1 = RMS(term1(1:nlongs,1:nlevs))
rms2 = RMS(term2(1:nlongs,1:nlevs))
norm = 1. / (rms1 + rms2)

! Compute the normalised geostrophic imbalance diagnostic
state % geost_imbal(1:nlongs,1:nlevs) = (term1(1:nlongs,1:nlevs) + term2(1:nlongs,1:nlevs)) * norm


END SUBROUTINE Calc_geost
!===================================================================================



!===================================================================================
SUBROUTINE Calc_hydro (state, Dims)

! Calculate hydrostatic imbalance in real space

USE DefConsTypes, ONLY :   &
    model_vars_type,       &
    Dimensions_type,       &
    ZREAL8,                &
    nlongs, nlevs, C

IMPLICIT NONE

! Declare parameters
!---------------------
TYPE(model_vars_type), INTENT(INOUT)  :: state
TYPE(Dimensions_type), INTENT(IN)     :: Dims

! Declare variables
!---------------------
INTEGER                               :: x, z
REAL(ZREAL8)                          :: rms1, rms2, norm
REAL(ZREAL8)                          :: term1(1:nlongs,1:nlevs), term2(1:nlongs,1:nlevs)

! Functions
! ---------
REAL(ZREAL8)                          :: RMS

! Calculate each term in the hydrostatic imbalance equation
DO x = 1, nlongs
  DO z = 1, nlevs
    term1(x,z) = C * (state % r(x,z+1) - state % r(x,z)) /         &
                     (Dims % half_levs(z+1) - Dims % half_levs(z))
    term2(x,z) = -1. * state % b(x,z)
  ENDDO
ENDDO

! Find the RMS of each
rms1 = RMS(term1(1:nlongs,1:nlevs))
rms2 = RMS(term2(1:nlongs,1:nlevs))
norm = 1. / (rms1 + rms2)

! Compute the normalised geostrophic imbalance diagnostic
state % hydro_imbal(1:nlongs,1:nlevs) = (term1(1:nlongs,1:nlevs) + term2(1:nlongs,1:nlevs)) * norm

END SUBROUTINE Calc_hydro
!===================================================================================



!===================================================================================
SUBROUTINE Calc_vert_mom_source (state, Dims)

! Calculate source of vertical momentum in real space

USE DefConsTypes, ONLY :   &
    model_vars_type,       &
    Dimensions_type,       &
    ZREAL8,                &
    nlongs, nlevs, C

IMPLICIT NONE

! Declare parameters
!---------------------
TYPE(model_vars_type), INTENT(INOUT)  :: state
TYPE(Dimensions_type), INTENT(IN)     :: Dims

! Declare variables
!---------------------
INTEGER                               :: x, z

! This has a very similiar form to the hydrostatic imbalance

DO x = 1, nlongs
  DO z = 1, nlevs
    state % vert_mom_source(x,z) = ( state % b(x,z) -                                  &
                                     C * (state % r(x,z+1) - state % r(x,z)) /         &
                                     (Dims % half_levs(z+1) - Dims % half_levs(z)) ) *                                &
                                   state % rho(x,z)
  ENDDO
ENDDO


END SUBROUTINE Calc_vert_mom_source
!===================================================================================



!===================================================================================
SUBROUTINE Lscales_from_fft (dataset, nxpoints, nypoints, xlength, ylength)

! Code to estimate the lengthscale of a field of data by the FFT

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  wp

USE nag_fft
USE nag_sym_fft

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: dataset(1:nxpoints,1:nypoints)
INTEGER,      INTENT(IN)    :: nxpoints, nypoints
REAL(ZREAL8), INTENT(INOUT) :: xlength, ylength

COMPLEX(wp)                 :: fftx(1:nxpoints), ffty(1:nypoints)
REAL(wp)                    :: trigx(1:2*nxpoints), trigy(1:2*nypoints)
REAL(wp)                    :: Dummy_x(1:nxpoints), Dummy_y(1:nypoints)
REAL(ZREAL8)                :: mean_len_x(1:nypoints), mean_len_y(1:nxpoints), mean_len
REAL(ZREAL8)                :: this_len, nxpoints_r, nypoints_r, norm, weight, modulus
INTEGER                     :: k, xx, yy


nxpoints_r = REAL(nxpoints)
nypoints_r = REAL(nypoints)


! Deal with the x-direction
! -------------------------

CALL nag_fft_trig(trigx)

DO yy = 1, nypoints
  ! Do fft in x for this y value
  Dummy_x(1:nxpoints) = dataset(1:nxpoints,yy)
  fftx(1:nxpoints) = nag_fft_1d_real(Dummy_x(1:nxpoints), inverse=.false., trig=trigx)

  ! Treat fft^2 as a PDF - find the mean lengthscale
  mean_len = 0.0
  norm     = 0.0
  DO k = 1, nxpoints/2
    this_len = nxpoints_r / REAL(k)
    modulus  = CABS(fftx(k))
    weight   = modulus * modulus
    mean_len = mean_len + weight * this_len
    norm     = norm + weight
  END DO

  IF (norm < 0.0000000001) PRINT *, 'WARNING FROM HORIZONTAL FFT ROUTINE - ZERO AT ', yy
  mean_len_x(yy) = mean_len / (2.0 * norm)    ! Half wavelength
END DO

xlength = SUM(mean_len_x(1:nypoints)) / nypoints_r


! Deal with the y-direction
! -------------------------

CALL nag_fft_trig(trigy)

DO xx = 1, nxpoints
  ! Do fft in y for this x value
  Dummy_y(1:nypoints) = dataset(xx,1:nypoints)
  ffty(1:nypoints) = nag_fft_1d_real(Dummy_y(1:nypoints), inverse=.false., trig=trigy)

  ! Treat fft^2 as a PDF - find the mean lengthscale
  mean_len = 0.0
  norm     = 0.0
  DO k = 1, nypoints/2
    this_len = nypoints_r / REAL(k)
    modulus  = CABS(ffty(k))
    weight   = modulus * modulus
    mean_len = mean_len + weight * this_len
    norm     = norm + weight
  END DO

  IF (norm < 0.0000000001) PRINT *, 'WARNING FROM VERTICAL  FFT ROUTINE - ZERO AT ', xx
  mean_len_y(xx) = mean_len / (2.0 * norm)    ! Half wavelength
END DO

ylength = SUM(mean_len_y(1:nxpoints)) / nxpoints_r

END SUBROUTINE Lscales_from_fft
!===================================================================================


!===================================================================================
SUBROUTINE Magnitude_rms (dataset, nxpoints, nypoints, rms)

! Code to calculate the RMS of a field

USE DefConsTypes, ONLY :  &
  ZREAL8

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: dataset(1:nxpoints,1:nypoints)
INTEGER,      INTENT(IN)    :: nxpoints, nypoints
REAL(ZREAL8), INTENT(INOUT) :: rms

INTEGER                     :: xx, yy
REAL(ZREAL8)                :: nxpoints_r, nypoints_r, ms, field


nxpoints_r = REAL(nxpoints)
nypoints_r = REAL(nypoints)

ms = 0.0
DO yy = 1, nypoints
  DO xx = 1, nxpoints
    field = dataset(xx,yy)
    ms    = ms + field * field
  END DO
END DO

ms  = ms / (nxpoints_r * nypoints_r)

rms = SQRT(ms)


END SUBROUTINE Magnitude_rms
!===================================================================================


