SUBROUTINE ModelDriver (StateIn, StateOut, Dims, ntimesteps, ndumps, fileout, &
                        diagnostics_file, convection_switch)

!***********************************************************************************
!*                                                                                 *
!*  Driver Routine for Non linear forward model                                    *
!*                                                                                 *
!*  StateIn                   - input state of model_vars_type                     *
!*  StateOut                  - output state of model_vars_type                    *
!*  Dims                      - dimension data                                     *
!*  ntimesteps                - total number of timesteps for model integration    *
!*  ndumps                    - number of output times                             *
!*  fileout                   - output file                                        *
!*  convection_switch         - convective forcing                                 *
!*                                                                                 *
!*                                                                                 *
!*   R. Petrie, 2.0:  10-6-2011                                                    *
!*   R. Petrie, 3.0:  30-7-2013                                                    *
!*                                                                                 *
!*                                                                                 *
!***********************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  model_vars_type,        &
  Dimensions_type,        &
  nlongs,                 &
  nlevs,                  &
  dt, deltat, dx, dz,     &
  Lengthscale_diagnostics,&
  A, B, C, f, g

IMPLICIT NONE

INCLUDE "Boundaries.interface"
!INCLUDE "/usr/include/netcdf.inc"

! Parameters
!-----------
TYPE(model_vars_type), INTENT(INOUT) :: StateIn
TYPE(model_vars_type), INTENT(INOUT) :: StateOut
TYPE(Dimensions_type), INTENT(INOUT) :: Dims
INTEGER,               INTENT(IN)    :: ntimesteps, ndumps
CHARACTER(*),          INTENT(IN)    :: fileout
CHARACTER(*),          INTENT(IN)    :: diagnostics_file
LOGICAL,               INTENT(IN)    :: convection_switch

! Local Variables
!----------------
INTEGER                              :: dumpspacing
INTEGER                              :: t, diagnostics_unit
REAL(ZREAL8)                         :: L_horiz_u, L_horiz_v, L_horiz_w, L_horiz_r, L_horiz_b, L_horiz_uv, L_horiz_rho
REAL(ZREAL8)                         :: L_vert_u, L_vert_v, L_vert_w, L_vert_r, L_vert_b, L_vert_uv, L_vert_rho
REAL(ZREAL8)                         :: rms_u, rms_v, rms_w, rms_r, rms_b, rms_uv, rms_rho, au
REAL(ZREAL8)                         :: ca, cg, Ro, M, Fr, RrA, Rrg, ar, sr, rr


! Functions
! ---------
REAL(ZREAL8)                         :: RMS

!PRINT*, 'Model driver called'

! CONSISTENCY CHECKS
!-------------------

dumpspacing = ntimesteps / ndumps

IF ( (dt/deltat /= 2.0 )) THEN
  PRINT*, '******************************'
  PRINT*, '********** Error *************'
  PRINT*, '  dt /deltat .NE. 2'
  STOP
ENDIF

PRINT*, 'Total no. of timesteps ',ntimesteps
PRINT*, 'Total run length ', dt*ntimesteps,' s'
PRINT*, 'Total run length ', dt*ntimesteps/3600,' hrs'
PRINT*, 'Number of output times ', ndumps
PRINT*, 'Output times every ', dt*ntimesteps/(3600*ndumps),' hrs'


! Diagnostics output file
OPEN (diagnostics_unit, file=diagnostics_file)


! Apply boundary conditions
!--------------------------
CALL Boundaries (StateIn, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_r=.TRUE., &
                          set_b=.TRUE., set_rho=.TRUE., set_tracer=.TRUE.)

! Calculate some diagnostics and write Initial Conditions
!-------------------------
CALL Calc_hydro(StateIn, Dims)
CALL Calc_geost(StateIn, Dims)
CALL Calc_vert_mom_source(StateIn, Dims)
CALL Energy(StateIn)
CALL Effective_buoyancy(StateIn)
CALL Boundaries (StateIn, set_beff=.TRUE.)
CALL Write_state_2d (fileout, StateIn, Dims, ndumps+1, 0, dumpspacing)
PRINT*, 'Output inital conditions:'
PRINT*, 'xconv -i ',fileout,' &'

PRINT*,'----------------------------'
PRINT*,'    Running prognostic model'
PRINT*,'----------------------------'


DO t = 1, ntimesteps
 ! PRINT*,' > ', t, ' of ', ntimesteps

  !-------------------------------------------
  IF (convection_switch) THEN
    IF (t <= 1350) THEN
      CALL Convection(StateIn, 180, 1, 50, 5 )
    ENDIF
  END IF
  !-------------------------------------------



  !*******************************************
  CALL Forward_model (StateIn, StateOut, Dims)
  !*******************************************

  !******************************************
  ! DATA DUMP
  !******************************************

  ! Calculate the energy contributions
  CALL Energy (StateOut)

  ! At output times
  !-----------------
  IF ((t/dumpspacing)*dumpspacing == t) THEN
    ! Calculate diagnostics
    !------------------------
    CALL Calc_hydro (StateOut, Dims)
    CALL Calc_geost (StateOut, Dims)
    CALL Calc_vert_mom_source(StateOut, Dims)
    CALL Effective_buoyancy (StateOut)
    CALL Write_state_2d (fileout, StateOut, Dims, ndumps+1, t/dumpspacing, dumpspacing)
  ENDIF

  WRITE (diagnostics_unit, *) t, t*dt, StateOut % Kinetic_Energy,      &
                                       StateOut % Buoyant_Energy,      &
                                       StateOut % Elastic_Energy,      &
                                       StateOut % Total_Energy

  ! Lengthscale diagnostics at the last timestep
  IF (Lengthscale_diagnostics .AND. (t == ntimesteps)) THEN
    WRITE (diagnostics_unit, *) '# -----------------------------------------'
    WRITE (diagnostics_unit, *) '# Diagnostics for last step (SI units)'
    WRITE (diagnostics_unit, *) '# -----------------------------------------'

    WRITE (*,*) 'Diagnosing lengthscales for u ...'
    CALL Lscales_from_fft(StateOut % u(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_u, L_vert_u)
    L_horiz_u = L_horiz_u * dx
    L_vert_u  = L_vert_u * dz
    WRITE (*,*) 'Diagnosing RMS for u ...'
    CALL Magnitude_rms(StateOut % u(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_u)

    WRITE (*,*) 'Diagnosing lengthscales for v ...'
    CALL Lscales_from_fft(StateOut % v(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_v, L_vert_v)
    L_horiz_v = L_horiz_v * dx
    L_vert_v  = L_vert_v * dz
    WRITE (*,*) 'Diagnosing RMS for v ...'
    CALL Magnitude_rms(StateOut % v(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_v)

    WRITE (*,*) 'Diagnosing lengthscales for w ...'
    CALL Lscales_from_fft(StateOut % w(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_w, L_vert_w)
    L_horiz_w = L_horiz_w * dx
    L_vert_w  = L_vert_w * dz
    WRITE (*,*) 'Diagnosing RMS for w ...'
    CALL Magnitude_rms(StateOut % w(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_w)

    WRITE (*,*) 'Diagnosing lengthscales for r ...'
    CALL Lscales_from_fft(StateOut % r(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_r, L_vert_r)
    L_horiz_r = L_horiz_r * dx
    L_vert_r  = L_vert_r * dz
    WRITE (*,*) 'Diagnosing RMS for r ...'
    CALL Magnitude_rms(StateOut % r(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_r)

    WRITE (*,*) 'Diagnosing lengthscales for rho ...'
    CALL Lscales_from_fft(StateOut % rho(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_rho, L_vert_rho)
    L_horiz_rho = L_horiz_rho * dx
    L_vert_rho  = L_vert_rho * dz
    WRITE (*,*) 'Diagnosing RMS for rho ...'
    CALL Magnitude_rms(StateOut % rho(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_rho)


    WRITE (*,*) 'Diagnosing lengthscales for b ...'
    CALL Lscales_from_fft(StateOut % b(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_b, L_vert_b)
    L_horiz_b = L_horiz_b * dx
    L_vert_b  = L_vert_b * dz
    WRITE (*,*) 'Diagnosing RMS for b ...'
    CALL Magnitude_rms(StateOut % b(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_b)

    ! Average data for horizontal wind
    L_horiz_uv = (L_horiz_u + L_horiz_v) / 2.0
    L_vert_uv  = (L_vert_u + L_vert_v) / 2.0
    rms_uv     = (rms_u + rms_v) / 2.0


    ! Acoustic wave speed
    ca  = SQRT(B*C)
    ! Gravity wave speed
    cg  = SQRT(g*L_vert_r)
    ! Rossby number
    Ro  = rms_u / (f * L_horiz_r)
    ! Mach number
    M   = rms_u / ca
    ! Froude number
    Fr  = rms_u / cg
    ! Rossby radius calculated with static stability
    RrA = A * L_vert_r / f
    ! Rossby radius calculated with gravity wave speed
    Rrg = cg / f
    ! Aspect ratio
    ar  = L_vert_u / L_horiz_u
    ! Speed ratio
    sr  = rms_w / rms_u
    ! Ratio of ratios
    rr  = ar / sr
    ! Ratio of horizontal wind speeds
    au = rms_v / rms_u

    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'A,                         ', A
    WRITE (diagnostics_unit, *) 'B,                         ', B
    WRITE (diagnostics_unit, *) 'C,                         ', C
    WRITE (diagnostics_unit, *) 'f,                         ', f
    WRITE (diagnostics_unit, *) 'g,                         ', g
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'u L_h,                     ', L_horiz_u
    WRITE (diagnostics_unit, *) 'u L_v,                     ', L_vert_u
    WRITE (diagnostics_unit, *) 'u RMS,                     ', rms_u
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'v L_h,                     ', L_horiz_v
    WRITE (diagnostics_unit, *) 'v L_v,                     ', L_vert_v
    WRITE (diagnostics_unit, *) 'v RMS,                     ', rms_v
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'hor-sp L_h,                ', L_horiz_uv
    WRITE (diagnostics_unit, *) 'hor-sp L_v,                ', L_vert_uv
    WRITE (diagnostics_unit, *) 'hor-sp RMS,                ', rms_uv
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'w L_h,                     ', L_horiz_w
    WRITE (diagnostics_unit, *) 'w L_v,                     ', L_vert_w
    WRITE (diagnostics_unit, *) 'w RMS,                     ', rms_w
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'r L_h,                     ', L_horiz_r
    WRITE (diagnostics_unit, *) 'r L_v,                     ', L_vert_r
    WRITE (diagnostics_unit, *) 'r RMS,                     ', rms_r
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'rho L_h,                   ', L_horiz_rho
    WRITE (diagnostics_unit, *) 'rho L_v,                   ', L_vert_rho
    WRITE (diagnostics_unit, *) 'rho RMS,                   ', rms_rho
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'b L_h,                     ', L_horiz_b
    WRITE (diagnostics_unit, *) 'b L_v,                     ', L_vert_b
    WRITE (diagnostics_unit, *) 'b RMS,                     ', rms_b
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'ca,                        ', ca
    WRITE (diagnostics_unit, *) 'cg,                        ', cg
    WRITE (diagnostics_unit, *) 'RossbyNo,                  ', Ro
    WRITE (diagnostics_unit, *) 'MachNo,                    ', M
    WRITE (diagnostics_unit, *) 'FroudeNo,                  ', Fr
    WRITE (diagnostics_unit, *) 'RbyRad (A),                ', RrA
    WRITE (diagnostics_unit, *) 'RbyRad (cg),               ', Rrg
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'Aspect ratio,              ', ar
    WRITE (diagnostics_unit, *) 'Speed ratio W/U,           ', sr
    WRITE (diagnostics_unit, *) 'Speed ratio V/U            ', au
    WRITE (diagnostics_unit, *) 'aspect ratio / speed ratio ', rr
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'MAGNITUDES OF TERMS IN EQUATIONS'
    WRITE (diagnostics_unit, *) 'For a key, see'
    WRITE (diagnostics_unit, *) 'www.met.reading.ac.uk/~ross/DARC/ConvScaleB/ABCmodel/ScaleAnalysis4Paper.html'
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'u1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'u2,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'u3,                        ', B * Ro * sr / ar
    WRITE (diagnostics_unit, *) 'u4,                        ', C * rms_r / (rms_u * f * L_horiz_r)
    WRITE (diagnostics_unit, *) 'u5,                        ', au
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'v1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'v2,                        ', B * Ro * L_horiz_u / L_horiz_v
    WRITE (diagnostics_unit, *) 'v3,                        ', B * Ro * L_horiz_u * sr / L_vert_v
    WRITE (diagnostics_unit, *) 'v4,                        ', 1.0 / au
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'w1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'w2,                        ', B * Ro * L_horiz_u / L_horiz_w
    WRITE (diagnostics_unit, *) 'w3,                        ', B * Ro * L_horiz_u * sr / L_vert_w
    WRITE (diagnostics_unit, *) 'w4,                        ', C * rms_r / (rms_w * f * L_vert_r)
    WRITE (diagnostics_unit, *) 'w5,                        ', rms_b / (rms_w * f)
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'rhop1,                     ', 1.0
    WRITE (diagnostics_unit, *) 'rho02,                     ', 1.0
    WRITE (diagnostics_unit, *) 'rhop3,                     ', L_horiz_u * sr / L_vert_w
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'bp1,                       ', B * Ro
    WRITE (diagnostics_unit, *) 'bp2,                       ', B * Ro * L_horiz_u / L_horiz_b
    WRITE (diagnostics_unit, *) 'bp3,                       ', B * Ro * L_horiz_u * sr / L_vert_b
    WRITE (diagnostics_unit, *) 'bp4,                       ', A * A * rms_w / (rms_b * f)
    WRITE (diagnostics_unit, *) '---'

    WRITE (*,*) 'Done diagnostics'
  END IF

  ! Re-assign
  StateIn = StateOut

END DO

CLOSE (diagnostics_unit)


PRINT*,'-----------------------------------------'
PRINT*,'  Output file'
PRINT*,'  xconv -i ', fileout, '&'
PRINT*,'-----------------------------------------'

END SUBROUTINE ModelDriver
!===================================================================================================




!===================================================================================================
SUBROUTINE Forward_model (input, stateout, Dims)

! Forward model integrated over dt
! adjust is the input state

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  model_vars_type,        &
  Dimensions_type,        &
  Averages_type,          &
  nlongs, nlevs,          &
  bdiva_f, deltat, recip_alpha_f, recipdx, f, recip2dx, bdiva_N, half, &
  recip_alpha_N, A, B, C, dt, &
  Adv_tracer

IMPLICIT NONE

INCLUDE "Boundaries.interface"

! Declare global parameters
!--------------------------
TYPE(model_vars_type), INTENT(IN)    :: input
TYPE(model_vars_type), INTENT(INOUT) :: stateout
TYPE(Dimensions_type), INTENT(IN)    :: Dims

! Declare local variables
!------------------------
TYPE(model_vars_type)                :: adjust, adjust_p1
TYPE(Averages_type)                  :: aves
REAL(ZREAL8)                         :: u_at_b(1:nlongs,1:nlevs), w_at_u(1:nlongs,1:nlevs)
REAL(ZREAL8)                         :: u_at_v(1:nlongs,1:nlevs), w_at_v(1:nlongs,1:nlevs)
REAL(ZREAL8)                         :: upstr_u, upstr_w, rhodivu
REAL(ZREAL8)                         :: deltat_recip_alpha_f, half_f, half_deltat_C
REAL(ZREAL8)                         :: deltat_recip_alpha_N, deltat_recip_alpha_N_N2, deltat_B
REAL(ZREAL8)                         :: deltat_recip_alpha_f_f, B_dt
REAL(ZREAL8)                         :: drdx1, drdx2, drdz, Nsq

! Declare functions
!------------------
REAL(ZREAL8)    :: INT_FH, INT_HF

! Declare Flags
!--------------
LOGICAL                              :: flag_u, flag_w, flag_u_at_b, flag_w_at_b
LOGICAL                              :: flag_w_at_u, flag_u_at_v, flag_w_at_v

! Declare Loop counters
!----------------------
INTEGER                              :: i, k, tsmall

! Initialise
CALL Initialise_Averages (aves)
Nsq                     = A * A
deltat_recip_alpha_f    = deltat * recip_alpha_f
half_f                  = half * f
half_deltat_C           = half * deltat * C
deltat_recip_alpha_N    = deltat * recip_alpha_N
deltat_recip_alpha_N_N2 = deltat * recip_alpha_N * Nsq
deltat_B                = deltat * B
deltat_recip_alpha_f_f  = deltat * recip_alpha_f * f
B_dt                    = B * dt

! copy input state into adjust array
adjust = input

!****************************
!***** ADJUSTMENT STAGE *****
!****************************

!Begin temporal loop, with small timestep
DO tsmall = 1, 2

!*** Forward-backward scheme ***

  DO k = 1, nlevs
    DO i = 1, nlongs
      drdx1              = ( adjust % r(i+1,k) - adjust % r(i,k) ) * recipdx
      drdx2              = ( adjust % r(i+1,k) - adjust % r(i-1,k) ) * recip2dx

      adjust_p1 % u(i,k) = bdiva_f * adjust % u(i,k) +                                                &
                           deltat_recip_alpha_f * ( half_f * (adjust % v(i,k) + adjust % v(i+1,k) ) - &
                                                    C * drdx1 )

      adjust_p1 % v(i,k) = bdiva_f * adjust % v(i,k) +                                                &
                           deltat_recip_alpha_f_f * ( half_deltat_C * drdx2 -                        &
                                                      half * ( adjust % u(i-1,k) + adjust % u(i,k) ) )

      drdz               = ( adjust % r(i,k+1) - adjust % r(i,k) ) * dims % recip_half_kp1_k(k)

      adjust_p1 % w(i,k) = bdiva_N * adjust % w(i,k) +                                                &
                           deltat_recip_alpha_N * ( adjust % b(i,k) - C * drdz )

      adjust_p1 % b(i,k) = bdiva_N * adjust % b(i,k) +                                                &
                           deltat_recip_alpha_N_N2 * ( half_deltat_C * drdz - adjust % w(i,k) )
    END DO
  END DO

  ! Reassign variables
  !-------------------
  adjust % u(1:nlongs,1:nlevs) = adjust_p1 % u(1:nlongs,1:nlevs)
  adjust % v(1:nlongs,1:nlevs) = adjust_p1 % v(1:nlongs,1:nlevs)
  adjust % w(1:nlongs,1:nlevs) = adjust_p1 % w(1:nlongs,1:nlevs)
  adjust % b(1:nlongs,1:nlevs) = adjust_p1 % b(1:nlongs,1:nlevs)
  adjust % w(1:nlongs,0)       = adjust_p1 % w(1:nlongs,0)

  CALL Boundaries (adjust, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_b=.TRUE.)

  !Store all adjustment variables
  !------------------------------
  IF (tsmall == 1) THEN
    aves % u_1(0:nlongs+1,0:nlevs+1) = adjust % u(0:nlongs+1,0:nlevs+1)
    aves % w_1(0:nlongs+1,0:nlevs+1) = adjust % w(0:nlongs+1,0:nlevs+1)
  ELSE IF (tsmall == 2 ) THEN
    aves % u_2(0:nlongs+1,0:nlevs+1) = adjust % u(0:nlongs+1,0:nlevs+1)
    aves % w_2(0:nlongs+1,0:nlevs+1) = adjust % w(0:nlongs+1,0:nlevs+1)
  END IF


!!!!!  adjust % rho(1:nlongs,1:nlevs) = 1.0 + adjust % r(1:nlongs,1:nlevs)
!!!!!
!!!!!  !Apply boundary conditions to full rho
!!!!!  CALL Boundaries (adjust, set_rho=.TRUE.)

  !Calculate rnext using a forward upstream
  !----------------------------------------

  !Calculate advecting interpolations
  DO k = 1, nlevs
    DO i = 1, nlongs
      u_at_v(i,k) = half * ( adjust % u(i-1, k) + adjust % u(i,k) )
      w_at_v(i,k) = INT_FH( adjust % w(i,k-1), adjust % w(i,k), k , Dims )
    END DO
  END DO

  DO k = 1, nlevs
    DO i = 1,nlongs

      !Calculate flags
      flag_u_at_v = ( u_at_v(i,k) .LE. 0.0 )
      flag_w_at_v = ( w_at_v(i,k) .LE. 0.0 )

      !Calculate upstream terms
      IF (flag_u_at_v) THEN
        upstr_u = u_at_v(i,k) * ( adjust % rho(i+1,k) - adjust % rho(i,k) ) * recipdx
      ELSE
        upstr_u = u_at_v(i,k) * ( adjust % rho(i,k) - adjust % rho(i-1,k) ) * recipdx
      END IF

      IF (flag_w_at_v) THEN
        upstr_w = w_at_v(i,k) * ( adjust % rho(i,k+1) - adjust % rho(i,k) ) *        &
                                Dims % recip_half_kp1_k(k)
      ELSE
        upstr_w = w_at_v(i,k) * ( adjust % rho(i,k) - adjust % rho(i,k-1) ) *        &
                                Dims % recip_half_k_km1(k)
      END IF

      rhodivu  =  adjust % rho(i,k) *                                                &
                  ( adjust % u(i,k) - adjust % u(i-1,k) ) * recipdx +                &
                  ( adjust % w(i,k) - adjust % w(i,k-1) ) * Dims % recip_full_k_km1(k)

      adjust_p1 % r(i,k) = adjust % r(i,k) - deltat_B *    &
                           ( rhodivu +                     &   ! rho div u
                             upstr_w + upstr_u )               ! u dot grad rho
    END DO
  END DO

  ! Update fields
  adjust % r(1:nlongs,1:nlevs)   = adjust_p1 % r(1:nlongs,1:nlevs)
  adjust % rho(1:nlongs,1:nlevs) = 1.0 + adjust % r(1:nlongs,1:nlevs)
  ! Apply boundary conditions
  CALL Boundaries (adjust, set_r=.TRUE., set_rho=.TRUE.)

  !End of temporal loop for tsmall
END DO

!*************************************
!***** Calculate advecting means *****
!*************************************

aves % u_m(1:nlongs,1:nlevs) = 0.5 * ( aves % u_1(1:nlongs,1:nlevs) + aves % u_2(1:nlongs,1:nlevs) )
aves % w_m(1:nlongs,1:nlevs) = 0.5 * ( aves % w_1(1:nlongs,1:nlevs) + aves % w_2(1:nlongs,1:nlevs) )

! Apply Boundary Conditions to advecting means
!---------------------------------------------
aves % u_m(0,0:nlevs+1)         = aves % u_m(nlongs, 0:nlevs+1)
aves % w_m(0,0:nlevs+1)         = aves % w_m(nlongs, 0:nlevs+1)
aves % u_m(nlongs+1,0:nlevs+1)  = aves % u_m(1, 0:nlevs+1)
aves % w_m(nlongs+1,0:nlevs+1)  = aves % w_m(1, 0:nlevs+1)
aves % u_m( 0:nlongs+1,nlevs+1) = aves % u_m( 0:nlongs+1, nlevs)
aves % u_m( 0:nlongs+1,0)       = -1.0 * aves % u_m( 0:nlongs+1,1)
aves % w_m(0:nlongs+1,0)        = 0.0
aves % w_m(0:nlongs+1,nlevs)    = 0.0
aves % w_m(0:nlongs+1,nlevs+1)  = 0.0

!Reassign r
stateout % r(1:nlongs,1:nlevs) = adjust % r(1:nlongs,1:nlevs)
CALL Boundaries (stateout, set_r=.TRUE.)


!***************************
!***** ADVECTION STAGE *****
!***************************

!Calculate advecting interpolations
DO k = 1, nlevs
  DO i = 1, nlongs
    u_at_b(i,k) = half * ( INT_HF(aves % u_m(i-1,k), aves % u_m(i-1,k+1), k, Dims)  + &
                           INT_HF(aves % u_m(i,k), aves % u_m(i,k+1), k, Dims) )
    u_at_v(i,k) = half * ( aves % u_m(i-1, k) + aves % u_m(i,k) )
    w_at_u(i,k) = half * ( INT_FH( aves % w_m(i,k-1), aves % w_m(i, k), k , Dims)   + &
                           INT_FH( aves % w_m(i+1,k-1), aves % w_m(i+1, k), k , Dims))
    w_at_v(i,k) = INT_FH( aves % w_m(i,k-1), aves % w_m(i, k), k , Dims)
  END DO
END DO

DO k = 1, nlevs
  DO i = 1,nlongs

    !Calculate flags
    !---------------
    flag_u      = ( aves % u_m(i,k) .LE. 0.0 )
    flag_w      = ( aves % w_m(i,k) .LE. 0.0 )
    flag_u_at_b = ( u_at_b(i,k) .LE. 0.0 )
    flag_w_at_u = ( w_at_u(i,k) .LE. 0.0 )
    flag_u_at_v = ( u_at_v(i,k) .LE. 0.0 )
    flag_w_at_v = ( w_at_v(i,k) .LE. 0.0 )

    ! u- component
    !-------------
    IF (flag_u) THEN
      upstr_u = aves % u_m(i,k) * ( input % u(i+1,k) - input % u(i,k) ) * recipdx
    ELSE
      upstr_u = aves % u_m(i,k) * ( input % u(i,k) - input % u(i-1,k) ) * recipdx
    END IF
    IF (flag_w_at_u) THEN
      upstr_w = w_at_u(i,k) * ( input % u(i,k+1) - input % u(i,k) ) *                 &
                              Dims % recip_half_kp1_k(k)
    ELSE
      upstr_w = w_at_u(i,k) * ( input % u(i,k) - input % u(i,k-1) ) *                 &
                              Dims % recip_half_k_km1(k)
    END IF
    stateout % u(i,k) = adjust % u(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad u

    ! v- component
    !-------------
    IF (flag_u_at_v) THEN
      upstr_u = u_at_v(i,k) * ( input % v(i+1,k) - input % v(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_v(i,k) * ( input % v(i,k) - input % v(i-1,k) ) * recipdx
    END IF
    IF (flag_w_at_v) THEN
      upstr_w = w_at_v(i,k) * ( input % v(i,k+1) - input % v(i,k) ) *                &
                              Dims % recip_half_kp1_k(k)
    ELSE
      upstr_w = w_at_v(i,k) * ( input % v(i,k) - input % v(i,k-1) ) *                &
                              Dims % recip_half_k_km1(k)
    END IF

    stateout % v(i,k) = adjust % v(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad v

    ! w- component
    !-------------
    IF (flag_u_at_b) THEN
      upstr_u = u_at_b(i,k) * ( input % w(i+1,k) - input % w(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_b(i,k) * ( input % w(i,k) - input % w(i-1,k) ) * recipdx
    END IF
    IF (flag_w) THEN
      upstr_w = aves % w_m(i,k) *                                                    &
                ( input % w(i,k+1) - input % w(i,k) ) *                              &
                Dims % recip_full_kp1_k(k)
    ELSE
      upstr_w = aves % w_m(i,k) *                                                    &
                ( input % w(i,k) - input % w(i,k-1) ) *                              &
                Dims % recip_full_k_km1(k)
    END IF

    stateout % w(i,k) = adjust % w(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad w


    ! b' component
    !-------------
    IF (flag_u_at_b) THEN
      upstr_u = u_at_b(i,k) * ( input % b(i+1,k) - input % b(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_b(i,k) * ( input % b(i,k) -  input % b(i-1,k) ) * recipdx
    END IF
    IF (flag_w) THEN
      upstr_w = aves % w_m(i,k) *                                                    &
                ( input % b(i,k+1) - input % b(i,k) ) *                              &
                Dims % recip_full_kp1_k(k)
    ELSE
      upstr_w = aves % w_m(i,k) *                                                    &
                ( input % b(i,k) - input % b(i,k-1) ) *                              &
                Dims % recip_full_k_km1(k)

    END IF

    stateout % b(i,k) = adjust % b(i,k) - B_dt * ( upstr_u + upstr_w )   ! u . grad b
  END DO
END DO

! Apply boundary conditions, u, v, w, b'
CALL Boundaries (stateout, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_b=.TRUE.)


!*************************
!***** Advect Tracer *****
!*************************

IF (Adv_tracer) THEN
  DO k = 1, nlevs
    DO i = 1,nlongs

      !Calculate flags
      flag_u_at_b = ( u_at_b(i,k) .LE. 0.0 )
      flag_w_at_b = ( aves % w_m(i,k) .LE. 0.0 )

      IF (flag_u_at_b) THEN
        upstr_u = u_at_b(i,k) * ( input % tracer(i+1,k) - input % tracer(i,k) ) * recipdx
      ELSE
        upstr_u = u_at_b(i,k) * ( input % tracer(i,k) - input % tracer(i-1,k) ) * recipdx
      END IF
      IF (flag_w_at_b) THEN
        upstr_w = aves % w_m(i,k) * ( input % tracer(i,k+1) - input % tracer(i,k) ) *    &
                  Dims % recip_full_kp1_k(k)
      ELSE
        upstr_w = aves % w_m(i,k) * ( input % tracer(i,k) - input % tracer(i,k-1) ) *    &
                  Dims % recip_full_k_km1(k)
      END IF

      stateout % tracer (i,k) = input % tracer(i,k) - dt * ( upstr_u + upstr_w )   ! vec(u) . grad tracer
    END DO
  END DO
  CALL Boundaries (stateout, set_tracer=.TRUE.)
END IF

! Update rho
stateout % rho(1:nlongs,1:nlevs) = 1.0 + stateout % r(1:nlongs,1:nlevs)
CALL Boundaries (stateout, set_rho=.TRUE.)

END SUBROUTINE Forward_model
!===================================================================================================
