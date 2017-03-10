!  File includes
!  MODULE DefConsTypes
!  SUBROUTINE SetParams

!===================================================================================================
MODULE DefConsTypes

!*****************************************************
!*   Definition of all global constants              *
!*   and derived variable types                      *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*                                                   *
!*****************************************************

IMPLICIT NONE

!Definition of some new data types
!Define an integer that describes a unified double precision
!-----------------------------------------------------------
INTEGER, PARAMETER        :: ZREAL8 = SELECTED_REAL_KIND(15,307)
INTEGER, PARAMETER        :: wp     = KIND(1.0D0)

! Model parameters
!-----------------
INTEGER, PARAMETER        :: nlongs = 360            ! total number of longitude points
INTEGER, PARAMETER        :: nlevs  = 60             ! total number of vertical levels
INTEGER                   :: ntimesteps              ! number of timesteps to integrate
REAL(ZREAL8)              :: dt     = 4.             ! full timestep
REAL(ZREAL8)              :: dx     = 1500.          ! dx grid resolution
REAL(ZREAL8)              :: H      = 14862.01       ! UM model height
REAL(ZREAL8)              :: deltat
REAL(ZREAL8)              :: Lx
REAL(ZREAL8)              :: dz

! Mathematical and physical constants
!------------------------------------
REAL(ZREAL8), PARAMETER   :: Rd    = 287.058         ! Gas constant for dry air
REAL(ZREAL8), PARAMETER   :: Re    = 6.371E6         ! Mean radius of the earth
REAL(ZREAL8), PARAMETER   :: Cp    = 1005.7          ! Specific heat at constant pressure
REAL(ZREAL8), PARAMETER   :: Cv    = 719.0           ! Specific heat at constant volume
REAL(ZREAL8), PARAMETER   :: g     = 9.81            ! Acceleration due to gravity
REAL(ZREAL8), PARAMETER   :: p00   = 1.0E5           ! Reference surface pressure 1000 hPa
REAL(ZREAL8), PARAMETER   :: kappa = 0.286           ! Rd/Cp
REAL(ZREAL8), PARAMETER   :: pi    = 3.141592654     ! pi

! Tuneable model parameters
!--------------------------
REAL(ZREAL8)              :: f     = 1.0E-4          ! Coriolis Parameter
REAL(ZREAL8)              :: A     = 0.02            ! A is the buoyancy frequency
REAL(ZREAL8)              :: B     = 5.0E-3          ! B accoustic wave speed modulator
REAL(ZREAL8)              :: C     = 1.0E5           ! Constant relating pressure and density
REAL(ZREAL8)              :: BoundSpread = 50.       ! No of grid points to spread boundary discontinuity info
REAL(ZREAL8)              :: theta_r     = 273.

! Useful constants
!-----------------
REAL(ZREAL8)              :: third, half
REAL(ZREAL8)              :: recippi
REAL(ZREAL8)              :: recipdx
REAL(ZREAL8)              :: recipdx2
REAL(ZREAL8)              :: recip2dx
REAL(ZREAL8)              :: alpha_f
REAL(ZREAL8)              :: alpha_N
REAL(ZREAL8)              :: beta_f
REAL(ZREAL8)              :: beta_N
REAL(ZREAL8)              :: recip_alpha_f
REAL(ZREAL8)              :: recip_alpha_N
REAL(ZREAL8)              :: bdiva_f
REAL(ZREAL8)              :: bdiva_N

! Variables to do with reading UM data
! ------------------------------------
! Set this to convert a long/level slice of a SUK-UM file to the domain of this simplified model
LOGICAL                   :: make_ics_from_um = .TRUE.    ! To read slice of UM data to make init conds
CHARACTER                 :: datadirUM*100                ! Directory to do with UM data
CHARACTER                 :: init_um_file*100             ! Input UM filename
CHARACTER                 :: model_ics_data_out_file*100  ! Simplified model ic file
INTEGER                   :: latitude = 144               ! The latitude to be extracted
LOGICAL                   :: Regular_vert_grid = .TRUE.   ! Set to use regular vertical levels
LOGICAL                   :: gravity_wave_switch = .FALSE.! To set u=0 (simulate gws)

! Variables to do with running the forward model
! ----------------------------------------------
LOGICAL                   :: run_toy_model = .TRUE.       ! To run simplified model
CHARACTER                 :: datadirMODEL*100             ! Directory to do with simplified model data
CHARACTER                 :: model_ics_data_read_file*36  ! Input file
CHARACTER                 :: model_output_file*36         ! Dump file
CHARACTER                 :: diagnostics_file*36          ! For timestep diagnostics (e.g. energy)
REAL(ZREAL8)              :: runlength = 60.0             ! Runlength in seconds
INTEGER                   :: ndumps = 10                  ! number of dump times
LOGICAL                   :: convection_switch = .FALSE.  ! Set to ?
LOGICAL                   :: pressure_pert = .FALSE.      ! To introduce a blob of pressure in the ics
INTEGER                   :: press_source_i = 180         ! Long pos of centre of pressure blob
INTEGER                   :: press_source_k = 30          ! Level pos of centre of pressure blob
INTEGER                   :: x_scale = 80                 ! Long size of pressure blob
INTEGER                   :: z_scale = 3                  ! Level size of pressure blob
REAL(ZREAL8)              :: press_amp = 0.01             ! Amplitude of pressure blob
LOGICAL                   :: Adv_tracer = .FALSE.         ! Set to advect tracer in calculations
INTEGER                   :: Tracer_level = 40            ! Location to put unit tracer at t=0
LOGICAL                   :: Lengthscale_diagnostics = .FALSE.
                                                          ! Set to do lengthscale diagnostics at final time

! Variables to do with linear analysis
! ------------------------------------
LOGICAL                   :: do_linear_analysis = .FALSE. ! To do linear analysis
CHARACTER                 :: linear_analysis_odir*100     ! Separate directory for the linear analysis data
CHARACTER                 :: wavespeed_experiment*36      ! Name of wave speed experiment (part of o/p file names)



!**************************************************************************************************
! Declare Compound Types
!**************************************************************************************************

!---------------------------------
TYPE um_data_type
  REAL(ZREAL8) :: longs_u(1:nlongs)
  REAL(ZREAL8) :: longs_v(1:nlongs)
  REAL(ZREAL8) :: half_levs(1:nlevs+1)
  REAL(ZREAL8) :: full_levs(0:nlevs)
  REAL(ZREAL8) :: u(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: v(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: w(1:nlongs,0:nlevs)
  REAL(ZREAL8) :: density(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: theta(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: exner_pressure(1:nlongs,1:nlevs+1)
  REAL(ZREAL8) :: orog_height(1:nlongs)
END TYPE um_data_type
!---------------------------------

!---------------------------------
TYPE Dimensions_type
  REAL(ZREAL8) :: longs_u(0:nlongs+1)
  REAL(ZREAL8) :: longs_v(0:nlongs+1)
  REAL(ZREAL8) :: half_levs(0:nlevs+1)
  REAL(ZREAL8) :: full_levs(0:nlevs+1)
  ! Variables required for vertical interpolation
  REAL(ZREAL8) :: a1(1:nlevs+1)
  REAL(ZREAL8) :: b1(1:nlevs+1)
  REAL(ZREAL8) :: a2(0:nlevs)
  REAL(ZREAL8) :: b2(0:nlevs)
  REAL(ZREAL8) :: recip_half_kp1_k(0:nlevs)
  REAL(ZREAL8) :: recip_half_k_km1(1:nlevs+1)
  REAL(ZREAL8) :: recip_full_kp1_k(0:nlevs)
  REAL(ZREAL8) :: recip_full_k_km1(1:nlevs+1)
END TYPE Dimensions_type
!---------------------------------

!---------------------------------
TYPE model_vars_type
  ! Horizontal grid is Awakara C grid
  ! Vertical grid is Charney-Phillips
  REAL(ZREAL8) :: u(0:nlongs+1,0:nlevs+1)               ! Zonal wind perturbation
  REAL(ZREAL8) :: v(0:nlongs+1,0:nlevs+1)               ! Meridional wind perturbation
  REAL(ZREAL8) :: w(0:nlongs+1,0:nlevs+1)               ! Vertical wind perturbation
  REAL(ZREAL8) :: r(0:nlongs+1,0:nlevs+1)               ! rho density perturbation
  REAL(ZREAL8) :: b(0:nlongs+1,0:nlevs+1)               ! buoyancy perturbation
  REAL(ZREAL8) :: rho(0:nlongs+1,0:nlevs+1)             ! rho full field
  REAL(ZREAL8) :: b_ef(0:nlongs+1,0:nlevs+1)            ! Effective buoyancy
  REAL(ZREAL8) :: tracer(0:nlongs+1,0:nlevs+1)          ! Tracer
  REAL(ZREAL8) :: hydro_imbal(0:nlongs+1,0:nlevs+1)     ! Hydrostatic imbalance
  REAL(ZREAL8) :: geost_imbal(0:nlongs+1,0:nlevs+1)     ! Geostrophic imbalance
  REAL(ZREAL8) :: vert_mom_source(0:nlongs+1,0:nlevs+1) ! Vertical momentum source
  REAL(ZREAL8) :: Kinetic_Energy
  REAL(ZREAL8) :: Buoyant_Energy
  REAL(ZREAL8) :: Elastic_Energy
  REAL(ZREAL8) :: Total_Energy
END TYPE model_vars_type
!---------------------------------


!---------------------------------
TYPE Averages_type
  REAL(ZREAL8) :: u_1(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: u_2(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: u_m(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_1(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_2(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_m(0:nlongs+1, 0:nlevs+1)
END TYPE Averages_type
!---------------------------------


! Define namelist
! ---------------
NAMELIST / UserOptions /                                                                                 &
! Reading and processing UM data
  make_ics_from_um, datadirUM, init_um_file, model_ics_data_out_file, latitude, Regular_vert_grid,       &
  gravity_wave_switch, f, A, B, C, BoundSpread, Tracer_level,                                            &
! Run the forward model
  run_toy_model, datadirMODEL, model_ics_data_read_file, model_output_file, diagnostics_file, dt, dx, H, &
  runlength, ndumps,                                                                                     &
  convection_switch, pressure_pert, press_source_i, press_source_k, x_scale, z_scale, press_amp,         &
  Adv_tracer, Lengthscale_diagnostics,                                                                   &
! Linear analysis
  do_linear_analysis, linear_analysis_odir, wavespeed_experiment



END MODULE DefConsTypes
!===================================================================================



!===================================================================================
SUBROUTINE SetParams

! Read namelists and set derived parameters

USE DefConsTypes, ONLY :     &
  nlongs, nlevs, pi,         &
  UserOptions,               &
  make_ics_from_UM,          &
  datadirUM,                 &
  init_um_file,              &
  model_ics_data_out_file,   &
  latitude,                  &
  Regular_vert_grid,         &
  gravity_wave_switch,       &
  run_toy_model,             &
  datadirMODEL,              &
  model_ics_data_read_file,  &
  model_output_file,         &
  diagnostics_file,          &
  dt, dx, H,                 &
  runlength,                 &
  ndumps, ntimesteps,        &
  f, A, B, C,                &
  BoundSpread,               &
  Adv_tracer, Tracer_level,  &
  Lengthscale_diagnostics,   &
  convection_switch,         &
  pressure_pert,             &
  press_source_i,            &
  press_source_k,            &
  x_scale, z_scale,          &
  press_amp,                 &
  do_linear_analysis,        &
  wavespeed_experiment,      &
  Lx, deltat, dz,            &
  recipdx2, recipdx,         &
  recip2dx,                  &
  third, half,               &
  recippi,                   &
  alpha_f, alpha_N,          &
  beta_f, beta_N,            &
  recip_alpha_f,             &
  recip_alpha_N,             &
  bdiva_f, bdiva_N



IMPLICIT NONE

INTEGER                   :: unitno, ErrStatus

! Read the namelist containing user's choice of parameters
! --------------------------------------------------------
ErrStatus = 0
unitno    = 12
OPEN (UNIT=unitno, FILE='UserOptions.nl', IOSTAT=ErrStatus)
IF (ErrStatus == 0) THEN
  READ (unitno, NML=UserOptions)
  CLOSE (unitno)
ELSE
  PRINT *, 'Expecting a namelist file called UserOptions.nl'
  STOP
END IF


! Set the derived parameters
! --------------------------
Lx            = dx * REAL(nlongs)
deltat        = dt / 2.
dz            = H / nlevs
ntimesteps    = runlength / dt
recipdx2      = 1. / (dx * dx)
recipdx       = 1. / dx
recip2dx      = 1. / ( 2. * dx )
third         = 1. / 3.
half          = 1. / 2.
recippi       = 1. / pi
alpha_f       = 1. + deltat*deltat * f*f / 4.
alpha_N       = 1. + deltat*deltat * A*A / 4.
beta_f        = 1. - deltat*deltat * f*f / 4.
beta_N        = 1. - deltat*deltat * A*A / 4.
recip_alpha_f = 1. / alpha_f
recip_alpha_N = 1. / alpha_N
bdiva_f       = beta_f / alpha_f
bdiva_N       = beta_N / alpha_N


PRINT *
PRINT *, 'User options'
PRINT *, '============================================================'
PRINT *, 'make_ics_from_UM         = ', make_ics_from_UM
PRINT *, 'datadirUM                = ', datadirUM
PRINT *, 'init_um_file             = ', init_um_file
PRINT *, 'model_ics_data_out_file  = ', model_ics_data_out_file
PRINT *, 'latitude                 = ', latitude
PRINT *, 'Regular_vert_grid        = ', Regular_vert_grid
PRINT *, 'gravity_wave_switch      = ', gravity_wave_switch
PRINT *, 'f                        = ', f
PRINT *, 'A                        = ', A
PRINT *, 'B                        = ', B
PRINT *, 'C                        = ', C
PRINT *, 'BoundSpread              = ', BoundSpread
PRINT *, 'Tracer_level             = ', Tracer_level
PRINT *, 'run_toy_model            = ', run_toy_model
PRINT *, 'datadirMODEL             = ', datadirMODEL
PRINT *, 'model_ics_data_read_file = ', model_ics_data_read_file
PRINT *, 'model_output_file        = ', model_output_file
PRINT *, 'diagnostics_file         = ', diagnostics_file
PRINT *, 'dt                       = ', dt
PRINT *, 'dx                       = ', dx
PRINT *, 'H                        = ', H
PRINT *, 'runlength                = ', runlength
PRINT *, 'ndumps                   = ', ndumps
PRINT *, 'convection_switch        = ', convection_switch
PRINT *, 'pressure_pert            = ', pressure_pert
PRINT *, 'press_source_i           = ', press_source_i
PRINT *, 'press_source_k           = ', press_source_k
PRINT *, 'x_scale                  = ', x_scale
PRINT *, 'z_scale                  = ', z_scale
PRINT *, 'press_amp                = ', press_amp
PRINT *, 'Adv_tracer               = ', Adv_tracer
PRINT *, 'Lengthscale_diagnostics  = ', Lengthscale_diagnostics
PRINT *, 'do_linear_analysis       = ', do_linear_analysis
PRINT *, 'wavespeed_experiment     = ', wavespeed_experiment
PRINT *
PRINT *, 'Derived variables'
PRINT *, '============================================================'
PRINT *, 'UM data file             = ', TRIM(datadirUM) // '/' // init_um_file
PRINT *, 'Output initial conds     = ', TRIM(datadirMODEL) // '/' // model_ics_data_out_file
PRINT *, 'Input initial conds      = ', TRIM(datadirMODEL) // '/' // model_ics_data_read_file
PRINT *, 'Output dumps             = ', TRIM(datadirMODEL) // '/' // model_output_file
PRINT *, 'Diagnostics file         = ', TRIM(datadirMODEL) // '/' // diagnostics_file
PRINT *, 'Lx                       = ', Lx
PRINT *, 'deltat                   = ', deltat
PRINT *, 'dz                       = ', dz
PRINT *, 'ntimesteps               = ', ntimesteps
PRINT *, 'recipdx2                 = ', recipdx2
PRINT *, 'recipdx                  = ', recipdx
PRINT *, 'third                    = ', third
PRINT *, 'recippi                  = ', recippi
PRINT *, 'alpha_f                  = ', alpha_f
PRINT *, 'beta_f                   = ', beta_f
PRINT *, 'beta_N                   = ', beta_N
PRINT *, 'recip_alpha_f            = ', recip_alpha_f
PRINT *, 'recip_alpha_N            = ', recip_alpha_N
PRINT *, 'bdiva_f                  = ', bdiva_f
PRINT *, 'bdiva_N                  = ', bdiva_N
PRINT *, '============================================================'

END SUBROUTINE SetParams
