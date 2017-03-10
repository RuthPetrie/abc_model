PROGRAM Main

!*****************************************************
!*   Main Driver program to run fowrard model        *
!*   and compute diagnostics, including              *
!*   numerical linear analysis                       *
!*                                                   *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*                                                   *
!*****************************************************

! Compile Statement
!===================

! Out of bounds checking: debug mode (carrot)
!--------------------------------------------
! f90 -C -dalign -I/opt/graphics/include Main.f90 DefConsType.f90 Linear_Analysis.f90 Write_data.f90 Model.f90 Functions.f90 Boundaries.f90 Read_data.f90 UM_data_proc.f90 -o Main.out -M/opt/graphics/include -L/opt/graphics/lib -lnetcdf -lhdf5_hl -lhdf5 -M/opt/local/lib/nag_mod_dir -L/opt/tools/lib -lnagfl90

! Standard (carrot)
!------------------
! f90 -C -dalign -I/opt/graphics/include Main.f90 DefConsType.f90  Linear_Analysis.f90 Write_data.f90 Model.f90 Functions.f90 Boundaries.f90  Read_data.f90 UM_data_proc.f90 Diagnostics.f90 Forcings.f90 -o Main.out -M/opt/graphics/include -L/opt/graphics/lib -lnetcdf -lhdf5_hl -lhdf5 -M/opt/local/lib/nag_mod_dir -L/opt/tools/lib -lnagfl90

! Optimized (carrot)
!-------------------
! f90  -fast -I/opt/graphics/include Main.f90 DefConsType.f90  Linear_Analysis.f90 Write_data.f90 Model.f90 Functions.f90 Boundaries.f90 Read_data.f90 UM_data_proc.f90 -o Main.out -M/opt/graphics/include -L/opt/graphics/lib -lnetcdf -lhdf5_hl -lhdf5 -M/opt/local/lib/nag_mod_dir -L/opt/tools/lib -lnagfl90


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    um_data_type,                &
    dimensions_type,             &
    model_vars_type,             &
    ZREAL8,                      &
    dt,                          &
    nlongs, nlevs,               &
    make_ics_from_UM,            &
    init_um_file,                &
    datadirUM,                   &
    model_ics_data_out_file,     &
    latitude,                    &
    Regular_vert_grid,           &
    gravity_wave_switch,         &
    run_toy_model,               &
    datadirMODEL,                &
    linear_analysis_odir,        &
    model_ics_data_read_file,    &
    model_output_file,           &
    diagnostics_file,            &
    runlength,                   &
    ndumps,                      &
    ntimesteps,                  &
    convection_switch,           &
    pressure_pert,               &
    press_source_i,              &
    press_source_k,              &
    x_scale, z_scale,            &
    press_amp,                   &
    do_linear_analysis,          &
    wavespeed_experiment


IMPLICIT NONE


! Declare variables
!==========================
TYPE(um_data_type)            :: init_um_data
TYPE(dimensions_type)         :: dims
TYPE(model_vars_type)         :: init_state, output_state
INTEGER, ALLOCATABLE          :: times(:)
REAL(ZREAL8)                  :: dist_x2, dist_z2, x_scale2, z_scale2
INTEGER                       :: dist_x, dist_z

LOGICAL                       :: read_processed_data
LOGICAL                       :: rescale
INTEGER                       :: balance_diags
INTEGER                       :: x_scale_length, z_scale_length

INTEGER                       :: j, x, z


PRINT*, '*************************************************************************'
PRINT*, 'RUNNING MAIN'
PRINT*, 'Simplified convective-scale model code'
PRINT*, 'R.E. Petrie and R.N.Bannister'
PRINT*, 'Last revised June 2014'
PRINT*, '*************************************************************************'


!*********************************************************************************
! INITIALISE
!*********************************************************************************

! Read namelist and set derived parameters
CALL SetParams
! Set state to zero
CALL Initialise_model_vars (output_state)


!*********************************************************************************
! Generate initial data (in variable init_state)
!*********************************************************************************
IF (.NOT.do_linear_analysis) THEN
  IF (make_ics_from_UM) THEN

    ! Read in raw UM data, store in init_um_data
    !-------------------------------------------
    CALL Initialise_um_data (init_um_data)
    PRINT*, 'Reading UM data ...'
    CALL Read_um_data_2d (init_um_data, TRIM(datadirUM) // '/' // TRIM(init_um_file), latitude)
    PRINT*, '--done'

    ! Process um data
    !----------------
    ! Set grid
    CALL Initialise_dims (dims)
    PRINT*, 'Setting Grid ...'

    CALL Set_grid (init_um_data, dims, Regular_vert_grid)
    PRINT*, '--done'
   
    ! Set some commonly-used constants
    PRINT*, 'Setting some constants ...'
    CALL Set_ht_dep_cons (dims)
    PRINT*, '-- done'

    ! Define variables for simplified model from UM data, store in init_state
    PRINT*, 'Processing UM data to make it compatible with simplified model ...'
    CALL Initialise_model_vars (init_state)
 
    CALL Process_um_data (init_um_data, init_state, dims, TRIM(datadirMODEL) // '/' // TRIM(model_ics_data_out_file))
  ELSE

    ! Read in preprocessed UM data store in init_state
    PRINT*, 'Reading in processed data ...'
    CALL Read_state_2d (TRIM(datadirMODEL) // '/' // TRIM(model_ics_data_read_file), init_state, dims, 1)

    ! Set some commonly-used constants
    PRINT*, 'Setting some constants ...'
    CALL Set_ht_dep_cons (dims)
    PRINT*, '-- done'

  END IF


  !*********************************************************************************
  ! Run forward model
  !*********************************************************************************

  IF (run_toy_model) THEN
    PRINT*, ' Run forward model'
    PRINT*, '------------------'

    IF (convection_switch) THEN
      ! What does this do?
      init_state % b (:,:)          = 0.0
      init_state % w (:,:)          = 0.0
      init_state % v (:,:)          = 0.0
      init_state % u (:,:)          = 0.0
      init_state % r (:,:)          = 0.0
      init_state % rho (:,:)        = 0.0
      init_state % geost_imbal(:,:) = 0.0
      init_state % hydro_imbal(:,:) = 0.0
    END IF


    IF (pressure_pert) THEN
      ! Initalise all fields to zero
      init_state % b (:,:)          = 0.0
      init_state % w (:,:)          = 0.0
      init_state % v (:,:)          = 0.0
      init_state % u (:,:)          = 0.0
      init_state % r (:,:)          = 0.0
      init_state % rho (:,:)        = 0.0
      init_state % geost_imbal(:,:) = 0.0
      init_state % hydro_imbal(:,:) = 0.0

      x_scale2 = REAL(x_scale * x_scale)
      z_scale2 = REAL(z_scale * z_scale)
      DO x = 1, nlongs
	DO z = 1, nlevs
          dist_x  = press_source_i - x
          dist_z  = press_source_k - z
          dist_x2 = REAL(dist_x * dist_x)
          dist_z2 = REAL(dist_z * dist_z)
          init_state % r(x,z) = press_amp * EXP ( -1. * (dist_x2/x_scale2 + dist_z2/z_scale2) )
	END DO
      END DO
    END IF

    CALL ModelDriver ( init_state, output_state, Dims, ntimesteps, ndumps,      &
                       TRIM(datadirMODEL) // '/' // TRIM(model_output_file),    &
                       TRIM(datadirMODEL) // '/' // TRIM(diagnostics_file),     &
                       convection_switch )
    PRINT*, '-- done'

  END IF

ELSE

  PRINT*, ' Doing linear analysis'
  PRINT*, '----------------------'
  CALL Linear_Analysis (TRIM(linear_analysis_odir), TRIM(wavespeed_experiment))
  PRINT*, '-- done'
END IF

END PROGRAM Main
