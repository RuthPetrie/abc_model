! Define namelist
! ---------------
&UserOptions

! Reading and processing UM data
! ------------------------------
  make_ics_from_um         = .TRUE.
  datadirUM                = '/net/quince/export/quince/data-03/wx019276/ABC_MODEL/INITCONDS'
  init_um_file             = 'UM_InitConds.nc'
  model_ics_data_out_file  = 'Processed_Initial_Data.nc'
  latitude                 = 144
  Regular_vert_grid        = .TRUE.
  gravity_wave_switch      = .FALSE.
  f                        = 1.0E-4
  A                        = 0.02
  B                        = 0.005
  C                        = 1.0E5
  BoundSpread              = 100.
  Tracer_level             = 20
! Run the forward model
! ------------------------------
  run_toy_model            = .TRUE.
  datadirMODEL             = '/net/quince/export/quince/data-03/wx019276/ABC_MODEL/OUTPUT'
  model_ics_data_read_file = 'Processed_Initial_Data.nc'
  model_output_file        = 'Test_B_is_0p005.nc'
  diagnostics_file         = 'Test_B_is_0p005.dat'
  dt                       = 0.1
  runlength                = 21600.0
  ndumps                   = 36
  convection_switch        = .FALSE.
  pressure_pert            = .FALSE.
  press_source_i           = 180
  press_source_k           = 30
  x_scale                  = 80
  z_scale                  = 3
  press_amp                = 0.01
  Adv_tracer               = .TRUE.
  Lengthscale_diagnostics  = .TRUE.
! Linear analysis
! ------------------------------
  do_linear_analysis       = .FALSE.
  linear_analysis_odir     = '/net/quince/export/quince/data-03/wx019276/ABC_MODEL/OUTPUT'
  wavespeed_experiment     = 'TestExp'
/
