!  File includes

!  SUBROUTINE Initialise_um_data
!  SUBROUTINE Initialise_dims
!  SUBROUTINE Initialise_model_vars
!  SUBROUTINE Initialise_Averages

!===================================================================================
SUBROUTINE Initialise_um_data (state)

! Initialise um-data variables

USE DefConsTypes, ONLY :   &
    um_data_type,          &
    nlongs, nlevs

IMPLICIT NONE

TYPE(um_data_type), INTENT(INOUT)   :: state

state % longs_u(1:nlongs)                  = 0.0
state % longs_v(1:nlongs)                  = 0.0
state % half_levs(1:nlevs+1)               = 0.0
state % full_levs(0:nlevs)                 = 0.0
state % u(1:nlongs,1:nlevs)                = 0.0
state % v(1:nlongs,1:nlevs)                = 0.0
state % w(1:nlongs,0:nlevs)                = 0.0
state % density(1:nlongs,1:nlevs)          = 0.0
state % theta(1:nlongs,1:nlevs)            = 0.0
state % exner_pressure(1:nlongs,1:nlevs+1) = 0.0
state % orog_height(1:nlongs)              = 0.0

END SUBROUTINE Initialise_um_data

!===================================================================================


!===================================================================================
SUBROUTINE Initialise_dims (state)

! Initialise dimension variables

USE DefConsTypes, ONLY :   &
    Dimensions_type,       &
    nlongs, nlevs

IMPLICIT NONE

TYPE(Dimensions_type), INTENT(INOUT)   :: state

state % longs_u(0:nlongs+1)  = 0.0
state % longs_v(0:nlongs+1)  = 0.0
state % half_levs(0:nlevs+1) = 0.0
state % full_levs(0:nlevs+1) = 0.0
state % a1(1:nlevs)          = 0.0
state % b1(1:nlevs)          = 0.0
state % a2(1:nlevs)          = 0.0
state % b2(1:nlevs)          = 0.0

END SUBROUTINE Initialise_dims

!===================================================================================


!===================================================================================
SUBROUTINE Initialise_model_vars(state)

! Initialise model variables

USE DefConsTypes, ONLY :     &
    model_vars_type,         &
    nlongs, nlevs

IMPLICIT NONE

TYPE(model_vars_type), INTENT(INOUT)   :: state

state % u(0:nlongs+1,0:nlevs+1)           = 0.0
state % v(0:nlongs+1,0:nlevs+1)           = 0.0
state % w(0:nlongs+1,0:nlevs+1)           = 0.0
state % r(0:nlongs+1,0:nlevs+1)           = 0.0     ! density perturbation
state % b(0:nlongs+1,0:nlevs+1)           = 0.0     ! buoyancy perturbation
state % rho(0:nlongs+1,0:nlevs+1)         = 0.0     ! density full field
state % b_ef(0:nlongs+1,0:nlevs+1)        = 0.0     ! Effective buoyancy
state % tracer(0:nlongs+1,0:nlevs+1)      = 0.0
state % hydro_imbal(0:nlongs+1,0:nlevs+1) = 0.0
state % geost_imbal(0:nlongs+1,0:nlevs+1) = 0.0
state % Kinetic_Energy                    = 0.0
state % Buoyant_Energy                    = 0.0
state % Elastic_Energy                    = 0.0
state % Total_Energy                      = 0.0

END SUBROUTINE Initialise_model_vars

!===================================================================================


!===================================================================================
SUBROUTINE Initialise_Averages(state)

! Initialise averages variables

USE DefConsTypes, ONLY :     &
    Averages_type,           &
    nlongs, nlevs

IMPLICIT NONE

TYPE(Averages_type), INTENT(INOUT)   :: state

state % u_1(0:nlongs+1, 0:nlevs+1) = 0.0
state % u_2(0:nlongs+1, 0:nlevs+1) = 0.0
state % u_m(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_1(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_2(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_m(0:nlongs+1, 0:nlevs+1) = 0.0

END SUBROUTINE Initialise_Averages

!===================================================================================

