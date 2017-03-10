!===================================================================================================
! This file contains all subroutines to process inital data taken from
! the um and convert it to be consistent with the requirements of
! the toy model.
!
! In this file are subroutines
! SUBROUTINE Process_um_data (um_data, state, dims, fileout, gravity_wave_switch, dump)
! SUBROUTINE Regularize_grid(init_um_data, dims,regularize)
! SUBROUTINE Read_um_data_2d ( umdata, filename, latitude )
!===================================================================================================
SUBROUTINE Process_um_data (um_data, state, dims, fileout)

!*******************************
!* Subroutine to process um    *
!* data for toy model          *
!*                             *
!*                             *
!*                             *
!* Version history             *
!* R. Petrie, 1.0: 6-6-2011    *
!* R. Petrie, 3.0: 30-7-2013   *
!* R. Bannister, Oct 2013      *
!*                             *
!*******************************

USE DefConsTypes, ONLY :   &
    um_data_type,          &
    dimensions_type,       &
    model_vars_type,       &
    ZREAL8,                &
    nlongs, nlevs,         &
    C, f,                  &
    Re, dx, recipdx,       &
    recip2dx, g, theta_r,  &
    Rd, kappa, p00,        &
    Adv_tracer,            &
    Tracer_level,          &
    gravity_wave_switch

IMPLICIT NONE

INCLUDE "Boundaries.interface"


! Declare Top level variables
TYPE(um_data_type),    INTENT(IN)     :: um_data
TYPE(model_vars_type), INTENT(INOUT)  :: state
TYPE(dimensions_type), INTENT(IN)     :: dims
CHARACTER(LEN=*),      INTENT(IN)     :: fileout

! Declare local variables
REAL(ZREAL8)                          :: cons, Re2
REAL(ZREAL8)                          :: rhopsumtemp, rhopsum(1:nlevs)
REAL(ZREAL8)                          :: rhobound(1:nlevs)
REAL(ZREAL8)                          :: diff(1:nlevs+1), ux, uxm1, px, px1, pxm1, pxp1, deltaz, dpdz
REAL(ZREAL8)                          :: theta_sum, theta_z(1:nlevs)
REAL(ZREAL8)                          :: pprimed(1:nlongs,1:nlevs), pav
REAL(ZREAL8)                          :: rhoprimed(1:nlongs,1:nlevs), rhoav(1:nlevs), rho0
REAL(ZREAL8)                          :: C_um
REAL(ZREAL8)                          :: ratio, ratioav, stddev
REAL(ZREAL8)                          :: INT_HF
INTEGER                               :: x, z, i, j, k
INTEGER, ALLOCATABLE                  :: times(:)


! Estimate C from the UM (method 1 - the formula given in the paper Sec 2.1.2, point 4)
! -------------------------------------------------------------------------------------
C_um = Rd * theta_r * SUM(um_data % exner_pressure(1:nlongs,1:nlevs)) / &
        ( (1.-kappa) * REAL(nlongs*nlevs) )
PRINT *, 'Value of C found from UM exner data is ', C_um

! Estimate C from the UM (method 2 - using p'/rho')
! -------------------------------------------------

! Calculate pprimed
pprimed(1:nlongs,1:nlevs) = p00 * (um_data % exner_pressure(1:nlongs,1:nlevs))**(1./kappa)
DO z = 1, nlevs
  pav                     = SUM(pprimed(1:nlongs,z)) / REAL(nlongs)
  pprimed(1:nlongs,z)     = pprimed(1:nlongs,z) - pav
END DO

! Calculate rhoprimed
rhoprimed(1:nlongs,1:nlevs) = um_data % density(1:nlongs,1:nlevs) / (Re * Re)
DO z = 1, nlevs
  rhoav(z)                  = SUM(rhoprimed(1:nlongs,z)) / REAL(nlongs)
  rhoprimed(1:nlongs,z)     = rhoprimed(1:nlongs,z) - rhoav(z)
END DO

! Calculate average ratio p'/rho'
ratioav = 0.
DO x = 1, nlongs
  DO z = 1, nlevs
    ratio   = pprimed(x,z) / rhoprimed(x,z)
    ratioav = ratioav + ABS(ratio)
  END DO
END DO
ratioav = ratioav / REAL(nlongs * nlevs)

! Calculate standard deviation of this
stddev = 0
DO x = 1, nlongs
  DO z = 1, nlevs
    ratio  = pprimed(x,z) / rhoprimed(x,z) - ratioav
    stddev = stddev + ratio * ratio
  END DO
END DO
stddev = SQRT(stddev / REAL(nlongs * nlevs))
PRINT *, 'Value of C found from p-primed/rho-primed is ', ratioav
PRINT *, '  with standard deviation ', stddev


! Caculate global average rho
rho0 = SUM(rhoav(1:nlevs)) / REAL(nlevs)
PRINT *, 'Global average rho (for rho0) = ', rho0

! Set u from UM
!--------------
PRINT*, 'Setting u ...'
state % u(1:nlongs, 1:nlevs) = um_data % u(1:nlongs, 1:nlevs)
CALL BoundaryMod (state % u(1:nlongs, 1:nlevs))
DO z =1, nlevs
  cons                  = SUM(state % u(1:nlongs, z)) / REAL(nlongs)
  state % u(1:nlongs,z) = state % u(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_u=.TRUE.)

! Set v from UM and enforce integral of v dx to be 0 on each vert level
! ----------------------------------------------------------------------
PRINT*, 'Setting v ...'
state % v(1:nlongs, 1:nlevs) = um_data % v(1:nlongs, 1:nlevs)
CALL BoundaryMod (state % v(1:nlongs, 1:nlevs))
DO z =1, nlevs
  cons                  = SUM(state % v(1:nlongs, z)) / REAL(nlongs)
  state % v(1:nlongs,z) = state % v(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_v=.TRUE.)


! Set p by integrating the geostrophic balance eq and enforce integral of p dx to be 0 on each vert level
! -------------------------------------------------------------------------------------------------------
PRINT*, 'Setting p by integrating geostrophic balance ...'
DO z = 1, nlevs
  state % r(1,z) = 0.
  DO x = 2, nlongs
    state % r(x,z) = state % r(x-1,z) + f * dx * (state % v(x,z) + state % v(x-1,z)) / (2. * C)
  END DO
  cons                  = SUM(state % r(1:nlongs,z)) / REAL(nlongs)
  state % r(1:nlongs,z) = state % r(1:nlongs,z) - cons
END DO
CALL Boundaries (state, set_r=.TRUE.)

! Set rho
! -------
PRINT*, 'Setting rho ...'
state % rho(1:nlongs,1:nlevs) = state % r(1:nlongs,1:nlevs) + 1.0
CALL Boundaries (state, set_rho=.TRUE.)

! Set bp from the hydrostatic relationship
! ----------------------------------------
PRINT*, 'Setting bp from hydrostatic balance ...'
DO z = 1, nlevs
  DO x = 1, nlongs
    state % b(x,z) = C * (state % r(x,z+1) - state % r(x,z)) /      &
                         (Dims % half_levs(z+1) - Dims % half_levs(z))
  END DO
END DO
DO z =1, nlevs
  cons                  = SUM(state % b(1:nlongs, z)) / REAL(nlongs)
  state % b(1:nlongs,z) = state % b(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_b=.TRUE.)

! Calculate w from the continuity equation
! ----------------------------------------
PRINT*, 'Setting w from the continuity equation ...'
DO x = 1, nlongs
  ! Calculate the d(rho u)/dx for a column at this longitude
  DO z = 1, nlevs
    ! u interpolated to full-level
    ux      = INT_HF(state % u(x,z), state % u(x,z+1), z, Dims)
    uxm1    = INT_HF(state % u(x-1,z), state % u(x-1,z+1), z, Dims)
    ! rho interpolated to full-level
    px    = INT_HF(state % rho(x,z), state % rho(x,z+1), z, Dims)
    pxm1  = INT_HF(state % rho(x-1,z), state % rho(x-1,z+1), z, Dims)
    pxp1  = INT_HF(state % rho(x+1,z), state % rho(x+1,z+1), z, Dims)
    ! d(p u)/dx
    diff(z) = recip2dx * ( (pxp1 + px) * ux - (px + pxm1) * uxm1 )
  END DO

  ! Linear extrapolation of diff to nlevs + 1
  diff(nlevs+1) = diff(nlevs) +                                              &
                  ( diff(nlevs) - diff(nlevs-1) ) *                          &
                  ( Dims % full_levs(nlevs+1) - Dims % full_levs(nlevs) ) /  &
                  ( Dims % full_levs(nlevs) - Dims % full_levs(nlevs-1) )

  ! Integrate the continuity equation from the ground upwards
  state % w(x,0) = 0.

  DO z = 1, nlevs - 1
    ! rho interpolated to full-level
    px = INT_HF(state % rho(x,z), state % rho(x,z+1), z, Dims)
    ! Level width
    deltaz  = Dims % half_levs(z+1) - Dims % half_levs(z)
    ! dp/dz
    dpdz  = (state % rho(x,z+1) - state % rho(x,z)) / deltaz
    ! Increment w
    state % w(x,z) = 0.0 * state % w(x,z-1) - (diff(z) + dpdz * state % w(x,z-1)) * deltaz / px
  END DO
  state % w(x,nlevs) = 0.0

END DO
DO z =1, nlevs
  cons                  = SUM(state % w(1:nlongs, z)) / REAL(nlongs)
  state % w(1:nlongs,z) = state % w(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_w=.TRUE.)

! Set the tracer
! --------------
IF (Adv_tracer) THEN
  PRINT*, 'Setting the tracer ...'
  DO i = 1, 4
    DO k = 1, 5
      state % tracer(INT(REAL(nlongs*i)/5.0)-2:INT(REAL(nlongs*i)/5.0)+2, k*10) = 1.0
      WRITE (*,*) 'Tracer at position',  i, INT(REAL(nlongs*i)/5.0)
    END DO
  END DO
  CALL Boundaries (state, set_tracer=.TRUE.)
END IF

! Add some gravity wave noise by setting u=0
! ------------------------------------------
IF (gravity_wave_switch) state % u(0:nlongs+1, 0:nlevs+1) = 0.



!***************************************************************************************************


! Calculate diagnostics
!------------------------
CALL Calc_hydro(state, Dims)
CALL Calc_geost(state, Dims)
CALL Energy(state)


!***************************************************************************************************
CALL Write_state_2d (fileout,state, Dims,1,0,0)
PRINT*,'-------------------------------------------------------------------------------'
PRINT*,'    Netcdf inital conds file: '
PRINT*,'    xconv -i ',fileout,'&'
PRINT*,'-------------------------------------------------------------------------------'!


!***************************************************************************************************


END SUBROUTINE Process_um_data


!===================================================================================================
SUBROUTINE Set_grid (init_um_data, dims, regular_grid)

!********************************
!* Subroutine to set model grid *
!*                              *
!*                              *
!* Set regular_grid = 1 to use a  *
!* regular vertical grid        *
!*                              *
!* R. Petrie                    *
!* version 1.0                  *
!* 6/6/2011                     *
!********************************

USE DefConsTypes, ONLY :   &
    ZREAL8,                &
    um_data_type,          &
    dimensions_type,       &
    nlongs, nlevs,         &
    dx

IMPLICIT NONE

! Declare Top level variables
TYPE(um_data_type),    INTENT(IN)     :: init_um_data
TYPE(dimensions_type), INTENT(INOUT)  :: dims
INTEGER,               INTENT(IN)     :: regular_grid

! Declare local variables
REAL(ZREAL8)                          :: level_reg_fact, half_lev_fact, halfdx
INTEGER                               :: i

! Set vertical dimensions
IF (regular_grid .EQ. 1) THEN
  PRINT*, 'Using regular vertical grid'
  level_reg_fact = init_um_data % full_levs(nlevs) / nlevs
  half_lev_fact  = level_reg_fact * 0.5
  halfdx         = dx * 0.5
  DO i = 0, nlevs
     dims % full_levs(i) = REAL(i) * level_reg_fact
  ENDDO

  DO i = 1, nlevs+1
     Dims % half_levs(i) = half_lev_fact + (REAL(i-1) * level_reg_fact)
  ENDDO
  Dims % half_levs(0) = -1. * Dims % half_levs(1)
ELSE
  PRINT*, 'Using MetO model vertical grid'
  Dims % full_levs(0:nlevs)   = init_um_data % full_levs(0:nlevs)
  Dims % half_levs(1:nlevs+1) = init_um_data % half_levs(1:nlevs+1)
  Dims % half_levs(0)         = -1. * Dims % half_levs(1)
END IF

! Set horizontal dimension
DO i = 0,nlongs+1
  Dims % longs_v(i) = REAL(i) * dx
ENDDO

DO i = 0,nlongs+1
  Dims % longs_u(i) = REAL(i) * dx + halfdx
ENDDO

END SUBROUTINE Set_grid


!===================================================================================================
SUBROUTINE Read_um_data_2d (umdata, filename, latitude)

!**********************************
!* Subroutine to read a latitude  *
!* slice of um netcdf output into *
!* a compound data type um_data   *
!*                                *
!* R. Petrie                      *
!* version 2                      *
!* 06/06/2011                     *
!**********************************

USE DefConsTypes, ONLY :  &
    ZREAL8,               &
    um_data_type,         &
    nlongs, nlevs

IMPLICIT NONE
!Include netcdf
!---------------
! Ubuntu:
!INCLUDE '/usr/include/netcdf.inc'
! Carrot:
!INCLUDE '/opt/local/include/netcdf.inc'
! Departmental linux:
INCLUDE '/opt/graphics/64/include/netcdf-4.0.inc'

!Declare parameters
!------------------
TYPE(um_data_type),   INTENT(INOUT)     :: umdata
CHARACTER (LEN=*),    INTENT(IN)        :: filename
INTEGER,              INTENT(IN)        :: latitude

! Declare local variables
!-------------------------
INTEGER             :: ncid, ierr
INTEGER             :: dimidLongs_u, dimidLongs_v, dimidhalf_levs, dimidfull_levs
INTEGER             :: varidLongs_u, varidLongs_v, varidhalf_levs, varidfull_levs
INTEGER             :: varidu, varidv, varidw, variddensity, varidtheta
INTEGER             :: varidorog, varidexpres
INTEGER             :: startA(1), countA(1), startB(4), countB(4),z
INTEGER             :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
REAL (ZREAL8)       :: temp(360)

! Open the netCDF file
!----------------------
ierr = NF_OPEN(filename, NF_NOWRITE, ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, ierr, NF_STRERROR(ierr)
  STOP
ENDIF

!Get the dimension ids
!-----------------------
ierr = NF_INQ_DIMID(ncid, 'x', dimidLongs_u)
ierr1 = ierr
ierr = NF_INQ_DIMID(ncid, 'x_1', dimidLongs_v)
ierr2 = ierr
ierr = NF_INQ_DIMID(ncid, 'hybrid_ht_3', dimidhalf_levs)
ierr3 = ierr
ierr = NF_INQ_DIMID(ncid, 'hybrid_ht_2', dimidfull_levs)
ierr4 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
! PRINT*, ' Dimension ids ok'
ENDIF

! Get the dimension variable ids
!---------------------------------
ierr = NF_INQ_VARID(ncid, 'x', varidLongs_u)
ierr1 = ierr
ierr = NF_INQ_VARID(ncid, 'x_1', varidLongs_v)
ierr2 = ierr
ierr = NF_INQ_VARID(ncid, 'hybrid_ht_3', varidhalf_levs)
ierr3 = ierr
ierr = NF_INQ_VARID(ncid, 'hybrid_ht_2', varidfull_levs)
ierr4 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
!ELSE
!  PRINT*, ' Dimension variable ids ok'
ENDIF


! Get the variable ids
!----------------------
ierr = NF_INQ_VARID(ncid, 'u', varidu)
ierr1 = ierr
ierr = NF_INQ_VARID(ncid, 'v', varidv)
ierr2 = ierr
ierr = NF_INQ_VARID(ncid, 'dz_dt', varidw)
ierr3 = ierr
ierr = NF_INQ_VARID(ncid, 'unspecified', variddensity)
ierr4 = ierr
ierr = NF_INQ_VARID(ncid, 'theta', varidtheta)
ierr5 = ierr
ierr = NF_INQ_VARID(ncid, 'ht', varidorog)
ierr6 = ierr
ierr = NF_INQ_VARID(ncid, 'field7', varidexpres)
ierr7 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) .OR. (ierr5 .NE. 0)&
      .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
  PRINT*, '***Error getting variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
  STOP
!ELSE
!  PRINT*, ' Variable ids ok'
ENDIF

! Get the dimension data from the file
!---------------------------------------
! Longitudes
!-------------
startA(1) = 1
countA(1) = nlongs

ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, umdata % longs_u(1:nlongs))
ierr2 = ierr
ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, umdata % longs_v(1:nlongs))
ierr3 = ierr

! Level Heights
!----------------
startA(1) = 1
countA(1) = nlevs+1
ierr = NF_GET_VARA_DOUBLE (ncid, varidfull_levs, startA, countA, umdata % full_levs(0:nlevs))
ierr4 = ierr
ierr = NF_GET_VARA_DOUBLE (ncid, varidhalf_levs, startA, countA, umdata % half_levs(1:nlevs+1))
ierr5 = ierr

IF ( (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
     (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
! ELSE
!  PRINT*, ' Dimension data ok'
ENDIF

! Get the main data from the file
!----------------------------------
startB(1) = 1
countB(1) = 360       ! All longitude points
startB(2) = latitude  ! Selected latitude slice
countB(2) = 1         !lats
startB(3) = 1
countB(3) = 1         !levs
startB(4) = 1
countB(4) = 1         !time

! u
!----
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidu, startB, countB,  temp)
  umdata % u(1:nlongs,z) = temp(1:nlongs)
END DO
ierr1 = ierr

! v
!----
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidv, startB, countB, temp)
  umdata % v(1:nlongs, z) = temp(1:nlongs)
END DO
ierr2 = ierr

! w
!----
DO z=1, nlevs+1
  startB(3) = z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidw, startB, countB, temp)
  umdata % w(1:nlongs, z-1) = temp(1:nlongs)
ENDDO
ierr3 = ierr

! density
!---------
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, variddensity, startB, countB, temp)
  umdata % density(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr4 = ierr

! theta
!-------
DO z=1, nlevs
  startB(3)= z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidtheta, startB, countB, temp)
  umdata % theta(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr5 = ierr

! exner presure
!----------------
DO z=1, nlevs + 1
  startB(3) = z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidexpres, startB, countB, temp)
  umdata % exner_pressure(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr6 = ierr

! orographic height
!--------------------
DO z=1,1
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidorog, startB, countB, umdata % orog_height(1:nlongs))
ENDDO
ierr7 = ierr

IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
     (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
  PRINT*, '***Error getting main variable data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
  STOP
! ELSE
!  PRINT*, ' Main variable data ok'
ENDIF

!Close the netCDF file
ierr = NF_CLOSE(ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, '***ERROR closing file***'
  STOP
ENDIF

END SUBROUTINE Read_um_data_2d



!===================================================================================================
SUBROUTINE Set_ht_dep_cons (dims)

!********************************
!* Subroutine to set some height*
!* dependent constants          *
!* R. Bannister                 *
!* version 1.0                  *
!* 20/02/14                     *
!********************************

USE DefConsTypes, ONLY :   &
    dimensions_type,       &
    nlevs

IMPLICIT NONE

! Declare Top level variables
TYPE(dimensions_type), INTENT(INOUT)  :: dims

! Declare local variables
INTEGER                               :: z

! For full-to-half level interpolation
DO z = 1, nlevs+1
   dims % a1(z) = (dims % half_levs(z) - dims % full_levs(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
   dims % b1(z) = (Dims % full_levs(z) - Dims % half_levs(z))   / (dims % full_levs(z) - dims % full_levs(z-1))
ENDDO

! For half-to-full level interpolation
DO z = 0, nlevs
   dims % a2(z) = (dims % full_levs(z)   - dims % half_levs(z)) / (dims % half_levs(z+1) - dims % half_levs(z))
   dims % b2(z) = (dims % half_levs(z+1) - dims % full_levs(z)) / (dims % half_levs(z+1) - dims % half_levs(z))
ENDDO

! Reciprocal of some level differences used commonly
DO z = 1, nlevs+1
  dims % recip_half_k_km1(z) = 1.0 / ( dims % half_levs(z) - dims % half_levs(z-1) )
  dims % recip_full_k_km1(z) = 1.0 / ( dims % full_levs(z) - dims % full_levs(z-1) )
END DO
DO z = 0, nlevs
  dims % recip_half_kp1_k(z) = 1.0 / ( dims % half_levs(z+1) - dims % half_levs(z) )
  dims % recip_full_kp1_k(z) = 1.0 / ( dims % full_levs(z+1) - dims % full_levs(z) )
END DO


END SUBROUTINE Set_ht_dep_cons


