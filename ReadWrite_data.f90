!  File includes

!  SUBROUTINE Read_state_2d
!  SUBROUTINE Write_state_2d


!==========================================================================
SUBROUTINE Read_state_2d (filename, state, Dims, t)

!********************************************************
!* Subroutine to read a field of 2d toy model data      *
!*        of model_vars_type                            *
!*                                                      *
!* t is the read time    t                              *
!*                                                      *
!* R. Petrie, 2.0 10-6-11                               *
!* R. Petrie, 3.0 30-7-13                               *
!*                                                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
  model_vars_type,         &
  Dimensions_type,         &
  nlongs,                  &
  nlevs

IMPLICIT NONE

!NetCDF library (file format used to read/write data)
!----------------------------------------------------
! Ubuntu:
!INCLUDE '/usr/include/netcdf.inc'
! Carrot:
!INCLUDE '/opt/local/include/netcdf.inc'
! Departmental linux:
INCLUDE '/opt/graphics/64/include/netcdf-4.0.inc'

!Declare parameters
!------------------
TYPE(model_vars_type),  INTENT(INOUT) :: state
TYPE(Dimensions_type),  INTENT(INOUT) :: Dims
INTEGER,                INTENT(IN)    :: t
CHARACTER(LEN=*)                      :: filename

!Declare local variables
!-----------------------
INTEGER         :: ncid, ierr, dd(3)
INTEGER         :: dimidLongs_u, varidLongs_u
INTEGER         :: dimidLongs_v, varidLongs_v
INTEGER         :: dimidHalf_levs, varidHalf_levs
INTEGER         :: dimidFull_levs, varidFull_levs
INTEGER         :: dimid_time, varid_time
INTEGER         :: dimidEnergy, varidEnergy
INTEGER         :: varid_u, varid_v, varid_w
INTEGER         :: varid_b, varid_beff
INTEGER         :: varid_r, varid_rho
INTEGER         :: varid_tracer, varid_hydro, varid_geost
INTEGER         :: varid_ke, varid_be, varid_ee, varid_te
INTEGER         :: startA(1), countA(1), startC(3), countC(3)
INTEGER         :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
INTEGER         :: ierr8, ierr9, ierr10, ierr11, ierr12, ierr13, ierr14
INTEGER         :: ierr15, ierr16
INTEGER         :: z
REAL(ZREAL8)    :: temp(360), tempvert(60)

!PRINT*,'read field called'

!Open the netCDF file
!---------------------
ierr = NF_OPEN(filename, NF_NOWRITE, ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, ierr, NF_STRERROR(ierr)
  PRINT*, 'FILE :: ', filename
  STOP
ENDIF

!Get the dimension ids
!---------------------
ierr1 = NF_INQ_DIMID(ncid, 'longs_u', dimidLongs_u)
ierr2 = NF_INQ_DIMID(ncid, 'longs_v', dimidLongs_v)
ierr3 = NF_INQ_DIMID(ncid, 'half_level', dimidHalf_levs)
ierr4 = NF_INQ_DIMID(ncid, 'full_level', dimidFull_levs)

IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
  PRINT*, '***Error getting dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
!  PRINT*, ' Dimension ids ok'
ENDIF

! Get the variable ids, for dimensions
!-------------------------------------
ierr1 = NF_INQ_VARID(ncid, 'longs_u', varidLongs_u)
ierr2 = NF_INQ_VARID(ncid, 'longs_v', varidLongs_v)
ierr3 = NF_INQ_VARID(ncid, 'half_level', varidHalf_levs)
ierr4 = NF_INQ_VARID(ncid, 'full_level', varidFull_levs)

IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
  PRINT*, '***Error getting variable dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
!  PRINT*, ' Variable dimension ids ok'
ENDIF


!Get the main variable ids
!-------------------------
ierr1  = NF_INQ_VARID(ncid, 'u', varid_u)
ierr2  = NF_INQ_VARID(ncid, 'v', varid_v)
ierr3  = NF_INQ_VARID(ncid, 'w', varid_w)
ierr4  = NF_INQ_VARID(ncid, 'r_prime', varid_r)
ierr5  = NF_INQ_VARID(ncid, 'b_prime', varid_b)
ierr6  = NF_INQ_VARID(ncid, 'rho', varid_rho)
ierr7  = NF_INQ_VARID(ncid, 'b_effective', varid_beff)
ierr8  = NF_INQ_VARID(ncid, 'tracer', varid_tracer)
ierr9  = NF_INQ_VARID(ncid, 'geo_imbal', varid_geost)
ierr10 = NF_INQ_VARID(ncid, 'hydro_imbal', varid_hydro)
ierr11 = NF_INQ_VARID(ncid, 'ke', varid_ke)
ierr12 = NF_INQ_VARID(ncid, 'be', varid_be)
ierr13 = NF_INQ_VARID(ncid, 'ee', varid_ee)
ierr14 = NF_INQ_VARID(ncid, 'te', varid_te)

IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
      ierr9 + ierr10 + ierr11 + ierr12 + ierr13 + ierr14 ) /= 0) THEN
  PRINT*, '***Error getting variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr7, NF_STRERROR(ierr7)
  PRINT*,'ierr8', ierr8, NF_STRERROR(ierr8)
  PRINT*,'ierr9', ierr9, NF_STRERROR(ierr9)
  PRINT*,'ierr10', ierr10, NF_STRERROR(ierr10)
  PRINT*,'ierr11', ierr11, NF_STRERROR(ierr11)
  PRINT*,'ierr12', ierr12, NF_STRERROR(ierr12)
  PRINT*,'ierr13', ierr13, NF_STRERROR(ierr13)
  PRINT*,'ierr14', ierr14, NF_STRERROR(ierr14)
  STOP
ELSE
!  PRINT*, ' Variable ids are ok'
ENDIF


!Get the dimension data from the file
!------------------------------------
startA(1) = 1
countA(1) = nlongs
ierr1 = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, Dims % longs_u(1:nlongs))
ierr2 = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, Dims % longs_v(1:nlongs))

countA(1) = nlevs
ierr3 = NF_GET_VARA_DOUBLE (ncid, varidHalf_levs, startA, countA, Dims % half_levs(1:nlevs))
ierr4 = NF_GET_VARA_DOUBLE (ncid, varidFull_levs, startA, countA, Dims % full_levs(1:nlevs))

Dims % half_levs(0) = - Dims % half_levs(1)
Dims % full_levs(0) = 0.0
Dims % half_levs(nlevs+1) = 2. * Dims % half_levs(nlevs) - Dims % half_levs(nlevs - 1)
Dims % full_levs(nlevs+1) = 2. * Dims % full_levs(nlevs) - Dims % full_levs(nlevs - 1)

IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
  PRINT*, '*** Error getting dimension data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
!  PRINT*, ' Dimension data ok'
ENDIF


!Get the main data from the file
!-------------------------------
startC(1) = 1
countC(1) = nlongs !longs
startC(2) = 1
countC(2) = 1      !levs
startC(3) = t
countC(3) = 1      !time

! Read 2D variables
DO z = 1, nlevs

  startC(2) = z

  ! Read u
  !-------
  ierr1 = NF_GET_VARA_DOUBLE (ncid, varid_u, startC, countC,  temp)
  state % u(1:nlongs,z) = temp(1:nlongs)

  ! Read v
  !-------
  ierr2 = NF_GET_VARA_DOUBLE (ncid, varid_v, startC, countC, temp)
  state % v(1:nlongs, z) = temp(1:nlongs)

  ! Read w
  !-------
  ierr3 = NF_GET_VARA_DOUBLE (ncid, varid_w, startC, countC, temp)
  state % w(1:nlongs, z) = temp(1:nlongs)

  ! Read r
  !--------
  ierr4 = NF_GET_VARA_DOUBLE (ncid, varid_r, startC, countC, temp)
  state % r(1:nlongs, z) = temp(1:nlongs)

  ! Read bp
  !--------
  ierr5 = NF_GET_VARA_DOUBLE (ncid, varid_b, startC, countC, temp)
  state % b(1:nlongs, z) = temp(1:nlongs)

  ! Read rho
  !---------
  ierr6 = NF_GET_VARA_DOUBLE (ncid, varid_rho, startC, countC, temp)
  state % rho (1:nlongs,z) = temp(1:nlongs)

  ! Read beffective
  !----------------
  ierr7 = NF_GET_VARA_DOUBLE (ncid, varid_beff, startC, countC, temp)
  state % b_ef (1:nlongs,z) = temp(1:nlongs)

  ! Read tracer
  !------------
  ierr8 = NF_GET_VARA_DOUBLE (ncid, varid_tracer, startC, countC, temp)
  state % tracer(1:nlongs, z) = temp(1:nlongs)

  ! Read geostrohpic imbalance
  !---------------------------
  ierr9 = NF_GET_VARA_DOUBLE (ncid, varid_geost, startC, countC, temp)
  state % geost_imbal(1:nlongs, z) = temp(1:nlongs)

  ! Read hydrostatic imbalance
  !---------------------------
  ierr10 = NF_GET_VARA_DOUBLE (ncid, varid_hydro, startC, countC, temp)
  state % hydro_imbal(1:nlongs, z) = temp(1:nlongs)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
       ierr9 + ierr10  ) /= 0) THEN
    PRINT*, '*** Error getting main variable data ***'
    PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
    PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
    PRINT*,'ierr7', ierr7, NF_STRERROR(ierr7)
    PRINT*,'ierr8', ierr8, NF_STRERROR(ierr8)
    PRINT*,'ierr9', ierr9, NF_STRERROR(ierr9)
    STOP
  ENDIF

ENDDO

! Read energy
!------------
startA(1) = t
countA(1) = 1
ierr1 = NF_GET_VARA_DOUBLE (ncid, varid_ke, startA, countA, state % Kinetic_Energy)
ierr2 = NF_GET_VARA_DOUBLE (ncid, varid_be, startA, countA, state % Buoyant_Energy)
ierr3 = NF_GET_VARA_DOUBLE (ncid, varid_ee, startA, countA, state % Elastic_Energy)
ierr4 = NF_GET_VARA_DOUBLE (ncid, varid_te, startA, countA, state % Total_Energy)

IF (( ierr1 + ierr2 + ierr3 + ierr4 ) /= 0) THEN
  PRINT*, '*** Error getting main variable data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
ENDIF

!Close the netCDF file
ierr = NF_CLOSE(ncid)
IF ( ierr .NE. 0 ) THEN
   PRINT*, '***ERROR closing file***'
   STOP
ENDIF

END SUBROUTINE Read_state_2d
!====================================================================================================




!================================================================================
SUBROUTINE Write_state_2d (filename, state, Dims, ntimes, t, wts)

!********************************************************
!* Subroutine to write a 2d field of model_vars_type    *
!*                                                      *
!* ntimes is the total number of output times           *
!*    should include one extra for inital conds         *
!* t is the time of output                              *
!*                                                      *
!* wts is the number of timesteps between output times  *
!*                                                      *
!*                                                      *
!* R. Petrie, 2.0 10-6-11                               *
!* R. Petrie, 3.0 30-7-13                               *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
  model_vars_type,         &
  Dimensions_type,         &
  nlongs,                  &
  nlevs,                   &
  dt

IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
! Ubuntu:
!INCLUDE '/usr/include/netcdf.inc'
! Carrot:
!INCLUDE '/opt/local/include/netcdf.inc'
! Departmental linux:
INCLUDE '/opt/graphics/64/include/netcdf-4.0.inc'

!Declare parameters
!------------------
TYPE(model_vars_type), INTENT(INOUT)   :: state
TYPE(Dimensions_type), INTENT(INOUT)   :: Dims
CHARACTER(LEN=*)                       :: filename
INTEGER                                :: ntimes, t, wts

!Declare local variables
!------------------------
INTEGER, ALLOCATABLE                   :: times(:)
INTEGER                                :: ncid, ierr, i, ddA(1), ddC(3)
INTEGER                                :: dimidLongs_u, varidLongs_u
INTEGER                                :: dimidLongs_v, varidLongs_v
INTEGER                                :: dimidHalf_levs, varidHalf_levs
INTEGER                                :: dimidFull_levs, varidFull_levs
INTEGER                                :: dimid_time, varid_time
INTEGER                                :: dimidEnergy, varidEnergy
INTEGER, SAVE                          :: varid_u, varid_v, varid_w
INTEGER, SAVE                          :: varid_b, varid_beff
INTEGER, SAVE                          :: varid_r, varid_rho
INTEGER, SAVE                          :: varid_tracer, varid_hydro, varid_geost, varid_wmomsource
INTEGER, SAVE                          :: varid_ke, varid_be, varid_ee, varid_te
INTEGER                                :: startA(1), countA(1), startC(3), countC(3)
INTEGER                                :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
INTEGER                                :: ierr8, ierr9, ierr10, ierr11, ierr12, ierr13, ierr14
INTEGER                                :: ierr15, ierr16, ierr10a

!*****************************************************************************************
!PRINT*, 'Write 2d data'
!*****************************************************************************************

! At initial time of run create a new NetCDF output file
IF (t .EQ. 0) THEN

  ! Create netCDF file
  !-------------------------------------
  ierr = NF_CREATE(filename, NF_CLOBBER, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error creating file ***'
    PRINT*, filename
    PRINT*, ierr, NF_STRERROR(ierr)
    STOP
    !ELSE
    !PRINT*, 'FILE CREATED'
  ENDIF

  ! Allocate array for time stamps
  IF (ntimes .EQ. 1) THEN
    ALLOCATE (times(1))  ! output times in seconds
    times(1) = 0
  ELSE
    ALLOCATE (times(0:ntimes-1))  ! output times in seconds
    DO i = 0, ntimes-1
      times(i) = i * wts * dt
    END DO
  END IF

  !Define the dimensions
  !---------------------
  ierr1 = NF_DEF_DIM(ncid, 'longs_u', nlongs, dimidLongs_u)
  ierr2 = NF_DEF_DIM(ncid, 'longs_v', nlongs, dimidLongs_v)
  ierr3 = NF_DEF_DIM(ncid, 'full_level', nlevs, dimidFull_levs)
  ierr4 = NF_DEF_DIM(ncid, 'half_level', nlevs, dimidHalf_levs)
  ierr5 = NF_DEF_DIM(ncid, 'time', ntimes, dimid_time)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 ) /= 0) THEN
     PRINT*, '***Error defining dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids defined'
  ENDIF

  !Define the variables (include variables giving the dim. values)
  !---------------------------------------------------------------

  ! Dimension variables

  ddA(1) = dimidLongs_u
  ierr1  = NF_DEF_VAR(ncid, 'longs_u', NF_DOUBLE, 1, ddA, varidLongs_u)

  ddA(1) = dimidLongs_v
  ierr2  = NF_DEF_VAR(ncid, 'longs_v', NF_DOUBLE, 1, ddA, varidLongs_v)

  ddA(1) = dimidFull_levs
  ierr3  = NF_DEF_VAR(ncid, 'full_level', NF_DOUBLE, 1, ddA, varidFull_levs)

  ddA(1) = dimidHalf_levs
  ierr4  = NF_DEF_VAR(ncid, 'half_level', NF_DOUBLE, 1, ddA, varidHalf_levs)

  ddA(1) = dimid_time
  ierr5  = NF_DEF_VAR(ncid, 'time', NF_DOUBLE, 1, ddA, varid_time)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5) /= 0) THEN
     PRINT*, '***Error defining dimension variable ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     STOP
  !ELSE
  !PRINT*, 'Dimension variable ids defined'
  ENDIF

  ! Main variables

  ddC(1)  = dimidLongs_u
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr1   = NF_DEF_VAR(ncid, 'u', NF_DOUBLE, 3, ddC, varid_u)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr2   = NF_DEF_VAR(ncid, 'v', NF_DOUBLE, 3, ddC, varid_v)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidFull_levs
  ddC(3)  = dimid_time
  ierr3   = NF_DEF_VAR(ncid, 'w', NF_DOUBLE, 3, ddC, varid_w)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr4   = NF_DEF_VAR(ncid, 'r_prime', NF_DOUBLE, 3, ddC, varid_r)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidFull_levs
  ddC(3)  = dimid_time
  ierr5   = NF_DEF_VAR(ncid, 'b_prime', NF_DOUBLE, 3, ddC, varid_b)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr6   = NF_DEF_VAR(ncid, 'rho', NF_DOUBLE, 3, ddC, varid_rho)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidFull_levs
  ddC(3)  = dimid_time
  ierr7   = NF_DEF_VAR(ncid, 'b_effective', NF_DOUBLE, 3, ddC, varid_beff)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr8   = NF_DEF_VAR(ncid, 'tracer', NF_DOUBLE, 3, ddC, varid_tracer)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidHalf_levs
  ddC(3)  = dimid_time
  ierr9   = NF_DEF_VAR(ncid, 'geo_imbal', NF_DOUBLE, 3, ddC, varid_geost)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidFull_levs
  ddC(3)  = dimid_time
  ierr10  = NF_DEF_VAR(ncid, 'hydro_imbal', NF_DOUBLE, 3, ddC, varid_hydro)

  ddC(1)  = dimidLongs_v
  ddC(2)  = dimidFull_levs
  ddC(3)  = dimid_time
  ierr10a = NF_DEF_VAR(ncid, 'wmom_source', NF_DOUBLE, 3, ddC, varid_wmomsource)

  ddA(1)  = dimid_time
  ierr11 = NF_DEF_VAR(ncid, 'ke', NF_DOUBLE, 1, ddA, varid_ke)
  ierr12 = NF_DEF_VAR(ncid, 'be', NF_DOUBLE, 1, ddA, varid_be)
  ierr13 = NF_DEF_VAR(ncid, 'ee', NF_DOUBLE, 1, ddA, varid_ee)
  ierr14 = NF_DEF_VAR(ncid, 'te', NF_DOUBLE, 1, ddA, varid_te)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
      ierr9 + ierr10 + ierr10a + ierr11 + ierr12 + ierr13 + ierr14 ) /= 0) THEN
    PRINT*, '***Error defining main variable ids ***'
    PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
    PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
    PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
    PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
    PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
    PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
    PRINT*, 'ierr7 ',  ierr7,  NF_STRERROR(ierr7)
    PRINT*, 'ierr8 ',  ierr8,  NF_STRERROR(ierr8)
    PRINT*, 'ierr9 ',  ierr9,  NF_STRERROR(ierr9)
    PRINT*, 'ierr10',  ierr10, NF_STRERROR(ierr10)
    PRINT*, 'ierr10a', ierr10, NF_STRERROR(ierr10a)
    PRINT*, 'ierr11',  ierr11, NF_STRERROR(ierr11)
    PRINT*, 'ierr12',  ierr12, NF_STRERROR(ierr12)
    PRINT*, 'ierr13',  ierr13, NF_STRERROR(ierr13)
    PRINT*, 'ierr14',  ierr14, NF_STRERROR(ierr14)
    STOP
  ENDIF

  !Change mode of netCDF operation from define to write
  !------------------------------------------------------
  ierr = NF_ENDDEF(ncid)

  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error changing mode of netCDF operation ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  !ELSE
  !PRINT*, 'Mode Changed'
  END IF

  !--------------------------------------------
  ! Output the values of the dimension variables
  ! --------------------------------------------
  startA(1) = 1
  countA(1) = nlongs
  ierr1     = NF_PUT_VARA_DOUBLE(ncid, varidLongs_u, startA, countA, Dims % longs_u(1:nlongs))
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varidLongs_v, startA, countA, Dims % longs_v(1:nlongs))
  countA(1) = nlevs
  ierr3     = NF_PUT_VARA_DOUBLE(ncid, varidHalf_levs, startA, countA, Dims % half_levs(1:nlevs))
  ierr4     = NF_PUT_VARA_DOUBLE(ncid, varidFull_levs, startA, countA, Dims % full_levs(1:nlevs))

  IF (ntimes .EQ. 1) THEN
    countA(1) = 1
    ierr5 = NF_PUT_VARA_INT(ncid, varid_time, startA, countA, times(1) )
  ELSE
    countA(1) = ntimes
    ierr5 = NF_PUT_VARA_INT(ncid, varid_time, startA, countA, times(0:ntimes-1) )
  ENDIF

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 ) /= 0) THEN
    PRINT*, '***Error writing dimension data ***'
    PRINT*, 'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*, 'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*, 'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*, 'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*, 'ierr5', ierr5, NF_STRERROR(ierr5)
    STOP
  !ELSE
    !PRINT*, 'Dimension Variables output ok'
  ENDIF

  DEALLOCATE (times)

ELSE  ! If time not equal to 0

  PRINT*, t
  !PRINT*, 'opening file'
  !Open the existing output file
  ierr = NF_OPEN(filename, NF_WRITE, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, '***ERROR opening file***', filename
    STOP
  ENDIF
ENDIF

PRINT*, 'beginning output'


!--------------------------------------------
! Output the values of the main variables
! -------------------------------------------

startC(1) = 1
countC(1) = nlongs
startC(2) = 1
countC(2) = nlevs
startC(3) = t+1 ! plus 1 to account for initial time step
countC(3) = 1

ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_u,          startC, countC, state % u(1:nlongs, 1:nlevs))
ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_v,          startC, countC, state % v(1:nlongs, 1:nlevs))
ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_w,          startC, countC, state % w(1:nlongs,1:nlevs))
ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_r,          startC, countC, state % r (1:nlongs, 1:nlevs))
ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_b,          startC, countC, state % b(1:nlongs, 1:nlevs))
ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_rho,        startC, countC, state % rho (1:nlongs, 1:nlevs))
ierr7     = NF_PUT_VARA_DOUBLE(ncid, varid_beff,       startC, countC, state % b_ef(1:nlongs, 1:nlevs))
ierr8     = NF_PUT_VARA_DOUBLE(ncid, varid_tracer,     startC, countC, state % tracer(1:nlongs, 1:nlevs))
ierr9     = NF_PUT_VARA_DOUBLE(ncid, varid_hydro,      startC, countC, state % hydro_imbal(1:nlongs, 1:nlevs))
ierr10    = NF_PUT_VARA_DOUBLE(ncid, varid_geost,      startC, countC, state % geost_imbal(1:nlongs, 1:nlevs))
ierr10a   = NF_PUT_VARA_DOUBLE(ncid, varid_wmomsource, startC, countC, state % vert_mom_source(1:nlongs, 1:nlevs))

IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
  ierr9 + ierr10 + ierr10a) /= 0) THEN
  PRINT*, '***Error writing main variable data ***'
  PRINT*,'ierr1',  ierr1,   NF_STRERROR(ierr1)
  PRINT*,'ierr2',  ierr2,   NF_STRERROR(ierr2)
  PRINT*,'ierr3',  ierr3,   NF_STRERROR(ierr3)
  PRINT*,'ierr4',  ierr4,   NF_STRERROR(ierr4)
  PRINT*,'ierr5',  ierr5,   NF_STRERROR(ierr5)
  PRINT*,'ierr6',  ierr6,   NF_STRERROR(ierr6)
  PRINT*,'ierr7',  ierr7,   NF_STRERROR(ierr7)
  PRINT*,'ierr8',  ierr8,   NF_STRERROR(ierr8)
  PRINT*,'ierr9',  ierr9,   NF_STRERROR(ierr9)
  PRINT*,'ierr10', ierr10,  NF_STRERROR(ierr10)
  PRINT*,'ierr10', ierr10a, NF_STRERROR(ierr10)
  STOP
!ELSE
  !PRINT*, 'Main data written'
ENDIF

startA(1) = t+1 ! plus 1 to account for initial time step
countA(1) = 1
ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_ke, startA, countA, state % Kinetic_energy)
ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_be, startA, countA, state % Buoyant_energy)
ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_ee, startA, countA, state % Elastic_energy)
ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_te, startA, countA, state % Total_energy)

IF ((ierr1 + ierr2 + ierr3 + ierr4 ) /= 0) THEN
  PRINT*, '***Error writing scalars and 1D main data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
  !PRINT*, 'Scalars and 1D main data data written'
ENDIF


!Close-up the file
!-----------------
ierr = NF_CLOSE(ncid)

IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error closing netCDF file ***'
  PRINT*,'ierr', ierr, NF_STRERROR(ierr)
  PRINT*,'xconv -i ', filename, ' &'
  STOP
!ELSE
  !PRINT*, 'File closed'
ENDIF


END SUBROUTINE Write_state_2d
!===================================================================================================
