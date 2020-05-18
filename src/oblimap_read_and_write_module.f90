! File name: oblimap_read_and_write_module.f90
!
! Copyright (C) 2016 Thomas Reerink.
!
! This file is distributed under the terms of the
! GNU General Public License.
!
! This file is part of OBLIMAP 2.0
!
! See Reerink et al. (2010,2016) for OBLIMAP's scientific documentation:
!  http://www.geosci-model-dev.net/3/13/2010/
!  http://www.geosci-model-dev.net/9/4111/2016/
!
! The OBLIMAP User Guide (Reerink, 2016) can be found at:
!  https://github.com/oblimap/oblimap-2.0/tree/master/documentation
!
! The OBLIMAP code can be downloaded by:
!  svn checkout https://svn.science.uu.nl/repos/project.oblimap
! or from OBLIMAP's Github by:
!  git clone https://github.com/oblimap/oblimap-2.0
!
! OBLIMAP is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! OBLIMAP is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with OBLIMAP. If not, see <http://www.gnu.org/licenses/>.
!
!
! OBLIMAP is maintained by:
!
! Thomas Reerink
! Institute for Marine and Atmospheric Research Utrecht (IMAU)
! Utrecht University
! Princetonplein 5
! 3584 CC Utrecht
! The Netherlands
!
! email: tjreerink@gmail.com
!

MODULE oblimap_read_and_write_module
  USE oblimap_configuration_module, ONLY: dp
  IMPLICIT NONE

  TYPE oblimap_netcdf_file_type
    ! This TYPE contains the meta data and the data of the dimensions and the variables:
    INTEGER                                       :: ncid                             ! The id of the netcdf file:  file_id
    INTEGER                                       :: number_of_fields                 ! The number of fields
    CHARACTER(LEN=256)                            :: file_name                        ! The name of the netcdf file
    LOGICAL                                       :: include_time_dimension           ! If .TRUE. the netcdf file contains time records
    LOGICAL                                       :: enable_ignore_option             ! If .TRUE. the ignore option (concerning the reading of the pre-mapped variables) will be enabled
    LOGICAL           , DIMENSION(:), ALLOCATABLE :: ignore_reading_pre_mapped_fields ! The choice whether the pre-mapped variable should be ignored while reading
    INTEGER           , DIMENSION(:), ALLOCATABLE :: spatial_dimension_of_field       ! The spatial dimension of the variable
    CHARACTER(LEN=128), DIMENSION(:), ALLOCATABLE :: field_name                       ! The name of the variable
    CHARACTER(LEN=128), DIMENSION(:), ALLOCATABLE :: field_unit                       ! The units of the variable
    CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE :: field_longname                   ! The long name of the variable (with units)
    INTEGER,            DIMENSION(:), ALLOCATABLE :: id                               ! The ID number of the variable
   !INTEGER,            DIMENSION(:), ALLOCATABLE :: type                             ! An integer determing the variable type: nf90_int, nf90_float or nf90_double
   !INTEGER,            DIMENSION(:), ALLOCATABLE :: case                             ! A number for each situation, for case selection

    INTEGER,            DIMENSION(:), ALLOCATABLE :: LEN_DIM                          ! Contains the dimensional size of each dimension
   !INTEGER                                       :: NDIM                             ! The number of dimensions
   !INTEGER                                       :: NVAR                             ! The number of different variables
   !INTEGER                                       :: N_loop                           ! The number of loops in which the dimensions + variables are written to the netcdf file
   !REAL(dp),           DIMENSION(:), ALLOCATABLE :: grid_size                        ! The grid size for each dimension (Only if it is equidistant)
    REAL(dp),           DIMENSION(:), ALLOCATABLE :: vertical_axis                    ! This array contains the vertical coordinate values
  END TYPE oblimap_netcdf_file_type



CONTAINS
  SUBROUTINE oblimap_open_netcdf_file(file_name, number_of_fields, field_name, LEN_DIM_1, LEN_DIM_2, coordinates_dimension_1, coordinates_dimension_2, enable_ignore_option, nc)
    ! This routine opens a netcdf file.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_open, nf90_nowrite, nf90_inquire, nf90_inq_varid, nf90_inquire_variable, nf90_get_var
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*)                                     , INTENT(IN)            :: file_name                  ! The netcdf file name
    INTEGER                                              , INTENT(IN)            :: number_of_fields           ! Number of Fields
    CHARACTER(LEN=128), DIMENSION(C%MND:number_of_fields), INTENT(IN)            :: field_name                 ! The netcdf field names in the netcdf file
    INTEGER                                              , INTENT(IN)            :: LEN_DIM_1                  ! Number of grid points in longitude direction
    INTEGER                                              , INTENT(IN)            :: LEN_DIM_2                  ! Number of grid points in latitude  direction
    LOGICAL                                              , INTENT(IN) , OPTIONAL :: enable_ignore_option       ! If .TRUE. the ignore option (concerning the reading of the pre-mapped variables) will be enabled

    ! Output variables:
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)   , INTENT(OUT), OPTIONAL :: coordinates_dimension_1
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)   , INTENT(OUT), OPTIONAL :: coordinates_dimension_2
    TYPE(oblimap_netcdf_file_type)                       , INTENT(OUT)           :: nc                         ! The struct which contains the meta data of the netcdf file

    ! Local variables:
    INTEGER                                                                      :: unlimited_dimension_id     ! The ID of the unlimited variable, if there is no unlimited variable this ID is -1
    INTEGER           , DIMENSION(4)                                             :: dimension_ids              ! The dimension ID's (time, NLAT, NLON, NVL) of a variable are gathered
    INTEGER                                                                      :: field_counter              ! The counter in the loop over the field numbers
    INTEGER                                                                      :: number_of_dimensions       ! The number of dimensions of a variable
    INTEGER                                                                      :: k                          ! The counter over the vertical axis

    INTEGER                                                                      :: i, j
    REAL(dp)          , DIMENSION(LEN_DIM_1          )                           :: coordinates_dimension_1_1D ! 1D coordinates for the first  spatial dimension
    REAL(dp)          , DIMENSION(          LEN_DIM_2)                           :: coordinates_dimension_2_1D ! 1D coordinates for the second spatial dimension

    ! Setting the number of fields:
    nc%number_of_fields = number_of_fields

    ALLOCATE(nc%ignore_reading_pre_mapped_fields(C%MND:nc%number_of_fields))
    ALLOCATE(nc%spatial_dimension_of_field      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_name                      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_unit                      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_longname                  (C%MND:nc%number_of_fields))
    ALLOCATE(nc%id                              (C%MND:nc%number_of_fields))
    ALLOCATE(nc%LEN_DIM(3))
    ALLOCATE(nc%vertical_axis(C%number_of_vertical_layers))

    ! The enable_ignore_option is stored in the struct:
    IF(PRESENT(enable_ignore_option)) THEN
     nc%enable_ignore_option = enable_ignore_option
    ELSE
     nc%enable_ignore_option = .FALSE.
    END IF

    ! For each variable the choice whether the pre-mapped variable should be ignored while reading is stored in the struct:
    nc%ignore_reading_pre_mapped_fields = C%ignore_reading_pre_mapped_fields

    ! Default initialization to 2D spatial fields:
    nc%spatial_dimension_of_field(:) = 2

    ! The file name is stored in the struct:
    nc%file_name = file_name

    ! The field names are stored in the struct:
    nc%field_name = field_name

    nc%LEN_DIM(1) = LEN_DIM_1
    nc%LEN_DIM(2) = LEN_DIM_2
    nc%LEN_DIM(3) = C%number_of_vertical_layers                                                ! Number of vertical layers, i.e. number of grid points of the vertical coordinate

    ! Open the netcdf file:
    CALL handle_error(nf90_open(nc%file_name, nf90_nowrite, nc%ncid),                                                       '. [ 1] From oblimap_open_netcdf_file(): The file '//TRIM(nc%file_name)//' is not found')

    CALL handle_error(nf90_inquire(ncid = nc%ncid, unlimitedDimID = unlimited_dimension_id),                                '. [ 2] From oblimap_open_netcdf_file(): it concerns the  file '//TRIM(nc%file_name))
    IF(unlimited_dimension_id == -1) THEN
     nc%include_time_dimension = .FALSE.
    !WRITE(UNIT=*, FMT='(/A/)') ' The netcdf file  '//TRIM(nc%file_name)//'  has no time dimension.'
    ELSE
     nc%include_time_dimension = .TRUE.
    !WRITE(UNIT=*, FMT='(/A/)') ' The netcdf file  '//TRIM(nc%file_name)//'  has a time dimension.'

     IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(0)) THEN
      ! Ignore reading the time axis name
     ELSE
      CALL handle_error(nf90_inq_varid(nc%ncid, nc%field_name(0), nc%id(0)),                                                '. [ 3] From oblimap_open_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(0))//'" in the file '//TRIM(nc%file_name))
     END IF
    END IF

    ! Read the other field variables:
    DO field_counter = 1, nc%number_of_fields
      IF(nc%enable_ignore_option .AND. nc%ignore_reading_pre_mapped_fields(field_counter)) THEN
       ! Ignore reading this pre mapped field.
      ELSE
       CALL handle_error(nf90_inq_varid(nc%ncid, nc%field_name(field_counter), nc%id(field_counter)),                       '. [ 4] From oblimap_open_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', with field number ', field_counter)

       ! Initialize all dimension_ids to a negative value which will be unequal to the time dimension ID:
       dimension_ids = -19062014
       CALL handle_error(nf90_inquire_variable(ncid = nc%ncid, varid = nc%id(field_counter), ndims = number_of_dimensions, &
                                                                                             dimids = dimension_ids),       '. [ 5] From oblimap_open_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name))
       IF(nc%include_time_dimension) THEN
        ! Check whether the considered field is realy time dependent:
        IF(ANY(dimension_ids(:) == unlimited_dimension_id)) THEN
         nc%spatial_dimension_of_field(field_counter) = number_of_dimensions - 1
        ELSE
         nc%spatial_dimension_of_field(field_counter) = number_of_dimensions
        END IF
       ELSE
        nc%spatial_dimension_of_field(field_counter) = number_of_dimensions
       END IF
      END IF
     !WRITE(UNIT=*, FMT='(A, A30, A, I2)') ' The netcdf field  ', TRIM(nc%field_name(field_counter)), '  has a spatial dimension of ', nc%spatial_dimension_of_field(field_counter)
    END DO

    IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(-6)) THEN
     ! Ignore reading the vertical axis
     nc%vertical_axis = (/ (k, k=1, nc%LEN_DIM(3)) /) ! Avoiding the case that nc%vertical_axis is uninitialized but is used by copying the entire nc struct.
    ELSE IF(ANY(nc%spatial_dimension_of_field(:) == 3)) THEN
     CALL handle_error(nf90_inq_varid(nc%ncid, nc%field_name(-6), nc%id(-6)             ),                                  '. [ 6] From oblimap_open_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(nf90_get_var(nc%ncid, nc%id(-6), nc%vertical_axis(:), start=(/1/)),                                  '. [ 7] From oblimap_open_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name)//', field number ', -6)
    ELSE
     nc%vertical_axis = (/ (k, k=1, nc%LEN_DIM(3)) /) ! Avoiding the case that nc%vertical_axis is uninitialized but is used by copying the entire nc struct.
    END IF

    IF(PRESENT(coordinates_dimension_1)) THEN
     CALL handle_error(nf90_inq_varid(nc%ncid, nc%field_name(-2), nc%id(-2)),                                               '. [ 8] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(nf90_inquire_variable(ncid = nc%ncid, varid = nc%id(-2), ndims = nc%spatial_dimension_of_field(-2)), '. [ 9] From oblimap_open_netcdf_file(): it concerns the field "'    //TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))

     ! Read the first spatial dimension variable:
     IF(nc%spatial_dimension_of_field(-2) == 1) THEN
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(-2), coordinates_dimension_1_1D),                                       '. [10] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))

      ! Creating from the one dimensional variable, the two dimensional one:
      DO i = 1, nc%LEN_DIM(1)
       coordinates_dimension_1(i,:) = coordinates_dimension_1_1D(i)
      END DO
     ELSE
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(-2), coordinates_dimension_1),                                          '. [11] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     END IF
     ! Unit conversion, e.g. from kilometers to meters:
     coordinates_dimension_1 = coordinates_dimension_1 * C%unit_conversion_x_ax
    END IF

    IF(PRESENT(coordinates_dimension_2)) THEN
     CALL handle_error(nf90_inq_varid(nc%ncid, nc%field_name(-4), nc%id(-4)),                                               '. [12] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(nf90_inquire_variable(ncid = nc%ncid, varid = nc%id(-4), ndims = nc%spatial_dimension_of_field(-4)), '. [13] From oblimap_open_netcdf_file(): it concerns the field "'    //TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))

     ! Read the second spatial dimension variable:
     IF(nc%spatial_dimension_of_field(-4) == 1) THEN
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(-4), coordinates_dimension_2_1D),                                       '. [14] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))

      ! Creating from the one dimensional variable, the two dimensional one:
      DO j = 1, nc%LEN_DIM(2)
       coordinates_dimension_2(:,j) = coordinates_dimension_2_1D(j)
      END DO
     ELSE
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(-4), coordinates_dimension_2),                                          '. [15] From oblimap_open_netcdf_file(): it concerns the dimension "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     END IF
     ! Unit conversion, e.g. from kilometers to meters:
     coordinates_dimension_2 = coordinates_dimension_2 * C%unit_conversion_y_ax
    END IF
  END SUBROUTINE oblimap_open_netcdf_file



  SUBROUTINE oblimap_create_netcdf_file(file_name, number_of_fields, include_time_dimension, spatial_dimension_of_field, &
                                        field_name, field_unit, field_longname, &
                                        LEN_DIM_1, LEN_DIM_2, coordinates_dimension_1, coordinates_dimension_2, vertical_axis, nc)
    ! This routine creates a netcdf file, the netcdf format is specified. In case the horizontal spatial coordinates are not specified by the optional
    ! arguments, a default 2D horizontal spatial C%NX by C%NY grid is assumed with grid sizes C%dx and C%dy.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_create, nf90_clobber, nf90_def_dim, nf90_def_var, nf90_unlimited, nf90_double, nf90_float, nf90_put_att, nf90_global, nf90_enddef, nf90_put_var
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*)                                           , INTENT(IN)           :: file_name
    INTEGER                                                    , INTENT(IN)           :: number_of_fields
    LOGICAL                                                    , INTENT(IN)           :: include_time_dimension
    INTEGER           , DIMENSION(C%MND:number_of_fields)      , INTENT(IN)           :: spatial_dimension_of_field
    CHARACTER(LEN=128), DIMENSION(C%MND:number_of_fields)      , INTENT(IN)           :: field_name
    CHARACTER(LEN=128), DIMENSION(C%MND:number_of_fields)      , INTENT(IN)           :: field_unit
    CHARACTER(LEN=256), DIMENSION(C%MND:number_of_fields)      , INTENT(IN)           :: field_longname
    INTEGER                                                    , INTENT(IN)           :: LEN_DIM_1
    INTEGER                                                    , INTENT(IN)           :: LEN_DIM_2
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)         , INTENT(IN), OPTIONAL :: coordinates_dimension_1        ! The optional input 2D coordinates for the first  spatial dimension
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)         , INTENT(IN), OPTIONAL :: coordinates_dimension_2        ! The optional input 2D coordinates for the second spatial dimension
    REAL(dp)          , DIMENSION(C%number_of_vertical_layers) , INTENT(IN), OPTIONAL :: vertical_axis                  ! Currently the vertical layers are only numbered, here the real ones could be parsed

    ! Output variables:
    TYPE(oblimap_netcdf_file_type)                             , INTENT(OUT)          :: nc

    ! Local variables:
    INTEGER                                                                           :: field_counter                  ! The counter in the loop over the field numbers
    INTEGER                                                                           :: xtype_in_netcdf                ! Writing the netcdf data in nf90_double or nf90_float precision
    INTEGER                                                                           :: k                              ! The counter over the vertical axis
    INTEGER                                                                           :: m, n                           ! The counter over the 'other' horrizontal spatial coordinate, to check for symmetry
    LOGICAL                                                                           :: include_two_spatial_dimensions
    LOGICAL                                                                           :: include_vertical_dimension
    LOGICAL                                                                           :: coordinates_dimension_1_are_1D
    LOGICAL                                                                           :: coordinates_dimension_2_are_1D
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)                                :: coordinates_dimension_1_2D     ! 2D coordinates for the first  spatial dimension
    REAL(dp)          , DIMENSION(LEN_DIM_1,LEN_DIM_2)                                :: coordinates_dimension_2_2D     ! 2D coordinates for the second spatial dimension
    REAL(dp)          , DIMENSION(LEN_DIM_1          )                                :: coordinates_dimension_1_1D     ! 1D coordinates for the first  spatial dimension
    REAL(dp)          , DIMENSION(          LEN_DIM_2)                                :: coordinates_dimension_2_1D     ! 1D coordinates for the second spatial dimension

    IF(C%use_double_instead_of_float_in_netcdf) THEN
     xtype_in_netcdf = nf90_double
    ELSE
     xtype_in_netcdf = nf90_float
    END IF

    ! Setting the number of fields:
    nc%number_of_fields = number_of_fields

    ALLOCATE(nc%ignore_reading_pre_mapped_fields(C%MND:nc%number_of_fields))
    ALLOCATE(nc%spatial_dimension_of_field      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_name                      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_unit                      (C%MND:nc%number_of_fields))
    ALLOCATE(nc%field_longname                  (C%MND:nc%number_of_fields))
    ALLOCATE(nc%id                              (C%MND:nc%number_of_fields))
    ALLOCATE(nc%LEN_DIM(3))
    ALLOCATE(nc%vertical_axis(C%number_of_vertical_layers))

    ! The file name is stored in the struct:
    nc%file_name = file_name

    ! Whether the fields contain a time dimension is stored in the struct:
    nc%include_time_dimension = include_time_dimension

    ! The spatial dimension of the variables are stored in the struct:
    nc%spatial_dimension_of_field = spatial_dimension_of_field

    ! The field names are stored in the struct:
    nc%field_name     = field_name
    nc%field_unit     = field_unit
    nc%field_longname = field_longname

    nc%LEN_DIM(1) = LEN_DIM_1
    nc%LEN_DIM(2) = LEN_DIM_2

    ! If both optional arguments coordinates_dimension_1 and coordinates_dimension_2 are specified, then this inititialization can be omitted:
    IF(.NOT.(PRESENT(coordinates_dimension_1) .AND. PRESENT(coordinates_dimension_2))) &
     CALL initialize_im_coordinates(C%NX, C%NY, C%dx, C%dy, coordinates_dimension_1_2D, coordinates_dimension_2_2D)

    IF(PRESENT(coordinates_dimension_1)) coordinates_dimension_1_2D = coordinates_dimension_1
    IF(PRESENT(coordinates_dimension_2)) coordinates_dimension_2_2D = coordinates_dimension_2

    coordinates_dimension_1_are_1D = .TRUE.
    DO n = 1, nc%LEN_DIM(2)
     IF(ALL(coordinates_dimension_1_2D(:,1) == coordinates_dimension_1_2D(:,n))) THEN
      ! Keep coordinates_dimension_1_are_1D = .TRUE.
     ELSE
      coordinates_dimension_1_are_1D = .FALSE.
      EXIT
     END IF
    END DO

    coordinates_dimension_2_are_1D = .TRUE.
    DO m = 1, nc%LEN_DIM(1)
     IF(ALL(coordinates_dimension_2_2D(1,:) == coordinates_dimension_2_2D(m,:))) THEN
      ! Keep coordinates_dimension_2_are_1D = .TRUE.
     ELSE
      coordinates_dimension_2_are_1D = .FALSE.
      EXIT
     END IF
    END DO

    ! Reduce the 2D coordinates_dimension_1_2D to 1D in case they are symmetric in coordinates_dimension_2 direction.
    IF(coordinates_dimension_1_are_1D) coordinates_dimension_1_1D = coordinates_dimension_1_2D(:,1)

    ! Reduce the 2D coordinates_dimension_2_2D to 1D in case they are symmetric in coordinates_dimension_1 direction.
    IF(coordinates_dimension_2_are_1D) coordinates_dimension_2_1D = coordinates_dimension_2_2D(1,:)

    IF(ANY(nc%spatial_dimension_of_field(:) == 2 .OR. nc%spatial_dimension_of_field(:) == 3)) THEN
     include_two_spatial_dimensions = .TRUE.
    ELSE
     include_two_spatial_dimensions = .FALSE.
    END IF

    IF(ANY(nc%spatial_dimension_of_field == 3)) THEN
     include_vertical_dimension = .TRUE.
     nc%LEN_DIM(3)              = C%number_of_vertical_layers   ! Number of vertical layers, i.e. number of grid points of the vertical coordinate
    ELSE
     include_vertical_dimension = .FALSE.
     nc%LEN_DIM(3)              = 1                             ! To save a waste of allocation and computational time in case C%number_of_vertical_layers is initiated to a higher number
    END IF

    CALL check_file_existence(nc%file_name)

    ! Create the netcdf file:
    CALL handle_error(nf90_create(nc%file_name, nf90_clobber, nc%ncid))

    ! Define the dimension dimensions:
    IF(nc%include_time_dimension)     &
     CALL handle_error(  nf90_def_dim(nc%ncid, nc%field_name(-1), nf90_unlimited,                             nc%id(-1)), '. [ 1] From oblimap_create_netcdf_file(): it concerns the dimension "'                        //TRIM(nc%field_name(-1))//'" in the file '//TRIM(nc%file_name))
    IF(include_vertical_dimension) &
     CALL handle_error(  nf90_def_dim(nc%ncid, nc%field_name(-7), nc%LEN_DIM(3) ,                             nc%id(-7)), '. [ 2] From oblimap_create_netcdf_file(): it concerns the dimension "'                        //TRIM(nc%field_name(-7))//'" in the file '//TRIM(nc%file_name))
    IF(include_two_spatial_dimensions) &
     CALL handle_error(  nf90_def_dim(nc%ncid, nc%field_name(-5), nc%LEN_DIM(2) ,                             nc%id(-5)), '. [ 3] From oblimap_create_netcdf_file(): it concerns the dimension "'                        //TRIM(nc%field_name(-5))//'" in the file '//TRIM(nc%file_name))
    CALL handle_error(   nf90_def_dim(nc%ncid, nc%field_name(-3), nc%LEN_DIM(1) ,                             nc%id(-3)), '. [ 4] From oblimap_create_netcdf_file(): it concerns the dimension "'                        //TRIM(nc%field_name(-3))//'" in the file '//TRIM(nc%file_name))

    ! Define the dimension variables:
    IF(coordinates_dimension_1_are_1D) THEN
     CALL handle_error(  nf90_def_var(nc%ncid, nc%field_name(-2), xtype_in_netcdf, (/nc%id(-3)           /) , nc%id(-2)), '. [ 5] From oblimap_create_netcdf_file(): it concerns the dimension variable "'               //TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
    ELSE
     CALL handle_error(  nf90_def_var(nc%ncid, nc%field_name(-2), xtype_in_netcdf, (/nc%id(-3), nc%id(-5)/) , nc%id(-2)), '. [ 6] From oblimap_create_netcdf_file(): it concerns the field "'                            //TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
    END IF
    CALL handle_error(   nf90_put_att(nc%ncid, nc%id(-2), 'long_name', nc%field_longname(-2))                           , '. [ 7] From oblimap_create_netcdf_file(): it concerns the long_name definition of the field "'//TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
    CALL handle_error(   nf90_put_att(nc%ncid, nc%id(-2), 'units'    , nc%field_unit(-2))                               , '. [ 8] From oblimap_create_netcdf_file(): it concerns the units definition of the field "'    //TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))

    IF(include_two_spatial_dimensions) THEN
     IF(coordinates_dimension_2_are_1D) THEN
      CALL handle_error( nf90_def_var(nc%ncid, nc%field_name(-4), xtype_in_netcdf, (/           nc%id(-5)/) , nc%id(-4)), '. [ 9] From oblimap_create_netcdf_file(): it concerns the dimension variable "'               //TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     ELSE
      CALL handle_error( nf90_def_var(nc%ncid, nc%field_name(-4), xtype_in_netcdf, (/nc%id(-3), nc%id(-5)/) , nc%id(-4)), '. [10] From oblimap_create_netcdf_file(): it concerns the field "'                            //TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     END IF
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(-4), 'long_name', nc%field_longname(-4))                           , '. [11] From oblimap_create_netcdf_file(): it concerns the long_name definition of the field "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(-4), 'units'    , nc%field_unit(-4))                               , '. [12] From oblimap_create_netcdf_file(): it concerns the units definition of the field "'    //TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
    END IF

    IF(include_vertical_dimension) THEN
     CALL handle_error(  nf90_def_var(nc%ncid, nc%field_name(-6), xtype_in_netcdf, (/nc%id(-7)           /) , nc%id(-6)), '. [13] From oblimap_create_netcdf_file(): it concerns the field "'                            //TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(-6), 'long_name', nc%field_longname(-6))                           , '. [14] From oblimap_create_netcdf_file(): it concerns the long_name definition of the field "'//TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(-6), 'units'    , nc%field_unit(-6))                               , '. [15] From oblimap_create_netcdf_file(): it concerns the units definition of the field "'    //TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
    END IF

    IF(nc%include_time_dimension) THEN
     CALL handle_error(  nf90_def_var(nc%ncid, nc%field_name( 0), xtype_in_netcdf, (/nc%id(-1)           /) , nc%id( 0)), '. [16] From oblimap_create_netcdf_file(): it concerns the dimension variable "'               //TRIM(nc%field_name( 0))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id( 0), 'long_name', nc%field_longname( 0))                           , '. [17] From oblimap_create_netcdf_file(): it concerns the long_name definition of the field "'//TRIM(nc%field_name( 0))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id( 0), 'units'    , nc%field_unit( 0))                               , '. [18] From oblimap_create_netcdf_file(): it concerns the units definition of the field "'    //TRIM(nc%field_name( 0))//'" in the file '//TRIM(nc%file_name))
    END IF

    ! Define the field variables:
    DO field_counter = 1, nc%number_of_fields
     IF(nc%include_time_dimension) THEN
      IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3), nc%id(-5), nc%id(-7), nc%id(-1)/), nc%id(field_counter)), &
        '. [19] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3)                      , nc%id(-1)/), nc%id(field_counter)), &
        '. [20] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3), nc%id(-5)           , nc%id(-1)/), nc%id(field_counter)), &
        '. [21] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      END IF
     ELSE
      IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3), nc%id(-5), nc%id(-7)           /), nc%id(field_counter)), &
        '. [22] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3)                                 /), nc%id(field_counter)), &
        '. [23] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE
       CALL handle_error(nf90_def_var(nc%ncid, nc%field_name(field_counter), xtype_in_netcdf, (/nc%id(-3), nc%id(-5)                      /), nc%id(field_counter)), &
        '. [24] From oblimap_create_netcdf_file(): it concerns the field variable "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      END IF
     END IF

     ! In case the description equals still the initial dummy default, a neater one is generated:
     IF(nc%field_longname(field_counter) == 'longname: ?') nc%field_longname(field_counter) = &
      TRIM(nc%field_name(field_counter))//' ('//TRIM(nc%field_unit(field_counter))//')'
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(field_counter), 'units'        , nc%field_unit         (field_counter)), '. [25] From oblimap_create_netcdf_file(): it concerns the units definition, of the field "'        //TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(field_counter), 'long_name'    , nc%field_longname     (field_counter)), '. [26] From oblimap_create_netcdf_file(): it concerns the long_name definition, of the field "'    //TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name))
     CALL handle_error(  nf90_put_att(nc%ncid, nc%id(field_counter), 'missing_value', C%invalid_output_value(field_counter)), '. [27] From oblimap_create_netcdf_file(): it concerns the missing_value definition, of the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name))
    END DO
    CALL handle_error(   nf90_put_att(nc%ncid, nf90_global         , 'acknowledgment', 'Created with OBLIMAP 2.0, copyright (C) 2016 Thomas Reerink under the GNU General Public License Copyright.'), '. [28] From oblimap_create_netcdf_file(): it concerns the global attribute definition, in the file '//TRIM(nc%file_name))
    CALL handle_error(   nf90_put_att(nc%ncid, nf90_global         , 'web-reference-1', 'OBLIMAP github: https://github.com/oblimap/oblimap-2.0'                      ), '. [29] From oblimap_create_netcdf_file(): it concerns the global attribute definition, in the file '//TRIM(nc%file_name))
    CALL handle_error(   nf90_put_att(nc%ncid, nf90_global         , 'web-reference-2', 'Reerink et al. (2010): http://www.geosci-model-dev.net/3/13/2010/'           ), '. [30] From oblimap_create_netcdf_file(): it concerns the global attribute definition, in the file '//TRIM(nc%file_name))
    CALL handle_error(   nf90_put_att(nc%ncid, nf90_global         , 'web-reference-3', 'Reerink et al. (2016): http://www.geosci-model-dev-discuss.net/gmd-2016-124/'), '. [31] From oblimap_create_netcdf_file(): it concerns the global attribute definition, in the file '//TRIM(nc%file_name))

    CALL handle_error(   nf90_enddef(nc%ncid))


    ! Put the coordinates of dimension 1:
    IF(coordinates_dimension_1_are_1D) THEN
     CALL handle_error(nf90_put_var( nc%ncid, nc%id(-2), coordinates_dimension_1_1D    ), '. [32] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
    ELSE
     CALL handle_error(nf90_put_var( nc%ncid, nc%id(-2), coordinates_dimension_1_2D    ), '. [33] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-2))//'" in the file '//TRIM(nc%file_name))
    END IF

    IF(include_two_spatial_dimensions) THEN
     ! Put the coordinates of dimension 2:
     IF(coordinates_dimension_2_are_1D) THEN
      CALL handle_error(nf90_put_var(nc%ncid, nc%id(-4), coordinates_dimension_2_1D    ), '. [34] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     ELSE
      CALL handle_error(nf90_put_var(nc%ncid, nc%id(-4), coordinates_dimension_2_2D    ), '. [35] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-4))//'" in the file '//TRIM(nc%file_name))
     END IF
    END IF

    ! Put the vertical coordinate variable:
    IF(include_vertical_dimension) THEN
     IF(PRESENT(vertical_axis)) THEN
      CALL handle_error(nf90_put_var(nc%ncid, nc%id(-6), vertical_axis                 ), '. [36] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
     ELSE
      CALL handle_error(nf90_put_var(nc%ncid, nc%id(-6), (/ (k, k=1, nc%LEN_DIM(3)) /) ), '. [37] From oblimap_create_netcdf_file(): it concerns the field "'//TRIM(nc%field_name(-6))//'" in the file '//TRIM(nc%file_name))
     END IF
    END IF
  END SUBROUTINE oblimap_create_netcdf_file



  SUBROUTINE create_netcdf_for_gcm_grid(longitude_coordinates_of_gcm_grid_points, latitude_coordinates_of_gcm_grid_points, nc_for_dimensional_shape, nc)
    ! This routine creates a gcm netcdf file, the netcdf format is specified, the dimensions are written.
    USE oblimap_configuration_module, ONLY: dp, C, check_directory_existence
    IMPLICIT NONE

    ! Input variables:
    REAL(dp)        , DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: longitude_coordinates_of_gcm_grid_points
    REAL(dp)        , DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: latitude_coordinates_of_gcm_grid_points
    TYPE(oblimap_netcdf_file_type)            , INTENT(IN)  :: nc_for_dimensional_shape

    ! Output variables:
    TYPE(oblimap_netcdf_file_type)            , INTENT(OUT) :: nc

    ! Check whether the directory in the path of C%im_created_filename exists:
    CALL check_directory_existence(C%im_created_filename)

    ! Output: nc
    CALL oblimap_create_netcdf_file(C%gcm_created_filename, C%number_of_mapped_fields, &
                                    nc_for_dimensional_shape%include_time_dimension, &
                                    nc_for_dimensional_shape%spatial_dimension_of_field, &
                                    C%gcm_field_name, C%gcm_field_unit, C%gcm_field_longname, &
                                    C%NLON, C%NLAT, longitude_coordinates_of_gcm_grid_points, latitude_coordinates_of_gcm_grid_points, &
                                    vertical_axis = nc_for_dimensional_shape%vertical_axis, nc = nc)
  END SUBROUTINE create_netcdf_for_gcm_grid



  SUBROUTINE create_netcdf_for_im_grid(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, nc_for_dimensional_shape, nc)
    ! This routine creates an im netcdf file, the netcdf format is specified, the dimensions are written.
    USE oblimap_configuration_module, ONLY: dp, C, check_directory_existence
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(IN)  :: x_coordinates_of_im_grid_points
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(IN)  :: y_coordinates_of_im_grid_points
    TYPE(oblimap_netcdf_file_type), INTENT(IN)  :: nc_for_dimensional_shape

    ! Output variables:
    TYPE(oblimap_netcdf_file_type), INTENT(OUT) :: nc

    ! Check whether the directory in the path of C%im_created_filename exists:
    CALL check_directory_existence(C%im_created_filename)

    ! Output: nc
    CALL oblimap_create_netcdf_file(C%im_created_filename, C%number_of_mapped_fields, &
                                    nc_for_dimensional_shape%include_time_dimension, &
                                    nc_for_dimensional_shape%spatial_dimension_of_field, &
                                    C%im_field_name, C%im_field_unit, C%im_field_longname, &
                                    C%NX, C%NY, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, &
                                    vertical_axis = nc_for_dimensional_shape%vertical_axis, nc = nc)
  END SUBROUTINE create_netcdf_for_im_grid



  SUBROUTINE oblimap_read_im_coordinates_from_netcdf_file(coordinates_dimension_1, coordinates_dimension_2)
    ! This routine reads the horizontal spatial coordinates from a netcdf file.
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Output variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: coordinates_dimension_1
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: coordinates_dimension_2

    ! Local variables:
    TYPE(oblimap_netcdf_file_type)              :: nc

    ! The spatial horizontal coordinates of the IM grid are read from a netcdf file:
    ! Output: nc
    CALL oblimap_open_netcdf_file(C%prefabricated_im_grid_filename, 0, C%prefabricated_im_grid_field_name, C%NX, C%NY, coordinates_dimension_1, coordinates_dimension_2, nc = nc)

    ! Output: -
    CALL oblimap_close_netcdf_file(nc)
  END SUBROUTINE oblimap_read_im_coordinates_from_netcdf_file



  SUBROUTINE initialize_im_coordinates(NX, NY, dx, dy, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)
    ! Initialize the x, y coordinate values of the IM grid points in plane S'. In total there are
    ! NX * NY grid points and default the grid cell distances are dx and dy (in meters). The central
    ! grid point is (0,0) and coincide with the longitude-latitude coordinates (lamda_M, phi_M), if not
    ! shifted by the routine shifting_center_im_grid (see oblimap_scan_contributions_module). Mind that
    ! the center (0,0) coincides with the middle point for an odd grid number, and with the lower one
    ! of the two middle points for an even grid number.
    ! Instead of using the default dx and dy, it is possible to read a prefabricated_im_grid (which might
    ! contain the coordinates in a 1D or 2D format). This generates the option to use 2D irregular IM grid
    ! coordinates.
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    INTEGER                   , INTENT(IN)  :: NX
    INTEGER                   , INTENT(IN)  :: NY
    REAL(dp)                  , INTENT(IN)  :: dx
    REAL(dp)                  , INTENT(IN)  :: dy

    ! Output variables:
    REAL(dp), DIMENSION(NX,NY), INTENT(OUT) :: x_coordinates_of_im_grid_points
    REAL(dp), DIMENSION(NX,NY), INTENT(OUT) :: y_coordinates_of_im_grid_points

    ! Local variables:
    INTEGER                                 :: m, n

    IF(C%use_prefabricated_im_grid_coordinates) THEN
     ! The x, y coordinates of the IM grid are read from a separate netcdf file:
     ! Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
     CALL oblimap_read_im_coordinates_from_netcdf_file(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)
    ELSE
     DO m = 1, NX
     DO n = 1, NY
       x_coordinates_of_im_grid_points(m,n) = dx * (m - ((NX+1) / 2.0))
       y_coordinates_of_im_grid_points(m,n) = dy * (n - ((NY+1) / 2.0))
     END DO
     END DO
    END IF
  END SUBROUTINE initialize_im_coordinates



  SUBROUTINE oblimap_read_netcdf_fields(nc, record_number, fields, time)
    ! This routine reads for a given record all the fields, it is used for reading both gcm and im netcdf files.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_get_var
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_netcdf_file_type)                                                    , INTENT(IN)            :: nc
    INTEGER                                                                           , INTENT(IN)            :: record_number    ! time record number

    ! Output variables:
    REAL(dp), DIMENSION(nc%number_of_fields,nc%LEN_DIM(1),nc%LEN_DIM(2),nc%LEN_DIM(3)), INTENT(OUT)           :: fields
    REAL(dp)                                                                          , INTENT(OUT), OPTIONAL :: time

    ! Local variables:
    INTEGER                                                                                                   :: field_counter    ! The counter in the loop over the field numbers

    ! In case time is in argument list, initialize in any case:
    IF(PRESENT(time)) time = 0.0_dp

    IF(PRESENT(time) .AND. nc%include_time_dimension) THEN
     IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(0)) THEN
      ! Ignore reading the time axis
     ELSE
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(0), time, start=(/record_number/)), '. [1] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(0))//'" in the file '//TRIM(nc%file_name))
     END IF
    END IF

    ! Read the field variables:
    DO field_counter = 1, nc%number_of_fields
      IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(field_counter)) THEN
       ! Ignore reading this pre mapped field.
       fields(field_counter,:,:,:) = C%invalid_input_value(field_counter)
      ELSE

       IF(nc%include_time_dimension) THEN
        IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,:), start=(/1, 1, 1, record_number/)), &
          '. [2] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,1,1), start=(/1      , record_number/)), &
          '. [3] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        ELSE
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,1), start=(/1, 1   , record_number/)), &
          '. [4] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        END IF
       ELSE
        IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,:), start=(/1, 1, 1/)), &
          '. [5] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,1,1), start=(/1      /)), &
          '. [6] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        ELSE
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,1), start=(/1, 1   /)), &
          '. [7] From oblimap_read_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
        END IF
       END IF

      END IF
    END DO
  END SUBROUTINE oblimap_read_netcdf_fields



  SUBROUTINE oblimap_read_netcdf_1D_fields(nc, record_number, fields, time)
    ! This routine reads for a given record all the 1D fields.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_get_var
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_netcdf_file_type)                        , INTENT(IN)            :: nc
    INTEGER                                               , INTENT(IN)            :: record_number    ! time record number

    ! Output variables:
    REAL(dp), DIMENSION(nc%number_of_fields,nc%LEN_DIM(1)), INTENT(OUT)           :: fields
    REAL(dp)                                              , INTENT(OUT), OPTIONAL :: time

    ! Local variables:
    INTEGER                                                                       :: field_counter    ! The counter in the loop over the field numbers

    ! In case time is in argument list, initialize in any case:
    IF(PRESENT(time)) time = 0.0_dp

    IF(PRESENT(time) .AND. nc%include_time_dimension) THEN
     IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(0)) THEN
      ! Ignore reading the time axis
     ELSE
      CALL handle_error(nf90_get_var(nc%ncid, nc%id(0), time, start=(/record_number/)), '. [1] From oblimap_read_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(0))//'" in the file '//TRIM(nc%file_name))
     END IF
    END IF

    ! Read the field variables:
    DO field_counter = 1, nc%number_of_fields
      IF(nc%enable_ignore_option .AND. C%ignore_reading_pre_mapped_fields(field_counter)) THEN
       ! Ignore reading this pre mapped field.
       fields(field_counter,:) = C%invalid_input_value(field_counter)
      ELSE

       IF(nc%include_time_dimension) THEN
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:), start=(/1, record_number/)), &
          '. [2] From oblimap_read_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
       ELSE
         CALL handle_error(nf90_get_var(nc%ncid, nc%id(field_counter), fields(field_counter,:), start=(/1               /)), &
          '. [3] From oblimap_read_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
       END IF
      END IF
    END DO
  END SUBROUTINE oblimap_read_netcdf_1D_fields



  SUBROUTINE oblimap_write_netcdf_fields(nc, record_number, fields, time)
    ! This routine writes for a given record all the fields, it is used for writing both gcm and im netcdf files.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_put_var, nf90_sync
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_netcdf_file_type)                                                    , INTENT(IN)           :: nc
    INTEGER                                                                           , INTENT(IN)           :: record_number   ! time record number
    REAL(dp), DIMENSION(nc%number_of_fields,nc%LEN_DIM(1),nc%LEN_DIM(2),nc%LEN_DIM(3)), INTENT(IN)           :: fields
    REAL(dp)                                                                          , INTENT(IN), OPTIONAL :: time

    ! Local variables:
    INTEGER                                                                                                  :: field_counter   ! The counter in the loop over the field numbers
    REAL(dp)                                                                                                 :: written_time

    IF(PRESENT(time)) THEN
     written_time = time
    ELSE
     written_time = record_number
    END IF

    ! Put the time variables:
    IF(nc%include_time_dimension) THEN
     field_counter = 0
     CALL handle_error(nf90_put_var( nc%ncid, nc%id(field_counter), written_time               , start=(/         record_number/)), &
      '. [1] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
    END IF

    ! Put the other field variables:
    DO field_counter = 1, nc%number_of_fields

     IF(nc%include_time_dimension) THEN
      IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,:), start=(/1, 1, 1, record_number/)), &
        '. [2] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,1,1), start=(/1,       record_number/)), &
        '. [3] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,1), start=(/1, 1,    record_number/)), &
        '. [4] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      END IF
     ELSE
      IF(nc%spatial_dimension_of_field(field_counter) == 3) THEN
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,:), start=(/1, 1, 1/)), &
        '. [5] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE IF(nc%spatial_dimension_of_field(field_counter) == 1) THEN
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,1,1), start=(/1      /)), &
        '. [6] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      ELSE
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:,:,1), start=(/1, 1   /)), &
        '. [7] From oblimap_write_netcdf_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
      END IF
     END IF

    END DO

    ! Synchronize the disk copy of a netcdf dataset: to minimize data loss in case of an abnormal termination:
    IF(C%synchronize_netcdf_writing) CALL handle_error(nf90_sync(nc%ncid), '. [8] From oblimap_write_netcdf_fields(): while synchronizing the file: '//TRIM(nc%file_name))
  END SUBROUTINE oblimap_write_netcdf_fields



  SUBROUTINE oblimap_write_netcdf_1D_fields(nc, record_number, fields, time)
    ! This routine writes for a given record all the 1D fields.
    USE oblimap_configuration_module, ONLY: dp, C
    USE netcdf, ONLY: nf90_put_var, nf90_sync
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_netcdf_file_type)                        , INTENT(IN)           :: nc
    INTEGER                                               , INTENT(IN)           :: record_number   ! time record number
    REAL(dp), DIMENSION(nc%number_of_fields,nc%LEN_DIM(1)), INTENT(IN)           :: fields
    REAL(dp)                                              , INTENT(IN), OPTIONAL :: time

    ! Local variables:
    INTEGER                                                                      :: field_counter   ! The counter in the loop over the field numbers
    REAL(dp)                                                                     :: written_time

    IF(PRESENT(time)) THEN
     written_time = time
    ELSE
     written_time = record_number
    END IF

    ! Put the time variables:
    IF(nc%include_time_dimension) THEN
     field_counter = 0
     CALL handle_error(nf90_put_var( nc%ncid, nc%id(field_counter), written_time, start=(/               record_number/)), &
        '. [1] From oblimap_write_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
    END IF

    ! Put the other field variables:
    DO field_counter = 1, nc%number_of_fields

     IF(nc%include_time_dimension) THEN
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:), start=(/1, record_number/)), &
        '. [2] From oblimap_write_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
     ELSE
       CALL handle_error(nf90_put_var(nc%ncid, nc%id(field_counter), fields(field_counter,:), start=(/1               /)), &
        '. [3] From oblimap_write_netcdf_1D_fields(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
     END IF
    END DO

    ! Synchronize the disk copy of a netcdf dataset: to minimize data loss in case of an abnormal termination:
    IF(C%synchronize_netcdf_writing) CALL handle_error(nf90_sync(nc%ncid), '. [4] From oblimap_write_netcdf_1D_fields(): while synchronizing the file: '//TRIM(nc%file_name))
  END SUBROUTINE oblimap_write_netcdf_1D_fields



  SUBROUTINE handle_error(stat, message, message_counter)
    USE oblimap_configuration_module, ONLY: C
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message
    INTEGER,          OPTIONAL, INTENT(IN) :: message_counter

    IF(stat /= nf90_noerr) THEN
     IF(PRESENT(message) .AND. PRESENT(message_counter)) THEN
      WRITE(UNIT=*,FMT='(4A, I3/, A/)') C%ERROR, ' netcdf failed because: ', TRIM(nf90_strerror(stat)), message, message_counter, '        The used config file is: '//TRIM(C%config_filename)
     ELSE IF(PRESENT(message)) THEN
      WRITE(UNIT=*,FMT='(4A/    , A/)') C%ERROR, ' netcdf failed because: ', TRIM(nf90_strerror(stat)), message                 , '        The used config file is: '//TRIM(C%config_filename)
     ELSE
      WRITE(UNIT=*,FMT='(3A/    , A/)') C%ERROR, ' netcdf failed because: ', TRIM(nf90_strerror(stat))                          , '        The used config file is: '//TRIM(C%config_filename)
     END IF
     STOP
    END IF
  END SUBROUTINE handle_error



  SUBROUTINE check_file_existence(file_name)
    USE oblimap_configuration_module, ONLY: C
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*), INTENT(IN) :: file_name

    ! Local variables:
    LOGICAL                      :: file_exists

    ! Continue or not if file exists:
    INQUIRE(EXIST = file_exists, FILE = file_name)
    IF(C%protect_file_overwriting .AND. file_exists) THEN
     WRITE(UNIT=*,FMT='(4A)') C%ERROR, ' The file "', TRIM(file_name), '" already exists! Remove it and run again!'
     STOP
    END IF
  END SUBROUTINE check_file_existence



  SUBROUTINE oblimap_close_netcdf_file(nc)
    USE netcdf, ONLY: nf90_close
    IMPLICIT NONE

    ! In/Output variables:
    TYPE(oblimap_netcdf_file_type), INTENT(INOUT) :: nc

    ! Close netcdf file:
    CALL handle_error(nf90_close(nc%ncid), '. From oblimap_close_netcdf_file(): it concerns the file '//TRIM(nc%file_name))

    DEALLOCATE(nc%ignore_reading_pre_mapped_fields)
    DEALLOCATE(nc%spatial_dimension_of_field      )
    DEALLOCATE(nc%field_name                      )
    DEALLOCATE(nc%field_unit                      )
    DEALLOCATE(nc%field_longname                  )
    DEALLOCATE(nc%id                              )
   !DEALLOCATE(nc%type                            )
   !DEALLOCATE(nc%case                            )
    DEALLOCATE(nc%LEN_DIM                         )
   !DEALLOCATE(nc%grid_size                       )
    DEALLOCATE(nc%vertical_axis                   )
  END SUBROUTINE oblimap_close_netcdf_file



  SUBROUTINE reduce_dummy_dimensions(file_name, number_of_fields, field_name, LEN_DIM_1, LEN_DIM_2)
    ! With this routine dummy dimensions which have size 1 are omitted by a call to nco command ncwa which
    ! does an dimension average over this dummy dimension. The ncwa omits a dimension with size 1.
    USE oblimap_configuration_module, ONLY: MND
    USE netcdf, ONLY: nf90_inq_dimid, nf90_inquire_dimension
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*)                                   , INTENT(IN) :: file_name                  ! The netcdf file name
    INTEGER                                            , INTENT(IN) :: number_of_fields           ! Number of Fields
    CHARACTER(LEN=128), DIMENSION(MND:number_of_fields), INTENT(IN) :: field_name                 ! The netcdf field names in the netcdf file
    INTEGER                                            , INTENT(IN) :: LEN_DIM_1                  ! Number of grid points in longitude direction
    INTEGER                                            , INTENT(IN) :: LEN_DIM_2                  ! Number of grid points in latitude  direction

    ! Local variables:
    INTEGER, PARAMETER                                              :: first_dimension = MND
    INTEGER, PARAMETER                                              :: last_dimension  = -1
    TYPE(oblimap_netcdf_file_type)                                  :: nc
    INTEGER                                                         :: field_counter
    INTEGER           , DIMENSION(first_dimension:last_dimension)   :: size_of_dimensions         ! The size of the dimensions
    LOGICAL                                                         :: ncwa_does_exist
    CHARACTER(LEN=256), DIMENSION(first_dimension:last_dimension)   :: dimension_name             ! The netcdf dimension name
    INTEGER                                                         :: first_dimension_to_loop
    INTEGER                                                         :: last_dimension_to_loop

    ncwa_does_exist = .FALSE.
    INQUIRE( file='/usr/bin/ncwa', EXIST=ncwa_does_exist)
    IF(.NOT. ncwa_does_exist) INQUIRE(file='/opt/local/bin/ncwa', EXIST=ncwa_does_exist)

    IF(ncwa_does_exist) THEN
     CALL oblimap_open_netcdf_file(file_name, number_of_fields, field_name, LEN_DIM_1, LEN_DIM_2, nc = nc)

     IF(ANY(nc%spatial_dimension_of_field(:) == 3)) THEN
      first_dimension_to_loop = first_dimension
     ELSE
      first_dimension_to_loop = first_dimension + 2
     END IF

     IF(nc%include_time_dimension) THEN
      last_dimension_to_loop = last_dimension
     ELSE
      last_dimension_to_loop = last_dimension - 2
     END IF

     DO field_counter = first_dimension_to_loop, last_dimension_to_loop, 2
       CALL handle_error(nf90_inq_dimid(nc%ncid, nc%field_name(field_counter), nc%id(field_counter)), &
        '. [ 1] From reduce_dummy_dimensions(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name))
       CALL handle_error(nf90_inquire_dimension(ncid = nc%ncid, dimid = nc%id(field_counter), name = dimension_name(field_counter), len = size_of_dimensions(field_counter)), &
        '. [ 2] From reduce_dummy_dimensions(): it concerns the field "'//TRIM(nc%field_name(field_counter))//'" in the file '//TRIM(nc%file_name)//', field number ', field_counter)
     END DO

     CALL oblimap_close_netcdf_file(nc)

     DO field_counter = first_dimension_to_loop, last_dimension_to_loop, 2
      IF(size_of_dimensions(field_counter) == 1) THEN
       ! Omit the time dimension if it is only one record long:
       WRITE(UNIT=*, FMT='(/A)') ' Dimension reduction: The dimension '//TRIM(dimension_name(field_counter))//' is omitted by an external call to ncwa:'
       CALL SYSTEM('echo \ \ ncwa  -Oh --average '//TRIM(dimension_name(field_counter))//' '//TRIM(file_name)//' '//TRIM(file_name))
       CALL SYSTEM('         ncwa  -Oh --average '//TRIM(dimension_name(field_counter))//' '//TRIM(file_name)//' '//TRIM(file_name))
      END IF
     END DO

    END IF
  END SUBROUTINE reduce_dummy_dimensions

END MODULE oblimap_read_and_write_module
