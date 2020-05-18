! File name: oblimap_configuration_module.f90
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

MODULE oblimap_configuration_module
      IMPLICIT NONE

      ! PRECISION
      ! =========
      ! The kind of real numbers used by default throughout the program.
      ! Reals should be declared as:
      !
      ! REAL(dp) :: example
      !  or
      ! REAL(KIND=dp) :: example
      !
      ! dp must be a PARAMETER
      INTEGER, PARAMETER :: dp        = KIND(1.0D0)     ! Kind of double precision numbers.

      INTEGER, PARAMETER :: N_ignore  = 2 * 100         ! The size of the array containing the ignored ranges
      INTEGER, PARAMETER :: MNF       = 100             ! The Maximum Number of Fields which can be mapped simultaneously
      INTEGER, PARAMETER :: MND       =  -7             ! The Maximum Number of Dimensions of the simultaneously mapped fields times minus one


   ! CONFIG VARIABLES:
   !==================
   ! Variables which are set eventually by the read_config_file subroutine:


      ! GRID SIZES AND GRID SPACING
      !============================
      ! GCM grid sizes:
      INTEGER                                    :: NLON_config                                      = 128                                                  ! config variable
      INTEGER                                    :: NLAT_config                                      = 64                                                   ! config variable

      ! IM grid sizes and grid spacing:
      INTEGER                                    :: NX_config                                        =   100                                                ! config variable
      INTEGER                                    :: NY_config                                        =     1                                                ! config variable
      REAL(dp)                                   :: dx_config                                        = 20000.0_dp                                           ! config variable
      REAL(dp)                                   :: dy_config                                        = 20000.0_dp                                           ! config variable

      ! OBLIMAP
      ! =======
      INTEGER                                    :: oblimap_message_level_config                     = 0                                                    ! config variable
      LOGICAL                                    :: suppress_check_on_scan_parameters_config         = .FALSE.                                              ! config variable
      REAL(dp)                                   :: oblimap_allocate_factor_config                   = 2.0_dp                                               ! config variable
      CHARACTER(LEN=256)                         :: choice_projection_method_config                  = 'oblique_stereographic_projection'                   ! config variable
      REAL(dp)                                   :: lambda_M_config                                  = 320.0_dp                                             ! config variable
      REAL(dp)                                   :: phi_M_config                                     = 72.0_dp                                              ! config variable
      REAL(dp)                                   :: alpha_stereographic_config                       = 7.5_dp                                               ! config variable
      REAL(dp)                                   :: theta_rotation_projection_config                 = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: shift_x_coordinate_rotation_projection_config    = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: shift_y_coordinate_rotation_projection_config    = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: shepard_exponent_config                          = 2.0_dp                                               ! config variable
      CHARACTER(LEN=256)                         :: gcm_input_filename_config                        = 'gcm-input-no-name.nc'                               ! config variable
      LOGICAL                                    :: enable_shift_im_grid_config                      = .FALSE.                                              ! config variable
      REAL(dp)                                   :: shift_x_coordinate_im_grid_config                = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: shift_y_coordinate_im_grid_config                = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: alternative_lambda_for_center_im_grid_config     = 0.0_dp                                               ! config variable
      REAL(dp)                                   :: alternative_phi_for_center_im_grid_config        = 0.0_dp                                               ! config variable
      INTEGER           , DIMENSION(2)           :: gcm_record_range_config                          = (/ 1, 1 /)                                           ! config variable
      INTEGER           , DIMENSION(2)           :: im_record_range_config                           = (/ 1, 1 /)                                           ! config variable
      INTEGER                                    :: number_of_vertical_layers_config                 = 0                                                    ! config variable
      INTEGER                                    :: number_of_mapped_fields_config                   = 3                                                    ! config variable
      LOGICAL           , DIMENSION(MND:MNF)     :: ignore_reading_pre_mapped_fields_config          = .FALSE.                                              ! config variable
      INTEGER           , DIMENSION(  1:MNF)     :: field_which_determines_invalid_value_mask_config = -1                                                   ! config variable
      INTEGER           , DIMENSION(  1:MNF)     :: invalid_value_mask_criterion_config              =  1                                                   ! config variable
      REAL(dp)          , DIMENSION(  1:MNF)     :: invalid_input_value_config                       = -9999._dp                                            ! config variable
      REAL(dp)          , DIMENSION(  1:MNF)     :: invalid_output_value_config                      = -9999._dp                                            ! config variable
      CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: gcm_field_name_config                                                                                   ! config variable
      CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: gcm_field_unit_config                                                                                   ! config variable
      CHARACTER(LEN=256), DIMENSION(MND:MNF)     :: gcm_field_longname_config                                                                               ! config variable
      CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: im_field_name_config                                                                                    ! config variable
      CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: im_field_unit_config                                                                                    ! config variable
      CHARACTER(LEN=256), DIMENSION(MND:MNF)     :: im_field_longname_config                                                                                ! config variable
      CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: prefabricated_im_grid_field_name_config                                                                 ! config variable
      REAL(dp)          , DIMENSION(  0:MNF)     :: field_factor_config                              = 1.0_dp                                               ! config variable
      REAL(dp)          , DIMENSION(  0:MNF)     :: field_shift_config                               = 0.0_dp                                               ! config variable
      CHARACTER(LEN=256)                         :: im_created_filename_config                       = 'im-created-no-name.nc'                              ! config variable
      LOGICAL                                    :: scanning_mode_config                             = .TRUE.                                               ! config variable
      LOGICAL                                    :: data_set_is_cyclic_in_longitude_config           = .FALSE.                                              ! config variable
      INTEGER                                    :: level_of_automatic_oblimap_scanning_config       = 3                                                    ! config variable
      LOGICAL                                    :: choice_quadrant_method_config                    = .TRUE.                                               ! config variable
      REAL(dp)                                   :: R_search_interpolation_config                    = 16000.0_dp                                           ! config variable
      INTEGER                                    :: scan_search_block_size_config                    = -3                                                   ! config variable
      INTEGER                                    :: scan_search_block_size_step_config               = 2                                                    ! config variable
      CHARACTER(LEN=256)                         :: sid_filename_config                              = 'sid-no-name.txt'                                    ! config variable 
      CHARACTER(LEN=256)                         :: backward_sid_filename_config                     = 'backward-sid-no-name.txt'                           ! config variable 
      CHARACTER(LEN=256)                         :: im_input_filename_config                         = 'im-input-no-name.nc'                                ! config variable
      LOGICAL                                    :: use_prefabricated_im_grid_coordinates_config     = .FALSE.                                              ! config variable
      CHARACTER(LEN=256)                         :: prefabricated_im_grid_filename_config            = 'prefabricated-im-grid-coordinates-no-name.nc'       ! config variable
      CHARACTER(LEN=256)                         :: gcm_created_filename_config                      = 'to_and_fro_mapped_GCM.nc'                           ! config variable
      LOGICAL                                    :: nearest_point_assignment_config                  = .FALSE.                                              ! config variable
      LOGICAL                                    :: vincenty_method_for_ellipsoid_config             = .FALSE.                                              ! config variable
      LOGICAL                                    :: reduce_dummy_dimensions_config                   = .FALSE.                                              ! config variable

      ! UNIT CONVERSION
      ! ===============
      REAL(dp)                                   :: unit_conversion_x_ax_config                      = 1.0_dp                                               ! config variable
      REAL(dp)                                   :: unit_conversion_y_ax_config                      = 1.0_dp                                               ! config variable


      ! PHYSICAL PARAMETERS
      ! ===================
      REAL(dp)                                   :: earth_radius_config                              = 6.371221E6_dp                                        ! config variable   ! Earth Radius [m], earth_radius mean IUGG (2a +b ) / 3 = 6.371009E6_dp
      REAL(dp)                                   :: ellipsoid_semi_major_axis_config                 = 6.378137E6_dp                                        ! config variable   ! The semi-major axis of the Earth Ellipsoid [m], a in Snyder (1987) at p. 160; default WGS84 value for a = 6378137.0           meter
      REAL(dp)                                   :: ellipsoid_eccentricity_config                    = 0.08181919084262149_dp                               ! config variable   ! The exentricity     of the Earth Ellipsoid [m], e in Snyder (1987) at p. 160, default WGS84 value for e = 0.08181919084262149 meter

      LOGICAL                                    :: use_double_instead_of_float_in_netcdf_config     = .FALSE.                                              ! config variable

      LOGICAL                                    :: synchronize_netcdf_writing_config                = .FALSE.                                              ! config variable
      LOGICAL                                    :: protect_file_overwriting_config                  = .TRUE.                                               ! config variable
      LOGICAL                                    :: enable_color_messaging_in_terminal_config        = .TRUE.                                               ! config variable


    ! TYPE DEFENITIONS
    !=================

      ! This TYPE contains all the information once the config file is read never will change during the run of the program
      TYPE constants_type
                                                 ! NAME OF THE CONFIG FILE
                                                 !========================
        CHARACTER(LEN=256)                         :: config_filename

        INTEGER                                    :: MND      = MND
        INTEGER                                    :: N_ignore = N_ignore

                                                 ! GRID SIZES AND GRID SPACING
                                                 !============================
        INTEGER                                    :: NLON                                ! Number of GCM grid points in the longitudinal-direction
        INTEGER                                    :: NLAT                                ! Number of GCM grid points in the latitudinal-direction
        INTEGER                                    :: NX                                  ! Number of IM grid points in the x-direction
        INTEGER                                    :: NY                                  ! Number of IM grid points in the y-direction
        REAL(dp)                                   :: dx                                  ! IM grid spacing in x-direction [meter]
        REAL(dp)                                   :: dy                                  ! IM grid spacing in y-direction [meter]

                                                 ! MATHEMATICAL CONSTANTS
                                                 !=======================
        REAL(dp)                                   :: pi
        REAL(dp)                                   :: degrees_to_radians                  ! Conversion factor between radians and degrees
        REAL(dp)                                   :: radians_to_degrees                  ! Conversion factor between degrees and radians

                                                 ! OBLIMAP
                                                 !========
        INTEGER                                    :: oblimap_message_level
        LOGICAL                                    :: suppress_check_on_scan_parameters
        REAL(dp)                                   :: oblimap_allocate_factor
        CHARACTER(LEN=256)                         :: choice_projection_method
        REAL(dp)                                   :: lambda_M
        REAL(dp)                                   :: phi_M
        LOGICAL                                    :: polar_projection
        REAL(dp)                                   :: alpha_stereographic
        REAL(dp)                                   :: theta_rotation_projection
        REAL(dp)                                   :: shift_x_coordinate_rotation_projection
        REAL(dp)                                   :: shift_y_coordinate_rotation_projection
        REAL(dp)                                   :: shepard_exponent
        CHARACTER(LEN=256)                         :: gcm_input_filename
        LOGICAL                                    :: enable_shift_im_grid
        REAL(dp)                                   :: shift_x_coordinate_im_grid
        REAL(dp)                                   :: shift_y_coordinate_im_grid
        REAL(dp)                                   :: alternative_lambda_for_center_im_grid
        REAL(dp)                                   :: alternative_phi_for_center_im_grid
        INTEGER           , DIMENSION(2)           :: gcm_record_range
        INTEGER           , DIMENSION(2)           :: im_record_range
        LOGICAL                                    :: include_vertical_dimension
        INTEGER                                    :: number_of_vertical_layers
        INTEGER                                    :: number_of_mapped_fields
        LOGICAL           , DIMENSION(MND:MNF)     :: ignore_reading_pre_mapped_fields           ! This array contains for each field the choice if the pre-mapped field should be read or not
        LOGICAL           , DIMENSION(  1:MNF)     :: masked_fields                              ! This array contains for each field the choice if it should be treated masked or not
        INTEGER           , DIMENSION(  1:MNF)     :: field_which_determines_invalid_value_mask  ! This array contains for each field the number of the field which determines the invalid value mask. Default this will equal the field number of the considered field itself
        INTEGER           , DIMENSION(  1:MNF)     :: invalid_value_mask_criterion               ! This array contains for each field the masking criterion. If 1: the destination point gets an invalid value in case the nearest projected departure point has an invalid value (default). If 2: as long there is any valid contribution available they will be used, if no valid contribution is detected then the point gets an invalid value
        REAL(dp)          , DIMENSION(  1:MNF)     :: invalid_input_value                        ! This array contains for each field the invalid value where the invalid value mask is based on for that field
        REAL(dp)          , DIMENSION(  1:MNF)     :: invalid_output_value                       ! This array contains for each field the invalid value as written to the output file. Default this equals the input invalid value
        CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: gcm_field_name                             ! This array contains for each gcm field the name of the field
        CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: gcm_field_unit                             ! This array contains for each gcm field the unit of the field
        CHARACTER(LEN=256), DIMENSION(MND:MNF)     :: gcm_field_longname                         ! This array contains for each gcm field the longname of the field
        CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: im_field_name                              ! This array contains for each im  field the name of the field
        CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: im_field_unit                              ! This array contains for each im  field the unit of the field
        CHARACTER(LEN=256), DIMENSION(MND:MNF)     :: im_field_longname                          ! This array contains for each im  field the longname of the field
        CHARACTER(LEN=128), DIMENSION(MND:MNF)     :: prefabricated_im_grid_field_name           ! This array contains for each im  field the name of the field, actually the x and y coordinate names are used only by reading the prefabricated im grid coordinates
        REAL(dp)          , DIMENSION(  0:MNF)     :: field_factor                               ! This array contains for each field the conversion factor
        REAL(dp)          , DIMENSION(  0:MNF)     :: field_shift                                ! This array contains for each field the conversion shift
        CHARACTER(LEN=256)                         :: im_created_filename
        LOGICAL                                    :: scanning_mode
        LOGICAL                                    :: data_set_is_cyclic_in_longitude
        INTEGER                                    :: level_of_automatic_oblimap_scanning
        INTEGER                                    :: scan_search_block_size
        INTEGER                                    :: scan_search_block_size_step
        LOGICAL                                    :: choice_quadrant_method
        REAL(dp)                                   :: R_search_interpolation
        CHARACTER(LEN=256)                         :: sid_filename
        CHARACTER(LEN=256)                         :: backward_sid_filename
        CHARACTER(LEN=256)                         :: im_input_filename
        LOGICAL                                    :: use_prefabricated_im_grid_coordinates
        CHARACTER(LEN=256)                         :: prefabricated_im_grid_filename
        CHARACTER(LEN=256)                         :: gcm_created_filename
        LOGICAL                                    :: nearest_point_assignment
        LOGICAL                                    :: vincenty_method_for_ellipsoid
        LOGICAL                                    :: reduce_dummy_dimensions
        LOGICAL                                    :: full_scanning_mode
        INTEGER                                    :: unit_scanning_file_content
        CHARACTER(LEN=256)                         :: filename_sid_content
        REAL(dp)                                   :: large_distance
        REAL(dp)                                   :: fls_latitude_border   = 85._dp  ! full longitude scan for high latitude grid rows near the pole
        INTEGER                                    :: fls_grid_range        =  2
        INTEGER                                    :: fls_limited_lat_range = 21


                                                 ! OBLIMAP ELLIPSOID
                                                 !==================
        REAL(dp)                                   :: a
        REAL(dp)                                   :: e
        REAL(dp)                                   :: ellipsoid_flattening
        REAL(dp)                                   :: ellipsoid_semi_minor_axis
        REAL(dp)                                   :: am
        REAL(dp)                                   :: akm
        REAL(dp)                                   :: chi_M

        REAL(dp)                                   :: q_M
        REAL(dp)                                   :: q_polar
        REAL(dp)                                   :: beta_M
        REAL(dp)                                   :: R_q_polar
        REAL(dp)                                   :: D

                                                 ! UNIT CONVERSION
                                                 ! ===============
        REAL(dp)                                   :: unit_conversion_x_ax
        REAL(dp)                                   :: unit_conversion_y_ax

        LOGICAL                                    :: use_double_instead_of_float_in_netcdf

        LOGICAL                                    :: synchronize_netcdf_writing
        LOGICAL                                    :: protect_file_overwriting
        LOGICAL                                    :: enable_color_messaging_in_terminal

                                                 ! PHYSICAL PARAMETERS
                                                 ! ===================
        REAL(dp)                                   :: earth_radius

        CHARACTER(LEN=16)                          :: ERROR                               ! The allocation is exactly, so it is possible to omit a TRIM on this string
        CHARACTER(LEN=18)                          :: WARNING                             ! The allocation is exactly, so it is possible to omit a TRIM on this string
        CHARACTER(LEN=24)                          :: OBLIMAP_ERROR                       ! The allocation is exactly, so it is possible to omit a TRIM on this string
        CHARACTER(LEN=26)                          :: OBLIMAP_WARNING                     ! The allocation is exactly, so it is possible to omit a TRIM on this string
        CHARACTER(LEN=25)                          :: OBLIMAP_ADVICE                      ! The allocation is exactly, so it is possible to omit a TRIM on this string

        INTEGER                                    :: processor_id_process_dependent
        INTEGER                                    :: number_of_processors
        INTEGER                                    :: max_nr_of_lines_per_partition_block ! The maximum numberr of lines per partition block, in a MPI parallel approach
        INTEGER                                    :: psi_process_dependent               ! Partition starting index, in a MPI parallel approach
      END TYPE constants_type

      ! C is the 'struct' containing all the Constants from the config file and/or the defaults
      TYPE(constants_type), SAVE :: C


      ! This TYPE contains variables which are related to the parallel OBLIMAP implementation using MPI.
      TYPE parallel_type
        INTEGER :: processor_id_process_dependent
        INTEGER :: number_of_processors
        INTEGER :: max_nr_of_lines_per_partition_block ! The maximum numberr of lines per partition block
        INTEGER :: psi_process_dependent               ! Partition starting index
      END TYPE parallel_type

      ! PAR is the 'struct' containing the parallel OBLIMAP implementation using MPI.
      TYPE(parallel_type), SAVE :: PAR


      TYPE oblimap_scan_parameter_type
        ! This struct contains some crucial scan parameters.
        LOGICAL  :: data_set_is_cyclic_in_longitude ! This should be TRUE for GCM to IM mapping if the gcm data set is cyclic in longitude, i.e. the gcm grid covers the entire 0-360 degrees longitude range
        INTEGER  :: search_block_size               ! The optimal search_block_size which is used in the fast scanning mode to find the nearest projected points of a point close to the previous handled point
        REAL(dp) :: alpha_stereographic             ! The optimal alpha based on the entire grid area, alpha determines the standard paralel in the stereographic projection
        LOGICAL  :: choice_quadrant_method          ! The best interpolation method is selected: a choice between the "quadrant method" (TRUE) or the  "radius method" (FALSE)
        REAL(dp) :: R_search_interpolation          ! The optimal size of the search radius in case the "radius onterpolation method" is used
      END TYPE oblimap_scan_parameter_type



CONTAINS
  SUBROUTINE default_pre_initialization_of_constants()
    ! This routine sets defaults to the gcm_field_name_config and im_field_name_config arrays.
    ! These defaults are set before reading the config file, and will be overwritten if they
    ! appear in the config file
    IMPLICIT NONE

    gcm_field_name_config(MND:0)   = 'gcm-field-dimension-no-name'
    gcm_field_name_config(  1:MNF) = 'gcm-field-no-name'
    gcm_field_name_config(-7)      = 'NVL'
    gcm_field_name_config(-6)      = 'vertical-coordinate'
    gcm_field_name_config(-5)      = 'NLAT'
    gcm_field_name_config(-4)      = 'latitude'
    gcm_field_name_config(-3)      = 'NLON'
    gcm_field_name_config(-2)      = 'longitude'
    gcm_field_name_config(-1)      = 'NTIME'
    gcm_field_name_config( 0)      = 'time'
    gcm_field_name_config( 1)      = 'gcm-field-1-no-name'
    gcm_field_name_config( 2)      = 'gcm-field-2-no-name'
    gcm_field_name_config( 3)      = 'gcm-field-3-no-name'
    gcm_field_name_config( 4)      = 'gcm-field-4-no-name'
    gcm_field_name_config( 5)      = 'gcm-field-5-no-name'
    gcm_field_name_config( 6)      = 'gcm-field-6-no-name'
    gcm_field_name_config( 7)      = 'gcm-field-7-no-name'
    gcm_field_name_config( 8)      = 'gcm-field-8-no-name'
    gcm_field_name_config( 9)      = 'gcm-field-9-no-name'

    gcm_field_unit_config(:) = 'unit: ?'
    gcm_field_longname_config(:) = 'longname: ?'


    im_field_name_config(MND:0)   = 'im-field-dimension-no-name'
    im_field_name_config(  1:MNF) = 'im-field-no-name'
    im_field_name_config(-7)      = 'NVL'
    im_field_name_config(-6)      = 'vertical-coordinate'
    im_field_name_config(-5)      = 'NY'
    im_field_name_config(-4)      = 'y'
    im_field_name_config(-3)      = 'NX'
    im_field_name_config(-2)      = 'x'
    im_field_name_config(-1)      = 'NTIME'
    im_field_name_config( 0)      = 'time'
    im_field_name_config( 1)      = 'im-field-1-no-name'
    im_field_name_config( 2)      = 'im-field-2-no-name'
    im_field_name_config( 3)      = 'im-field-3-no-name'
    im_field_name_config( 4)      = 'im-field-4-no-name'
    im_field_name_config( 5)      = 'im-field-5-no-name'
    im_field_name_config( 6)      = 'im-field-6-no-name'
    im_field_name_config( 7)      = 'im-field-7-no-name'
    im_field_name_config( 8)      = 'im-field-8-no-name'
    im_field_name_config( 9)      = 'im-field-9-no-name'

    im_field_unit_config(:) = 'unit: ?'
    im_field_longname_config(:) = 'longname: ?'


    prefabricated_im_grid_field_name_config(MND:0)   = 'prefabricated-im-grid-field-dimension-no-name'
    prefabricated_im_grid_field_name_config(  1:MNF) = 'prefabricated-im-grid-field-no-name'
    prefabricated_im_grid_field_name_config(-7)      = 'NVL'
    prefabricated_im_grid_field_name_config(-6)      = 'vertical-coordinate'
    prefabricated_im_grid_field_name_config(-5)      = 'NY'
    prefabricated_im_grid_field_name_config(-4)      = 'y'
    prefabricated_im_grid_field_name_config(-3)      = 'NX'
    prefabricated_im_grid_field_name_config(-2)      = 'x'
    prefabricated_im_grid_field_name_config(-1)      = 'NTIME'
    prefabricated_im_grid_field_name_config( 0)      = 'time'


    C%enable_color_messaging_in_terminal = .TRUE.
    C%ERROR                              = coloring(' ERROR:', 'red')
    C%WARNING                            = coloring(' WARNING:', 'red')
    C%OBLIMAP_ERROR                      = coloring(' OBLIMAP ERROR:', 'red')
    C%OBLIMAP_WARNING                    = coloring(' OBLIMAP WARNING:', 'red')
    C%OBLIMAP_ADVICE                     = coloring(' OBLIMAP ADVICE:', 'yellow')
  END SUBROUTINE default_pre_initialization_of_constants



  SUBROUTINE read_config_file(config_file_number, config_filename)
    ! This subroutine reads the config variables from a configuration file. These config variables
    ! are defined in this module. The name of the configuration file should be specified on the
    ! command line. If no name is specified on the command line, then the default values are used.
    IMPLICIT NONE

    ! Input variables:
    INTEGER           , INTENT(IN) :: config_file_number
    CHARACTER(LEN=256), INTENT(IN) :: config_filename

    ! Local variables:
    INTEGER, PARAMETER             :: config_unit = 28              ! Unit number which is used for the configuration file.
    INTEGER                        :: ios
    CHARACTER(LEN=256)             :: filename_config_variable_list

    ! List of items in the configuration file:
    NAMELIST /CONFIG/NLON_config                                               , &
                     NLAT_config                                               , &
                     NX_config                                                 , &
                     NY_config                                                 , &
                     dx_config                                                 , &
                     dy_config                                                 , &
                     oblimap_message_level_config                              , &
                     suppress_check_on_scan_parameters_config                  , &
                     oblimap_allocate_factor_config                            , &
                     choice_projection_method_config                           , &
                     lambda_M_config                                           , &
                     phi_M_config                                              , &
                     alpha_stereographic_config                                , &
                     theta_rotation_projection_config                          , &
                     shift_x_coordinate_rotation_projection_config             , &
                     shift_y_coordinate_rotation_projection_config             , &
                     shepard_exponent_config                                   , &
                     gcm_input_filename_config                                 , &
                     enable_shift_im_grid_config                               , &
                     shift_x_coordinate_im_grid_config                         , &
                     shift_y_coordinate_im_grid_config                         , &
                     alternative_lambda_for_center_im_grid_config              , &
                     alternative_phi_for_center_im_grid_config                 , &
                     gcm_record_range_config                                   , &
                     im_record_range_config                                    , &
                     number_of_vertical_layers_config                          , &
                     number_of_mapped_fields_config                            , &
                     ignore_reading_pre_mapped_fields_config                   , &
                     field_which_determines_invalid_value_mask_config          , &
                     invalid_value_mask_criterion_config                       , &
                     invalid_input_value_config                                , &
                     invalid_output_value_config                               , &
                     gcm_field_name_config                                     , &
                     gcm_field_unit_config                                     , &
                     gcm_field_longname_config                                 , &
                     im_field_name_config                                      , &
                     im_field_unit_config                                      , &
                     im_field_longname_config                                  , &
                     prefabricated_im_grid_field_name_config                   , &
                     field_factor_config                                       , &
                     field_shift_config                                        , &
                     im_created_filename_config                                , &
                     scanning_mode_config                                      , &
                     data_set_is_cyclic_in_longitude_config                    , &
                     level_of_automatic_oblimap_scanning_config                , &
                     scan_search_block_size_config                             , &
                     scan_search_block_size_step_config                        , &
                     choice_quadrant_method_config                             , &
                     R_search_interpolation_config                             , &
                     sid_filename_config                                       , &
                     backward_sid_filename_config                              , &
                     im_input_filename_config                                  , &
                     use_prefabricated_im_grid_coordinates_config              , &
                     prefabricated_im_grid_filename_config                     , &
                     gcm_created_filename_config                               , &
                     nearest_point_assignment_config                           , &
                     vincenty_method_for_ellipsoid_config                      , &
                     reduce_dummy_dimensions_config                            , &
                     unit_conversion_x_ax_config                               , &
                     unit_conversion_y_ax_config                               , &
                     earth_radius_config                                       , &
                     ellipsoid_semi_major_axis_config                          , &
                     ellipsoid_eccentricity_config                             , &
                     use_double_instead_of_float_in_netcdf_config              , &
                     synchronize_netcdf_writing_config                         , &
                     protect_file_overwriting_config                           , &
                     enable_color_messaging_in_terminal_config

    ! Open the configuration file and read it:
    OPEN(UNIT=config_unit, FILE=TRIM(config_filename), STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF(ios /= 0) THEN
     WRITE(UNIT=*, FMT='(/3A/)') C%ERROR, ' Could not open the configuration file: ', TRIM(config_filename)
     IF(config_file_number > 1) WRITE(UNIT=*, FMT='(2A/)') C%WARNING, ' Be aware that you are trying to read more than one configuration file.'
     STOP
    END IF

    ! In the following statement the entire configuration file is read, using the namelist (NML=CONFIG)
    READ(UNIT=config_unit, NML=CONFIG, IOSTAT=ios)
    CLOSE(UNIT=config_unit)

    ! Compose the filename of the file which contains the log of the config variable list:
    IF(config_file_number < 10) THEN
     WRITE(filename_config_variable_list, FMT='(A, I1, A)') 'log-config-variable-list-0', config_file_number, '.txt'
    ELSE IF(config_file_number < 100) THEN
     WRITE(filename_config_variable_list, FMT='(A, I2, A)') 'log-config-variable-list-' , config_file_number, '.txt'
    ELSE
     WRITE(UNIT=*, FMT='(/2A/)') C%ERROR, ' The number of config files is currently limited to 100, this can be extended easily in the code.'
     STOP
    END IF

    ! Writing the log of all (initialized) config variables:
    OPEN( UNIT=31082015, FILE=filename_config_variable_list)
    WRITE(UNIT=31082015, NML=CONFIG)
    CLOSE(UNIT=31082015)

    IF(ios /= 0) THEN
     WRITE(UNIT=*, FMT='(/3A)') C%ERROR, ' while reading configuration file: ', TRIM(config_filename)
     CALL checking_the_config_variable_names(config_filename, filename_config_variable_list)
     WRITE(UNIT=*, FMT='(A)') ''
     STOP
    END IF
  END SUBROUTINE read_config_file



  SUBROUTINE checking_the_config_variable_names(config_filename, namelist_filename)
    ! This routine reads one by one the config variable names from the config file and compares
    ! them with the CONFIG NAMELIST. In case a config variable name is not present in the
    ! NAMELIST list this config variable name will be messaged and the program will be stopped.
    ! A few other syntax errors are detected and messaged as well. However, not all errors are
    ! fully explained, like e.g. a wrong array index of a config variable or syntax errors at
    ! the right side of the '=' sign.
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*), INTENT(IN) :: config_filename      ! The name of the config file which is given as an argument to the program
    CHARACTER(LEN=*), INTENT(IN) :: namelist_filename    ! The name of the file which contains the printed NAMELIST variables

    ! Local variables:
    CHARACTER(256) :: text_per_line
    INTEGER        :: line_counter
    CHARACTER(256) :: config_variable_name               ! The name of the config variable in the config file which is given as an argument to the program
    CHARACTER(256) :: namelist_variable_name             ! The name of the namelist varibale in the file which contains the printed NAMELIST variables
    INTEGER        :: status
    LOGICAL        :: config_variable_has_been_found
    INTEGER        :: string_index_1
    INTEGER        :: string_index_2
    INTEGER        :: string_index_3

    ! Opening the config file:
    OPEN(UNIT=1188, FILE=config_filename)

    line_counter = 1
    DO
     READ(UNIT=1188, FMT='(A)', IOSTAT=status) text_per_line
     IF (status < 0) exit
     ! Copying the first string of each text line (thus the name of namelist variable), the left part of the string is selected until a '(', 'space' or an 'equal sign' is encountered:
     config_variable_name = capitalize_string(text_per_line(1:SCAN(text_per_line, '( =')-1))

     IF(config_variable_name == '/') THEN
      ! No message for the part of the config file which is not considered
      exit
     ELSE IF(config_variable_name == '' .OR. config_variable_name == '&CONFIG' .OR. config_variable_name(1:1) == '!') THEN
      ! Deselect all lines without config variables.
      ! WRITE(UNIT=*, FMT='(2A)') ' deselected line: ', TRIM(config_variable_name)
     ELSE
        ! Scanning until a ')', 'space' or an 'equal sign' is encountered:
        string_index_1 = SCAN(text_per_line(2:), ') =') + 1

        ! Messaging about other invalid config syntax:
        string_index_2 = SCAN(text_per_line(string_index_1+1:), '=') + string_index_1
        string_index_3 = LEN(TRIM(text_per_line(string_index_1+1:string_index_2-1))) + string_index_1
        IF((text_per_line(string_index_3:string_index_3) /= ' ' .AND. string_index_1 /= string_index_3) .OR. &
           (SCAN(text_per_line(string_index_2 + 1:), '=') > 0) .OR. (SCAN(text_per_line(string_index_1 + 1:), '=') == 0)) THEN
         WRITE(UNIT=*, FMT='(A, I5, 4A)') ' The invalid config syntax at line ', line_counter, ' is: "', TRIM(coloring(text_per_line(string_index_1+1:string_index_3), 'red')), '", the whole line is cited below: '
         WRITE(UNIT=*, FMT='(2A)') ' ', TRIM(coloring(TRIM(text_per_line), 'red'))
        END IF

        ! Opening the config file:
        OPEN(UNIT=1088, FILE=namelist_filename)

        config_variable_has_been_found = .FALSE.
        DO
         READ(UNIT=1088, FMT='(A)', IOSTAT=status) text_per_line
         IF (status < 0) exit
         ! Copy left part of the string (thus the name of namelist variable), the left part of the string is selected until a 'space' or an 'equal sign' is encountered:
         namelist_variable_name = text_per_line(2:SCAN(text_per_line(2:), ' ='))
         IF(namelist_variable_name == '' .OR. namelist_variable_name == 'CONFIG' .OR. namelist_variable_name == '/' .OR. namelist_variable_name(1:1) == '!') THEN
          ! Deselect all lines without namelist variables.
          ! WRITE(UNIT=*, FMT='(2A)') ' deselected line: ', TRIM(namelist_variable_name)
         ELSE IF(namelist_variable_name == config_variable_name) THEN
          config_variable_has_been_found = .TRUE.
          exit
         ELSE
         END IF
        END DO
        IF(.NOT. config_variable_has_been_found) THEN
         WRITE(UNIT=*, FMT='(3A, I5, A)') ' The variable:  ', TRIM(coloring(TRIM(config_variable_name), 'red')), '  at line ', line_counter, ' is not a valid config variable.'
        END IF

        ! Closing the config file:
        CLOSE(UNIT=1088)

     END IF
     line_counter = line_counter + 1
    END DO

    ! Closing the config file:
    CLOSE(UNIT=1188)
  END SUBROUTINE checking_the_config_variable_names



  PURE FUNCTION capitalize_string(string_with_lower_cases) RESULT (string_with_upper_cases)
    ! Changes a string which contains lower case letters to a string with upper case letters
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*)                       , INTENT(IN) :: string_with_lower_cases

    ! Result variables:
    CHARACTER(LEN(string_with_lower_cases))             :: string_with_upper_cases

    ! Local variables:
    INTEGER                                             :: i
    INTEGER                                             :: index_cap

    CHARACTER(26), PARAMETER                            :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER                            :: low = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    string_with_upper_cases = string_with_lower_cases
    DO i = 1, LEN_TRIM(string_with_lower_cases)
     index_cap = INDEX(low, string_with_lower_cases(i:i))
    IF(index_cap > 0) string_with_upper_cases(i:i) = cap(index_cap:index_cap)
    END DO
  END FUNCTION capitalize_string



  FUNCTION coloring(string, color_choice) RESULT(color_string)
    ! This function will create a red color on terminal output for the string argument.
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*), INTENT(IN)           :: string
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: color_choice

    ! Result variables:
    CHARACTER(LEN=:), ALLOCATABLE          :: color_string

    ALLOCATE(CHARACTER(LEN(string)+9)      :: color_string)       ! The +9 is just enough to store the color characters

    IF(C%enable_color_messaging_in_terminal) THEN
     ! The 91m gives red, 0m sets the default back [Available colors: 90:gray, 91:red, 92:green, 93:yellow, 94:blue, 95:pink, 96:light blue]
     IF(PRESENT(color_choice)) THEN
      SELECT CASE(color_choice)
      CASE('default')
       DEALLOCATE(color_string)
       ALLOCATE(CHARACTER(LEN(string)) :: color_string)
       color_string = string
      CASE('gray')
       color_string = achar(27)//'[90m'//string//achar(27)//'[0m'
      CASE('red')
       color_string = achar(27)//'[91m'//string//achar(27)//'[0m'
      CASE('green')
       color_string = achar(27)//'[92m'//string//achar(27)//'[0m'
      CASE('yellow')
       color_string = achar(27)//'[93m'//string//achar(27)//'[0m'
      CASE('blue')
       color_string = achar(27)//'[94m'//string//achar(27)//'[0m'
      CASE('pink')
       color_string = achar(27)//'[95m'//string//achar(27)//'[0m'
      CASE('light blue')
       color_string = achar(27)//'[96m'//string//achar(27)//'[0m'
      CASE DEFAULT
       WRITE(UNIT=*, FMT='(3A)') ' The function "coloring" needs one of the following keywords: "default", "gray", "red", ', &
                                 ' "green", "yellow", "blue", "pink", "light blue" instead of: ', TRIM(color_choice)
       STOP
      END SELECT
     ELSE
      color_string = achar(27)//'[91m'//string//achar(27)//'[0m'
     END IF
    ELSE
     DEALLOCATE(color_string)
     ALLOCATE(CHARACTER(LEN(string)) :: color_string)
     color_string = string
    END IF

    RETURN

    DEALLOCATE(color_string)
  END FUNCTION coloring



  SUBROUTINE initialize_constants(config_filename)
    ! This routine puts all the constants which will never change during the run after the config file
    ! has been read, into a special constant 'struct'
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=256), INTENT(IN) :: config_filename

    C%config_filename                                = TRIM(config_filename)                                 ! The name of the config file, for error messaging

    C%enable_color_messaging_in_terminal             = enable_color_messaging_in_terminal_config
    C%ERROR                                          = coloring(' ERROR:', 'red')
    C%WARNING                                        = coloring(' WARNING:', 'red')
    C%OBLIMAP_ERROR                                  = coloring(' OBLIMAP ERROR:', 'red')
    C%OBLIMAP_WARNING                                = coloring(' OBLIMAP WARNING:', 'red')
    C%OBLIMAP_ADVICE                                 = coloring(' OBLIMAP ADVICE:', 'yellow')

    ! IM grid sizes:
    C%NX                                             = NX_config                                             ! Number of IM grid points in the x-direction
    C%NY                                             = NY_config                                             ! Number of IM grid points in the y-direction
    C%dx                                             = dx_config                                             ! IM grid spacing in x-direction [meter]
    C%dy                                             = dy_config                                             ! IM grid spacing in y-direction [meter]

    ! GCM grid sizes:
    C%NLON                                           = NLON_config                                           ! Number of GCM grid points in the longitudinal-direction
    C%NLAT                                           = NLAT_config                                           ! Number of GCM grid points in the latitudinal-direction

    C%pi                                             = 2.0_dp * ACOS(0.0_dp)                                 ! Just pi=3.14159... exactly
    C%degrees_to_radians                             = C%pi / 180.0_dp                                       ! Conversion factor between radians and degrees
    C%radians_to_degrees                             = 180.0_dp / C%pi                                       ! Conversion factor between degrees and radians

    ! OBLIMAP:
    C%oblimap_message_level                          = oblimap_message_level_config
    C%suppress_check_on_scan_parameters              = suppress_check_on_scan_parameters_config
    C%oblimap_allocate_factor                        = oblimap_allocate_factor_config
    C%choice_projection_method                       = choice_projection_method_config

    SELECT CASE(C%choice_projection_method)
    CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder','oblique_stereographic_projection_ellipsoid_snyder', &
    'oblique_lambert_equal-area_projection_snyder','oblique_lambert_equal-area_projection_ellipsoid_snyder','rotation_projection')
    CASE DEFAULT
     WRITE(UNIT=*, FMT='(/3A )') C%OBLIMAP_ERROR, ' In the config file: ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='( 2A/)') '                Invalid value for:  choice_projection_method_config = ', TRIM(C%choice_projection_method)
     STOP
    END SELECT

    ! Assign a C%lambda_M value for North and South Pole which generate the correct IM coordinate system,
    ! see Reerink et al. (2010) equation (2.3) or equation (A.53):
   !IF(phi_M_config == -90.0_dp .OR. phi_M_config == 90.0_dp) lambda_M_config = 0.0_dp                       ! This gives the default polar stereographic projection.
    C%lambda_M                                       = C%degrees_to_radians * lambda_M_config                ! Coordinate (C%lamda_M, C%phi_M) is the middle of the GCM's longitude-latitude region
    C%phi_M                                          = C%degrees_to_radians * phi_M_config                   ! of interest that will be projected to the IM, convert the degrees to radians

   !IF(ABS(phi_M_config) == 90.0_dp) THEN
   ! SELECT CASE(C%choice_projection_method)
   ! ! A work around to avoid the polar equations: in this way the oblique ones are used instead of the
   ! ! polar ones by a tiny adjustment of C%phi_M. Mind that closer to 1 the results become worse again.
   ! CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder','oblique_stereographic_projection_ellipsoid_snyder')
   !  C%phi_M = 0.999999_dp * C%phi_M
   ! END SELECT
   !END IF

    IF(ABS(phi_M_config) == 90.0_dp) THEN
     C%polar_projection = .TRUE.
    ELSE
     C%polar_projection = .FALSE.
    END IF

    ! The exact projection plane can be adjusted by specifying a certain angle alpha_stereographic
    ! This projection plane below and parallel to the tangent plane in (C%lamda_M, C%phi_M).
    C%alpha_stereographic                            = C%degrees_to_radians * alpha_stereographic_config

    ! The 2D rotation projection is defined by the angle below:
    C%theta_rotation_projection                      = C%degrees_to_radians * theta_rotation_projection_config
    C%shift_x_coordinate_rotation_projection         = shift_x_coordinate_rotation_projection_config
    C%shift_y_coordinate_rotation_projection         = shift_y_coordinate_rotation_projection_config

    C%shepard_exponent                               = shepard_exponent_config
    C%gcm_input_filename                             = gcm_input_filename_config
    C%enable_shift_im_grid                           = enable_shift_im_grid_config
    C%shift_x_coordinate_im_grid                     = shift_x_coordinate_im_grid_config
    C%shift_y_coordinate_im_grid                     = shift_y_coordinate_im_grid_config
    C%alternative_lambda_for_center_im_grid          = alternative_lambda_for_center_im_grid_config
    C%alternative_phi_for_center_im_grid             = alternative_phi_for_center_im_grid_config
    IF(choice_projection_method_config /= 'rotation_projection') THEN
     IF(C%alternative_lambda_for_center_im_grid      == 0.0_dp) C%alternative_lambda_for_center_im_grid = C%radians_to_degrees * C%lambda_M
     IF(C%alternative_phi_for_center_im_grid         == 0.0_dp) C%alternative_phi_for_center_im_grid    = C%radians_to_degrees * C%phi_M
    ELSE
     IF(C%enable_shift_im_grid) THEN
      WRITE(UNIT=*, FMT='(/2A, /2A/)') C%OBLIMAP_WARNING, '  In case "choice_projection_method_config = rotation_projection" the option "enable_shift_im_grid_config" has to be FALSE.', &
      ' OBLIMAP continues with  "enable_shift_im_grid_config = .FALSE." ignoring thus its setting in the config file: ', TRIM(C%config_filename)
      C%enable_shift_im_grid = .FALSE.
      C%shift_x_coordinate_im_grid                   = 0.0_dp
      C%shift_y_coordinate_im_grid                   = 0.0_dp
      C%alternative_lambda_for_center_im_grid        = 0.0_dp
      C%alternative_phi_for_center_im_grid           = 0.0_dp
     END IF
    END IF
    C%gcm_record_range                               = gcm_record_range_config
    C%im_record_range                                = im_record_range_config
    IF(C%gcm_record_range(2) < C%gcm_record_range(1)) THEN
     C%gcm_record_range(2)                           = C%gcm_record_range(1)
     WRITE(UNIT=*, FMT='(/2A, I3, A, I3/)') C%WARNING, ' The GCM record range is adjusted to:', C%gcm_record_range(1), '  --', C%gcm_record_range(2)
    END IF
    IF(C%im_record_range(2)  < C%im_record_range(1))  THEN
     C%im_record_range(2)                            = C%im_record_range(1)
     WRITE(UNIT=*, FMT='(/2A, I3, A, I3/)') C%WARNING, ' The IM record range is adjusted to:' , C%im_record_range(1),  '  --', C%im_record_range(2)
    END IF

    C%number_of_vertical_layers                      = number_of_vertical_layers_config                      ! Number of vertcal layers, i.e. number of grid points of the vertical coordinate
    IF(C%number_of_vertical_layers < 1) C%number_of_vertical_layers = 1
    C%number_of_mapped_fields                        = number_of_mapped_fields_config
    IF(C%number_of_mapped_fields < 1 .OR. C%number_of_mapped_fields > MNF) THEN
     WRITE(UNIT=*, FMT='(/A, I5, 2A)') '  The "number_of_mapped_fields" should be between 1 and ', MNF, ', adapt this in your: ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='(A)')          '  Or you have to higher the Maximum Number of Fields (MNF) in the configuration_module, followed by:'
     WRITE(UNIT=*, FMT='(A)')          '   make clean'
     WRITE(UNIT=*, FMT='(A/)')         '   make all'
     C%number_of_mapped_fields = 1
    END IF

    C%ignore_reading_pre_mapped_fields               = ignore_reading_pre_mapped_fields_config

    C%field_which_determines_invalid_value_mask      = field_which_determines_invalid_value_mask_config
    IF(ANY(C%field_which_determines_invalid_value_mask(1:C%number_of_mapped_fields) > C%number_of_mapped_fields)) THEN
     WRITE(UNIT=*, FMT='(/2A, I5, /2A/)') C%OBLIMAP_ERROR, ' One or more values of the "field_which_determines_invalid_value_mask_config" exceed the number_of_mapped_fields_config = ', C%number_of_mapped_fields, ' Adjust this in your config file: ', TRIM(C%config_filename)
     STOP
    END IF
    WHERE(C%field_which_determines_invalid_value_mask <= 0)
     C%masked_fields = .FALSE.   ! Default situation
    ELSEWHERE
     C%masked_fields = .TRUE.
    END WHERE

    C%invalid_value_mask_criterion                   = invalid_value_mask_criterion_config
    IF(ANY(C%invalid_value_mask_criterion(1:C%number_of_mapped_fields) < 1) .OR. ANY(C%invalid_value_mask_criterion(1:C%number_of_mapped_fields) > 2)) THEN
     WRITE(UNIT=*, FMT='(/3A/)') C%OBLIMAP_ERROR, ' The "invalid_value_mask_criterion" should be between 1 and 2. Adjust this in your config file: ', TRIM(C%config_filename)
     STOP
    END IF

    C%invalid_input_value                            = invalid_input_value_config
    C%invalid_output_value                           = invalid_output_value_config

    C%gcm_field_name                                 = gcm_field_name_config
    C%gcm_field_unit                                 = gcm_field_unit_config
    C%gcm_field_longname                             = gcm_field_longname_config
    C%im_field_name                                  = im_field_name_config
    C%im_field_unit                                  = im_field_unit_config
    C%im_field_longname                              = im_field_longname_config
    C%prefabricated_im_grid_field_name               = prefabricated_im_grid_field_name_config
    C%field_factor                                   = field_factor_config
    C%field_shift                                    = field_shift_config

    C%im_created_filename                            = im_created_filename_config
    C%scanning_mode                                  = scanning_mode_config
    C%data_set_is_cyclic_in_longitude                = data_set_is_cyclic_in_longitude_config
    C%level_of_automatic_oblimap_scanning            = level_of_automatic_oblimap_scanning_config
    C%scan_search_block_size                         = scan_search_block_size_config
    C%scan_search_block_size_step                    = scan_search_block_size_step_config
    C%choice_quadrant_method                         = choice_quadrant_method_config
    C%R_search_interpolation                         = R_search_interpolation_config
    C%sid_filename                                   = sid_filename_config
    C%backward_sid_filename                          = backward_sid_filename_config
    C%im_input_filename                              = im_input_filename_config
    C%use_prefabricated_im_grid_coordinates          = use_prefabricated_im_grid_coordinates_config
    C%prefabricated_im_grid_filename                 = prefabricated_im_grid_filename_config

    C%gcm_created_filename                           = gcm_created_filename_config
    C%nearest_point_assignment                       = nearest_point_assignment_config
    C%vincenty_method_for_ellipsoid                  = vincenty_method_for_ellipsoid_config
    C%reduce_dummy_dimensions                        = reduce_dummy_dimensions_config

    ! With scan_search_block_size_config = -1 the much slower full scanning mode is available for testing and benchmarking:
    IF(C%scan_search_block_size == -1) THEN
     C%full_scanning_mode                            = .TRUE.
    ELSE
     C%full_scanning_mode                            = .FALSE.
    END IF

    ! Unit used for the temporal file containing the fast input file content:
    C%unit_scanning_file_content                     = 14111984
    C%filename_sid_content                           = 'content-sid-file.txt'

    ! A predefined large distance (more then the earth circumference) used to initialize distances when searching the nearest projected points:
    C%large_distance                                 = 1.0E8_dp

    ! Description of the numbers which are used for the WGS84 ellipsoid:
    ! The relation between the flattening f and the eccentricity is given in Snyder (1987) p. 13:
    !  e^2 = 2f - f^2  or  f = 1 - (1 - e^2)^(0.5)
    ! The relation between the flattening f and the equatorial axis a and the polar axis b is
    ! given in Snyder (1987) p. 12 in the table description:
    !  b = a(1 - f)  or  f = 1 - b/a
    ! Where a and b repectively might have alternative names: semi-major axis a and semi-minor
    ! axis b, below given in meters.
    ! Litteral numbers taken from:                           Ellipsoid reference  Semi-major axis a  Semi-minor axis b  Inverse flattening (1/f)
    !  http://en.wikipedia.org/wiki/Earth_ellipsoid        :      WGS 1984          6378137          6356752.314245179    298.257223563
    !  https://en.wikipedia.org/wiki/World_Geodetic_System :      WGS 1984          6378137          6356752.314245       298.257223563
    ! Note that:
    !  the eccentricity C%e below is calculted with python by e = (2*f - f**2)**0.5 = 0.08181919084262149 with f = 1 / 298.257223563
    !  the semi-minor axis b      is calculted with python by b = a*(1 - f)         = 6356752.314245179   with f = 1 / 298.257223563 and a = 6378137

    C%a                                              = ellipsoid_semi_major_axis_config       ! The semi-major axis or the equatorial radius of the ellipsoid (in case of the Earth), a in Snyder (1987) at p. 160; default WGS84 value for a = 6.378137E6
    C%e                                              = ellipsoid_eccentricity_config          ! The eccentricity of the ellipsoid, e in Snyder (1987) at p. 160, WGS84, see Snyder (1987) p. 13, default WGS84 value for e = 0.08181919084262149
    C%ellipsoid_flattening                           = 1._dp - (1._dp - C%e**2)**0.5_dp       ! Flattening of the ellipsoid, f = 1 - (1 - e**2)**0.5 in Snyder (1987) at p. 13, WGS84 value for f = 0.0033528106647474805
    C%ellipsoid_semi_minor_axis                      = C%a * (1._dp - C%ellipsoid_flattening) ! The semi-minor axis or the polar radius of the ellipsoid (in case of the Earth) b in Snyder (1987) at p. 160, b = a(1-f) given in Snyder (1987) p. 12 in the table description; WGS84 value for b = 6356752.314245179 meter

    ! See equations (14-15) and (21-27) on page 160 in Snyder (1987), am corresponds with a*m_1 and akm corresponds with 2a*k0*m_1 in Snyder (1987):
    C%am                                             = C%a * (COS(C%phi_M) / SQRT(1.0_dp - (C%e * SIN(C%phi_M))**2))
    C%akm                                            = (1.0_dp + COS(C%alpha_stereographic)) * C%am
    ! See equations (3-1a) on page 160 in Snyder (1987),  chi_M corresponds with chi_1 in Snyder (1987):
    C%chi_M                                          = 2.0_dp * ATAN(SQRT(((1.0_dp + SIN(C%phi_M)) / (1.0_dp - SIN(C%phi_M))) * &
                                                         ((1.0_dp - C%e * SIN(C%phi_M)) / (1.0_dp + C%e * SIN(C%phi_M)))**(C%e))) - 0.5_dp * C%pi

    ! See equation (3-12) on page 187 in Snyder (1987):
    C%q_M                                            = (1.0_dp - C%e**2) * ((SIN(C%phi_M) / (1.0_dp - (C%e * SIN(C%phi_M))**2)) - (1.0_dp / (2.0_dp * C%e)) * LOG((1.0_dp - C%e * SIN(C%phi_M)) / (1.0_dp + C%e * SIN(C%phi_M))))
    ! See equation (3-12) on page 187 in Snyder (1987):
    C%q_polar                                        = (1.0_dp - C%e**2) * ((1.0_dp / (1.0_dp - C%e**2)) - (1.0_dp / (2.0_dp * C%e)) * LOG((1.0_dp - C%e) / (1.0_dp + C%e)))
    ! See equation (3-11) on page 187 in Snyder (1987):
    C%beta_M                                         = ASIN(C%q_M / C%q_polar)
    ! See equation (3-13) on page 187 in Snyder (1987):
    C%R_q_polar                                      = C%a * SQRT(0.5_dp * C%q_polar)
    ! See equation (24-20) on page 187 in Snyder (1987):
    C%D                                              = C%am / (C%R_q_polar * COS(C%beta_M))

    C%unit_conversion_x_ax                           = unit_conversion_x_ax_config
    C%unit_conversion_y_ax                           = unit_conversion_y_ax_config

    C%use_double_instead_of_float_in_netcdf          = use_double_instead_of_float_in_netcdf_config

    C%synchronize_netcdf_writing                     = synchronize_netcdf_writing_config                     ! If True the netcdf writing will be synchronized after each record which is an advantage in case the program is aborted or crashed because all fields up to then are written. But it can be signifcantly slower.
    C%protect_file_overwriting                       = protect_file_overwriting_config

    ! Physical parameters:
    C%earth_radius                                   = earth_radius_config                                   ! Earth Radius [m]
  END SUBROUTINE initialize_constants



  SUBROUTINE initialize_config_variables_for_one_config_file(config_file_number, config_filename)
    IMPLICIT NONE

    ! Input variables:
    INTEGER           , INTENT(IN) :: config_file_number
    CHARACTER(LEN=256), INTENT(IN) :: config_filename

    ! Pre initialization of *_field_*_config arrays like e.g.: gcm_field_name_config:
    CALL default_pre_initialization_of_constants()

    ! Read the configuration file:
    CALL read_config_file(config_file_number, config_filename)

    ! Initialization of the struckt  C% :
    CALL initialize_constants(config_filename)
  END SUBROUTINE initialize_config_variables_for_one_config_file



  SUBROUTINE initialize_config_variables
    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256) :: config_filename
    INTEGER            :: total_number_of_config_files

    total_number_of_config_files = COMMAND_ARGUMENT_COUNT()

    SELECT CASE(total_number_of_config_files)
    CASE(0)
     ! In this case no configuration file is read, instead the default defined values in the configuration module are used.
    CASE DEFAULT
     ! Get the name of the configuration file:
     CALL getarg(1, config_filename)

     ! Initialize the C% struct for each config file:
     CALL initialize_config_variables_for_one_config_file(1, config_filename)
    END SELECT
  END SUBROUTINE initialize_config_variables



  SUBROUTINE check_directory_existence(full_string)
    ! This subroutine checks whether the directory exists if a directory path is part of an filename.
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*), INTENT(IN) :: full_string

    ! Local variables:
    INTEGER                      :: index_of_last_slash
    LOGICAL                      :: file_exists

    index_of_last_slash = INDEX(trim(full_string), '/', .TRUE.)

    IF(index_of_last_slash /= 0) THEN
     ! Abort in case the directry in the path of the filename does not exist:
     INQUIRE(EXIST = file_exists, FILE = full_string(1:index_of_last_slash))
     IF(.NOT. file_exists) THEN
      WRITE(UNIT=*,FMT='(/6A/)') C%ERROR,' The directory "', TRIM(full_string(1:index_of_last_slash)), '" for the file "', TRIM(full_string), '" does not exists.'
      STOP
     END IF
    END IF
  END SUBROUTINE check_directory_existence


  SUBROUTINE oblimap_licence(program_name)
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*), INTENT(IN) :: program_name

    WRITE(UNIT=*, FMT='(/A/, A/, A/)') ' OBLIMAP 2.0', ' Copyright (C) 2016 Thomas Reerink', ' Free under GNU GPL version 3'
    WRITE(UNIT=*, FMT='( A/, 4A/)') ' The following run is executed by OBLIMAP:', '  ', './src/'//TRIM(program_name), ' ', TRIM(C%config_filename)
  END SUBROUTINE oblimap_licence



  PURE FUNCTION rounding(x, n) RESULT (x_rounded)
    ! This routine simply rounds a real(dp) value x at the n-th decimal.
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) :: x          ! x is the value which has to be rounded
    INTEGER , INTENT(IN) :: n          ! x will be rounded at the n-th decimal

    ! Result variables:
    REAL(dp)             :: x_rounded

    x_rounded = NINT(x * 10.0_dp ** n, dp) / 10.0_dp ** n
  END FUNCTION rounding

END MODULE oblimap_configuration_module
