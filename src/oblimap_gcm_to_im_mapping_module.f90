! File name: oblimap_gcm_to_im_mapping_module.f90
!
! Copyright (C) 2016 Thomas Reerink.
!
! This file is distributed under the terms of the
! GNU General Public License.
!
! This file is part of OBLIMAP 2.0
!
! The scientific documentation of OBLIMAP is published at:
!  http://www.geosci-model-dev.net/3/13/2010/gmd-3-13-2010.html
!  http://www.geosci-model-dev-discuss.net/gmd-2016-124/#discussion
!
! The OBLIMAP User Guide can be found at:
!  https://github.com/oblimap/oblimap-2.0/tree/master/documentation
!
! The OBLIMAP code can be downloaded by:
!  svn checkout https://svn.science.uu.nl/repos/project.oblimap
! Or from OBLIMAP's Github:
!  https://github.com/oblimap/oblimap-2.0
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

MODULE oblimap_gcm_to_im_mapping_module

CONTAINS
  SUBROUTINE oblimap_gcm_to_im_mapping()
    ! This program reads several field variables from a Global Circulation Model (GCM) netcdf file,
    ! which are projected by an oblique stereographic projection on a rectangular Ice Model (IM)
    ! grid. The GCM grid point coordinates are projected on the IM grid, where in general they fall
    ! irregular and between the IM grid points. Therefore the field values defined on the GCM grid
    ! points have to be interpolated at each IM grid point, with help of the nearby projected GCM
    ! grid points. Two interpolation methods are available in this program both based on a distance
    ! weigthing Shepard technique. One method, the 'quadrant method', searches within each quadrant
    ! around each IM grid point the nearest projected GCM point and interpolates with help of this
    ! four points to obtain the field value at such a IM grid point. Another method, the 'radius
    ! method', searches all projected GCM points within a certain radius and interpolates with help
    ! of these points to obtain the field value at such a IM grid point. Finally the projected and
    ! interpolated IM fields are written to an IM netcdf file.
    !
    ! The GCM geographical coordinates are in longitude (lon) and latitude (lat) in degrees, while the
    ! IM rectangular coordinates are in x and y in meters. The (lon,lat) coordinates are defined in the
    ! curved spherical surface S, and the (x,y) coordinates in the flat surface S'. For a more extended
    ! description of the projection and the interpolation method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    !
    ! In the NAMELIST or config file, several options can be changed without compiling the program.
    ! For example the center of the area of interest (the center of the projected area) can be specified
    ! by setting the:
    !   lamda_M_config = 319  (value for the case of Greenland)
    !   phi_M_config   =  72  (value for the case of Greenland)
    ! in the config file. The exact (inverse) oblique stereographic projection can be chosen with
    !   alpha_stereographic_config = 7.1 (value for the case of Greenland)
    ! The center of the IM grid will coincide with (lamda_M_config,phi_M_config), and the extensions of
    ! the IM grid are determined by the IM grid spacings C%dx and C%dy and the IM grid sizes C%NX and C%NY.
    !
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    USE oblimap_read_and_write_module, ONLY: oblimap_netcdf_file_type, oblimap_open_netcdf_file, initialize_im_coordinates, create_netcdf_for_im_grid, &
          oblimap_read_netcdf_fields, oblimap_write_netcdf_fields, oblimap_close_netcdf_file, reduce_dummy_dimensions
    USE oblimap_scan_contributions_module, ONLY: scan_with_quadrant_method_gcm_to_im, scan_with_radius_method_gcm_to_im, check_for_GCM_points_at_the_point_of_projection, &
          shifting_center_im_grid, determining_scan_parameters
    USE oblimap_mapping_module, ONLY: oblimap_ddo_type, oblimap_read_sid_file, oblimap_mapping, oblimap_deallocate_ddo
    IMPLICIT NONE

    ! Local variables:
    REAL(dp), DIMENSION(                          C%NLON,C%NLAT                            ) :: lon_gcm
    REAL(dp), DIMENSION(                          C%NLON,C%NLAT                            ) :: lat_gcm
    REAL(dp), DIMENSION(                          C%NX  ,C%NY                              ) :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    REAL(dp), DIMENSION(                          C%NX  ,C%NY                              ) :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'
    LOGICAL,  DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers) :: mask_of_invalid_contributions    ! For each field and for each layer a mask represents the invalid contributions (like e.g. missing values) of the GCM grid points
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers) :: gcm_field
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NX  ,C%NY  ,C%number_of_vertical_layers) :: im_field
    REAL(dp)                                                                                 :: time                             ! The time value for the considered record
    INTEGER                                                                                  :: field_counter                    ! The counter in the loop over the field numbers
    INTEGER                                                                                  :: record_counter                   ! The counter over the time records
    INTEGER                                                                                  :: layer_counter                    ! The counter over the vertical layers
    TYPE(oblimap_netcdf_file_type)                                                           :: im_netcdf_file
    TYPE(oblimap_netcdf_file_type)                                                           :: gcm_netcdf_file
    TYPE(oblimap_ddo_type)                                                                   :: oblimap_ddo                      ! The DDO containing all the scanned contributions
    TYPE(oblimap_scan_parameter_type)                                                        :: advised_scan_parameter


    ! Opening the GCM netcdf file, and reading the longitude and latitude coordinates of the GCM grid:
    ! Output: gcm_netcdf_file, lon_gcm, lat_gcm
    CALL oblimap_open_netcdf_file(C%gcm_input_filename, C%number_of_mapped_fields, C%gcm_field_name, C%NLON, C%NLAT, lon_gcm, lat_gcm, nc = gcm_netcdf_file)

    IF(C%choice_projection_method == 'rotation_projection') THEN
     lon_gcm = lon_gcm - C%shift_x_coordinate_rotation_projection
     lat_gcm = lat_gcm - C%shift_y_coordinate_rotation_projection
    ELSE
     ! It is required that all angles are in the 0 - 360 degree range:
     WHERE(lon_gcm <    0._dp) lon_gcm = lon_gcm + 360._dp
     WHERE(lon_gcm >= 360._dp) lon_gcm = lon_gcm - 360._dp

     ! In/Output: lon_gcm, lat_gcm
     CALL check_for_GCM_points_at_the_point_of_projection(lon_gcm, lat_gcm)
    END IF

    ! Determine the x, y coordinates of the IM grid points in plane S'
    ! Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    CALL initialize_im_coordinates(C%NX, C%NY, C%dx, C%dy, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)

    ! In/Output: x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    IF(C%enable_shift_im_grid) CALL shifting_center_im_grid(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)


    ! Ahead of the mapping, the scanning has to be done. In case that the file C%sid_filename exists,
    ! the scanning phase can be omitted because this file contains the scanned projection data. This file contains for each
    ! target grid point the coordinates and the distances between the contributing points and the target point. These
    ! coordinates and distances are used for the field interpolation during the mapping phase.
    ! The points which are projected most nearby the target point contribute relative to the sheppard weighing depending on
    ! their distance. There are two methods: the 'radius method' and the 'quadrant method' which can be used to select the
    ! contributing points. This 'scan' part in fact contains the most technical and CPU consuming part of OBLIMAP, covering
    ! the projection and the selection of points for the interpolation (and keeping the relative distances of each projected
    ! point relative to the target point).
    IF(C%scanning_mode) THEN
     ! Output: advised_scan_parameter
     CALL determining_scan_parameters('gcm-to-im', lon_gcm, lat_gcm, advised_scan_parameter)

     IF(C%choice_quadrant_method) THEN
      ! Output: the C%sid_filename file is created
      CALL scan_with_quadrant_method_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
     ELSE
      ! Output: the C%sid_filename file is created
      CALL scan_with_radius_method_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
     END IF
    END IF


    ! Reading the contributions of the scanned projection data into the Dynamic Data Object (DDO):
    ! Output: oblimap_ddo
    CALL oblimap_read_sid_file(C%sid_filename, oblimap_ddo)

    ! The IM netcdf file is created, this file contains the IM fields which are mapped on the IM grid:
    ! Output: im_netcdf_file
    CALL create_netcdf_for_im_grid(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, gcm_netcdf_file, im_netcdf_file)

    DO record_counter = 0, C%gcm_record_range(2) - C%gcm_record_range(1)
     ! Output: gcm_field, time
     CALL oblimap_read_netcdf_fields(gcm_netcdf_file, C%gcm_record_range(1) + record_counter, gcm_field, time)

     ! Determine the mask_of_invalid_contributions based on the invalid values for the specified field:
     mask_of_invalid_contributions = .FALSE.
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
       IF(C%masked_fields(field_counter)) WHERE(gcm_field(C%field_which_determines_invalid_value_mask(field_counter),:,:,layer_counter) == C%invalid_input_value(field_counter)) &
        mask_of_invalid_contributions(field_counter,:,:,layer_counter) = .TRUE.
     END DO
     END DO

     ! The GCM fields are mapped (= projected + interpolated) on the IM grid. For each target grid point the coordinates and
     ! the relative distances of the nearest projected points are stored in the C%sid_filename file, these are
     ! used here to map the fields:
     ! Output: im_field
     CALL oblimap_mapping(oblimap_ddo, C%NLON, C%NLAT, C%NX, C%NY, mask_of_invalid_contributions, gcm_field, im_field)

     ! Rescaling each field by multiplication with a gcm_to_im_factor and by adding a gcm_to_im_shift (in case the units differ):
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
       WHERE(im_field(field_counter,:,:,layer_counter) /= C%invalid_input_value(field_counter))
        im_field(field_counter,:,:,layer_counter) = C%field_factor(field_counter) * im_field(field_counter,:,:,layer_counter) + C%field_shift(field_counter)
       ELSEWHERE
        im_field(field_counter,:,:,layer_counter) = C%invalid_output_value(field_counter)
       END WHERE
     END DO
     END DO

     IF(C%gcm_record_range(2) - C%gcm_record_range(1) > 0) &
      WRITE(UNIT=*, FMT='(A, I4, A)') ' Time record ', 1 + record_counter, ' is written by OBLIMAP.'

     ! Finally the mapped fields im_field are written to an output file:
     ! Output: -
     CALL oblimap_write_netcdf_fields(im_netcdf_file, 1 + record_counter, im_field, time)
    END DO

    ! Output: -
    CALL oblimap_close_netcdf_file(gcm_netcdf_file)
    CALL oblimap_close_netcdf_file(im_netcdf_file)

    IF(C%reduce_dummy_dimensions) CALL reduce_dummy_dimensions(C%im_created_filename, C%number_of_mapped_fields, C%im_field_name, C%NX, C%NY)

    IF(C%oblimap_message_level > 0) THEN
     WRITE(UNIT=*, FMT='(A)')
     DO field_counter = 1, C%number_of_mapped_fields
     DO layer_counter = 1, C%number_of_vertical_layers
      IF(C%masked_fields(field_counter)) WRITE(UNIT=*, FMT='(A, E24.16, A, A20, A, I3)') &
       ' Using masked mapping for the invalid value = ', C%invalid_input_value(field_counter), ' in the field ', TRIM(C%gcm_field_name(field_counter)), '   for layer = ', layer_counter
     END DO
     END DO
    END IF

    ! Finishing message:
    WRITE(UNIT=*, FMT='(/3A/2A/)') ' Finished! The file  ', TRIM(C%im_created_filename), '  is created. Which can be viewed by:', '  ncview ', TRIM(C%im_created_filename)

    ! Output: -
    CALL oblimap_deallocate_ddo(oblimap_ddo)

  END SUBROUTINE oblimap_gcm_to_im_mapping

END MODULE oblimap_gcm_to_im_mapping_module
