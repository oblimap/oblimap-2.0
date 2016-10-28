! File name: oblimap_convert_module.f90
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

MODULE oblimap_convert_module

CONTAINS
  SUBROUTINE oblimap_convert
    ! This program reads several field variables from a rectangular Ice Model (IM) and calcultes the
    ! longitude and latitude (lon,lat)-coordinates of all the (x,y)-grid points. The IM data set
    ! combined with these (lon,lat)-coordinates can be considered as being an irrugalar distributed
    ! data set of an Global Circulation Model (GCM).
    ! This program only convert the (x,y)-coordinates to their (lon,lat)-coordinate equivalents.
    ! Actually the fields are maintained, only the coordinates belonging to each grid cell are
    ! converted by an inverse oblique stereographic or lambert equal area projection. This means
    ! that the grid cells are irrecgularly distributed according to the (lon,lat)-coordinate system
    ! because they are not interpolated at a regular (lon,lat)-grid. The (lon,lat)-coordinate fields
    ! must be 2D due to the irregular distribution of the grid cells. The irregular gridded GCM field
    ! is written to a netcdf file. This GCM file can be used to remap/reproject these data on a IM grid
    ! which is the actual purpose.
    !
    ! The GCM geographical coordinates are in longitude (lon) and latitude (lat) in degrees, while the
    ! IM rectangular coordinates are in x and y in meters. The (lon,lat) coordinates are defined in the
    ! curved spherical surface S, and the (x,y)-coordinates in the flat surface S'. For a more extended
    ! description of the projection and the interpolation method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    !
    ! In the NAMELIST or config file, several options can be changed without compiling the program.
    ! For example the center of the area of interest (the center of the projected area) can be specified
    ! by setting the:
    !   lamda_M_config = 320  (value for the case of Greenland)
    !   phi_M_config   = 72   (value for the case of Greenland)
    ! in the config file. The exact (inverse) oblique stereographic projection can be chosen with
    !   alpha_stereographic_config = 7.5  (value for the case of Greenland)
    ! The center of the IM grid will coincide with (lamda_M_config,phi_M_config), and the extensions of
    ! the IM grid are determined by the IM grid spacings C%dx and C%dy and the IM grid sizes C%NX and C%NY.
    !
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type, initialize_config_variables, oblimap_licence
    USE oblimap_read_and_write_module, ONLY: oblimap_netcdf_file_type, oblimap_open_netcdf_file, create_netcdf_for_gcm_grid, oblimap_read_netcdf_fields, oblimap_write_netcdf_fields, oblimap_close_netcdf_file, reduce_dummy_dimensions
    USE oblimap_scan_contributions_module, ONLY: projecting_the_im_xy_coordinates_to_lonlat, determining_scan_parameters
    IMPLICIT NONE

    REAL(dp), DIMENSION(                          C%NX,C%NY                            ) :: x_coordinates_of_im_grid_points   ! The x-coordinates of the IM points in S'
    REAL(dp), DIMENSION(                          C%NX,C%NY                            ) :: y_coordinates_of_im_grid_points   ! The y-coordinates of the IM points in S'
    REAL(dp), DIMENSION(                          C%NX,C%NY                            ) :: lon_gcm
    REAL(dp), DIMENSION(                          C%NX,C%NY                            ) :: lat_gcm
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NX,C%NY,C%number_of_vertical_layers) :: im_field
    REAL(dp)                                                                             :: time                              ! The time value for the considered record
    INTEGER                                                                              :: field_counter                     ! The counter in the loop over the field numbers
    INTEGER                                                                              :: record_counter
    TYPE(oblimap_netcdf_file_type)                                                       :: im_netcdf_file
    TYPE(oblimap_netcdf_file_type)                                                       :: gcm_netcdf_file
    TYPE(oblimap_scan_parameter_type)                                                    :: advised_scan_parameter

    ! Opening the IM netcdf file:
    ! Output: im_netcdf_file, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points
    CALL oblimap_open_netcdf_file(C%im_input_filename, C%number_of_mapped_fields, C%im_field_name, C%NX, C%NY, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, nc = im_netcdf_file)

    ! Projection of the IM coordinates to the GCM coordinates with the inverse oblique stereographic projection:
    ! Output: lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points
    CALL projecting_the_im_xy_coordinates_to_lonlat(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm)

    ! It is required that all angles are in the 0 - 360 degree range:
    WHERE(lon_gcm <    0._dp) lon_gcm = lon_gcm + 360._dp
    WHERE(lon_gcm >= 360._dp) lon_gcm = lon_gcm - 360._dp

    ! An option to print the corner coordinates:
    IF(C%oblimap_message_level > 0) THEN
     WRITE(UNIT=*, FMT='(2(A, F21.16) )') '  longitude(1,NY) = ', lon_gcm(1,C%NY), '   longitude(NX,NY) = ', lon_gcm(C%NX,C%NY)
     WRITE(UNIT=*, FMT='(2(A, F21.16)/)') '  longitude(1, 1) = ', lon_gcm(1,   1), '   longitude(NX, 1) = ', lon_gcm(C%NX,   1)

     WRITE(UNIT=*, FMT='(2(A, F21.16) )') '  latitude (1,NY) = ', lat_gcm(1,C%NY), '   latitude (NX,NY) = ', lat_gcm(C%NX,C%NY)
     WRITE(UNIT=*, FMT='(2(A, F21.16)/)') '  latitude (1, 1) = ', lat_gcm(1,   1), '   latitude (NX, 1) = ', lat_gcm(C%NX,   1)

     WRITE(UNIT=*, FMT='(2(A, E16.8) )') '  x(1,NY) = ', x_coordinates_of_im_grid_points(1,C%NY), '   x(NX,NY) = ', x_coordinates_of_im_grid_points(C%NX,C%NY)
     WRITE(UNIT=*, FMT='(2(A, E16.8)/)') '  x(1, 1) = ', x_coordinates_of_im_grid_points(1,   1), '   x(NX, 1) = ', x_coordinates_of_im_grid_points(C%NX,   1)

     WRITE(UNIT=*, FMT='(2(A, E16.8) )') '  y(1,NY) = ', y_coordinates_of_im_grid_points(1,C%NY), '   y(NX,NY) = ', y_coordinates_of_im_grid_points(C%NX,C%NY)
     WRITE(UNIT=*, FMT='(2(A, E16.8) )') '  y(1, 1) = ', y_coordinates_of_im_grid_points(1,   1), '   y(NX, 1) = ', y_coordinates_of_im_grid_points(C%NX,   1)
    END IF

    ! Output: advised_scan_parameter
    IF(C%oblimap_message_level > 1) CALL determining_scan_parameters('gcm-to-im', lon_gcm, lat_gcm, advised_scan_parameter)

    ! Output: -
    CALL create_netcdf_for_gcm_grid(lon_gcm, lat_gcm, im_netcdf_file, gcm_netcdf_file)

    ! From the IM netcdf file we read the IM variables which will be mapped:
    DO record_counter = 0, C%im_record_range(2) - C%im_record_range(1)
     ! Output: im_field
     CALL oblimap_read_netcdf_fields(im_netcdf_file, C%im_record_range(1) + record_counter, im_field, time)

     DO field_counter = 1, C%number_of_mapped_fields
      im_field(field_counter,:,:,:) = C%field_factor(field_counter) * im_field(field_counter,:,:,:) + C%field_shift(field_counter)
     END DO

     IF(C%im_record_range(2) - C%im_record_range(1) > 0) &
      WRITE(UNIT=*, FMT='(A, I4, A)') ' Time record ', 1 + record_counter, ' is written by OBLIMAP.'

     ! Without any interpolation the IM fields are written to a GCM output file:
     ! Output: -
     CALL oblimap_write_netcdf_fields(gcm_netcdf_file, 1 + record_counter, im_field, time)
    END DO

    ! Output: -
    CALL oblimap_close_netcdf_file(gcm_netcdf_file)
    CALL oblimap_close_netcdf_file(im_netcdf_file)

    IF(C%reduce_dummy_dimensions) CALL reduce_dummy_dimensions(C%gcm_created_filename, C%number_of_mapped_fields, C%gcm_field_name, C%NLON, C%NLAT)

    ! Finishing message:
    WRITE(UNIT=*, FMT='(/3A/2A/)') ' Finished! The file  ', TRIM(C%gcm_created_filename), '  is created. Which can be viewed by:', '  ncview ', TRIM(C%gcm_created_filename)
  END SUBROUTINE oblimap_convert

END MODULE oblimap_convert_module
