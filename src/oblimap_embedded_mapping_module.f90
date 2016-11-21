! File name: oblimap_embedded_mapping_module.f90
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

MODULE oblimap_embedded_mapping_module

CONTAINS

  SUBROUTINE oblimap_initialize_embedded_mapping(oblimap_ddo_gcm_to_im, oblimap_ddo_im_to_gcm)
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_mapping_module, ONLY: oblimap_ddo_type, oblimap_read_sid_file
    IMPLICIT NONE

    ! Output variables:
    TYPE(oblimap_ddo_type), INTENT(OUT) :: oblimap_ddo_gcm_to_im   ! The DDO containing all the scanned contributions
    TYPE(oblimap_ddo_type), INTENT(OUT) :: oblimap_ddo_im_to_gcm   ! The DDO containing all the scanned contributions for the backward mapping

    ! Reading the contributions of the scanned projection data into the Dynamic Data Object (DDO):
    ! Output: oblimap_ddo_gcm_to_im
    CALL oblimap_read_sid_file(C%sid_filename, oblimap_ddo_gcm_to_im)

    ! Reading the contributions of the scanned projection data into the Dynamic Data Object (DDO):
    ! Output: oblimap_ddo_im_to_gcm
    CALL oblimap_read_sid_file(C%backward_sid_filename, oblimap_ddo_im_to_gcm)


  END SUBROUTINE oblimap_initialize_embedded_mapping



  SUBROUTINE oblimap_embedded_gcm_to_im_mapping(oblimap_ddo, gcm_field, im_field)
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_mapping_module, ONLY: oblimap_ddo_type, oblimap_mapping
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_ddo_type)                                                                  , INTENT(IN)  :: oblimap_ddo                   ! The DDO containing all the scanned contributions
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers), INTENT(IN)  :: gcm_field

    ! Output variables:
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NX  ,C%NY  ,C%number_of_vertical_layers), INTENT(OUT) :: im_field

    ! Local variables:
    INTEGER                                                                                               :: field_counter                 ! The counter in the loop over the field numbers
    INTEGER                                                                                               :: layer_counter                 ! The counter over the vertical layers
    LOGICAL,  DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers)              :: mask_of_invalid_contributions ! For each field and for each layer a mask represents the invalid contributions (like e.g. missing values) of the GCM grid points

    ! Determine the mask_of_invalid_contributions based on the invalid values for the specified field and layer:
    mask_of_invalid_contributions = .FALSE.
    DO field_counter = 1, C%number_of_mapped_fields
    DO layer_counter = 1, C%number_of_vertical_layers
      IF(C%masked_fields(field_counter)) WHERE(gcm_field(C%field_which_determines_invalid_value_mask(field_counter),:,:,layer_counter) == C%invalid_input_value(field_counter)) &
       mask_of_invalid_contributions(field_counter,:,:,layer_counter) = .TRUE.
    END DO
    END DO

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
  END SUBROUTINE oblimap_embedded_gcm_to_im_mapping



  SUBROUTINE oblimap_embedded_im_to_gcm_mapping(oblimap_ddo, im_field, input_gcm_field, gcm_field)
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_mapping_module, ONLY: oblimap_ddo_type, oblimap_mapping
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_ddo_type)                                                                  , INTENT(IN)  :: oblimap_ddo                   ! The DDO containing all the scanned contributions
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NX  ,C%NY  ,C%number_of_vertical_layers), INTENT(IN)  :: im_field                      ! The cluster of IM fields which will be mapped
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers), INTENT(IN)  :: input_gcm_field               ! The cluster of GCM fields in the stage before the mapping as conducted by this routine

    ! Output variables:
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers), INTENT(OUT) :: gcm_field                     ! The resulting cluster of GCM fields after the mapping as conducted by this routine

    ! Local variables:
    INTEGER                                                                                               :: field_counter                 ! The counter in the loop over the field numbers
    INTEGER                                                                                               :: layer_counter                 ! The counter over the vertical layers
    LOGICAL,  DIMENSION(C%number_of_mapped_fields,C%NX  ,C%NY  ,C%number_of_vertical_layers)              :: mask_of_invalid_contributions ! For each field and for each layer a mask represents the invalid contributions (like e.g. missing values) of the IM grid points
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NX  ,C%NY  ,C%number_of_vertical_layers)              :: im_field_converted            ! The rescaled and/or shifted input im_field
    REAL(dp), DIMENSION(C%number_of_mapped_fields,C%NLON,C%NLAT,C%number_of_vertical_layers)              :: mapped_gcm_field
    INTEGER,  DIMENSION(                          C%NLON,C%NLAT                            )              :: mapping_participation_mask    ! A mask for points which participate (mask = 1) in the mapping, so which fall within the mapped area.

    ! Determine the mask_of_invalid_contributions based on the invalid values for the specified field and layer:
    mask_of_invalid_contributions = .FALSE.
    DO field_counter = 1, C%number_of_mapped_fields
    DO layer_counter = 1, C%number_of_vertical_layers
      IF(C%masked_fields(field_counter)) WHERE(im_field(C%field_which_determines_invalid_value_mask(field_counter),:,:,layer_counter) == C%invalid_input_value(field_counter)) &
       mask_of_invalid_contributions(field_counter,:,:,layer_counter) = .TRUE.
    END DO
    END DO

    ! Rescaling each field by dividing by a gcm_to_im_factor and by subtracting a gcm_to_im_shift (in case the units differ):
    DO field_counter = 1, C%number_of_mapped_fields
    DO layer_counter = 1, C%number_of_vertical_layers
      WHERE(im_field(field_counter,:,:,layer_counter) /= C%invalid_input_value(field_counter))
       im_field_converted(field_counter,:,:,layer_counter) = im_field(field_counter,:,:,layer_counter) / C%field_factor(field_counter) - C%field_shift(field_counter)
      ELSEWHERE
       im_field_converted(field_counter,:,:,layer_counter) = C%invalid_input_value(field_counter)
      END WHERE
    END DO
    END DO

    ! Output: mapped_gcm_field
    CALL oblimap_mapping(oblimap_ddo, C%NX, C%NY, C%NLON, C%NLAT, mask_of_invalid_contributions, im_field_converted, mapped_gcm_field, mapping_participation_mask)

    ! Compose the fields which are given back to the GCM. Part of these fields (the coordinates with mapping_participation_mask == 0) are
    ! the same as the original GCM fields: input_gcm_field_*, and the rest of them equal the mapped values as in mapped_gcm_field_*.
    DO field_counter = 1, C%number_of_mapped_fields
    DO layer_counter = 1, C%number_of_vertical_layers
      gcm_field(field_counter,:,:,layer_counter) =        mapping_participation_mask(:,:)  * mapped_gcm_field(field_counter,:,:,layer_counter) &
                                                   + (1 - mapping_participation_mask(:,:)) * input_gcm_field (field_counter,:,:,layer_counter)
      WHERE(gcm_field(field_counter,:,:,layer_counter) == C%invalid_input_value(field_counter)) gcm_field(field_counter,:,:,layer_counter) = C%invalid_output_value(field_counter)
    END DO
    END DO
  END SUBROUTINE oblimap_embedded_im_to_gcm_mapping

END MODULE oblimap_embedded_mapping_module
