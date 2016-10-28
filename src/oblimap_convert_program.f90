! File name: oblimap_convert_program.f90
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

PROGRAM oblimap_convert_program
  USE oblimap_configuration_module, ONLY: initialize_config_variables, oblimap_licence
  USE oblimap_convert_module, ONLY: oblimap_convert
  IMPLICIT NONE

  ! Read the configuration file and initialization of the struckt C%:
  CALL initialize_config_variables()

  ! Output: -
  CALL oblimap_licence('oblimap_convert_program')

  ! Calling the oblimap_gcm_to_im_mapping :
  CALL oblimap_convert()

END PROGRAM oblimap_convert_program
