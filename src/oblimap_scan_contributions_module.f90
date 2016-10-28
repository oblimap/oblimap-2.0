! File name: oblimap_scan_contributions_module.f90
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

MODULE oblimap_scan_contributions_module
  USE oblimap_configuration_module, ONLY: dp
  IMPLICIT NONE

  TYPE triplet
    INTEGER  :: row_index     ! row index of nearest point within one quadrant
    INTEGER  :: column_index  ! column index of nearest point within one quadrant
    REAL(dp) :: distance      ! distance of this nearest point relative to the IM point (m,n)
  END TYPE triplet

  ! In case there are no contributions, the distance is set to a huge number: C%large_distance = 1.0E8_dp, and the indices to -999999
 !TYPE(triplet), PARAMETER :: no_contribution = triplet(-999999, -999999, 1.0E8_dp)



CONTAINS

  ! -----------------------------------------------------------------------------
  ! ROUTINES WHICH SCAN THE CONTRIBUTING POINTS FOR INTERPOLATION OF GCM TO IM
  ! -----------------------------------------------------------------------------

  SUBROUTINE scan_with_quadrant_method_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
    ! This routine selects the contributing points for each target grid point, by searching with the quadrant method. First
    ! the coordinates of the GCM grid points are projected with the oblique stereographic projection to the IM coordinates.
    ! Thereafter with these projected coordinates the distances of the projected points relative to each target grid point
    ! are calculated and used to select the nearest contributing grid points. The GCM-grid indices of the contributing points
    ! and the relative distance to 'their' target grid point are stored by writing them to the C%sid_filename
    ! file. With the indices and the distances of the contributing points the GCM fields can be mapped fast and simultaneously
    ! on to the IM grid.
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    IMPLICIT NONE

    ! Input variables:
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: x_coordinates_of_im_grid_points   ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: y_coordinates_of_im_grid_points   ! The y-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lon_gcm                           ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lat_gcm                           ! The latitude  coordinates (degrees) of the GCM grid points
    TYPE(oblimap_scan_parameter_type)        , INTENT(IN) :: advised_scan_parameter            ! The struct containing the crucial scan parameters.

    ! Local variables:
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT)             :: x_coordinates_of_gcm_grid_points  ! The x-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT)             :: y_coordinates_of_gcm_grid_points  ! The y-coordinates of the GCM points projected on S'
    LOGICAL                                               :: latitude_parallel_to_grid_numbers ! True if the latitudes increase in the same direction as the grid numbers increase, this depends on the input dataset
    INTEGER                                               :: i, j
    INTEGER                                               :: m, n
    REAL(dp)                                              :: m_message = 0._dp
    INTEGER                                               :: situation
    INTEGER                                               :: count_iteratations                ! Counting the DO WHILE iterations
    INTEGER                                               :: counter                           ! This counter counts each time a nearer point is found, if it stays zero no more points are found by extending the search block
    INTEGER                                               :: scan_search_block_size
    INTEGER                                               :: highest_scan_search_block_size = 0
    INTEGER                                               :: number_of_situations
    INTEGER                                               :: amount_of_mapped_points
    INTEGER                                               :: number_points_no_contribution = 0 ! Number of points for which no any contribution is found
    INTEGER                                               :: count_contributions
    INTEGER                                               :: quadrant                          ! The quadrant I, II, III or IV relative to an IM grid point
    INTEGER                                               :: loop
    TYPE(triplet)                                         :: projected_gcm                     ! Projected GCM point on S'
    TYPE(triplet), DIMENSION(4,C%NX  ,C%NY  )             :: contribution                      ! Nearest projected GCM point in quadrant (DIM=I,II,III or IV) in S', relative to the IM grid point
    TYPE(triplet)                                         :: no_contribution                   ! In case there are no contributions, the nearest contribution elements are set to some specific values: the distance to a huge number, and the indices are set to -1
    TYPE(triplet)                                         :: pivot_contribution                ! The selected pivot contribution, this pivot determines the position of the scan block
    TYPE(triplet), DIMENSION(  C%NX  ,C%NY  )             :: nearest_contribution              ! Keep the nearest projected GCM contribution for each IM grid point
    LOGICAL                                               :: do_full_scan                      ! Do a full scan off all projected departing grid points for this destination grid point
    INTEGER                                               :: i_start                           ! The starting i indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: i_end                             ! The ending   i indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: j_start                           ! The starting j indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: j_end                             ! The ending   j indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: i_start_previous_iteration        ! The i_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: i_end_previous_iteration          ! The i_end   of the previous iteration in the WHILE-loop
    INTEGER                                               :: j_start_previous_iteration        ! The j_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: j_end_previous_iteration          ! The j_end   of the previous iteration in the WHILE-loop
    LOGICAL,       DIMENSION(  C%NLON,C%NLAT)             :: mask
    REAL(dp)                                              :: local_gcm_grid_distance           ! The distance between two local GCM grid neighbour points

    IF(lat_gcm(1,1) < lat_gcm(1,C%NLAT)) THEN
     latitude_parallel_to_grid_numbers = .TRUE.
    ELSE
     latitude_parallel_to_grid_numbers = .FALSE.
    END IF

    ! Projection of the GCM coordinates to the IM coordinates with the oblique stereographic projection:
    ! Output: x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points
    CALL projecting_the_gcm_lonlat_coordinates_to_xy(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)

    no_contribution%distance     = C%large_distance
    no_contribution%row_index    = -1
    no_contribution%column_index = -1

    pivot_contribution           = no_contribution

    amount_of_mapped_points      = 0

    IF(C%data_set_is_cyclic_in_longitude) THEN
     number_of_situations = 3
    ELSE
     number_of_situations = 1
    END IF

    ! Opening the file to which the coordinates of the nearest projected points are written, which will be the content of the SID file:
    OPEN(UNIT=C%unit_scanning_file_content, FILE=TRIM(C%filename_sid_content))

    ! For each IM grid point the four nearest projected GCM points are determined:
    WRITE(UNIT=*,FMT='(A)') '  The progress of the OBLIMAP scanning phase is at:'
    DO m = 1, C%NX
      IF(m >= m_message) THEN
       IF(C%oblimap_message_level == 0) THEN
        WRITE(UNIT=*,FMT='(F9.1, A    )') 100._dp * REAL(m, dp) / REAL(C%NX, dp), ' %'
        m_message = m_message + 0.10_dp * C%NX
       ELSE
        WRITE(UNIT=*,FMT='(F9.1, A, I5)') 100._dp * REAL(m, dp) / REAL(C%NX, dp), ' %,  at  m = ', m
        m_message = m_message + 0.05_dp * C%NX
       END IF
      END IF
    DO n = 1, C%NY

      IF(C%full_scanning_mode .OR. (m == 1 .AND. n == 1)) THEN
       ! For the very first point always a full scan is conducted. In case the full_scanning_mode = TRUE, the full scan is conducted at any point
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      ELSE IF(n == 1) THEN
       ! Low frequent situation (At n == 1, the starting and most left column, only a contribution of the previous row can be used if it exits)
       IF(nearest_contribution(m-1,n)%distance /= C%large_distance) THEN
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(m-1,n)
       ELSE
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       END IF
      ELSE IF(nearest_contribution(m,n-1)%distance /= C%large_distance) THEN
       ! Most frequent situation (continuing at the same row): take the contribution of that neighbour point which is located at the previous column
       do_full_scan = .FALSE.
       pivot_contribution = nearest_contribution(m,n-1)
      ELSE IF(m == 1) THEN
       ! Low frequent situation (If no neighbour contribution at the same row is found, it is not possible to try the previous row, because m == 1 is the lowest and first scanned row)
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      ELSE IF(nearest_contribution(m-1,n)%distance /= C%large_distance) THEN
       ! Second frequent situation: take the contribution of that neighbour point which is located at the previous row
       do_full_scan = .FALSE.
       pivot_contribution = nearest_contribution(m-1,n)
      ELSE
       ! Low frequent situation (no contributing neighbour point is found for this point which can serve as a pivot)
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      END IF

      ! Initialize the contributions to inappropriate values:
      contribution(:,m,n) = no_contribution

      IF(do_full_scan) THEN
       scan_search_block_size = C%scan_search_block_size ! Actually a dummy value to avoid the compiler warning: scan_search_block_size may be used uninitialized
      ELSE
       IF(C%scan_search_block_size < -1) THEN
        ! Estimating the distance of two projected gcm points in the flat IM plane: to obtain a local estimate the distance between the pivot and its neighbour at the next row is calculated:
        IF(pivot_contribution%row_index == C%NLON) THEN
         local_gcm_grid_distance = distance_in_flat_surface( &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index - 1, pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index - 1, pivot_contribution%column_index) )
        ELSE
         local_gcm_grid_distance = distance_in_flat_surface( &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index + 1, pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index + 1, pivot_contribution%column_index) )
        END IF
        ! There are situations around the pole for which different grid points might almost coincide. In that case the local_gcm_grid_distance will be very small.
        ! To avoid in such case an extreme blow up of the scan_search_block_size (or even a NaN due to a zero division) we set:
        IF(local_gcm_grid_distance < 0.1_dp) local_gcm_grid_distance = C%large_distance

        ! Calculate the local scan_search_block_size based on the estimated distance of two local projected neighbour points:
        scan_search_block_size = INT(MAX(C%dx, C%dy) / local_gcm_grid_distance) + 3
       ELSE
        scan_search_block_size = C%scan_search_block_size
       END IF
      END IF

      count_iteratations = 0
      counter = -1
      ! The dynamic scan_search_block_size iteration: each iteration the scan_search_block_size is increased with the C%scan_search_block_size_step:
      iterate: DO WHILE(counter /= 0)
       counter = 0
       count_iteratations = count_iteratations + 1

       IF(.NOT. do_full_scan) THEN
        IF(C%scan_search_block_size == -3) THEN
         IF(count_iteratations > 1) scan_search_block_size = scan_search_block_size + C%scan_search_block_size_step
        END IF
        highest_scan_search_block_size = MAX(scan_search_block_size, highest_scan_search_block_size)
       END IF

       IF(do_full_scan) THEN
        ! In case there is no clue where to start, a full search is done:
        IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, 2(I5, A))') ' Full scan for (m, n) = (', m, ',', n, ')'
        j_start =      1
        j_end   = C%NLAT
       ELSE
        j_start = MAX(pivot_contribution%column_index - scan_search_block_size,      1)
        j_end   = MIN(pivot_contribution%column_index + scan_search_block_size, C%NLAT)

        IF(latitude_parallel_to_grid_numbers) THEN
         ! Due to the spreading of points close to the South Pole in longitudinal direction at the low latitude edge of the grid, a full longitude scan is done:
         IF(j_start <= C%fls_grid_range) THEN
          IF(lat_gcm(1,j_start) < - C%fls_latitude_border) THEN
           j_start = 1
           j_end   = MIN(C%fls_limited_lat_range, C%NLAT)
           do_full_scan = .TRUE.
          END IF
         END IF
         ! Due to the spreading of points close to the North Pole in longitudinal direction at the high latitude edge of the grid, a full longitude scan is done:
         IF(j_end > C%NLAT - C%fls_grid_range) THEN
          IF(lat_gcm(1,j_end)   >   C%fls_latitude_border) THEN
           j_start = MAX(C%NLAT - C%fls_limited_lat_range + 1, 1)
           j_end   = C%NLAT
           do_full_scan = .TRUE.
          END IF
         END IF
        ELSE
         ! Due to the spreading of points close to the South Pole in longitudinal direction at the low latitude edge of the grid, a full longitude scan is done:
         IF(j_end > C%NLAT - C%fls_grid_range) THEN
          IF(lat_gcm(1,j_end) < - C%fls_latitude_border) THEN
           j_start = MAX(C%NLAT - C%fls_limited_lat_range + 1, 1)
           j_end   = C%NLAT
           do_full_scan = .TRUE.
          END IF
         END IF
         ! Due to the spreading of points close to the North Pole in longitudinal direction at the high latitude edge of the grid, a full longitude scan is done:
         IF(j_start <= C%fls_grid_range) THEN
          IF(lat_gcm(1,j_start) >  C%fls_latitude_border) THEN
           j_start = 1
           j_end   = MIN(C%fls_limited_lat_range, C%NLAT)
           do_full_scan = .TRUE.
          END IF
         END IF
        END IF
       END IF

       ! For cyclic cases all three situations are passed. While passing for non-cyclic cases only situation 1.
       ! A data set which is periodical in the longitude direction is called a cyclic case.
       DO situation = 1, number_of_situations
        IF(situation == 1) THEN
         IF(do_full_scan) THEN
          ! In case there is no clue where to start, a full search is done:
          i_start =      1
          i_end   = C%NLON
         ELSE
          ! A quick search within a local block will be done:
          i_start = MAX(pivot_contribution%row_index    - scan_search_block_size,      1)
          i_end   = MIN(pivot_contribution%row_index    + scan_search_block_size, C%NLON)
         END IF
         mask(i_start:i_end,j_start:j_end) = .TRUE.
         IF(count_iteratations > 1) THEN
          ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
          i_start_previous_iteration = MAX(pivot_contribution%row_index    - (scan_search_block_size - C%scan_search_block_size_step),      1)
          i_end_previous_iteration   = MIN(pivot_contribution%row_index    + (scan_search_block_size - C%scan_search_block_size_step), C%NLON)
         END IF
        ELSE IF(situation == 2) THEN
         IF(pivot_contribution%row_index + scan_search_block_size > C%NLON) THEN
          ! Search for contributions at the west side of the grid if the east side of the grid has been reached:
          i_start = 1
          i_end   = pivot_contribution%row_index + scan_search_block_size - C%NLON
          mask(i_start:i_end,j_start:j_end) = .TRUE.
          IF(count_iteratations > 1) THEN
           ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
           i_start_previous_iteration = 1
           i_end_previous_iteration   = pivot_contribution%row_index + (scan_search_block_size - C%scan_search_block_size_step) - C%NLON
          END IF
         ELSE
          i_start = 1
          i_end   = 0  ! This will immediately stop the ij-loop (note that therefore the mask will not be evaluated in the ij-loop with an i_end = 0 which would be out of range)
         END IF
        ELSE IF(situation == 3) THEN
         IF(pivot_contribution%row_index - scan_search_block_size <      1) THEN
          ! Search for contributions at the east side of the grid if the west side of the grid has been reached:
          i_start = pivot_contribution%row_index - scan_search_block_size + C%NLON
          i_end   = C%NLON
          mask(i_start:i_end,j_start:j_end) = .TRUE.
          IF(count_iteratations > 1) THEN
           ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
           i_start_previous_iteration = pivot_contribution%row_index - (scan_search_block_size - C%scan_search_block_size_step) + C%NLON
           i_end_previous_iteration   = C%NLON
          END IF
         ELSE
          i_start = 1
          i_end   = 0  ! This will immediately stop the ij-loop (note that therefore the mask will not be evaluated in the ij-loop with an i_end = 0 which would be out of range)
         END IF
        END IF

        IF(count_iteratations > 1) THEN
         ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
         j_start_previous_iteration = MAX(pivot_contribution%column_index - (scan_search_block_size - C%scan_search_block_size_step),      1)
         j_end_previous_iteration   = MIN(pivot_contribution%column_index + (scan_search_block_size - C%scan_search_block_size_step), C%NLAT)
         mask(i_start_previous_iteration:i_end_previous_iteration,j_start_previous_iteration:j_end_previous_iteration) = .FALSE.
        END IF

        ! See equation (2.17) in Reerink et al. (2010):
        DO i = i_start, i_end
        DO j = j_start, j_end
          IF(mask(i,j)) THEN

           ! Determine the quadrant in which the projected point lies relative to the considered grid point:
           ! Output: quadrant
           CALL find_quadrant_around_IM_grid_point(x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j), &
                                                   x_coordinates_of_im_grid_points(m,n),  y_coordinates_of_im_grid_points(m,n), quadrant)

           ! Determine in the flat plane S' the distance between the projected GCM coordinates relative to the considered IM grid point:
           projected_gcm%row_index    = i
           projected_gcm%column_index = j
           projected_gcm%distance     = distance_in_flat_surface(x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j), &
                                                                 x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n))
           ! In case the projected point coincides with the grid point we put it at the very close distance of 1 centimeter, preventing devision by zero:
           IF(projected_gcm%distance == 0._dp) projected_gcm%distance = 0.01_dp

           ! Select the in S' projected GCM point with the shortest distance to the considered IM grid point in this quadrant,
           ! and keep this distance and the GCM-grid indices of this GCM point in S:
           IF(projected_gcm%distance < contribution(quadrant,m,n)%distance) THEN
            contribution(quadrant,m,n) = projected_gcm
            counter = counter + 1
           END IF

          END IF
        END DO
        END DO
        ! Leaving the do loop and the do while loop immediately if a full search was carried out:
        IF(do_full_scan) EXIT iterate
       END DO
       ! Leaving the do while loop immediately if a full search was carried out or a fixed scan_search_block_size is used:
       IF(do_full_scan .OR. (C%scan_search_block_size /= -3)) EXIT iterate
      END DO iterate

      count_contributions = 4
      DO loop = 1, 4
       IF(contribution(loop,m,n)%distance == C%large_distance) count_contributions = count_contributions - 1
      END DO

      IF(count_contributions == 0) THEN
       IF(C%oblimap_message_level > 2) WRITE(UNIT=*, FMT='(2A, 2(I5, A))') TRIM(C%OBLIMAP_WARNING), ' from scan_with_quadrant_method_gcm_to_im():  In four quadrants no single point is found for point (m, n) = (', m, ',', n, ')'

       nearest_contribution(m,n) = no_contribution

       number_points_no_contribution = number_points_no_contribution + 1
      ELSE
       ! The nearest contribution is selected:
       nearest_contribution(m,n) = contribution(MINLOC(contribution(:,m,n)%distance, 1),m,n)

       WRITE(UNIT=C%unit_scanning_file_content, FMT='(3I6)', ADVANCE='NO') m, n, count_contributions
       DO loop = 1, 4
        ! Filter the appropriate contributions (leave out the quadrants in which no contributing point is found, e.g. at the grid border):
        IF(contribution(loop,m,n)%distance /= C%large_distance) THEN
         WRITE(UNIT=C%unit_scanning_file_content, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop,m,n)%row_index, contribution(loop,m,n)%column_index, contribution(loop,m,n)%distance
        END IF
       END DO
       WRITE(UNIT=C%unit_scanning_file_content, FMT='(A)') ''
       amount_of_mapped_points = amount_of_mapped_points + 1
      END IF

    END DO
    END DO

    IF(C%scan_search_block_size == -3) highest_scan_search_block_size = highest_scan_search_block_size - 2
    IF(C%oblimap_message_level > 0) WRITE(UNIT=*,FMT='(/A, I6/)') ' The highest dynamic scan_search_block_size was: ', highest_scan_search_block_size

    ! Closing the the SID file:
    CLOSE(UNIT=C%unit_scanning_file_content)

    ! Output: -
    CALL write_sid_file(advised_scan_parameter, highest_scan_search_block_size, amount_of_mapped_points, number_points_no_contribution, maximum_contributions = 4, gcm_to_im_direction = .TRUE.)
  END SUBROUTINE scan_with_quadrant_method_gcm_to_im



  SUBROUTINE scan_with_radius_method_gcm_to_im(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
    ! This routine selects the contributing points for each target grid point, by searching with the radius method. First
    ! the coordinates of the GCM grid points are projected with the oblique stereographic projection to the IM coordinates.
    ! Thereafter with these projected coordinates the distances of the projected points relative to each target grid point
    ! are calculated and used to select the nearest contributing grid points. The GCM-grid indices of the contributing points
    ! and the relative distance to 'their' target grid point are stored by writing them to the C%sid_filename
    ! file. With the indices and the distances of the contributing points the GCM fields can be mapped fast and simultaneously
    ! on to the IM grid.
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    IMPLICIT NONE

    ! Input variables:
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: x_coordinates_of_im_grid_points   ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: y_coordinates_of_im_grid_points   ! The y-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lon_gcm                           ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lat_gcm                           ! The latitude  coordinates (degrees) of the GCM grid points
    TYPE(oblimap_scan_parameter_type)        , INTENT(IN) :: advised_scan_parameter            ! The struct containing the crucial scan parameters.

    ! Local variables:
    INTEGER                                               :: max_size
    INTEGER                                               :: status
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT)             :: x_coordinates_of_gcm_grid_points  ! The x-coordinates of the GCM points projected on S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT)             :: y_coordinates_of_gcm_grid_points  ! The y-coordinates of the GCM points projected on S'
    LOGICAL                                               :: latitude_parallel_to_grid_numbers ! True if the latitudes increase in the same direction as the grid numbers increase, this depends on the input dataset
    INTEGER                                               :: i, j
    INTEGER                                               :: m, n
    REAL(dp)                                              :: m_message = 0._dp
    INTEGER                                               :: situation
    INTEGER                                               :: count_iteratations                ! Counting the DO WHILE iterations
    INTEGER                                               :: counter                           ! This counter counts each time a nearer point is found, if it stays zero no more points are found by extending the search block
    INTEGER                                               :: scan_search_block_size
    INTEGER                                               :: highest_scan_search_block_size = 0
    INTEGER                                               :: number_of_situations
    INTEGER                                               :: amount_of_mapped_points
    INTEGER                                               :: number_points_no_contribution = 0 ! Number of points for which no any contribution is found
    INTEGER                                               :: count_contributions
    INTEGER                                               :: maximum_contributions = 0
    INTEGER                                               :: loop
    TYPE(triplet)                                         :: projected_gcm                     ! Projected GCM point on S'
    TYPE(triplet), DIMENSION(:,:,:), ALLOCATABLE          :: contribution                      ! Nearest projected GCM point in quadrant (DIM=I,II,III or IV) in S', relative to the IM grid point
    TYPE(triplet)                                         :: no_contribution                   ! In case there are no contributions, the nearest contribution elements are set to some specific values: the distance to a huge number, and the indices are set to -1
    TYPE(triplet)                                         :: pivot_contribution                ! The selected pivot contribution, this pivot determines the position of the scan block
    TYPE(triplet), DIMENSION(  C%NX  ,C%NY  )             :: nearest_contribution              ! Keep the nearest projected GCM contribution for each IM grid point
    LOGICAL                                               :: do_full_scan                      ! Do a full scan off all projected departing grid points for this destination grid point
    INTEGER                                               :: i_start                           ! The starting i indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: i_end                             ! The ending   i indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: j_start                           ! The starting j indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: j_end                             ! The ending   j indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: i_start_previous_iteration        ! The i_start of the previous iteration in the DO WHILE-loop
    INTEGER                                               :: i_end_previous_iteration          ! The i_end   of the previous iteration in the DO WHILE-loop
    INTEGER                                               :: j_start_previous_iteration        ! The j_start of the previous iteration in the DO WHILE-loop
    INTEGER                                               :: j_end_previous_iteration          ! The j_end   of the previous iteration in the DO WHILE-loop
    LOGICAL,       DIMENSION(  C%NLON,C%NLAT)             :: mask
    REAL(dp)                                              :: local_gcm_grid_distance           ! The distance between two local GCM grid neighbour points

    IF(lat_gcm(1,1) < lat_gcm(1,C%NLAT)) THEN
     latitude_parallel_to_grid_numbers = .TRUE.
    ELSE
     latitude_parallel_to_grid_numbers = .FALSE.
    END IF

    ! The devision by 1000 is to prevent the failure of CEILING with large numbers:
    max_size = CEILING(MAX(4._dp * C%pi * (C%R_search_interpolation / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)), &
                                   C%pi * (C%R_search_interpolation / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)))  * C%oblimap_allocate_factor)
    ALLOCATE(contribution(max_size,C%NX,C%NY), STAT=status)
    IF(status /= 0) THEN
     WRITE(UNIT=*, FMT='(/2A, /2(A, I8), A, F16.3, A/)') &
      C%OBLIMAP_ERROR, ' message from: scan_with_radius_method_gcm_to_im():  The allocation of the "contribution struct" exceeds your system allocation capacity.', &
      '                The combination of NLON_config = ', C%NX, ', and NLAT_config = ', C%NY, ', with R_search_interpolation_config = ', C%R_search_interpolation, ' is to large. Reduce e.g. the size of R_search_interpolation_config'
     STOP
    END IF

    ! Projection of the GCM coordinates to the IM coordinates with the oblique stereographic projection:
    ! Output: x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points
    CALL projecting_the_gcm_lonlat_coordinates_to_xy(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)

    no_contribution%distance     = C%large_distance
    no_contribution%row_index    = -1
    no_contribution%column_index = -1

    pivot_contribution           = no_contribution

    amount_of_mapped_points      = 0

    IF(C%data_set_is_cyclic_in_longitude) THEN
     number_of_situations = 3
    ELSE
     number_of_situations = 1
    END IF

    ! Opening the file to which the coordinates of the nearest projected points are written, which will be the content of the SID file:
    OPEN(UNIT=C%unit_scanning_file_content, FILE=TRIM(C%filename_sid_content))

    ! For each IM grid point the nearest projected GCM points are determined:
    WRITE(UNIT=*,FMT='(A)') '  The progress of the OBLIMAP scanning phase is at:'
    DO m = 1, C%NX
      IF(m >= m_message) THEN
       IF(C%oblimap_message_level == 0) THEN
        WRITE(UNIT=*,FMT='(F9.1, A    )') 100._dp * REAL(m, dp) / REAL(C%NX, dp), ' %'
        m_message = m_message + 0.10_dp * C%NX
       ELSE
        WRITE(UNIT=*,FMT='(F9.1, A, I5)') 100._dp * REAL(m, dp) / REAL(C%NX, dp), ' %,  at  m = ', m
        m_message = m_message + 0.05_dp * C%NX
       END IF
      END IF
    DO n = 1, C%NY

      IF(C%full_scanning_mode .OR. (m == 1 .AND. n == 1)) THEN
       ! For the very first point always a full scan is conducted. In case the full_scanning_mode = TRUE, the full scan is conducted at any point
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      ELSE IF(n == 1) THEN
       ! Low frequent situation (At n == 1, the starting and most left column, only a contribution of the previous row can be used if it exits)
       IF(nearest_contribution(m-1,n)%distance /= C%large_distance) THEN
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(m-1,n)
       ELSE
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       END IF
      ELSE IF(nearest_contribution(m,n-1)%distance /= C%large_distance) THEN
       ! Most frequent situation (continuing at the same row): take the contribution of that neighbour point which is located at the previous column
       do_full_scan = .FALSE.
       pivot_contribution = nearest_contribution(m,n-1)
      ELSE IF(m == 1) THEN
       ! Low frequent situation (If no neighbour contribution at the same row is found, it is not possible to try the previous row, because m == 1 is the lowest and first scanned row)
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      ELSE IF(nearest_contribution(m-1,n)%distance /= C%large_distance) THEN
       ! Second frequent situation: take the contribution of that neighbour point which is located at the previous row
       do_full_scan = .FALSE.
       pivot_contribution = nearest_contribution(m-1,n)
      ELSE
       ! Low frequent situation (no contributing neighbour point is found for this point which can serve as a pivot)
       do_full_scan = .TRUE.
       pivot_contribution = no_contribution
      END IF

      ! Initialize the contributions to inappropriate values:
      contribution(:,m,n) = no_contribution

      IF(do_full_scan) THEN
       scan_search_block_size = C%scan_search_block_size ! Actually a dummy value to avoid the compiler warning: scan_search_block_size may be used uninitialized
      ELSE
       IF(C%scan_search_block_size < -1) THEN
        ! Estimating the distance of two projected gcm points in the flat IM plane: to obtain a local estimate the distance between the pivot and its neighbour at the next row is calculated:
        IF(pivot_contribution%row_index == C%NLON) THEN
         local_gcm_grid_distance = distance_in_flat_surface( &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index - 1, pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index - 1, pivot_contribution%column_index) )
        ELSE
         local_gcm_grid_distance = distance_in_flat_surface( &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index    , pivot_contribution%column_index), &
          x_coordinates_of_gcm_grid_points(pivot_contribution%row_index + 1, pivot_contribution%column_index), &
          y_coordinates_of_gcm_grid_points(pivot_contribution%row_index + 1, pivot_contribution%column_index) )
        END IF
        ! There are situations around the pole for which different grid points might almost coincide. In that case the local_gcm_grid_distance will be very small.
        ! To avoid in such case an extreme blow up of the scan_search_block_size (or even a NaN due to a zero division) we set:
        IF(local_gcm_grid_distance < 0.1_dp) local_gcm_grid_distance = C%large_distance

        ! Calculate the local scan_search_block_size based on the estimated distance of two local projected neighbour points:
        scan_search_block_size = INT(MAX(C%dx, C%dy) / local_gcm_grid_distance) + INT(C%R_search_interpolation / local_gcm_grid_distance) + 2
       ELSE
        scan_search_block_size = C%scan_search_block_size
       END IF
      END IF

      count_contributions = 0
      count_iteratations = 0
      counter = -1
      ! The dynamic scan_search_block_size iteration: each iteration the scan_search_block_size is increased with the C%scan_search_block_size_step:
      iterate: DO WHILE(counter /= 0)
       counter = 0
       count_iteratations = count_iteratations + 1

       IF(.NOT. do_full_scan) THEN
        IF(C%scan_search_block_size == -3) THEN
         IF(count_iteratations > 1) scan_search_block_size = scan_search_block_size + C%scan_search_block_size_step
        END IF
        highest_scan_search_block_size = MAX(scan_search_block_size, highest_scan_search_block_size)
       END IF

       IF(do_full_scan) THEN
        ! In case there is no clue where to start, a full search is done:
        IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, 2(I5, A))') ' Full scan for (m, n) = (', m, ',', n, ')'
        j_start =      1
        j_end   = C%NLAT
       ELSE
        j_start = MAX(pivot_contribution%column_index - scan_search_block_size,      1)
        j_end   = MIN(pivot_contribution%column_index + scan_search_block_size, C%NLAT)

        IF(latitude_parallel_to_grid_numbers) THEN
         ! Due to the spreading of points close to the South Pole in longitudinal direction at the low latitude edge of the grid, a full longitude scan is done:
         IF(j_start <= C%fls_grid_range) THEN
          IF(lat_gcm(1,j_start) < - C%fls_latitude_border) THEN
           j_start = 1
           j_end   = MIN(C%fls_limited_lat_range, C%NLAT)
           do_full_scan = .TRUE.
          END IF
         END IF
         ! Due to the spreading of points close to the North Pole in longitudinal direction at the high latitude edge of the grid, a full longitude scan is done:
         IF(j_end > C%NLAT - C%fls_grid_range) THEN
          IF(lat_gcm(1,j_end)   >   C%fls_latitude_border) THEN
           j_start = MAX(C%NLAT - C%fls_limited_lat_range + 1, 1)
           j_end   = C%NLAT
           do_full_scan = .TRUE.
          END IF
         END IF
        ELSE
         ! Due to the spreading of points close to the South Pole in longitudinal direction at the low latitude edge of the grid, a full longitude scan is done:
         IF(j_end > C%NLAT - C%fls_grid_range) THEN
          IF(lat_gcm(1,j_end) < - C%fls_latitude_border) THEN
           j_start = MAX(C%NLAT - C%fls_limited_lat_range + 1, 1)
           j_end   = C%NLAT
           do_full_scan = .TRUE.
          END IF
         END IF
         ! Due to the spreading of points close to the North Pole in longitudinal direction at the high latitude edge of the grid, a full longitude scan is done:
         IF(j_start <= C%fls_grid_range) THEN
          IF(lat_gcm(1,j_start) >  C%fls_latitude_border) THEN
           j_start = 1
           j_end   = MIN(C%fls_limited_lat_range, C%NLAT)
           do_full_scan = .TRUE.
          END IF
         END IF
        END IF
       END IF

       ! For cyclic cases all three situations are passed. While passing for non-cyclic cases only situation 1.
       ! A data set which is periodical in the longitude direction is called a cyclic case.
       DO situation = 1, number_of_situations
        IF(situation == 1) THEN
         IF(do_full_scan) THEN
          ! In case there is no clue where to start, a full search is done:
          i_start =      1
          i_end   = C%NLON
         ELSE
          ! A quick search within a local block will be done:
          i_start = MAX(pivot_contribution%row_index    - scan_search_block_size,      1)
          i_end   = MIN(pivot_contribution%row_index    + scan_search_block_size, C%NLON)
         END IF
         mask(i_start:i_end,j_start:j_end) = .TRUE.
         IF(count_iteratations > 1) THEN
          ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
          i_start_previous_iteration = MAX(pivot_contribution%row_index    - (scan_search_block_size - C%scan_search_block_size_step),      1)
          i_end_previous_iteration   = MIN(pivot_contribution%row_index    + (scan_search_block_size - C%scan_search_block_size_step), C%NLON)
         END IF
        ELSE IF(situation == 2) THEN
         IF(pivot_contribution%row_index + scan_search_block_size > C%NLON) THEN
          ! Search for contributions at the west side of the grid if the east side of the grid has been reached:
          i_start = 1
          i_end   = pivot_contribution%row_index + scan_search_block_size - C%NLON
          mask(i_start:i_end,j_start:j_end) = .TRUE.
          IF(count_iteratations > 1) THEN
           ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
           i_start_previous_iteration = 1
           i_end_previous_iteration   = pivot_contribution%row_index + (scan_search_block_size - C%scan_search_block_size_step) - C%NLON
          END IF
         ELSE
          i_start = 1
          i_end   = 0  ! This will immediately stop the ij-loop (note that therefore the mask will not be evaluated in the ij-loop with an i_end = 0 which would be out of range)
         END IF
        ELSE IF(situation == 3) THEN
         IF(pivot_contribution%row_index - scan_search_block_size <      1) THEN
          ! Search for contributions at the east side of the grid if the west side of the grid has been reached:
          i_start = pivot_contribution%row_index - scan_search_block_size + C%NLON
          i_end   = C%NLON
          mask(i_start:i_end,j_start:j_end) = .TRUE.
          IF(count_iteratations > 1) THEN
           ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
           i_start_previous_iteration = pivot_contribution%row_index - (scan_search_block_size - C%scan_search_block_size_step) + C%NLON
           i_end_previous_iteration   = C%NLON
          END IF
         ELSE
          i_start = 1
          i_end   = 0  ! This will immediately stop the ij-loop (note that therefore the mask will not be evaluated in the ij-loop with an i_end = 0 which would be out of range)
         END IF
        END IF

        IF(count_iteratations > 1) THEN
         ! Mask the part which has been scanned in the previous iterations to FALSE so this area can be skipped this iteration:
         j_start_previous_iteration = MAX(pivot_contribution%column_index - (scan_search_block_size - C%scan_search_block_size_step),      1)
         j_end_previous_iteration   = MIN(pivot_contribution%column_index + (scan_search_block_size - C%scan_search_block_size_step), C%NLAT)
         mask(i_start_previous_iteration:i_end_previous_iteration,j_start_previous_iteration:j_end_previous_iteration) = .FALSE.
        END IF

        ! See equation (2.19) in Reerink et al. (2010):
        DO i = i_start, i_end
        DO j = j_start, j_end
          IF(mask(i,j)) THEN

           projected_gcm%row_index    = i
           projected_gcm%column_index = j
           projected_gcm%distance     = distance_in_flat_surface(x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j), &
                                                                 x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n))

           IF(projected_gcm%distance <= C%R_search_interpolation .AND. projected_gcm%distance > 0._dp) THEN
            count_contributions = count_contributions + 1
            maximum_contributions = MAX(maximum_contributions, count_contributions)
            IF(count_contributions > max_size) THEN
             WRITE(UNIT=*,FMT='(/A, 2(A, I10)/, 2A/)') C%OBLIMAP_ERROR, ' the array contribution is not allocated properly [scan_with_radius_method_gcm_to_im], number of contributions =', count_contributions, ', max size = ', max_size, &
                                                       '                Increasing the "oblimap_allocate_factor_config" significantly, by a factor of 10 to 1000 or more, might solve this issue. Adjust your config file: ', TRIM(C%config_filename)
             STOP
            END IF
            contribution(count_contributions,m,n) = projected_gcm
            counter = counter + 1
           END IF

          END IF
        END DO
        END DO
        ! Leaving the do loop and the do while loop immediately if a full search was carried out:
        IF(do_full_scan) EXIT iterate
       END DO
       ! Leaving the do while loop immediately if a full search was carried out or a fixed scan_search_block_size is used:
       IF(do_full_scan .OR. (C%scan_search_block_size /= -3)) EXIT iterate
      END DO iterate

      IF(count_contributions == 0) THEN
       IF(C%oblimap_message_level > 2) WRITE(UNIT=*,FMT='(2A, F12.2, A, 2(I5, A))') TRIM(C%OBLIMAP_WARNING), ' from scan_with_radius_method_gcm_to_im(): There are 0 points within C%R_search_interpolation = ', C%R_search_interpolation, ' for point (m, n) = (', m, ',', n, ')'

       nearest_contribution(m,n) = no_contribution

       number_points_no_contribution = number_points_no_contribution + 1
      ELSE
       ! The nearest contribution is selected:
       nearest_contribution(m,n) = contribution(MINLOC(contribution(:,m,n)%distance, 1),m,n)

       WRITE(UNIT=C%unit_scanning_file_content, FMT='(3I6)', ADVANCE='NO') m, n, count_contributions
       DO loop = 1, count_contributions
        WRITE(UNIT=C%unit_scanning_file_content, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop,m,n)%row_index, contribution(loop,m,n)%column_index, contribution(loop,m,n)%distance
       END DO
       WRITE(UNIT=C%unit_scanning_file_content, FMT='(A)') ''
       amount_of_mapped_points = amount_of_mapped_points + 1
      END IF

    END DO
    END DO

    IF(C%scan_search_block_size == -3) highest_scan_search_block_size = highest_scan_search_block_size - 2
    IF(C%oblimap_message_level > 0) WRITE(UNIT=*,FMT='(/A, I6/)') ' The highest dynamic scan_search_block_size was: ', highest_scan_search_block_size

    DEALLOCATE(contribution)

    ! Closing the the SID file:
    CLOSE(UNIT=C%unit_scanning_file_content)

    ! Output: -
    CALL write_sid_file(advised_scan_parameter, highest_scan_search_block_size, amount_of_mapped_points, number_points_no_contribution, maximum_contributions, gcm_to_im_direction = .TRUE.)

    IF(maximum_contributions > max_size) THEN
     WRITE(UNIT=*, FMT='(/3A       )') C%OBLIMAP_ERROR, ' scan_with_radius_method_gcm_to_im(): in the config file: ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='(  A, F5.1/)') '                The oblimap_allocate_factor_config should be increased to ', 1.1_dp * maximum_contributions / REAL(max_size)
     STOP
    END IF
  END SUBROUTINE scan_with_radius_method_gcm_to_im



  ! -----------------------------------------------------------------------------
  ! ROUTINES WHICH SCAN THE CONTRIBUTING POINTS FOR INTERPOLATION OF IM TO GCM
  ! -----------------------------------------------------------------------------


  SUBROUTINE scan_with_quadrant_method_im_to_gcm(mapping_participation_mask, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
    ! This routine selects the contributing points for each target grid point, by searching with the quadrant method. First
    ! the coordinates of the IM grid points are projected with the inverse oblique stereographic projection to the GCM
    ! coordinates. Thereafter with these projected coordinates the distances of the projected points relative to each target
    ! grid point are calculated and used to select the nearest contributing grid points. The IM-grid indices of the
    ! contributing points and the relative distance to 'their' target grid point are stored by writing them to the
    ! C%sid_filename file. With the indices and the distances of the contributing points the IM fields can be
    ! mapped fast and simultaneously on to the GCM grid.
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    IMPLICIT NONE

    ! Input variables:
    INTEGER,       DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: mapping_participation_mask        ! A mask for points which participate (mask = 1) in the mapping, so which fall within the mapped area.
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: x_coordinates_of_im_grid_points   ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: y_coordinates_of_im_grid_points   ! The y-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lon_gcm                           ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lat_gcm                           ! The latitude  coordinates (degrees) of the GCM grid points
    TYPE(oblimap_scan_parameter_type)        , INTENT(IN) :: advised_scan_parameter            ! The struct containing the crucial scan parameters.

    ! Local variables:
    REAL(dp)                                              :: minimum_im_grid_distance          ! Just MIN(C%dx, C%dy)
    REAL(dp),      DIMENSION(  C%NX  ,C%NY  )             :: lon_coordinates_of_im_grid_points
    REAL(dp),      DIMENSION(  C%NX  ,C%NY  )             :: lat_coordinates_of_im_grid_points
    INTEGER                                               :: i, j
    INTEGER                                               :: m, n
    REAL(dp)                                              :: i_message = 0._dp
    INTEGER                                               :: loop
    INTEGER                                               :: count_iteratations                ! Counting the DO WHILE iterations
    INTEGER                                               :: counter                           ! This counter counts each time a nearer point is found, if it stays zero no more points are found by extending the search block
    INTEGER                                               :: scan_search_block_size
    INTEGER                                               :: highest_scan_search_block_size = 0
    INTEGER                                               :: amount_of_mapped_points
    INTEGER                                               :: number_points_no_contribution = 0 ! Number of points for which no any contribution is found
    INTEGER                                               :: count_contributions
    INTEGER                                               :: quadrant                          ! The quadrant I, II, III or IV relative to a GCM grid point
    TYPE(triplet)                                         :: projected_im                      ! Projected GCM point on S'
    TYPE(triplet), DIMENSION(4,C%NLON,C%NLAT)             :: contribution                      ! Nearest projected IM point in quadrant (DIM=I,II,III or IV) in S, relative to the GCM grid point
    TYPE(triplet)                                         :: no_contribution                   ! In case there are no contributions, the nearest contribution elements are set to some specific values: the distance to a huge number, and the indices are set to -1
    TYPE(triplet)                                         :: pivot_contribution                ! The selected pivot contribution, this pivot determines the position of the scan block
    TYPE(triplet), DIMENSION(  C%NLON,C%NLAT)             :: nearest_contribution              ! Keep the nearest projected IM contribution for each GCM grid point
    LOGICAL                                               :: do_full_scan                      ! Do a full scan off all projected departing grid points for this destination grid point
    INTEGER                                               :: m_start                           ! The starting m indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: m_end                             ! The ending   m indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: n_start                           ! The starting n indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: n_end                             ! The ending   n indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: m_start_previous_iteration        ! The m_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: m_end_previous_iteration          ! The m_end   of the previous iteration in the WHILE-loop
    INTEGER                                               :: n_start_previous_iteration        ! The n_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: n_end_previous_iteration          ! The n_end   of the previous iteration in the WHILE-loop
    LOGICAL,       DIMENSION(  C%NX  ,C%NY  )             :: mask

    minimum_im_grid_distance = MIN(C%dx, C%dy)

    ! Projection of the IM coordinates to the GCM coordinates with the inverse oblique stereographic projection:
    ! Output: lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points
    CALL projecting_the_im_xy_coordinates_to_lonlat(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, &
                                                    lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)

    no_contribution%distance     = C%large_distance
    no_contribution%row_index    = -1
    no_contribution%column_index = -1

    pivot_contribution           = no_contribution

    amount_of_mapped_points      = 0

    ! Opening the file to which the coordinates of the nearest projected points are written, which will be the content of the SID file:
    OPEN(UNIT=C%unit_scanning_file_content, FILE=TRIM(C%filename_sid_content))

    ! For each GCM grid point the four nearest projected IM points are determined:
    WRITE(UNIT=*,FMT='(A)') '  The progress of the OBLIMAP scanning phase is at:'
    DO i = 1, C%NLON
      IF(i >= i_message) THEN
       IF(C%oblimap_message_level == 0) THEN
        WRITE(UNIT=*,FMT='(F9.1, A    )') 100._dp * REAL(i, dp) / REAL(C%NLON, dp), ' %'
        i_message = i_message + 0.10_dp * C%NLON
       ELSE
        WRITE(UNIT=*,FMT='(F9.1, A, I5)') 100._dp * REAL(i, dp) / REAL(C%NLON, dp), ' %,  at  i = ', i
        i_message = i_message + 0.05_dp * C%NLON
       END IF
      END IF
    DO j = 1, C%NLAT
      IF(mapping_participation_mask(i,j) == 1) THEN

       ! Because of the mapping_participation_mask only for a part of the GCM grid points there are contributions available, however earlier scanned
       ! grid neighbour contributions at the same row or same column are used as much as possible to save cpu:
       IF(C%full_scanning_mode .OR. (i == 1 .AND. j == 1)) THEN
        ! For the very first point always a full scan is conducted. In case the full_scanning_mode = TRUE, the full scan is conducted at any point
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       ELSE IF(j == 1) THEN
        ! Low frequent situation (At j == 1, the starting and most left column, only a contribution of the previous row can be used if it exits)
        IF(nearest_contribution(i-1,j)%distance /= C%large_distance) THEN
         do_full_scan = .FALSE.
         pivot_contribution = nearest_contribution(i-1,j)
        ELSE
         do_full_scan = .TRUE.
         pivot_contribution = no_contribution
        END IF
       ELSE IF(nearest_contribution(i,j-1)%distance /= C%large_distance) THEN
        ! Most frequent situation (continuing at the same row): take the contribution of that neighbour point which is located at the previous column
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(i,j-1)
       ELSE IF(i == 1) THEN
        ! Low frequent situation (If no neighbour contribution at the same row is found, it is not possible to try the previous row, because i == 1 is the lowest and first scanned row)
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       ELSE IF(nearest_contribution(i-1,j)%distance /= C%large_distance) THEN
        ! Second frequent situation: take the contribution of that neighbour point which is located at the previous row
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(i-1,j)
       ELSE
        ! Low frequent situation (no contributing neighbour point is found for this point which can serve as a pivot)
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       END IF

       ! Initialize the contributions to inappropriate values:
       contribution(:,i,j) = no_contribution

       IF(C%scan_search_block_size < -1) THEN
        ! Calculate the local scan_search_block_size based on the estimated distance of two local projected neighbour points:
        IF(j == C%NLAT) THEN
         scan_search_block_size = INT(distance_in_im_plane_between_two_gcm_points(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i,j-1), lat_gcm(i,j-1)) / minimum_im_grid_distance) + 2
        ELSE
         scan_search_block_size = INT(distance_in_im_plane_between_two_gcm_points(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i,j+1), lat_gcm(i,j+1)) / minimum_im_grid_distance) + 2
        END IF
       ELSE
        scan_search_block_size = C%scan_search_block_size
       END IF

       count_iteratations = 0
       counter = -1
       ! The dynamic scan_search_block_size iteration: each iteration the scan_search_block_size is increased with the C%scan_search_block_size_step:
      iterate: DO WHILE(counter /= 0)
        counter = 0
        count_iteratations = count_iteratations + 1

        IF(do_full_scan) THEN
         ! In case there is no clue where to start, a full search is done:
         IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, 2(I5, A))') ' Full scan for (i, j) = (', i, ',', j, ')'
         m_start =    1
         m_end   = C%NX
         n_start =    1
         n_end   = C%NY
        ELSE
         IF(C%scan_search_block_size == -3) THEN
          IF(count_iteratations > 1) scan_search_block_size = scan_search_block_size + C%scan_search_block_size_step
         END IF
         highest_scan_search_block_size = MAX(scan_search_block_size, highest_scan_search_block_size)

         ! A quick search within a local block will be done:
         m_start = MAX(pivot_contribution%row_index    - scan_search_block_size,     1)
         m_end   = MIN(pivot_contribution%row_index    + scan_search_block_size,  C%NX)
         n_start = MAX(pivot_contribution%column_index - scan_search_block_size,     1)
         n_end   = MIN(pivot_contribution%column_index + scan_search_block_size,  C%NY)
        END IF
        mask(m_start:m_end,n_start:n_end) = .TRUE.

        ! Mask the part which is already scanned in the previous iterations to FALSE so this area can be skipped in the next iterations:
        IF(count_iteratations > 1) THEN
         m_start_previous_iteration = MAX(pivot_contribution%row_index    - (scan_search_block_size - C%scan_search_block_size_step),    1)
         m_end_previous_iteration   = MIN(pivot_contribution%row_index    + (scan_search_block_size - C%scan_search_block_size_step), C%NX)
         n_start_previous_iteration = MAX(pivot_contribution%column_index - (scan_search_block_size - C%scan_search_block_size_step),    1)
         n_end_previous_iteration   = MIN(pivot_contribution%column_index + (scan_search_block_size - C%scan_search_block_size_step), C%NY)
         mask(m_start_previous_iteration:m_end_previous_iteration,n_start_previous_iteration:n_end_previous_iteration) = .FALSE.
        END IF

        ! See equation (2.17) in Reerink et al. (2010):
        DO m = m_start, m_end
        DO n = n_start, n_end
          IF(mask(m,n)) THEN

           ! Determine the quadrant in which the projected point lies relative to the considered grid point:
           ! Output: quadrant
           CALL find_quadrant_around_GCM_grid_point(lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n), &
                                                    lon_gcm(i,j), lat_gcm(i,j), quadrant)

           ! Determine in the curved plane S the distance between the projected IM coordinates relative to the considered GCM grid point:
           projected_im%row_index    = m
           projected_im%column_index = n
           projected_im%distance     = distance_over_earth_surface(lon_gcm(i,j), lat_gcm(i,j), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
           ! In case the projected point coincides with the grid point we put it at the very close distance of 1 centimeter, preventing devision by zero:
           IF(projected_im%distance == 0._dp) projected_im%distance = 0.01_dp

           ! Select the in S projected IM point with the shortest distance to the considered GCM grid point in this quadrant,
           ! and keep this distance and the IM-grid indices of this IM point in S':
           IF(projected_im%distance < contribution(quadrant,i,j)%distance) THEN
            contribution(quadrant,i,j) = projected_im
            counter = counter + 1
           END IF

          END IF
        END DO
        END DO
        ! Leaving the loop immediately if a full search was carried out or a fixed scan_search_block_size is used:
        IF(do_full_scan .OR. (C%scan_search_block_size /= -3)) EXIT iterate
       END DO iterate

       count_contributions = 4
       DO loop = 1, 4
        IF(contribution(loop,i,j)%distance == C%large_distance) count_contributions = count_contributions - 1
       END DO

       IF(count_contributions == 0) THEN
        IF(C%oblimap_message_level > 2) WRITE(UNIT=*, FMT='(2A, 2(I5, A))') TRIM(C%OBLIMAP_WARNING), ' from scan_with_quadrant_method_im_to_gcm():  In four quadrants no single point is found for point (i, j) = (', i, ',', j, ')'

        nearest_contribution(i,j) = no_contribution

        number_points_no_contribution = number_points_no_contribution + 1
       ELSE
        ! The nearest contribution is selected:
        nearest_contribution(i,j) = contribution(MINLOC(contribution(:,i,j)%distance, 1),i,j)

        WRITE(UNIT=C%unit_scanning_file_content, FMT='(3I6)', ADVANCE='NO') i, j, count_contributions
        DO loop = 1, 4
         ! Filter the appropriate contributions (leave out the quadrants in which no contributing point is found, e.g. at the grid border):
         IF(contribution(loop,i,j)%distance /= C%large_distance) THEN
          WRITE(UNIT=C%unit_scanning_file_content, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop,i,j)%row_index, contribution(loop,i,j)%column_index, contribution(loop,i,j)%distance
         END IF
        END DO
        WRITE(UNIT=C%unit_scanning_file_content, FMT='(A)') ''
        amount_of_mapped_points = amount_of_mapped_points + 1
       END IF

      ELSE
       nearest_contribution(i,j) = no_contribution
      END IF
    END DO
    END DO

    IF(C%scan_search_block_size == -3) highest_scan_search_block_size = highest_scan_search_block_size - 2
    IF(C%oblimap_message_level > 0) WRITE(UNIT=*,FMT='(/A, I6/)') ' The highest dynamic scan_search_block_size was: ', highest_scan_search_block_size

    ! Closing the the SID file:
    CLOSE(UNIT=C%unit_scanning_file_content)

    ! Output: -
    CALL write_sid_file(advised_scan_parameter, highest_scan_search_block_size, amount_of_mapped_points, number_points_no_contribution, maximum_contributions = 4, gcm_to_im_direction = .FALSE.)
  END SUBROUTINE scan_with_quadrant_method_im_to_gcm



  SUBROUTINE scan_with_radius_method_im_to_gcm(mapping_participation_mask, x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, advised_scan_parameter)
    ! This routine selects the contributing points for each target grid point, by searching with the radius method. First
    ! the coordinates of the IM grid points are projected with the inverse oblique stereographic projection to the GCM
    ! coordinates. Thereafter with these projected coordinates the distances of the projected points relative to each target
    ! grid point are calculated and used to select the nearest contributing grid points. The IM-grid indices of the
    ! contributing points and the relative distance to 'their' target grid point are stored by writing them to the
    ! C%sid_filename file. With the indices and the distances of the contributing points the IM fields can be
    ! mapped fast and simultaneously on to the GCM grid.
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    IMPLICIT NONE

    ! Input variables:
    INTEGER,       DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: mapping_participation_mask        ! A mask for points which participate (mask = 1) in the mapping, so which fall within the mapped area.
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: x_coordinates_of_im_grid_points   ! The x-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NX,  C%NY  ), INTENT(IN) :: y_coordinates_of_im_grid_points   ! The y-coordinates of the IM points in S'
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lon_gcm                           ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp),      DIMENSION(  C%NLON,C%NLAT), INTENT(IN) :: lat_gcm                           ! The latitude  coordinates (degrees) of the GCM grid points
    TYPE(oblimap_scan_parameter_type)        , INTENT(IN) :: advised_scan_parameter            ! The struct containing the crucial scan parameters.

    ! Local variables:
    INTEGER                                               :: max_size
    INTEGER                                               :: status
    REAL(dp)                                              :: minimum_im_grid_distance          ! Just MIN(C%dx, C%dy)
    REAL(dp),      DIMENSION(  C%NX,  C%NY  )             :: lon_coordinates_of_im_grid_points
    REAL(dp),      DIMENSION(  C%NX,  C%NY  )             :: lat_coordinates_of_im_grid_points
    INTEGER                                               :: i, j
    INTEGER                                               :: m, n
    REAL(dp)                                              :: i_message = 0._dp
    INTEGER                                               :: amount_of_mapped_points
    INTEGER                                               :: number_points_no_contribution = 0 ! Number of points for which no any contribution is found
    INTEGER                                               :: count_contributions
    INTEGER                                               :: maximum_contributions = 0
    INTEGER                                               :: loop
    INTEGER                                               :: count_iteratations                ! Counting the DO WHILE iterations
    INTEGER                                               :: counter                           ! This counter counts each time a nearer point is found, if it stays zero no more points are found by extending the search block
    INTEGER                                               :: scan_search_block_size
    INTEGER                                               :: highest_scan_search_block_size = 0
    TYPE(triplet)                                         :: projected_im                      ! Projected IM point on S
    TYPE(triplet), DIMENSION(:,:,:), ALLOCATABLE          :: contribution                      ! Projected IM points on S within C%R_search_interpolation, relative to the GCM grid point
    TYPE(triplet)                                         :: no_contribution                   ! In case there are no contributions, the nearest contribution elements are set to some specific values: the distance to a huge number, and the indices are set to -1
    TYPE(triplet)                                         :: pivot_contribution                ! The selected pivot contribution, this pivot determines the position of the scan block
    TYPE(triplet), DIMENSION(  C%NLON,C%NLAT)             :: nearest_contribution              ! Keep the nearest projected IM contribution for each GCM grid point
    LOGICAL                                               :: do_full_scan                      ! Do a full scan off all projected departing grid points for this destination grid point
    INTEGER                                               :: m_start                           ! The starting m indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: m_end                             ! The ending   m indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: n_start                           ! The starting n indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: n_end                             ! The ending   n indices to walk through the local block to find the nearest contributions, avoiding a search through the whole domain for each point
    INTEGER                                               :: m_start_previous_iteration        ! The m_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: m_end_previous_iteration          ! The m_end   of the previous iteration in the WHILE-loop
    INTEGER                                               :: n_start_previous_iteration        ! The n_start of the previous iteration in the WHILE-loop
    INTEGER                                               :: n_end_previous_iteration          ! The n_end   of the previous iteration in the WHILE-loop
    LOGICAL,       DIMENSION(  C%NX  ,C%NY  )             :: mask                              ! The mask which is used in the dynamic scanning, setting the interior of the already scanned block to FALSE

    max_size = CEILING(MAX(4._dp * C%pi * (C%R_search_interpolation / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)), &
                                   C%pi * (C%R_search_interpolation / 1000._dp)**2 / ((C%dx / 1000._dp) * (C%dy / 1000._dp)))  * C%oblimap_allocate_factor)
    ALLOCATE(contribution(max_size,C%NLON,C%NLAT), STAT=status)
    IF(status /= 0) THEN
     WRITE(UNIT=*, FMT='(/2A, /2(A, I8), A, F16.3, A/)') &
      C%OBLIMAP_ERROR, ' message from: scan_with_radius_method_im_to_gcm():  The allocation of the "contribution struct" exceeds your system allocation capacity.', &
      '                The combination of NLON_config = ', C%NLON, ', and NLAT_config = ', C%NLAT, ', with R_search_interpolation_config = ', C%R_search_interpolation, ' is to large. Reduce e.g. the size of R_search_interpolation_config'
     STOP
    END IF

    minimum_im_grid_distance = MIN(C%dx, C%dy)

    ! Projection of the IM coordinates to the GCM coordinates with the inverse oblique stereographic projection:
    ! Output: lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points
    CALL projecting_the_im_xy_coordinates_to_lonlat(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, &
                                                    lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)

    no_contribution%distance     = C%large_distance
    no_contribution%row_index    = -1
    no_contribution%column_index = -1

    pivot_contribution           = no_contribution

    amount_of_mapped_points      = 0

    ! Opening the file to which the coordinates of the nearest projected points are written, which will be the content of the SID file:
    OPEN(UNIT=C%unit_scanning_file_content, FILE=TRIM(C%filename_sid_content))

    ! For each GCM grid point the nearest projected IM points are determined:
    WRITE(UNIT=*,FMT='(A)') '  The progress of the OBLIMAP scanning phase is at:'
    DO i = 1, C%NLON
      IF(i >= i_message) THEN
       IF(C%oblimap_message_level == 0) THEN
        WRITE(UNIT=*,FMT='(F9.1, A    )') 100._dp * REAL(i, dp) / REAL(C%NLON, dp), ' %'
        i_message = i_message + 0.10_dp * C%NLON
       ELSE
        WRITE(UNIT=*,FMT='(F9.1, A, I5)') 100._dp * REAL(i, dp) / REAL(C%NLON, dp), ' %,  at  i = ', i
        i_message = i_message + 0.05_dp * C%NLON
       END IF
      END IF
    DO j = 1, C%NLAT
      IF(mapping_participation_mask(i,j) == 1) THEN

       ! Because of the mapping_participation_mask only for a part of the GCM grid points there are contributions available, however earlier scanned
       ! grid neighbour contributions at the same row or same column are used as much as possible to save cpu:
       IF(C%full_scanning_mode .OR. (i == 1 .AND. j == 1)) THEN
        ! For the very first point always a full scan is conducted. In case the full_scanning_mode = TRUE, the full scan is conducted at any point
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       ELSE IF(j == 1) THEN
        ! Low frequent situation (At j == 1, the starting and most left column, only a contribution of the previous row can be used if it exits)
        IF(nearest_contribution(i-1,j)%distance /= C%large_distance) THEN
         do_full_scan = .FALSE.
         pivot_contribution = nearest_contribution(i-1,j)
        ELSE
         do_full_scan = .TRUE.
         pivot_contribution = no_contribution
        END IF
       ELSE IF(nearest_contribution(i,j-1)%distance /= C%large_distance) THEN
        ! Most frequent situation (continuing at the same row): take the contribution of that neighbour point which is located at the previous column
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(i,j-1)
       ELSE IF(i == 1) THEN
        ! Low frequent situation (If no neighbour contribution at the same row is found, it is not possible to try the previous row, because i == 1 is the lowest and first scanned row)
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       ELSE IF(nearest_contribution(i-1,j)%distance /= C%large_distance) THEN
        ! Second frequent situation: take the contribution of that neighbour point which is located at the previous row
        do_full_scan = .FALSE.
        pivot_contribution = nearest_contribution(i-1,j)
       ELSE
        ! Low frequent situation (no contributing neighbour point is found for this point which can serve as a pivot)
        do_full_scan = .TRUE.
        pivot_contribution = no_contribution
       END IF

       ! Initialize the contributions to inappropriate values:
       contribution(:,i,j) = no_contribution

       IF(C%scan_search_block_size < -1) THEN
        ! Calculate the local scan_search_block_size based on the estimated distance of two local projected neighbour points:
        IF(j == C%NLAT) THEN
         scan_search_block_size = INT(distance_in_im_plane_between_two_gcm_points(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i,j-1), lat_gcm(i,j-1)) / minimum_im_grid_distance) &
                                   + INT(C%R_search_interpolation / minimum_im_grid_distance) + 2
        ELSE
         scan_search_block_size = INT(distance_in_im_plane_between_two_gcm_points(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i,j+1), lat_gcm(i,j+1)) / minimum_im_grid_distance) &
                                   + INT(C%R_search_interpolation / minimum_im_grid_distance) + 2
        END IF
       ELSE
        scan_search_block_size = C%scan_search_block_size
       END IF

       count_contributions = 0
       count_iteratations = 0
       counter = -1
       ! The dynamic scan_search_block_size iteration: each iteration the scan_search_block_size is increased with the C%scan_search_block_size_step:
      iterate: DO WHILE(counter /= 0)
        counter = 0
        count_iteratations = count_iteratations + 1

        IF(do_full_scan) THEN
         ! In case there is no clue where to start, a full search is done:
         IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, 2(I5, A))') ' Full scan for (i, j) = (', i, ',', j, ')'
         m_start =    1
         m_end   = C%NX
         n_start =    1
         n_end   = C%NY
        ELSE
         IF(C%scan_search_block_size == -3) THEN
          IF(count_iteratations > 1) scan_search_block_size = scan_search_block_size + C%scan_search_block_size_step
         END IF
         highest_scan_search_block_size = MAX(scan_search_block_size, highest_scan_search_block_size)

         ! A quick search within a local block will be done:
         m_start = MAX(pivot_contribution%row_index    - scan_search_block_size,     1)
         m_end   = MIN(pivot_contribution%row_index    + scan_search_block_size,  C%NX)
         n_start = MAX(pivot_contribution%column_index - scan_search_block_size,     1)
         n_end   = MIN(pivot_contribution%column_index + scan_search_block_size,  C%NY)
        END IF
        mask(m_start:m_end,n_start:n_end) = .TRUE.

        ! Mask the part which is already scanned in the previous iterations to FALSE so this area can be skipped in the next iterations:
        IF(count_iteratations > 1) THEN
         m_start_previous_iteration = MAX(pivot_contribution%row_index    - (scan_search_block_size - C%scan_search_block_size_step),    1)
         m_end_previous_iteration   = MIN(pivot_contribution%row_index    + (scan_search_block_size - C%scan_search_block_size_step), C%NX)
         n_start_previous_iteration = MAX(pivot_contribution%column_index - (scan_search_block_size - C%scan_search_block_size_step),    1)
         n_end_previous_iteration   = MIN(pivot_contribution%column_index + (scan_search_block_size - C%scan_search_block_size_step), C%NY)
         mask(m_start_previous_iteration:m_end_previous_iteration,n_start_previous_iteration:n_end_previous_iteration) = .FALSE.
        END IF

        ! See equation (2.19) in Reerink et al. (2010):
        DO m = m_start, m_end
        DO n = n_start, n_end
          IF(mask(m,n)) THEN

           projected_im%row_index    = m
           projected_im%column_index = n
           projected_im%distance     = distance_over_earth_surface(lon_gcm(i,j), lat_gcm(i,j), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))

           IF(projected_im%distance <= C%R_search_interpolation .AND. projected_im%distance > 0._dp) THEN
            count_contributions = count_contributions + 1
            maximum_contributions = MAX(maximum_contributions, count_contributions)
            IF(count_contributions > max_size) THEN
             WRITE(UNIT=*,FMT='(/A, 2(A, I10)/, 2A/)') C%OBLIMAP_ERROR, ' the array contribution is not allocated properly [scan_with_radius_method_im_to_gcm], number of contributions =', count_contributions, ', max size = ', max_size, &
                                                       '                Increasing the "oblimap_allocate_factor_config" significantly, by a factor of 10 to 1000 or more, might solve this issue. Adjust your config file: ', TRIM(C%config_filename)
             STOP
            END IF
            contribution(count_contributions,i,j) = projected_im
            counter = counter + 1
           END IF

          END IF
        END DO
        END DO
        ! Leaving the loop immediately if a full search was carried out or a fixed scan_search_block_size is used:
        IF(do_full_scan .OR. (C%scan_search_block_size /= -3)) EXIT iterate
       END DO iterate

       IF(count_contributions == 0) THEN
        IF(C%oblimap_message_level > 2) WRITE(UNIT=*,FMT='(2A, F12.2, A, 2(I5, A))') TRIM(C%OBLIMAP_WARNING), ' from scan_with_radius_method_im_to_gcm: There are 0 points within C%R_search_interpolation = ', C%R_search_interpolation, ' for point (i, j) = (', i, ',', j, ')'

        nearest_contribution(i,j) = no_contribution

        number_points_no_contribution = number_points_no_contribution + 1
       ELSE
        ! The nearest contribution is selected:
        nearest_contribution(i,j) = contribution(MINLOC(contribution(:,i,j)%distance, 1),i,j)

        WRITE(UNIT=C%unit_scanning_file_content, FMT='(3I6)', ADVANCE='NO') i, j, count_contributions
        DO loop = 1, count_contributions
         WRITE(UNIT=C%unit_scanning_file_content, FMT='(2I6,E23.15)', ADVANCE='NO') contribution(loop,i,j)%row_index, contribution(loop,i,j)%column_index, contribution(loop,i,j)%distance
        END DO
        WRITE(UNIT=C%unit_scanning_file_content, FMT='(A)') ''
        amount_of_mapped_points = amount_of_mapped_points + 1
       END IF

      ELSE
       nearest_contribution(i,j) = no_contribution
      END IF
    END DO
    END DO

    IF(C%scan_search_block_size == -3) highest_scan_search_block_size = highest_scan_search_block_size - 2
    IF(C%oblimap_message_level > 0) WRITE(UNIT=*,FMT='(/A, I6/)') ' The highest dynamic scan_search_block_size was: ', highest_scan_search_block_size

    DEALLOCATE(contribution)

    ! Closing the the SID file:
    CLOSE(UNIT=C%unit_scanning_file_content)

    ! Output: -
    CALL write_sid_file(advised_scan_parameter, highest_scan_search_block_size, amount_of_mapped_points, number_points_no_contribution, maximum_contributions, gcm_to_im_direction = .FALSE.)

    IF(maximum_contributions > max_size) THEN
     WRITE(UNIT=*, FMT='(/3A       )') C%OBLIMAP_ERROR, ' scan_with_radius_method_im_to_gcm(): in the config file: ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='(  A, F5.1/)') '                The oblimap_allocate_factor_config should be increased to ', 1.1_dp * maximum_contributions / REAL(max_size)
     STOP
    END IF
  END SUBROUTINE scan_with_radius_method_im_to_gcm



  ! --------------------
  ! SUPPORTING ROUTINES
  ! --------------------

  SUBROUTINE write_sid_file(advised_scan_parameter, highest_scan_search_block_size, amount_of_mapped_points, number_points_no_contribution, maximum_contributions, gcm_to_im_direction)
    ! This routine writes the SID file (the file which contains the scanned indices and distances).
    USE oblimap_configuration_module, ONLY: dp, C, oblimap_scan_parameter_type
    IMPLICIT NONE

    ! Input variables:
    TYPE(oblimap_scan_parameter_type) , INTENT(IN) :: advised_scan_parameter
    INTEGER                           , INTENT(IN) :: highest_scan_search_block_size
    INTEGER                           , INTENT(IN) :: amount_of_mapped_points
    INTEGER                           , INTENT(IN) :: number_points_no_contribution
    INTEGER                           , INTENT(IN) :: maximum_contributions
    LOGICAL                           , INTENT(IN) :: gcm_to_im_direction                   ! This variable has to be TRUE for the GCM -> IM mapping direction, and FALSE vice versa.

    ! Local variables:
    INTEGER                                        :: unit_number = 107

    ! Opening the SID file:
    OPEN(UNIT=unit_number, FILE=TRIM(C%sid_filename))

    ! Writing the header of the C%sid_filename file:
    WRITE(UNIT=unit_number,   FMT='( A        )') '# Do not remove this header. The data format of this file is:'
    IF(gcm_to_im_direction) THEN
     WRITE(UNIT=unit_number,  FMT='( A        )') '#  m  n  N  N(i  j  distance)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with m = the x-axis IM grid counter (of the considered destination point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with n = the y-axis IM grid counter (of the considered destination point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with N = the number of contributions for this destination point'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with i = the longitudinal GCM grid counter (of the projetced departure point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with j = the latitudinal  GCM grid counter (of the projetced departure point)'
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A        )') '#  i  j  N  N(m  n  distance)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with i = the longitudinal GCM grid counter (of the considered destination point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with j = the latitudinal  GCM grid counter (of the considered destination point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with N = the number of contributions for this destination point'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with m = the x-axis IM grid counter (of the projetced departure point)'
     WRITE(UNIT=unit_number,  FMT='( A        )') '# with n = the y-axis IM grid counter (of the projetced departure point)'
    END IF
    WRITE(UNIT=unit_number,   FMT='( A        )') '# where the distance is the distance between the projected departure point and the destination point'
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    IF(gcm_to_im_direction) THEN
     WRITE(UNIT=unit_number,  FMT='(2A        )') '# This file is created by:  ../src/oblimap_gcm_to_im_program ', TRIM(C%config_filename)
    ELSE
     WRITE(UNIT=unit_number,  FMT='(2A        )') '# This file is created by:  ../src/oblimap_im_to_gcm_program ', TRIM(C%config_filename)
    END IF
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A        )') '# Summary of the OBLIMAP scan parameters:'
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    IF(gcm_to_im_direction) THEN
     WRITE(UNIT=unit_number,  FMT='( A, A     )') '#  gcm_input_filename_config                                 = ', TRIM(C%gcm_input_filename)
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A, A     )') '#  im_input_filename_config                                  = ', TRIM(C%im_input_filename)
    END IF
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  NLON_config                                               = ', C%NLON
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  NLAT_config                                               = ', C%NLAT
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  NX_config                                                 = ', C%NX
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  NY_config                                                 = ', C%NY
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  dx_config                                                 = ', C%dx
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  dy_config                                                 = ', C%dy
    IF(C%choice_projection_method == 'rotation_projection') THEN
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  shift_x_coordinate_rotation_projection_config             = ', C%shift_x_coordinate_rotation_projection
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  shift_y_coordinate_rotation_projection_config             = ', C%shift_y_coordinate_rotation_projection
     WRITE(UNIT=unit_number,  FMT='( A        )') '# '
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  theta_rotation_projection_config                          = ', C%radians_to_degrees * C%theta_rotation_projection
     WRITE(UNIT=unit_number,  FMT='( A        )') '# '
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  lambda_M_config                                           = ', C%radians_to_degrees * C%lambda_M
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  phi_M_config                                              = ', C%radians_to_degrees * C%phi_M
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  alpha_stereographic_config                                = ', C%radians_to_degrees * C%alpha_stereographic
     SELECT CASE(C%choice_projection_method)
     CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder','oblique_lambert_equal-area_projection_snyder')
      WRITE(UNIT=unit_number, FMT='( A, E24.16)') '#  earth_radius_config                                       = ', C%earth_radius
      WRITE(UNIT=unit_number, FMT='( A        )') '# '
     CASE('oblique_stereographic_projection_ellipsoid_snyder','oblique_lambert_equal-area_projection_ellipsoid_snyder')
      WRITE(UNIT=unit_number, FMT='( A, E24.16)') '#  ellipsoid_semi_major_axis_config                          = ', C%a
      WRITE(UNIT=unit_number, FMT='( A, E24.16)') '#  ellipsoid_eccentricity_config                             = ', C%e
     END SELECT
    END IF
    WRITE(UNIT=unit_number,   FMT='( A, A     )') '#  choice_projection_method_config                           = ', TRIM(C%choice_projection_method)
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  enable_shift_im_grid_config                               = ', C%enable_shift_im_grid
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  shift_x_coordinate_im_grid_config                         = ', C%shift_x_coordinate_im_grid
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  shift_y_coordinate_im_grid_config                         = ', C%shift_y_coordinate_im_grid
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  alternative_lambda_for_center_im_grid_config              = ', C%alternative_lambda_for_center_im_grid
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  alternative_phi_for_center_im_grid_config                 = ', C%alternative_phi_for_center_im_grid
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  unit_conversion_x_ax_config                               = ', C%unit_conversion_x_ax
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  unit_conversion_y_ax_config                               = ', C%unit_conversion_y_ax
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  use_prefabricated_im_grid_coordinates_config              = ', C%use_prefabricated_im_grid_coordinates
    IF(C%use_prefabricated_im_grid_coordinates) THEN
     WRITE(UNIT=unit_number,  FMT='( A, A     )') '#  prefabricated_im_grid_filename_config                     = ', TRIM(C%prefabricated_im_grid_filename)
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A        )') '#  prefabricated_im_grid_filename_config                     = -'
    END IF
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  level_of_automatic_oblimap_scanning_config                = ', C%level_of_automatic_oblimap_scanning
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  data_set_is_cyclic_in_longitude_config                    = ', C%data_set_is_cyclic_in_longitude
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  choice_quadrant_method_config                             = ', C%choice_quadrant_method
    IF(C%choice_quadrant_method) THEN
     WRITE(UNIT=unit_number,  FMT='( A        )') '#  R_search_interpolation_config                             = -'
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A, E24.16)') '#  R_search_interpolation_config                             = ', C%R_search_interpolation
    END IF
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  scan_search_block_size_config                             = ', C%scan_search_block_size
    WRITE(UNIT=unit_number,   FMT='( A, I9    )') '#  scan_search_block_size_step_config                        = ', C%scan_search_block_size_step
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  vincenty_method_for_ellipsoid_config                      = ', C%vincenty_method_for_ellipsoid
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    WRITE(UNIT=unit_number,   FMT='( A        )') '# Below the best estimates by OBLIMAP itself are listed. Consider to use them if they differ from the used values above.'
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  optimal data_set_is_cyclic_in_longitude_config  = ', advised_scan_parameter%data_set_is_cyclic_in_longitude
    IF(highest_scan_search_block_size == 0) THEN
     WRITE(UNIT=unit_number,  FMT='( A        )') '#  optimal search_block_size_config                = -'
    ELSE
     WRITE(UNIT=unit_number,  FMT='( A, I9    )') '#  optimal search_block_size_config                = ', highest_scan_search_block_size
    END IF
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  optimal alpha_stereographic_config              = ', advised_scan_parameter%alpha_stereographic
    WRITE(UNIT=unit_number,   FMT='( A, L     )') '#  optimal choice_quadrant_method_config           = ', advised_scan_parameter%choice_quadrant_method
    WRITE(UNIT=unit_number,   FMT='( A, E24.16)') '#  optimal R_search_interpolation_config           = ', advised_scan_parameter%R_search_interpolation
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '
    IF(gcm_to_im_direction) THEN
     WRITE(UNIT=unit_number,  FMT='( I20, A       )') maximum_contributions  , '     # The maximum number of GCM contributions used to obtain the interpolated field values for the IM points.'
     WRITE(UNIT=unit_number,  FMT='( I20, A, I8, A)') amount_of_mapped_points, '     # The number of mapped IM points from which ' , number_points_no_contribution, ' points have no contribution.'
    ELSE
     WRITE(UNIT=unit_number,  FMT='( I20, A       )') maximum_contributions  , '     # The maximum number of IM contributions used to obtain the interpolated field values for the GCM points.'
     WRITE(UNIT=unit_number,  FMT='( I20, A, I8, A)') amount_of_mapped_points, '     # The number of mapped GCM points from which ', number_points_no_contribution, ' points have no contribution.'
    END IF
    WRITE(UNIT=unit_number,   FMT='( A        )') '# '

    ! Closing the the SID file:
    CLOSE(UNIT=unit_number)

    ! Appending the content to the header:
    CALL SYSTEM('cat '//TRIM(C%filename_sid_content)//' >> '//TRIM(C%sid_filename))
    CALL SYSTEM('rm -f '//TRIM(C%filename_sid_content))
  END SUBROUTINE write_sid_file



  SUBROUTINE projecting_the_gcm_lonlat_coordinates_to_xy(lon_gcm, lat_gcm, x_coordinates_of_gcm_grid_points, y_coordinates_of_gcm_grid_points)
    ! This routine projects the GCM coordinates on the requested plane S' which coincides with the IM grid,
    ! with an oblique stereographic or an oblique lambert equal area projection.
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_sg_projection_ellipsoid_snyder, &
      oblique_laea_projection_snyder, oblique_laea_projection_ellipsoid_snyder, rotation_projection
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm                          ! longitude coordinates (degrees) of GCM grid
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm                          ! latitude coordinates  (degrees) of GCM grid

    ! Output variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: x_coordinates_of_gcm_grid_points ! The x-coordinates of the GCM points projected on S'
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: y_coordinates_of_gcm_grid_points ! The y-coordinates of the GCM points projected on S'

    ! Local variables:
    INTEGER                                         :: i, j

    ! Determine the x,y coordinates of each GCM longitude-latitude coordinate after
    ! The oblique stereographic projection on the projection plane S'
    DO i = 1, C%NLON
    DO j = 1, C%NLAT
      ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
      SELECT CASE(C%choice_projection_method)
      CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder')
       CALL oblique_sg_projection(                                  lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      CASE('oblique_stereographic_projection_ellipsoid_snyder')
       CALL oblique_sg_projection_ellipsoid_snyder(                 lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      CASE('oblique_lambert_equal-area_projection_snyder')
       CALL oblique_laea_projection_snyder(                         lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder')
       CALL oblique_laea_projection_ellipsoid_snyder(               lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      CASE('rotation_projection')
       CALL rotation_projection(                                    lon_gcm(i,j), lat_gcm(i,j), x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j))
      END SELECT
    END DO
    END DO
  END SUBROUTINE projecting_the_gcm_lonlat_coordinates_to_xy



  SUBROUTINE projecting_the_im_xy_coordinates_to_lonlat(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, &
                                                        lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    ! This routine projects the IM coordinates on the requested plane S' which coincides with the GCM grid,
    ! with an inverse oblique stereographic or an inverse oblique lambert equal area projection.
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: inverse_oblique_sg_projection, inverse_oblique_sg_projection_snyder, inverse_oblique_sg_projection_ellipsoid_snyder, &
      inverse_oblique_laea_projection_snyder, inverse_oblique_laea_projection_ellipsoid_snyder, inverse_rotation_projection
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(IN)  :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(IN)  :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'

    ! Output variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: lon_coordinates_of_im_grid_points
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(OUT) :: lat_coordinates_of_im_grid_points

    ! Local variables:
    INTEGER                                     :: m, n

    ! Inverse projection of the IM (x,y)-coordinates to the GCM-(longitude,latitude) coordinates:
    DO m = 1, C%NX
    DO n = 1, C%NY
      ! Output: lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n)
      SELECT CASE(C%choice_projection_method)
      CASE('oblique_stereographic_projection')
       CALL inverse_oblique_sg_projection(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      CASE('oblique_stereographic_projection_snyder')
       CALL inverse_oblique_sg_projection_snyder(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      CASE('oblique_stereographic_projection_ellipsoid_snyder')
       CALL inverse_oblique_sg_projection_ellipsoid_snyder(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      CASE('oblique_lambert_equal-area_projection_snyder')
       CALL inverse_oblique_laea_projection_snyder(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder')
       CALL inverse_oblique_laea_projection_ellipsoid_snyder(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      CASE('rotation_projection')
       CALL inverse_rotation_projection(x_coordinates_of_im_grid_points(m,n), y_coordinates_of_im_grid_points(m,n), lon_coordinates_of_im_grid_points(m,n), lat_coordinates_of_im_grid_points(m,n))
      END SELECT
    END DO
    END DO
  END SUBROUTINE projecting_the_im_xy_coordinates_to_lonlat



  SUBROUTINE shifting_center_im_grid(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points)
    ! A shift is usually not necessary and not recommended in an oblique case, but can be very handy in a testing phase, or
    ! in case of a non optimal projection such as a polar projection of a non-polar region (i.e. a non-centered polar projection).
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_laea_projection_snyder, oblique_sg_projection_ellipsoid_snyder, oblique_laea_projection_ellipsoid_snyder
    IMPLICIT NONE

    ! In/Output variables:
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(INOUT) :: x_coordinates_of_im_grid_points
    REAL(dp), DIMENSION(C%NX,C%NY), INTENT(INOUT) :: y_coordinates_of_im_grid_points

    ! Local variables:
    REAL(dp)                                      :: shift_x_coordinate
    REAL(dp)                                      :: shift_y_coordinate

    ! Output: shift_x_coordinate, shift_y_coordinate
    SELECT CASE(C%choice_projection_method)
    CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder')
     CALL oblique_sg_projection(C%alternative_lambda_for_center_im_grid, C%alternative_phi_for_center_im_grid, shift_x_coordinate, shift_y_coordinate)
    CASE('oblique_stereographic_projection_ellipsoid_snyder')
     CALL oblique_sg_projection_ellipsoid_snyder(C%alternative_lambda_for_center_im_grid, C%alternative_phi_for_center_im_grid, shift_x_coordinate, shift_y_coordinate)
    CASE('oblique_lambert_equal-area_projection_snyder')
     CALL oblique_laea_projection_snyder(C%alternative_lambda_for_center_im_grid, C%alternative_phi_for_center_im_grid, shift_x_coordinate, shift_y_coordinate)
    CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder')
     CALL oblique_laea_projection_ellipsoid_snyder(C%alternative_lambda_for_center_im_grid, C%alternative_phi_for_center_im_grid, shift_x_coordinate, shift_y_coordinate)
    CASE('rotation_projection')
     ! No alternative center should be specified for the rotational projection via this option:
     shift_x_coordinate = 0._dp
     shift_y_coordinate = 0._dp
    CASE DEFAULT
     WRITE(UNIT=*, FMT='(/3A )') C%OBLIMAP_ERROR, ' shifting_center_im_grid(): in the config file: ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='( 2A/)') '                Invalid value for:  choice_projection_method_config = ', TRIM(C%choice_projection_method)
     STOP
    END SELECT

    ! There are three ways to apply a shift:
    !  1. Specifying an alternative lon-lat im center (angles in degrees)
    !  2. Specifying a shift in the IM coordinates (shift in meter)
    !  3. A combination of these two ways
    x_coordinates_of_im_grid_points = x_coordinates_of_im_grid_points + shift_x_coordinate + C%shift_x_coordinate_im_grid
    y_coordinates_of_im_grid_points = y_coordinates_of_im_grid_points + shift_y_coordinate + C%shift_y_coordinate_im_grid
  END SUBROUTINE shifting_center_im_grid



  SUBROUTINE make_mapping_participation_mask(x_coordinates_of_im_grid_points, y_coordinates_of_im_grid_points, lon_gcm, lat_gcm, mapping_participation_mask)
    ! This routine determines the GCM longitude, latitude grid points that participate
    ! in the mapping. It distinguishes with a mask between points that are projected within the
    ! IM grid domain and points which are projected outside that domain.
    !
    ! With the mapping from the IM to the GCM (the inverse oblique stereographic projection
    ! followed by the interpolation) it is necessary to know which GCM points are affected
    ! by the mapping, the rest of the GCM points keep their initial value.
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_sg_projection_ellipsoid_snyder, &
      oblique_laea_projection_snyder, oblique_laea_projection_ellipsoid_snyder, rotation_projection
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NX,  C%NY  ), INTENT(IN)  :: x_coordinates_of_im_grid_points  ! The x-coordinates of the IM points in S'
    REAL(dp), DIMENSION(C%NX,  C%NY  ), INTENT(IN)  :: y_coordinates_of_im_grid_points  ! The y-coordinates of the IM points in S'
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lon_gcm
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)  :: lat_gcm

    ! Output variables:
    INTEGER,  DIMENSION(C%NLON,C%NLAT), INTENT(OUT) :: mapping_participation_mask

    ! Local variables:
    INTEGER                                         :: i, j
    REAL(dp)                                        :: x_gcm
    REAL(dp)                                        :: y_gcm

    mapping_participation_mask = 0
    DO i = 1, C%NLON
    DO j = 1, C%NLAT
      ! Output: x_gcm, y_gcm
      SELECT CASE(C%choice_projection_method)
      CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder')
       CALL oblique_sg_projection(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      CASE('oblique_stereographic_projection_ellipsoid_snyder')
       CALL oblique_sg_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      CASE('oblique_lambert_equal-area_projection_snyder')
       CALL oblique_laea_projection_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder')
       CALL oblique_laea_projection_ellipsoid_snyder(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      CASE('rotation_projection')
       CALL rotation_projection(lon_gcm(i,j), lat_gcm(i,j), x_gcm, y_gcm)
      END SELECT

      IF(x_gcm > x_coordinates_of_im_grid_points(1,1) .AND. x_gcm < x_coordinates_of_im_grid_points(C%NX,1) .AND. &
         y_gcm > y_coordinates_of_im_grid_points(1,1) .AND. y_gcm < y_coordinates_of_im_grid_points(1,C%NY)        ) mapping_participation_mask(i,j) = 1
    END DO
    END DO
  END SUBROUTINE make_mapping_participation_mask



  SUBROUTINE find_quadrant_around_IM_grid_point(x_value_projected_point,  y_value_projected_point, &
                                                x_value_IM_grid_point, y_value_IM_grid_point, quadrant)
    ! Determing the quadrant in which a 'projected point' is situated relative to an IM grid point.
    !   quadrants:
    !   II  |   I
    !       |
    !  -----|-----
    !       |
    !  III  |  IV
    !
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_value_projected_point
    REAL(dp), INTENT(IN)  :: y_value_projected_point
    REAL(dp), INTENT(IN)  :: x_value_IM_grid_point
    REAL(dp), INTENT(IN)  :: y_value_IM_grid_point

    ! Output variables:
    INTEGER,  INTENT(OUT) :: quadrant

    IF(      x_value_projected_point >  x_value_IM_grid_point .AND. y_value_projected_point >= y_value_IM_grid_point) THEN
     ! Check for quadrant I
     quadrant = 1
    ELSE IF( x_value_projected_point <= x_value_IM_grid_point .AND. y_value_projected_point >  y_value_IM_grid_point) THEN
     ! Check for quadrant II
     quadrant = 2
    ELSE IF( x_value_projected_point <  x_value_IM_grid_point .AND. y_value_projected_point <= y_value_IM_grid_point) THEN
     ! Check for quadrant III
     quadrant = 3
    ELSE IF( x_value_projected_point >= x_value_IM_grid_point .AND. y_value_projected_point <= y_value_IM_grid_point) THEN
     ! Check for quadrant IV
     quadrant = 4
    ELSE
     STOP ' quadrant not found [find_quadrant_around_IM_grid_point, oblimap_scan_contributions_module], --STOPPED'
    END IF
  END SUBROUTINE find_quadrant_around_IM_grid_point



  SUBROUTINE find_quadrant_around_GCM_grid_point(lon_value_projected_point, lat_value_projected_point, &
                                                 lon_value_GCM_grid_point,  lat_value_GCM_grid_point, quadrant)
    ! Determing the quadrant in which a 'projected point' is situated relative to a GCM grid point.
    !   quadrants:
    !   II  |   I
    !       |
    !  -----|-----
    !       |
    !  III  |  IV
    !
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lon_value_projected_point
    REAL(dp), INTENT(IN)  :: lat_value_projected_point
    REAL(dp), INTENT(IN)  :: lon_value_GCM_grid_point
    REAL(dp), INTENT(IN)  :: lat_value_GCM_grid_point

    ! Output variables:
    INTEGER,  INTENT(OUT) :: quadrant

    IF(((lon_value_projected_point - lon_value_GCM_grid_point) <=  180._dp .AND. (lon_value_projected_point - lon_value_GCM_grid_point) > 0._dp) .OR. &
       ((lon_value_projected_point - lon_value_GCM_grid_point) <= -180._dp)) THEN
      IF(lat_value_projected_point > lat_value_GCM_grid_point) THEN
        quadrant = 1
      ELSE
        quadrant = 4
      END IF
    ELSE
      IF(lat_value_projected_point <= lat_value_GCM_grid_point) THEN
        quadrant = 3
      ELSE
        quadrant = 2
      END IF
    END IF
  END SUBROUTINE find_quadrant_around_GCM_grid_point



  REAL(dp) FUNCTION distance_over_earth_surface(lon_gcm, lat_gcm, lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points) RESULT (distance)
    ! Calculation of the distance over the Earth surface of two points which are situated at this earth surface.
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) :: lon_gcm
    REAL(dp), INTENT(IN) :: lat_gcm
    REAL(dp), INTENT(IN) :: lon_coordinates_of_im_grid_points
    REAL(dp), INTENT(IN) :: lat_coordinates_of_im_grid_points

    IF(C%choice_projection_method == 'rotation_projection') THEN
     ! In this case the distances are measured in the rectangular plane:
     distance = distance_in_flat_surface(lon_gcm, lat_gcm, lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    ELSE IF(C%vincenty_method_for_ellipsoid .AND. &
       C%choice_projection_method == 'oblique_stereographic_projection_ellipsoid_snyder' .OR. &
       C%choice_projection_method == 'oblique_lambert_equal-area_projection_ellipsoid_snyder') THEN
     distance = distance_over_curved_surface_ellipsoid(lon_gcm, lat_gcm, lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    ELSE
     distance = distance_over_curved_surface(lon_gcm, lat_gcm, lon_coordinates_of_im_grid_points, lat_coordinates_of_im_grid_points)
    END IF

    RETURN
  END FUNCTION distance_over_earth_surface



  REAL(dp) FUNCTION distance_in_flat_surface(x_coordinates_of_point_1, y_coordinates_of_point_1, x_coordinates_of_point_2, y_coordinates_of_point_2) RESULT (distance)
    ! Calculation of the distance between two points which are situated in the flat IM surface.
    USE oblimap_configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) :: x_coordinates_of_point_1
    REAL(dp), INTENT(IN) :: y_coordinates_of_point_1
    REAL(dp), INTENT(IN) :: x_coordinates_of_point_2
    REAL(dp), INTENT(IN) :: y_coordinates_of_point_2

    distance = SQRT( (x_coordinates_of_point_1 - x_coordinates_of_point_2)**2 + (y_coordinates_of_point_1 - y_coordinates_of_point_2)**2 )

    RETURN
  END FUNCTION distance_in_flat_surface



  REAL(dp) FUNCTION distance_over_curved_surface(point_1_lon, point_1_lat, point_2_lon, point_2_lat) RESULT (distance)
    ! Calculation of the distance over the Earth surface of two points which are situated at this earth surface.
    ! Note that in this function the earth surface is considered as the surface of a sphere. Currently also used for the
    ! ellipsoid projection.
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) :: point_1_lon
    REAL(dp), INTENT(IN) :: point_1_lat
    REAL(dp), INTENT(IN) :: point_2_lon
    REAL(dp), INTENT(IN) :: point_2_lat

    ! See equation (2.18) in Reerink et al. (2010):
    distance = C%earth_radius * ACOS(COS(C%degrees_to_radians * point_1_lat) * COS(C%degrees_to_radians * point_2_lat) * COS(C%degrees_to_radians * (point_1_lon - point_2_lon)) + &
                                     SIN(C%degrees_to_radians * point_1_lat) * SIN(C%degrees_to_radians * point_2_lat))

    RETURN
  END FUNCTION distance_over_curved_surface



  REAL(dp) FUNCTION distance_over_curved_surface_ellipsoid(point_1_lon, point_1_lat, point_2_lon, point_2_lat) RESULT (distance)
    ! Vincenty Inverse Solution of Geodesics on the WGS84 Ellipsoid
    !
    ! from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the
    !       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975
    !       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
    !       http://www.movable-type.co.uk/scripts/latlong-vincenty.html Chris Veness 2002-2012
    !
    ! Calculates geodetic distance between two points specified by latitude/longitude using
    ! Vincenty inverse formula for ellipsoids
    !
    ! The output variable 'distance' given in meters is:
    !  the distance between the two given points measured over the WGS84 ellipsoid
    ! or 
    !  the length of the geodesic on the WGS84 ellipsoid between the two given points
    !
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN) :: point_1_lat                                  ! Latitude  of first  point [degrees]
    REAL(dp), INTENT(IN) :: point_1_lon                                  ! Longitude of first  point [degrees]
    REAL(dp), INTENT(IN) :: point_2_lat                                  ! Latitude  of second point [degrees]
    REAL(dp), INTENT(IN) :: point_2_lon                                  ! Longitude of second point [degrees]

    ! Local variables:
    REAL(dp)             :: a                                            ! Semi-major axis of the WGS84 ellipsoid
    REAL(dp)             :: b                                            ! Semi-minor axis of the WGS84 ellipsoid
    REAL(dp)             :: f                                            ! Flattening      of the WGS84 ellipsoid
    REAL(dp)             :: L                                            ! First approximation of lambda
    REAL(dp)             :: U1                                           ! Reduced latitude of point 1
    REAL(dp)             :: U2                                           ! Reduced latitude of point 2
    REAL(dp)             :: sin_U1
    REAL(dp)             :: cos_U1
    REAL(dp)             :: sin_U2
    REAL(dp)             :: cos_U2
    REAL(dp)             :: lambda                                       ! Difference in longitude on an auxiliary sphere, positive east
    REAL(dp)             :: lambda_P
    REAL(dp)             :: sin_lambda
    REAL(dp)             :: cos_lambda
    REAL(dp)             :: sin_sigma
    REAL(dp)             :: cos_sigma          = 0._dp                   ! To prevent 'uninitialized' compiler warnings
    REAL(dp)             :: sigma              = 0._dp                   ! To prevent 'uninitialized' compiler warnings
    REAL(dp)             :: sin_alpha
    REAL(dp)             :: cos_alpha_squared  = 0._dp                   ! To prevent 'uninitialized' compiler warnings
    REAL(dp)             :: cos_2sigma_m       = 0._dp                   ! To prevent 'uninitialized' compiler warnings
    REAL(dp)             :: C_capital
    REAL(dp)             :: delta_sigma
    REAL(dp)             :: s                                            ! The length of the geodesic on the WGS84 ellipsoid between the two given points
    REAL(dp)             :: u_squared
    REAL(dp)             :: A_capital
    REAL(dp)             :: B_capital
    INTEGER              :: counter                                      ! Iteration counter

    ! WGS84 ellipsoid parameters (see WGS84 description in the oblimap_configuration_module):
    a     = C%a                                                          ! Semi-major axis of the WGS84 ellipsoid: a = 6378137           meter
    f     = C%ellipsoid_flattening                                       ! Flattening of the ellipsoid, f = 1 - (1 - e**2)**0.5 in Snyder (1987) at p. 13, f = (a-b)/a = 1 / 298.257223563, WGS84 value for f = 0.0033528106647474805
    b     = C%ellipsoid_semi_minor_axis                                  ! The semi-minor axis or the polar radius of the ellipsoid (in case of the Earth), b in Snyder (1987) at p. 160; WGS84 value for b = 6356752.314245179 meter
    L     = (point_2_lon - point_1_lon) * C%degrees_to_radians           ! Difference in longitude on an auxiliary sphere, positive east
    U1    = ATAN((1._dp - f) * TAN(point_1_lat * C%degrees_to_radians))  ! Reduced latitude, defined by tan(U) = tan(1-f) * tan(phi) 
    U2    = ATAN((1._dp - f) * TAN(point_2_lat * C%degrees_to_radians))  !  with phi the geodetic latitude, positive north of the equator
    sin_U1 = SIN(U1)
    cos_U1 = COS(U1)
    sin_U2 = SIN(U2)
    cos_U2 = COS(U2)

    ! Change the difference in longitude a tiny bit to overcome numerical problems:
    IF(ABS(point_1_lon - point_2_lon) <= 1.0E-5_dp) L = 1.0E-5_dp * C%degrees_to_radians

    lambda = L
    DO counter = 1, 101
      sin_lambda = SIN(lambda)
      cos_lambda = COS(lambda)
      sin_sigma  = SQRT((cos_U2 * sin_lambda) * (cos_U2 * sin_lambda) + (cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda) * (cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda))
      IF(sin_sigma == 0._dp) THEN
       ! Co-incident points:
       EXIT
      END IF
      cos_sigma         = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda
      sigma             = ATAN2(sin_sigma, cos_sigma)
      sin_alpha         = cos_U1 * cos_U2 * sin_lambda / sin_sigma
      cos_alpha_squared = 1._dp - sin_alpha * sin_alpha
      cos_2sigma_m      = cos_sigma - 2._dp * sin_U1 * sin_U2 / cos_alpha_squared
      IF(cos_alpha_squared == 0._dp) THEN
       ! Equatorial line: cos_alpha_squared = 0:
       cos_2sigma_m = 0._dp
      END IF
      C_capital = f / 16._dp * cos_alpha_squared * (4._dp + f * (4._dp - 3._dp * cos_alpha_squared))
      lambda_P  = lambda
      lambda    = L + (1._dp - C_capital) * f * sin_alpha * (sigma + C_capital * sin_sigma * (cos_2sigma_m + C_capital * cos_sigma*(-1._dp + 2._dp * cos_2sigma_m * cos_2sigma_m)))

      IF(ABS(lambda - lambda_P) > 1.0E-12_dp) EXIT
    END DO

    ! The formula failed to converge:
    IF(counter >= 100) THEN
     WRITE(UNIT=*, FMT='(A, 2(2F12.5, A))') ' No converge in Vincenty method for ellipsoid after 100 iterations, for (', point_1_lon, point_1_lat,') and (', point_2_lon, point_2_lat, ')'
    END IF

    IF(sin_sigma == 0._dp) THEN
     ! Co-incident points:
     distance = 0._dp
    ELSE
     u_squared   = cos_alpha_squared * (a * a - b * b) / (b * b)
     A_capital   = 1._dp + u_squared / 16384._dp * (4096._dp + u_squared * (-768._dp + u_squared * (320._dp - 175._dp * u_squared)))
     B_capital   = u_squared / 1024._dp * (256._dp + u_squared * (-128._dp + u_squared * (74._dp - 47._dp * u_squared)))
     delta_sigma = B_capital * sin_sigma * (cos_2sigma_m + B_capital / 4._dp * (cos_sigma * (-1._dp + 2._dp * cos_2sigma_m * cos_2sigma_m) - B_capital / 6._dp * cos_2sigma_m*(-3._dp + 4._dp * sin_sigma*sin_sigma)*(-3._dp + 4._dp * cos_2sigma_m * cos_2sigma_m)))
     s           = b * A_capital * (sigma - delta_sigma)

     distance = s
    END IF

    ! The initial and final bearing are not used in OBLIMAP, they could be calculated with:
    ! initial_bearing = C%radians_to_degrees * ATAN2(cos_U2 * sin_lambda,   cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda)
    ! final_bearing   = C%radians_to_degrees * ATAN2(cos_U1 * sin_lambda, - sin_U1 * cos_U2 + cos_U1 * sin_U2 * cos_lambda)

    RETURN
  END FUNCTION distance_over_curved_surface_ellipsoid



  REAL(dp) FUNCTION distance_in_im_plane_between_two_gcm_points(lon_gcm_point_1, lat_gcm_point_1, lon_gcm_point_2, lat_gcm_point_2) RESULT (distance)
    ! This routine projects two GCM points on the requested plane S' which coincides with the IM grid,
    ! with an oblique projection, and calculates the distance between them in the IM plane.
    USE oblimap_configuration_module, ONLY: dp, C
    USE oblimap_projection_module, ONLY: oblique_sg_projection, oblique_sg_projection_ellipsoid_snyder, &
      oblique_laea_projection_snyder, oblique_laea_projection_ellipsoid_snyder, rotation_projection
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lon_gcm_point_1         ! longitude coordinate (degrees) of GCM grid point 1
    REAL(dp), INTENT(IN)  :: lat_gcm_point_1         ! latitude  coordinate (degrees) of GCM grid point 1
    REAL(dp), INTENT(IN)  :: lon_gcm_point_2         ! longitude coordinate (degrees) of GCM grid point 2
    REAL(dp), INTENT(IN)  :: lat_gcm_point_2         ! latitude  coordinate (degrees) of GCM grid point 2

    ! Local variables:
    REAL(dp)              :: x_coordinate_point_1    ! The x-coordinate of GCM grid point 1 projected on S'
    REAL(dp)              :: y_coordinate_point_1    ! The y-coordinate of GCM grid point 1 projected on S'
    REAL(dp)              :: x_coordinate_point_2    ! The x-coordinate of GCM grid point 2 projected on S'
    REAL(dp)              :: y_coordinate_point_2    ! The y-coordinate of GCM grid point 2 projected on S'

    ! Determine the x,y coordinates of each GCM longitude-latitude coordinate after
    ! The oblique stereographic projection on the projection plane S'
    ! Output: x_coordinates_of_gcm_grid_points(i,j), y_coordinates_of_gcm_grid_points(i,j)
    SELECT CASE(C%choice_projection_method)
    CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder')
     CALL oblique_sg_projection(                                  lon_gcm_point_1, lat_gcm_point_1, x_coordinate_point_1, y_coordinate_point_1)
     CALL oblique_sg_projection(                                  lon_gcm_point_2, lat_gcm_point_2, x_coordinate_point_2, y_coordinate_point_2)
    CASE('oblique_stereographic_projection_ellipsoid_snyder')
     CALL oblique_sg_projection_ellipsoid_snyder(                 lon_gcm_point_1, lat_gcm_point_1, x_coordinate_point_1, y_coordinate_point_1)
     CALL oblique_sg_projection_ellipsoid_snyder(                 lon_gcm_point_2, lat_gcm_point_2, x_coordinate_point_2, y_coordinate_point_2)
    CASE('oblique_lambert_equal-area_projection_snyder')
     CALL oblique_laea_projection_snyder(                         lon_gcm_point_1, lat_gcm_point_1, x_coordinate_point_1, y_coordinate_point_1)
     CALL oblique_laea_projection_snyder(                         lon_gcm_point_2, lat_gcm_point_2, x_coordinate_point_2, y_coordinate_point_2)
    CASE('oblique_lambert_equal-area_projection_ellipsoid_snyder')
     CALL oblique_laea_projection_ellipsoid_snyder(               lon_gcm_point_1, lat_gcm_point_1, x_coordinate_point_1, y_coordinate_point_1)
     CALL oblique_laea_projection_ellipsoid_snyder(               lon_gcm_point_2, lat_gcm_point_2, x_coordinate_point_2, y_coordinate_point_2)
    CASE('rotation_projection')
     CALL rotation_projection(                                    lon_gcm_point_1, lat_gcm_point_1, x_coordinate_point_1, y_coordinate_point_1)
     CALL rotation_projection(                                    lon_gcm_point_2, lat_gcm_point_2, x_coordinate_point_2, y_coordinate_point_2)
    END SELECT

    ! The distances between the two GCM points measured after their projection in the rectangular IM plane:
    distance = distance_in_flat_surface(x_coordinate_point_1, y_coordinate_point_1, x_coordinate_point_2, y_coordinate_point_2)

    RETURN
  END FUNCTION distance_in_im_plane_between_two_gcm_points



  SUBROUTINE check_for_GCM_points_at_the_point_of_projection(lon_gcm, lat_gcm)
    ! This routine checks if a GCM point coincides with the point of projection itself:
    ! point C in Fig. 2a of Reerink et al. (2010). In fact that is a rare situation,
    ! because it is not likely to project points near point C, but in the mode of scanning
    ! all GCM points (if C%scan_search_block_size doesn't reduce the search area) this can
    ! occur. In that case it is desirable that zero division is omitted. In such case the
    ! GCM lat lon coordinates of the considered point are minimal shifted. Of course this
    ! minimal shifted point lies very very close to C, and such points are projected almost
    ! at an infinity large distance in the plane of projection, and will thus never be
    ! selected.
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! In/Output variables:
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(INOUT) :: lon_gcm   ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(INOUT) :: lat_gcm   ! The latitude  coordinates (degrees) of the GCM grid points

    ! Local variables:
    REAL(dp)                                          :: lambda_C
    REAL(dp)                                          :: phi_C
    INTEGER                                           :: i, j
    REAL(dp), PARAMETER                               :: tiny_shift = 1.0E-14_dp    ! in degrees

    ! Calculation of the coordinates of anti-pole C:
    phi_C    = - C%phi_M
    IF(C%lambda_M >= 0._dp .AND. C%lambda_M < C%Pi) THEN
     lambda_C = C%lambda_M + C%Pi
    ELSE
     lambda_C = C%lambda_M - C%Pi
    END IF
    ! Convert from radians to degree for comparison with the GCM input coordinates:
    phi_C    = phi_C    * C%radians_to_degrees
    lambda_C = lambda_C * C%radians_to_degrees

    ! Checking if any points coincide with C, and shifting these coordinates if that is the case:
    DO i = 1, C%NLON
    DO j = 1, C%NLAT
      IF(lat_gcm(i,j) == phi_C) THEN
       ! For the north and south pole any lambda is allowed because the longitude coordinate at the poles is often poorly defined:
       IF(lon_gcm(i,j) == lambda_C .OR. phi_C == 90._dp .OR. phi_C == -90._dp) THEN

        ! A little arbitrary split up, main thing is to prevent final angles being out of the -90 to +90 degree range:
        IF(phi_C > 0._dp) THEN
         lat_gcm(i,j) = phi_C - tiny_shift
        ELSE
         lat_gcm(i,j) = phi_C + tiny_shift
        END IF

        ! A little arbitrary split up, main thing is to prevent final angles being out of the 0 to 360 degree range:
        IF(lambda_C > 180._dp) THEN
         lon_gcm(i,j) = lambda_C - tiny_shift
        ELSE
         lon_gcm(i,j) = lambda_C + tiny_shift
        END IF

       END IF
      END IF
    END DO
    END DO
  END SUBROUTINE check_for_GCM_points_at_the_point_of_projection



  SUBROUTINE determining_scan_parameters(check_direction, lon_gcm, lat_gcm, advised_scan_parameter)
    ! This routine determines the scan parameters and informs the user over the scan parameters as chosen in the config file.
    USE oblimap_configuration_module, ONLY: C, dp, oblimap_scan_parameter_type, rounding
    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*)                  , INTENT(IN)           :: check_direction                                     ! This string distinguishes in the map direction
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lon_gcm                                             ! The longitude coordinates (degrees) of the GCM grid points
    REAL(dp), DIMENSION(C%NLON,C%NLAT), INTENT(IN)           :: lat_gcm                                             ! The latitude  coordinates (degrees) of the GCM grid points

    ! Output variables:
    TYPE(oblimap_scan_parameter_type) , INTENT(OUT)          :: advised_scan_parameter

    ! Local variables:
    INTEGER                                                  :: i, j
    INTEGER                                                  :: p                     =  0
    INTEGER                                                  :: zero_distance_counter =  0                          ! Counter which counts zero distances
    INTEGER                                                  :: counter               =  0                          ! Counter which counts all no-zero distances occurences
    REAL(dp)                                                 :: estimated_distance    =  0._dp
    REAL(dp)                                                 :: average_grid_distance =  0._dp
    REAL(dp)                                                 :: minimum_grid_distance =  1.0E30_dp
    REAL(dp)                                                 :: maximum_grid_distance = -1.0_dp
    REAL(dp)                                                 :: grid_distance_in_lon_direction
    REAL(dp)                                                 :: grid_distance_in_lat_direction
    REAL(dp)                                                 :: decisive_ratio
    LOGICAL                                                  :: selected_choice_quadrant_method                     ! Selecting between the advised choice_quadrant_method or the specified choice_quadrant_method in the config file

    ! The grid distance in the center of the GCM grid is estimated. This is an average of a few
    ! central neighbour points situated on the earth surface.
    SELECT CASE(C%choice_projection_method)
    CASE('rotation_projection')
     estimated_distance = MIN(ABS(lon_gcm(2,1) - lon_gcm(1,1)), ABS(lat_gcm(1,2) - lat_gcm(1,1)))
    CASE DEFAULT
     DO i = MAX(1, (C%NLON / 2) - 3), MIN(C%NLON - 1, (C%NLON / 2) + 2)
     DO j = MAX(1, (C%NLAT / 2) - 3), MIN(C%NLAT - 1, (C%NLAT / 2) + 2)
       p = p + 1
       estimated_distance = estimated_distance + distance_over_earth_surface(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i+1,j), lat_gcm(i+1,j))
     END DO
     END DO
     estimated_distance = estimated_distance / p

     IF(C%oblimap_message_level > 0) THEN
      ! Searching the minimum, maximum and average GCM grid resolution:
      DO i = 2, C%NLON-1
      DO j = 2, C%NLAT-1
        grid_distance_in_lon_direction = distance_over_earth_surface(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i+1,j  ), lat_gcm(i+1,j  ))
        grid_distance_in_lat_direction = distance_over_earth_surface(lon_gcm(i,j), lat_gcm(i,j), lon_gcm(i  ,j+1), lat_gcm(i  ,j+1))
        IF(grid_distance_in_lon_direction == 0._dp) THEN
         IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(2A, 2I4)') TRIM(C%OBLIMAP_WARNING), ' Grid distance equals zero in longitude direction for coordinates: ', i, j
         zero_distance_counter = zero_distance_counter + 1
        ELSE
         minimum_grid_distance = MIN(minimum_grid_distance, grid_distance_in_lon_direction)
         average_grid_distance = average_grid_distance + grid_distance_in_lon_direction
         counter = counter + 1
        END IF
        IF(grid_distance_in_lat_direction == 0._dp) THEN
         IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(2A, 2I4)') TRIM(C%OBLIMAP_WARNING), ' Grid distance equals zero in latitude direction for coordinates: ', i, j
         zero_distance_counter = zero_distance_counter + 1
        ELSE
         minimum_grid_distance = MIN(minimum_grid_distance, grid_distance_in_lat_direction)
         average_grid_distance = average_grid_distance + grid_distance_in_lat_direction
         counter = counter + 1
        END IF
        maximum_grid_distance = MAX(maximum_grid_distance, grid_distance_in_lon_direction)
        maximum_grid_distance = MAX(maximum_grid_distance, grid_distance_in_lat_direction)
      END DO
      END DO
      average_grid_distance = average_grid_distance / REAL(counter, dp)
      WRITE(UNIT=*, FMT='(4(A, E24.16/), A, I10)')                     &
       ' Minimum grid resolution:       '     , minimum_grid_distance, &
       ' Maximum grid resolution:       '     , maximum_grid_distance, &
       ' Average grid resolution:       '     , average_grid_distance, &
       ' Grid resolution at grid center:'     , estimated_distance   , &
       ' Number of points with zero distance:', zero_distance_counter
     END IF
    END SELECT


    ! In case the departing grid is a GCM grid, it is determined if this grid is cyclic (= periodical) in the longitude direction:
    advised_scan_parameter%data_set_is_cyclic_in_longitude = .FALSE.
    IF(check_direction == 'gcm-to-im' .AND. C%choice_projection_method /= 'rotation_projection') THEN
     IF(C%NLAT > 1) THEN
      ! Check if the grid is cyclic in logitude direction at C%NLAT / 2 which is halfway the grid in the latitude direction:
      IF(2._dp * lon_gcm(C%NLON,C%NLAT / 2) - lon_gcm(C%NLON-1,C%NLAT / 2) >= 360._dp  .AND. &
         2._dp * lon_gcm(     1,C%NLAT / 2) - lon_gcm(       2,C%NLAT / 2) <=   0._dp)       &
       advised_scan_parameter%data_set_is_cyclic_in_longitude = .TRUE.
     ELSE
      IF(2._dp * lon_gcm(C%NLON,         1) - lon_gcm(C%NLON-1,         1) >= 360._dp  .AND. &
         2._dp * lon_gcm(     1,         1) - lon_gcm(       2,         1) <=   0._dp)       &
       advised_scan_parameter%data_set_is_cyclic_in_longitude = .TRUE.
     END IF
    END IF

    IF(C%level_of_automatic_oblimap_scanning < 1) THEN
     IF(.NOT. C%full_scanning_mode) THEN
      IF(C%data_set_is_cyclic_in_longitude .NEQV. advised_scan_parameter%data_set_is_cyclic_in_longitude) WRITE(UNIT=*, FMT='(/2A, L, A, L, 2A/, A/)') &
       TRIM(C%OBLIMAP_WARNING), ' data_set_is_cyclic_in_longitude_config = ', C%data_set_is_cyclic_in_longitude, ' is advised to be ', advised_scan_parameter%data_set_is_cyclic_in_longitude, ' in ', TRIM(C%config_filename), ' Reason of warning: OBLIMAP detected a grid with longitude coordinates which cover the entire 0-360 degrees longitude range.'
     END IF
     IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(/2(A, L)/)') ' The "data_set_is_cyclic_in_longitude_config" is advised to be: ', advised_scan_parameter%data_set_is_cyclic_in_longitude, '. This run uses data_set_is_cyclic_in_longitude_config = ', C%data_set_is_cyclic_in_longitude
    END IF


    ! Determining the optimal intersection angle alpha_stereographic which determines the standard paralel in the stereographic projection:
    SELECT CASE(C%choice_projection_method)
    CASE('oblique_stereographic_projection','oblique_stereographic_projection_snyder','oblique_stereographic_projection_ellipsoid_snyder')
     IF(REAL(C%NX * C%NY * C%dx * C%dy) > 2._dp * C%pi * C%earth_radius**2) THEN
      WRITE(UNIT=*, FMT='(/A, /A, F8.3, A/)') &
       ' No optimal alpha_stereographic_config could be calcultated for this grid configuration. The best estimate might ', &
       ' be alpha_stereographic_config = 90 degrees, or take a smaller grid configuration. This run uses alpha_stereographic_config =', C%radians_to_degrees * C%alpha_stereographic, ' degrees.'
      advised_scan_parameter%alpha_stereographic        = 90._dp
     ELSE
      advised_scan_parameter%alpha_stereographic = C%radians_to_degrees * ASIN(SQRT((C%NX * C%dx * C%NY * C%dy) / (2._dp * C%pi)) / C%earth_radius)

      IF(C%oblimap_message_level > 0) WRITE(UNIT=*, FMT='(/2(A, F8.3, A/))') &
       ' An optimal alpha_stereographic_config = ', advised_scan_parameter%alpha_stereographic,   ' degrees, based on the entire IM grid area.', &
       ' Here an    alpha_stereographic_config = ', C%radians_to_degrees * C%alpha_stereographic, ' degrees is used.'
     END IF
    CASE DEFAULT
     ! In other non stereographic projections alpha is not used and default set to zero:
     advised_scan_parameter%alpha_stereographic        = 0._dp
    END SELECT
    advised_scan_parameter%alpha_stereographic = rounding(advised_scan_parameter%alpha_stereographic, 4)


    ! Determining which interpolation method is optimal: the quadrant or the radius method. In case the target grid
    ! is about four or more times coarser, then the radius method should be used (Reerink et al. 2010, p. 32 last alinea)
    SELECT CASE(check_direction)
    CASE('gcm-to-im')
     decisive_ratio = MIN(C%dx, C%dy) / (4._dp * estimated_distance)
    CASE('im-to-gcm')
     decisive_ratio = estimated_distance / (4._dp * MIN(C%dx, C%dy))
    CASE DEFAULT
     STOP 'PROGRAMMER ERROR: The first argument of the determining_scan_parameters() should be either "gcm-to-im" or "im-to-gcm".'
    END SELECT

    IF(decisive_ratio > 1._dp) THEN
     advised_scan_parameter%choice_quadrant_method = .FALSE.
    ELSE
     advised_scan_parameter%choice_quadrant_method = .TRUE.
    END IF

    IF(advised_scan_parameter%choice_quadrant_method) THEN
     IF(.NOT. C%choice_quadrant_method) THEN
      ! Advising the quadrant interpolation method (which is advised in contrast to the actual used method):
      IF(C%level_of_automatic_oblimap_scanning < 2) &
       WRITE(UNIT=*, FMT='(2A/, 2A/)') C%OBLIMAP_ADVICE, ' switch from interpolation method due to the detected ratio of the used grids: use the quadrant method for interpolation:', ' "choice_quadrant_method_config = .TRUE." in your config file: ', TRIM(C%config_filename)
     ELSE IF(C%oblimap_message_level > 0) THEN
      WRITE(UNIT=*, FMT='(A/                 )') ' OBLIMAP advices to use the quadrant method for interpolation:   choice_quadrant_method_config = .TRUE.'
     END IF
    ELSE
     IF(C%choice_quadrant_method) THEN
      ! Advising the radius interpolation method (which is advised in contrast to the actual used method):
      ! Note a first step equals the quadrant method but an additional second step is required here:
      IF(C%level_of_automatic_oblimap_scanning < 2) &
       WRITE(UNIT=*, FMT='(2A/, 2A/)') C%OBLIMAP_ADVICE, ' switch from interpolation method due to the detected ratio of the used grids: use the radius method for interpolation:',  ' "choice_quadrant_method_config = .FALSE." in your config file: ', TRIM(C%config_filename)
     ELSE IF(C%oblimap_message_level > 0) THEN
      WRITE(UNIT=*, FMT='(A/                 )') ' OBLIMAP advices to use the radius method for interpolation:   choice_quadrant_method_config = .FALSE.'
     END IF
    END IF

    SELECT CASE(C%choice_projection_method)
    CASE('rotation_projection')
     IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, F16.6, A/, 2(A, F16.1, A/))') ' The decisive ratio           =', decisive_ratio, ' which distinguises between the quadrant or the radius interpolation method.', ' The       IM grid resolution =', estimated_distance, ' meter', ' The local IM grid resolution =', MIN(C%dx, C%dy), ' meter.'
    CASE DEFAULT
     IF(C%oblimap_message_level > 1) WRITE(UNIT=*, FMT='(A, F16.6, A/, 2(A, F16.1, A/))') ' The decisive ratio      =', decisive_ratio, ' which distinguises between the quadrant or the radius interpolation method.', ' The  IM grid resolution =', MIN(C%dx, C%dy),    ' meter', ' The GCM grid resolution =', estimated_distance, ' meter.'
    END SELECT

    IF(C%level_of_automatic_oblimap_scanning >= 2) THEN
     selected_choice_quadrant_method = advised_scan_parameter%choice_quadrant_method
    ELSE
     selected_choice_quadrant_method = C%choice_quadrant_method
    END IF


    ! Providing an error message in case an invalid scan_search_block_size_config is specified:
    IF(C%scan_search_block_size < -3) THEN
     WRITE(UNIT=*, FMT='(/2A, I5, 2A/)') C%OBLIMAP_ERROR, ' The scan_search_block_size_config = ', C%scan_search_block_size, ' while it should be an integer of -3 or higher,  in ', TRIM(C%config_filename)
     STOP
    END IF


    ! In case the radius interpolation method is used, a best estimate is given for the search radius:
    SELECT CASE(check_direction)
    CASE('gcm-to-im')
     IF(.NOT. selected_choice_quadrant_method .AND. advised_scan_parameter%choice_quadrant_method .AND. C%R_search_interpolation < 0.6_dp * estimated_distance) THEN
      WRITE(UNIT=*, FMT='(2A, I12, A/)') C%OBLIMAP_ADVICE, ' If on pupose the advice to use the quadrant interpolation method is negclected, the search radius should be of the order: R_search_interpolation_config =', INT(estimated_distance), ' meter to avoid missing vales.'
      advised_scan_parameter%R_search_interpolation = estimated_distance
     ELSE
      advised_scan_parameter%R_search_interpolation = 0.4_dp * MIN(C%dx, C%dy)
     END IF
    CASE('im-to-gcm')
     advised_scan_parameter%R_search_interpolation = 0.4_dp * estimated_distance
    END SELECT
    ! The advised_scan_parameter%R_search_interpolation is rounded to a significance of 3 digits or to meters for low values:
    advised_scan_parameter%R_search_interpolation = rounding(advised_scan_parameter%R_search_interpolation, MIN(-(INT(LOG10(advised_scan_parameter%R_search_interpolation)) - (3 - 1)), 0))

    IF((.NOT. selected_choice_quadrant_method) .AND. (C%oblimap_message_level > 0)) &
     WRITE(UNIT=*, FMT='(2(A, F16.3, A/))') &
      ' A good estimate of R_search_interpolation_config =', advised_scan_parameter%R_search_interpolation, ' meter, based on 0.8 times half the grid size of the target grid.', &
      ' Here an            R_search_interpolation_config =', C%R_search_interpolation, ' meter is used.'


    IF(C%level_of_automatic_oblimap_scanning >= 0 .AND. C%level_of_automatic_oblimap_scanning <= 4) THEN
     ! Here four values of the C% struct are overwritten. Overwriting elements of the C% struct after their initialization is
     ! actually against the policy of using the C% struct in the OBLIMAP code. This is the only place violating this policy.
     IF((C%level_of_automatic_oblimap_scanning >= 1)                                       ) C%data_set_is_cyclic_in_longitude = advised_scan_parameter%data_set_is_cyclic_in_longitude
     IF((C%level_of_automatic_oblimap_scanning >= 2)                                       ) C%choice_quadrant_method          = advised_scan_parameter%choice_quadrant_method
     IF((C%level_of_automatic_oblimap_scanning >= 3) .AND. (.NOT. C%choice_quadrant_method)) C%R_search_interpolation          = advised_scan_parameter%R_search_interpolation
     IF((C%level_of_automatic_oblimap_scanning >= 4)                                       ) C%alpha_stereographic             = advised_scan_parameter%alpha_stereographic * C%degrees_to_radians
     ! The update of C%alpha_stereographic requires the reinitialization of C%akm:
     IF((C%level_of_automatic_oblimap_scanning >= 4)                                       ) C%akm                             = (1.0_dp + COS(C%alpha_stereographic)) * C%am
     IF(C%oblimap_message_level > 1) THEN
      WRITE(UNIT=*, FMT='(A       )') ' The following scan parameters are estimated by OBLIMAP in this run:'
      IF((C%level_of_automatic_oblimap_scanning >= 1)                                       ) WRITE(UNIT=*, FMT='( A, L        )') '   data_set_is_cyclic_in_longitude_config = ', C%data_set_is_cyclic_in_longitude 
      IF((C%level_of_automatic_oblimap_scanning >= 2)                                       ) WRITE(UNIT=*, FMT='( A, L        )') '   choice_quadrant_method_config = '         , C%choice_quadrant_method          
      IF((C%level_of_automatic_oblimap_scanning >= 3) .AND. (.NOT. C%choice_quadrant_method)) WRITE(UNIT=*, FMT='( A, E24.16   )') '   R_search_interpolation_config = '         , C%R_search_interpolation          
      IF((C%level_of_automatic_oblimap_scanning >= 4)                                       ) WRITE(UNIT=*, FMT='( A, F8.3  , A)') '   alpha_stereographic_config = '            , C%alpha_stereographic * C%radians_to_degrees, ' degrees'
      WRITE(UNIT=*, FMT='( A       )') ''
     END IF
    ELSE
     WRITE(UNIT=*, FMT='(/2A, I5, 2A/)') C%OBLIMAP_ERROR, ' Invalid value for  level_of_automatic_oblimap_scanning_config = ', C%level_of_automatic_oblimap_scanning, '  in ', TRIM(C%config_filename)
     WRITE(UNIT=*, FMT='(A )') '                The level_of_automatic_oblimap_scanning_config should have one of the following levels:'
     WRITE(UNIT=*, FMT='(A )') '                 level = 0: No scan parameter is determined automatically by OBLIMAP.'
     WRITE(UNIT=*, FMT='(A )') '                 level = 1: OBLIMAP determines the data_set_is_cyclic_in_longitude_config automatically.'
     WRITE(UNIT=*, FMT='(A )') '                 level = 2: OBLIMAP determines the choice_quadrant_method_config          automatically, further like level 1.'
     WRITE(UNIT=*, FMT='(A )') '                 level = 3: OBLIMAP determines the R_search_interpolation_config          automatically, further like level 2.'
     WRITE(UNIT=*, FMT='(A/)') '                 level = 4: OBLIMAP determines the alpha_stereographic_config             automatically, further like level 3.'
     STOP
    END IF
  END SUBROUTINE determining_scan_parameters

END MODULE oblimap_scan_contributions_module
