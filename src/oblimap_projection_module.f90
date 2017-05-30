! File name: oblimap_projection_module.f90
!
! Copyright (C) 2016 Thomas Reerink & Michael Kliphuis.
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

MODULE oblimap_projection_module

CONTAINS
  SUBROUTINE oblique_sg_projection(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate
    ! system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: t_P_prime

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%degrees_to_radians * phi
    lambda_P = C%degrees_to_radians * lambda

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(C%alpha_stereographic)) / (1._dp + COS(phi_P) * COS(C%phi_M) * COS(lambda_P - C%lambda_M) + SIN(phi_P) * SIN(C%phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  C%earth_radius * (COS(phi_P) * SIN(lambda_P - C%lambda_M)) * t_P_prime
    y_IM_P_prime =  C%earth_radius * (SIN(phi_P) * COS(C%phi_M) - (COS(phi_P) * SIN(C%phi_M)) * COS(lambda_P - C%lambda_M)) * t_P_prime

    ! See equation (21-4) on page 157 in Snyder (1987):
    IF(PRESENT(k_P)) k_P = (1._dp + COS(C%alpha_stereographic)) / (1._dp + SIN(C%phi_M) * SIN(phi_P) + COS(C%phi_M) * COS(phi_P) * COS(lambda_P - C%lambda_M))
  END SUBROUTINE oblique_sg_projection



  SUBROUTINE inverse_oblique_sg_projection(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime  ! in meter
    REAL(dp)              :: y_3D_P_prime  ! in meter
    REAL(dp)              :: z_3D_P_prime  ! in meter
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P        ! in meter
    REAL(dp)              :: y_3D_P        ! in meter
    REAL(dp)              :: z_3D_P        ! in meter

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = C%earth_radius * COS(C%alpha_stereographic) * COS(C%lambda_M) * COS(C%phi_M) - SIN(C%lambda_M) * x_IM_P_prime - COS(C%lambda_M) * SIN(C%phi_M) * y_IM_P_prime
    y_3D_P_prime = C%earth_radius * COS(C%alpha_stereographic) * SIN(C%lambda_M) * COS(C%phi_M) + COS(C%lambda_M) * x_IM_P_prime - SIN(C%lambda_M) * SIN(C%phi_M) * y_IM_P_prime
    z_3D_P_prime = C%earth_radius * COS(C%alpha_stereographic) *                   SIN(C%phi_M)                                  +                   COS(C%phi_M) * y_IM_P_prime

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(C%lambda_M) * COS(C%phi_M) * x_3D_P_prime  +  SIN(C%lambda_M) * COS(C%phi_M) * y_3D_P_prime  +  SIN(C%phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * C%earth_radius**2 + 2._dp * C%earth_radius * a) / (C%earth_radius**2 + 2._dp * C%earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  C%earth_radius * COS(C%lambda_M) * COS(C%phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  C%earth_radius * SIN(C%lambda_M) * COS(C%phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  C%earth_radius *                   SIN(C%phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
     lambda_P = 180._dp + C%radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
     lambda_P =           C%radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 360._dp + C%radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
     lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
     lambda_P =   0._dp
    END IF

   ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
   IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
    phi_P = C%radians_to_degrees * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
   ELSE IF(z_3D_P >  0._dp) THEN
    phi_P =   90._dp
   ELSE IF(z_3D_P <  0._dp) THEN
    phi_P =  -90._dp
   END IF
  END SUBROUTINE inverse_oblique_sg_projection



  SUBROUTINE inverse_oblique_sg_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection the
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C       ! in radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator

    ! See equation (20-18) on page 159 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
    ! See equation (21-15) on page 159 Snyder (1987), because the denumerator is always positive this ATAN doesn't
    ! need a correction like note 2 on page ix in Snyder (1987):
    angle_C  = 2._dp * ATAN(rho / ((1._dp + COS(C%alpha_stereographic)) * C%earth_radius))

    ! See equation (20-14) on page 158 Snyder (1987):
    phi_P    = C%radians_to_degrees * ( ASIN(COS(angle_C) * SIN(C%phi_M) + ((y_IM_P_prime * SIN(angle_C) * COS(C%phi_M)) / rho)) )

    ! See equation (20-15) on page 159 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = rho * COS(C%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%phi_M) * SIN(angle_C)
    lambda_P    = C%radians_to_degrees * (C%lambda_M + arctangens_quotient(numerator, denumerator))

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

    ! In case point P coincides with M (see condition at the first line of page  159 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%radians_to_degrees * C%lambda_M
     phi_P    = C%radians_to_degrees * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_sg_projection_snyder



  SUBROUTINE oblique_laea_projection_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection the
    ! longitude-latitude coordinates which coincide with the GCM grid points to the rectangular IM
    ! coordinate system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda        ! in degrees
    REAL(dp), INTENT(IN)  :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT) :: y_IM_P_prime  ! in meter

    ! Local variables:
    REAL(dp)              :: phi_P         ! in radians
    REAL(dp)              :: lambda_P      ! in radians
    REAL(dp)              :: t_P_prime

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%degrees_to_radians * phi
    lambda_P = C%degrees_to_radians * lambda

    ! See equation (21-4) on page 185 of Snyder (1987):
    t_P_prime = SQRT(2._dp / (1._dp + COS(phi_P) * COS(C%phi_M) * COS(lambda_P - C%lambda_M) + SIN(phi_P) * SIN(C%phi_M)))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010), page 185 of Snyder (1987):
    x_IM_P_prime =  C%earth_radius * (COS(phi_P) * SIN(lambda_P - C%lambda_M)) * t_P_prime
    y_IM_P_prime =  C%earth_radius * (SIN(phi_P) * COS(C%phi_M) - (COS(phi_P) * SIN(C%phi_M)) * COS(lambda_P - C%lambda_M)) * t_P_prime
  END SUBROUTINE oblique_laea_projection_snyder



  SUBROUTINE inverse_oblique_laea_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C       ! in radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator

    ! See equation (20-18) on page 187 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
    ! See equation (24-16) on page 187 Snyder (1987):
    angle_C  = 2._dp * ASIN(rho / (2._dp * C%earth_radius))

    ! See equation (20-14) on page 186 Snyder (1987):
    phi_P    = C%radians_to_degrees * ( ASIN(COS(angle_C) * SIN(C%phi_M) + ((y_IM_P_prime * SIN(angle_C) * COS(C%phi_M)) / rho)) )

    ! See equation (20-15) on page 186 Snyder (1987):
    numerator   = x_IM_P_prime * SIN(angle_C)
    denumerator = rho * COS(C%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%phi_M) * SIN(angle_C)
    lambda_P    = C%radians_to_degrees * (C%lambda_M + arctangens_quotient(numerator, denumerator))

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%radians_to_degrees * C%lambda_M
     phi_P    = C%radians_to_degrees * C%phi_M
    END IF
  END SUBROUTINE inverse_oblique_laea_projection_snyder



  SUBROUTINE oblique_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with Snyder's oblique stereographic projection for the ellipsoid
    ! the the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians,  phi    in Snyder (1987)
    REAL(dp)                        :: lambda_P      ! in radians,  lambda in Snyder (1987)
    REAL(dp)                        :: chi_P         ! in radians,  chi    in Snyder (1987)
    REAL(dp)                        :: A

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IF(C%polar_projection) THEN
     ! The polar case is excepted from the oblique formula's, see page 161 Snyder (1987)

     IF(PRESENT(k_P)) THEN
      ! Output: x_IM_P_prime, y_IM_P_prime, k_P
      CALL polar_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
     ELSE
      ! Output: x_IM_P_prime, y_IM_P_prime
      CALL polar_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
     END IF
    ELSE
     ! The oblique case, see page 160 Snyder (1987)

     ! Convert longitude-latitude coordinates to radians:
     phi_P    = C%degrees_to_radians * phi
     lambda_P = C%degrees_to_radians * lambda

     ! See equations (3-1a) and (21-27) on page 160 in Snyder (1987):
     chi_P = 2._dp * ATAN(SQRT(((1._dp +       SIN(phi_P)) / (1._dp -       SIN(phi_P))) * &
                               ((1._dp - C%e * SIN(phi_P)) / (1._dp + C%e * SIN(phi_P)))**(C%e))) - 0.5_dp * C%pi
     A     = C%akm / (COS(C%chi_M) * (1._dp + SIN(C%chi_M) * SIN(chi_P) + COS(C%chi_M) * COS(chi_P) * COS(lambda_P - C%lambda_M)))


     ! See equations (21-24) and (21-25) on page 160 in Snyder (1987):
     x_IM_P_prime =  A * COS(chi_P) * SIN(lambda_P - C%lambda_M)
     y_IM_P_prime =  A * (COS(C%chi_M) * SIN(chi_P) - SIN(C%chi_M) * COS(chi_P) * COS(lambda_P - C%lambda_M))

     ! See equation (21-26) on page 160 in Snyder (1987):
     IF(PRESENT(k_P)) k_P = (A * COS(chi_P)) / (C%a * (COS(phi_P) / SQRT(1.0_dp - (C%e * SIN(phi_P))**2)))
    END IF
  END SUBROUTINE oblique_sg_projection_ellipsoid_snyder



  SUBROUTINE polar_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with Snyder's polar stereographic projection for the ellipsoid
    ! the the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y). See Snyder (1987) p. 161.
    !
    ! The examples of Snyder (1987) at p. 314-315 with the international ellipsoid are used to
    ! validate this forward SG projection on the ellipsoid for the polar aspect.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: phi_C         ! in radians,  phi_c  in Snyder (1987), the standard parallel
    REAL(dp)                        :: t_P           ! 
    REAL(dp)                        :: t_C           ! 
    REAL(dp)                        :: m_C           ! 
    REAL(dp)                        :: rho           ! 
    REAL(dp)                        :: pf            ! polar factor: -1.0 for SP and +1.0 for NP
    REAL(dp)                        :: k0            ! Length scale factor at center M [-]

    IF(C%phi_M == - 90.0_dp * C%degrees_to_radians) THEN
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    phi_C    = (C%degrees_to_radians * 90.0_dp - C%alpha_stereographic) * pf

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%degrees_to_radians * phi
    lambda_P = C%degrees_to_radians * lambda

    IF(C%alpha_stereographic == 0.0_dp) THEN
     ! This variant is only considered for the case that C%alpha_stereographic = 0, i.e. k0 = 1, for which the else-option does not work).
     ! The POLAR ASPECT WITH KNOWN k0 case:
    !t_P = TAN(pi / 4.0_dp - phi_P * pf / 2.0_dp) / ( (1.0_dp - C%e * SIN(phi_P * pf)) / (1.0_dp + C%e * SIN(phi_P * pf)) )**(C%e / 2.0_dp)                      ! (15-9)  on page 161 in Snyder (1987)
     t_P = ( ((1.0_dp - SIN(phi_P * pf)) / (1.0_dp + SIN(phi_P * pf))) * ((1.0_dp + C%e * SIN(phi_P * pf)) / (1.0_dp - C%e * SIN(phi_P * pf)))**C%e )**0.5_dp    ! (15-9a) on page 161 in Snyder (1987)
     k0  = 0.5_dp * (1.0_dp + COS(C%alpha_stereographic))                                                                                                        ! in fact it will be always 1 because it is only used for C%alpha_stereographic = 0
     rho = 2.0_dp * C%a * k0 * t_P / ( (1.0_dp + C%e)**(1.0_dp + C%e) * (1.0_dp - C%e)**(1.0_dp - C%e) )**0.5_dp                                                 ! (21-33) on page 161 in Snyder (1987)
    ELSE
     ! The POLAR ASPECT WITH KNOWN STANDARD PARALLEL NOT AT POLE case:
    !t_P   = TAN(pi / 4.0_dp - phi_P * pf / 2.0_dp) / ( (1.0_dp - C%e * SIN(phi_P * pf)) / (1.0_dp + C%e * SIN(phi_P * pf)) )**(C%e / 2.0_dp)                    ! (15-9)  on page 161 in Snyder (1987)
     t_P   = ( ((1.0_dp - SIN(phi_P * pf)) / (1.0_dp + SIN(phi_P * pf))) * ((1.0_dp + C%e * SIN(phi_P * pf)) / (1.0_dp - C%e * SIN(phi_P * pf)))**C%e )**0.5_dp  ! (15-9a) on page 161 in Snyder (1987)
    !t_C   = TAN(pi / 4.0_dp - phi_C * pf / 2.0_dp) / ( (1.0_dp - C%e * SIN(phi_C * pf)) / (1.0_dp + C%e * SIN(phi_C * pf)) )**(C%e / 2.0_dp)                    ! (15-9)  on page 161 in Snyder (1987)
     t_C   = ( ((1.0_dp - SIN(phi_C * pf)) / (1.0_dp + SIN(phi_C * pf))) * ((1.0_dp + C%e * SIN(phi_C * pf)) / (1.0_dp - C%e * SIN(phi_C * pf)))**C%e )**0.5_dp  ! (15-9a) on page 161 in Snyder (1987)
     m_C   = COS(phi_C * pf) / (1.0_dp - C%e**2 * SIN(phi_C * pf)**2)**0.5_dp                                                                                    ! (14-15) on page 160 in Snyder (1987)
     rho   = C%a * m_C * t_P / t_C                                                                                                                               ! (21-34) on page 161 in Snyder (1987)
    END IF

    x_IM_P_prime   =   rho * SIN(lambda_P * pf - C%lambda_M * pf) * pf                                                                                           ! (21-30) on page 161 in Snyder (1987)
    y_IM_P_prime   = - rho * COS(lambda_P * pf - C%lambda_M * pf) * pf                                                                                           ! (21-31) on page 161 in Snyder (1987)

    IF(PRESENT(k_P)) k_P = rho / (C%a * (COS(phi_P * pf) / SQRT(1.0_dp - (C%e * SIN(phi_P * pf))**2)))                                                           ! (21-32) on page 161 in Snyder (1987)
  END SUBROUTINE polar_sg_projection_ellipsoid_snyder



  SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection for the ellipsoid
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C       ! in radians
    REAL(dp)              :: chi_P         ! in radians,  chi in Snyder (1987)
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator

    IF(C%polar_projection) THEN
     ! The inverse polar case is excepted from the inverse oblique formula's, see page 161 Snyder (1987)

     ! Output: lambda_P, phi_P
     CALL inverse_polar_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ELSE
     ! The inverse oblique case, see page 160 Snyder (1987)

     ! See equation (20-18) on page 162 Snyder (1987):
     rho     = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)
     ! See equation (21-38) on page 162 Snyder (1987):
     angle_C = 2._dp * ATAN(rho * COS(C%chi_M) / C%akm)

     ! See equations (21-37) on page 161 in Snyder (1987):
     chi_P   = ASIN(COS(angle_C) * SIN(C%chi_M) + y_IM_P_prime * SIN(angle_C) * COS(C%chi_M) / rho)

     ! See equation (3-5) on page 162 instead of equation (3-4) on page 161 Snyder (1987):
     phi_P = C%radians_to_degrees * (chi_P + &
             (C%e**2 / 2._dp + 5._dp * C%e**4 / 24._dp +          C%e**6 / 12._dp  +   13._dp * C%e**8 /    360._dp) * SIN(2._dp * chi_P) + &
             (                 7._dp * C%e**4 / 48._dp + 29._dp * C%e**6 / 240._dp +  811._dp * C%e**8 /  11520._dp) * SIN(4._dp * chi_P) + &
             (                                            7._dp * C%e**6 / 120._dp +   81._dp * C%e**8 /   1120._dp) * SIN(6._dp * chi_P) + &
             (                                                                       4279._dp * C%e**8 / 161280._dp) * SIN(8._dp * chi_P))

     ! See equation (21-36) on page 161 Snyder (1987):
     numerator   = x_IM_P_prime * SIN(angle_C)
     denumerator = rho * COS(C%chi_M) * COS(angle_C) - y_IM_P_prime * SIN(C%chi_M) * SIN(angle_C)
     lambda_P    = C%radians_to_degrees * (C%lambda_M + arctangens_quotient(numerator, denumerator))

     ! Our choice is to return lambda in the 0-360 degree range:
     IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

     ! In case point P coincides with M (see condition at the first line of page  162 Snyder (1987):
     IF(rho == 0._dp) THEN
      lambda_P = C%radians_to_degrees * C%lambda_M
      phi_P    = C%radians_to_degrees * C%phi_M
     END IF

    END IF
  END SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder



  SUBROUTINE inverse_polar_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse polar stereographic projection for the ellipsoid
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees. See Snyder (1987) p. 162.
    !
    ! The examples of Snyder (1987) at p. 317-318 with the international ellipsoid are used to
    ! validate this inverse SG projection on the ellipsoid for the polar aspect.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees,  lambda in Snyder (1987)
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees,  phi    in Snyder (1987)

    ! Local variables:
    REAL(dp)              :: chi_P         ! in radians,  chi    in Snyder (1987)
    REAL(dp)              :: phi_C         ! in radians,  phi_c  in Snyder (1987), the standard parallel
    REAL(dp)              :: t_P
    REAL(dp)              :: t_C
    REAL(dp)              :: m_C
    REAL(dp)              :: rho
    REAL(dp)              :: pf            ! polar factor: -1.0 for SP and +1.0 for NP
    REAL(dp)              :: k0

    IF(C%phi_M == - 90.0_dp * C%degrees_to_radians) THEN
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    phi_C = (C%degrees_to_radians * 90.0_dp - C%alpha_stereographic) * pf

    rho   = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)                                                                                                           ! (20-18) on page 162 in Snyder (1987)

    IF(C%alpha_stereographic == 0.0_dp) THEN
     ! This variant is only considered for the case that C%alpha_stereographic = 0, i.e. k0 = 1, for which the else-option does not work).
     ! The POLAR ASPECT WITH KNOWN k0 case:
     k0  = 0.5_dp * (1.0_dp + COS(C%alpha_stereographic))                                                                                                     ! in fact it will be always 1 because it is only used for C%alpha_stereographic = 0
     t_P = rho * ( (1.0_dp + C%e)**(1.0_dp + C%e) * (1.0_dp - C%e)**(1.0_dp - C%e) )**0.5_dp / (2.0_dp * C%a * k0)                                            ! (21-39) on page 162 in Snyder (1987)
    ELSE
     ! The POLAR ASPECT WITH KNOWN STANDARD PARALLEL NOT AT POLE case:
    !t_C = TAN(pi / 4.0_dp - phi_C * pf / 2.0_dp) / ( (1.0_dp - C%e * SIN(phi_C * pf)) / (1.0_dp + C%e * SIN(phi_C * pf)) )**(C%e / 2.0_dp)                   ! (15-9)  on page 161 in Snyder (1987)
     t_C = ( ((1.0_dp - SIN(phi_C * pf)) / (1.0_dp + SIN(phi_C * pf))) * ((1.0_dp + C%e * SIN(phi_C * pf)) / (1.0_dp - C%e * SIN(phi_C * pf)))**C%e )**0.5_dp ! (15-9a) on page 161 in Snyder (1987)
     m_C = COS(phi_C * pf) / (1.0_dp - C%e**2 * SIN(phi_C * pf)**2)**0.5_dp                                                                                   ! (14-15) on page 160 in Snyder (1987)
     t_P = rho * t_C / (C%a * m_C)                                                                                                                            ! (21-40) on page 162 in Snyder (1987)
    END IF

    ! Note: Eventually replace ATAN with arctangens_quotient(numerator = rho * t_C, denumerator = C%a * m_C) :
    chi_P = C%pi / 2.0_dp - 2.0_dp * ATAN(t_P)                                                                                                                ! (7-13)  on page 162 in Snyder (1987)

    ! See equation (3-5) on page 162 instead of equation (7-9) on page 162 Snyder (1987):
    phi_P = C%radians_to_degrees * (chi_P + &
            (C%e**2 / 2._dp + 5._dp * C%e**4 / 24._dp +          C%e**6 / 12._dp  +   13._dp * C%e**8 /    360._dp) * SIN(2._dp * chi_P) + &
            (                 7._dp * C%e**4 / 48._dp + 29._dp * C%e**6 / 240._dp +  811._dp * C%e**8 /  11520._dp) * SIN(4._dp * chi_P) + &
            (                                            7._dp * C%e**6 / 120._dp +   81._dp * C%e**8 /   1120._dp) * SIN(6._dp * chi_P) + &
            (                                                                       4279._dp * C%e**8 / 161280._dp) * SIN(8._dp * chi_P)) *pf                 ! (3-5)   on page 162 in Snyder (1987)
    lambda_P = C%radians_to_degrees * (C%lambda_M * pf + arctangens_quotient(x_IM_P_prime * pf, - y_IM_P_prime * pf)) *pf                                     ! (20-16) on page 162 in Snyder (1987)

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
  END SUBROUTINE inverse_polar_sg_projection_ellipsoid_snyder



  SUBROUTINE oblique_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection for
    ! the ellipsoid the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: q_P           ! in radians,  q in Snyder (1987)
    REAL(dp)                        :: beta
    REAL(dp)                        :: B

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IF(C%polar_projection) THEN
     ! The polar case is excepted from the oblique formula's, see page 187-188 Snyder (1987)

     IF(PRESENT(k_P)) THEN
      ! Output: x_IM_P_prime, y_IM_P_prime, k_P
      CALL polar_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
     ELSE
      ! Output: x_IM_P_prime, y_IM_P_prime
      CALL polar_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime)
     END IF
    ELSE
     ! The oblique case, see page 160 Snyder (1987)

     ! Convert longitude-latitude coordinates to radians:
     phi_P    = C%degrees_to_radians * phi
     lambda_P = C%degrees_to_radians * lambda

     ! See equation (3-12) on page 187 in Snyder (1987):
     q_P = (1._dp - C%e**2) * ((SIN(phi_P) / (1._dp - (C%e * SIN(phi_P))**2)) - (1._dp / (2._dp * C%e)) * LOG((1._dp - C%e * SIN(phi_P)) / (1._dp + C%e * SIN(phi_P))))
     ! See equation (3-11) on page 187 in Snyder (1987):
     beta = ASIN(q_P / C%q_polar)
     ! See equation (24-19) on page 187 in Snyder (1987):
     B = C%R_q_polar * SQRT(2._dp / (1._dp + SIN(C%beta_M) * SIN(beta) + COS(C%beta_M) * COS(beta) * COS(lambda_P - C%lambda_M)))

     ! See equation (24-17) and (24-18) on page 187 in Snyder (1987):
     x_IM_P_prime = B * C%D * COS(beta) * SIN(lambda_P - C%lambda_M)
     y_IM_P_prime = (B / C%D) * (COS(C%beta_M) * SIN(beta) - SIN(C%beta_M) * COS(beta) * COS(lambda_P - C%lambda_M))
    END IF
  END SUBROUTINE oblique_laea_projection_ellipsoid_snyder



  SUBROUTINE polar_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with Snyder's polar Lambert azimuthal equal-area projection for
    ! the ellipsoid the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y). See Snyder (1987) p. 188.
    !
    ! The examples of Snyder (1987) at p. 334-345 with the international ellipsoid are used to
    ! validate this forward LAEA projection on the ellipsoid for the polar aspect.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-], k in Snyder

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: q_P           ! in radians,  q in Snyder (1987)
    REAL(dp)                        :: rho
    REAL(dp)                        :: pf            ! polar factor: -1.0 for SP and +1.0 for NP

    IF(C%phi_M == - 90.0_dp * C%degrees_to_radians) THEN
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = C%degrees_to_radians * phi
    lambda_P = C%degrees_to_radians * lambda

    ! See equation (3-12) on page 187 in Snyder (1987):
    q_P = (1._dp - C%e**2) * ((SIN(phi_P) / (1._dp - (C%e * SIN(phi_P))**2)) - (1._dp / (2._dp * C%e)) * LOG((1._dp - C%e * SIN(phi_P)) / (1._dp + C%e * SIN(phi_P))))
    ! See equation (24-23) and (24-25) on page 188 in Snyder (1987):
    rho = C%a * SQRT(C%q_polar - q_P * pf)

    ! See equation (21-30), (21-31) and (24-24) on page 188 in Snyder (1987):
    x_IM_P_prime =   rho * SIN(lambda_P - C%lambda_M)
    y_IM_P_prime = - rho * COS(lambda_P - C%lambda_M) * pf

    ! See equation (21-32) on page 188 in Snyder (1987):
    IF(PRESENT(k_P)) k_P = rho / (C%a * (COS(phi_P     ) / SQRT(1.0_dp - (C%e * SIN(phi_P     ))**2)))
   !IF(PRESENT(k_P)) k_P = rho / (C%a * (COS(phi_P * pf) / SQRT(1.0_dp - (C%e * SIN(phi_P * pf))**2))) ! Check if perhaps this is intended by Snyder, but makes no difference.
  END SUBROUTINE polar_laea_projection_ellipsoid_snyder



  SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection for
    ! the ellipsoid the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C       ! in radians
    REAL(dp)              :: beta          ! in radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator

    IF(C%polar_projection) THEN
     ! The inverse polar case is excepted from the inverse oblique formula's, see page 190 Snyder (1987)

     ! Output: lambda_P, phi_P
     CALL inverse_polar_laea_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ELSE
     ! The inverse oblique case, see page 188-189 Snyder (1987)

     ! See equation (24-28) on page 189 Snyder (1987):
     rho      = SQRT((x_IM_P_prime / C%D)**2 + (C%D * y_IM_P_prime)**2)
     ! See equation (24-29) on page 189 Snyder (1987):
     angle_C  = 2._dp * ASIN(rho / (2._dp * C%R_q_polar))

     ! See equation (24-30) on page 189 Snyder (1987):
     beta = ASIN(COS(angle_C) * SIN(C%beta_M) + (C%D * y_IM_P_prime * SIN(angle_C) * COS(C%beta_M) / rho))

     ! See equation (3-18) on page 189 instead of equation (3-16) on page 188 Snyder (1987):
     phi_P = C%radians_to_degrees * (beta + &
             (C%e**2 / 3._dp + 31._dp * C%e**4 / 180._dp + 517._dp * C%e**6 /  5040._dp) * SIN(2._dp * beta) + &
             (                 23._dp * C%e**4 / 360._dp + 251._dp * C%e**6 /  3780._dp) * SIN(4._dp * beta) + &
             (                                             761._dp * C%e**6 / 45360._dp) * SIN(6._dp * beta))

     ! See equation (20-26) on page 188 Snyder (1987):
     numerator   = x_IM_P_prime * SIN(angle_C)
     denumerator = C%D * rho * COS(C%beta_M) * COS(angle_C) - C%D**2 * y_IM_P_prime * SIN(C%beta_M) * SIN(angle_C)
     lambda_P    = C%radians_to_degrees * (C%lambda_M + arctangens_quotient(numerator, denumerator))

     ! Our choice is to return lambda in the 0-360 degree range:
     IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

     ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
     IF(rho == 0._dp) THEN
      lambda_P = C%radians_to_degrees * C%lambda_M
      phi_P    = C%radians_to_degrees * C%phi_M
     END IF
    END IF
  END SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder



  SUBROUTINE inverse_polar_laea_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection for
    ! the ellipsoid the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees. See Snyder (1987) p. 190.
    !
    ! The examples of Snyder (1987) at p. 336-337 with the international ellipsoid are used to
    ! validate this inverse SG projection on the ellipsoid for the polar aspect.

    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: beta          ! in radians
    REAL(dp)              :: pf            ! polar factor: -1.0 for SP and +1.0 for NP

    IF(C%phi_M == - 90.0_dp * C%degrees_to_radians) THEN
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    ! See equation (20-18) on page 190 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)

    ! See equation (24-32) on page 190 Snyder (1987):
    beta = pf * ASIN(1.0_dp - rho**2 / (C%a**2 * ( 1.0_dp - ((1.0_dp - C%e**2) / (2.0_dp * C%e)) * LOG((1.0_dp - C%e) / (1.0_dp + C%e)) )))

    ! See equation (3-18) on page 189 instead of equation (3-16) on page 188 as described on page 190 Snyder (1987):
    phi_P = C%radians_to_degrees * (beta + &
            (C%e**2 / 3._dp + 31._dp * C%e**4 / 180._dp + 517._dp * C%e**6 /  5040._dp) * SIN(2._dp * beta) + &
            (                 23._dp * C%e**4 / 360._dp + 251._dp * C%e**6 /  3780._dp) * SIN(4._dp * beta) + &
            (                                             761._dp * C%e**6 / 45360._dp) * SIN(6._dp * beta))

    ! See equation (20-16) on page 190 Snyder (1987):
    lambda_P    = C%radians_to_degrees * (C%lambda_M + arctangens_quotient(x_IM_P_prime, - y_IM_P_prime *pf))

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = C%radians_to_degrees * C%lambda_M
     phi_P    = C%radians_to_degrees * C%phi_M
    END IF
  END SUBROUTINE inverse_polar_laea_projection_ellipsoid_snyder



  FUNCTION arctangens_quotient(numerator, denumerator) RESULT(angle)
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: numerator
    REAL(dp), INTENT(IN)  :: denumerator

    ! Result variables:
    REAL(dp)              :: angle       ! in radians

    ! Local variables:
    REAL(dp)              :: quadrant_correction

    ! See note 2 on page ix in Snyder (1987), to distinguish between the quadrants:
    quadrant_correction = 0._dp
    IF(denumerator <  0._dp) quadrant_correction =          C%pi
    IF(denumerator == 0._dp) quadrant_correction = 0.5_dp * C%pi
    IF(numerator   <  0._dp) quadrant_correction = - quadrant_correction

    IF(denumerator == 0._dp) THEN
     IF(numerator == 0._dp) THEN
      ! The angle is indetermined, usually zero is taken:
      angle = 0._dp
     ELSE
      angle = 0.5_dp * C%pi * (numerator / ABS(numerator))
     END IF
    ELSE
     angle = ATAN(numerator / denumerator) + quadrant_correction
    END IF
  END FUNCTION arctangens_quotient



  SUBROUTINE rotation_projection(x_IM, y_IM, x_IM_prime, y_IM_prime)
    ! This subroutine transforms the 2D coordinates which are projected with a rotation from
    ! a rectangular IM coordinate system (x,y) to another rectangular IM prime coordinate
    ! system (x',y')
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM        ! in meter
    REAL(dp), INTENT(IN)  :: y_IM        ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_prime  ! in meter
    REAL(dp), INTENT(OUT) :: y_IM_prime  ! in meter

    ! See for a derivation: http://www.youtube.com/watch?v=h11ljFJeaLo
    x_IM_prime =   x_IM * COS(C%theta_rotation_projection) + y_IM * SIN(C%theta_rotation_projection)
    y_IM_prime = - x_IM * SIN(C%theta_rotation_projection) + y_IM * COS(C%theta_rotation_projection)
  END SUBROUTINE rotation_projection



  SUBROUTINE inverse_rotation_projection(x_IM_prime, y_IM_prime, x_IM, y_IM)
    ! This subroutine transforms the 2D coordinates which are projected with a inverse rotation
    ! from a rectangular IM prime coordinate system (x',y') to another rectangular IM coordinate
    ! system (x,y)
    USE oblimap_configuration_module, ONLY: dp, C
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM        ! in meter
    REAL(dp), INTENT(OUT) :: y_IM        ! in meter

    ! This equations are derived by taking a linear combination of the rotation projection equations.
    x_IM = x_IM_prime * COS(C%theta_rotation_projection) - y_IM_prime * SIN(C%theta_rotation_projection)
    y_IM = x_IM_prime * SIN(C%theta_rotation_projection) + y_IM_prime * COS(C%theta_rotation_projection)
  END SUBROUTINE inverse_rotation_projection

END MODULE oblimap_projection_module
