! ############################################################
! #                                                          #
! #             THE TWOSTREAM LIDORT MODEL                   #
! #                                                          #
! #      (LInearized Discrete Ordinate Radiative Transfer)   #
! #       --         -        -        -         -           #
! #                                                          #
! ############################################################

! ############################################################
! #                                                          #
! #  Authors :      Robert. J. D. Spurr (1)                  #
! #                 Vijay Natraj        (2)                  #
! #                                                          #
! #  Address (1) :     RT Solutions, Inc.                    #
! #                    9 Channing Street                     #
! #                    Cambridge, MA 02138, USA              #
! #  Tel:             (617) 492 1183                         #
! #  Email :           rtsolutions@verizon.net               #
! #                                                          #
! #  Address (2) :     CalTech                               #
! #                    Department of Planetary Sciences      #
! #                    1200 East California Boulevard        #
! #                    Pasadena, CA 91125                    #
! #  Tel:             (626) 395 6962                         #
! #  Email :           vijay@gps.caltech.edu                 #
! #                                                          #
! #  Version 1.0-1.3 :                                       #
! #     Mark 1: October  2010                                #
! #     Mark 2: May      2011, with BRDFs                    #
! #     Mark 3: October  2011, with Thermal sources          #
! #                                                          #
! #  Version 2.0-2.1 :                                       #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays   #
! #     Mark 5: December 2012, Observation Geometry option   #
! #                                                          #
! #  Version 2.2-2.3 :                                       #
! #     Mark 6: July     2013, Level outputs + control       #
! #     Mark 7: December 2013, Flux outputs  + control       #
! #     Mark 8: January  2014, Surface Leaving + control     #
! #     Mark 9: June     2014, Inverse Pentadiagonal         #
! #                                                          #
! #  Version 2.4 :                                           #
! #     Mark 10: August  2014, Green's function Regular      #
! #     Mark 11: January 2015, Green's function Linearized   #
! #                            Taylor, dethreaded, OpenMP    #
! #                                                          #
! #  Version 2.5 :                                           #
! #     Mark 12: December 2021, Code optimized for comp. eff.#
! #                             Separate routines for OBS/LAT#
! #                             Wavelength-independent ops   #
! #                             separated                    #
! #                                                          #
! ############################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              TWOSTREAM_BEAM_GEOMETRY_PREPARE                #
! #                                                             #
! ###############################################################

module twostream_geometry_m


PUBLIC

contains

SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE &
           ( MAXLAYERS,                                  & ! Dimension
             NLAYERS, SZA_GEOM_TRUE, REARTH, HEIGHTS,    & ! Inputs
             CHAPMAN_FACTORS )                             ! Output

      IMPLICIT NONE

!  Precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Generate path CHAPMAN_FACTORS for a curved ray-traced beam through a multilayer atmosphere.

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)  :: MAXLAYERS

!  Number of layers

      INTEGER, INTENT(IN)  :: NLAYERS

!  True solar zenith angle (degrees)

      REAL(kind=dp), INTENT(IN)  :: SZA_GEOM_TRUE

!  Earth radius (km)

      REAL(kind=dp), INTENT(IN)  :: REARTH

!  Coarse grids of heights, pressures and temperatures

      REAL(kind=dp), INTENT(IN)  :: HEIGHTS (0:MAXLAYERS)

!  Output arguments
!  ================

!  Path segments distances (km)

      REAL(kind=dp), INTENT(OUT) :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS)

!  Local variables
!  ===============

!  Local height array

      REAL(kind=dp) :: DELZ(MAXLAYERS)

!  Help variables

      INTEGER       :: N, K
      REAL(kind=dp) :: GM_TOA, TH_TOA, MU_TOA, DEG_TO_RAD
      REAL(kind=dp) :: STH1, SINTH1, STH2, SINTH2, PHI, SINPHI, &
                       RE_LOWER, RE_UPPER, DIST, H0, HN1, HN, DELZK

!  Some setup operations
!  =====================

!  Earth radii and heights differences

      H0 = HEIGHTS(0)
      HN1 = H0
      DO N = 1, NLAYERS
        HN = HEIGHTS(N)
        DELZ(N) = HN1-HN
        HN1 = HN
      ENDDO

!  TOA values

      DEG_TO_RAD = ATAN(1.0_DP) / 45.0_dp
      TH_TOA = SZA_GEOM_TRUE * DEG_TO_RAD
      MU_TOA = COS(TH_TOA)
      GM_TOA = SQRT ( 1.0_DP - MU_TOA * MU_TOA )

!  Straight line geometry
!  ======================

      H0 = H0 + REARTH
      DO N = 1, NLAYERS

!  Start values

        HN = HEIGHTS(N) + REARTH
        SINTH1 = GM_TOA * HN / H0
        STH1   = DASIN(SINTH1)
        RE_UPPER = H0

!  Loop over layers K from 1 to layer N

        DO K = 1, N

!  Sine-rule; PHI = earth-centered angle

          DELZK = DELZ(K)
          RE_LOWER = RE_UPPER - DELZK
          SINTH2 = RE_UPPER * SINTH1 / RE_LOWER
          STH2   = ASIN(SINTH2)
          PHI    = STH2 - STH1
          SINPHI = SIN(PHI)
          DIST = RE_UPPER * SINPHI / SINTH2
          CHAPMAN_FACTORS(K,N) = DIST / DELZK

!  Re-set

          RE_UPPER = RE_LOWER
          SINTH1 = SINTH2
          STH1   = STH2

        ENDDO

        DO K = N+1, NLAYERS
          CHAPMAN_FACTORS(K,N) = zero
        ENDDO

!  Finish main layer loop

      ENDDO

!  End of routine

      RETURN
END SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE

end module twostream_geometry_m
