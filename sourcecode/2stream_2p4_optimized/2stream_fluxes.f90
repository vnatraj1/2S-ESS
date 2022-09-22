! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #          TWOSTREAM_FLUXES_SOLAR                         #
! #          TWOSTREAM_FLUXES_THERMAL                       #
! #                                                         #
! ###########################################################

module twostream_fluxes_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_FLUXES_SOLAR &
            ( DO_UPWELLING, DO_DNWELLING,                          & ! Flags
              DO_DIRECTBEAMB, PI4, STREAM_VALUE,                   & ! Inputs
              FLUXFAC, FLUXMULT, X0B, TRANS_SOLAR_BEAMB,           & ! Inputs
              LCON_XVEC21, MCON_XVEC21, EIGENTRANS1, WUPPER21,     & ! Inputs
              LCON_XVEC1NL, MCON_XVEC1NL, EIGENTRANSNL, WLOWER1NL, & ! Inputs
              FLUXES_TOA_SOLARB, FLUXES_BOA_SOLARB )                 ! Outputs

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA

      IMPLICIT NONE

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_DIRECTBEAMB

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXFAC, FLUXMULT, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Cosines and transmittance of solar beam (for the direct fluxes)

      REAL(kind=dp), INTENT(IN) :: X0B
      REAL(kind=dp), INTENT(IN) :: TRANS_SOLAR_BEAMB

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS1, EIGENTRANSNL

!  Solution constants of integration multiplied by eigensolutions

      REAL(kind=dp), INTENT(IN) :: LCON_XVEC21, LCON_XVEC1NL
      REAL(kind=dp), INTENT(IN) :: MCON_XVEC21, MCON_XVEC1NL

!  Particular solutions at layer boundaries

      REAL(kind=dp), INTENT(IN) :: WUPPER21
      REAL(kind=dp), INTENT(IN) :: WLOWER1NL

!  Flux output

      REAL(kind=dp), INTENT(OUT) :: FLUXES_TOA_SOLARB(2)
      REAL(kind=dp), INTENT(OUT) :: FLUXES_BOA_SOLARB(2)

!  Local variables

      REAL(kind=dp) :: PI2, SHOM, SPAR, QUADINTENS_TOA, QUADINTENS_BOA, DMEAN, DFLUX

!  Upwelling Flux at TOA

      PI2 = 0.5_dp * PI4
      IF ( DO_UPWELLING ) THEN
         SHOM = LCON_XVEC21 + MCON_XVEC21 * EIGENTRANS1
         SPAR = WUPPER21
         QUADINTENS_TOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_TOA_SOLARB(1) = 0.5_dp * QUADINTENS_TOA
         FLUXES_TOA_SOLARB(2) = PI2 * STREAM_VALUE * QUADINTENS_TOA
      ENDIF

!  Downwelling Flux at BOA

      IF ( DO_DNWELLING ) THEN
         SHOM = LCON_XVEC1NL * EIGENTRANSNL + MCON_XVEC1NL
         SPAR = WLOWER1NL
         QUADINTENS_BOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_BOA_SOLARB(1) = 0.5_dp * QUADINTENS_BOA             !actinic flux (diffuse)
         FLUXES_BOA_SOLARB(2) = PI2 * STREAM_VALUE * QUADINTENS_BOA !regular flux (diffuse)
         IF (DO_DIRECTBEAMB) THEN
            DMEAN = FLUXFAC * TRANS_SOLAR_BEAMB / PI4
            DFLUX = FLUXFAC * TRANS_SOLAR_BEAMB * X0B
            FLUXES_BOA_SOLARB(1) = FLUXES_BOA_SOLARB(1) + DMEAN         !actinic flux (difffuse + direct)
            FLUXES_BOA_SOLARB(2) = FLUXES_BOA_SOLARB(2) + DFLUX         !regular flux (difffuse + direct)
         ENDIF
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES_SOLAR

!

SUBROUTINE TWOSTREAM_FLUXES_THERMAL &
            ( DO_UPWELLING, DO_DNWELLING,                          & ! Flags
              PI4, STREAM_VALUE, FLUXMULT,                         & ! Inputs
              LCON_XVEC21, MCON_XVEC21, EIGENTRANS1, WUPPER21,     & ! Inputs
              LCON_XVEC1NL, MCON_XVEC1NL, EIGENTRANSNL, WLOWER1NL, & ! Inputs
              FLUXES_TOA_THERMAL, FLUXES_BOA_THERMAL )               ! Outputs

!  New routine 11/5/13. Diffuse Fluxes at TOA and BOA for thermal-only scenario

      IMPLICIT NONE

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT, PI4

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Eigen-Transmittance

      REAL(kind=dp), INTENT(IN) :: EIGENTRANS1, EIGENTRANSNL

!  Solution constants of integration multiplied by eigensolutions

      REAL(kind=dp), INTENT(IN) :: LCON_XVEC21, LCON_XVEC1NL
      REAL(kind=dp), INTENT(IN) :: MCON_XVEC21, MCON_XVEC1NL

!  Particular solutions at layer boundaries

      REAL(kind=dp), INTENT(IN) :: WUPPER21
      REAL(kind=dp), INTENT(IN) :: WLOWER1NL

!  Flux output

      REAL(kind=dp), INTENT(OUT) :: FLUXES_TOA_THERMAL(2)
      REAL(kind=dp), INTENT(OUT) :: FLUXES_BOA_THERMAL(2)

!  Local variables

      REAL(kind=dp) :: PI2, SHOM, SPAR, QUADINTENS_TOA, QUADINTENS_BOA

!  Upwelling Flux at TOA

      PI2 = 0.5_dp * PI4
      IF ( DO_UPWELLING ) THEN
         SHOM = LCON_XVEC21 + MCON_XVEC21 * EIGENTRANS1
         SPAR = WUPPER21
         QUADINTENS_TOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_TOA_THERMAL(1) = 0.5_dp * QUADINTENS_TOA
         FLUXES_TOA_THERMAL(2) = PI2 * STREAM_VALUE * QUADINTENS_TOA
      ENDIF

!  Downwelling Flux at BOA

      IF ( DO_DNWELLING ) THEN
         SHOM = LCON_XVEC1NL * EIGENTRANSNL + MCON_XVEC1NL
         SPAR = WLOWER1NL
         QUADINTENS_BOA = FLUXMULT * ( SPAR + SHOM )
         FLUXES_BOA_THERMAL(1) = 0.5_dp * QUADINTENS_BOA             !actinic flux (diffuse)
         FLUXES_BOA_THERMAL(2) = PI2 * STREAM_VALUE * QUADINTENS_BOA !regular flux (diffuse)
      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_FLUXES_THERMAL

end module twostream_fluxes_m
