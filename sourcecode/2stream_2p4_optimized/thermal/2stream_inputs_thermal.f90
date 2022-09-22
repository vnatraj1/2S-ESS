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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            TWOSTREAM_INPUTS_THERMAL                         #
! #            TWOSTREAM_CHECK_INPUT_DIMS_THERMAL               #
! #            TWOSTREAM_CHECK_INPUTS_BASIC_THERMAL             #
! #                                                             #
! ###############################################################

module twostream_inputs_thermal_m

USE twostream_brdf_supplement_thermal_m
USE twostream_geometry_thermal_m


      PRIVATE
      PUBLIC :: TWOSTREAM_INPUTS_THERMAL, TWOSTREAM_CHECK_INPUT_DIMS_THERMAL, & 
                TWOSTREAM_CHECK_INPUTS_BASIC_THERMAL

      CONTAINS

SUBROUTINE TWOSTREAM_INPUTS_THERMAL &
        ( MAXLAYERS, MAXTOTAL, MAX_USER_STREAMS, MAXMESSAGES,             & ! Dimensions
          DO_UPWELLING, DO_DNWELLING,                                     & ! Inputs
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_SURFACE_EMISSION,        & ! Inputs
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_STREAMS, USER_VZANGLES,   & ! Inputs
          MAXSTREAMS_BRDF, NSTREAMS_BRDF, DO_BRDF_SURFACE,                & ! BRDF Inputs
          DEG_TO_RAD, PI4, DO_POSTPROCESSING, DO_INCLUDE_MVOUT,           & ! Outputs
          USER_STREAMS, USER_SECANTS, PXSQ, SURFACE_FACTOR, DELTA_FACTOR, & ! Outputs
          PIE, STREAM_SINE, USER_SINES, NBRDF_HALF,                       & ! BRDF Outputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF, & ! BRDF Outputs
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS )           ! Exception handling

      IMPLICIT NONE

!  Precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: one = 1.0_dp

!  Dimensions

      INTEGER, INTENT(IN)          :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)          :: MAX_USER_STREAMS, MAXMESSAGES

!  Directional Flags

      LOGICAL, INTENT(IN)          :: DO_UPWELLING, DO_DNWELLING

!  Flux option flags

      LOGICAL, INTENT(IN)          :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)          :: DO_ADDITIONAL_MVOUT

!  Source control

      LOGICAL, INTENT(IN)          :: DO_SURFACE_EMISSION

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)          :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)    :: STREAM_VALUE

!  VZA input

      INTEGER, INTENT(IN)          :: N_USER_STREAMS
      REAL(kind=dp), INTENT(IN)    :: USER_VZANGLES(MAX_USER_STREAMS)

!  Constants

      REAL(kind=dp), INTENT(OUT)   :: DEG_TO_RAD, PI4

!  Post-processing flag

      LOGICAL, INTENT(OUT)         :: DO_POSTPROCESSING

!  Flux output flag

      LOGICAL, INTENT(OUT)         :: DO_INCLUDE_MVOUT

!  Cosines and sines

      REAL(kind=dp), INTENT(OUT)   :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT)   :: USER_SECANTS ( MAX_USER_STREAMS )

!  Polynomial expansion

      REAL(kind=dp), INTENT(OUT )  :: PXSQ

!  Surface reflectance factors

      REAL(kind=dp), INTENT(OUT)   :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(OUT)   :: DELTA_FACTOR

!  BRDF variables
!  ==============

!  Dimensions

      INTEGER, INTENT(IN)          :: MAXSTREAMS_BRDF

!  Number of azimuth quadrature streams for BRDF   

      INTEGER, INTENT(IN)          :: NSTREAMS_BRDF

!  Flag

      LOGICAL, INTENT(IN)          :: DO_BRDF_SURFACE

!  Constants
            
      REAL(kind=dp), INTENT(OUT)   :: PIE

!  2-Stream angle sine

      REAL(kind=dp), INTENT(OUT)   :: STREAM_SINE

!  Viewing zenith sine

      REAL(kind=dp), INTENT(OUT)   :: USER_SINES(MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER, INTENT(OUT)         :: NBRDF_HALF
      REAL(kind=dp), INTENT(OUT)   :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(kind=dp), INTENT(OUT)   :: BAX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: CXE_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Exception handling
!  ------------------

!    Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER          :: STATUS_SUB
      INTEGER          :: UM, I, I1, K
      REAL(kind=dp)    :: X_BRDFI, USI

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Check input dimensions
!  ----------------------

      CALL TWOSTREAM_CHECK_INPUT_DIMS_THERMAL &
      ( MAXMESSAGES, MAXLAYERS, MAXTOTAL, MAX_USER_STREAMS, & ! Dimensions
        NLAYERS, NTOTAL, N_USER_STREAMS,                    & ! Inputs
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )      ! Outputs

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Constants
!  ---------

      DEG_TO_RAD = ACOS( - one ) / 180.0_dp
      PI4 = DEG_TO_RAD * 720.0_dp

!  Input checking
!  ==============

!  Check input Basic.
!  SS inputs are omitted in this version........
!  Includes check on Flux output flags, and setting of Post-Processing flag

      CALL TWOSTREAM_CHECK_INPUTS_BASIC_THERMAL &
        ( MAXMESSAGES, MAX_USER_STREAMS,      & ! Dimensions
          DO_UPWELLING, DO_DNWELLING,         & ! Inputs
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, & ! Inputs
          N_USER_STREAMS, USER_VZANGLES,      & ! Inputs
          DO_POSTPROCESSING,                  & ! Output
          STATUS_SUB, C_NMESSAGES,            & ! Outputs
          C_MESSAGES, C_ACTIONS )               ! Outputs

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Get derived inputs
!  ==================

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
         DO I = 1, N_USER_STREAMS
            USI = COS(DEG_TO_RAD * USER_VZANGLES(I))
            USER_STREAMS(I) = USI
            USER_SECANTS(I) = ONE / USI
         ENDDO
      ENDIF

!  Auxiliary geometry quantities

      CALL TWOSTREAM_AUXGEOM_THERMAL &
            ( STREAM_VALUE, & ! Input
              PXSQ )          ! Output

!  Set flags
!  ---------

! Inclusion of mean value calculation
! Control for the Flux calculation

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        DO_INCLUDE_MVOUT = .TRUE.
      ENDIF

!  Surface reflectance factors

      SURFACE_FACTOR = 2.0_dp
      DELTA_FACTOR   = ONE

!  BRDF set up

      IF ( DO_BRDF_SURFACE ) THEN

!  Constants

         PIE = DEG_TO_RAD * 180.0_dp
         STREAM_SINE = SQRT(ONE - STREAM_VALUE * STREAM_VALUE)

!  Half number of moments

         NBRDF_HALF = NSTREAMS_BRDF / 2

!  Viewing angles

         DO UM = 1, N_USER_STREAMS
            USER_SINES(UM)   = SQRT(1.0_dp-USER_STREAMS(UM)*USER_STREAMS(UM))
         ENDDO

!  BRDF quadrature

         CALL TWOSTREAM_GAULEG ( 0.0_dp, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
         DO I = 1, NBRDF_HALF
            I1 = I + NBRDF_HALF
            X_BRDFI = X_BRDF(I)
            X_BRDF(I1) = - X_BRDFI
            A_BRDF(I1) =   A_BRDF(I)
            CXE_BRDF(I) = X_BRDFI
            SXE_BRDF(I) = SQRT(ONE-X_BRDFI*X_BRDFI)
         ENDDO
         DO I = 1, NSTREAMS_BRDF
            X_BRDFI = X_BRDF(I)
            X_BRDF(I) = PIE * X_BRDFI
            CX_BRDF(I) = COS (X_BRDFI)
            SX_BRDF(I) = SIN (X_BRDFI)
         ENDDO

!  Half space cosine-weight arrays (emission only, non-Lambertian)

         IF ( DO_SURFACE_EMISSION ) THEN
            DO K = 1, NBRDF_HALF
               BAX_BRDF(K) = X_BRDF(K) * A_BRDF(K) / PIE
            ENDDO
         ENDIF

      ENDIF

!  End BRDF set up

!  Finish

      RETURN

END SUBROUTINE TWOSTREAM_INPUTS_THERMAL

SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS_THERMAL &
      ( MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAX_USER_STREAMS, &
        NLAYERS,   NTOTAL, N_USER_STREAMS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES
      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXTOTAL
      INTEGER, INTENT(IN)  ::  MAX_USER_STREAMS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  NLAYERS
      INTEGER, INTENT(IN)  ::  NTOTAL
      INTEGER, INTENT(IN)  ::  N_USER_STREAMS

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGES(0:MAXMESSAGES)
      CHARACTER*(*), INTENT(INOUT) :: ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = 0
      NM = NMESSAGES

!  Check active input dimensions against maximum dimensions
!  ========================================================

!  1. Basic dimensions - always checked

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension'
        STATUS       = 1
      ENDIF

      IF ( NTOTAL .GT. MAXTOTAL ) THEN
        NM = NM + 1
        MESSAGES(NM) = '2*Number of layers NTOTAL > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXTOTAL dimension'
        STATUS       = 1
      ENDIF

!  2. Geometry dimensions - always checked

      IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Observation angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
        STATUS       = 1
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish
      RETURN

END SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS_THERMAL

SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC_THERMAL &
         ( MAXMESSAGES, MAX_USER_STREAMS,       & ! Dimensions
           DO_UPWELLING, DO_DNWELLING,          & ! Inputs
           DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,  & ! Inputs
           N_USER_STREAMS, USER_VZANGLES,       & ! Inputs
           DO_POSTPROCESSING,                   & ! Output
           STATUS, NMESSAGES, MESSAGE, ACTION )   ! Outputs, Input/Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ------

!  Dimensions :

      INTEGER, INTENT(IN) :: MAXMESSAGES, MAX_USER_STREAMS

!  Flags

      LOGICAL, INTENT(IN) :: DO_UPWELLING, DO_DNWELLING

!  Flux output flags

      LOGICAL, INTENT(IN)    :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)    :: DO_ADDITIONAL_MVOUT

!  VZA input

      INTEGER, INTENT(IN) :: N_USER_STREAMS
      REAL(kind=dp), INTENT(IN) :: USER_VZANGLES(MAX_USER_STREAMS)

!  Processing flag

      LOGICAL, INTENT(OUT) :: DO_POSTPROCESSING

!  Module output

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGE(MAXMESSAGES)
      CHARACTER*(*), INTENT(INOUT) :: ACTION(MAXMESSAGES)

!  local variables

      INTEGER           :: I, NM
      CHARACTER(LEN=2)  :: C2
      LOGICAL           :: LOOP

!  Initialize output status

      STATUS = 0
      NM = NMESSAGES

!  Initialize post-processing flag

      DO_POSTPROCESSING = .false.

!  Check MAX number of messages should be at least 8 (covers all errors below as well as those checking 
!  input dimensions )

      IF ( MAXMESSAGES .LT. 8 ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: Not enough possible messages in Checks'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAXMESSAGES to 13'
        STATUS = 1
      ENDIF

!  No point in going on if dimension checks have failed

      IF ( STATUS .EQ. 1 .AND. NM .GT. 0 ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: At least one dimensioning check failed'
        ACTION(NM)  = 'Read previous messages to determine actions'
        NMESSAGES = NM
        RETURN
      ENDIF

!  Check directional input

      IF ( .NOT. DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: no directional input is set'
        ACTION(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: set one!'
        STATUS = 1
      ENDIF

!  Check FLUX flags, set post-processing

      IF ( DO_MVOUT_ONLY .AND. DO_ADDITIONAL_MVOUT ) then
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: both Flux-output flags are set'
        ACTION(NM)  = 'Abort: Turn off 1 of the MVOUT flags'
        STATUS = 1
      ENDIF

      IF ( DO_ADDITIONAL_MVOUT .OR. .NOT. DO_MVOUT_ONLY ) then
        DO_POSTPROCESSING = .TRUE.
      ENDIF

!  Check user-defined stream angles (should always be [0,90])
!  Avoid this section if MVOUT_ONLY

      IF ( .NOT. DO_MVOUT_ONLY ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE ( LOOP .AND. I .LT. N_USER_STREAMS )
          I = I + 1
          IF ( USER_VZANGLES(I) .GT. 90.0D0 .or.USER_VZANGLES(I) .LT. 0.0D0 ) THEN
            NM = NM + 1
            WRITE(C2,'(I2)')I  
            MESSAGE(NM)='Bad input: out-of-range viewing angle, no. '//C2
            ACTION(NM) = 'Look at viewing angle input, range [0,90]'
            LOOP = .FALSE.
            STATUS = 1
          ENDIF
        ENDDO
      ENDIF

!  Set number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC_THERMAL

end module twostream_inputs_thermal_m
