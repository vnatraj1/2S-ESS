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
! #            TWOSTREAM_CHECK_INPUT_DIMS                       #
! #            TWOSTREAM_CHECK_INPUTS_BASIC                     #
! #                                                             #
! ###############################################################

module twostream_inputs_m

      PRIVATE
      PUBLIC :: TWOSTREAM_CHECK_INPUT_DIMS, TWOSTREAM_CHECK_INPUTS_BASIC

      CONTAINS

SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS &
      ( DO_MVOUT_ONLY, DO_USER_OBSGEOMS, MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL,   NBEAMS,   N_USER_STREAMS,   N_USER_RELAZMS,   N_USER_OBSGEOMS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Control variables

      LOGICAL, INTENT(IN)  ::  DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  ::  DO_USER_OBSGEOMS

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES

      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXTOTAL
      INTEGER, INTENT(IN)  ::  MAXBEAMS
      INTEGER, INTENT(IN)  ::  MAX_USER_STREAMS
      INTEGER, INTENT(IN)  ::  MAX_USER_RELAZMS
      INTEGER, INTENT(IN)  ::  MAX_USER_OBSGEOMS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  NLAYERS
      INTEGER, INTENT(IN)  ::  NTOTAL
      INTEGER, INTENT(IN)  ::  NBEAMS
      INTEGER, INTENT(IN)  ::  N_USER_STREAMS
      INTEGER, INTENT(IN)  ::  N_USER_RELAZMS
      INTEGER, INTENT(IN)  ::  N_USER_OBSGEOMS

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

!  2a. Geometry dimensions - always checked

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles NBEAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS       = 1
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of relative azimuths N_USER_RELAZMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = 1
      ENDIF

!  2b. Geometry dimensions - conditionally checked

      IF ( .NOT. DO_MVOUT_ONLY ) THEN
        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of user streams N_USER_STREAMS > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS       = 1
        ENDIF
      ENDIF

      IF ( DO_USER_OBSGEOMS ) THEN
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Observation Geometries > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = 1
        ENDIF
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

END SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS

SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC &
         ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,              & ! Dimensions !@@
           MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,           & ! Dimensions
           DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,          & ! inputs
           DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                  & ! inputs
           DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,  & ! Input !@@ New line, 2p3
           DO_USER_OBSGEOMS, N_USER_OBSGEOMS, USER_OBSGEOMS,       & ! Input !@@ New
           NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,        & ! inputs
           BEAM_SZAS, USER_ANGLES, USER_RELAZMS,                   & ! inputs
           EARTH_RADIUS, HEIGHT_GRID,                              & ! inputs
           STATUS, NMESSAGES, MESSAGE, ACTION )                      ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  Notes 05 November 2013. Flux output flags. Version 2p3

!  Inputs
!  ------

!  Dimensions :

      INTEGER, INTENT(IN) :: MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS  !@@
      INTEGER, INTENT(IN) :: MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS

!  Flags

      LOGICAL, INTENT(IN) :: DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN) :: DO_SOLAR_SOURCES, DO_THERMAL_EMISSION

!  Single scatter flags omitted from this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSCORR_OUTGOING, DO_SSCORR_NADIR, DO_SSFULL

!  !@@ Version 2p3, 11/5/13. Flux output flags, processing flag  

      LOGICAL, INTENT(IN)    :: DO_MVOUT_ONLY        !@@
      LOGICAL, INTENT(IN)    :: DO_ADDITIONAL_MVOUT  !@@
      LOGICAL, INTENT(INOUT) :: DO_POSTPROCESSING    !@@ Will always be set here.

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@

      LOGICAL, INTENT(IN) :: DO_USER_OBSGEOMS !@@
      INTEGER, INTENT(IN) :: N_USER_OBSGEOMS  !@@
      REAL(kind=dp), INTENT(IN) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@

!  Number of layers

      INTEGER, INTENT(IN) :: NLAYERS

!  Angle Numbers. [Now Intent(inout), thanks to option for ObsGeom !@@]

      INTEGER, INTENT(INOUT) :: NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

!  Geometry. [Now Intent(inout), thanks to option for ObsGeom !@@]

      REAL(kind=dp), INTENT(INOUT) :: BEAM_SZAS    ( MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: USER_ANGLES  ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(INOUT) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  height and earth radius

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )

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

!  !@@ 2p3 Initialize post-processing flag

      DO_POSTPROCESSING = .false.

!  Single scattering stuff omitted
!      IF ( DO_SSCORR_OUTGOING ) THEN
!        IF ( NFINELAYERS .GT. 10 ) THEN
!          MAIL   = 'Number of fine layers > 10'
!          ACTION = 'Re-set input value'
!          STATUS = 1
!          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!          RETURN
!        ENDIF
!      ENDIF

!  May 2010, MARK II. Dimensioning checks removed

!      IF ( NMOMENTS_INPUT .GT. 300 ) THEN
!      IF ( NBEAMS .GT. 10 ) THEN
!      IF ( N_USER_STREAMS .GT. 10 ) THEN
!      IF ( N_USER_RELAZMS .GT. 10 ) THEN

!  Check MAX number of messages should be at least 13 (covers all errors below)

      IF ( MAXMESSAGES .LT. 14 ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: Not enough possible messages in Checks'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAXMESSAGES to 13'
        STATUS = 1
      ENDIF

!  November 2012. Dimensioning check re-introduced (Properly this time)

      IF ( NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: Number of layers NLAYERS > Maximum dimension MAXLAYERS'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAXLAYERS'
        STATUS = 1
      ENDIF

! !@@ New 12/21/12. Observational Geometry check

      if ( DO_USER_OBSGEOMS ) THEN
        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: Number of User ObsGeoms N_USER_OBSGEOMS > Maximum dimension'
          ACTION(NM)  = 'Increase value of symbolic Dimension MAX_USER_OBSGEOMS'
          STATUS = 1 ; go to 5665
        ENDIF
      ENDIF

!  !@@ Skip Next 3 checks if using observational geometry

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: Number of beams NBEAMS > Maximum dimension MAXBEAMS'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAXBEAMS'
        STATUS = 1  
      ENDIF

      IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: Number of User streams N_USER_STREAMS > Maximum dimension'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAX_USER_STREAMS'
        STATUS = 1
      ENDIF

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1  
        MESSAGE(NM) = 'Bad input: Number of User azimuths N_USER_RELAZMS > Maximum dimension'
        ACTION(NM)  = 'Increase value of symbolic Dimension MAX_USER_RELAZMS'
        STATUS = 1
      ENDIF

!  !@@ Continuation point for skipping normal geometry control checks

5665  continue

!  !@@ No point in going on if dimemnsion checks have failed

      if ( status .eq. 1 .and. NM.gt.0) then
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: At least one dimensioning check failed'
        ACTION(NM)  = 'Read previous messages to determine actions'
        nmessages = nm
        return
      endif

!  !@@ Reset angle input in Observational Geometry mode
!  ====================================================

! @@ Note differing treatment for Thermal-emission-only

      IF ( DO_USER_OBSGEOMS ) THEN
         IF ( DO_SOLAR_SOURCES ) THEN
            NBEAMS          = N_USER_OBSGEOMS
            N_USER_STREAMS  = N_USER_OBSGEOMS
            N_USER_RELAZMS  = N_USER_OBSGEOMS
            BEAM_SZAS    (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
            USER_ANGLES  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
            USER_RELAZMS (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
         ELSE
            NBEAMS          = 1
            N_USER_STREAMS  = N_USER_OBSGEOMS
            N_USER_RELAZMS  = 1
            USER_ANGLES  (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         ENDIF
      ENDIF  

!  check directional input

      IF ( .NOT.DO_UPWELLING .AND. .NOT. DO_DNWELLING ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: no directional input is set'
        ACTION(NM)  = 'Check DO_UPWELLING & DO_DNWELLING: set one!'
        STATUS = 1
      ENDIF

!  Thermal related checks
!  ======================

      IF ( .NOT.DO_SOLAR_SOURCES.AND..NOT.DO_THERMAL_EMISSION ) THEN
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: No solar or thermal sources'
        ACTION(NM)  = 'Abort: must set one of the source flags!'
        STATUS = 1
      ENDIF

!  Set number of solar sources NBEAMS to 1 for the thermal-only default

      IF ( .NOT.DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( NBEAMS .NE. 1 ) THEN
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: NBEAMS > 1 for thermal-only run'
          ACTION(NM)  = 'Abort: must set NBEAMS = 1 for thermal-only run'
          STATUS = 1
        ENDIF
      ENDIF

!  Set number of Azimuth angles to 1 for the thermal-only default

      IF ( .NOT. DO_SOLAR_SOURCES.AND.DO_THERMAL_EMISSION ) THEN
        IF ( N_USER_RELAZMS .NE. 1 ) THEN
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: N_USER_RELAZMS > 1 for thermal-only run'
          ACTION(NM)  = 'Abort: N_USER_RELAZMS = 1 for thermal-only run'
          STATUS = 1
        ENDIF
      ENDIF

!  11/5/13. Version 2p3, Check FLUX flags, set post-processing
!  -----------------------------------------------------------

      IF ( DO_MVOUT_ONLY .and. DO_ADDITIONAL_MVOUT ) then
        NM = NM + 1
        MESSAGE(NM) = 'Bad input: both Flux-output flags are set'
        ACTION(NM)  = 'Abort: Turn off 1 of the MVOUT flags'
        STATUS = 1
      ENDIF

      IF ( DO_ADDITIONAL_MVOUT .or. .not.DO_MVOUT_ONLY ) then
        DO_POSTPROCESSING = .true.
      ENDIF

!  check viewing geometry input
!  ============================

!  Check earth radius (Chapman function only)
!    ---WARNING. Default value of 6371.0 will be set

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( EARTH_RADIUS.LT.6320.0D0 .OR. EARTH_RADIUS.GT.6420.0D0 ) THEN
          NM = NM + 1
          MESSAGE(NM)= 'Bad input: Earth radius outside 6320-6420 km'
          ACTION(NM) = 'Warning: default value of 6371.0 km was set'
          EARTH_RADIUS = 6371.0D0
        ENDIF
      ENDIF

!  Check solar zenith angle input

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.NBEAMS)
        I = I + 1 
        IF ( BEAM_SZAS(I) .LT. 0.0d0 .OR. BEAM_SZAS(I).GE.90.0D0 ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM)= 'Bad input: out-of-range beam angle, no. '//C2
          ACTION(NM) = 'Look at BEAM_SZAS input, should be < 90 & > 0'
          LOOP = .FALSE.
          STATUS = 1 
        ENDIF
      ENDDO

!  Check relative azimuths
!  @@ 2p3 Avoid this section if MVOUT_ONLY

      IF ( .not.DO_MVOUT_ONLY ) THEN
        LOOP = .TRUE. ; I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_RELAZMS)
          I = I + 1
          IF ( USER_RELAZMS(I).GT.360.0D0 .OR. USER_RELAZMS(I).LT.0.0d0 ) THEN
            NM = NM + 1
            WRITE(C2,'(I2)')I
            MESSAGE(NM)='Bad input: out-of-range azimuth angle, no. '//C2
            ACTION(NM) = 'Look at azimuth angle input, range [0,360]'
            LOOP = .FALSE.
            STATUS = 1
          ENDIF
        ENDDO
      ENDIF

!  check user-defined stream angles (should always be [0,90])
!  @@ 2p3 Avoid this section if MVOUT_ONLY

      IF ( .not.DO_MVOUT_ONLY ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE (LOOP .AND. I.LT.N_USER_STREAMS)
          I = I + 1
          IF ( USER_ANGLES(I) .GT. 90.0 .or.USER_ANGLES(I) .LT. 0.0d0 ) THEN
            NM = NM + 1
            WRITE(C2,'(I2)')I  
            MESSAGE(NM)='Bad input: out-of-range viewing angle, no. '//C2
            ACTION(NM) = 'Look at viewing angle input, range [0,90]'
            LOOP = .FALSE.
            STATUS = 1
          ENDIF
        ENDDO
      ENDIF

!  Check height grid input (Chapman function only)

      LOOP = .TRUE.
      I = 0  
      DO WHILE (LOOP .AND. I.LT.NLAYERS)
        I = I + 1
        IF ( HEIGHT_GRID(I-1).LE.HEIGHT_GRID(I) ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM) = 'Bad input: Height-grid not monotonic decreasing; Layer '//C2
          ACTION(NM) = 'Look at Height-grid input'
          LOOP = .FALSE.
          STATUS = 1
        ENDIF
      ENDDO

!  set number of messages

      nmessages = nm 

!  Finish  

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC

end module twostream_inputs_m
