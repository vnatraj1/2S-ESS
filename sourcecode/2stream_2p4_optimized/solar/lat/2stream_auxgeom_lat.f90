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
! #            TWOSTREAM_CHECK_INPUTS_BASIC                     #
! #            TWOSTREAM_CHECK_INPUTS_OPTICAL                   #
! #                                                             #
! #              TWOSTREAM_AUXGEOM                              #
! #              TWOSTREAM_QSPREP                               #
! #              TWOSTREAM_PREPTRANS                            #
! #              TWOSTREAM_DIRECTBEAM                           #
! #              TWOSTREAM_EMULTMASTER                          #
! #                                                             #
! #              TWOSTREAM_BEAM_GEOMETRY_PREPARE                #
! #                                                             #
! ###############################################################

module twostream_miscsetups_m

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains

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

!  check SS inputs (both file-read and derived)
!  --------------------------------------------

!  This section commented out in the streamlined version

!  Check sphericity corrections....Cannot both be turned on
!      IF ( DO_SSCORR_NADIR .and. DO_SSCORR_OUTGOING ) THEN
!        MAIL   = 'Cannot have both single scatter corrections on'
!        ACTION = 'Turn off DO_SSCORR_NADIR and/or DO_SSCORR_OUTGOING'
!        STATUS = 1
!        CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!      ENDIF

!   Single scatter corrections must be turned on
!      IF ( DO_SSFULL ) THEN
!        IF ( .not.DO_SSCORR_NADIR.and..not.DO_SSCORR_OUTGOING ) THEN
!          MAIL   = 'Bad input: Full SS, must have one SSCORR flag set'
!          ACTION = 'Full SS: default to use outgoing SS correction'
!          DO_SSCORR_NADIR    = .FALSE.
!          DO_SSCORR_OUTGOING = .TRUE.
!          CALL TWOSTREAM_ERROR_TRACE ( INIT, MAIL, ACTION, STATUS )
!        ENDIF
!      ENDIF

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

!

SUBROUTINE TWOSTREAM_CHECK_INPUTS_OPTICAL &
           ( MAXLAYERS, MAXMESSAGES, NLAYERS,        & ! input
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,  & ! input
             STATUS, NMESSAGES, MESSAGE, ACTION )      ! output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Module inputs
!  -------------

!  Dimensions :
     
      INTEGER, INTENT(IN) :: MAXLAYERS, MAXMESSAGES

!  Numbers

      INTEGER, INTENT(IN) :: NLAYERS

!  Optical properties

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT (MAXLAYERS )
      REAL(kind=dp), INTENT(IN) :: OMEGA_TOTAL (MAXLAYERS )
      REAL(kind=dp), INTENT(IN) :: ASYMM_TOTAL (MAXLAYERS )

!  Module output
!  -------------

      INTEGER      , INTENT(OUT)   :: STATUS
      INTEGER      , INTENT(INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT(INOUT) :: MESSAGE (MAXMESSAGES)
      CHARACTER*(*), INTENT(INOUT) :: ACTION (MAXMESSAGES)

!  local variables

      INTEGER           :: L, NM
      CHARACTER(LEN=3)  :: C3

!  Initialize output status

      STATUS = 0
      NM = NMESSAGES

!  check optical property inputs
!  -----------------------------

!  Check non-negative optical thickness values

      DO L = 1, NLAYERS
        IF ( DELTAU_VERT(L).LE.0.0_dp ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTION(NM)  = 'Check opt-thickness input '
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for conservative scattering limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L).GT.0.999999d0 ) THEN !original
        IF ( OMEGA_TOTAL(L).GT.0.999999999d0 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo close to 1, layer '//C3
          ACTION (NM) = 'Check single scattering albedo input'
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for smallness limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L).LT.1.0d-06 ) THEN !original
        IF ( OMEGA_TOTAL(L).LT.1.0d-9 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo too small, layer '//C3
          ACTION (NM) = 'Check single scattering albedo input'
          STATUS = 1
        ENDIF
      ENDDO

!  check asymmetry parameter, between -1 and 1

      DO L = 1, NLAYERS
        IF ( ASYMM_TOTAL(L).LE.-1.0d0 .OR. &
             ASYMM_TOTAL(L).GE. 1.0d0  ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: Asymm parameter outside [-1,1], layer'//C3
          ACTION(NM)  = 'Check Asymmetry parameter input'
          STATUS = 1
        ENDIF
      ENDDO

!  set number of messages

      nmessages = nm

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_OPTICAL

!

SUBROUTINE TWOSTREAM_AUXGEOM &
      ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING,  & ! Dimensions, Flag
        N_USER_STREAMS, NBEAMS, FOURIER, & ! inputs
        X0, USER_STREAMS, STREAM_VALUE,  & ! inputs
        PX11, PXSQ, POX, PX0X, ULP )       ! outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        ::  MAX_USER_STREAMS, MAXBEAMS

!  Flag for post-processing, @@ 2p3, 11/5/13

      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Numbers

      INTEGER, INTENT(IN)        ::  N_USER_STREAMS, NBEAMS

!  Fourier

      INTEGER, INTENT(IN)        ::  FOURIER

!  stream directions

      REAL(kind=dp), INTENT(IN)  :: X0 ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

      REAL(kind=dp), INTENT(OUT) :: ULP ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT) :: POX  ( MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: PX0X ( MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: PXSQ, PX11

!  Local variables

      INTEGER       :: UM, IBEAM
      REAL(kind=dp) :: MU, MU0, POXIB

!  Saved quantities

      IF ( DO_POSTPROCESSING ) THEN
        DO UM = 1, N_USER_STREAMS
          MU = USER_STREAMS(UM)
          ULP(UM) =  -DSQRT(0.5d0*(1.0d0-MU*MU))
        ENDDO
      ELSE
        ULP = 0.0d0
      ENDIF

      if ( fourier.eq.0) then
        PXSQ = STREAM_VALUE * STREAM_VALUE
      else if ( fourier.eq.1 ) then
        PX11 = DSQRT(0.5d0*(1.0d0-STREAM_VALUE*STREAM_VALUE))
        PXSQ = PX11 * PX11
      endif

      DO IBEAM = 1, NBEAMS
        MU0 = X0(IBEAM)
        POXIB = DSQRT(0.5d0*(1.0d0-MU0*MU0))
        POX(IBEAM) = POXIB
        if ( fourier.eq.0) then
          PX0X(IBEAM) = MU0 * STREAM_VALUE
        else
          PX0X(IBEAM) = POXIB * PX11
        endif
      ENDDO

!  FInish

      RETURN
END SUBROUTINE TWOSTREAM_AUXGEOM

!

SUBROUTINE TWOSTREAM_QSPREP &
       ( MAXLAYERS, MAXBEAMS,                              & ! Dimensions
         NLAYERS, NBEAMS, DO_PLANE_PARALLEL,               & ! Input
         DELTAU_VERT, CHAPMAN_FACTORS, X0,                 & ! Input
         DO_REFLECTED_DIRECTBEAM,                          & ! In/Out
         INITIAL_TRANS, AVERAGE_SECANT, AVERAGE_SECANT_PP, & ! Output
         LAYER_PIS_CUTOFF, TRANS_SOLAR_BEAM )                ! Output

      implicit none

!   precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Inputs
!  ------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXBEAMS

!  Control

      INTEGER, INTENT(IN)        :: NLAYERS, NBEAMS

!  optical thickness input

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT ( MAXLAYERS )

!  Path segments distances (km)

      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS,MAXBEAMS)

!  Plane parallel control

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  Beam SZA cosines

      REAL(kind=dp), INTENT(IN)  :: X0(NBEAMS)

!  Output
!  ------

!  Reflectance flags
!mick fix 1/24/12 - Set before subroutine.  Possibly reset below.
      LOGICAL, INTENT(INOUT)     :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(OUT)       :: LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(OUT) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: AVERAGE_SECANT_PP ( MAXBEAMS )

!  Solar beam attenuations

      REAL(kind=dp), INTENT(OUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )

!  Local variables
!  ---------------

!  Derived optical thickness inputs

      INTEGER       ::  N, K, IB
      REAL(kind=dp) :: S_T_0, S_T_1, SEC0, DELS
      REAL(kind=dp) :: TAUGRID, TAUGRIDN1, TAU_SOLAR(MAXBEAMS)
      REAL(kind=dp) :: TAUSLANT, TAUSLANTN1
      
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

!  plane-parallel case
!  -------------------

      IF ( DO_PLANE_PARALLEL ) THEN

       DO IB = 1, NBEAMS

        SEC0 = one / X0(IB)
        LAYER_PIS_CUTOFF(IB) = NLAYERS
        TAUGRID = zero
        DO N = 1, NLAYERS
          TAUGRIDN1 = TAUGRID
          TAUGRID = TAUGRID + DELTAU_VERT(N)
          TAUSLANT = TAUGRID * SEC0
          IF ( N .LE. LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ENDIF
            INITIAL_TRANS(N,IB)  = EXP ( - TAUGRIDN1 * SEC0 )
          ENDIF
        ENDDO
        AVERAGE_SECANT_PP(IB) = SEC0
        TAU_SOLAR(IB) = TAUSLANT

!  Set Direct Beam Flag and solar beam total attenuation to surface

        IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
          TRANS_SOLAR_BEAM(IB)        = zero
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          TRANS_SOLAR_BEAM(IB) = EXP( - TAUSLANT )
!mick fix 1/24/12 - Set before subroutine
          !DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF

       ENDDO

      ELSE

!  pseudo-spherical case
!  ---------------------

       DO IB = 1, NBEAMS

        S_T_0 = 1.0_dp
        S_T_1 = zero
        LAYER_PIS_CUTOFF(IB) = NLAYERS

!  Get the total spherical attenuation from layer thickness sums and set up the average secant 
!  formulation

        TAUSLANT = zero
        DO N = 1, NLAYERS
          TAUSLANTN1 = TAUSLANT
          TAUSLANT = DOT_PRODUCT(DELTAU_VERT(1:N),CHAPMAN_FACTORS(1:N,N,IB))
          IF  (N.LE.LAYER_PIS_CUTOFF(IB) ) THEN
            IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
              LAYER_PIS_CUTOFF(IB) = N
            ELSE
              S_T_1 = EXP ( - TAUSLANT )
            ENDIF
            AVERAGE_SECANT(N,IB) = &
                (TAUSLANT-TAUSLANTN1) / DELTAU_VERT(N)
            INITIAL_TRANS(N,IB)  = S_T_0
            S_T_0             = S_T_1
          ELSE
            AVERAGE_SECANT(N,IB) = zero
            INITIAL_TRANS(N,IB)  = zero
          ENDIF
        ENDDO
        TAU_SOLAR(IB) = TAUSLANT

!  Set Direct Beam Flag and solar beam total attenuation to surface

        IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
          TRANS_SOLAR_BEAM(IB)        = zero
          DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
        ELSE
          TRANS_SOLAR_BEAM(IB) = EXP( - TAUSLANT )
!mick fix 1/24/12 - Set before subroutine
          !DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
        ENDIF

       ENDDO

      ENDIF

!  debug
!      do N = 1, nlayers
!      write(*,'(i3,1p2e20.10)')n,initial_trans(n,1),average_secant(n,1)
!      enddo

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_QSPREP

!

SUBROUTINE TWOSTREAM_PREPTRANS &
    ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS,                       & ! Dimensions
      DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_POSTPROCESSING,       & ! Input flags (2p1,2p3)
      NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_SECANTS,  & ! Input
      INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,             & ! Input
      T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM  )                    ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Prepare transmittances and transmittance factors

!  Inputs
!  ------

!  Dimensions

      INTEGER, INTENT(IN)       :: MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS

!  Control Flags.
!    !@@ Add observational geometry option 12/21/12
!    !@@ Add Post-processing flag, 11/5/13
!    !@@ Add solar sources flag 6/29/20 VN

      LOGICAL, intent(in)       :: DO_SOLAR_SOURCES
      LOGICAL, intent(in)       :: DO_USER_OBSGEOMS   !@@
      LOGICAL, intent(in)       :: DO_POSTPROCESSING  !@@

! Numbers

      INTEGER, INTENT(IN)       :: NLAYERS, N_USER_STREAMS, NBEAMS

!  Input optical depths after delta-M scaling and Chapman function

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT    ( MAXLAYERS )

!  User streams

      REAL(kind=dp), INTENT(IN) :: USER_SECANTS ( MAX_USER_STREAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)       ::  LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

!  Outputs
!  -------

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: T_DELT_USERM ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp), INTENT(OUT) :: ITRANS_USERM ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

!  local variables
!  ---------------

      INTEGER       :: N, UM, IB, LUM
      REAL(kind=dp) :: SPHER, USIB, ITNIB, DVN
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

!  Local user index. !@@ For the Observational Geometry option.

      LUM = 1

!  Transmittance factors for average secant streams
!  ================================================

      IF (DO_SOLAR_SOURCES) THEN

        DO IB = 1, NBEAMS
          DO N = 1, NLAYERS
            IF ( N .GT. LAYER_PIS_CUTOFF(IB) ) THEN
              T_DELT_MUBAR(N,IB) = zero
            ELSE
              SPHER = DELTAU_VERT(N) * AVERAGE_SECANT(N,IB)
              IF ( SPHER .GT. MAX_TAU_PATH ) THEN
                T_DELT_MUBAR(N,IB) = zero
              ELSE
                T_DELT_MUBAR(N,IB) = EXP ( - SPHER )
              ENDIF
            ENDIF
          ENDDO
        ENDDO

      ENDIF

!  Transmittances for User Streams
!  ===============================

!  Return if no post processing

      IF ( .not. DO_POSTPROCESSING ) RETURN

!  Initial transmittances divided by user streams
!  !@@ Option for Observational Geometry, 12/21/12.

      IF (DO_SOLAR_SOURCES) THEN

        IF ( DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            USIB = USER_SECANTS(IB)
            DO N = 1, NLAYERS
              ITRANS_USERM(LUM,N,IB) = INITIAL_TRANS(N,IB) * USIB
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO N = 1, NLAYERS
              ITNIB = INITIAL_TRANS(N,IB)
              DO UM = 1, N_USER_STREAMS
                ITRANS_USERM(UM,N,IB) = ITNIB * USER_SECANTS(UM)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Whole Layer transmittances

      DO UM = 1, N_USER_STREAMS
        DO N = 1, NLAYERS
          DVN = DELTAU_VERT(N)
          SPHER = DVN * USER_SECANTS(UM)
          IF ( SPHER.GT.MAX_TAU_PATH ) THEN
            T_DELT_USERM(UM,N) = zero
          ELSE
            T_DELT_USERM(UM,N) = EXP ( - SPHER )
          ENDIF
        ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_PREPTRANS

!

SUBROUTINE TWOSTREAM_DIRECTBEAM & 
          ( MAXBEAMS,                                   & ! Dimensions
            DO_BRDF_SURFACE,                            & ! Input
            DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,        & ! input @@ 2p3
            NBEAMS, FOURIER, FLUX_FACTOR, X0, PI4,      & ! Input
            DELTA_FACTOR, ALBEDO, BRDF_F_0M,            & ! Input 
            SLTERM_ISOTROPIC, SLTERM_F_0M,              & ! input @@ 2p3
            TRANS_SOLAR_BEAM,                           & ! Input
            ATMOS_ATTN, DIRECT_BEAM )                     ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXBEAMS

!  Surface Control

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  !@@ Version 2p3, 1/23/14. Surface leaving control

      LOGICAL, INTENT(IN)  :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  Numbers

      INTEGER, INTENT(IN)        :: NBEAMS

!  Fourier component

      INTEGER      , INTENT(IN)  :: FOURIER

!  4*pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  FLux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar beams, cosines

      REAL(kind=dp), INTENT(IN)  :: X0(NBEAMS)

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTOR
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0M  ( MAXBEAMS )

!  Do not need this, MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0M ( MAX_USER_STREAMS, MAXBEAMS )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!    ** Isotropic Surface leaving term (if flag set)
!    ** Fourier components of Surface-leaving terms:
!          Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) :: SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: SLTERM_F_0M ( MAXBEAMS )

!  Solar beam attenuations

      REAL(kind=dp), INTENT(IN)  :: TRANS_SOLAR_BEAM        ( MAXBEAMS )

!  output arguments
!  ----------------

!  Atmospheric attenuation

      REAL(kind=dp), INTENT(OUT) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solution, do not need USER_DIRECT_BEAM value

      REAL(kind=dp), INTENT(OUT) :: DIRECT_BEAM      ( MAXBEAMS )
!      REAL(kind=dp), INTENT(OUT) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(kind=dp) :: X0_FLUX, ATTN, HELP
      INTEGER       :: IB

!  Attenuation of solar beam
!  -------------------------

      DO IB = 1, NBEAMS

!  There should be no flux factor here.
!    Bug fixed 18 November 2005. Earlier Italian Job!!
!    Flux Factor put back, 1 March 2007. Using 1 / pi4

        X0_FLUX        = 4.0_dp * X0(IB) * FLUX_FACTOR / DELTA_FACTOR / PI4
        ATTN           = X0_FLUX * TRANS_SOLAR_BEAM(IB)
        ATMOS_ATTN(IB) = ATTN

!  Total contributions, BRDF or Lambertian

        IF ( DO_BRDF_SURFACE ) THEN
          DIRECT_BEAM(IB) = ATTN * BRDF_F_0M(IB)
        ELSE
          DIRECT_BEAM(IB) = ATTN * ALBEDO
        ENDIF

!  New Surface-Leaving stuff 1/23/14, Version 2.4
!    Normalized to Flux-factor / DELTA_Factor
!    Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

        IF ( DO_SURFACE_LEAVING ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTOR
          IF ( DO_SL_ISOTROPIC .and. FOURIER .EQ. 0 ) THEN
            DIRECT_BEAM(IB) = DIRECT_BEAM(IB) + SLTERM_ISOTROPIC(IB) * HELP
          ELSE
            DIRECT_BEAM(IB) = DIRECT_BEAM(IB) + SLTERM_F_0M(IB) * HELP
          ENDIF
        ENDIF

!  end direct beam calculation

      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_DIRECTBEAM

!

SUBROUTINE TWOSTREAM_EMULTMASTER &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                  & ! Dimensions
             DO_UPWELLING, DO_DNWELLING, DO_USER_OBSGEOMS, NLAYERS, NBEAMS, & ! Input
             N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, TAYLOR_SMALL, & ! Input
             USER_SECANTS, DELTAU_VERT, T_DELT_MUBAR, T_DELT_USERM,  & ! Input
             LAYER_PIS_CUTOFF, ITRANS_USERM, ITRANS_USERM_OBS,AVERAGE_SECANT,         & ! Input
             SIGMA_M, SIGMA_P, SIGMA_M_OBS, SIGMA_P_OBS,             & ! Output
             EMULT_UP, EMULT_DN, EMULT_UP_OBS, EMULT_DN_OBS )          ! Output

!  Version 2.4 Overhaul-----
!     Rob  Fix 8/15/14  - Small numbers analysis using Taylor parameters
!     Rob  Fix 8/15/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice choices
!     Rob  Fix 8/15/14  - Compact code in a single subroutine

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Prepare multipliers for the Beam source terms

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)       ::  MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS

!  Control

      LOGICAL, INTENT(IN)       :: DO_UPWELLING, DO_DNWELLING, DO_USER_OBSGEOMS

!  Numbers

      INTEGER, INTENT(IN)       ::  NLAYERS, NBEAMS

!  Post-processing control mask

      INTEGER, INTENT(IN)       :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Version 2p4 - Taylor series control
      
      INTEGER, intent(in)       :: TAYLOR_ORDER
      REAL(kind=dp), intent(in) :: TAYLOR_SMALL

!  User streams

      REAL(kind=dp), INTENT(IN) :: USER_SECANTS ( MAX_USER_STREAMS )

!  Layer optical thickness

      REAL(kind=dp), INTENT(IN) :: DELTAU_VERT ( MAXLAYERS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp), INTENT(IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(IN) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)       ::  LAYER_PIS_CUTOFF(MAXBEAMS)

!  Average-secant and initial tramsittance factors for solar beams.

      REAL(kind=dp), INTENT(IN) :: ITRANS_USERM    ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: ITRANS_USERM_OBS( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT  ( MAXLAYERS, MAXBEAMS )

!  Output = Global multipliers
!  ===========================

!  Coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(OUT) :: SIGMA_M(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: SIGMA_P(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: SIGMA_M_OBS(MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: SIGMA_P_OBS(MAXLAYERS,MAXBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: EMULT_UP_OBS (MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: EMULT_DN_OBS (MAXLAYERS,MAXBEAMS)

!  Local variables
!  ---------------

!  Other variables

      INTEGER         :: N, UM, IB, LUM
      LOGICAL         :: EMULT_HOPRULE
      REAL(kind=dp)   :: WDEL, WUDEL, UDEL, ITUDEL
      REAL(kind=dp)   :: DIFF, SB, SU, SD, SM
      REAL(kind=dp)   :: SIGMA_PN, SIGMA_MN, DELTAU_VERTN

!  sigma functions (SIGMA_P/N) and upwelling/downwelling external source function multipliers (EMULT_UP/DN)
!  --------------------------------------------------------------------------------------------------------

      DO IB = 1, NBEAMS

         LAYER_PIS_CUTOFFB = LAYER_PIS_CUTOFF(IB)
         USER_SECANTSB = USER_SECANTS(IB)
         DO N = 1, LAYER_PIS_CUTOFFB
            SB = AVERAGE_SECANT(N,IB)
            WDEL = T_DELT_MUBAR(N,IB)
            UDEL = T_DELT_USERM(N,IB)
            ITUDEL = ITRANS_USERM_OBS(N,IB)
            IF (DO_USER_OBSGEOMS) THEN
               IF (DO_UPWELLING) THEN
                  SIGMA_PN = SB + USER_SECANTSB
                  SIGMA_P_OBS(N,IB) = SIGMA_PN
                  WUDEL = WDEL * UDEL
                  SU = ( ONE - WUDEL ) / SIGMA_PN
                  EMULT_UP_OBS(N,IB) = ITUDEL * SU
               ENDIF
               IF (DO_DNWELLING) THEN
                  DIFF = ABS ( USER_SECANTSB - SB )
                  IF ( DIFF .LT. TAYLOR_SMALL ) EMULT_HOPRULE = .TRUE.
                  SIGMA_MN = SB - USER_SECANTSB
                  SIGMA_M_OBS(N,IB) = SIGMA_MN
                  IF ( EMULT_HOPRULE ) THEN
                    CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, SIGMA_MN, DELTAU_VERTN, WDEL, ONE, SD )
                  ELSE
                    SD = ( UDEL - WDEL ) / SIGMA_MN
                  ENDIF
                  EMULT_DN_OBS(N,IB) = ITUDEL * SD
               ENDIF
            ELSE
               DELTAU_VERTN = DELTAU_VERT(N)
               DO LUM = 1, N_PPSTREAMS
                  UM = PPSTREAM_MASKB(LUM,IB)
                  SM = USER_SECANTS(UM)
                  UDEL = T_DELT_USERM(N,UM)
                  ITUDEL = ITRANS_USERM(LUM,N,IB)
                  IF (DO_UPWELLING) THEN
                     SIGMA_PN = SB + SM
                     SIGMA_P(LUM,N,IB) = SIGMA_PN
                     WUDEL = WDEL * UDEL
                     SU = ( ONE - WUDEL ) / SIGMA_PN
                     EMULT_UP(LUM,N,IB) = ITUDEL * SU
                  ENDIF
                  IF (DO_DNWELLING) THEN
                     DIFF = ABS ( SM - SB )
                     IF ( DIFF .LT. TAYLOR_SMALL ) EMULT_HOPRULE = .TRUE.
                     SIGMA_MN = SB - SM
                     SIGMA_M(LUM,N,IB) = SIGMA_MN
                     IF ( EMULT_HOPRULE ) THEN
                       CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, SIGMA_MN, DELTAU_VERTN, WDEL, ONE, SD )
                     ELSE
                       SD = ( UDEL - WDEL ) / SIGMA_MN
                     ENDIF
                     EMULT_DN(LUM,N,IB) = ITUDEL * SD
               ENDDO
            ENDIF
         ENDDO

         DO N = 1, LAYER_PIS_CUTOFFB+1, NLAYERS
            IF (DO_USER_OBSGEOMS) THEN
               IF (DO_UPWELLING) THEN
                  SIGMA_P_OBS(N,IB) = ZERO
                  EMULT_UP_OBS(N,IB) = ZERO
               ENDIF
               IF (DO_DNWELLING) THEN
                  SIGMA_M_OBS(N,IB) = ZERO
                  EMULT_DN_OBS(N,IB) = ZERO
               ENDIF
            ELSE
               DO LUM = 1, N_PPSTREAMS
                  IF (DO_UPWELLING) THEN
                     SIGMA_P(LUM,N,IB) = ZERO
                     EMULT_UP(LUM,N,IB) = ZERO
                  ENDIF
                  IF (DO_DNWELLING) THEN
                     SIGMA_M(LUM,N,IB) = ZERO
                     EMULT_DN(LUM,N,IB) = ZERO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

      ENDDO

!  debug
!      do n = 1, nlayers
!        write(57,'(i3,1pe24.12)')n,EMULT_UP(1,N,1)
!        write(57,'(i3,1pe24.12)')n,T_DELT_MUBAR(N,1)
!      enddo

!  debug
!      do n = 1, nlayers
!c        write(57,'(i3,1pe24.12)')n,EMULT_DN(1,N,1)
!c      enddo

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_EMULTMASTER

!

SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE &
           ( MAXLAYERS,                                  & ! Dimensions
             NLAYERS, SZA_GEOM_TRUE, REARTH, HEIGHTS,    & ! Input
             CHAPMAN_FACTORS )                             ! In/Out

      implicit none

!  precision and parameters

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

      REAL(kind=dp), INTENT(INOUT)  :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS)

!  Local variables
!  ===============

!  local height array

      REAL(kind=dp) :: DELZ(MAXLAYERS)

!  help variables

      INTEGER       :: N, K, IB
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

!  finish main layer loop

      ENDDO

!  end of routine

      RETURN
END SUBROUTINE TWOSTREAM_BEAM_GEOMETRY_PREPARE

end module twostream_miscsetups_m

