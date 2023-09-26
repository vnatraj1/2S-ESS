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
! #            TWOSTREAM_INPUTS_SOLAR_OBS                       #
! #            TWOSTREAM_CHECK_INPUT_DIMS_SOLAR_OBS             #
! #            TWOSTREAM_CHECK_INPUTS_BASIC_SOLAR_OBS           #
! #                                                             #
! ###############################################################

module twostream_inputs_solar_obs_m

USE twostream_brdf_supplement_solar_obs_m
USE twostream_sleave_supplement_m
USE twostream_geometry_m
USE twostream_geometry_solar_obs_m


      PRIVATE
      PUBLIC :: TWOSTREAM_INPUTS_SOLAR_OBS, TWOSTREAM_CHECK_INPUT_DIMS_SOLAR_OBS, &
                TWOSTREAM_CHECK_INPUTS_BASIC_OBS

      CONTAINS

SUBROUTINE TWOSTREAM_INPUTS_SOLAR_OBS &
        ( MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, MAXMESSAGES,              & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                    & ! Inputs
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                               & ! Inputs
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,    & ! Inputs
          EARTH_RADIUS, HEIGHT_GRID,                                        & ! Inputs
          MAXSTREAMS_BRDF, NSTREAMS_BRDF, DO_BRDF_SURFACE,                  & ! BRDF Inputs
          DO_SURFACE_LEAVING, DO_FLUORESCENCE, FL_DO_DataGaussian, DO_ISOTROPIC, & ! SLEAVE Inputs
          FL_Latitude, FL_Longitude, FL_Epoch, FL_InputGAUSSIANS,           & ! SLEAVE Inputs
          N_FOURIERS, DEG_TO_RAD, PI4, DO_POSTPROCESSING, DO_INCLUDE_MVOUT, & ! Outputs
          AVERAGE_SECANT_PP, CHAPMAN_FACTORS, X0,                           & ! Outputs
          USER_STREAMS, USER_SECANTS, AZMFAC,                               & ! Outputs
          PX11, PXSQ, PX0X, ULP,                                            & ! Outputs
          SURFACE_FACTOR, DELTA_FACTOR,                                     & ! Outputs
          PIE, STREAM_SINE, SX0, USER_SINES, NBRDF_HALF,                    & ! BRDF Outputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF,                                 & ! BRDF Outputs
          FL_GAUSSIANS, Fs755,                                              & ! SLEAVE Outputs
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS )             ! Exception handling

      IMPLICIT NONE

!  Precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: one = 1.0_dp

!  Dimensions

      INTEGER, INTENT(IN)          :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)          :: MAX_USER_OBSGEOMS, MAXMESSAGES

!  Directional Flags

      LOGICAL, INTENT(IN)          :: DO_UPWELLING, DO_DNWELLING

!  Plane parallel flag

      LOGICAL, INTENT(IN)          :: DO_PLANE_PARALLEL

!  Flux option flags

      LOGICAL, INTENT(IN)          :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)          :: DO_ADDITIONAL_MVOUT

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)          :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)    :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]

      INTEGER, INTENT(IN)          :: N_USER_OBSGEOMS
      REAL(kind=dp), INTENT(IN)    :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS

!  Number of Fourier terms

      INTEGER, INTENT(OUT)         :: N_FOURIERS

!  Constants

      REAL(kind=dp), INTENT(OUT)   :: DEG_TO_RAD, PI4

!  Post-processing flag

      LOGICAL, INTENT(OUT)         :: DO_POSTPROCESSING

!  Flux flag

      LOGICAL, INTENT(OUT)         :: DO_INCLUDE_MVOUT(0:1)

!  Average secant (plane-parallel case)

      REAL(kind=dp), INTENT(OUT)   :: AVERAGE_SECANT_PP ( MAX_USER_OBSGEOMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(OUT)   :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAX_USER_OBSGEOMS )

!  Cosines and sines

      REAL(kind=dp), INTENT(OUT)   :: X0 ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(OUT)   :: USER_STREAMS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(OUT)   :: USER_SECANTS ( MAX_USER_OBSGEOMS )

!  Azimuth cosine factor

      REAL(kind=dp), INTENT(OUT)   :: AZMFAC ( MAX_USER_OBSGEOMS )

!  Polynomial expansions

      REAL(kind=dp), INTENT(OUT )  :: PX11, PXSQ(0:1)
      REAL(kind=dp), INTENT(OUT)   :: PX0X ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp), INTENT(OUT)   :: ULP ( MAX_USER_OBSGEOMS ) 

!  Surface reflectance factors

      REAL(kind=dp), INTENT(OUT)   :: SURFACE_FACTOR(0:1)
      REAL(kind=dp), INTENT(OUT)   :: DELTA_FACTOR(0:1)

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

!  Solar zenith sine

      REAL(kind=dp), INTENT(OUT)   :: SX0(MAX_USER_OBSGEOMS)

!  Viewing zenith sine

      REAL(kind=dp), INTENT(OUT)   :: USER_SINES(MAX_USER_OBSGEOMS)

!  BRDF azimuth quadrature streams

      INTEGER, INTENT(OUT)         :: NBRDF_HALF
      REAL(kind=dp), INTENT(OUT)   :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), INTENT(OUT)   :: A_BRDF  ( MAXSTREAMS_BRDF )

!  SLEAVE variables
!  ================

!  Surface leaving flag

      LOGICAL, INTENT(IN) :: DO_SURFACE_LEAVING

!  Isotropic flag  

      LOGICAL, INTENT(IN) :: DO_ISOTROPIC

!  Flo flag

      LOGICAL, INTENT(IN) :: DO_FLUORESCENCE

!  Fluorescence variables
!  ----------------------

!  Input Latitude/Longitude in [degs]

      REAL(kind=dp), INTENT(IN) :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER, INTENT(IN)       :: FL_Epoch(6)

!  Flag for using Data Gaussians

      LOGICAL, INTENT(IN)       :: FL_DO_DataGaussian

!  Input Gaussians (alternative to Data Gaussians)

      REAL(kind=dp), INTENT(IN) ::  FL_InputGAUSSIANS(3,2)

!  Output Gaussians

      REAL(kind=dp), INTENT(OUT) ::  FL_GAUSSIANS(3,2)

!  Fluorescence at 755 nm

      REAL(kind=dp), INTENT(OUT) :: Fs755(MAX_USER_OBSGEOMS)

!  Exception handling
!  ------------------

!    Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER          :: STATUS_SUB
      INTEGER          :: IB, UM, I, I1
      REAL(kind=dp)    :: X_BRDFI

!  Fluorescence file

      CHARACTER*60     :: Fluofile

!  Fluorescence Gaussian parameters
!  Parameters of the fluorescence Gaussian spectral shape model
!           Gaussian    A (Wm−2 μm−1 sr−1) Lambda(nm) Sigma(nm)
!              1           1.445           736.8        21.2
!              2           0.868           685.2        9.55

      REAL(kind=dp)    :: FL_DataGAUSSIANS(3,2)
      DATA FL_DataGAUSSIANS(1,1) / 1.445d0 /
      DATA FL_DataGAUSSIANS(2,1) / 736.8d0 /
      DATA FL_DataGAUSSIANS(3,1) / 21.2d0  /
      DATA FL_DataGAUSSIANS(1,2) / 0.868d0 /
      DATA FL_DataGAUSSIANS(2,2) / 685.2d0 /
      DATA FL_DataGAUSSIANS(3,2) / 9.55d0  /

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Check input dimensions
!  ----------------------

      CALL TWOSTREAM_CHECK_INPUT_DIMS_SOLAR_OBS &
      ( MAXMESSAGES, MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, & ! Dimensions
        NLAYERS, NTOTAL, N_USER_OBSGEOMS,                    & ! Inputs
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )       ! Outputs

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

      CALL TWOSTREAM_CHECK_INPUTS_BASIC_SOLAR_OBS  &
        ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,             & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,         & ! Inputs
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                    & ! Inputs
          NLAYERS, N_USER_OBSGEOMS, USER_OBSGEOMS,               & ! Inputs
          HEIGHT_GRID, EARTH_RADIUS, DO_POSTPROCESSING,          & ! Input, Input/Output, Output
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )         ! Outputs

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Chapman function calculation
!  ----------------------------

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
         DO IB = 1, N_USER_OBSGEOMS
            CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE &
               ( MAXLAYERS,                    & ! Dimension
                 NLAYERS, USER_OBSGEOMS(IB,1), & ! Inputs
                 EARTH_RADIUS, HEIGHT_GRID,    & ! Inputs
                 CHAPMAN_FACTORS(:,:,IB) )       ! In/Out
         ENDDO
      ENDIF

!  Get derived inputs
!  ==================

!  Solar zenith angle cosine and average secant for PP case

      DO IB = 1, N_USER_OBSGEOMS
         X0(IB) = COS ( USER_OBSGEOMS(IB,1) * DEG_TO_RAD )
         IF (DO_PLANE_PARALLEL) AVERAGE_SECANT_PP(IB) = ONE / X0(IB)   
      ENDDO

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN
         DO I = 1, N_USER_OBSGEOMS
            USER_STREAMS(I) = COS(DEG_TO_RAD * USER_OBSGEOMS(I,2))
            USER_SECANTS(I) = ONE / USER_STREAMS(I)
         ENDDO
      ENDIF

!  Auxiliary geometry quantities

      CALL TWOSTREAM_AUXGEOM_SOLAR_OBS &
         ( MAX_USER_OBSGEOMS,              & ! Dimension
           DO_POSTPROCESSING,              & ! Flag
           N_USER_OBSGEOMS,                & ! Input
           X0, USER_STREAMS, STREAM_VALUE, & ! Inputs
           PX11, PXSQ, PX0X, ULP )           ! Outputs

!  Fourier loop set-up
!  ===================

!  Set Fourier number, Nominally 1 in absence of SS-only flag
!  Zero if no solar sources (Thermal-only run)
!  Set NFOURIERS equal to zero for MVOUT_ONLY

      N_FOURIERS = 1
      IF (  DO_MVOUT_ONLY )         N_FOURIERS = 0

!  Mick fix 1/7/2012 - (test - make this permanent?)

      IF ( (N_USER_OBSGEOMS == 1) .AND. (USER_OBSGEOMS(1,1) < 1.0D-8) ) &
        N_FOURIERS = 0

!  Azimuth cosine factor (Fourier = 1)
!  Not required for Flux-only output

      AZMFAC = 0.0_dp
      IF ( DO_POSTPROCESSING ) THEN
         DO IB = 1, N_USER_OBSGEOMS
            AZMFAC(IB) = COS(DEG_TO_RAD*USER_OBSGEOMS(IB,3))
         ENDDO
      ENDIF

!  Set flags
!  ---------

! Inclusion of mean value calculation
! Control for the Flux calculation

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        DO_INCLUDE_MVOUT(0) = .TRUE.
      ENDIF

!  Surface reflectance factors

      SURFACE_FACTOR(0) = 2.0_dp
      DELTA_FACTOR(0)   = ONE
      SURFACE_FACTOR(1) = ONE
      DELTA_FACTOR(1)   = 2.0_dp

!  BRDF set up

      IF ( DO_BRDF_SURFACE ) THEN

!  Constants

         PIE = DEG_TO_RAD * 180.0_dp
         STREAM_SINE = SQRT(ONE - STREAM_VALUE * STREAM_VALUE)

!  Half number of moments

         NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams. Optionality, added 12/31/12
!  Warning, this should be the BOA angle. OK for the non-refractive case.

         DO IB = 1, N_USER_OBSGEOMS
            SX0(IB) = SQRT(ONE-X0(IB)*X0(IB))
         ENDDO

!  Viewing angles

         DO UM = 1, N_USER_OBSGEOMS
            USER_SINES(UM)   = SQRT(1.0_dp-USER_STREAMS(UM)*USER_STREAMS(UM))
         ENDDO

!  BRDF quadrature

         CALL TWOSTREAM_GAULEG ( 0.0_dp, ONE, X_BRDF, A_BRDF, NBRDF_HALF )
         DO I = 1, NBRDF_HALF
            I1 = I + NBRDF_HALF
            X_BRDFI = X_BRDF(I)
            X_BRDF(I1) = - X_BRDFI
            A_BRDF(I1) =   A_BRDF(I)
         ENDDO
         DO I = 1, NSTREAMS_BRDF
            X_BRDFI = X_BRDF(I)
            X_BRDF(I) = PIE * X_BRDFI
            CX_BRDF(I) = COS (X_BRDFI)
            SX_BRDF(I) = SIN (X_BRDFI)
         ENDDO

      ENDIF

!  End BRDF set up

!---------------
!  SLEAVE set up
!  =============
!---------------

!  Fluorescence
!  ============

      IF ( DO_SURFACE_LEAVING .AND. DO_FLUORESCENCE ) THEN

!  Set FL Gaussians

         IF ( FL_DO_DataGaussian ) then
            FL_GAUSSIANS(1:3,1) = FL_DataGAUSSIANS(1:3,1)
            FL_GAUSSIANS(1:3,2) = FL_DataGAUSSIANS(1:3,2)
         ELSE
            FL_GAUSSIANS(1:3,1) = FL_InputGAUSSIANS(1:3,1)
            FL_GAUSSIANS(1:3,2) = FL_InputGAUSSIANS(1:3,2)
         ENDIF

!  Temporary - Only Isotropic yet

         IF ( .NOT. DO_ISOTROPIC ) &
            Stop 'Non-isotropic not allowed yet if doing fluorescence'

!  F_755 data file

         Fluofile = 'lidort_test/data/fluor_data_2009_fortran.dat'

!  For each Solar zenith angle

         DO IB = 1, N_USER_OBSGEOMS

!  Get the F_755 data from the subroutine

            CALL get_fluorescence_755 &
               ( FL_Latitude, FL_Longitude, FL_Epoch, USER_OBSGEOMS(IB,1), FluoFile, Fs755(IB) )

!  End Beam loop
        
         ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( DO_SURFACE_LEAVING .AND. .NOT. DO_FLUORESCENCE ) THEN

!  Temporary - Only Isotropic yet

         IF ( .NOT. DO_ISOTROPIC ) &
            STOP 'Non-isotropic not allowed yet if not doing fluorescence'

!  PLACEHOLDERS for other Water-leaving options

      ENDIF  

!  Finish

      RETURN

END SUBROUTINE TWOSTREAM_INPUTS_SOLAR_OBS

SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS_SOLAR_OBS &
      ( MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL, N_USER_OBSGEOMS, &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  Max dimension input variables

      INTEGER, INTENT(IN)  ::  MAXMESSAGES
      INTEGER, INTENT(IN)  ::  MAXLAYERS
      INTEGER, INTENT(IN)  ::  MAXTOTAL
      INTEGER, INTENT(IN)  ::  MAX_USER_OBSGEOMS

!  Active dimension input variables

      INTEGER, INTENT(IN)  ::  NLAYERS
      INTEGER, INTENT(IN)  ::  NTOTAL
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

!  2. Geometry dimensions - always checked

      IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Observation Geometries > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_OBSGEOMS dimension'
        STATUS       = 1
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish
      RETURN

END SUBROUTINE TWOSTREAM_CHECK_INPUT_DIMS_SOLAR_OBS

SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC_SOLAR_OBS &
         ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,     & ! Dimensions
           DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, & ! Inputs
           DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,            & ! Inputs
           NLAYERS, N_USER_OBSGEOMS, USER_OBSGEOMS,       & ! Inputs
           HEIGHT_GRID, EARTH_RADIUS, DO_POSTPROCESSING,  & ! Input, Input/Output, Output
           STATUS, NMESSAGES, MESSAGE, ACTION )             ! Outputs, Input/Output

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Inputs
!  ------

!  Dimensions :

      INTEGER, INTENT(IN) :: MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS

!  Flags

      LOGICAL, INTENT(IN) :: DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL

!  Flux output flags

      LOGICAL, INTENT(IN)    :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)    :: DO_ADDITIONAL_MVOUT

!  Number of layers

      INTEGER, INTENT(IN) :: NLAYERS

!  Observational geometry input. [Same as LIDORT]

      INTEGER, INTENT(IN) :: N_USER_OBSGEOMS
      REAL(kind=dp), INTENT(IN) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3)

!  Height and Earth radius

      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )
      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS

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

!  Check MAX number of messages should be at least 12 (covers all errors below as well as those checking 
!  input dimensions )

      IF ( MAXMESSAGES .LT. 12 ) THEN
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

!  Check Earth radius (Chapman function only)
!  ---WARNING. Default value of 6371.0 will be set

      IF ( .NOT. DO_PLANE_PARALLEL ) THEN
        IF ( EARTH_RADIUS .LT. 6320.0D0 .OR. EARTH_RADIUS .GT. 6420.0D0 ) THEN
          NM = NM + 1
          MESSAGE(NM)= 'Bad input: Earth radius outside 6320-6420 km'
          ACTION(NM) = 'Warning: default value of 6371.0 km was set'
          EARTH_RADIUS = 6371.0D0
        ENDIF
      ENDIF

!  Check solar zenith angle input

      LOOP = .TRUE.
      I = 0
      DO WHILE ( LOOP .AND. I .LT. N_USER_OBSGEOMS )
        I = I + 1 
        IF ( USER_OBSGEOMS(I,1) .LT. 0.0D0 .OR. USER_OBSGEOMS(I,1) .GE. 90.0D0 ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM)= 'Bad input: out-of-range beam angle, no. '//C2
          ACTION(NM) = 'Look at SZA input, should be < 90 & > 0'
          LOOP = .FALSE.
          STATUS = 1 
        ENDIF
      ENDDO

!  Check relative azimuths
!  Avoid this section if MVOUT_ONLY

      IF ( .NOT. DO_MVOUT_ONLY ) THEN
        LOOP = .TRUE. ; I = 0
        DO WHILE ( LOOP .AND. I .LT. N_USER_OBSGEOMS )
          I = I + 1
          IF ( USER_OBSGEOMS(I,3) .GT. 360.0D0 .OR. USER_OBSGEOMS(I,3) .LT. 0.0D0 ) THEN
            NM = NM + 1
            WRITE(C2,'(I2)')I
            MESSAGE(NM)='Bad input: out-of-range azimuth angle, no. '//C2
            ACTION(NM) = 'Look at azimuth angle input, range [0,360]'
            LOOP = .FALSE.
            STATUS = 1
          ENDIF
        ENDDO
      ENDIF

!  Check user-defined stream angles (should always be [0,90])
!  Avoid this section if MVOUT_ONLY

      IF ( .NOT. DO_MVOUT_ONLY ) THEN
        LOOP = .TRUE.
        I = 0
        DO WHILE ( LOOP .AND. I .LT. N_USER_OBSGEOMS )
          I = I + 1
          IF ( USER_OBSGEOMS(I,2) .GT. 90.0D0 .or.USER_OBSGEOMS(I,2) .LT. 0.0D0 ) THEN
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
      DO WHILE ( LOOP .AND. I .LT. NLAYERS )
        I = I + 1
        IF ( HEIGHT_GRID(I-1) .LE. HEIGHT_GRID(I) ) THEN
          NM = NM + 1
          WRITE(C2,'(I2)')I
          MESSAGE(NM) = 'Bad input: Height-grid not monotonic decreasing; Layer '//C2
          ACTION(NM) = 'Look at Height-grid input'
          LOOP = .FALSE.
          STATUS = 1
        ENDIF
      ENDDO

!  Set number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_BASIC_SOLAR_OBS

end module twostream_inputs_solar_obs_m
