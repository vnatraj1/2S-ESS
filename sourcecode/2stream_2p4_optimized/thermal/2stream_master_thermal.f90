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
! #            TWOSTREAM_MASTER_THERMAL (top-level master)      #
! #            TWOSTREAM_FOURIER_MASTER_THERMAL                 #
! #                                                             #
! ###############################################################

module twostream_master_thermal_m

Use twostream_miscsetups_m
Use twostream_thermalsup_m
Use twostream_writemodules_thermal_m
Use twostream_solutions_m
Use twostream_solutions_thermal_m
Use twostream_bvproblem_m
Use twostream_intensity_thermal_m
Use twostream_fluxes_m

!Use twostream_geometry_m
!Use twostream_geometry_obs_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_MASTER_THERMAL &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAX_USER_STREAMS,        & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,                & ! Inputs     !@@ 2p2
          DO_POSTPROCESSING, DO_INCLUDE_MVOUT,                       & ! Inputs     !@@ 2p3
          DO_SURFACE_EMISSION, DO_D2S_SCALING, DO_BRDF_SURFACE,      & ! Inputs     !@@ 2p1
          DO_PENTADIAG_INVERSE, BVPSCALEFACTOR,                      & ! Input !@@ 2p3 6/25/14, 8/15/14 
          TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,                       & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_STREAMS,             & ! Inputs     !@@ 2p1
          SURFACE_FACTOR, DELTA_FACTOR,                              & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,       & ! Inputs
          THERMAL_BB_INPUT, ALBEDO, BRDF_F, UBRDF_F,                 & ! Inputs
          EMISSIVITY, SURFBB, PI4, USER_STREAMS, USER_SECANTS, PXSQ, & ! Inputs (thermal,geometry)
          INTENSITY_TOA, INTENSITY_BOA,                              & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,                                  & ! Outputs !@@ 2p2
          FLUXES_TOA_THERMAL, FLUXES_BOA_THERMAL,                    & ! Outputs !@@ 2p3 (Fluxes)
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,     & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )         ! Exception handling

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS

!  @@ Rob Spurr, 15 August 2014, Version 2.4, MAXTHREADS dimension removed

      INTEGER, INTENT(IN)        :: MAXMESSAGES, MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT     ! @@ 2p2

!  Post-processing flag (new for Version 2p3)

      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Flux output flag

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_MVOUT

!  ** New **. October 2011, Sources control

      LOGICAL, INTENT(IN)        :: DO_SURFACE_EMISSION

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  BVP control --- New 6/25/14, Version 2.3 and higher
!  * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!  * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!  * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Version 2.4, August 2014. Order of Taylor series (N) with Smallness number ( EPS)
!                            (Series including terms up to EPS^N)
      
      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!  Rob, introduced 14 May 2015, Version 2.4, following 2p3 implementation (2014)
!  Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Number of user angles

      INTEGER, INTENT(IN)        :: N_USER_STREAMS                     !@@ 2p1

!  Surface reflectance factors

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTOR

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Surface properties

      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  BRDF fourier component 0
!  incident quadrature stream, reflected quadrature stream
!  incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( MAX_USER_STREAMS )

!  Surface thermal sources

      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Constants

      REAL(kind=dp), INTENT(IN) :: PI4

!  Cosines and sines

      REAL(kind=dp), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN) :: USER_SECANTS ( MAX_USER_STREAMS )

!  Polynomial expansion

      REAL(kind=dp), INTENT(IN) :: PXSQ

!  Output
!  ------

!  Radiance Results

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_TOA(MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_BOA(MAX_USER_STREAMS)

!  output solutions at ALL levels
!  ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_UP (0:MAXLAYERS,MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_DN (0:MAXLAYERS,MAX_USER_STREAMS)

!  Flux output
!  ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(OUT) :: FLUXES_TOA_THERMAL(2)
     REAL(kind=dp), INTENT(OUT) :: FLUXES_BOA_THERMAL(2)

!  Exception handling
!  ------------------

!    1. Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!    2. Execution message and 2 Traces

      INTEGER      , INTENT(OUT) :: STATUS_EXECUTION
      CHARACTER*100, INTENT(OUT) :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!  Local definitions
!  =================

!  Local Atmospheric Optical properties
!  ------------------------------------

!  After application of delta-m scaling

      REAL(kind=dp) :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp) :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp) :: ASYMM_TOTAL(MAXLAYERS)

!  Miscsetup operations
!  ====================

!  Transmittance Setups
!  --------------------

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Thermal help variables

      REAL(kind=dp) :: THERMCOEFFS ( 2, MAXLAYERS )

!  Other local variables
!  =====================

!  Local error handling

      INTEGER          :: STATUS_SUB
      INTEGER          :: N
      REAL(kind=dp)    :: OMFAC, M1FAC, GDIFF

!  Singularity buster output

      LOGICAL          :: SBUST(6)

!  Test variables

      LOGICAL          :: DO_DEBUG_INPUT=.FALSE.
      !LOGICAL          :: DO_DEBUG_INPUT=.TRUE.

!  Other local variables

      REAL(kind=dp)    :: ASYMM_TOTALN, D2S_SCALINGN, DELTAU_INPUTN, DELTAU_VERTN
      REAL(kind=dp)    :: OMEGA_INPUTN, OMEGA_TOTALN

!real e1, e2, ftime, ctime

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0

!  Execution status and message/traces

      STATUS_EXECUTION  = 0
      E_MESSAGE = ' '
      E_TRACE_1 = ' '
      E_TRACE_2 = ' '

!  TWOSTREAM input debug
!  ---------------------

      IF (DO_DEBUG_INPUT) THEN
        CALL TWOSTREAM_DEBUG_INPUT_MASTER_THERMAL()
      END IF

!  Input checking
!  ==============

!  Check input optical values (IOPs in the atmosphere)

!call cpu_time(e1)

     CALL TWOSTREAM_CHECK_INPUTS_OPTICAL &
       ( MAXLAYERS, MAXMESSAGES, NLAYERS,                  & ! input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Get derived inputs
!  ==================

!  Set local atmospheric optical properties (Apply delta 2s scaling)
!  Just copy inputs, if not required

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if necessary

      SBUST = .false.
      IF ( DO_D2S_SCALING ) THEN  
         DO N = 1, NLAYERS
            OMEGA_INPUTN = OMEGA_INPUT(N)
            D2S_SCALINGN = D2S_SCALING(N)
            DELTAU_INPUTN = DELTAU_INPUT(N)
            OMFAC = one - OMEGA_INPUTN * D2S_SCALINGN
            M1FAC = one - D2S_SCALINGN
            GDIFF = ASYMM_INPUT(N) - D2S_SCALINGN
            DELTAU_VERTN = OMFAC * DELTAU_INPUTN
            DELTAU_VERT(N) = DELTAU_VERTN
            OMEGA_TOTALN = M1FAC * OMEGA_INPUTN / OMFAC
            !Singularity buster for single scatter albedo
            IF (OMEGA_TOTALN > 0.999999999D0) THEN
               OMEGA_TOTAL(N) = 0.999999999D0
               SBUST(1) = .true.
            ELSE IF (OMEGA_TOTALN < 1.0D-9) THEN
               OMEGA_TOTAL(N) = 1.0D-9
               SBUST(2) = .true.
            ELSE
               OMEGA_TOTAL(N) = OMEGA_TOTALN
            END IF   
            ASYMM_TOTALN = GDIFF / M1FAC
            !Singularity buster for asymmetry parameter
            IF (ASYMM_TOTALN > 0.999999999D0) THEN
               ASYMM_TOTAL(N) = 0.999999999D0
               SBUST(3) = .true.
            ELSE IF (ASYMM_TOTALN < -0.999999999D0) THEN
               ASYMM_TOTAL(N) = -0.999999999D0  
               SBUST(4) = .true.
            ELSE IF ((ASYMM_TOTALN >= ZERO) .AND. &
                     (ASYMM_TOTALN < 1.0D-9)) THEN
               ASYMM_TOTAL(N) = 1.0D-9
               SBUST(5) = .true.
            ELSE IF ((ASYMM_TOTALN < ZERO) .AND. &
                     (ASYMM_TOTALN > -1.0D-9)) THEN
               ASYMM_TOTAL(N) = -1.0D-9
               SBUST(6) = .true.
            ELSE
               ASYMM_TOTAL(N) = ASYMM_TOTALN
            END IF
            !WRITE(*,*)
            !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
            !DO I=1,6
            !   WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
            !ENDDO
         ENDDO
      ELSE
         DO N = 1, NLAYERS
            DELTAU_VERT(N) = DELTAU_INPUT(N)
            OMEGA_TOTAL(N) = OMEGA_INPUT(N)
            ASYMM_TOTAL(N) = ASYMM_INPUT(N)
         ENDDO
      ENDIF

!call cpu_time(e2)

!write(0,*) 'input manipulation time = ', e2-e1

!  SETUP OPERATIONS (moved from Fourier, Version 2p3, 15 August 2014)
!  ==================================================================

!  MISCSETUPS (ONE subroutine, which does the following):
!  transmittance setup (T_DELT_USERM: postprocessing only)

!call cpu_time(e1)

      IF (DO_POSTPROCESSING) THEN
         CALL TWOSTREAM_PREPTRANS &
            ( MAXLAYERS, MAX_USER_STREAMS,  & ! Dimensions
              NLAYERS, N_USER_STREAMS,      & ! Input
              DELTAU_VERT, USER_SECANTS,    & ! Input
              T_DELT_USERM )                  ! Output
      ENDIF

!call cpu_time(e2)

!write(0,*) 'preptrans time = ', e2-e1

!  ** New, October 2011 **. Thermal setup

!call cpu_time(e1)

      CALL TWOSTREAM_THERMALSETUP &
         ( MAXLAYERS, NLAYERS,               & ! Input
           DELTAU_VERT, THERMAL_BB_INPUT,    & ! Input
           THERMCOEFFS )                       ! input/Output

!call cpu_time(e2)

!write(0,*) 'thermal set up time = ', e2-e1

!  Fourier loop
!  ============

!ftime = 0.0
!ctime = 0.0

!  Main call to TWOSTREAM Fourier module
!  -------------------------------------

!call cpu_time(e1)

      CALL TWOSTREAM_FOURIER_MASTER_THERMAL &
           ( MAXLAYERS, MAXTOTAL, MAX_USER_STREAMS, DO_2S_LEVELOUT,            & ! Dimensions + flag
             DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                      & ! Input flags control
             DO_POSTPROCESSING, DO_PENTADIAG_INVERSE,                          & ! Inputs
             DO_SURFACE_EMISSION, DO_INCLUDE_MVOUT,                            & ! Inputs
             BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,              & ! Input !@@ 2p3 6/25/14, 8/15/14
             NLAYERS, NTOTAL, N_USER_STREAMS,                                  & ! Input integer control
             PI4, SURFACE_FACTOR, DELTA_FACTOR,                                & ! Input real control
             STREAM_VALUE, USER_STREAMS, USER_SECANTS, PXSQ,                   & ! Inputs
             ALBEDO, BRDF_F, UBRDF_F, SURFBB, EMISSIVITY,                      & ! Input real surface
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS, T_DELT_USERM, & ! Input real optical
             INTENSITY_TOA, INTENSITY_BOA, RADLEVEL_UP, RADLEVEL_DN,           & ! Output
             FLUXES_TOA_THERMAL, FLUXES_BOA_THERMAL,                           & ! Output (modified 2p3, Fluxes)
             STATUS_SUB, E_MESSAGE, E_TRACE_1 )                                  ! Outputs

!  Exception handling

         IF ( STATUS_SUB .NE. 0 ) THEN
            STATUS_EXECUTION = 1
            E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # 0'
            RETURN
         ENDIF

!call cpu_time(e2)

!ftime = ftime+e2-e1

!write(0,*) 'fourier subroutine time = ', ftime

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER_THERMAL()

      CALL TWOSTREAM_WRITE_STD_INPUT_THERMAL ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAX_USER_STREAMS,          &
        DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,                  &
        DO_INCLUDE_MVOUT, DO_POSTPROCESSING, DO_PENTADIAG_INVERSE,   &
        DO_SURFACE_EMISSION, DO_D2S_SCALING, DO_BRDF_SURFACE,        &
        BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,         &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_STREAMS, USER_STREAMS, &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,         &
        THERMAL_BB_INPUT, ALBEDO, SURFBB )

      IF (DO_BRDF_SURFACE) THEN
        CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT_THERMAL (   &
          MAX_USER_STREAMS, N_USER_STREAMS, &
          BRDF_F, UBRDF_F, EMISSIVITY )
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER_THERMAL

END SUBROUTINE TWOSTREAM_MASTER_THERMAL

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER_THERMAL &
        ( MAXLAYERS, MAXTOTAL, MAX_USER_STREAMS, DO_2S_LEVELOUT,            & ! Dimensions + flag
          DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE,                      & ! Input flags control
          DO_POSTPROCESSING, DO_PENTADIAG_INVERSE,                          & ! Input flags !@@ 2p3 6/25/14
          DO_INCLUDE_SURFEMISS, DO_INCLUDE_MVOUT,                           & ! Inputs
          BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,              & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, N_USER_STREAMS,                                  & ! Input integer control
          PI4, SURFACE_FACTOR, DELTA_FACTOR,                                & ! Input real control
          STREAM_VALUE, USER_STREAMS, USER_SECANTS, PXSQ,                   & ! Inputs
          ALBEDO, BRDF_F, UBRDF_F, SURFBB, EMISSIVITY,                      & ! Input real surface
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS, T_DELT_USERM, & ! Input real optical
          INTENSITY_TOA, INTENSITY_BOA, RADLEVEL_UP, RADLEVEL_DN,           & ! Output
          FLUXES_TOA_THERMAL, FLUXES_BOA_THERMAL,                           & ! Output
          STATUS, MESSAGE, TRACE )                                            ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions :
!      MAXTOTAL  = 2 * MAXLAYERS
     
      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)  :: MAX_USER_STREAMS

!  Flags

      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE

!  Post processing flag

      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!  Other flags

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_SURFEMISS
      LOGICAL, INTENT(IN)  :: DO_INCLUDE_MVOUT

!  BVP control --- New 6/25/14, Version 2.3 and higher
!  * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!  * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!  * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)        :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!  Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!  Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: N_USER_STREAMS

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Surface reflectance factors

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR
      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTOR

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Polynomial expansion

      REAL(kind=dp), INTENT(IN) :: PXSQ

!  Surface variables
!  ------------------

!  Lambertian Albedo

     REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident quadrature stream, reflected quadrature stream
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F    ( MAX_USER_STREAMS )

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_TOA (MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_BOA (MAX_USER_STREAMS)

!  Flux output

     REAL(kind=dp), INTENT(OUT) :: FLUXES_TOA_THERMAL(2)
     REAL(kind=dp), INTENT(OUT) :: FLUXES_BOA_THERMAL(2)

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_UP (0:MAXLAYERS,MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_DN (0:MAXLAYERS,MAX_USER_STREAMS)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Miscsetups Arrays
!  =================

!  Thermal help variables
!  ----------------------

      REAL(kind=dp), INTENT(IN) :: THERMCOEFFS ( 2, MAXLAYERS )

!  Transmittance for user-defined stream angles
!  --------------------------------------------

      REAL(kind=dp), INTENT(IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Local Arrays for Use in Subroutines
!  ===================================

!  Multiplier arrays (Homogeneous solutions)
!  -----------------

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Solutions to the homogeneous RT equations
!  -----------------------------------------

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANS(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANSNL

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp) :: NORM_SAVED(MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!  EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp) :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Boundary Value Problem
!  Original and Elimination matrices (Pentadiagonal, 2x2)
!  Split into multiple vectors

      REAL(kind=dp) :: SELM   (2,2)
      REAL(kind=dp) :: ELM1   (MAXTOTAL-1)
      REAL(kind=dp) :: ELM2   (MAXTOTAL-2)
      REAL(kind=dp) :: ELM3   (MAXTOTAL)
      REAL(kind=dp) :: ELM4   (3:MAXTOTAL)
      REAL(kind=dp) :: MAT1   (3:MAXTOTAL)
      REAL(kind=dp) :: MAT22

!  BVP solution
!  ------------

!  Solution constants of integration and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC1NL, MCON_XVEC1NL, LCON_XVEC21, MCON_XVEC21

!  Thermal solutions
!  -----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp) :: T_C_MINUS (0:2,MAXLAYERS)
      REAL(kind=dp) :: T_C_PLUS  (0:2,MAXLAYERS)

!  Thermal solution at layer boundaries and related quantities

      REAL(kind=dp) :: T_WUPPER (2,MAXLAYERS)
      REAL(kind=dp) :: T_WLOWER (2,MAXLAYERS)
      REAL(kind=dp) :: WUPPER21
      REAL(kind=dp) :: WLOWER1NL


!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Local help variables
!  --------------------

!  local inclusion flags. ** New October 2011 **, thermal flags
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      LOGICAL :: DO_INCLUDE_SURFACE

!  Other local variables

      REAL(kind=dp) :: EIGENTRANS1

!  Error tracing

      INTEGER       :: STATUS_SUB

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Surface flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
        IF ( ALBEDO .NE. ZERO) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Get Discrete ordinate solutions for all layers
!  Version 2p4. Green's function output = NORM_SAVED

       CALL TWOSTREAM_HOM_SOLUTION_THERMAL &
          ( MAXLAYERS, NLAYERS, STREAM_VALUE, PXSQ,    & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,     & ! Input
            EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED )   ! Output

!  Get Post-processing ("user") solutions for this layer
!  !@@ 2p3. 11/5/13. Post-processing control

       IF ( DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_HOM_USERSOLUTION_THERMAL &
             ( MAXLAYERS, MAX_USER_STREAMS,                  & ! Dimensions
               NLAYERS, N_USER_STREAMS, STREAM_VALUE,        & ! Input
               USER_STREAMS, XPOS, OMEGA_TOTAL, ASYMM_TOTAL, & ! Input
               U_XPOS, U_XNEG )                                ! Output
       ENDIF

!  Prepare homogeneous solution multipliers
!  !@@ 2p3. 11/5/13. Post-processing control
!  !@@ 2p4. 8/15/14. User secants, Taylor-series control

      IF ( DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_HMULT_MASTER &
            ( MAXLAYERS, MAX_USER_STREAMS,             & ! Dimensions
              TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, & ! Inputs 
              NLAYERS, N_USER_STREAMS, USER_SECANTS,   & ! Input
              EIGENVALUE, EIGENTRANS, T_DELT_USERM,    & ! Input
              HMULT_1, HMULT_2 )                         ! Output
      ENDIF

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)
!  Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
          ( MAXLAYERS, MAXTOTAL, BVPSCALEFACTOR, DO_PENTADIAG_INVERSE, & ! Dimensions
            DO_INCLUDE_SURFACE, NLAYERS, NTOTAL,                       & ! Input
            DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F,           & ! Input
            XPOS, EIGENTRANS, STREAM_VALUE,                            & ! Input
            MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM,                 & ! Output
            STATUS_SUB, MESSAGE )                                        ! Output

!  Exception handling for Pentadiagonal Matrix setup

      IF ( STATUS_SUB .NE. 0 ) THEN
         TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
         STATUS = 1 ; RETURN
      ENDIF

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. Greens function solution for Regular code
!   !@@ 2p4. 5/14/15. Thermal Cutoff variable introduced.

      CALL TWOSTREAM_THERMALGFSOLUTION &
         ( MAXLAYERS, NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMCOEFFS,  & ! Inputs
           TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED,          & ! Input
           T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )         ! Outputs
      IF ( DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_THERMALSTERMS_THERMAL &
            ( MAXLAYERS, MAX_USER_STREAMS,            & ! Dimensions
              DO_UPWELLING, DO_DNWELLING,             & ! Input
              NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Input
              TCUTOFF, T_DELT_USERM, DELTAU_VERT,     & ! Input
              U_XPOS, U_XNEG, HMULT_1, HMULT_2,       & ! Inputs
              T_C_PLUS, T_C_MINUS, TTERM_SAVE,        & ! Inputs
              LAYER_TSUP_UP, LAYER_TSUP_DN  )           ! Output
      ENDIF

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Solve boundary value problem (Pentadiagonal solution)

         CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
            ( MAXLAYERS, MAXTOTAL,                                & ! Dimensions
              BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,               & ! BVP control, 6/24/14
              DO_INCLUDE_SURFACE, .FALSE.,                        & ! Input
              DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,              & ! Input
              NLAYERS, NTOTAL,                                    & ! Input
              SURFACE_FACTOR, ALBEDO, BRDF_F, EMISSIVITY, SURFBB, & ! Input
              ZERO, T_WUPPER, T_WLOWER, STREAM_VALUE,             & ! Input
              MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM,          & ! Input
              LCON, MCON )                                          ! Output

         IF (( DO_UPWELLING .and. DO_POSTPROCESSING ) .OR. &
             ( DO_DNWELLING .and. DO_INCLUDE_MVOUT )) THEN
            EIGENTRANSNL = EIGENTRANS(NLAYERS) 
            LCON_XVEC1NL = LCON(NLAYERS)*XPOS(1,NLAYERS)
            MCON_XVEC1NL = MCON(NLAYERS)*XPOS(2,NLAYERS)
            WLOWER1NL    = T_WLOWER(1,NLAYERS)
         ENDIF

!  Upwelling, MSMODE only, no Direct Beam inclusion

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_INTENSITY_THERMAL &
               ( MAXLAYERS, MAX_USER_STREAMS,                         & ! Dimensions
                 DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_2S_LEVELOUT, & ! Flags
                 NLAYERS, N_USER_STREAMS,                             & ! Inputs
                 SURFACE_FACTOR, ALBEDO, UBRDF_F,                     & ! Inputs
                 DELTA_FACTOR, STREAM_VALUE,                          & ! Inputs
                 T_DELT_USERM, EIGENTRANSNL,                          & ! Inputs
                 LCON, LCON_XVEC1NL, MCON, MCON_XVEC1NL, WLOWER1NL,   & ! Inputs
                 U_XPOS, U_XNEG, HMULT_1, HMULT_2, LAYER_TSUP_UP,     & ! Inputs
                 INTENSITY_TOA, RADLEVEL_UP )                           ! Outputs
         ENDIF

!  Downwelling, MSMODE only

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_DNUSER_INTENSITY_THERMAL &
               ( MAXLAYERS, MAX_USER_STREAMS,     & ! Dimensions
                 DO_2S_LEVELOUT,                  & ! Flag
                 NLAYERS, N_USER_STREAMS,         & ! Inputs
                 DELTA_FACTOR, T_DELT_USERM,      & ! Inputs
                 LCON, MCON, U_XPOS, U_XNEG,      & ! Inputs
                 HMULT_1, HMULT_2, LAYER_TSUP_DN, & ! Inputs
                 INTENSITY_BOA, RADLEVEL_DN )       ! Outputs
         ENDIF

!  Flux output

         IF ( DO_INCLUDE_MVOUT ) THEN
            IF (DO_UPWELLING) THEN
               EIGENTRANS1 = EIGENTRANS(1) 
               LCON_XVEC21 = LCON(1)*XPOS(2,1)
               MCON_XVEC21 = MCON(1)*XPOS(1,1)
               WUPPER21    = T_WUPPER(2,1)
            ENDIF
            CALL TWOSTREAM_FLUXES_THERMAL &
               ( DO_UPWELLING, DO_DNWELLING,                          & ! Flags
                 PI4, STREAM_VALUE, DELTA_FACTOR,                     & ! Inputs
                 LCON_XVEC21, MCON_XVEC21, EIGENTRANS1, WUPPER21,     & ! Inputs
                 LCON_XVEC1NL, MCON_XVEC1NL, EIGENTRANSNL, WLOWER1NL, & ! Inputs
                 FLUXES_TOA_THERMAL, FLUXES_BOA_THERMAL )               ! Outputs
         ENDIF

!  Finish Thermal only.

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER_THERMAL

end module twostream_master_thermal_m
