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
! #            TWOSTREAM_MASTER (top-level master)              #
! #            TWOSTREAM_FOURIER_MASTER                         #
! #                                                             #
! ###############################################################

module twostream_master_m

Use twostream_inputs_m
Use twostream_writemodules_m
Use twostream_miscsetups_m
Use twostream_geometry_m
Use twostream_solutions_m
Use twostream_bvproblem_m
Use twostream_intensity_m
Use twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & ! Dimensions
          MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          & ! Dimensions !@@ 2p1
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     & ! Inputs
          DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      & ! Input !@@ 2p3 6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs     !@@ 2p1
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            & ! Inputs
          THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, & ! Inputs
          EMISSIVITY, SURFBB, SLTERM_ISOTROPIC, SLTERM_F_0,               & ! Inputs  !@@ 2p3 (Sleave)
          N_FOURIERS, DEG_TO_RAD, PI4, UMOFF, DO_POSTPROCESSING, N_PPSTREAMS, & ! Inputs (geometry)
          PPSTREAM_MASK, CHAPMAN_FACTORS, X0, USER_STREAMS, USER_SECANTS, & ! Inputs (geometry)
          INTENSITY_TOA, INTENSITY_BOA, FLUXES_TOA, FLUXES_BOA,           & ! Outputs !@@ 2p3 (Fluxes)
          RADLEVEL_UP, RADLEVEL_DN,                                       & ! Outputs !@@ 2p2
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,          & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )              ! Exception handling

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Notes 21 december 2012. Observational Geometry Inputs. Marked with !@@ 2p1

!     Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!     Observation-Geometry input control.       DO_USER_OBSGEOMS
!     Observation-Geometry input control.       N_USER_OBSGEOMS
!     User-defined Observation Geometry angles. USER_OBSGEOMS

!  Notes 17 July 2013, Optional output at all levels. Marked with !@@ 2p2
!      New flag for input : DO_2S_LEVELOUT

!  Notes 05 November 2013. Flux output options. Two New Flags
!   DO_MVOUT_ONLY
!   DO_ADDITIONAL_MVOUT

!  Notes 23 January 2014. surface Leaving options. Two New Flags
!   DO_SURFACE_LEAVING
!   DO_SL_ISOTROPIC

!  Notes 25 June 2014. BVProblem control
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

!  Version 2.4 Notes. 15 August 2014. Greens Function implementation.
!  ------------------------------------------------------------------

!   * Greens function particular integral and postprocessing
!   * Use of PPSTREAM maks to reduce coding for observational geometry
!   * Dethreading (removal of MAXTHREADS dimension)
!   * Use of Taylor series routines for limiting cases

!  Subroutine input arguments
!  --------------------------

!  Dimensions :
!      MAXTOTAL       = 2 * MAXLAYERS
!      MAX_GEOMETRIES = MAXBEAMS * MAX_USER_STREAMS * MAX_USER_RELAZMS
!      !@@ MAX_USER_OBSGEOMS >/= MAXBEAMS

!  @@ Rob Spurr, 15 August 2014, Version 2.4, MAXTHREADS dimension removed

      INTEGER, INTENT(IN)        :: MAXMESSAGES, MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_GEOMETRIES
      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  MS-only flag (Mode of operation). NOT REQUIRED
!    IF set, only calculating  MS field
!      LOGICAL, INTENT(IN)        :: DO_MSMODE_2STREAM

!  Plane parallel flag

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL, INTENT(IN)        :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL, INTENT(IN)        :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)        :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)        :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)        :: DO_SURFACE_EMISSION

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)        :: DO_USER_OBSGEOMS !@@ 2p1

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE option flags

      LOGICAL, INTENT(IN)        :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)        :: DO_SL_ISOTROPIC

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Version 2.4, August 2014. Order of Taylor series (N) with Smallness number ( EPS)
!                            (Series including terms up to EPS^N)
      
      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, Version 2.4, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@ 2p1

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS                    !@@ 2p1
      REAL(kind=dp), INTENT(IN)  :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@ 2p1

!  Viewing geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: N_USER_STREAMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_ANGLES  ( MAX_USER_STREAMS )
      INTEGER, INTENT(INOUT)        :: N_USER_RELAZMS
      REAL(kind=dp), INTENT(INOUT)  :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Flux factor

      REAL(kind=dp), INTENT(IN)     :: FLUX_FACTOR

!  Solar geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)        :: NBEAMS
      REAL(kind=dp), INTENT(INOUT)  :: BEAM_SZAS ( MAXBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS )

!  Geometry specification height
!      REAL(kind=dp), INTENT(IN)  :: GEOMETRY_SPECHEIGHT

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS)

!  Atmospheric thermal sources

      REAL(kind=dp), INTENT(IN)  :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Surface properties

      REAL(kind=dp), INTENT(IN)  :: LAMBERTIAN_ALBEDO
      REAL(kind=dp), INTENT(IN)  :: SURFBB

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams    !  NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( MAXBEAMS, 0:1 )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( MAX_USER_STREAMS, 0:1 )

!  Surface thermal sources

      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( MAXBEAMS, 0:1 )

!  Exact Surface-Leaving term
!      REAL(kind=dp) ::  SLTERM_USERANGLES ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted user streams. First order truncated
!      REAL(kind=dp) ::  USER_SLTERM_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )

!  Number of Fourier terms

      INTEGER, INTENT(IN)       :: N_FOURIERS

!  Constants

      REAL(kind=dp), INTENT(IN) :: DEG_TO_RAD, PI4

!  Geometry offset array

      INTEGER, INTENT(IN)       :: UMOFF ( MAX_USER_STREAMS, MAXBEAMS )

!  Post-processing flag (new for Version 2p3)

      LOGICAL, INTENT(IN)       :: DO_POSTPROCESSING

!  Post-processing control mask

      INTEGER, INTENT(IN)       :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(IN) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Cosines and sines

      REAL(kind=dp), INTENT(IN) :: X0  ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN) :: USER_SECANTS ( MAX_USER_STREAMS )

!  Output
!  ------

!  Radiance Results

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

!  Flux output
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(2,MAXBEAMS)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(2,MAXBEAMS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (MAX_GEOMETRIES,0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (MAX_GEOMETRIES,0:MAXLAYERS)

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

!  After application of deltam scaling

      REAL(kind=dp) :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp) :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp) :: ASYMM_TOTAL(MAXLAYERS)

!  Miscsetup operations
!  ====================

!  Pseudo-spherical preparation
!  ----------------------------

!     Last layer to include Particular integral solution
!     Average-secant and initial tramsittance factors for solar beams.
!     Solar beam attenuation

      INTEGER       :: LAYER_PIS_CUTOFF ( MAXBEAMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp) :: TRANS_SOLAR_BEAM( MAXBEAMS )

!  Reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( MAXBEAMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Transmittance factors for user-defined stream angles
!    Computed in the initial setup stage for Fourier m = 0

      REAL(kind=dp) :: T_DELT_USERM ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp) :: ITRANS_USERM ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )

!  Forcing term multiplieres
!  -------------------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp) :: SIGMA_M(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  L'Hopital's rule logical variables

      LOGICAL       :: EMULT_HOPRULE (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Thermal help variables

      REAL(kind=dp) :: THERMCOEFFS ( 2, MAXLAYERS )

!  Fourier-component solutions
!  ===========================

      REAL(kind=dp) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp) :: RADLEVEL_F_UP (0:MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp) :: RADLEVEL_F_DN (0:MAXLAYERS,MAX_USER_STREAMS,MAXBEAMS)

!  Single scatter solutions, commented out. THIS IS MS CODE !!!!
!      REAL(kind=dp) :: INTENSITY_SS_UP(N_GEOMETRIES)
!      REAL(kind=dp) :: INTENSITY_SS_DN(N_GEOMETRIES)

!  Other local variables
!  =====================

!  Local error handling

      CHARACTER(LEN=3) :: CF
      INTEGER          :: FOURIER, STATUS_SUB
      INTEGER          :: N, UA, UM, IB, IBEAM, I, LUM, LUA
      REAL(kind=dp)    :: AZM_ARGUMENT, DFC
      REAL(kind=dp)    :: OMFAC, M1FAC, GDIFF, ALBEDO

!  Local azimuth factors

      REAL(kind=dp)    :: AZMFAC (MAX_USER_RELAZMS,MAX_USER_STREAMS,MAXBEAMS)

!mick - singularity buster output
      LOGICAL          :: SBUST(6)

!  Test variables

      LOGICAL          :: DO_DEBUG_INPUT=.FALSE.
      !LOGICAL          :: DO_DEBUG_INPUT=.TRUE.

real :: e1, e2

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

!  Local user indices; !@@ Only required for OBSGEOM option

      LUM = 1
      LUA = 1

!  TWOSTREAM input debug
!  ---------------------

      IF (DO_DEBUG_INPUT) THEN
        CALL TWOSTREAM_DEBUG_INPUT_MASTER()
      END IF

!  Input checking
!  ==============

call cpu_time(e1)

!  Check input optical values (IOPs in the atmosphere)

     CALL TWOSTREAM_CHECK_INPUTS_OPTICAL &
       ( MAXLAYERS, MAXMESSAGES, NLAYERS,                  & ! input
         DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT,           & ! Input
         STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )    ! Output

call cpu_time(e2)

write(0,*) 'check inp opt', e2-e1

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Get derived inputs
!  ==================

!  Set local atmospheric optical properties (Apply delta 2s scaling)
!  Just copy inputs, if not required

call cpu_time(e1)

      IF ( DO_D2S_SCALING ) THEN
        DO N = 1, NLAYERS
          OMFAC = one - OMEGA_INPUT(N) * D2S_SCALING(N)
          M1FAC = one - D2S_SCALING(N)
          GDIFF = ASYMM_INPUT(N) - D2S_SCALING(N)
          DELTAU_VERT(N) = OMFAC * DELTAU_INPUT(N)
          OMEGA_TOTAL(N) = M1FAC * OMEGA_INPUT(N) / OMFAC
          ASYMM_TOTAL(N) = GDIFF / M1FAC
        ENDDO
      ELSE
        DO N = 1, NLAYERS
          DELTAU_VERT(N) = DELTAU_INPUT(N)
          OMEGA_TOTAL(N) = OMEGA_INPUT(N)
          ASYMM_TOTAL(N) = ASYMM_INPUT(N)
        ENDDO
      ENDIF

!mick fix 1/7/2012 - singularity busters added

!  Note: If running a case close to optical property numerical limits,
!        delta-m scaling may modify omega and/or g in such a way as to make
!        them unphysical or introduce instability; therefore, we recheck
!        omega and g AFTER delta-m scaling and slightly adjust them if necessary

      DO N = 1, NLAYERS
        SBUST = .false.

        !Singularity buster for single scatter albedo
        IF (OMEGA_TOTAL(N) > 0.999999999D0) THEN
          OMEGA_TOTAL(N) = 0.999999999D0
          SBUST(1) = .true.
        ELSE IF (OMEGA_TOTAL(N) < 1.0D-9) THEN
          OMEGA_TOTAL(N) = 1.0E-9_dp
          SBUST(2) = .true.
        END IF

        !Singularity buster for asymmetry parameter
        IF (ASYMM_TOTAL(N) > 0.999999999D0) THEN
          ASYMM_TOTAL(N) = 0.999999999D0
          SBUST(3) = .true.
        ELSE IF (ASYMM_TOTAL(N) < -0.999999999D0) THEN
          ASYMM_TOTAL(N) = -0.999999999D0
          SBUST(4) = .true.
        ELSE IF ((ASYMM_TOTAL(N) >= ZERO) .AND. &
                 (ASYMM_TOTAL(N) < 1.0D-9)) THEN
          ASYMM_TOTAL(N) = 1.0D-9
          SBUST(5) = .true.
        ELSE IF ((ASYMM_TOTAL(N) < ZERO) .AND. &
                 (ASYMM_TOTAL(N) > -1.0D-9)) THEN
          ASYMM_TOTAL(N) = -1.0D-9
          SBUST(6) = .true.
        END IF

        !WRITE(*,*)
        !WRITE(*,'(A,I2)') 'FOR LAYER: ',N
        !DO I=1,6
        !  WRITE(*,'(A,I1,A,L1)') '  SBUST(',I,') = ',SBUST(I)
        !ENDDO

      ENDDO

call cpu_time(e2)

write(0,*) 'derived inp', e2-e1

!  SETUP OPERATIONS (moved from Fourier, Version 2p3, 15 August 2014)
!  ==================================================================

!   MISCSETUPS (4 subroutines)  :
!       average-secant formulation,
!       transmittance setup
!       Thermal setup
!       Beam solution multipliers

!  Prepare quasi-spherical attenuation

call cpu_time(e1)

      if (do_solar_sources) then

      CALL TWOSTREAM_QSPREP &
        ( MAXLAYERS, MAXBEAMS, NLAYERS, NBEAMS, DO_PLANE_PARALLEL,     & ! Input
          DELTAU_VERT, CHAPMAN_FACTORS, X0, DO_DIRECTBEAM,             & ! In/Out
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF, & ! Output
          TRANS_SOLAR_BEAM )                     ! Output

      endif

call cpu_time(e2)

write(0,*) 'qsprep', e2-e1

!  Transmittances and Transmittance factors. !@@ Add flag for Observation Geometry
!    !@@ Add Post-processing flag, 11/5/13
!    !@@ Add solar sources flag 6/29/20 VN

call cpu_time(e1)

      CALL TWOSTREAM_PREPTRANS &
        ( MAXLAYERS, MAX_USER_STREAMS, MAXBEAMS,                      & ! Dimensions
          DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_POSTPROCESSING,      & ! Input flags (2p1,2p3)
          NLAYERS, N_USER_STREAMS, NBEAMS, DELTAU_VERT, USER_SECANTS, & ! Input
          INITIAL_TRANS, AVERAGE_SECANT, LAYER_PIS_CUTOFF,            & ! Input
          T_DELT_MUBAR, T_DELT_USERM, ITRANS_USERM )                    ! Output

call cpu_time(e2)

write(0,*) 'preptrans', e2-e1

!  ** New, October 2011 **. Thermal setup

call cpu_time(e1)

      IF ( DO_THERMAL_EMISSION ) THEN
         CALL TWOSTREAM_THERMALSETUP &
           ( MAXLAYERS, NLAYERS, OMEGA_TOTAL,  & ! Input
             DELTAU_VERT, THERMAL_BB_INPUT,    & ! Input
             THERMCOEFFS )                       ! input/Output
      ENDIF

call cpu_time(e2)

write(0,*) 'therm setup', e2-e1

!   EMULT_MASTER  : Beam source function multipliers.
!      !@@ Add alternative for Observational geometry, 2p1
!      !@@ Avoid altogether if no post-processing

!  Version 2.4 Overhaul-----
!     Rob  Fix 8/15/14  - Small numbers analysis using Taylor parameters
!     Rob  Fix 8/15/14  - Use of PPSTREAM and mask to deal with Obsgeom/Lattice choices
!     Rob  Fix 8/15/14  - Compact code in a single subroutine

call cpu_time(e1)

      IF ( DO_SOLAR_SOURCES.and.DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_EMULTMASTER &
           ( MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,                  & ! Dimensions
             DO_UPWELLING, DO_DNWELLING, NLAYERS, NBEAMS,            & ! Input
             N_PPSTREAMS, PPSTREAM_MASK, TAYLOR_ORDER, TAYLOR_SMALL, & ! Input
             USER_SECANTS, DELTAU_VERT, T_DELT_MUBAR, T_DELT_USERM,  & ! Input
             LAYER_PIS_CUTOFF, ITRANS_USERM, AVERAGE_SECANT,         & ! Input
             SIGMA_M, SIGMA_P, EMULT_HOPRULE, EMULT_UP, EMULT_DN )     ! Output
      ENDIF

call cpu_time(e2)

write(0,*) 'emult', e2-e1

!      write(*,*)'EMULT',EMULT_UP(1,15,1), EMULT_DN(1,16,2)
!      write(*,*)'EMULT',EMULT_UP(2,23,2), EMULT_DN(2,1,1)

!  Albedo

      ALBEDO = LAMBERTIAN_ALBEDO

!  Fourier loop
!  ============

      DO FOURIER = 0, N_FOURIERS

!  Azimuth cosine factors (Fourier = 1). !@@ 2p1, Notice OBSGEOM option
!  !@@ 2p3, not required for FLux-only output

         AZMFAC = zero
         IF ( DO_POSTPROCESSING ) THEN
            IF ( FOURIER .GT. 0 ) THEN
               DFC = DBLE(FOURIER)
               IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
                  DO IB = 1, NBEAMS
                     AZM_ARGUMENT = USER_RELAZMS(IB) * DFC
                     AZMFAC(LUA,LUM,IB) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                  ENDDO
               ELSE
                  DO IB = 1, NBEAMS
                     DO UM = 1, N_USER_STREAMS
                        DO UA = 1, N_USER_RELAZMS
                           AZM_ARGUMENT = USER_RELAZMS(UA) * DFC
                           AZMFAC(UA,UM,IB) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

!  Main call to Lidort Fourier module.
!  ----------------------------------

!  !@@ Add Observational Geometry dimension and control variables   !@@ 2p1
!  !@@ Call statement expanded to include ALL-LEVEL outputs         !@@ 2p2
!  !@@ Call statement expanded to include Flux outputs and control  !@@ 2p3  11/5/13
!  !@@ Call statement expanded to include Sleave inputs and control !@@ 2p3  01/23/14
!  !@@ Call statement expanded to include BVProblem control         !@@ 2p3  06/25/14
!  !@@ Call statement expanded to include Taylor-series control     !@@ 2p4  08/15/14
!  !@@ Call statement changed  to exclude Threading                 !@@ 2p4  08/15/14
!  !@@ Call statement changed  to include TCutoff                   !@@ 2p4  05/14/15

call cpu_time(e1)

         CALL TWOSTREAM_FOURIER_MASTER &
           ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, DO_2S_LEVELOUT, & ! Dimensions + flag
             DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,   & ! Input flags control
             DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Input flags sources
             DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,           & ! Input flags !@@ New line, 2p3
             N_PPSTREAMS, PPSTREAM_MASK,                                      & ! PostProc Inputs
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input flags !@@ 2p3 6/25/14
             BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,   & ! Input !@@ 2p3 6/25/14, 8/15/14
             NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, PI4,           & ! Input integer control
             FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, USER_SECANTS,       & ! Input real control
             ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                               & ! Input real surface
             SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                & ! Input real sleave and thermal
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,              & ! Input real optical
             LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! In/Out Miscsetups
             T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,     & ! In/Out Miscsetups
             SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                            & ! In/Out Miscsetups
             INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,    & ! Output
             FLUXES_TOA, FLUXES_BOA, STATUS_SUB, E_MESSAGE, E_TRACE_1 )         ! Output (modified 2p3, Fluxes)

call cpu_time(e2)

write(0,*) 'fourier master', e2-e1

!  Exception handling

         IF ( STATUS_SUB .NE. 0 ) THEN
            STATUS_EXECUTION = 1
            WRITE(CF,'(I2)')FOURIER
            E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # ' //CF
            RETURN
         ENDIF

!  Fourier summation and Convergence examination
!  SS code not included in this version---------------
!  !@@ Alternative Call for Observationsl Geometry case      !@@ 2p1
!  !@@ Call statements expanded to include ALL-LEVEL outputs !@@ 2p2
!  !@@ Convergence skipped for MVOUT_ONLY option             !@@ 2p3 !mick fix 12/17/2013 - fixed logic

call cpu_time(e1)

         IF ( .not. DO_MVOUT_ONLY ) then
            DO IBEAM = 1, NBEAMS
               IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
                  CALL TWOSTREAM_CONVERGE_OBSGEO &
                    ( MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXLAYERS,        & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,           & ! Inputs ! @@ 2p2
                      NLAYERS, FOURIER, AZMFAC(:,:,IBEAM),                  & ! Inputs ! @@ 2p2
                      INTENSITY_F_UP(:,IBEAM), INTENSITY_F_DN(:,IBEAM),     & ! Inputs
                      RADLEVEL_F_UP(0:,:,IBEAM), RADLEVEL_F_DN(0:,:,IBEAM), & ! Inputs ! @@ 2p2
                      INTENSITY_TOA(IBEAM),   INTENSITY_BOA(IBEAM),         & ! In/Out
                      RADLEVEL_UP(0:,IBEAM),  RADLEVEL_DN(0:,IBEAM)   )       ! In/Out ! @@ 2p2
               ELSE
                  CALL TWOSTREAM_CONVERGE &
                    ( MAX_USER_STREAMS, MAX_USER_RELAZMS,             & ! Dimensions
                      MAX_GEOMETRIES, MAXLAYERS,                      & ! Dimensions ! @@ 2p2
                      DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
                      NLAYERS, FOURIER, N_USER_STREAMS,               & ! Inputs     ! @@ 2p2
                      N_USER_RELAZMS, AZMFAC(:,:,IBEAM), UMOFF(:,IBEAM), & ! Inputs
                      INTENSITY_F_UP(:,IBEAM), INTENSITY_F_DN(:,IBEAM),  & ! Inputs
                      RADLEVEL_F_UP(0:,:,IBEAM), RADLEVEL_F_DN(0:,:,IBEAM), & ! Inputs ! @@ 2p2
                      INTENSITY_TOA,   INTENSITY_BOA,                 & ! In/Out
                      RADLEVEL_UP,     RADLEVEL_DN   )                  ! In/Out     ! @@ 2p2
               ENDIF
            END DO
         ENDIF

call cpu_time(e2)

write(0,*) 'conv', e2-e1

!  End Fourier loop

      ENDDO

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER()

      CALL TWOSTREAM_WRITE_STD_INPUT ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS, MAX_GEOMETRIES,     & 
        MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,  &
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             &
        DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,     &
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,              &
        DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,      &
        BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,  &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  &
        N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      &
        FLUX_FACTOR, NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,      &
        DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,            &
        THERMAL_BB_INPUT, LAMBERTIAN_ALBEDO, SURFBB )

      IF (DO_BRDF_SURFACE) THEN
        CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT (   &
          MAXBEAMS, MAX_USER_STREAMS, &
          NBEAMS, N_USER_STREAMS,     &
          BRDF_F_0, BRDF_F, UBRDF_F, EMISSIVITY )
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL TWOSTREAM_WRITE_SUP_SLEAVE_INPUT ( &
          MAXBEAMS,NBEAMS,&
          SLTERM_ISOTROPIC,SLTERM_F_0 )
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER

END SUBROUTINE TWOSTREAM_MASTER

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER &
        ( MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, DO_2S_LEVELOUT, & ! Dimensions + flag
          DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,   & ! Input flags control
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION, DO_SURFACE_EMISSION,      & ! Input flags sources
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING,           & ! Input flags !@@ New line, 2p3
          N_PPSTREAMS, PPSTREAM_MASK,                                      & ! PostProc Inputs
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input flags !@@ 2p3 6/25/14
          BVPINDEX, BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL, TCUTOFF,   & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, NBEAMS, N_USER_STREAMS, FOURIER, PI4,           & ! Input integer control
          FLUX_FACTOR, STREAM_VALUE, X0, USER_STREAMS, USER_SECANTS,       & ! Input real control
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                               & ! Input real surface
          SLTERM_ISOTROPIC, SLTERM_F_0, SURFBB, EMISSIVITY,                & ! Input real sleave and thermal
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL, THERMCOEFFS,              & ! Input real optical
          LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! In/Out Miscsetups
          T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,     & ! In/Out Miscsetups
          SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                            & ! In/Out Miscsetups
          INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,    & ! Output
          FLUXES_TOA, FLUXES_BOA, STATUS, MESSAGE, TRACE )                   ! Output (modified 2p3, Fluxes)

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions :
!      MAXTOTAL  = 2 * MAXLAYERS
     
      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)        :: MAXBEAMS, MAX_USER_STREAMS

!  Flags
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)  :: DO_THERMAL_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SURFACE_EMISSION
      LOGICAL, INTENT(IN)  :: DO_SOLAR_SOURCES

!   !@@ Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)  :: DO_USER_OBSGEOMS !@@ 2p1

!  !@@ Version 2p3, 11/5/13. Flux output flags, processing flag

      LOGICAL, INTENT(IN)  :: DO_MVOUT_ONLY
      LOGICAL, INTENT(IN)  :: DO_ADDITIONAL_MVOUT
      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!  Post-processing control mask

      INTEGER, INTENT(IN)       :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  !@@ Version 2p3, 1/23/14. Surface leaving control

      LOGICAL, INTENT(IN)  :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  BVP control --- New 6/25/14, Version 2.3 and higher
!   * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!   * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!   * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      INTEGER      , INTENT(IN)  :: BVPINDEX
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)        :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: NBEAMS, N_USER_STREAMS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Surface variables
!  ------------------

!  Lambertian Albedo

     REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!    incident solar directions,  reflected quadrature stream
!    incident quadrature stream, reflected quadrature stream
!    incident solar directions,  reflected user streams -- NOT REQUIRED
!    incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( MAXBEAMS, 0:1 )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( 0:1, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( MAX_USER_STREAMS, 0:1 )

!  ** New **. October 2011. Thermal variables
!  ------------------------------------------

      REAL(kind=dp), INTENT(IN)  :: SURFBB
      REAL(kind=dp), INTENT(IN)  :: EMISSIVITY

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAXBEAMS )

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( MAXBEAMS, 0:1 )

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (MAX_USER_STREAMS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (MAX_USER_STREAMS,MAXBEAMS)

!  Flux output (already initialized here)
!     ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA(2,MAXBEAMS)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA(2,MAXBEAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (MAX_USER_STREAMS,0:MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (MAX_USER_STREAMS,0:MAXLAYERS,MAXBEAMS)

!  Single scatter solutions, commented out in this streamlined version
!      REAL(kind=dp) INTENSITY_SS_UP (N_GEOMETRIES)
!      REAL(kind=dp) INTENSITY_SS_DN (N_GEOMETRIES)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Miscsetups Arrays required as In/OUT
!  ====================================

!  Thermal help variables
!  ----------------------

      REAL(kind=dp), INTENT(INOUT) :: THERMCOEFFS ( 2, MAXLAYERS )

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!     Last layer to include Particular integral solution
!     Average-secant and initial/layer tramsittance factors for solar beams.

      INTEGER      , INTENT(INOUT) :: LAYER_PIS_CUTOFF ( MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: INITIAL_TRANS    ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: AVERAGE_SECANT   ( MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )

!  Solar beam attenuation, reflectance flag
!  ----------------------------------------

      REAL(kind=dp), INTENT(INOUT) :: TRANS_SOLAR_BEAM ( MAXBEAMS )
      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( MAXBEAMS )

!  Transmittance for user-defined stream angles
!  --------------------------------------------

      REAL(kind=dp), INTENT(INOUT) :: ITRANS_USERM ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
      REAL(kind=dp), INTENT(INOUT) :: T_DELT_USERM ( MAX_USER_STREAMS, MAXLAYERS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(INOUT) :: SIGMA_P(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: SIGMA_M(MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(INOUT) :: EMULT_UP (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)
      REAL(kind=dp), INTENT(INOUT) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

!  Local Arrays for Use in Subroutines
!  ===================================

!  Geometry arrays
!  ---------------

!  These just save some Polynomial expansions

      REAL(kind=dp) :: ULP  ( MAX_USER_STREAMS )
      REAL(kind=dp) :: POX  ( MAXBEAMS )
      REAL(kind=dp) :: PX0X ( MAXBEAMS )
      REAL(kind=dp) :: PX11, PXSQ

!  Solar beam Attenuation
!  ----------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solutions. No USER-term required, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( MAXBEAMS )

!  Multiplier arrays (Homogeneous solutions)
!  -----------------

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Solutions to the homogeneous RT equations
!  -----------------------------------------

!  local matrices for eigenvalue computation

      REAL(kind=dp) :: SAB(MAXLAYERS), DAB(MAXLAYERS)

!  Eigensolutions

      REAL(kind=dp) :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp) :: EIGENTRANS(MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp) :: XPOS(2,MAXLAYERS)

!  Green;s function normalization factors
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp) :: NORM_SAVED(MAXLAYERS)

!  Saved help variables

      REAL(kind=dp) :: U_HELP_P(0:1)
      REAL(kind=dp) :: U_HELP_M(0:1)

!  Eigenvectors defined at user-defined stream angles
!     EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp) :: U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_HOMP, H_HOMM

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

!  particular integrals and BVP solution
!  -------------------------------------

!  Solutions at layer boundaries

      REAL(kind=dp) :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp) :: WLOWER(2,MAXLAYERS)

!  Downwelling BOA solution, before reflectance

      REAL(kind=dp) :: H_PARTIC

!  Single-scatter Particular beam solutions at user-defined angles
!    ****** NOT REQUIRED for MS-mode only
!      REAL(kind=dp) :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)
!      REAL(kind=dp) :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)

!  Solution constants of integration, and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC(2,MAXLAYERS)
      REAL(kind=dp) :: MCON_XVEC(2,MAXLAYERS)

!  Beam Solutions (Greens function)
!  --------------------------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: BTERM_SAVE(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp) :: GAMMA_P(MAXLAYERS)

!  Thermal solutions
!  -----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp) :: T_C_MINUS (0:2,MAXLAYERS)
      REAL(kind=dp) :: T_C_PLUS  (0:2,MAXLAYERS)

!  Thermal solution at layer boundaries

      REAL(kind=dp) :: T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) :: T_WLOWER ( 2, MAXLAYERS )

!  Complete layer term solutions

      REAL(kind=dp) :: LAYER_TSUP_UP(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp) :: LAYER_TSUP_DN(MAX_USER_STREAMS,MAXLAYERS)

!  Local help variables
!  --------------------

      INTEGER :: N, IBEAM, i

!  local inclusion flags. ** New October 2011 **, thermal flags
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      LOGICAL :: DO_INCLUDE_SURFACE
      LOGICAL :: DO_INCLUDE_SURFEMISS
      LOGICAL :: DO_INCLUDE_THERMEMISS
      LOGICAL :: DO_INCLUDE_DIRECTBEAM
      LOGICAL :: DO_INCLUDE_MVOUT

!  Flux multiplier and Fourier component numbers

      REAL(kind=dp) :: FLUXMULT
      REAL(kind=dp) :: DELTA_FACTOR
      REAL(kind=dp) :: SURFACE_FACTOR

!  Other local variables

      INTEGER       :: LAYER_PIS_CUTOFFB
      REAL(kind=dp) :: PX0XB
      REAL(kind=dp) :: AVERAGE_SECANTB
      REAL(kind=dp) :: INITIAL_TRANSB
      REAL(kind=dp) :: T_DELT_MUBARB
      REAL(kind=dp) :: OMEGA_TOTALN
      REAL(kind=dp) :: ASYMM_TOTALN
      REAL(kind=dp) :: XPOSN(2)

!  Error tracing

      INTEGER       :: STATUS_SUB

real :: e1, e2

!  ##############
!  initialization
!  ##############

!  Exception handling initialize

      STATUS = 0
      MESSAGE = ' '
      TRACE   = ' '

!  Set local flags
!  ---------------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE.
      IF ( DO_SURFACE_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_SURFEMISS = .TRUE.
        ENDIF
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_THERMEMISS = .TRUE.
        ENDIF
      ENDIF

!  Surface flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .FALSE.
      IF ( DO_BRDF_SURFACE ) THEN
        DO_INCLUDE_SURFACE = .TRUE.
      ELSE
!mick fix 1/30/2015 - refined control logic
        IF ( FOURIER .EQ. 0 .AND. ALBEDO .NE. ZERO) DO_INCLUDE_SURFACE = .TRUE.
        !IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above surface flag has been set)

!mick fix 1/30/2015 - refined control logic
      IF ( DO_SOLAR_SOURCES .and. DO_INCLUDE_SURFACE ) THEN
      !IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, NBEAMS
          DO_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

!  Inclusion of mean value calculation
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        IF ( FOURIER .EQ. 0 ) THEN
          DO_INCLUDE_MVOUT = .TRUE.
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = 2.0_dp
        DELTA_FACTOR   = one
      ELSE
        SURFACE_FACTOR = one
        DELTA_FACTOR   = 2.0_dp
      ENDIF

!  Flux multiplier
!   = 1 / 4.pi with beam sources, 1.0 for thermal

      FLUXMULT   = DELTA_FACTOR

call cpu_time(e1)

!  Reflected Direct beam attenuation.
!  ! @@2p3, 1/23/14 add SLEAVE inputs

      if (do_solar_sources) then

      CALL TWOSTREAM_DIRECTBEAM & 
        ( MAXBEAMS,                             & ! Dimensions
          DO_INCLUDE_SURFACE, DO_BRDF_SURFACE,  & ! Input
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC,  & ! input @@ 2p3
          NBEAMS, FOURIER, FLUX_FACTOR, X0,     & ! Input
          DELTA_FACTOR, ALBEDO, BRDF_F_0(:,FOURIER), & ! Input
          SLTERM_ISOTROPIC, SLTERM_F_0(:,FOURIER),   & ! input @@ 2p3
          TRANS_SOLAR_BEAM, DO_DIRECTBEAM,      & ! Input
          ATMOS_ATTN, DIRECT_BEAM )               ! Output

       endif

call cpu_time(e2)

write(0,*) 'direct beam', e2-e1

call cpu_time(e1)

!  Auxiliary Geometry
!  ! @@2p3, 11/5/13 add Post-processing flag

      CALL TWOSTREAM_AUXGEOM &
        ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING, & ! Dimensions, Flag (2p3
          N_USER_STREAMS, NBEAMS, FOURIER, & ! Input
          X0, USER_STREAMS, STREAM_VALUE,  & ! Input
          PX11, PXSQ, POX, PX0X, ULP )       ! Output

call cpu_time(e2)

write(0,*) 'auxgeom', e2-e1

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

call cpu_time(e1)

!  Start layer loop

      DO N = 1, NLAYERS

         OMEGA_TOTALN = OMEGA_TOTAL(N)
         ASYMM_TOTALN = ASYMM_TOTAL(N)

!  Get Discrete ordinate solutions for this layer
!    Version 2p4. Green's function output = NORM_SAVED

         CALL TWOSTREAM_HOM_SOLUTION &
           ( FOURIER, STREAM_VALUE, PXSQ,             & ! Input
             OMEGA_TOTALN, ASYMM_TOTALN, DELTAU_VERT(N),   & ! Input
             SAB(N), DAB(N), EIGENVALUE(N), EIGENTRANS(N), & ! In/Out
             XPOSN(1:2), NORM_SAVED(N) )                     ! In/Out

!  Get Post-processing ("user") solutions for this layer
!   !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_HOM_USERSOLUTION &
              ( MAX_USER_STREAMS,                                          & ! Dimensions
                N_USER_STREAMS, FOURIER, STREAM_VALUE, PX11,               & ! Input
                USER_STREAMS, ULP, XPOSN(1:2), OMEGA_TOTALN, ASYMM_TOTALN, & ! Input
                U_XPOS(:,N), U_XNEG(:,N), U_HELP_P, U_HELP_M )               ! Output
         ENDIF

         XPOS(1:2,N) = XPOSN(1:2)

!  end layer loop

      ENDDO

call cpu_time(e2)

write(0,*) 'hom solution and user sol', e2-e1

call cpu_time(e1)

!  Prepare homogeneous solution multipliers
!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. User secants, Taylor-series control

      IF ( DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS,             & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS,   & ! Input
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,    & ! Input
             HMULT_1, HMULT_2 )         ! Output
      ENDIF

call cpu_time(e2)

write(0,*) 'hmult', e2-e1

call cpu_time(e1)

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
          ( MAXLAYERS, MAXTOTAL, BVPSCALEFACTOR, DO_PENTADIAG_INVERSE, & ! Dimensions
            DO_INCLUDE_SURFACE, NLAYERS, NTOTAL,                       & ! Input
            DO_BRDF_SURFACE, SURFACE_FACTOR, ALBEDO, BRDF_F(FOURIER),  & ! Input
            XPOS, EIGENTRANS, STREAM_VALUE,                            & ! Input
            H_HOMP, H_HOMM, MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM, & ! Output
            STATUS_SUB, MESSAGE )                                        ! Output

call cpu_time(e2)

write(0,*) 'bvp setup', e2-e1

!  Exception handling for Pentadiagonal Matrix setup

      IF ( STATUS_SUB .NE. 0 ) THEN
         TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
         STATUS = 1 ; RETURN
      ENDIF

call cpu_time(e1)

!  Thermal solutions
!     1. Find the Particular solution (NOT FOR transmittance only)
!     2. Compute thermal layer source terms. (Upwelling and Downwelling)
!       These will be scaled up by factor 4.pi if solar beams as well

!   !@@ 2p3. 11/5/13. Post-processing control
!   !@@ 2p4. 8/15/14. Greens function solution for Regular code
!   !@@ 2p4. 5/14/15. Thermal Cutoff variable introduced.

      IF ( DO_INCLUDE_THERMEMISS ) THEN
         CALL TWOSTREAM_THERMALGFSOLUTION &
           ( MAXLAYERS, NLAYERS, OMEGA_TOTAL, DELTAU_VERT, THERMCOEFFS,  & ! Inputs
             TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED,          & !input
             T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )         ! Outputs
         IF ( DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_THERMALSTERMS &
              ( MAXLAYERS, MAX_USER_STREAMS,                    & ! Dimensions
                DO_UPWELLING, DO_DNWELLING, DO_SOLAR_SOURCES,   & ! Input
                NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS,     & ! Input
                TCUTOFF, T_DELT_USERM, DELTAU_VERT,             & ! Input
                U_XPOS, U_XNEG, HMULT_1, HMULT_2,               & ! Inputs
                T_C_PLUS, T_C_MINUS, TTERM_SAVE,                & ! Inputs
                LAYER_TSUP_UP, LAYER_TSUP_DN  )                   ! Output
         ENDIF
      ENDIF

call cpu_time(e2)

write(0,*) 'thermal gf and sterms', e2-e1

!  Skip the thermal-only section if there are solar sources

      IF ( DO_SOLAR_SOURCES ) GO TO 455

!  ####################################################
!  Complete Radiation Field with Thermal-only solutions
!  ####################################################

!  Only one solution, local direct_beam flag NOT set

      IBEAM = 1
      DO_INCLUDE_DIRECTBEAM = .FALSE.

!  set the BVP PI solution at the lower/upper boundaries

      DO N = 1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = T_WUPPER(I,N)
            WLOWER(I,N) = T_WLOWER(I,N)
         ENDDO
      ENDDO

call cpu_time(e1)

!  Solve boundary value problem (Pentadiagonal solution)
!     Version 2p3. Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
        ( MAXLAYERS, MAXTOTAL,                                 & ! Dimensions
          BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,                & ! BVP control, 6/24/14
          DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAM,           & ! Input
          DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
          NLAYERS, NTOTAL,                                     & ! Input
          SURFACE_FACTOR, ALBEDO, BRDF_F(FOURIER), EMISSIVITY, SURFBB, & ! Input
          DIRECT_BEAM(IBEAM), XPOS, WUPPER, WLOWER,                & ! Input
          STREAM_VALUE, MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM, & ! Input
          H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

call cpu_time(e2)

write(0,*) 'bvp sol', e2-e1

call cpu_time(e1)

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_UPUSER_INTENSITY &
           ( MAXLAYERS, MAX_USER_STREAMS,                                   & ! Dimensions
             DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,         & ! Input !@@ 2p1
             DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,       & ! Input !@@ 2p2
             NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,                         & ! inputs !@@ 2p3 Greens
             N_PPSTREAMS, PPSTREAM_MASK(:,IBEAM),                           & ! inputs
             LAYER_PIS_CUTOFF(IBEAM), PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F(:,FOURIER), & ! inputs
             FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,             & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_P(:,:,IBEAM), ATERM_SAVE, BTERM_SAVE,  & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS(:,IBEAM), ITRANS_USERM(:,:,IBEAM), T_DELT_USERM, & ! Inputs !@@ 2p3 Greens
             T_DELT_MUBAR(:,IBEAM), EIGENTRANS, LCON, LCON_XVEC,            & ! inputs
             MCON, MCON_XVEC, WLOWER, U_XPOS, U_XNEG, HMULT_1, HMULT_2,     & ! inputs
             EMULT_UP(:,:,IBEAM), LAYER_TSUP_UP,                            & ! inputs
             INTENSITY_F_UP(:,IBEAM), RADLEVEL_F_UP(:,0:,IBEAM) )             ! Output !@@ 2p2
      ENDIF

call cpu_time(e2)

write(0,*) 'upuser int', e2-e1

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

      IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_DNUSER_INTENSITY &
           ( MAXLAYERS, MAX_USER_STREAMS,                                   & ! Dimensions
             DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                       & ! Dimensions
             DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                              & ! Inputs !@@ 2p1, 2p2
             NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,                         & ! inputs !@@ 2p3 Greens
             N_PPSTREAMS, PPSTREAM_MASK(:,IBEAM),                           & ! inputs
             LAYER_PIS_CUTOFF(IBEAM), PI4, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, & ! inputs
             GAMMA_P, GAMMA_M, SIGMA_M(:,:,IBEAM), ATERM_SAVE, BTERM_SAVE,  & ! Inputs !@@ 2p3 Greens
             INITIAL_TRANS(:,IBEAM), ITRANS_USERM(:,:,IBEAM), T_DELT_USERM, & ! Inputs !@@ 2p3 Greens
             T_DELT_MUBAR(:,IBEAM), LCON, MCON, U_XPOS, U_XNEG,             & ! Inputs
             HMULT_1, HMULT_2, EMULT_DN(:,:,IBEAM), LAYER_TSUP_DN,          & ! Inputs
             INTENSITY_F_DN(:,IBEAM), RADLEVEL_F_DN(:,0:,IBEAM) )             ! Output !@@ 2p2
      ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

      IF ( DO_INCLUDE_MVOUT ) THEN
         CALL TWOSTREAM_FLUXES &
           ( MAXLAYERS, DO_UPWELLING, DO_DNWELLING,                     & ! Input Dimensions, flags
             DO_INCLUDE_DIRECTBEAM, NLAYERS, PI4, STREAM_VALUE,         & ! Input Control
             FLUX_FACTOR, FLUXMULT, X0(IBEAM), TRANS_SOLAR_BEAM(IBEAM), & ! Input Control
             LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,          & ! Input 2-stream solution
             FLUXES_TOA(:,IBEAM), FLUXES_BOA(:,IBEAM) )                   ! Output
      ENDIF

!  Finish Thermal only.

      RETURN

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Continuation point

 455  CONTINUE

!  Start loop over various solar beams

      DO IBEAM = 1, NBEAMS

!  Solar beam Particular solution
!  ------------------------------

         LAYER_PIS_CUTOFFB = LAYER_PIS_CUTOFF(IBEAM)
         PX0XB = PX0X(IBEAM)

         AVERAGE_SECANTB(:) = AVERAGE_SECANT(:,IBEAM)
         INITIAL_TRANSB(:)  = INITIAL_TRANS(:,IBEAM)
         T_DELT_MUBARB(:)   = T_DELT_MUBAR(:,IBEAM)

!  Version 2p4 Greens function solution

         CALL TWOSTREAM_GBEAM_SOLUTION &
            ( TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT,            & ! Inputs 
              FOURIER, PI4, FLUX_FACTOR,                          & ! Inputs
              LAYER_PIS_CUTOFFB, PX0XB, OMEGA_TOTAL, ASYMM_TOTAL, & ! Inputs
              AVERAGE_SECANTB, INITIAL_TRANSB, T_DELT_MUBARB,     & ! Inputs
              XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Output
              WUPPER, WLOWER )                                      ! Output

!  Add thermal solutions if flagged. NO modulus on the thermal contribution.

         IF ( DO_INCLUDE_THERMEMISS ) THEN
            DO N = 1, NLAYERS
               DO I = 1, 2
                  WUPPER(I,N) = WUPPER(I,N) + T_WUPPER(I,N)
                  WLOWER(I,N) = WLOWER(I,N) + T_WLOWER(I,N)
               ENDDO
            ENDDO
         ENDIF

!  Solve boundary value problem (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

         CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
           ( MAXLAYERS, MAXTOTAL,                                 & ! Dimensions
             BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,                & ! BVP control, 6/24/14
             DO_INCLUDE_SURFACE, DO_DIRECTBEAM(IBEAM),            & ! Input
             DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,               & ! Input
             NLAYERS, NTOTAL,                                     & ! Input
             SURFACE_FACTOR, ALBEDO, BRDF_F(FOURIER), EMISSIVITY, SURFBB, & ! Input
             DIRECT_BEAM(IBEAM), XPOS, WUPPER, WLOWER,                & ! Input
             STREAM_VALUE, MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM, & ! Input
             H_PARTIC, LCON, MCON, LCON_XVEC, MCON_XVEC )           ! Output

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion.
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_INTENSITY &
              ( MAXLAYERS, MAX_USER_STREAMS,                                   & ! Dimensions
                DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_USER_OBSGEOMS,         & ! Input !@@ 2p1
                DO_SOLAR_SOURCES, DO_INCLUDE_THERMEMISS, DO_2S_LEVELOUT,       & ! Input !@@ 2p2
                NLAYERS, N_USER_STREAMS, TAYLOR_ORDER,                         & ! inputs !@@ 2p3 Greens
                N_PPSTREAMS, PPSTREAM_MASK(:,IBEAM),                           & ! inputs
                LAYER_PIS_CUTOFF(IBEAM), PI4, SURFACE_FACTOR, ALBEDO, UBRDF_F(:,FOURIER), & ! inputs
                FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,             & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_P, ATERM_SAVE, BTERM_SAVE,             & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS(:,IBEAM), ITRANS_USERM(:,:,IBEAM), T_DELT_USERM, & ! Inputs !@@ 2p3 Greens
                T_DELT_MUBAR(:,IBEAM), EIGENTRANS, LCON, LCON_XVEC,            & ! inputs
                MCON, MCON_XVEC, WLOWER, U_XPOS, U_XNEG, HMULT_1, HMULT_2,     & ! inputs
                EMULT_UP(:,:,IBEAM), LAYER_TSUP_UP,                            & ! inputs
                INTENSITY_F_UP(:,IBEAM), RADLEVEL_F_UP(:,0:,IBEAM) )             ! Output !@@ 2p2
         ENDIF

!  Downwelling, MSMODE only,
!         !@@ 2p1 New OBSGEOM     option 12/21/12
!         !@@ 2p2 New 2S_LEVELOUT option 07/17/13
!         !@@ 2p3. 11/5/13. Post-processing control

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_DNUSER_INTENSITY &
              ( MAXLAYERS, MAX_USER_STREAMS,                                       & ! Dimensions
                DO_INCLUDE_THERMEMISS, DO_SOLAR_SOURCES,                           & ! Dimensions
                DO_USER_OBSGEOMS, DO_2S_LEVELOUT,                                  & ! Inputs !@@ 2p1, 2p2
                NLAYERS, N_USER_STREAMS,  TAYLOR_ORDER,                            & ! inputs !@@ 2p3 Greens
                N_PPSTREAMS, PPSTREAM_MASK(:,IBEAM),                               & ! inputs
                LAYER_PIS_CUTOFF(IBEAM), PI4, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT, & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_M(:,:,IBEAM), ATERM_SAVE, BTERM_SAVE,      & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANS(:,IBEAM), ITRANS_USERM(:,:,IBEAM), T_DELT_USERM,     & ! Inputs !@@ 2p3 Greens
                T_DELT_MUBAR(:,IBEAM), LCON, MCON, U_XPOS, U_XNEG,                 & ! Inputs
                HMULT_1, HMULT_2, EMULT_DN(:,:,IBEAM), LAYER_TSUP_DN,              & ! Inputs
                INTENSITY_F_DN(:,IBEAM), RADLEVEL_F_DN(:,0:,IBEAM) )                 ! Output !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

         IF ( DO_INCLUDE_MVOUT ) THEN
            CALL TWOSTREAM_FLUXES &
              ( MAXLAYERS, DO_UPWELLING, DO_DNWELLING,                     & ! Input Dimensions, flags
                DO_DIRECTBEAM(IBEAM), NLAYERS, PI4, STREAM_VALUE,          & ! Input Control
                FLUX_FACTOR, FLUXMULT, X0(IBEAM), TRANS_SOLAR_BEAM(IBEAM), & ! Input Control
                LCON_XVEC, MCON_XVEC, EIGENTRANS, WUPPER, WLOWER,          & ! Input 2-stream solution
                FLUXES_TOA(:,IBEAM), FLUXES_BOA(:,IBEAM) )                   ! Output
         ENDIF

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER

end module twostream_master_m
