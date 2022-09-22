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
! #            TWOSTREAM_MASTER_SOLAR_OBS (top-level master)    #
! #            TWOSTREAM_FOURIER_MASTER_SOLAR_OBS               #
! #                                                             #
! ###############################################################

module twostream_master_solar_obs_m

Use twostream_miscsetups_m
Use twostream_miscsetups_solar_obs_m
Use twostream_miscsetups_solar_obs_pp_m
Use twostream_converge_obs_m
Use twostream_writemodules_solar_obs_m
Use twostream_solutions_m
Use twostream_solutions_solar_m
Use twostream_bvproblem_m
Use twostream_intensity_solar_obs_m
Use twostream_fluxes_m

!Use twostream_geometry_m
!Use twostream_geometry_obs_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_MASTER_SOLAR_OBS &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAX_USER_OBSGEOMS,             & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,   & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_INCLUDE_MVOUT,                                 & ! Inputs     !@@ 2p3
          DO_POSTPROCESSING, DO_D2S_SCALING, DO_BRDF_SURFACE,              & ! Inputs     !@@ 2p1
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input !@@ 2p3 6/25/14
          BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,                      & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,   & ! Inputs     !@@ 2p1
          FLUX_FACTOR, SURFACE_FACTOR, DELTA_FACTOR,                       & ! Inputs
          DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,             & ! Inputs
          ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0, & ! Inputs  !@@ 2p3 (Sleave)
          N_FOURIERS, PI4, AVERAGE_SECANT_PP, CHAPMAN_FACTORS, X0,         & ! Inputs (geometry)
          USER_STREAMS, USER_SECANTS, AZMFAC, PX11, PXSQ, PX0X, ULP,       & ! Inputs (geometry)
          INTENSITY_TOA, INTENSITY_BOA,                                    & ! Outputs
          RADLEVEL_UP, RADLEVEL_DN,                                        & ! Outputs !@@ 2p2
          FLUXES_TOA_SOLAR, FLUXES_BOA_SOLAR,                              & ! Outputs !@@ 2p3 (Fluxes)
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,           & ! Exception handling
          STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )               ! Exception handling

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
      INTEGER, INTENT(IN)        :: MAX_USER_OBSGEOMS

!  Directional Flags

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING

!  Plane parallel flag

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL, INTENT(IN)        :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_MVOUT(0:1)

!  Post-processing flag (new for Version 2p3)

      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Deltam-2stream scaling flag

      LOGICAL, INTENT(IN)        :: DO_D2S_SCALING

!  BRDF surface flag

      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE option flags

      LOGICAL, INTENT(IN)        :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)        :: DO_SL_ISOTROPIC

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

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@ 2p1

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS                    !@@ 2p1
      REAL(kind=dp), INTENT(IN)  :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@ 2p1

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Surface reflectance factors

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR(0:1)
      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTOR(0:1)

!  Atmospheric optical properties

      REAL(kind=dp), INTENT(IN)  :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: D2S_SCALING (MAXLAYERS)

!  Surface albedo

      REAL(kind=dp), INTENT(IN)  :: ALBEDO

!  BRDF fourier components
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!  incident solar directions,  reflected quadrature stream
!  incident quadrature stream, reflected quadrature stream
!  incident solar directions,  reflected user streams    !  NOT REQUIRED
!  incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0  ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp), INTENT(IN)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0 ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_F   ( MAX_USER_OBSGEOMS, 0:1 )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!  Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAX_USER_OBSGEOMS )

!  Fourier components of Surface-leaving terms:
!  Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0 ( MAX_USER_OBSGEOMS, 0:1 )

!  Exact Surface-Leaving term
!      REAL(kind=dp) ::  SLTERM_USERANGLES ( MAX_USER_OBSGEOMS )
!  Fourier components of Surface-leaving terms:
!  Every solar direction, SL-transmitted user streams. First order truncated
!      REAL(kind=dp) ::  USER_SLTERM_F_0 ( MAX_USER_OBSGEOMS, 0:1 )

!  Number of Fourier terms

      INTEGER, INTENT(IN)       :: N_FOURIERS

!  Constants

      REAL(kind=dp), INTENT(IN) :: PI4

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(IN) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAX_USER_OBSGEOMS )

!  Average secant (plane-parallel case)

      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT_PP ( MAX_USER_OBSGEOMS )

!  Cosines and sines

      REAL(kind=dp), INTENT(IN) :: X0  ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: USER_STREAMS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: USER_SECANTS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: AZMFAC ( MAX_USER_OBSGEOMS )

!  Polynomial expansions

      REAL(kind=dp), INTENT(IN) :: PX11, PXSQ(0:1)
      REAL(kind=dp), INTENT(IN) :: PX0X ( 0:1, MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: ULP ( MAX_USER_OBSGEOMS )

!  Output
!  ------

!  Radiance Results

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_USER_OBSGEOMS)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (0:MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (0:MAXLAYERS,MAX_USER_OBSGEOMS)

!  Flux output
!  ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp), INTENT(INOUT) :: FLUXES_TOA_SOLAR(2,MAX_USER_OBSGEOMS)
     REAL(kind=dp), INTENT(INOUT) :: FLUXES_BOA_SOLAR(2,MAX_USER_OBSGEOMS)

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

      INTEGER       :: LAYER_PIS_CUTOFF ( MAX_USER_OBSGEOMS )
      REAL(kind=dp) :: INITIAL_TRANS    ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp) :: AVERAGE_SECANT   ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp) :: TRANS_SOLAR_BEAM( MAX_USER_OBSGEOMS )

!  Reflectance flags

      LOGICAL       :: DO_DIRECTBEAM ( MAX_USER_OBSGEOMS )

!  Transmittance Setups
!  --------------------

!  Transmittance factors for average secant stream

      REAL(kind=dp) :: T_DELT_MUBAR ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Forcing term multiplieres
!  -------------------------

!  coefficient functions for user-defined angles

      REAL(kind=dp) :: SIGMA_P(MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp) :: SIGMA_M(MAXLAYERS,MAX_USER_OBSGEOMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp) :: EMULT_UP (MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp) :: EMULT_DN (MAXLAYERS,MAX_USER_OBSGEOMS)

!  Fourier-component solutions
!  ===========================

      REAL(kind=dp) :: INTENSITY_F_UP (MAX_USER_OBSGEOMS)
      REAL(kind=dp) :: INTENSITY_F_DN (MAX_USER_OBSGEOMS)

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp) :: RADLEVEL_F_UP (0:MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp) :: RADLEVEL_F_DN (0:MAXLAYERS,MAX_USER_OBSGEOMS)

!  Other local variables
!  =====================

!  Local error handling

      CHARACTER(LEN=3) :: CF
      INTEGER          :: FOURIER, STATUS_SUB
      INTEGER          :: N, IBEAM
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
        CALL TWOSTREAM_DEBUG_INPUT_MASTER_SOLAR_OBS()
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

!  MISCSETUPS (ONE integrated subroutine, which does the following)  :
!       1. average-secant formulation (solar sources only),
!       2. transmittance setup (T_DELT_MUBAR: solar sources only; T_DELT_USERM, ITRANS_USERM: postprocessing only)
!       3. Beam solution multipliers (solar sources and post processing only)

!  1. Solar setup

!call cpu_time(e1)

         DO_DIRECTBEAM = .TRUE.
         IF (DO_PLANE_PARALLEL) THEN
            CALL TWOSTREAM_QSPREP_PREPTRANS_EMULT_OBS_PP &
               ( MAXLAYERS, MAX_USER_OBSGEOMS,                         & ! Dimensions
                 DO_UPWELLING, DO_DNWELLING, DO_POSTPROCESSING,        & ! Flags
                 NLAYERS, N_USER_OBSGEOMS, TAYLOR_SMALL, TAYLOR_ORDER, & ! Input
                 DELTAU_VERT, AVERAGE_SECANT_PP, USER_SECANTS,         & ! Input
                 DO_DIRECTBEAM,                                        & ! In/Out
                 LAYER_PIS_CUTOFF, INITIAL_TRANS, TRANS_SOLAR_BEAM,    & ! Output
                 T_DELT_MUBAR, ITRANS_USERM, T_DELT_USERM,             & ! Output
                 SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN )                  ! Output
         ELSE
            CALL TWOSTREAM_QSPREP_PREPTRANS_EMULT_OBS &
               ( MAXLAYERS, MAX_USER_OBSGEOMS,                               & ! Dimensions
                 DO_UPWELLING, DO_DNWELLING, DO_POSTPROCESSING,              & ! Flags
                 NLAYERS, N_USER_OBSGEOMS, TAYLOR_SMALL, TAYLOR_ORDER,       & ! Input
                 DELTAU_VERT, CHAPMAN_FACTORS, USER_SECANTS,                 & ! Input
                 DO_DIRECTBEAM,                                              & ! In/Out
                 LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,            & ! Output
                 TRANS_SOLAR_BEAM, T_DELT_MUBAR, ITRANS_USERM, T_DELT_USERM, & ! Output
                 SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN )                        ! Output
         ENDIF

!call cpu_time(e2)

!write(0,*) 'qsprep time = ', e2-e1

!      write(*,*)'EMULT',EMULT_UP(15,1), EMULT_DN(16,1)
!      write(*,*)'EMULT',EMULT_UP(23,1), EMULT_DN(2,1)

!  Fourier loop
!  ============

!ftime = 0.0
!ctime = 0.0

      DO FOURIER = 0, N_FOURIERS

!  Main call to TWOSTREAM Fourier module
!  -------------------------------------

!call cpu_time(e1)

         CALL TWOSTREAM_FOURIER_MASTER_SOLAR_OBS &
           ( MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, DO_2S_LEVELOUT,           & ! Dimensions + flag
             DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_PLANE_PARALLEL,   & ! Input flags control
             DO_POSTPROCESSING, DO_INCLUDE_MVOUT(FOURIER),                     & ! Input flag
             DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,        & ! Input flags !@@ 2p3 6/25/14
             BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,                       & ! Input !@@ 2p3 6/25/14, 8/15/14
             NLAYERS, NTOTAL, N_USER_OBSGEOMS, FOURIER, PI4,                   & ! Input integer control
             FLUX_FACTOR, SURFACE_FACTOR(FOURIER), DELTA_FACTOR(FOURIER),      & ! Input real control
             STREAM_VALUE, X0, AVERAGE_SECANT_PP,                              & ! Input real control
             USER_STREAMS, USER_SECANTS,                                       & ! Input real control
             PX11, PXSQ(FOURIER), PX0X(FOURIER,:), ULP,                        & ! Inputs
             ALBEDO, BRDF_F_0(:,FOURIER), BRDF_F(FOURIER), UBRDF_F(:,FOURIER), & ! Input real surface
             SLTERM_ISOTROPIC, SLTERM_F_0(:,FOURIER),                          & ! Input real sleave and thermal
             DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,                            & ! Input real optical
             LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,    & ! In/Out Miscsetups
             T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,      & ! In/Out Miscsetups
             SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                             & ! In/Out Miscsetups
             INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,     & ! Outputs (Radiances)
             FLUXES_TOA_SOLAR, FLUXES_BOA_SOLAR,                               & ! Outputs (Fluxes)
             STATUS_SUB, E_MESSAGE, E_TRACE_1 )                                  ! Outputs (Errors)

!  Exception handling

         IF ( STATUS_SUB .NE. 0 ) THEN
            STATUS_EXECUTION = 1
            WRITE(CF,'(I2)')FOURIER
            E_TRACE_2 = 'Error from 2S_FOURIER_MASTER, Fourier # ' //CF
            RETURN
         ENDIF

!call cpu_time(e2)

!ftime = ftime+e2-e1

!  Fourier summation and Convergence examination

!call cpu_time(e1)

         IF ( .NOT. DO_MVOUT_ONLY ) THEN
            DO IBEAM = 1, N_USER_OBSGEOMS
               CALL TWOSTREAM_CONVERGE_OBSGEO_SOLAR &
                  ( MAXLAYERS,        & ! Dimensions ! @@ 2p2
                    DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,           & ! Inputs ! @@ 2p2
                    NLAYERS, FOURIER, AZMFAC(IBEAM),                      & ! Inputs ! @@ 2p2
                    INTENSITY_F_UP(IBEAM), INTENSITY_F_DN(IBEAM),         & ! Inputs
                    RADLEVEL_F_UP(0:,IBEAM), RADLEVEL_F_DN(0:,IBEAM),     & ! Inputs ! @@ 2p2
                    INTENSITY_TOA(IBEAM),   INTENSITY_BOA(IBEAM),         & ! In/Out
                    RADLEVEL_UP(0:,IBEAM),  RADLEVEL_DN(0:,IBEAM)   )       ! In/Out ! @@ 2p2
            ENDDO
         ENDIF

!call cpu_time(e2)

!ctime = ctime+e2-e1

!  End Fourier loop

      ENDDO

!write(0,*) 'fourier subroutine time = ', ftime
!write(0,*) 'converge routine time = ', ctime

!  Finish

      RETURN

      CONTAINS

      SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER_SOLAR_OBS()

      CALL TWOSTREAM_WRITE_STD_INPUT_SOLAR_OBS ( &
        MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAX_USER_OBSGEOMS,                  &
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT,        &
        DO_INCLUDE_MVOUT, DO_POSTPROCESSING, DO_PENTADIAG_INVERSE,            &
        DO_D2S_SCALING, DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, &
        BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,                           &
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,        &
        FLUX_FACTOR, DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING, &
        ALBEDO )

      IF (DO_BRDF_SURFACE) THEN
        CALL TWOSTREAM_WRITE_SUP_BRDF_INPUT_SOLAR_OBS (   &
          MAX_USER_OBSGEOMS, N_USER_OBSGEOMS, &
          BRDF_F_0, BRDF_F, UBRDF_F )
      END IF

      IF (DO_SURFACE_LEAVING) THEN
        CALL TWOSTREAM_WRITE_SUP_SLEAVE_INPUT_OBS ( &
          MAX_USER_OBSGEOMS, N_USER_OBSGEOMS,&
          SLTERM_ISOTROPIC, SLTERM_F_0 )
      END IF

      END SUBROUTINE TWOSTREAM_DEBUG_INPUT_MASTER_SOLAR_OBS

END SUBROUTINE TWOSTREAM_MASTER_SOLAR_OBS

!

SUBROUTINE TWOSTREAM_FOURIER_MASTER_SOLAR_OBS &
        ( MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, DO_2S_LEVELOUT,          & ! Dimensions + flag
          DO_UPWELLING, DO_DNWELLING, DO_BRDF_SURFACE, DO_PLANE_PARALLEL,  & ! Input flags control
          DO_POSTPROCESSING, DO_INCLUDE_MVOUTM,                            & ! Input flag
          DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PENTADIAG_INVERSE,       & ! Input flags !@@ 2p3 6/25/14
          BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,                      & ! Input !@@ 2p3 6/25/14, 8/15/14
          NLAYERS, NTOTAL, N_USER_OBSGEOMS, FOURIER, PI4,                  & ! Input integer control
          FLUX_FACTOR, SURFACE_FACTORM, DELTA_FACTORM,                     & ! Input real control
          STREAM_VALUE, X0, AVERAGE_SECANT_PP,                             & ! Input real control
          USER_STREAMS, USER_SECANTS,                                      & ! Input real control
          PX11, PXSQM, PX0XM, ULP,                                         & ! Inputs
          ALBEDO, BRDF_F_0M, BRDF_FM, UBRDF_FM,                            & ! Input real surface
          SLTERM_ISOTROPIC, SLTERM_F_0M,                                   & ! Input real sleave and thermal
          DELTAU_VERT, OMEGA_TOTAL, ASYMM_TOTAL,                           & ! Input real optical
          LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,   & ! Inputs Miscsetups
          T_DELT_USERM, ITRANS_USERM, TRANS_SOLAR_BEAM, DO_DIRECTBEAM,     & ! Inputs Miscsetups
          SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN,                            & ! Inputs Miscsetups
          INTENSITY_F_UP, INTENSITY_F_DN, RADLEVEL_F_UP, RADLEVEL_F_DN,    & ! Output
          FLUXES_TOA_SOLAR, FLUXES_BOA_SOLAR,                              & ! Output
          STATUS, MESSAGE, TRACE )                                           ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions :
!      MAXTOTAL  = 2 * MAXLAYERS
     
      INTEGER, INTENT(IN)  :: MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)  :: MAX_USER_OBSGEOMS

!  Flags

      LOGICAL, INTENT(IN)  :: DO_2S_LEVELOUT
      LOGICAL, INTENT(IN)  :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)  :: DO_BRDF_SURFACE
      LOGICAL, INTENT(IN)  :: DO_PLANE_PARALLEL

!  Post processing flag

      LOGICAL, INTENT(IN)  :: DO_POSTPROCESSING

!  Surface leaving control

      LOGICAL, INTENT(IN)  :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)  :: DO_SL_ISOTROPIC

!  Other flag

      LOGICAL, INTENT(IN)  :: DO_INCLUDE_MVOUTM

!  BVP control --- New 6/25/14, Version 2.3 and higher
!  * PentaDiagonal Inverse flag (BVP solved from bottom to top). Only for BVPIndex = 1
!  * BVP Index : 0 = LAPACK, 1 = Penta # 1 (original), 2 = Penta # 2 (new, 2012 Kanal paper)
!  * BVP Scale Factor. Debug only. Set this to 1.0 on input

      LOGICAL      , INTENT(IN)  :: DO_PENTADIAG_INVERSE
      REAL(kind=dp), INTENT(IN)  :: BVPSCALEFACTOR

!  Order of Taylor series (including terms up to EPS^n)
      
      INTEGER, intent(in)        :: TAYLOR_ORDER
      REAL(kind=dp), intent(in)  :: TAYLOR_SMALL

!  Numbers

      INTEGER, INTENT(IN)  :: NLAYERS, NTOTAL
      INTEGER, INTENT(IN)  :: N_USER_OBSGEOMS

!  Input Fourier component number

      INTEGER, INTENT(IN)        :: FOURIER

!  4pi

      REAL(kind=dp), INTENT(IN)  :: PI4

!  Flux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Surface reflectance factors

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTORM
      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTORM

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Geometry

      REAL(kind=dp), INTENT(IN)  :: X0           ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_OBSGEOMS )

!  Average secant (plane-parallel case)

      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT_PP ( MAX_USER_OBSGEOMS )

!  Polynomial expansions

      REAL(kind=dp), INTENT(IN) :: PX11, PXSQM
      REAL(kind=dp), INTENT(IN) :: PX0XM ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: ULP ( MAX_USER_OBSGEOMS )

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

      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0M  ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: BRDF_FM
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0M ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: UBRDF_FM   ( MAX_USER_OBSGEOMS )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  --------------------------------------------

!    Do not require any first-order inputs (exact or Fourier)

!  Isotropic Surface leaving term (if flag set)

      REAL(kind=dp), INTENT(IN) ::  SLTERM_ISOTROPIC ( MAX_USER_OBSGEOMS )

!  Fourier components of Surface-leaving terms:
!  Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) ::  SLTERM_F_0M ( MAX_USER_OBSGEOMS )

!  Optical properties
!  ------------------

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: OMEGA_TOTAL(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM_TOTAL(MAXLAYERS)

!  Output
!  ------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP (MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN (MAX_USER_OBSGEOMS)

!  Flux output

     REAL(kind=dp), INTENT(OUT) :: FLUXES_TOA_SOLAR(2,MAX_USER_OBSGEOMS)
     REAL(kind=dp), INTENT(OUT) :: FLUXES_BOA_SOLAR(2,MAX_USER_OBSGEOMS)

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (0:MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (0:MAXLAYERS,MAX_USER_OBSGEOMS)

!  Exception handling

      INTEGER      , INTENT(OUT)  :: STATUS
      CHARACTER*(*), INTENT(OUT)  :: MESSAGE, TRACE

!  Miscsetups Arrays
!  =================

!  Solar beam pseudo-spherical setup
!  ---------------------------------

!  Last layer to include Particular integral solution
!  Average-secant and initial/layer transmittance factors for solar beams.

      INTEGER      , INTENT(IN) :: LAYER_PIS_CUTOFF ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: INITIAL_TRANS    ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: AVERAGE_SECANT   ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: T_DELT_MUBAR     ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Solar beam attenuation, reflectance flag
!  ----------------------------------------

      REAL(kind=dp), INTENT(IN) :: TRANS_SOLAR_BEAM ( MAX_USER_OBSGEOMS )
      LOGICAL      , INTENT(INOUT) :: DO_DIRECTBEAM ( MAX_USER_OBSGEOMS )

!  Transmittance for user-defined stream angles
!  --------------------------------------------

      REAL(kind=dp), INTENT(IN) :: ITRANS_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Multiplier arrays
!  -----------------

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(IN) :: SIGMA_P(MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(IN) :: SIGMA_M(MAXLAYERS,MAX_USER_OBSGEOMS)

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(IN) :: EMULT_UP (MAXLAYERS,MAX_USER_OBSGEOMS)
      REAL(kind=dp), INTENT(IN) :: EMULT_DN (MAXLAYERS,MAX_USER_OBSGEOMS)

!  Local Arrays for Use in Subroutines
!  ===================================

!  Solar beam Attenuation
!  ----------------------

!  Atmospheric attenuation

      REAL(kind=dp) :: ATMOS_ATTN ( MAX_USER_OBSGEOMS )

!  Direct beam solutions. No USER-term required, MS-mode only

      REAL(kind=dp) :: DIRECT_BEAM ( MAX_USER_OBSGEOMS )

!  Multiplier arrays (Homogeneous solutions)
!  -----------------

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp) :: HMULT_1(MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp) :: HMULT_2(MAX_USER_OBSGEOMS,MAXLAYERS)

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

      REAL(kind=dp) :: U_XPOS(MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp) :: U_XNEG(MAX_USER_OBSGEOMS,MAXLAYERS)

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
      REAL(kind=dp) :: WUPPER21
      REAL(kind=dp) :: WLOWER1NL

!  Solution constants of integration and related quantities

      REAL(kind=dp) :: LCON(MAXLAYERS)
      REAL(kind=dp) :: MCON(MAXLAYERS)
      REAL(kind=dp) :: LCON_XVEC1NL, MCON_XVEC1NL, LCON_XVEC21, MCON_XVEC21

!  Beam Solutions (Greens function)
!  --------------------------------

!  Saved quantities for the Green function solution

      REAL(kind=dp) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp) :: BTERM_SAVE(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp) :: GAMMA_P(MAXLAYERS)

!  Local help variables
!  --------------------

      INTEGER :: IBEAM

!  local inclusion flags. ** New October 2011 **, thermal flags
! !@@ 2p3 11/5/13. Control for the Flux calculation  

      LOGICAL :: DO_INCLUDE_SURFACE

!  Other local variables

      INTEGER       :: LAYER_PIS_CUTOFFB
      REAL(kind=dp) :: EIGENTRANS1
      REAL(kind=dp) :: INITIAL_TRANSB(MAXLAYERS)
      REAL(kind=dp) :: T_DELT_MUBARB(MAXLAYERS)

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
!mick fix 1/30/2015 - refined control logic
        IF ( FOURIER .EQ. 0 .AND. ALBEDO .NE. ZERO) DO_INCLUDE_SURFACE = .TRUE.
        !IF ( FOURIER .EQ. 0 ) DO_INCLUDE_SURFACE = .TRUE.
      ENDIF

!  Direct beam flag (only if above surface flag has been set)

      IF ( DO_INCLUDE_SURFACE ) THEN
        DO IBEAM = 1, N_USER_OBSGEOMS
          DO_DIRECTBEAM(IBEAM) = .TRUE.
        ENDDO
      ELSE
        DO IBEAM = 1, N_USER_OBSGEOMS
          DO_DIRECTBEAM(IBEAM) = .FALSE.
        ENDDO
      ENDIF

!  Reflected Direct beam attenuation.
!  ! @@2p3, 1/23/14 add SLEAVE inputs

         IF (.NOT. DO_INCLUDE_SURFACE) THEN
            ATMOS_ATTN  = ZERO
            DIRECT_BEAM = ZERO
         ELSE
            CALL TWOSTREAM_DIRECTBEAM & 
          ( MAX_USER_OBSGEOMS,                                    & ! Dimension
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, & ! Flags
            N_USER_OBSGEOMS, FOURIER, FLUX_FACTOR, X0,            & ! Inputs
            DELTA_FACTORM, ALBEDO, BRDF_F_0M,     & ! Inputs
            SLTERM_ISOTROPIC, SLTERM_F_0M,        & ! Inputs
            TRANS_SOLAR_BEAM, DO_DIRECTBEAM,      & ! Inputs
            ATMOS_ATTN, DIRECT_BEAM )               ! Outputs
          ENDIF

!  #########################################
!   RTE HOMOGENEOUS SOLUTIONS and BVP SETUP
!  #########################################

!  Get Discrete ordinate solutions for all layers
!  Version 2p4. Green's function output = NORM_SAVED

       CALL TWOSTREAM_HOM_SOLUTION_SOLAR &
          ( MAXLAYERS, NLAYERS, FOURIER, STREAM_VALUE, PXSQM, & ! Input
            OMEGA_TOTAL, ASYMM_TOTAL, DELTAU_VERT,   & ! Input
            EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED ) ! Output

!  Get Post-processing ("user") solutions for this layer
!  !@@ 2p3. 11/5/13. Post-processing control

       IF ( DO_POSTPROCESSING ) THEN
          CALL TWOSTREAM_HOM_USERSOLUTION_SOLAR &
             ( MAXLAYERS, MAX_USER_OBSGEOMS,                          & ! Dimensions
               NLAYERS, N_USER_OBSGEOMS, FOURIER, STREAM_VALUE, PX11, & ! Input
               USER_STREAMS, ULP, XPOS, OMEGA_TOTAL, ASYMM_TOTAL,     & ! Input
               U_XPOS, U_XNEG )                                         ! Output
       ENDIF

!  Prepare homogeneous solution multipliers
!  !@@ 2p3. 11/5/13. Post-processing control
!  !@@ 2p4. 8/15/14. User secants, Taylor-series control

      IF ( DO_POSTPROCESSING ) THEN
         CALL TWOSTREAM_HMULT_MASTER &
            ( MAXLAYERS, MAX_USER_OBSGEOMS,            & ! Dimensions
              TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT, & ! Inputs 
              NLAYERS, N_USER_OBSGEOMS, USER_SECANTS,  & ! Input
              EIGENVALUE, EIGENTRANS, T_DELT_USERM,    & ! Input
              HMULT_1, HMULT_2 )                         ! Output
      ENDIF

!  Boundary value problem - MATRIX PREPARATION (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

      CALL TWOSTREAM_BVP_MATSETUP_PENTADIAG &
          ( MAXLAYERS, MAXTOTAL, BVPSCALEFACTOR, DO_PENTADIAG_INVERSE, & ! Dimensions
            DO_INCLUDE_SURFACE, NLAYERS, NTOTAL,                       & ! Input
            DO_BRDF_SURFACE, SURFACE_FACTORM, ALBEDO, BRDF_FM,         & ! Input
            XPOS, EIGENTRANS, STREAM_VALUE,                            & ! Input
            MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM,                 & ! Output
            STATUS_SUB, MESSAGE )                                        ! Output

!  Exception handling for Pentadiagonal Matrix setup

      IF ( STATUS_SUB .NE. 0 ) THEN
         TRACE  = 'Call BVP_MATSETUP_PENTADIAG in 2S_FOURIER_MASTER'
         STATUS = 1 ; RETURN
      ENDIF

!  ##################################################
!  Complete Radiation Field with Solar Beam solutions
!  ##################################################

!  Start loop over various solar beams

      DO IBEAM = 1, N_USER_OBSGEOMS

!  Solar beam Particular solution
!  ------------------------------

         LAYER_PIS_CUTOFFB = LAYER_PIS_CUTOFF(IBEAM)
         INITIAL_TRANSB(:)  = INITIAL_TRANS(:,IBEAM)
         T_DELT_MUBARB(:)   = T_DELT_MUBAR(:,IBEAM)

!  Version 2p4 Greens function solution

         CALL TWOSTREAM_GBEAM_SOLUTION &
            ( MAXLAYERS, NLAYERS,                                       & ! Inputs
              DO_PLANE_PARALLEL, DO_POSTPROCESSING,                     & ! Inputs
              TAYLOR_ORDER, TAYLOR_SMALL, DELTAU_VERT,                  & ! Inputs 
              FOURIER, PI4, FLUX_FACTOR,                                & ! Inputs
              LAYER_PIS_CUTOFFB, PX0XM(IBEAM), OMEGA_TOTAL, ASYMM_TOTAL, & ! Inputs
              AVERAGE_SECANT_PP(IBEAM), AVERAGE_SECANT(:,IBEAM),        & ! Inputs
              INITIAL_TRANSB, T_DELT_MUBARB,                            & ! Inputs
              XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,                 & ! Inputs
              GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,                 & ! Output
              WUPPER, WLOWER )                                            ! Output

!  Solve boundary value problem (Pentadiagonal solution)
!     Pentadiagonal inverse option introduced, 25 June 2014

         CALL TWOSTREAM_BVP_SOLUTION_PENTADIAG &
           ( MAXLAYERS, MAXTOTAL,                                 & ! Dimensions
             BVPSCALEFACTOR, DO_PENTADIAG_INVERSE,                & ! BVP control, 6/24/14
             DO_INCLUDE_SURFACE, DO_INCLUDE_SURFACE,              & ! Input
             .FALSE., DO_BRDF_SURFACE,                            & ! Input
             NLAYERS, NTOTAL,                                     & ! Input
             SURFACE_FACTORM, ALBEDO, BRDF_FM, ONE, ZERO,         & ! Input
             DIRECT_BEAM(IBEAM), WUPPER, WLOWER,                  & ! Input
             STREAM_VALUE, MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM, & ! Input
             LCON, MCON )                                           ! Output

         IF (( DO_UPWELLING .and. DO_POSTPROCESSING ) .OR. &
             ( DO_DNWELLING .and. DO_INCLUDE_MVOUTM )) THEN
            EIGENTRANSNL = EIGENTRANS(NLAYERS) 
            LCON_XVEC1NL = LCON(NLAYERS)*XPOS(1,NLAYERS)
            MCON_XVEC1NL = MCON(NLAYERS)*XPOS(2,NLAYERS)
            WLOWER1NL    = WLOWER(1,NLAYERS)
         ENDIF

! ##################################
!   Radiance Field Post Processing
! ##################################

!  upwelling, MSMODE only, no Direct Beam inclusion.

         IF ( DO_UPWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_UPUSER_INTENSITY_SOLAR_OBS &
              ( MAXLAYERS,                                                       & ! Dimension
                DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_2S_LEVELOUT,             & ! Input !@@ 2p2
                NLAYERS, TAYLOR_ORDER,                                           & ! inputs !@@ 2p3 Greens
                LAYER_PIS_CUTOFFB, SURFACE_FACTORM, ALBEDO, UBRDF_FM(IBEAM),     & ! inputs
                DELTA_FACTORM, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,          & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_P(:,IBEAM), ATERM_SAVE, BTERM_SAVE,      & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANSB, ITRANS_USERM(:,IBEAM), T_DELT_USERM(:,IBEAM),    & ! Inputs !@@ 2p3 Greens
                T_DELT_MUBARB, EIGENTRANSNL, LCON, LCON_XVEC1NL,                 & ! inputs
                MCON, MCON_XVEC1NL, WLOWER1NL,                                   & ! inputs
                U_XPOS(IBEAM,:), U_XNEG(IBEAM,:), HMULT_1(IBEAM,:), HMULT_2(IBEAM,:), & ! inputs
                EMULT_UP(:,IBEAM),                                               & ! inputs
                INTENSITY_F_UP(IBEAM), RADLEVEL_F_UP(0:,IBEAM) )                   ! Output !@@ 2p2
         ENDIF

!  Downwelling, MSMODE only,

         IF ( DO_DNWELLING .and. DO_POSTPROCESSING ) THEN
            CALL TWOSTREAM_DNUSER_INTENSITY_SOLAR_OBS &
              ( MAXLAYERS,                                                          & ! Dimension
                DO_2S_LEVELOUT,                                                     & ! Inputs !@@ 2p1, 2p2
                NLAYERS, TAYLOR_ORDER,                                              & ! inputs !@@ 2p3 Greens
                LAYER_PIS_CUTOFFB, DELTA_FACTORM, TAYLOR_SMALL, DELTAU_VERT,        & ! inputs
                GAMMA_P, GAMMA_M, SIGMA_M(:,IBEAM), ATERM_SAVE, BTERM_SAVE,         & ! Inputs !@@ 2p3 Greens
                INITIAL_TRANSB, ITRANS_USERM(:,IBEAM), T_DELT_USERM(:,IBEAM),       & ! Inputs !@@ 2p3 Greens
                T_DELT_MUBARB, LCON, MCON, U_XPOS(IBEAM,:), U_XNEG(IBEAM,:),        & ! Inputs
                HMULT_1(IBEAM,:), HMULT_2(IBEAM,:), EMULT_DN(:,IBEAM),              & ! Inputs
                INTENSITY_F_DN(IBEAM), RADLEVEL_F_DN(0:,IBEAM) )                      ! Output !@@ 2p2
         ENDIF

!  Flux output. New Subroutine, 11/5/13 Version 2.3

         IF ( DO_INCLUDE_MVOUTM ) THEN
           IF (DO_UPWELLING) THEN
               EIGENTRANS1 = EIGENTRANS(1)
               LCON_XVEC21 = LCON(1)*XPOS(2,1)
               MCON_XVEC21 = MCON(1)*XPOS(1,1)
               WUPPER21    = WUPPER(2,1)
            ENDIF
            CALL TWOSTREAM_FLUXES_SOLAR &
              ( DO_UPWELLING, DO_DNWELLING,                                    & ! Flags
                DO_DIRECTBEAM(IBEAM), PI4, STREAM_VALUE,                       & ! Inputs
                FLUX_FACTOR, DELTA_FACTORM, X0(IBEAM), TRANS_SOLAR_BEAM(IBEAM), & ! Inputs
                LCON_XVEC21, MCON_XVEC21, EIGENTRANS1, WUPPER21,               & ! Inputs
                LCON_XVEC1NL, MCON_XVEC1NL, EIGENTRANSNL, WLOWER1NL,           & ! Inputs
                FLUXES_TOA_SOLAR(:,IBEAM), FLUXES_BOA_SOLAR(:,IBEAM) )           ! Outputs
         ENDIF

!  End loop over beam solutions

      END DO

!  ######
!  finish
!  ######

      RETURN
END SUBROUTINE TWOSTREAM_FOURIER_MASTER_SOLAR_OBS

end module twostream_master_solar_obs_m
