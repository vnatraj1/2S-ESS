module twostream_geometry_m

use twostream_inputs_m
use twostream_miscsetups_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_GEOMETRY &
        ( MAXLAYERS, MAXTOTAL, MAXMESSAGES, MAXBEAMS,                     & ! Dimensions
          MAX_USER_RELAZMS, MAX_USER_STREAMS, MAX_USER_OBSGEOMS,          & ! Dimensions !@@ 2p1
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                  & ! Inputs     !@@ 2p2
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs     !@@ 2p3
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                          & ! Inputs
          DO_USER_OBSGEOMS,                                               & ! Inputs     !@@ 2p1   
          NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs     !@@ 2p1
          N_USER_STREAMS, USER_ANGLES, N_USER_RELAZMS, USER_RELAZMS,      & ! Inputs
          NBEAMS, BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,                   & ! Inputs
          N_FOURIERS, DEG_TO_RAD, PI4, UMOFF, DO_POSTPROCESSING, N_PPSTREAMS, & ! Outputs
          PPSTREAM_MASK, CHAPMAN_FACTORS, X0, USER_STREAMS, USER_SECANTS, & ! Outputs
          STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS )           ! Exception handling

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: one = 1.0_dp

!  @@ Rob Spurr, 15 August 2014, Version 2.4, MAXTHREADS dimension removed

      INTEGER, INTENT(IN)          :: MAXMESSAGES, MAXLAYERS, MAXTOTAL
      INTEGER, INTENT(IN)          :: MAXBEAMS
      INTEGER, INTENT(IN)          :: MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS

!  Directional Flags

      LOGICAL, INTENT(IN)          :: DO_UPWELLING, DO_DNWELLING

!  Plane parallel flag

      LOGICAL, INTENT(IN)          :: DO_PLANE_PARALLEL 

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL, INTENT(IN)          :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL, INTENT(IN)          :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  ** New **. October 2011, Sources control, including thermal

      LOGICAL, INTENT(IN)          :: DO_SOLAR_SOURCES
      LOGICAL, INTENT(IN)          :: DO_THERMAL_EMISSION

!  Observational Geometry flag !@@ 2p1

      LOGICAL, INTENT(IN)          :: DO_USER_OBSGEOMS !@@ 2p1

!  Numbers (basic), NTOTAL = 2 * NLAYERS

      INTEGER, INTENT(IN)          :: NLAYERS, NTOTAL

!  Stream value

      REAL(kind=dp), INTENT(IN)    :: STREAM_VALUE

!  Observational geometry input. [Same as LIDORT]. New 12/21/12 !@@ 2p1

      INTEGER, INTENT(IN)          :: N_USER_OBSGEOMS                    !@@ 2p1
      REAL(kind=dp), INTENT(IN)    :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@ 2p1

!  Viewing geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)       :: N_USER_STREAMS
      REAL(kind=dp), INTENT(INOUT) :: USER_ANGLES  ( MAX_USER_STREAMS )
      INTEGER, INTENT(INOUT)       :: N_USER_RELAZMS  
      REAL(kind=dp), INTENT(INOUT) :: USER_RELAZMS ( MAX_USER_RELAZMS )

!  Solar geometry. [Now Intent(inout), thanks to option for ObsGeom !@@ 2p1

      INTEGER, INTENT(INOUT)       :: NBEAMS
      REAL(kind=dp), INTENT(INOUT) :: BEAM_SZAS ( MAXBEAMS )

!  Height and earth radius (latter could be re-set internally)

      REAL(kind=dp), INTENT(INOUT) :: EARTH_RADIUS
      REAL(kind=dp), INTENT(IN)    :: HEIGHT_GRID ( 0:MAXLAYERS ) 

!  Number of Fourier terms

      INTEGER, INTENT(OUT)         :: N_FOURIERS

!  Constants

      REAL(kind=dp), INTENT(OUT)   :: DEG_TO_RAD, PI4

!  Geometry offset array

      INTEGER, INTENT(OUT)         :: UMOFF ( MAX_USER_STREAMS, MAXBEAMS )

!  Post-processing flag (new for Version 2p3)

      LOGICAL, INTENT(OUT)         :: DO_POSTPROCESSING

!  Post-processing control mask

      INTEGER, INTENT(OUT)         :: N_PPSTREAMS, PPSTREAM_MASK ( MAX_USER_STREAMS, MAXBEAMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp), INTENT(OUT)   :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Cosines and sines

      REAL(kind=dp), INTENT(OUT)   :: X0  ( MAXBEAMS )
      REAL(kind=dp), INTENT(OUT)   :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT)   :: USER_SECANTS ( MAX_USER_STREAMS )

!  Exception handling
!  ------------------

!    Check Messages and actions

      INTEGER      , INTENT(OUT) :: STATUS_INPUTCHECK
      INTEGER      , INTENT(OUT) :: C_NMESSAGES
      CHARACTER*100, INTENT(OUT) :: C_MESSAGES(0:MAXMESSAGES)
      CHARACTER*100, INTENT(OUT) :: C_ACTIONS (0:MAXMESSAGES)

!  Local variables

      INTEGER          :: STATUS_SUB
      INTEGER          :: UM, IB, N_VIEWING, IBEAM, IBOFF, I

!  Initialize some variables
!  -------------------------

!  Input check

      STATUS_INPUTCHECK = 0
      C_NMESSAGES       = 0   

!  Check input dimensions
!  ----------------------

      CALL TWOSTREAM_CHECK_INPUT_DIMS &
      ( DO_MVOUT_ONLY, DO_USER_OBSGEOMS, MAXMESSAGES, &
        MAXLAYERS, MAXTOTAL, MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
        NLAYERS,   NTOTAL,   NBEAMS,   N_USER_STREAMS,   N_USER_RELAZMS,   N_USER_OBSGEOMS, &
        STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )

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
!    SS inputs are omitted in this version........
!    !@@ 2p1, Observational Geometry inputs are included (New 12/21/12)
!    !@@ 2p3 Includes check on Flux output flags, and setting of Post-Processing flag

      CALL TWOSTREAM_CHECK_INPUTS_BASIC  &
        ( MAXLAYERS, MAXMESSAGES, MAX_USER_OBSGEOMS,             & ! Dimensions !@@
          MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS,          & ! Dimensions
          DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,         & ! Input
          DO_SOLAR_SOURCES, DO_THERMAL_EMISSION,                 & ! Input
          DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT, DO_POSTPROCESSING, & ! Input !@@ New line, 2p3
          DO_USER_OBSGEOMS, N_USER_OBSGEOMS, USER_OBSGEOMS,      & ! Input !@@ New line
          NLAYERS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,       & ! Input
          BEAM_SZAS, USER_ANGLES, USER_RELAZMS,                  & ! Input
          EARTH_RADIUS, HEIGHT_GRID,                             & ! Input
          STATUS_SUB, C_NMESSAGES, C_MESSAGES, C_ACTIONS )         ! Output

      IF ( STATUS_SUB .EQ. 1 ) THEN
        STATUS_INPUTCHECK = 1
        RETURN
      ENDIF

!  Geometry offsets/masks
!  ======================

!  Save some offsets for indexing geometries

!   !@@ 2p1, This section revised for the Observational Geometry option   
!   !@@ N_GEOMETRIES = NBEAMS * N_USER_STREAMS * N_USER_RELAZMS
!   !@@ 2p3, This section revised for the post-processing flag
!   !@@ 2p4, PPSTREAM masking for Lattice/Obsgeom choice

      N_VIEWING = 0
      UMOFF     = 0
      IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
         N_VIEWING    = N_USER_OBSGEOMS
      ELSE
         if ( DO_POSTPROCESSING ) THEN
            N_VIEWING    = N_USER_STREAMS * N_USER_RELAZMS
            DO IBEAM = 1, NBEAMS
               IBOFF = N_VIEWING * ( IBEAM - 1 )
               DO UM = 1, N_USER_STREAMS
                  UMOFF(UM,IBEAM) = IBOFF +  N_USER_RELAZMS * (UM - 1)
               END DO
            END DO
         ENDIF
      ENDIF

!  Local post-processing control

      PPSTREAM_MASK = 0
      DO IB = 1, NBEAMS
         IF ( DO_USER_OBSGEOMS .AND. DO_SOLAR_SOURCES ) THEN
            N_PPSTREAMS = 1; PPSTREAM_MASK(1,IB) = IB
         else
            N_PPSTREAMS = N_USER_STREAMS
            do UM = 1, N_PPSTREAMS
               PPSTREAM_MASK(UM,IB) = UM
            enddo
         endif
      enddo

!  Chapman function calculation
!  ----------------------------

      if (do_solar_sources) then

      DO IB = 1, NBEAMS
         CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE &
            ( MAXLAYERS,                                & ! Dimensions
              NLAYERS, DO_PLANE_PARALLEL,               & ! Input
              BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID, & ! Input
              CHAPMAN_FACTORS(:,:,IB) )                   ! In/Out
      ENDDO

      endif

!  Get derived inputs
!  ==================  

      if (do_solar_sources) then

!  Solar zenith angle cosine

      DO IB = 1, NBEAMS
         X0(IB) = COS ( BEAM_SZAS(IB) * DEG_TO_RAD )
         IF (DO_PLANE_PARALLEL) AVERAGE_SECANT_PP(IB) = ONE / X0(IB)
      ENDDO

      endif

!  User stream cosines. 11/5/13 2p3 Post-processing control

      IF ( DO_POSTPROCESSING ) THEN  
         DO I = 1, N_USER_STREAMS
            USER_STREAMS(I) = COS(DEG_TO_RAD * USER_ANGLES(I))
            USER_SECANTS(I) = ONE / USER_STREAMS(I)
         ENDDO
      ENDIF

!  Initialise Fourier loop
!  =======================

!  Set Fourier number, Nominally 1 in absence of SS-only flag
!  Zero if no solar sources (Thermal-only run)
!  !@@ 2p3, Set NFOURIERS equal to zero for MVOUT_ONLY

      N_FOURIERS = 1
      IF (  DO_MVOUT_ONLY )         N_FOURIERS = 0   
      IF ( .NOT. DO_SOLAR_SOURCES ) N_FOURIERS = 0

!mick fix 1/7/2012 - (test - make this permanent?)
      IF ( (NBEAMS == 1) .AND. (BEAM_SZAS(1) < 1.0D-8) ) &
        N_FOURIERS = 0

!  Azimuth cosine factors (Fourier = 1). !@@ 2p1, Notice OBSGEOM option
!  !@@ 2p3, not required for FLux-only output

      AZMFAC = zero
      IF ( DO_POSTPROCESSING ) THEN
         IF ( DO_USER_OBSGEOMS.and.DO_SOLAR_SOURCES ) THEN
            DO IB = 1, NBEAMS
               AZM_ARGUMENT = USER_RELAZMS(IB)
               AZMFAC(IB) = COS(DEG_TO_RAD*AZM_ARGUMENT)
            ENDDO
         ELSE
            DO IB = 1, NBEAMS
               DO UM = 1, N_USER_STREAMS
                  DO UA = 1, N_USER_RELAZMS
                     AZM_ARGUMENT = USER_RELAZMS(UA)
                     AZMFAC(UA,UM,IB) = COS(DEG_TO_RAD*AZM_ARGUMENT)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

!  Set flags
!  ---------

!  inclusion of thermal surface emission term, only for Fourier = 0

      DO_INCLUDE_SURFEMISS = .FALSE. 
      IF ( DO_SURFACE_EMISSION ) THEN
        DO_INCLUDE_SURFEMISS(0) = .TRUE.
      ENDIF

!  inclusion of thermal emission term, only for Fourier = 0

      DO_INCLUDE_THERMEMISS = .FALSE.
      IF ( DO_THERMAL_EMISSION ) THEN 
        DO_INCLUDE_THERMEMISS(0) = .TRUE.
      ENDIF

!  Inclusion of mean value calculation
! !@@ 2p3 11/5/13. Control for the Flux calculation

      DO_INCLUDE_MVOUT = .FALSE.
      IF ( DO_ADDITIONAL_MVOUT .OR. DO_MVOUT_ONLY ) THEN
        DO_INCLUDE_MVOUT(0) = .TRUE.
      ENDIF

!  surface reflectance factors

      SURFACE_FACTOR(0) = 2.0_dp
      DELTA_FACTOR(0)   = one
      SURFACE_FACTOR(1) = one
      DELTA_FACTOR(1)   = 2.0_dp

!  Finish

      RETURN

END SUBROUTINE TWOSTREAM_GEOMETRY

SUBROUTINE TWOSTREAM_AUXGEOM &
      ( MAX_USER_STREAMS, MAXBEAMS, DO_POSTPROCESSING,  & ! Dimensions, Flag
        N_USER_STREAMS, NBEAMS,          & ! inputs
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

!  stream directions

      REAL(kind=dp), INTENT(IN)  :: X0 ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

      REAL(kind=dp), INTENT(OUT) :: ULP ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(OUT) :: POX  ( MAXBEAMS )
      REAL(kind=dp), INTENT(OUT) :: PX0X ( MAXBEAMS, 0:1 )
      REAL(kind=dp), INTENT(OUT) :: PXSQ(0:1), PX11

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

      PXSQ(0) = STREAM_VALUE * STREAM_VALUE
      PX11 = DSQRT(0.5d0*(1.0d0-STREAM_VALUE*STREAM_VALUE))
      PXSQ(1) = PX11 * PX11

      DO IBEAM = 1, NBEAMS
        MU0 = X0(IBEAM)
        POXIB = DSQRT(0.5d0*(1.0d0-MU0*MU0))
        POX(IBEAM) = POXIB
        PX0X(IBEAM,0) = MU0 * STREAM_VALUE
        PX0X(IBEAM,1) = POXIB * PX11
      ENDDO

!  FInish

      RETURN
END SUBROUTINE TWOSTREAM_AUXGEOM

end module twostream_geometry_m
