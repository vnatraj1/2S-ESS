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
! #            TWOSTREAM_BRDFMASTER_SOLAR_OBS (master), calling #
! #                                                             #
! #              TWOSTREAM_BRDF_MAKER_SOLAR_OBS, calling        #
! #                TWOSTREAM_BRDF_FUNCTION                      #
! #              TWOSTREAM_BRDF_FOURIER_SOLAR_OBS               #
! #              TWOSTREAM_GAULEG                               #
! #                                                             #
! ###############################################################

module twostream_brdf_supplement_solar_obs_m

!  Module for calling the 2S brdf supplement, BRDFs only

!  This contains all the necessary pieces

!  Version 2.0.
!     Construction October 21, 2011
!     R. Spurr, RT SOLUTIONS Inc.
!     Upgraded by M. Christi, November 2012

!  Version 2.1. Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Exception handling (dimension checks)
!     User Relazimuths are not required, so have been removed

!  Version 2.4
!    TWOSTREAM_BRDF_FUNCTION --> Upgrade to BPDF kernels (SOIL, VEGN, NDVI)

USE twostream_brdfkernels_m

CONTAINS

!

SUBROUTINE TWOSTREAM_BRDFMASTER_SOLAR_OBS &
     ( MAX_USER_OBSGEOMS, MAXSTREAMS_BRDF,            & ! Dimensions
       MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS,         & ! Dimensions
       LAMBERTIAN_KERNEL_FLAG, DO_SHADOW_EFFECT,      & ! Inputs
       N_USER_OBSGEOMS,                               & ! Inputs
       STREAM_VALUE, NSTREAMS_BRDF,                   & ! Inputs
       N_BRDF_KERNELS, WHICH_BRDF, BRDF_FACTORS,      & ! Inputs
       N_BRDF_PARAMETERS, BRDF_PARAMETERS,            & ! Inputs
       BRDF_F_0, BRDF_F, UBRDF_F,                     & ! Outputs
       STATUS_BRDFSUP, MESSAGE, ACTION )                ! Outputs

!  Prepares the bidirectional reflectance functions
!  necessary for 2S code (MULTIPLE SCATTERING ONLY)

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)       :: MAX_USER_OBSGEOMS
      INTEGER, INTENT(IN)       :: MAXSTREAMS_BRDF, MAX_BRDF_KERNELS, &
                                   MAX_BRDF_PARAMETERS

!  Lambertian surface control

      LOGICAL, INTENT(IN)       :: LAMBERTIAN_KERNEL_FLAG (MAX_BRDF_KERNELS)

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, INTENT(IN)       :: DO_SHADOW_EFFECT

!  Observational geometry input. [Same as LIDORT].

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS

!  2-Stream angle cosine

      REAL(kind=dp), INTENT(IN) :: STREAM_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER, INTENT(IN)       :: NSTREAMS_BRDF

!  Number and index-list of bidirectional functions
!  Lambertian Surface control

      INTEGER, INTENT(IN)       :: N_BRDF_KERNELS
      INTEGER, INTENT(IN)       :: WHICH_BRDF (MAX_BRDF_KERNELS)

!  kernel amplitude factors and Parameters required for Kernel families

      REAL(kind=dp), INTENT(IN) :: BRDF_FACTORS (MAX_BRDF_KERNELS)
      INTEGER, INTENT(IN)       :: N_BRDF_PARAMETERS (MAX_BRDF_KERNELS)
      REAL(kind=dp), INTENT(IN) :: BRDF_PARAMETERS (MAX_BRDF_PARAMETERS,MAX_BRDF_KERNELS)

!  Output arguments
!  ================

!  BRDF Fourier components (NOT threaded)
!  0 and 1 Fourier components of BRDF, following order (same all threads)
!  incident solar directions,  reflected quadrature stream
!  incident quadrature stream, reflected quadrature stream
!  incident solar directions,  reflected user streams -- NOT REQUIRED
!  incident quadrature stream, reflected user streams

      REAL(kind=dp), INTENT(OUT)  :: BRDF_F_0  ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp), INTENT(OUT)  :: BRDF_F    ( 0:1 )
!      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F_0 ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp), INTENT(OUT)  :: UBRDF_F   ( MAX_USER_OBSGEOMS, 0:1 )

!  Exception handling

      INTEGER      , intent(out) :: status_brdfsup
      CHARACTER*(*), intent(out) :: message, action

!  Local BRDF functions
!  ====================

!  At quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: BRDFUNC_0 ( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  At user-defined stream directions

      REAL(kind=dp)  :: USER_BRDFUNC   ( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  Local angles, and cosine/sines/weights
!  ======================================

!  SZAs

      REAL(kind=dp)  :: SZASURCOS(MAX_USER_OBSGEOMS)
      REAL(kind=dp)  :: SZASURSIN(MAX_USER_OBSGEOMS)

!  Viewing zenith streams

      REAL(kind=dp)  :: USER_STREAMS(MAX_USER_OBSGEOMS)
      REAL(kind=dp)  :: USER_SINES  (MAX_USER_OBSGEOMS)

!  BRDF azimuth quadrature streams

      REAL(kind=dp)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  Azimuth factors

      REAL(kind=dp)  :: BRDF_AZMFAC(MAXSTREAMS_BRDF)

!  Local kernel Fourier components
!  ===============================

!  At quadrature (discrete ordinate) angles

      REAL(kind=dp)  :: LOCAL_BRDF_F
      REAL(kind=dp)  :: LOCAL_BRDF_F_0 ( MAX_USER_OBSGEOMS )

!  At user-defined stream directions

      REAL(kind=dp)  :: LOCAL_USER_BRDF_F (MAX_USER_OBSGEOMS )

!  Other local variables
!  =====================

!  Spherical albedo

!      REAL(kind=dp)  :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  Help

      INTEGER        :: K, B, IB, UM, M, I
      INTEGER        :: LOCAL_BRDF_NPARS
      REAL(kind=dp)  :: LOCAL_BRDF_PARS(MAX_BRDF_KERNELS), LFAC
      REAL(kind=dp)  :: DELFAC, STREAM_SINE
      LOGICAL        :: ADD_FOURIER

      INTEGER, PARAMETER :: COXMUNK_IDX = 9

!  Initialize
!  ----------

!  Exception handling

      STATUS_BRDFSUP = 0
      MESSAGE = ' '
      ACTION  = ' '

!  Initialise BRDF arrays (IMPORTANT)
!  ----------------------------------

      BRDF_F_0        = 0.0_dp
      BRDF_F          = 0.0_dp
      UBRDF_F         = 0.0_dp

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Local variables

        LOCAL_BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, LOCAL_BRDF_NPARS
          LOCAL_BRDF_PARS(B) = BRDF_PARAMETERS(B,K)
        ENDDO
        LFAC = BRDF_FACTORS(K)

!  Coxmunk shadow flag

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) LOCAL_BRDF_PARS(3) = 1.0_dp
        ENDIF

!  Get the kernels. Solar sources optionality, added 12/31/12

        CALL twostream_brdfmaker_solar_obs &
           ( MAX_USER_OBSGEOMS, MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS,   & ! Dimensions
             WHICH_BRDF(K), LOCAL_BRDF_NPARS, LOCAL_BRDF_PARS,          & ! Inputs
             NSTREAMS_BRDF, N_USER_OBSGEOMS, STREAM_VALUE, STREAM_SINE, & ! Inputs
             USER_STREAMS, USER_SINES, SZASURCOS, SZASURSIN,            & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF,                                  & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0  )                          ! Outputs

!  Two Fourier components

        DO M = 0, 1

!  Fourier addition flag

          ADD_FOURIER = (  .NOT. LAMBERTIAN_KERNEL_FLAG(K) .OR. &
                          (LAMBERTIAN_KERNEL_FLAG(K) .AND. M .EQ. 0) )

!  Surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = 1.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I)
            ENDDO
          ELSE
            DELFAC   = 2.0_dp
            DO I = 1, NSTREAMS_BRDF
              BRDF_AZMFAC(I) = A_BRDF(I) * COS ( M * X_BRDF(I) )
            ENDDO
          ENDIF

!  Call. Solar sources optionality, added 12/31/12

          CALL twostream_brdf_fourier_solar_obs &
           ( MAX_USER_OBSGEOMS, MAXSTREAMS_BRDF,        & ! Dimensions
             LAMBERTIAN_KERNEL_FLAG(K),                 & ! Inputs
             M, DELFAC, N_USER_OBSGEOMS, NSTREAMS_BRDF, & ! Inputs
             BRDFUNC, USER_BRDFUNC, BRDFUNC_0,          & ! Inputs
             BRDF_AZMFAC,                               & ! Inputs
             LOCAL_BRDF_F, LOCAL_BRDF_F_0,              & ! Outputs
             LOCAL_USER_BRDF_F )                          ! Outputs

!  Spherical albedo (debug only)

!          IF ( M .EQ. 0 ) THEN
!            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
!              HELP_A  = 4.0_dp * LOCAL_BRDF_F * STREAM_VALUE * STREAM_VALUE
!              SPHERICAL_ALBEDO(K) = HELP_A
!            ENDIF
!          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature-quadrature reflectance)

            BRDF_F(M) = BRDF_F(M) + LFAC * LOCAL_BRDF_F

!  Kernel combinations (for Solar-quadrature reflectance)

            DO IB = 1, N_USER_OBSGEOMS
              BRDF_F_0(IB,M) = BRDF_F_0(IB,M) + LFAC * LOCAL_BRDF_F_0(IB)
            ENDDO

!  Kernel combinations (for quadrature-userstream reflectance)

            DO UM = 1, N_USER_OBSGEOMS
              UBRDF_F(UM,M) = UBRDF_F(UM,M) + LFAC * LOCAL_USER_BRDF_F(UM)
            ENDDO

!  End Fourier addition clause and Fourier loop

          ENDIF
        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      return
end subroutine twostream_brdfmaster_solar_obs

!

subroutine twostream_brdfmaker_solar_obs &
    ( MAX_USER_OBSGEOMS, MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS, & ! Dimensions
      WHICH_BRDF, BRDF_NPARS, BRDF_PARS,                       & ! Inputs
      NSTREAMS_BRDF, N_USER_OBSGEOMS,                          & ! Inputs
      STREAM_VALUE, STREAM_SINE, USER_STREAMS, USER_SINES,     & ! Inputs
      SZAC, SZAS, X_BRDF, CX_BRDF, SX_BRDF,                    & ! Inputs
      BRDFUNC, USER_BRDFUNC, BRDFUNC_0 )                         ! Outputs

!  Prepares the bidirectional reflectance scatter matrices
!  Solar sources optionality, added 12/31/12

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, intent(in)  :: MAX_USER_OBSGEOMS
      INTEGER, intent(in)  :: MAXSTREAMS_BRDF, MAX_BRDF_PARAMETERS

!  Which BRDF index

      INTEGER, intent(in)  :: WHICH_BRDF

!  Local number of parameters and local parameter array

      INTEGER, intent(in)        :: BRDF_NPARS
      REAL(kind=dp), intent(in)  :: BRDF_PARS (MAX_BRDF_PARAMETERS)

!  Local angle control

      INTEGER  , intent(in)      :: N_USER_OBSGEOMS

!  Local angles

      REAL(kind=dp), intent(in)  :: SZAC(MAX_USER_OBSGEOMS)
      REAL(kind=dp), intent(in)  :: SZAS(MAX_USER_OBSGEOMS)

      REAL(kind=dp), intent(in)  :: STREAM_VALUE
      REAL(kind=dp), intent(in)  :: STREAM_SINE

      REAL(kind=dp), intent(in)  :: USER_STREAMS(MAX_USER_OBSGEOMS)
      REAL(kind=dp), intent(in)  :: USER_SINES  (MAX_USER_OBSGEOMS)

!  Azimuth quadrature streams for BRDF

      INTEGER  , intent(in)  :: NSTREAMS_BRDF
      REAL(kind=dp), intent(in)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  At quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(out) :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(out) :: BRDFUNC_0 ( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  At user-defined stream directions

      REAL(kind=dp), intent(out) :: USER_BRDFUNC( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  Local variables
!  ---------------

      INTEGER   :: UI, K, IB
      REAL(kind=dp) :: SZACIB, SZASIB, USER_STREAMSUI, USER_SINESUI

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, N_USER_OBSGEOMS
        SZACIB = SZAC(IB)
        SZASIB = SZAS(IB)
        DO K = 1, NSTREAMS_BRDF
          call twostream_brdf_function &
           ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
             SZACIB, SZASIB, STREAM_VALUE, STREAM_SINE,              & ! Inputs
             X_BRDF(K), CX_BRDF(K), SX_BRDF(K),                      & ! Inputs
             BRDFUNC_0(K,IB) )
        ENDDO
      ENDDO

!  Incident quadrature directions

      DO K = 1, NSTREAMS_BRDF
         call twostream_brdf_function &
           ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
             STREAM_VALUE, STREAM_SINE, STREAM_VALUE,                & ! Inputs
             STREAM_SINE, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         & ! Inputs
             BRDFUNC(K) )
      ENDDO

!  User-streams outgoing direction

        DO UI = 1, N_USER_OBSGEOMS
           USER_STREAMSUI = USER_STREAMS(UI)
           USER_SINESUI =  USER_SINES(UI)
           DO K = 1, NSTREAMS_BRDF
              call twostream_brdf_function &
                 ( MAX_BRDF_PARAMETERS, WHICH_BRDF, BRDF_NPARS, BRDF_PARS, & ! Inputs
                   STREAM_VALUE, STREAM_SINE, USER_STREAMSUI,              & ! Inputs
                   USER_SINESUI, X_BRDF(K), CX_BRDF(K), SX_BRDF(K),        & ! Inputs
                   USER_BRDFUNC(K,UI) )
           ENDDO
        ENDDO

!  Finish

      RETURN
end subroutine twostream_brdfmaker_solar_obs

!

subroutine twostream_brdf_function &
   ( MAXPARS, WHICH_BRDF, NPARS, PARS,   & ! Inputs
     XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI, & ! Inputs
     KERNEL )

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  indices

!  These refer to the BRDF kernel functions currently included.
!    Updated list of kernels (Version 2.4)

      INTEGER, PARAMETER :: LAMBERTIAN_IDX  = 1
      INTEGER, PARAMETER :: ROSSTHIN_IDX    = 2
      INTEGER, PARAMETER :: ROSSTHICK_IDX   = 3
      INTEGER, PARAMETER :: LISPARSE_IDX    = 4
      INTEGER, PARAMETER :: LIDENSE_IDX     = 5
      INTEGER, PARAMETER :: HAPKE_IDX       = 6
      INTEGER, PARAMETER :: ROUJEAN_IDX     = 7
      INTEGER, PARAMETER :: RAHMAN_IDX      = 8
      INTEGER, PARAMETER :: COXMUNK_IDX     = 9

      INTEGER, PARAMETER :: BPDFVEGN_IDX    = 10
      INTEGER, PARAMETER :: BPDFSOIL_IDX    = 11
      INTEGER, PARAMETER :: BPDFNDVI_IDX    = 12

      INTEGER, PARAMETER :: MAXBRDF_IDX = BPDFNDVI_IDX

!  Subroutine arguments

      INTEGER      , intent(in)  :: MAXPARS
      INTEGER      , intent(in)  :: WHICH_BRDF
      INTEGER      , intent(in)  :: NPARS
      REAL(kind=dp), intent(in)  :: PARS (MAXPARS)
      REAL(kind=dp), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
      REAL(kind=dp), intent(out) :: KERNEL

!  Trawl through

      IF ( WHICH_BRDF .EQ. LAMBERTIAN_IDX ) THEN
        CALL TWOSTREAM_LAMBERTIAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHIN_IDX ) THEN
        CALL TWOSTREAM_ROSSTHIN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROSSTHICK_IDX ) THEN
        CALL TWOSTREAM_ROSSTHICK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LISPARSE_IDX ) THEN
        CALL TWOSTREAM_LISPARSE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. LIDENSE_IDX ) THEN
        CALL TWOSTREAM_LIDENSE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. RAHMAN_IDX ) THEN
        CALL TWOSTREAM_RAHMAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. ROUJEAN_IDX ) THEN
        CALL TWOSTREAM_ROUJEAN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. HAPKE_IDX ) THEN
        CALL TWOSTREAM_HAPKE_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. COXMUNK_IDX ) THEN
        CALL TWOSTREAM_COXMUNK_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFVEGN_IDX ) THEN
        CALL TWOSTREAM_BPDFVEGN_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFSOIL_IDX ) THEN
        CALL TWOSTREAM_BPDFSOIL_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ELSE IF ( WHICH_BRDF .EQ. BPDFNDVI_IDX ) THEN
        CALL TWOSTREAM_BPDFNDVI_FUNCTION &
        ( MAXPARS, NPARS, PARS, &
          XJ, SXJ, XI, SXI, PHI, CPHI, &
          SKPHI, KERNEL )
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdf_function

!

subroutine twostream_brdf_fourier_solar_obs &
         ( MAX_USER_OBSGEOMS, MAXSTREAMS_BRDF, & ! Dimensions
           LAMBERTIAN_FLAG, M, DELFAC,         & ! Inputs
           N_USER_OBSGEOMS, NSTREAMS_BRDF,     & ! Inputs
           BRDFUNC, USER_BRDFUNC, BRDFUNC_0,   & ! Inputs
           BRDF_AZMFAC,                        & ! Inputs
           LOCAL_BRDF_F, LOCAL_BRDF_F_0,       & ! Outputs
           LOCAL_USER_BRDF_F )                   ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions
!  Solar sources optionality, added 12/31/12

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimensions

      INTEGER  , intent(in)  :: MAX_USER_OBSGEOMS
      INTEGER  , intent(in)  :: MAXSTREAMS_BRDF

!  Control

      LOGICAL      , intent(in)  :: LAMBERTIAN_FLAG
      REAL(kind=dp), intent(in)  :: DELFAC
      INTEGER      , intent(in)  :: M

!  Local numbers

      INTEGER  , intent(in)  :: N_USER_OBSGEOMS
      INTEGER  , intent(in)  :: NSTREAMS_BRDF

!  Azimuth cosines

      REAL(kind=dp), intent(in)  :: BRDF_AZMFAC ( MAXSTREAMS_BRDF )

!  at quadrature (discrete ordinate) angles

      REAL(kind=dp), intent(in)  :: BRDFUNC   ( MAXSTREAMS_BRDF )
      REAL(kind=dp), intent(in)  :: BRDFUNC_0 ( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  at user-defined stream directions

      REAL(kind=dp), intent(in)  :: USER_BRDFUNC ( MAXSTREAMS_BRDF, MAX_USER_OBSGEOMS )

!  Output: Local kernel Fourier components
!  =======================================

      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F
      REAL(kind=dp), intent(out) :: LOCAL_BRDF_F_0    ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), intent(out) :: LOCAL_USER_BRDF_F ( MAX_USER_OBSGEOMS )

!  local variables
!  ===============

      INTEGER        :: UI, IB
      REAL(kind=dp)  :: SUM1, HELP

!  surface factor

      HELP = 0.5_dp * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO IB = 1, N_USER_OBSGEOMS
          SUM1 = DOT_PRODUCT(BRDFUNC_0(1:NSTREAMS_BRDF,IB),BRDF_AZMFAC(1:NSTREAMS_BRDF))
          LOCAL_BRDF_F_0(IB) = SUM1 * HELP
        ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
        DO IB = 1, N_USER_OBSGEOMS
          LOCAL_BRDF_F_0(IB) = 1.0_dp
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         SUM1 = DOT_PRODUCT(BRDFUNC(1:NSTREAMS_BRDF),BRDF_AZMFAC(1:NSTREAMS_BRDF))
         LOCAL_BRDF_F = SUM1 * HELP
      ELSE IF ( M .EQ. 0 ) THEN
         LOCAL_BRDF_F = 1.0_dp
      ENDIF

!  User-streams incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
         DO UI = 1, N_USER_OBSGEOMS
            SUM1 = DOT_PRODUCT(USER_BRDFUNC(1:NSTREAMS_BRDF,UI),BRDF_AZMFAC(1:NSTREAMS_BRDF))
            LOCAL_USER_BRDF_F(UI) = SUM1 * HELP
         ENDDO
      ELSE IF ( M .EQ. 0 ) THEN
         DO UI = 1, N_USER_OBSGEOMS
            LOCAL_USER_BRDF_F(UI) = 1.0_dp
         ENDDO
      ENDIF

!  Finish

      RETURN
end subroutine twostream_brdf_fourier_solar_obs

!  end

SUBROUTINE TWOSTREAM_GAULEG(X1,X2,X,W,N)

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Arguments

      INTEGER       :: N
      REAL(kind=dp) :: X1,X2,X(N),W(N)

!  Local variables

      INTEGER       :: I, M, J 
      REAL(kind=dp) :: XM,XL,P1,P2,P3,PP,Z,Z1

      REAL(kind=dp), PARAMETER :: EPS = 3.D-14

      M =  (N+1)/2
      XM = 0.5D0*(X2+X1)
      XL = 0.5D0*(X2-X1)

      DO I=1,M
        Z = DCOS ( 3.141592654D0 * (I-.25D0) / (N+.5D0) )
 1      CONTINUE
        P1 = 1.D0
        P2 = 0.D0
        DO J = 1, N
          P3 = P2
          P2 = P1
          P1 = ( (2.D0*J-1.D0) * Z *P2 - (J-1.D0) * P3 ) / J
        ENDDO
        PP = N * (Z*P1-P2) / (Z*Z-1.D0)
        Z1 = Z
        Z  = Z1 - P1 / PP
        IF ( DABS(Z-Z1) .GT. EPS ) GOTO 1
        X(I)     = XM - XL * Z
        X(N+1-I) = XM + XL * Z
        W(I)     = 2.D0 * XL / ( (1.D0-Z*Z) * PP * PP )
        W(N+1-I) = W(I)
      ENDDO

      RETURN
END SUBROUTINE TWOSTREAM_GAULEG

end module twostream_brdf_supplement_solar_obs_m
