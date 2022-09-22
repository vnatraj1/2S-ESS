! ############################################################
! #                                                          #
! #             THE TWOSTREAM LIDORT MODEL                   #
! #                                                          #
! #      (LInearized Discrete Ordinate Radiative Transfer)   #
! #       --         -        -        -         -           #
! #                                                          #
! ############################################################

! ############################################################
! #                                                          #
! #  Authors :      Robert. J. D. Spurr (1)                  #
! #                 Vijay Natraj        (2)                  #
! #                                                          #
! #  Address (1) :     RT Solutions, Inc.                    #
! #                    9 Channing Street                     #
! #                    Cambridge, MA 02138, USA              #
! #  Tel:             (617) 492 1183                         #
! #  Email :           rtsolutions@verizon.net               #
! #                                                          #
! #  Address (2) :     CalTech                               #
! #                    Department of Planetary Sciences      #
! #                    1200 East California Boulevard        #
! #                    Pasadena, CA 91125                    #
! #  Tel:             (626) 395 6962                         #
! #  Email :           vijay@gps.caltech.edu                 #
! #                                                          #
! #  Version 1.0-1.3 :                                       #
! #     Mark 1: October  2010                                #
! #     Mark 2: May      2011, with BRDFs                    #
! #     Mark 3: October  2011, with Thermal sources          #
! #                                                          #
! #  Version 2.0-2.1 :                                       #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays   #
! #     Mark 5: December 2012, Observation Geometry option   #
! #                                                          #
! #  Version 2.2-2.3 :                                       #
! #     Mark 6: July     2013, Level outputs + control       #
! #     Mark 7: December 2013, Flux outputs  + control       #
! #     Mark 8: January  2014, Surface Leaving + control     #
! #     Mark 9: June     2014, Inverse Pentadiagonal         #
! #                                                          #
! #  Version 2.4 :                                           #
! #     Mark 10: August  2014, Green's function Regular      #
! #     Mark 11: January 2015, Green's function Linearized   #
! #                            Taylor, dethreaded, OpenMP    #
! #                                                          #
! #  Version 2.5 :                                           #
! #     Mark 12: December 2021, Code optimized for comp. eff.#
! #                             Separate routines for OBS/LAT#
! #                             Wavelength-independent ops   #
! #                             separated                    #
! #                                                          #
! ############################################################

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
! #            TWOSTREAM_CHECK_INPUTS_OPTICAL                   #
! #            TWOSTREAM_DIRECTBEAM                             #
! #            TWOSTREAM_PREPTRANS                              #
! #                                                             #
! ###############################################################

module twostream_miscsetups_m

PUBLIC

contains

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
        IF ( DELTAU_VERT(L) .LE. 0.0_dp ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTION(NM)  = 'Check opt-thickness input '
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for conservative scattering limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L) .GT. 0.999999d0 ) THEN !original
        IF ( OMEGA_TOTAL(L) .GT. 0.999999999d0 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo close to 1, layer '//C3
          ACTION (NM) = 'Check single scattering albedo input'
          STATUS = 1
        ENDIF
      ENDDO

!  check single scatter albedos, for smallness limit

      DO L = 1, NLAYERS
        !IF ( OMEGA_TOTAL(L) .LT. 1.0d-06 ) THEN !original
        IF ( OMEGA_TOTAL(L) .LT. 1.0d-9 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: SS-albedo too small, layer '//C3
          ACTION (NM) = 'Check single scattering albedo input'
          STATUS = 1
        ENDIF
      ENDDO

!  check asymmetry parameter, between -1 and 1

      DO L = 1, NLAYERS
        IF ( ASYMM_TOTAL(L) .LE. -1.0d0 .OR. &
             ASYMM_TOTAL(L) .GE. 1.0d0 ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGE(NM) = 'Bad input: Asymm parameter outside [-1,1], layer'//C3
          ACTION(NM)  = 'Check Asymmetry parameter input'
          STATUS = 1
        ENDIF
      ENDDO

!  set number of messages

      NMESSAGES = NM

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CHECK_INPUTS_OPTICAL

!

SUBROUTINE TWOSTREAM_DIRECTBEAM & 
          ( MAXBEAMS,                                   & ! Dimension
            DO_BRDF_SURFACE, DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, & ! Flags
            NBEAMS, FOURIER, FLUX_FACTOR, X0,           & ! Input
            DELTA_FACTORM, ALBEDO, BRDF_F_0M,           & ! Inputs
            SLTERM_ISOTROPIC, SLTERM_F_0M,              & ! Inputs
            TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,  & ! Inputs
            ATMOS_ATTN, DIRECT_BEAM )                     ! Outputs

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

!  Surface leaving control

      LOGICAL, INTENT(IN)        :: DO_SURFACE_LEAVING
      LOGICAL, INTENT(IN)        :: DO_SL_ISOTROPIC

!  Numbers

      INTEGER, INTENT(IN)        :: NBEAMS

!  Fourier component

      INTEGER      , INTENT(IN)  :: FOURIER

!  FLux factor

      REAL(kind=dp), INTENT(IN)  :: FLUX_FACTOR

!  Solar beams, cosines

      REAL(kind=dp), INTENT(IN)  :: X0(MAXBEAMS)

!  Surface inputs

      REAL(kind=dp), INTENT(IN)  :: DELTA_FACTORM
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_F_0M  ( MAXBEAMS )

!  Do not need this, MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: UBRDF_F_0M ( MAXBEAMS )

!  SLEAVE stuff
!  ** Isotropic Surface leaving term (if flag set)
!  ** Fourier components of Surface-leaving terms:
!     Every solar direction, SL-transmitted quadrature streams

      REAL(kind=dp), INTENT(IN) :: SLTERM_ISOTROPIC ( MAXBEAMS )
      REAL(kind=dp), INTENT(IN) :: SLTERM_F_0M ( MAXBEAMS )

!  Solar beam attenuations and reflectance flags

      REAL(kind=dp), INTENT(IN)  :: TRANS_SOLAR_BEAM        ( MAXBEAMS )
      LOGICAL, INTENT(IN)        :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )

!  output arguments
!  ----------------

!  Atmospheric attenuation

      REAL(kind=dp), INTENT(OUT) :: ATMOS_ATTN ( MAXBEAMS )

!  Direct beam solution, do not need USER_DIRECT_BEAM value

      REAL(kind=dp), INTENT(OUT) :: DIRECT_BEAM      ( MAXBEAMS )
!      REAL(kind=dp), INTENT(OUT) :: USER_DIRECT_BEAM ( MAX_USER_STREAMS, MAXBEAMS )

!  Local variables
!  ---------------

      REAL(kind=dp) :: PI, ATTN, HELP
      INTEGER       :: IB

!  Attenuation of solar beam
!  -------------------------

      PI = ACOS(-ONE)
      DO IB = 1, NBEAMS
       IF ( DO_REFLECTED_DIRECTBEAM(IB) ) THEN

        ATTN           = FLUX_FACTOR * X0(IB) / DELTA_FACTORM / PI * TRANS_SOLAR_BEAM(IB)
        ATMOS_ATTN(IB) = ATTN

!  Total contributions, BRDF or Lambertian

        IF ( DO_BRDF_SURFACE ) THEN
          DIRECT_BEAM(IB) = ATTN * BRDF_F_0M(IB)
        ELSE
          DIRECT_BEAM(IB) = ATTN * ALBEDO
        ENDIF

!  Surface-Leaving stuff
!  Normalized to Flux-factor / DELTA_Factor
!  Delta_Factor = 1.0 for the Isotropic or non-iso Fourier = 0 cases

        IF ( DO_SURFACE_LEAVING ) THEN
          HELP = FLUX_FACTOR / DELTA_FACTORM
          IF ( DO_SL_ISOTROPIC .AND. FOURIER .EQ. 0 ) THEN
            DIRECT_BEAM(IB) = DIRECT_BEAM(IB) + SLTERM_ISOTROPIC(IB) * HELP
          ELSE
            DIRECT_BEAM(IB) = DIRECT_BEAM(IB) + SLTERM_F_0M(IB) * HELP
          ENDIF
        ENDIF

!  end direct beam calculation

       ENDIF
      ENDDO

!  finish

      RETURN
END SUBROUTINE TWOSTREAM_DIRECTBEAM

!

SUBROUTINE TWOSTREAM_PREPTRANS &
       ( MAXLAYERS, MAX_USER_STREAMS, & ! Dimensions
         NLAYERS, N_USER_STREAMS,     & ! Input
         DELTAU_VERT, USER_SECANTS,   & ! Input
         T_DELT_USERM )                 ! Output

      implicit none

!   precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp 

!  Inputs
!  ------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Control

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  optical thickness input

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT ( MAXLAYERS )  

!  User secants

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Output
!  ------

!  Transmittance factors for user-defined stream angles
      
      REAL(kind=dp), INTENT(OUT) :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Local variables
!  ---------------
      
!  Derived optical thickness inputs

      INTEGER       :: N, UM
      REAL(kind=dp) :: USUM, SPHER, UDEL

      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

    DO UM = 1, N_USER_STREAMS
      
      USUM = USER_SECANTS(UM)

      DO N = 1, NLAYERS

!  Whole Layer transmittances
!  ==========================

         SPHER = DELTAU_VERT(N) * USUM
         IF ( SPHER .GT. MAX_TAU_PATH ) THEN
            UDEL = ZERO
         ELSE
            UDEL = EXP ( - SPHER )
         ENDIF
         T_DELT_USERM(N,UM) = UDEL

      ENDDO

!  End geometry loop
    ENDDO

      RETURN
END SUBROUTINE TWOSTREAM_PREPTRANS

!

end module twostream_miscsetups_m
