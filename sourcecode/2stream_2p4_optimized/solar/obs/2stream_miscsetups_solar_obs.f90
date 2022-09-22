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
! #              TWOSTREAM_QSPREP_PREPTRANS_EMULT_OBS           #
! #                                                             #
! ###############################################################

module twostream_miscsetups_solar_obs_m

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_QSPREP_PREPTRANS_EMULT_OBS &
       ( MAXLAYERS, MAXGEOMS,                                        & ! Dimensions
         DO_UPWELLING, DO_DNWELLING, DO_POSTPROCESSING,              & ! Flags
         NLAYERS, NGEOMS, TAYLOR_SMALL, TAYLOR_ORDER,                & ! Input
         DELTAU_VERT, CHAPMAN_FACTORS, USER_SECANTS,                 & ! Input
         DO_REFLECTED_DIRECTBEAM,                                    & ! In/Out
         LAYER_PIS_CUTOFF, INITIAL_TRANS, AVERAGE_SECANT,            & ! Output
         TRANS_SOLAR_BEAM, T_DELT_MUBAR, ITRANS_USERM, T_DELT_USERM, & ! Output
         SIGMA_P, SIGMA_M, EMULT_UP, EMULT_DN )                        ! Output

      implicit none

!   precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  Inputs
!  ------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXGEOMS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Control

      INTEGER, INTENT(IN)        :: NLAYERS, NGEOMS

!  Order of Taylor series (including terms up to EPS^n)

      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL
      INTEGER, INTENT(IN)        :: TAYLOR_ORDER

!  optical thickness input

      REAL(kind=dp), INTENT(IN)  :: DELTAU_VERT ( MAXLAYERS )

!  Path segments distances (km)

      REAL(kind=dp), INTENT(IN)  :: CHAPMAN_FACTORS(MAXLAYERS,MAXLAYERS,MAXGEOMS)

!  User secants

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAXGEOMS )

!  Output
!  ------

!  Reflectance flags

      LOGICAL, INTENT(INOUT)     :: DO_REFLECTED_DIRECTBEAM ( MAXGEOMS )

!  Last layer to include Particular integral solution

      INTEGER, INTENT(OUT)       :: LAYER_PIS_CUTOFF(MAXGEOMS)

!  Initial transmittance factors for solar beams and average secant

      REAL(kind=dp), INTENT(OUT) :: INITIAL_TRANS  ( MAXLAYERS, MAXGEOMS )
      REAL(kind=dp), INTENT(OUT) :: AVERAGE_SECANT ( MAXLAYERS, MAXGEOMS )

!  Solar beam attenuations

      REAL(kind=dp), INTENT(OUT) :: TRANS_SOLAR_BEAM ( MAXGEOMS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), INTENT(OUT) :: T_DELT_MUBAR ( MAXLAYERS, MAXGEOMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(OUT) :: ITRANS_USERM ( MAXLAYERS, MAXGEOMS )
      REAL(kind=dp), INTENT(OUT) :: T_DELT_USERM ( MAXLAYERS, MAXGEOMS )

!  coefficient functions for user-defined angles

      REAL(kind=dp), INTENT(OUT) :: SIGMA_P ( MAXLAYERS, MAXGEOMS )
      REAL(kind=dp), INTENT(OUT) :: SIGMA_M ( MAXLAYERS, MAXGEOMS )

!  Forcing term multipliers (saved for whole atmosphere)

      REAL(kind=dp), INTENT(OUT) :: EMULT_UP ( MAXLAYERS, MAXGEOMS )
      REAL(kind=dp), INTENT(OUT) :: EMULT_DN ( MAXLAYERS, MAXGEOMS )   

!  Local variables
!  ---------------

!  Derived optical thickness inputs

      INTEGER       :: N, IB, LAYER_PIS_CUTOFFB
      REAL(kind=dp) :: SB, USIB, S_T_0, S_T_1, TAUSLANT,TAUSLANTN1, DELTAU_VERTN
      REAL(kind=dp) :: INITIAL_TRANS_NIB, SPHER, WDEL, UDEL, WUDEL, ITUDEL
      REAL(kind=dp) :: SIGMA_PN, SIGMA_MN, SU, SD, DIFF

      LOGICAL       :: EMULT_HOPRULE
      
      REAL(kind=dp), PARAMETER :: MAX_TAU_PATH = 88.0d0

    DO IB = 1, NGEOMS

      IF (DO_POSTPROCESSING) USIB = USER_SECANTS(IB)

      S_T_0 = ONE
      S_T_1 = ZERO
      LAYER_PIS_CUTOFFB = NLAYERS
      TAUSLANT = ZERO

!  Average secant, initial transmittance and transmittance factors for average secant streams
!  ==========================================================================================

      DO N = 1, NLAYERS

        DELTAU_VERTN = DELTAU_VERT(N)

!  Get the total spherical attenuation from layer thickness sums and set up the average secant
!  formulation

        TAUSLANTN1 = TAUSLANT
        TAUSLANT = DOT_PRODUCT(DELTAU_VERT(1:N),CHAPMAN_FACTORS(1:N,N,IB))

        IF  ( N .LE. LAYER_PIS_CUTOFFB ) THEN

          IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
            LAYER_PIS_CUTOFFB = N
          ELSE
            S_T_1 = EXP ( - TAUSLANT )
          ENDIF
          SB = (TAUSLANT-TAUSLANTN1) / DELTAU_VERTN
          INITIAL_TRANS_NIB = S_T_0
          S_T_0 = S_T_1

!  Transmittance factors for average secant streams

          SPHER = DELTAU_VERTN * SB
          IF ( SPHER .GT. MAX_TAU_PATH ) THEN
            WDEL = ZERO
          ELSE
            WDEL = EXP ( - SPHER )
          ENDIF

        ELSE

          SB                = ZERO
          INITIAL_TRANS_NIB = ZERO
          WDEL              = ZERO

        ENDIF

        IF (DO_POSTPROCESSING) THEN
          ITUDEL = INITIAL_TRANS_NIB * USIB   
        ENDIF

        AVERAGE_SECANT(N,IB) = SB
        INITIAL_TRANS(N,IB) = INITIAL_TRANS_NIB
        T_DELT_MUBAR(N,IB) = WDEL
        IF (DO_POSTPROCESSING) THEN
          ITRANS_USERM(N,IB) = ITUDEL
        ENDIF

!  Whole Layer transmittances
!  ==========================

        IF (DO_POSTPROCESSING) THEN
          SPHER = DELTAU_VERTN * USIB
          IF ( SPHER.GT.MAX_TAU_PATH ) THEN
            UDEL = ZERO
          ELSE
            UDEL = EXP ( - SPHER )
          ENDIF
          T_DELT_USERM(N,IB) = UDEL
        ENDIF

!  Sigma functions (SIGMA_P/N) and upwelling/downwelling external source function multipliers (EMULT_UP/DN)
!  ========================================================================================================

        IF (DO_POSTPROCESSING) THEN
          IF ( N .LE. LAYER_PIS_CUTOFFB ) THEN
            IF (DO_UPWELLING) THEN
              SIGMA_PN = SB + USIB
              SIGMA_P(N,IB) = SIGMA_PN
              WUDEL = WDEL * UDEL
              SU = ( ONE - WUDEL ) / SIGMA_PN
              EMULT_UP(N,IB) = ITUDEL * SU
            ENDIF
            IF (DO_DNWELLING) THEN
              DIFF = ABS ( USIB - SB )
              IF ( DIFF .LT. TAYLOR_SMALL ) EMULT_HOPRULE = .TRUE.
              SIGMA_MN = SB - USIB
              SIGMA_M(N,IB) = SIGMA_MN
              IF ( EMULT_HOPRULE ) THEN
                CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, SIGMA_MN, DELTAU_VERTN, WDEL, ONE, SD )
              ELSE
                SD = ( UDEL - WDEL ) / SIGMA_MN
              ENDIF
              EMULT_DN(N,IB) = ITUDEL * SD
            ENDIF
          ELSE
            IF (DO_UPWELLING) THEN
              SIGMA_P(N,IB) = ZERO
              EMULT_UP(N,IB) = ZERO
            ENDIF
            IF (DO_DNWELLING) THEN
              SIGMA_M(N,IB) = ZERO
              EMULT_DN(N,IB) = ZERO
            ENDIF
          ENDIF
        ENDIF

      ENDDO


      LAYER_PIS_CUTOFF(IB) = LAYER_PIS_CUTOFFB

!  Set Direct Beam Flag and solar beam total attenuation to surface

      IF ( TAUSLANT .GT. MAX_TAU_PATH ) THEN
        TRANS_SOLAR_BEAM(IB)        = zero
        DO_REFLECTED_DIRECTBEAM(IB) = .FALSE.
      ELSE
        TRANS_SOLAR_BEAM(IB) = EXP( - TAUSLANT )
!mick fix 1/24/12 - Set before subroutine
        !DO_REFLECTED_DIRECTBEAM(IB) = .TRUE.
      ENDIF

!  End geometry loop

    ENDDO

!  debug
!      do N = 1, nlayers
!      write(*,'(i3,1p2e20.10)')n,initial_trans(n,1),average_secant(n,1)
!      enddo

!  debug
!      do n = 1, nlayers
!        write(57,'(i3,1pe24.12)')n,EMULT_UP(N,1)
!        write(57,'(i3,1pe24.12)')n,T_DELT_MUBAR(N,1)
!      enddo

!  debug
!      do n = 1, nlayers
!c        write(57,'(i3,1pe24.12)')n,EMULT_DN(N,1)
!c      enddo


!  finish

      RETURN
END SUBROUTINE TWOSTREAM_QSPREP_PREPTRANS_EMULT_OBS

end module twostream_miscsetups_solar_obs_m
