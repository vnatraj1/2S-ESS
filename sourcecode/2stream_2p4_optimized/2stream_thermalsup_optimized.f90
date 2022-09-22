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
! #   setups                                                    #
! #                                                             #
! #          TWOSTREAM_THERMALSETUP                             #
! #                                                             #
! #   discrete ordinate particular integral                     #
! #                                                             #
! #          TWOSTREAM_THERMALGFSOLUTION                        #
! #                                                             #
! #   postprocessing source terms                               #
! #                                                             #
! #          TWOSTREAM_THERMALSTERMS                            #
! #                                                             #
! ###############################################################

module twostream_thermalsup_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_THERMALSETUP &
          ( MAXLAYERS, NLAYERS,            & ! Inputs
            DELTAU_VERT, THERMAL_BB_INPUT, & ! inputs
            THERMCOEFFS )                    ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Number of layers

      INTEGER, INTENT(IN)        :: NLAYERS

!  Optical properties

      REAL(kind=dp), INTENT (IN) :: DELTAU_VERT ( MAXLAYERS )

!  Thermal input

      REAL(kind=dp), INTENT (IN) :: THERMAL_BB_INPUT ( 0:MAXLAYERS )

!  Output Help variables

      REAL (kind=dp), INTENT (OUT) :: THERMCOEFFS  ( 2, MAXLAYERS )

!  Local variables

      INTEGER        :: N
      REAL (kind=dp) :: TCN1, TCN, TC2

!  thermal coefficients and powers and auxiliary output

      TCN1 = THERMAL_BB_INPUT(0)
      DO N = 1, NLAYERS
        TCN = THERMAL_BB_INPUT(N)
        THERMCOEFFS(1,N)  = TCN1
        TC2 = (TCN-TCN1)/DELTAU_VERT(N)
        THERMCOEFFS(2,N)  = TC2
        TCN1 = TCN
      END DO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSETUP


SUBROUTINE TWOSTREAM_THERMALGFSOLUTION &
      ( MAXLAYERS, NLAYERS, OMEGA, DELTAUS, THERMCOEFFS,      & ! Inputs
        TCUTOFF, EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED,    & ! Input
        T_C_PLUS, T_C_MINUS, TTERM_SAVE, T_WUPPER, T_WLOWER )   ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS

!  OMEGA and DELTAUS

      REAL(kind=dp), INTENT(IN)   :: OMEGA ( MAXLAYERS )
      REAL(kind=dp), INTENT (IN)  :: DELTAUS ( MAXLAYERS )

!  Help variable

      REAL(kind=dp), INTENT (IN) :: THERMCOEFFS ( 2, MAXLAYERS )

!  Thermal Cutoff (actually a layer optical thickness minimum)
!     Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!    Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Transmittance factors and eigenvalues

      REAL(kind=dp), intent(in)  :: EIGENVALUE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: EIGENTRANS (MAXLAYERS)

!  Eigenvector solutions

      REAL(kind=dp), intent(in)  :: XPOS (2,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: NORM_SAVED (MAXLAYERS)

!  Output variables
!  ----------------

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(out) :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(out) :: T_C_MINUS (0:2,MAXLAYERS)
      REAL(kind=dp), intent(out) :: T_C_PLUS  (0:2,MAXLAYERS)

!  Thermal solution at layer boundaries

      REAL(kind=dp) , INTENT (OUT) ::  T_WUPPER ( 2, MAXLAYERS )
      REAL(kind=dp) , INTENT (OUT) ::  T_WLOWER ( 2, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: N
      REAL(kind=dp) :: TK, K1, TT, t_gmult_dn, t_gmult_up
      REAL(kind=dp) :: omega1, help, sum_m, sum_p
      REAL(kind=dp) :: deltausn, xpos1n, xpos2n
      REAL(kind=dp) :: tc1n, tc2n, tcm2n, tcp2n, tcm1n, tcp1n, tcm0n, tcp0n

! ---------------------------------------
! Green function solutions for all layers
! ---------------------------------------

   DO n = 1, nlayers

      deltausn = deltaus(n)

!  Only do this if cutoff is not in play

      if ( deltausn .gt. TCUTOFF )  then

        xpos1n = xpos(1,n)
        xpos2n = xpos(2,n)
        tc1n   = thermcoeffs(1,n)
        tc2n   = thermcoeffs(2,n)

!  get the constant term

        omega1 = one - omega(n)
        help   = xpos1n+xpos2n
        tterm_save(n) = omega1 * help / norm_saved(n)
        tt = tterm_save(n)

! Green function multipliers

        k1 = one / eigenvalue(n)
        tk = eigentrans(n)
        tcm2n = k1 * tc2n
        tcp2n = tcm2n
        t_c_minus(2,n)  = tcm2n
        t_c_plus(2,n)   = tcp2n
        tcm1n = k1 * ( tc1n - tcm2n )
        tcp1n = k1 * ( tc1n + tcp2n )
        t_c_minus(1,n) = tcm1n
        t_c_plus(1,n)  = tcp1n
        sum_m = tcm1n + tcm2n * deltausn
        sum_p = tcp1n + tcp2n * deltausn
        tcm0n = - tcm1n
        tcp0n = - sum_p
        t_c_minus(0,n) = tcm0n
        t_c_plus(0,n)  = tcp0n
        t_gmult_dn = tt * ( tk * tcm0n + sum_m )
        t_gmult_up = tt * ( tk * tcp0n + tcp1n )

! Set particular integral from Green function expansion

        t_wupper(1,n) = t_gmult_up*xpos2n
        t_wupper(2,n) = t_gmult_up*xpos1n
        t_wlower(2,n) = t_gmult_dn*xpos2n
        t_wlower(1,n) = t_gmult_dn*xpos1n

! Zero boundary layer values if deltau below cutoff

      else

        t_wupper(1:2,n) = zero
        t_wlower(1:2,n) = zero

!  End active layer condition

      ENDIF

!  End layer loop

   END DO

!  finish

   RETURN
END SUBROUTINE TWOSTREAM_THERMALGFSOLUTION

!

SUBROUTINE TWOSTREAM_THERMALSTERMS_THERMAL &
      ( MAXLAYERS, MAX_USER_STREAMS,           & ! Dimensions
        DO_UPWELLING, DO_DNWELLING,            & ! Inputs
        NLAYERS, N_USER_STREAMS, USER_STREAMS, & ! Inputs
        TCUTOFF, T_DELT_USERM, DELTAUS,        & ! Inputs
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,      & ! Inputs
        T_C_PLUS, T_C_MINUS, TTERM_SAVE,       & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN  )          ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAX_USER_STREAMS

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT (IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Thermal Cutoff (actually a layer optical thickness minimum)
!  Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!  Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Optical thickness

      REAL(kind=dp), INTENT (IN)  :: DELTAUS ( MAXLAYERS )

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)   :: U_XPOS (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)   :: U_XNEG (MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (global, whole layer)

      REAL(kind=dp), intent(in)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_MINUS (0:2,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_PLUS  (0:2,MAXLAYERS)

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: SUM_M, SUM_P, SU, SD, COSMUM, DELTAUSN
      REAL(kind=dp) :: TSGM_DU0N, TSGM_DU1N, TSGM_DU2N
      REAL(kind=dp) :: TSGM_DD0N, TSGM_DD1N, TSGM_DD2N
      REAL(kind=dp) :: TSGM_UU0N, TSGM_UU1N, TSGM_UU2N
      REAL(kind=dp) :: TSGM_UD0N, TSGM_UD1N, TSGM_UD2N
      REAL(kind=dp) :: TCP0N, TCP1N, TCM0N, TCM1N
      REAL(kind=dp) :: TTERM_SAVEN, T_DELT_USERM_NUM

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  Rob  Fix, 5.14.15. Only get solutions for valid non-Cutoff layers.

!  Upwelling

    IF ( DO_UPWELLING ) THEN
      DO N = 1, NLAYERS
        deltausn = deltaus(n)
        IF ( deltausn .gt. TCUTOFF ) THEN
          tsgm_uu2n = t_c_plus (2,n)
          tsgm_ud2n = t_c_minus(2,n)
          tcp1n = t_c_plus (1,n)
          tcm1n = t_c_minus(1,n)
          tcp0n = t_c_plus (0,n)
          tcm0n = t_c_minus(0,n)
          tterm_saven = tterm_save(n)
          DO UM = 1, N_USER_STREAMS
            COSMUM = USER_STREAMS(UM)
            t_delt_userm_num = t_delt_userm(n,um)
            tsgm_uu1n = tcp1n + cosmum * tsgm_uu2n
            tsgm_ud1n = tcm1n + cosmum * tsgm_ud2n
            tsgm_uu0n = - tsgm_uu1n - tsgm_uu2n * deltausn
            tsgm_ud0n = - tsgm_ud1n - tsgm_ud2n * deltausn
            su = tcp0n  * hmult_1(um,n) + &
                      tsgm_uu0n * t_delt_userm_num + tsgm_uu1n
            sd = tcm0n * hmult_2(um,n) + &
                      tsgm_ud0n * t_delt_userm_num + tsgm_ud1n
            layer_tsup_up(um,n) = tterm_saven * ( u_xpos(um,n)*sd + u_xneg(um,n)*su )
          ENDDO
        ELSE
          layer_tsup_up(:,n) = zero
        ENDIF
      ENDDO
    ENDIF

    IF ( DO_DNWELLING ) THEN
      DO N = 1, NLAYERS
        deltausn = deltaus(n)
        IF ( deltaus(n) .gt. TCUTOFF ) THEN
          tsgm_du2n = t_c_plus (2,n)
          tsgm_dd2n = t_c_minus(2,n)
          tcp1n = t_c_plus (1,n)
          tcm1n = t_c_minus(1,n)
          tcp0n = t_c_plus (0,n)
          tcm0n = t_c_minus(0,n)
          tterm_saven = tterm_save(n)
          DO UM = 1, N_USER_STREAMS
            COSMUM = USER_STREAMS(UM)
            t_delt_userm_num = t_delt_userm(n,um)
            tsgm_du1n = tcp1n - cosmum * tsgm_du2n
            tsgm_dd1n = tcm1n - cosmum * tsgm_dd2n
            tsgm_du0n = - tsgm_du1n
            tsgm_dd0n = - tsgm_dd1n
            sum_p = tsgm_du1n  + tsgm_du2n * deltausn
            sum_m = tsgm_dd1n  + tsgm_dd2n * deltausn
            su = tcp0n  * hmult_2(um,n) + &
                    tsgm_du0n * t_delt_userm_num + sum_p
            sd = tcm0n * hmult_1(um,n) + &
                    tsgm_dd0n * t_delt_userm_num + sum_m
            layer_tsup_dn(um,n) = tterm_saven * (  u_xneg(um,n)*sd + u_xpos(um,n)*su )
          ENDDO
        ELSE
          layer_tsup_dn(:,n) = zero
        ENDIF
      ENDDO
    ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSTERMS_THERMAL

SUBROUTINE TWOSTREAM_THERMALSTERMS_SOLARTHERMAL &
      ( MAXLAYERS, MAX_USER_STREAMS,                & ! Dimensions
        DO_UPWELLING, DO_DNWELLING,                 & ! Inputs
        NLAYERS, N_USER_STREAMS, PI4, USER_STREAMS, & ! Inputs
        TCUTOFF, T_DELT_USERM, DELTAUS,             & ! Inputs
        U_XPOS, U_XNEG, HMULT_1, HMULT_2,           & ! Inputs
        T_C_PLUS, T_C_MINUS, TTERM_SAVE,            & ! Inputs
        LAYER_TSUP_UP, LAYER_TSUP_DN  )               ! Outputs

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS, MAX_USER_STREAMS

!  Flags

      LOGICAL, INTENT (IN)        :: DO_UPWELLING
      LOGICAL, INTENT (IN)        :: DO_DNWELLING

!  Constant

      REAL(kind=dp), INTENT (IN)  :: PI4

!  Numbers

      INTEGER, INTENT (IN)        :: NLAYERS
      INTEGER, INTENT (IN)        :: N_USER_STREAMS

!  User streams and layer transmittances

      REAL(kind=dp), INTENT (IN)  :: USER_STREAMS  ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT (IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Thermal Cutoff (actually a layer optical thickness minimum)
!  Rob, introduced 14 May 2015, following 2p3 implementation (2014)
!  Solutions are avoided for optically thin layers

      REAL(kind=dp), INTENT (IN) :: TCUTOFF

!  Optical thickness

      REAL(kind=dp), INTENT (IN)  :: DELTAUS ( MAXLAYERS )

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), intent(in)   :: U_XPOS (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)   :: U_XNEG (MAX_USER_STREAMS,MAXLAYERS)

!  Integrated homogeneous solution multipliers (global, whole layer)

      REAL(kind=dp), intent(in)  :: HMULT_1 (MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: HMULT_2 (MAX_USER_STREAMS,MAXLAYERS)

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: TTERM_SAVE (MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_MINUS (0:2,MAXLAYERS)
      REAL(kind=dp), intent(in)  :: T_C_PLUS  (0:2,MAXLAYERS)

!  Output variables
!  ----------------

      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_UP ( MAX_USER_STREAMS, MAXLAYERS )
      REAL(kind=dp), INTENT (OUT) :: LAYER_TSUP_DN ( MAX_USER_STREAMS, MAXLAYERS )

!  Local variables
! ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: SUM_M, SUM_P, SU, SD, SPAR, COSMUM, DELTAUSN
      REAL(kind=dp) :: TSGM_DU0N, TSGM_DU1N, TSGM_DU2N
      REAL(kind=dp) :: TSGM_DD0N, TSGM_DD1N, TSGM_DD2N
      REAL(kind=dp) :: TSGM_UU0N, TSGM_UU1N, TSGM_UU2N
      REAL(kind=dp) :: TSGM_UD0N, TSGM_UD1N, TSGM_UD2N
      REAL(kind=dp) :: TCP0N, TCP1N, TCM0N, TCM1N
      REAL(kind=dp) :: TTERM_SAVEN, T_DELT_USERM_NUM

!  PARTICULAR SOLUTION LAYER SOURCE TERMS ( GREEN'S FUNCTION SOLUTION )
!  --------------------------------------------------------------------

!  INITIAL MODULUS = 4.PI SINCE SOLAR SOURCES ARE INCLUDED

!  UPWELLING and DOWNWELLING WHOLE LAYER SOURCE TERMS
!  Rob  Fix, 5.14.15. Only get solutions for valid non-Cutoff layers.

!  Upwelling

    IF ( DO_UPWELLING ) THEN
      DO N = 1, NLAYERS
        deltausn = deltaus(n)
        IF ( deltausn .gt. TCUTOFF ) THEN
          tsgm_uu2n = t_c_plus (2,n)
          tsgm_ud2n = t_c_minus(2,n)
          tcp1n = t_c_plus (1,n)
          tcm1n = t_c_minus(1,n)
          tcp0n = t_c_plus (0,n)
          tcm0n = t_c_minus(0,n)
          tterm_saven = tterm_save(n)
          DO UM = 1, N_USER_STREAMS
            COSMUM = USER_STREAMS(UM)
            t_delt_userm_num = t_delt_userm(n,um)
            tsgm_uu1n = tcp1n + cosmum * tsgm_uu2n
            tsgm_ud1n = tcm1n + cosmum * tsgm_ud2n
            tsgm_uu0n = - tsgm_uu1n - tsgm_uu2n * deltausn
            tsgm_ud0n = - tsgm_ud1n - tsgm_ud2n * deltausn
            su = tcp0n  * hmult_1(um,n) + &
                      tsgm_uu0n * t_delt_userm_num + tsgm_uu1n
            sd = tcm0n * hmult_2(um,n) + &
                      tsgm_ud0n * t_delt_userm_num + tsgm_ud1n
            spar = tterm_saven * ( u_xpos(um,n)*sd + u_xneg(um,n)*su )
            layer_tsup_up(um,n) = PI4 * spar
          ENDDO
        ELSE
          layer_tsup_up(:,n) = zero
        ENDIF
      ENDDO
    ENDIF

    IF ( DO_DNWELLING ) THEN
      DO N = 1, NLAYERS
        deltausn = deltaus(n)
        IF ( deltaus(n) .gt. TCUTOFF ) THEN
          tsgm_du2n = t_c_plus (2,n)
          tsgm_dd2n = t_c_minus(2,n)
          tcp1n = t_c_plus (1,n)
          tcm1n = t_c_minus(1,n)
          tcp0n = t_c_plus (0,n)
          tcm0n = t_c_minus(0,n)
          tterm_saven = tterm_save(n)
          DO UM = 1, N_USER_STREAMS
            COSMUM = USER_STREAMS(UM)
            t_delt_userm_num = t_delt_userm(n,um)
            tsgm_du1n = tcp1n - cosmum * tsgm_du2n
            tsgm_dd1n = tcm1n - cosmum * tsgm_dd2n
            tsgm_du0n = - tsgm_du1n
            tsgm_dd0n = - tsgm_dd1n
            sum_p = tsgm_du1n  + tsgm_du2n * deltausn
            sum_m = tsgm_dd1n  + tsgm_dd2n * deltausn
            su = tcp0n  * hmult_2(um,n) + &
                    tsgm_du0n * t_delt_userm_num + sum_p
            sd = tcm0n * hmult_1(um,n) + &
                    tsgm_dd0n * t_delt_userm_num + sum_m
            spar = tterm_saven * (  u_xneg(um,n)*sd + u_xpos(um,n)*su )
            layer_tsup_dn(um,n) = PI4 * spar
          ENDDO
        ELSE
          layer_tsup_dn(:,n) = zero
        ENDIF
      ENDDO
    ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_THERMALSTERMS_SOLARTHERMAL

!  End module

end module twostream_thermalsup_m

