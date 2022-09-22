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

! ###########################################################
! #                                                         #
! #   Contains the following Master subroutines             #
! #                                                         #
! #     TWOSTREAM_UPUSER_INTENSITY_SOLAR_OBS   (master)     #
! #     TWOSTREAM_DNUSER_INTENSITY_SOLAR_OBS   (master)     #
! #                                                         #
! ###########################################################

module twostream_intensity_solar_obs_m

  use Twostream_Taylor_m, only : Twostream_Taylor_Series_2

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_INTENSITY_SOLAR_OBS &
  ( MAXLAYERS,                                                   & ! Dimension
    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_2S_LEVELOUT,         & ! inputs !@@ 2p2
    NLAYERS, TAYLOR_ORDER,                                       & ! inputs !@@ 2p4 Greens
    LAYER_PIS_CUTOFFB, SURFACE_FACTOR, ALBEDO, UBRDF_FMB,        & ! inputs
    FLUXMULT, STREAM_VALUE, TAYLOR_SMALL, DELTAU_VERT,           & ! inputs
    GAMMA_P, GAMMA_M, SIGMA_PB, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p4 Greens
    INITIAL_TRANSB, ITRANS_USERMB, T_DELT_USERMB, T_DELT_MUBARB, & ! Inputs !@@ 2p4 Greens
    T_DELT_EIGENNL, LCON, LCON_XVEC1NL, MCON, MCON_XVEC1NL, WLOWER1NL, & ! inputs
    U_XPOSB, U_XNEGB, HMULT_1B, HMULT_2B, EMULT_UPB,             & ! inputs
    INTENSITY_F_UPB, RADLEVEL_F_UPB )                                ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimension

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  local surface control flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Number

      INTEGER, INTENT(IN)        :: NLAYERS

!  Taylor control. Greens function solution

      INTEGER, INTENT(IN)        :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_FMB

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT

!  Transmittance factors for user-defined stream angle

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERMB ( MAXLAYERS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)      :: LAYER_PIS_CUTOFFB

!  Quantities for Greens function solution
!  ---------------------------------------

!  Average secant/eigenvalue coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_P      ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_M      ( MAXLAYERS )

!  Average secant/user secant coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_PB     ( MAXLAYERS )

!  Input optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!  Initial transmittance factors for solar beam, and divided by user-cosine

      REAL(kind=dp), intent(in)  :: INITIAL_TRANSB ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: ITRANS_USERMB  ( MAXLAYERS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), intent(in)  :: T_DELT_MUBARB  ( MAXLAYERS )

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Older variables (Version 2p2)
!  -----------------------------

!  No USER_DIRECT_BEAM (MSMODE only ===> No Direct BOA source term)
!      DOUBLE PRECISION USER_DIRECT_BEAMB

!  Transmittance factors for +/- eigenvalues for lowest layer

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGENNL

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions for downwelling at BOA

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC1NL
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC1NL

!  Beam particular downwelling solution at BOA

      REAL(kind=dp), INTENT(IN)  :: WLOWER1NL

!  Eigenvectors defined at user-defined stream angle
!  EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOSB(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEGB(MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angle
!  @@@ NOT REQUIRED For MS-mode only
!     REAL(kind=dp), INTENT(IN)  :: U_WPOS1(MAXLAYERS)

!  Solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1B (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2B (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_UPB(MAXLAYERS)

!  Outputs
!  -------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UPB

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UPB (0:MAXLAYERS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER       :: N, N1
      REAL(kind=dp) :: LCONN, MCONN
      REAL(kind=dp) :: LAYERSOURCE
      REAL(kind=dp) :: BOA_SOURCE, IDOWNSURF
      REAL(kind=dp) :: SHOM, SPAR, PAR, HOM, SU, SD
      REAL(kind=dp) :: FAC2, ITRANS, WDEL, MULT
      REAL(kind=dp) :: GAMMAMN, EMULT_UPN
      REAL(kind=dp) :: U_XPOSN, U_XNEGN, HMULT_1N, HMULT_2N
      REAL(kind=dp) :: PMULT_UU, PMULT_UD
      REAL(kind=dp) :: CUMSOURCE_UP_OLD, CUMSOURCE_UP_NEW

!  BOA source terms
!  ----------------
!  MSMODE only ===> No Direct BOA source term

      IF ( DO_INCLUDE_SURFACE )THEN

!  Full solution: Downward intensity at computational angles (beam/homog)
!  --> Develop reflectance integrand  a(j).x(j).I(-j)

        PAR = WLOWER1NL
        HOM = LCON_XVEC1NL*T_DELT_EIGENNL + MCON_XVEC1NL
        IDOWNSURF = ( PAR + HOM ) * STREAM_VALUE

!  Reflected multiple scatter intensity at user defined-angles

        IF ( DO_BRDF_SURFACE  ) THEN
          BOA_SOURCE = SURFACE_FACTOR * IDOWNSURF * UBRDF_FMB
        ELSE
          BOA_SOURCE = SURFACE_FACTOR * ALBEDO * IDOWNSURF
        ENDIF

!  MSMODE only ===> No Direct BOA source term
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!            DIRECT_BOA_SOURCE = USER_DIRECT_BEAMB
!        ENDIF

      ELSE

        BOA_SOURCE = ZERO

      ENDIF

!  Add surface emission term if flagged
!  ********  Direct surface emission not included

!  Initialize post-processing recursion
!  ====================================

!  Set the cumulative source term equal to BOA values
!  MSMODE only ===> No Direct BOA source term

      CUMSOURCE_UP_OLD = BOA_SOURCE
      IF ( DO_2S_LEVELOUT ) RADLEVEL_F_UPB(NLAYERS) = FLUXMULT * CUMSOURCE_UP_OLD

!  Recursion Loop in Source function integration
!  =============================================

      DO N = NLAYERS, 1, -1
        N1 = N - 1

        LCONN = LCON(N)
        MCONN = MCON(N)

!write(0,*) 'LM',LCONN, MCONN

!  Only present if not beyond the cut-off layer

        IF ( N .LE. LAYER_PIS_CUTOFFB ) THEN

!  Layer quantities

          WDEL       = T_DELT_MUBARB(N)
          ITRANS     = INITIAL_TRANSB(N)  
          GAMMAMN    = GAMMA_M(N)
          EMULT_UPN = EMULT_UPB(N)

        ENDIF

        U_XPOSN =  U_XPOSB(N)
        U_XNEGN =  U_XNEGB(N)
        HMULT_1N = HMULT_1B(N)
        HMULT_2N = HMULT_2B(N)

!  Homogeneous solutions.

        SHOM = LCON(N) * U_XPOSN * HMULT_2N + &
               MCON(N) * U_XNEGN * HMULT_1N
        LAYERSOURCE = SHOM

!write(0,*) 'SHOM',N,SHOM

!  Add solar source term
!  Only present if not beyond the cut-off layer

        IF ( N .LE. LAYER_PIS_CUTOFFB ) THEN

!  Multipliers PMULT_UD, PMULT_UU; SPAR = Greens function addition

          IF ( ABS(GAMMAMN) .LT. TAYLOR_SMALL ) THEN
             FAC2 = WDEL * T_DELT_USERMB(N)
             CALL Twostream_Taylor_Series_2 ( TAYLOR_ORDER, TAYLOR_SMALL, GAMMAMN, SIGMA_PB(N), &
                                              DELTAU_VERT(N), ONE, FAC2, ONE, MULT )
             SD = ITRANS_USERMB(N) * MULT
          ELSE
             SD = ( ITRANS * HMULT_2N - EMULT_UPN ) / GAMMAMN
          ENDIF
          SU = ( - ITRANS * WDEL * HMULT_1N + EMULT_UPN ) / GAMMA_P(N)
          PMULT_UD = SD * ATERM_SAVE(N)
          PMULT_UU = SU * BTERM_SAVE(N)
          SPAR = U_XPOSN*PMULT_UD + U_XNEGN*PMULT_UU
          LAYERSOURCE = LAYERSOURCE + SPAR

        ENDIF

!write(0,*) 'SHOM+SPAR',N,LAYERSOURCE

!  Upward recursion for source function.

!write(0,*) N,LAYERSOURCE,T_DELT_USERMB(N),CUMSOURCE_UP_OLD
        CUMSOURCE_UP_NEW = LAYERSOURCE + T_DELT_USERMB(N)*CUMSOURCE_UP_OLD
        IF (DO_2S_LEVELOUT)RADLEVEL_F_UPB(N1) = FLUXMULT * CUMSOURCE_UP_NEW
        CUMSOURCE_UP_OLD = CUMSOURCE_UP_NEW

!  End recursion loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term at TOA

      INTENSITY_F_UPB = FLUXMULT * CUMSOURCE_UP_OLD

!  Debug
!      write(*,*)INTENSITY_F_UPB

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_INTENSITY_SOLAR_OBS

!

SUBROUTINE TWOSTREAM_DNUSER_INTENSITY_SOLAR_OBS &
   ( MAXLAYERS,                                                   & ! Dimensions
     DO_2S_LEVELOUT,                                              & ! Inputs !@@ 2p1, 2p2
     NLAYERS, TAYLOR_ORDER,                                       & ! Inputs !@@ 2p3 Greens
     LAYER_PIS_CUTOFFB, FLUXMULT, TAYLOR_SMALL, DELTAU_VERT,      & ! Inputs
     GAMMA_P, GAMMA_M, SIGMA_MB, ATERM_SAVE, BTERM_SAVE,          & ! Inputs !@@ 2p3 Greens
     INITIAL_TRANSB, ITRANS_USERMB, T_DELT_USERMB, T_DELT_MUBARB, & ! Inputs !@@ 2p3 Greens
     LCON, MCON, U_XPOSB, U_XNEGB,                                & ! Inputs
     HMULT_1B, HMULT_2B, EMULT_DNB,                               & ! Inputs
     INTENSITY_F_DNB, RADLEVEL_F_DNB )                              ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimension

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Number

      INTEGER, INTENT(IN)        :: NLAYERS

!  Taylor control. Greens function solution

      INTEGER, INTENT(IN)        :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Multiplier, 4pi

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT

!  Transmittance factors for user-defined stream angle

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERMB ( MAXLAYERS )

!  Last layer to include Particular integral solution

      INTEGER  , intent(in)      :: LAYER_PIS_CUTOFFB

!  New quantities for Greens function solution
!  -------------------------------------------

!  Average secant/eigenvalue coefficients

      REAL(kind=dp), intent(in)  :: GAMMA_P      ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: GAMMA_M      ( MAXLAYERS )

!  Average secant/user secant coefficients

      REAL(kind=dp), intent(in)  :: SIGMA_MB     ( MAXLAYERS )

!  Input optical depths

      REAL(kind=dp), intent(in)  :: DELTAU_VERT  ( MAXLAYERS )

!  Initial transmittance factors for solar beam, and divided by user-cosine

      REAL(kind=dp), intent(in)  :: INITIAL_TRANSB ( MAXLAYERS )
      REAL(kind=dp), intent(in)  :: ITRANS_USERMB  ( MAXLAYERS )

!  Transmittance factors for average secant stream

      REAL(kind=dp), intent(in)  :: T_DELT_MUBARB  ( MAXLAYERS )

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(in)  :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(in)  :: BTERM_SAVE(MAXLAYERS)

!  Older variables
!  ---------------

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!  EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOSB(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEGB(MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!  @@@ NOT REQUIRED For MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1B(MAXLAYERS)

!  Solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1B (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2B (MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EMULT_DNB(MAXLAYERS)

!  Outputs
!  -------

!  User-defined solution

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DNB

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DNB (0:MAXLAYERS)

!  Local variables
!  ---------------

!  Help variables

      INTEGER       :: N, NC
      REAL(kind=dp) :: LCONN, MCONN
      REAL(kind=dp) :: LAYERSOURCE
      REAL(kind=dp) :: SHOM, SPAR, SU, SD
      REAL(kind=dp) :: ITRANS, WDEL, MULT
      REAL(kind=dp) :: GAMMAMN, EMULT_DNN
      REAL(kind=dp) :: U_XPOSN, U_XNEGN, HMULT_1N, HMULT_2N
      REAL(kind=dp) :: PMULT_DU, PMULT_DD
      REAL(kind=dp) :: CUMSOURCE_DN_OLD, CUMSOURCE_DN_NEW

!  Recursion Loop in Source function integration
!  =============================================

!  Initialize recursion for user-defined stream angles only

      CUMSOURCE_DN_OLD = ZERO

      DO N = 1, NLAYERS

        NC = N

        LCONN = LCON(N)
        MCONN = MCON(N)

!  Only present if not beyond the cut-off layer

        IF ( N .LE. LAYER_PIS_CUTOFFB ) THEN

!  Layer quantities

          WDEL       = T_DELT_MUBARB(N)
          ITRANS     = INITIAL_TRANSB(N)
          GAMMAMN    = GAMMA_M(N)
          EMULT_DNN  = EMULT_DNB(N)

        ENDIF

        U_XPOSN =  U_XPOSB(N)
        U_XNEGN =  U_XNEGB(N)
        HMULT_1N = HMULT_1B(N)
        HMULT_2N = HMULT_2B(N)

!  Homogeneous solutions.

        SHOM = LCON(N) * U_XNEGN * HMULT_1N + &
               MCON(N) * U_XPOSN * HMULT_2N
        LAYERSOURCE = SHOM

!  Add solar source term
!  Only present if not beyond the cut-off layer

        IF ( N .LE. LAYER_PIS_CUTOFFB ) THEN

!  Multipliers PMULT_UD, PMULT_UU; SPAR = Greens function addition

          IF ( ABS(GAMMAMN) .lt. TAYLOR_SMALL ) THEN
             CALL Twostream_Taylor_Series_2 ( TAYLOR_ORDER, TAYLOR_SMALL, GAMMAMN, SIGMA_MB(N), &
                                                DELTAU_VERT(N), T_DELT_USERMB(N), WDEL, ONE, MULT )
             SD = ITRANS_USERMB(N) * MULT 
          ELSE
             SD = ( ITRANS * HMULT_1N - EMULT_DNN ) / GAMMAMN
          ENDIF
          SU = ( - ITRANS * WDEL * HMULT_2N + EMULT_DNN ) / GAMMA_P(N)
          PMULT_DD = SD * ATERM_SAVE(N)
          PMULT_DU = SU * BTERM_SAVE(N)
          SPAR = U_XNEGN*PMULT_DD + U_XPOSN*PMULT_DU
          LAYERSOURCE = LAYERSOURCE + SPAR

        ENDIF

!  Upward recursion for source function.

        CUMSOURCE_DN_NEW = LAYERSOURCE + T_DELT_USERMB(N)*CUMSOURCE_DN_OLD
        IF (DO_2S_LEVELOUT)RADLEVEL_F_DNB(N) = FLUXMULT * CUMSOURCE_DN_NEW
        CUMSOURCE_DN_OLD = CUMSOURCE_DN_NEW

!  End recursion loop

      ENDDO

!  User-defined stream output, just set to the cumulative source term at BOA

      INTENSITY_F_DNB = FLUXMULT * CUMSOURCE_DN_OLD

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_INTENSITY_SOLAR_OBS

end module twostream_intensity_solar_obs_m
