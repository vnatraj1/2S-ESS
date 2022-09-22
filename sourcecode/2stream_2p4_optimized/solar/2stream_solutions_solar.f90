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
! #     Mark 9  : June   2014, Inverse Pentadiagonal        #
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
! #     Homogeneous solution                                    #
! #                                                             #
! #              TWOSTREAM_HOM_SOLUTION_SOLAR                   #
! #              TWOSTREAM_HOM_USERSOLUTION_SOLAR               #
! #                                                             #
! #     Particular integrals                                    #
! #                                                             #
! #              TWOSTREAM_GBEAM_SOLUTION                       #
! #                                                             #
! ###############################################################

module twostream_solutions_solar_m

!    Introduced for V2p4, Mark 10

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_HOM_SOLUTION_SOLAR &
          ( MAXLAYERS, NLAYERS, FOURIER, STREAM_VALUE, PXSQ, & ! Inputs
            OMEGA, ASYMM, DELTAU_VERT,                       & ! Inputs
            EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED )         ! Out

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)         :: MAXLAYERS

!  Number of layers

      INTEGER, INTENT(IN)         :: NLAYERS

!  Fourier number (input)

      INTEGER, INTENT(IN)         :: FOURIER

!  Stream value

      REAL(kind=dp), INTENT(IN)   :: STREAM_VALUE

!  Polynomials

      REAL(kind=dp), INTENT(IN)   :: PXSQ
      
!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: ASYMM(MAXLAYERS)

!  optical thickness

      REAL(kind=dp), INTENT(IN)   :: DELTAU_VERT(MAXLAYERS)

!  Solutions to the homogeneous RT equations 
!  -----------------------------------------

!  Eigensolutions

      REAL(kind=dp), INTENT(OUT)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT)  :: EIGENTRANS(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(OUT)  :: XPOS(2,MAXLAYERS)

!  Green's function normalization factors
!  Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(OUT)  :: NORM_SAVED(MAXLAYERS)

!  Local variables
!  ---------------

!  local variable for eigenvalue computation

      REAL(kind=dp) :: SAB, DAB

!  parameter Introduced for [V2p3, Mark 10]

      INTEGER       :: N
      REAL(kind=dp) :: EP, EM, XINV, OMEGAN, OMEGA_ASYMM_3, DIFVEC, HELP
      REAL(kind=dp) :: EIGENVALUEN, XPOS1N, XPOS2N
      REAL(kind=dp), parameter :: MAX_TAU_QPATH = 88.0_dp

!  Develop Sum and Difference matrices, set Eigenvalue

      XINV = one / STREAM_VALUE

      DO N = 1, NLAYERS

         OMEGAN = OMEGA(N)
         OMEGA_ASYMM_3 = 3.0_dp * OMEGAN * ASYMM(N)
         if ( fourier .eq. 0 ) then
           EP = OMEGAN + PXSQ * OMEGA_ASYMM_3
           EM = OMEGAN - PXSQ * OMEGA_ASYMM_3
         Else if ( fourier .eq. 1 ) then
           EP = OMEGA_ASYMM_3 * PXSQ
           EM = OMEGA_ASYMM_3 * PXSQ
         ENDIF
         SAB = XINV * ( ( EP + EM ) * 0.5_dp - one )
         DAB = XINV * ( ( EP - EM ) * 0.5_dp - one )
         EIGENVALUEN = SQRT(SAB*DAB)
         EIGENVALUE(N) = EIGENVALUEN

!  Eigentrans, defined properly. [V2p3, Mark 10]

         HELP = EIGENVALUEN*DELTAU_VERT(N)
         IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            EIGENTRANS(N) = zero
         ELSE
            EIGENTRANS(N) = EXP(-HELP)
         ENDIF

!  Auxiliary equation to get up and down solutions

         DIFVEC = - SAB / EIGENVALUEN
         XPOS1N = 0.5d0 * ( one + DIFVEC )
         XPOS2N = 0.5d0 * ( one - DIFVEC )
         XPOS(1,N) = XPOS1N
         XPOS(2,N) = XPOS2N

!  Green's function norm
!  Introduced for Version 2p4, Mark 10

         NORM_SAVED(N) = STREAM_VALUE * ( XPOS1N*XPOS1N - XPOS2N*XPOS2N )

!  debug
!      if (fourier.eq.0)write(*,'(i4,1p2e24.12)')n,EIGENTRANS(N),norm_saved(n)
!      if ( n.eq.23)pause

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_SOLUTION_SOLAR

!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_SOLAR &
         ( MAXLAYERS, MAX_USER_STREAMS,                 & ! Dimensions
           NLAYERS, N_USER_STREAMS, FOURIER, STREAM_VALUE, PX11, & ! Input
           USER_STREAMS, ULP, XPOS, OMEGA, ASYMM,       & ! Input
           U_XPOS, U_XNEG )                               ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  Fourier number (inputs)

      INTEGER, INTENT(IN)        :: FOURIER

!  Stream value and polynomial

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE, PX11

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: ULP          ( MAX_USER_STREAMS )

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)  :: OMEGA(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: ASYMM(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Subroutine output arguments
!  ---------------------------

!  Eigenvectors defined at user-defined stream angles

      REAL(kind=dp), INTENT(OUT) ::  U_XPOS(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) ::  U_XNEG(MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: N, UM
      REAL(kind=dp) :: SUM_NEG, SUM_POS
      REAL(kind=dp) :: OMEGAN, OMEGA_MOM, HMU_STREAM
      REAL(kind=dp) :: XPOS1N, XPOS2N
      REAL(kind=dp) :: U_HELP_P0, U_HELP_P1
      REAL(kind=dp) :: U_HELP_M0, U_HELP_M1
      REAL(kind=dp) :: user_streams_um, ulp_um

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE

      DO N = 1, NLAYERS

         XPOS1N = XPOS(1,N)
         XPOS2N = XPOS(2,N)
         if ( fourier .eq. 0) then
           u_help_p0 = ( XPOS2N + XPOS1N ) * 0.5d0
           u_help_p1 = ( XPOS2N - XPOS1N ) * HMU_STREAM
           u_help_M0 =   u_help_p0
           u_help_M1 = - u_help_p1
         else
           u_help_p1 = - ( XPOS2N + XPOS1N ) * PX11 * 0.5d0
           u_help_M1 = u_help_p1
         endif

!  Now sum over harmonic contributions at each user-defined stream

         OMEGAN = OMEGA(N)
         OMEGA_MOM = 3.0d0 * OMEGAN * ASYMM(N)
         DO UM = 1, N_USER_STREAMS
           if (fourier.eq.0 ) then
             user_streams_um = user_streams(um)
             sum_pos = u_help_p0 * omegaN &
                    +  u_help_p1 * omega_mom * user_streams_um
             sum_neg = u_help_m0 * omegaN &
                    +  u_help_m1 * omega_mom * user_streams_um
           else
             ulp_um = ulp(um)
             sum_pos = u_help_p1 * omega_mom * ulp_um
             sum_neg = u_help_m1 * omega_mom * ulp_um
           endif
           U_XPOS(UM,N) = SUM_POS
           U_XNEG(UM,N) = SUM_NEG
         ENDDO

!  debug
!      if (fourier.eq.1)
!     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,n),u_xneg(1,n)

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_SOLAR

!

SUBROUTINE TWOSTREAM_GBEAM_SOLUTION &
         ( MAXLAYERS, NLAYERS,                                 & ! Inputs
           DO_PLANE_PARALLEL, DO_POSTPROCESSING,               & ! Inputs
           TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS,                & ! Inputs 
           FOURIER, PI4, FLUX_FACTOR,                          & ! Inputs
           LAYER_PIS_CUTOFFB, PX0XB, OMEGA, ASYMM,             & ! Inputs
           AVERAGE_SECANT_PPB, AVERAGE_SECANTB, INITIAL_TRANSB, T_DELT_MUBARB, & ! Inputs
           XPOS, EIGENVALUE, EIGENTRANS, NORM_SAVED,           & ! Inputs
           GAMMA_M, GAMMA_P, ATERM_SAVE, BTERM_SAVE,           & ! Output
           WUPPER, WLOWER )                                      ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  subroutine arguments
!  --------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS

!  Order of Taylor series (including terms up to EPS^n).

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)

!  Flags

      LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL
      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Fourier number

      INTEGER, INTENT(IN)         :: FOURIER

!  Flux factor

      REAL(kind=dp), INTENT(IN)   :: FLUX_FACTOR, PI4

!  Last layer to include Particular integral solution

      INTEGER, INTENT(IN)         :: LAYER_PIS_CUTOFFB

!  Beam SZA polynomial factors

      REAL(kind=dp), INTENT(IN)   :: PX0XB

!  OMEGA and ASYMM

      REAL(kind=dp), INTENT(IN)   :: OMEGA(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: ASYMM(MAXLAYERS)

!  Average-secant and initial transmittance factors for solar beams.

      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANT_PPB
      REAL(kind=dp), INTENT(IN)   :: AVERAGE_SECANTB(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: INITIAL_TRANSB(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: T_DELT_MUBARB(MAXLAYERS)

!  UP and down solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  Eigenvalues and eigentransmittance

      REAL(kind=dp), INTENT(IN)   :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)   :: EIGENTRANS(MAXLAYERS)

!  Green's function normalization factors
!  Introduced for [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: NORM_SAVED(MAXLAYERS)

!  subroutine output arguments
!  ===========================

!  Saved quantities for the Green function solution

      REAL(kind=dp), intent(out) :: ATERM_SAVE(MAXLAYERS)
      REAL(kind=dp), intent(out) :: BTERM_SAVE(MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      REAL(kind=dp), intent(out) :: GAMMA_M(MAXLAYERS)
      REAL(kind=dp), intent(out) :: GAMMA_P(MAXLAYERS)

!  Solutions at layer boundaries

      REAL(kind=dp), INTENT(OUT)  :: WUPPER(2,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT)  :: WLOWER(2,MAXLAYERS)

!  Help variables
!  --------------

!  Layer C and D functions

      REAL(kind=dp) :: CFUNC
      REAL(kind=dp) :: DFUNC

!  Green function Multipliers for solution

      REAL(kind=dp) :: GFUNC_UPN
      REAL(kind=dp) :: GFUNC_DNN

      INTEGER       :: I, N
      REAL(kind=dp) :: TP, TM, EPS, CONST, SECBAR, NORM_SAVEDN
      REAL(kind=dp) :: GAMMA_PN, GAMMA_MN, EIGENVALUEN, OMEGAN, XPOS1N, XPOS2N
      REAL(kind=dp) :: WDEL, ZDEL, ZWDEL, F1, SUM_LA, SUM_LB, OMEGA_ASYMM
      REAL(kind=dp) :: ATERM_SAVEN, BTERM_SAVEN
      REAL(kind=dp) :: DMIN, DPIN

!  Flux factor

      F1 = FLUX_FACTOR / PI4

!  Start loops over layers
!  No particular solution beyond the cutoff layer
!  Or no scattering in this layer...
!  ... Zero the boundary layer values beyond cutoff layer


      DO N = 1, LAYER_PIS_CUTOFFB

!  constants for the layer

         IF (.NOT. DO_PLANE_PARALLEL) THEN
            SECBAR = AVERAGE_SECANTB(N)
         ENDIF
         CONST  = INITIAL_TRANSB(N)
         WDEL   = T_DELT_MUBARB(N)
         EIGENVALUEN = EIGENVALUE(N)

!  Optical depth integrations for the discrete ordinate solution
!  =============================================================

         IF (.NOT. DO_PLANE_PARALLEL) THEN
            GAMMA_PN = SECBAR + EIGENVALUEN
            GAMMA_MN = SECBAR - EIGENVALUEN
         ELSE
            GAMMA_PN = AVERAGE_SECANT_PPB + EIGENVALUEN
            GAMMA_MN = AVERAGE_SECANT_PPB - EIGENVALUEN
         ENDIF
         IF (DO_POSTPROCESSING) THEN
            GAMMA_P(N) = GAMMA_PN
            GAMMA_M(N) = GAMMA_MN
         ENDIF
         ZDEL  = EIGENTRANS(N)
         ZWDEL = ZDEL * WDEL
         IF ( ABS(GAMMA_MN) .LT. TAYLOR_SMALL ) THEN
            EPS = GAMMA_MN
            CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DELTAUS(N), WDEL, ONE, CFUNC )
         ELSE
            CFUNC =  ( ZDEL - WDEL ) / GAMMA_MN
         ENDIF
         DFUNC  = ( one - ZWDEL ) / GAMMA_PN

!  Help quantities for Green's function

         OMEGAN = OMEGA(N)
         OMEGA_ASYMM = OMEGAN * ASYMM(N) * 3.0_dp
         if ( fourier.eq.0) then
            TP = OMEGAN + PX0XB * OMEGA_ASYMM
            TM = OMEGAN - PX0XB * OMEGA_ASYMM
         else if ( fourier .eq. 1 ) then
            TP = PX0XB * OMEGA_ASYMM
            TM = PX0XB * OMEGA_ASYMM
         endif
         DPIN = TP * F1
         DMIN = TM * F1

         XPOS1N = XPOS(1,N) 
         XPOS2N = XPOS(2,N)   
         NORM_SAVEDN = NORM_SAVED(N)

!  Green function multipliers GFUNC

         SUM_LA  = DPIN*XPOS1N+DMIN*XPOS2N
         SUM_LB  = DMIN*XPOS1N+DPIN*XPOS2N
         ATERM_SAVEN = SUM_LA / NORM_SAVEDN
         BTERM_SAVEN = SUM_LB / NORM_SAVEDN
         IF (DO_POSTPROCESSING) THEN
            ATERM_SAVE(N) = ATERM_SAVEN
            BTERM_SAVE(N) = BTERM_SAVEN
         ENDIF
         GFUNC_DNN = CFUNC * ATERM_SAVEN * CONST
         GFUNC_UPN = DFUNC * BTERM_SAVEN * CONST

!  particular integrals at lower and upper boundaries

         WUPPER(1,N) = GFUNC_UPN*XPOS2N
         WUPPER(2,N) = GFUNC_UPN*XPOS1N
         WLOWER(1,N) = GFUNC_DNN*XPOS1N
         WLOWER(2,N) = GFUNC_DNN*XPOS2N

!      if ( FOURIER.eq.1)write(*,*)Fourier, n, WUPPER(1:2,N),WLOWER(1:2,N)

      ENDDO

      DO N = LAYER_PIS_CUTOFFB+1, NLAYERS
         DO I = 1, 2
            WUPPER(I,N) = zero
            WLOWER(I,N) = zero
         ENDDO
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_GBEAM_SOLUTION

end module twostream_solutions_solar_m
