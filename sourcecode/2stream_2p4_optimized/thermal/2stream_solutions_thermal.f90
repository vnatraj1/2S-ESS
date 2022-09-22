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
! #              TWOSTREAM_HOM_SOLUTION_THERMAL                 #
! #              TWOSTREAM_HOM_USERSOLUTION_THERMAL             #
! #                                                             #
! ###############################################################

module twostream_solutions_thermal_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_HOM_SOLUTION_THERMAL &
          ( MAXLAYERS, NLAYERS, STREAM_VALUE, PXSQ,    & ! Inputs
            OMEGA, ASYMM, DELTAU_VERT,                 & ! Inputs
            EIGENVALUE, EIGENTRANS, XPOS, NORM_SAVED )   ! Out

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
         EP = OMEGAN + PXSQ * OMEGA_ASYMM_3
         EM = OMEGAN - PXSQ * OMEGA_ASYMM_3
         SAB = XINV * ( ( EP + EM ) * 0.5_dp - one )
         DAB = XINV * ( ( EP - EM ) * 0.5_dp - one )
         EIGENVALUEN = SQRT(SAB*DAB)
         EIGENVALUE(N) = EIGENVALUEN

!  Eigentrans, defined properly. [V2p3, Mark 10]

         HELP = EIGENVALUEN*DELTAU_VERT(N)
         IF ( HELP .GT. MAX_TAU_QPATH ) THEN
            EIGENTRANS(N) = ZERO
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

         NORM_SAVED(N) = STREAM_VALUE * ( XPOS1N*XPOS1N - XPOS2N*XPOS2N )

!  debug
!      if (fourier.eq.0)write(*,'(i4,1p2e24.12)')n,EIGENTRANS(N),norm_saved(n)
!      if ( n.eq.23)pause

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_SOLUTION_THERMAL

!

SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_THERMAL &
         ( MAXLAYERS, MAX_USER_STREAMS,           & ! Dimensions
           NLAYERS, N_USER_STREAMS, STREAM_VALUE, & ! Input
           USER_STREAMS, XPOS, OMEGA, ASYMM,      & ! Input
           U_XPOS, U_XNEG )                         ! Output

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

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  User-defined post-processing stream directions

      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_STREAMS )

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
      REAL(kind=dp) :: user_streams_um

!  Eigenvector interpolation to user-defined angles
!  ------------------------------------------------

!  For each moment, do inner sum over computational angles
!  for the positive and negative eigenvectors

      HMU_STREAM = 0.5d0 * STREAM_VALUE

      DO N = 1, NLAYERS

         XPOS1N = XPOS(1,N)
         XPOS2N = XPOS(2,N)
         u_help_p0 = ( XPOS2N + XPOS1N ) * 0.5d0
         u_help_p1 = ( XPOS2N - XPOS1N ) * HMU_STREAM
         u_help_M0 =   u_help_p0
         u_help_M1 = - u_help_p1

!  Now sum over harmonic contributions at each user-defined stream

         OMEGAN = OMEGA(N)
         OMEGA_MOM = 3.0d0 * OMEGAN * ASYMM(N)
         DO UM = 1, N_USER_STREAMS
           user_streams_um = user_streams(um)
           sum_pos = u_help_p0 * omegaN &
                  +  u_help_p1 * omega_mom * user_streams_um
           sum_neg = u_help_m0 * omegaN &
                  +  u_help_m1 * omega_mom * user_streams_um
           U_XPOS(UM,N) = SUM_POS
           U_XNEG(UM,N) = SUM_NEG
         ENDDO

!  debug
!      if (fourier.eq.1)
!     &   write(57,'(i4,1p2e24.12)')n,u_xpos(1,n),u_xneg(1,n)

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HOM_USERSOLUTION_THERMAL

end module twostream_solutions_thermal_m
