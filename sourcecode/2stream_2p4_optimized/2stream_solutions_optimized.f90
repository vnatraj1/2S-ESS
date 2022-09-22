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
! #              TWOSTREAM_HMULT_MASTER                         #
! #                                                             #
! ###############################################################

module twostream_solutions_m

!    Introduced for V2p4, Mark 10

   use Twostream_Taylor_m, only : Twostream_Taylor_Series_1

PUBLIC

contains

SUBROUTINE TWOSTREAM_HMULT_MASTER &
           ( MAXLAYERS, MAX_USER_STREAMS,            & ! Dimensions
             TAYLOR_ORDER, TAYLOR_SMALL, DELTAUS,    & ! Inputs 
             NLAYERS, N_USER_STREAMS, USER_SECANTS,  & ! Inputs
             EIGENVALUE, EIGENTRANS, T_DELT_USERM,   & ! Inputs
             HMULT_1, HMULT_2 )                        ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Input arguments
!  ===============

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_STREAMS

!  Order of Taylor series (including terms up to EPS^n).
!    Introduced for [V2p3, Mark 10]

      INTEGER      , intent(in)  :: TAYLOR_ORDER
      REAL(kind=dp), INTENT(IN)  :: TAYLOR_SMALL

!  Input Optical depths required for Taylor-series limiting cases
!    Introduced for [V2p3, Mark 10]

      REAL(kind=dp), intent(in)  :: DELTAUS(MAXLAYERS)

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_STREAMS

!  User secants (formerly streams). [V2p3, Mark 10]

      REAL(kind=dp), INTENT(IN)  :: USER_SECANTS ( MAX_USER_STREAMS )

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS )

!  Eigensolutions

      REAL(kind=dp), INTENT(IN)  :: EIGENVALUE(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: EIGENTRANS(MAXLAYERS)

!  Output = Global multipliers
!  ===========================

!  Integrated homogeneous solution multipliers, whole layer

      REAL(kind=dp), INTENT(OUT) :: HMULT_1(MAX_USER_STREAMS,MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: HMULT_2(MAX_USER_STREAMS,MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: UM, N
      REAL(kind=dp) :: UDEL, SM, EPS, ZDEL, ZUDEL
      REAL(kind=dp) :: EIGVN, EIGTN, DN, ZP, ZM

!  whole layer multipliers
!  -----------------------

!  Start loops over layers and user-streams

      DO N = 1, NLAYERS
        EIGVN = EIGENVALUE(N)
        EIGTN = EIGENTRANS(N)
        DN = DELTAUS(N)
        DO UM = 1, N_USER_STREAMS
          UDEL = T_DELT_USERM(N,UM)
          SM   = USER_SECANTS(UM)
          ZP = SM + EIGVN
          ZM = SM - EIGVN
          ZDEL    = EIGTN
          ZUDEL   = ZDEL * UDEL
          HMULT_2(UM,N) = SM * ( ONE - ZUDEL ) / ZP
          IF ( ABS(ZM) .LT. TAYLOR_SMALL ) THEN
            EPS = ZM
            CALL TWOSTREAM_TAYLOR_SERIES_1 ( TAYLOR_ORDER, EPS, DN, UDEL, SM, HMULT_1(UM,N) )
          ELSE
            HMULT_1(UM,N) = SM * ( ZDEL - UDEL ) / ZM
          ENDIF
        ENDDO
      ENDDO

!  debug
!      do n = 1, 3
!        write(*,*)HMULT_1(1,N),HMULT_2(2,N)
!      enddo
!      pause

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_HMULT_MASTER

end module twostream_solutions_m
