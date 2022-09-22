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
! #          TWOSTREAM_CONVERGE_OBSGEO_SOLAR                #
! #          TWOSTREAM_CONVERGE_OBSGEO_THERMAL              #
! #                                                         #
! ###########################################################

module twostream_converge_obs_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO_SOLAR &
       ( MAXLAYERS,                                           & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,          & ! Inputs     ! @@ 2p2
         NLAYERS, FOURIER_COMPONENT, AZMFACB,                 & ! Inputs     ! @@ 2p2
         INTENSITY_F_UPB,  INTENSITY_F_DNB,                   & ! Inputs
         RADLEVEL_F_UPB,   RADLEVEL_F_DNB,                    & ! Inputs     ! @@ 2p2
         INTENSITY_TOAB,   INTENSITY_BOAB,                    & ! In/Out
         RADLEVEL_UPB,     RADLEVEL_DNB   )                     ! In/Out     ! @@ 2p2

      IMPLICIT NONE

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Fourier component, nlayers

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS

!  Local  azimuth factors

      REAL(kind=dp), INTENT(IN)  :: AZMFACB

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UPB
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DNB

!  Fourier-component solutions at ALL levels

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UPB (0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DNB (0:MAXLAYERS)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOAB
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOAB

!  Output solutions at ALL levels

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UPB (0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DNB (0:MAXLAYERS)

!  Local variables
!  ---------------

      REAL(kind=dp) :: TOLD, TAZM

!  Fourier 0 component

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

        IF ( DO_UPWELLING ) THEN
          INTENSITY_TOAB = INTENSITY_F_UPB
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UPB(0:NLAYERS) = RADLEVEL_F_UPB(0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          INTENSITY_BOAB = INTENSITY_F_DNB
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DNB(0:NLAYERS) = RADLEVEL_F_DNB(0:NLAYERS)
        ENDIF

!  Fourier component = 1

      ELSE
        IF ( DO_UPWELLING ) THEN
          TOLD = INTENSITY_TOAB
          TAZM = AZMFACB*INTENSITY_F_UPB
          INTENSITY_TOAB = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UPB(0:NLAYERS) = &
                   RADLEVEL_UPB(0:NLAYERS) + AZMFACB * RADLEVEL_F_UPB(0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          TOLD = INTENSITY_BOAB
          TAZM = AZMFACB*INTENSITY_F_DNB
          INTENSITY_BOAB = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DNB(0:NLAYERS) = &
                   RADLEVEL_DNB(0:NLAYERS) + AZMFACB * RADLEVEL_F_DNB(0:NLAYERS)
        ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO_SOLAR

!

SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO_THERMAL &
       ( MAX_USER_STREAMS, MAXLAYERS,                    & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
         NLAYERS, N_USER_STREAMS,                        & ! Inputs
         INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
         INTENSITY_TOA, INTENSITY_BOA,                   & ! In/Out
         RADLEVEL_UP,   RADLEVEL_DN   )                    ! In/Out     ! @@ 2p2

      IMPLICIT NONE

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS
      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Control

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Number

      INTEGER, INTENT(IN)        :: N_USER_STREAMS

!  nlayers

      INTEGER, INTENT(IN)        :: NLAYERS 

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS)

!  Fourier-component solutions at ALL levels

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,0:MAXLAYERS)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_USER_STREAMS)

!  Output solutions at ALL levels

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (0:MAXLAYERS,MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (0:MAXLAYERS,MAX_USER_STREAMS)

!  Local variable
!  --------------

      INTEGER       :: I

      DO I = 1, N_USER_STREAMS
        IF ( DO_UPWELLING ) THEN
          INTENSITY_TOA(I) = INTENSITY_F_UP(I)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(0:NLAYERS,I) = RADLEVEL_F_UP(0:NLAYERS,I)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          INTENSITY_BOA(I) = INTENSITY_F_DN(I)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(0:NLAYERS,I) = RADLEVEL_F_DN(0:NLAYERS,I)
        ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO_THERMAL

end module twostream_converge_obs_m
