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
! #          TWOSTREAM_CONVERGE (master)                    #
! #          TWOSTREAM_CONVERGE_OBSGEO (master) !@@ 2p1     #
! #                                                         #
! ###########################################################

module twostream_intensity_obs_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_CONVERGE &
       ( MAX_USER_STREAMS, MAX_USER_RELAZMS,             & ! Dimensions
         MAX_GEOMETRIES, MAXLAYERS,                      & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,     & ! Inputs     ! @@ 2p2
         NLAYERS, FOURIER_COMPONENT,                     & ! Inputs     ! @@ 2p2
         N_USER_STREAMS, N_USER_RELAZMS, AZMFAC, UMOFF,  & ! Inputs
         INTENSITY_F_UP,  INTENSITY_F_DN,                & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                 & ! Inputs     ! @@ 2p2
         INTENSITY_TOA, INTENSITY_BOA,                   & ! In/Out
         RADLEVEL_UP,   RADLEVEL_DN   )                    ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAX_GEOMETRIES, MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_STREAMS, N_USER_RELAZMS

!  Fourier component and  beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS

!  Local  azimuth factors

      INTEGER, INTENT(IN)        :: UMOFF ( MAX_USER_STREAMS )
      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_RELAZMS,MAX_USER_STREAMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA(MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA(MAX_GEOMETRIES)

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (0:MAXLAYERS,MAX_GEOMETRIES)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (0:MAXLAYERS,MAX_GEOMETRIES)

!  Local variables
!  ---------------

      INTEGER       :: I, UA, V, UMOFFI
      REAL(kind=dp) :: TOLD, TAZM

!  ###################
!  Fourier 0 component
!  ###################

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

!  Copy DIFFUSE Fourier component at all output angles and optical depths
!    If no SSCORR and no DBCORR, then two options apply:
!     (a) Convergence on RADIANCE = DIFFUSE + SSTRUNCATED + DBTRUNCATED
!              (full radiance, no SS correction, no DB correction)
!     (b) Convergence on RADIANCE = DIFFUSE alone (MS only mode)
!              (SSTRUNCATED + DBTRUNCATED do not get calculated)

!  Code only for the NON-SSFULL case
!        IF ( .not. DO_SSFULL ) THEN

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

          DO I = 1, N_USER_STREAMS
            UMOFFI = UMOFF(I)
            DO UA = 1, N_USER_RELAZMS
              V = UMOFFI + UA
              IF ( DO_UPWELLING ) THEN
                INTENSITY_TOA(V) = INTENSITY_F_UP(I)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(0:NLAYERS,V) = RADLEVEL_F_UP(I,0:NLAYERS)
              ENDIF
              IF ( DO_DNWELLING ) THEN
                INTENSITY_BOA(V) = INTENSITY_F_DN(I)
                IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(0:NLAYERS,V) = RADLEVEL_F_DN(I,0:NLAYERS)
              ENDIF
            ENDDO
          ENDDO

!  Commented out in the streamlined version - NO SS OUTPUT
!        IF ( DO_SSFULL ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V) = 0.0d0
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V) = 0.0d0
!              ENDIF
!            ENDDO
!          ENDDO
!        ENDIF

!    Add the single scatter component if flagged
!  Commented out in the streamlined version - NO SS OUTPUT HERE
!        IF ( DO_SSFULL.OR.DO_SSCORR_NADIR.OR.DO_SSCORR_OUTGOING ) THEN
!           DO I = 1, N_USER_STREAMS
!            DO UA = 1, N_USER_RELAZMS
!              V = UMOFF(IBEAM,I) + UA
!              IF ( DO_UPWELLING ) THEN
!                INTENSITY_TOA(V) = 
!     &             INTENSITY_TOA(V) + INTENSITY_SS_UP(V)
!              ENDIF
!              IF ( DO_DNWELLING ) THEN
!                INTENSITY_BOA(V) = 
!     &             INTENSITY_BOA(V) + INTENSITY_SS_DN(V)
!              ENDIF
!            ENDDO
!          ENDDO
!      ENDIF

!  ######################
!  Fourier component = 1
!  ######################

      ELSE

!  No examination of convergence
!  -----------------------------

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

        DO I = 1, N_USER_STREAMS
          UMOFFI = UMOFF(I)
          DO UA = 1, N_USER_RELAZMS
            V = UMOFFI + UA
            IF ( DO_UPWELLING ) THEN
              TOLD = INTENSITY_TOA(V)
              TAZM = AZMFAC(UA,I)*INTENSITY_F_UP(I)
              INTENSITY_TOA(V) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(0:NLAYERS,V) = &
                   RADLEVEL_UP(0:NLAYERS,V) + AZMFAC(UA,I) * RADLEVEL_F_UP(I,0:NLAYERS)
            ENDIF
            IF ( DO_DNWELLING ) THEN
              TOLD = INTENSITY_BOA(V)
              TAZM = AZMFAC(UA,I)*INTENSITY_F_DN(I)
              INTENSITY_BOA(V) = TOLD + TAZM
              IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(0:NLAYERS,V) = &
                   RADLEVEL_DN(0:NLAYERS,V) + AZMFAC(UA,I) * RADLEVEL_F_DN(I,0:NLAYERS)
            ENDIF
          ENDDO
        ENDDO

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE

!

SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO &
       ( MAX_USER_STREAMS, MAX_USER_RELAZMS,                & ! Dimensions
         MAXLAYERS,                                         & ! Dimensions ! @@ 2p2
         DO_UPWELLING, DO_DNWELLING, DO_2S_LEVELOUT,        & ! Inputs     ! @@ 2p2
         NLAYERS, FOURIER_COMPONENT, AZMFAC,                & ! Inputs     ! @@ 2p2
         INTENSITY_F_UP,  INTENSITY_F_DN,                   & ! Inputs
         RADLEVEL_F_UP,   RADLEVEL_F_DN,                    & ! Inputs     ! @@ 2p2
         INTENSITY_TOA,   INTENSITY_BOA,                    & ! In/Out
         RADLEVEL_UP,     RADLEVEL_DN   )                     ! In/Out     ! @@ 2p2

!  Alterations for version 2.2, 17 July 2013

      implicit none

!  Precision

      INTEGER, PARAMETER :: dp = KIND( 1.0D0 )

!  Input variables
!  ---------------

!  Dimensions (2p2, add MAXLAYERS)

      INTEGER, INTENT(IN)        :: MAX_USER_STREAMS, MAX_USER_RELAZMS
      INTEGER, INTENT(IN)        :: MAXLAYERS

!  Control
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_UPWELLING, DO_DNWELLING
      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  SS control, not required in this streamlined version
!      LOGICAL, INTENT(IN) :: DO_SSFULL, DO_SSCORR_OUTGOING, DO_SSCORR_NADIR

!  Fourier component and beam, nlayers (2p2, added)

      INTEGER, INTENT(IN)        :: FOURIER_COMPONENT, NLAYERS

!  Local  azimuth factors

      REAL(kind=dp), INTENT(IN)  :: AZMFAC(MAX_USER_RELAZMS,MAX_USER_STREAMS)

!  User-defined solutions

      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_UP(MAX_USER_STREAMS)
      REAL(kind=dp), INTENT(IN)  :: INTENSITY_F_DN(MAX_USER_STREAMS)

!  Fourier-component solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_UP (MAX_USER_STREAMS,0:MAXLAYERS)
      REAL(kind=dp), INTENT(IN) :: RADLEVEL_F_DN (MAX_USER_STREAMS,0:MAXLAYERS)

!  Single scatter solutions, Not required here
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_UP(MAX_GEOMETRIES)
!      REAL(kind=dp), INTENT(IN)  :: INTENSITY_SS_DN(MAX_GEOMETRIES)

!  Output
!  ------

!  TOA and BOA output

      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_TOA
      REAL(kind=dp), INTENT(INOUT) :: INTENSITY_BOA

!  output solutions at ALL levels
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_UP (0:MAXLAYERS)
      REAL(kind=dp), INTENT(INOUT) :: RADLEVEL_DN (0:MAXLAYERS)

!  Local variables
!  ---------------

      INTEGER       :: LUM, LUA
      REAL(kind=dp) :: TOLD, TAZM

!  Local user indices

      LUM = 1
      LUA = 1

!  Fourier 0 component
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      IF ( FOURIER_COMPONENT.EQ.0 ) THEN

        IF ( DO_UPWELLING ) THEN
          INTENSITY_TOA = INTENSITY_F_UP(LUM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(0:NLAYERS) = RADLEVEL_F_UP(LUM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          INTENSITY_BOA = INTENSITY_F_DN(LUM)
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(0:NLAYERS) = RADLEVEL_F_DN(LUM,0:NLAYERS)
        ENDIF

!  Fourier component = 1
!     ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      ELSE
        IF ( DO_UPWELLING ) THEN
          TOLD = INTENSITY_TOA
          TAZM = AZMFAC(LUA,LUM)*INTENSITY_F_UP(LUM)
          INTENSITY_TOA = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_UP(0:NLAYERS) = &
                   RADLEVEL_UP(0:NLAYERS) + AZMFAC(LUA,LUM) * RADLEVEL_F_UP(LUM,0:NLAYERS)
        ENDIF
        IF ( DO_DNWELLING ) THEN
          TOLD = INTENSITY_BOA
          TAZM = AZMFAC(LUA,LUM)*INTENSITY_F_DN(LUM)
          INTENSITY_BOA = TOLD + TAZM
          IF ( DO_2S_LEVELOUT ) RADLEVEL_DN(0:NLAYERS) = &
                   RADLEVEL_DN(0:NLAYERS) + AZMFAC(LUA,LUM) * RADLEVEL_F_DN(LUM,0:NLAYERS)
        ENDIF

!  Finish Fourier

      ENDIF

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_CONVERGE_OBSGEO

end module twostream_intensity_obs_m
