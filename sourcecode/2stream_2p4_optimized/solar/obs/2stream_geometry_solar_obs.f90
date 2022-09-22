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
! #              TWOSTREAM_AUXGEOM_SOLAR_OBS                    #
! #                                                             #
! ###############################################################

module twostream_geometry_solar_obs_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_AUXGEOM_SOLAR_OBS &
      ( MAX_USER_OBSGEOMS,              & ! Dimension
        DO_POSTPROCESSING,              & ! Flag
        N_USER_OBSGEOMS,                & ! Input
        X0, USER_STREAMS, STREAM_VALUE, & ! Inputs
        PX11, PXSQ, PX0X, ULP )           ! Outputs

      implicit none

!  precision

      INTEGER, PARAMETER :: dp     = KIND( 1.0D0 )

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAX_USER_OBSGEOMS

!  Flag for post-processing

      LOGICAL, INTENT(IN)        :: DO_POSTPROCESSING

!  Numbers

      INTEGER, INTENT(IN)        :: N_USER_OBSGEOMS

!  Stream directions

      REAL(kind=dp), INTENT(IN)  :: X0 ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: USER_STREAMS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

      REAL(kind=dp), INTENT(OUT) :: PX11, PXSQ(0:1)
      REAL(kind=dp), INTENT(OUT) :: PX0X ( 0:1, MAX_USER_OBSGEOMS )
      REAL(kind=dp), INTENT(OUT) :: ULP ( MAX_USER_OBSGEOMS )

!  Local variables

      INTEGER       :: UM, IBEAM
      REAL(kind=dp) :: MU, MU0, POX, SVSQ

!  Saved quantities

      IF ( DO_POSTPROCESSING ) THEN
        DO UM = 1, N_USER_OBSGEOMS
          MU = USER_STREAMS(UM)
          ULP(UM) =  -DSQRT(0.5D0*(1.0d0-MU*MU))
        ENDDO
      ELSE
        ULP = 0.0D0
      ENDIF

      SVSQ = STREAM_VALUE * STREAM_VALUE
      PXSQ(0) = SVSQ
      PX11 = DSQRT(0.5D0*(1.0D0-SVSQ))
      PXSQ(1) = PX11 * PX11

      DO IBEAM = 1, N_USER_OBSGEOMS
         MU0 = X0(IBEAM)
         POX = DSQRT(0.5D0*(1.0D0-MU0*MU0))
         PX0X(0,IBEAM) = MU0 * STREAM_VALUE
         PX0X(1,IBEAM) = POX * PX11
      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_AUXGEOM_SOLAR_OBS

end module twostream_geometry_solar_obs_m
