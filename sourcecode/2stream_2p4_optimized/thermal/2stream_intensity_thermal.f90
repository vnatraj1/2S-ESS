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
! #       TWOSTREAM_UPUSER_INTENSITY_THERMAL (master)       #
! #       TWOSTREAM_DNUSER_INTENSITY_THERMAL (master)       #
! #                                                         #
! ###########################################################

module twostream_intensity_thermal_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_UPUSER_INTENSITY_THERMAL &
  ( MAXLAYERS, MAX_USER_OBSGEOMS,                                & ! Dimensions
    DO_INCLUDE_SURFACE, DO_BRDF_SURFACE, DO_2S_LEVELOUT,         & ! inputs
    NLAYERS, N_USER_OBSGEOMS,                                    & ! inputs !@@ 2p4 Greens
    SURFACE_FACTOR, ALBEDO, UBRDF_FM,                            & ! inputs
    FLUXMULT, STREAM_VALUE,                                      & ! inputs
    T_DELT_USERM,                                                & ! Inputs !@@ 2p4 Greens
    T_DELT_EIGENNL, LCON, LCON_XVEC1NL, MCON, MCON_XVEC1NL, WLOWER1NL, & ! inputs
    U_XPOS, U_XNEG, HMULT_1, HMULT_2, LAYER_TSUP_UP,             & ! inputs
    INTENSITY_F_UP, RADLEVEL_F_UP )                                ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_OBSGEOMS

!  local surface control flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_BRDF_SURFACE

!  ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_OBSGEOMS

!  Surface stuff

      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTOR, ALBEDO
      REAL(kind=dp), INTENT(IN)  :: UBRDF_FM ( MAX_USER_OBSGEOMS )

!  Multiplier

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Stream value

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Older variables
!  ---------------

!  transmittance factors for +/- eigenvalues for lowest layer

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGENNL

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Solution constants of integration multiplied by homogeneous solutions for downwelling at BOA

      REAL(kind=dp), INTENT(IN)  :: LCON_XVEC1NL
      REAL(kind=dp), INTENT(IN)  :: MCON_XVEC1NL

!  Particular downwelling solution at BOA

      REAL(kind=dp), INTENT(IN)  :: WLOWER1NL

!  Eigenvectors defined at user-defined stream angles
!  EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_OBSGEOMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!  @@@ NOT REQUIRED For MS-mode only
!  REAL(kind=dp), INTENT(IN)  :: U_WPOS1(MAX_USER_STREAMS,MAXLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (MAX_USER_OBSGEOMS,MAXLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_UP(MAX_USER_OBSGEOMS,MAXLAYERS)

!  Outputs
!  -------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_UP(MAX_USER_OBSGEOMS)

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_UP (0:MAXLAYERS,MAX_USER_OBSGEOMS)

!  local variables
!  ---------------

!  Reflectance integrand  a(j).x(j).I(-j)

      REAL(kind=dp) :: IDOWNSURF

!  help variables

      INTEGER       :: UM, N, N1
      REAL(kind=dp) :: LAYERSOURCE
      REAL(kind=dp) :: BOA_SOURCE  ( MAX_USER_OBSGEOMS )
      REAL(kind=dp) :: SHOM, PAR, HOM, KMULT
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

!  reflected multiple scatter intensity at user defined-angles
!  BRDF code added 4 May 2009

        IF ( DO_BRDF_SURFACE  ) THEN
          KMULT = SURFACE_FACTOR * IDOWNSURF
          DO UM = 1, N_USER_OBSGEOMS
            BOA_SOURCE(UM) = KMULT * UBRDF_FM(UM)
          ENDDO
        ELSE
          KMULT  = SURFACE_FACTOR * ALBEDO * IDOWNSURF
          DO UM = 1, N_USER_OBSGEOMS
            BOA_SOURCE(UM) = KMULT
          ENDDO
        ENDIF

!    MSMODE only ===> No Direct BOA source term
!        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
!          DO UM = 1, N_USER_OBSGEOMS
!            DIRECT_BOA_SOURCE(UM) = USER_DIRECT_BEAM(UM,IPARTIC)
!          ENDDO
!        ENDIF

      ENDIF

!  Add surface emission term if flagged
!  ********  Direct surface emission not included

!  Recursion Loop in Source function integration
!  =============================================

      DO UM = 1, N_USER_OBSGEOMS

!  Set the cumulative source term equal to BOA values
!  MSMODE only ===> No Direct BOA source term

         CUMSOURCE_UP_OLD = BOA_SOURCE(UM)
         IF ( DO_2S_LEVELOUT ) RADLEVEL_F_UP(NLAYERS,UM) = FLUXMULT * CUMSOURCE_UP_OLD

         DO N = NLAYERS, 1, -1

            N1 = N - 1

!  Homogeneous solutions.

            SHOM = LCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N) + &
                   MCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N)
            LAYERSOURCE = SHOM

!  Add thermal emission term (direct and diffuse)

            LAYERSOURCE = LAYERSOURCE + LAYER_TSUP_UP(UM,N)

!  Upward recursion for source function.

            CUMSOURCE_UP_NEW = LAYERSOURCE + T_DELT_USERM(N,UM)*CUMSOURCE_UP_OLD
            IF (DO_2S_LEVELOUT)RADLEVEL_F_UP(N1,UM) = FLUXMULT * CUMSOURCE_UP_NEW
            CUMSOURCE_UP_OLD = CUMSOURCE_UP_NEW

!  End layer loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term at TOA

         INTENSITY_F_UP(UM) = FLUXMULT * CUMSOURCE_UP_OLD

!  End recursion loop

      ENDDO

!  debug 28
!      write(*,*)INTENSITY_F_UP(1)

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_UPUSER_INTENSITY_THERMAL

!

SUBROUTINE TWOSTREAM_DNUSER_INTENSITY_THERMAL &
   ( MAXLAYERS, MAX_USER_OBSGEOMS,          & ! Dimensions
     DO_2S_LEVELOUT,                        & ! Inputs !@@ 2p1, 2p2
     NLAYERS, N_USER_OBSGEOMS,              & ! inputs !@@ 2p3 Greens
     FLUXMULT,                              & ! inputs
     T_DELT_USERM,                          & ! Inputs !@@ 2p3 Greens
     LCON, MCON, U_XPOS, U_XNEG,            & ! Inputs
     HMULT_1, HMULT_2, LAYER_TSUP_DN,       & ! Inputs
     INTENSITY_F_DN, RADLEVEL_F_DN )          ! Output !@@ 2p2

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: ZERO = 0.0_dp, ONE = 1.0_dp

!  Subroutine input arguments
!  --------------------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAX_USER_OBSGEOMS

!     ! @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL, INTENT(IN)        :: DO_2S_LEVELOUT

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, N_USER_OBSGEOMS

!  Multiplier

      REAL(kind=dp), INTENT(IN)  :: FLUXMULT

!  Transmittance factors for user-defined stream angles

      REAL(kind=dp), INTENT(IN)  :: T_DELT_USERM ( MAXLAYERS, MAX_USER_OBSGEOMS )

!  Older variables
!  ---------------

!  Solution constants of integration

      REAL(kind=dp), INTENT(IN)  :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: MCON(MAXLAYERS)

!  Eigenvectors defined at user-defined stream angles
!  EP for the positive KEIGEN values, EM for -ve KEIGEN

      REAL(kind=dp), INTENT(IN)  :: U_XPOS(MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: U_XNEG(MAX_USER_OBSGEOMS,MAXLAYERS)

!  Single-scatter Particular beam solution at user-defined angles
!  @@@ NOT REQUIRED For MS-mode only
!      REAL(kind=dp), INTENT(IN)  :: U_WNEG1(MAX_USER_OBSGEOMS,MAXLAYERS)

!  solution multipliers

      REAL(kind=dp), INTENT(IN)  :: HMULT_1 (MAX_USER_OBSGEOMS,MAXLAYERS)
      REAL(kind=dp), INTENT(IN)  :: HMULT_2 (MAX_USER_OBSGEOMS,MAXLAYERS)

!  Thermal layer source term

      REAL(kind=dp), INTENT(IN)  :: LAYER_TSUP_DN(MAX_USER_OBSGEOMS,MAXLAYERS)

!  Outputs
!  -------

!  User-defined solutions

      REAL(kind=dp), INTENT(OUT) :: INTENSITY_F_DN(MAX_USER_OBSGEOMS)

!  Fourier-component solutions at ALL levels (Optional Output)

      REAL(kind=dp), INTENT(OUT) :: RADLEVEL_F_DN (0:MAXLAYERS,MAX_USER_OBSGEOMS)

!  local variables
!  ---------------

!  help variables

      INTEGER       :: UM, N
      REAL(kind=dp) :: LAYERSOURCE
      REAL(kind=dp) :: SHOM
      REAL(kind=dp) :: CUMSOURCE_DN_OLD, CUMSOURCE_DN_NEW

!  Recursion Loop in Source function integration
!  =============================================

!  Initialize recursion for user-defined stream angles only

      DO UM = 1, N_USER_OBSGEOMS

         CUMSOURCE_DN_OLD = ZERO

         DO N = 1, NLAYERS

!  Homogeneous solutions.

            SHOM = LCON(N) * U_XNEG(UM,N) * HMULT_1(UM,N) + &
                   MCON(N) * U_XPOS(UM,N) * HMULT_2(UM,N)
            LAYERSOURCE = SHOM

!  Add thermal emission term (direct and diffuse)

            LAYERSOURCE = LAYERSOURCE + LAYER_TSUP_DN(UM,N)

!  Upward recursion for source function.

            CUMSOURCE_DN_NEW = LAYERSOURCE + T_DELT_USERM(N,UM)*CUMSOURCE_DN_OLD
            IF (DO_2S_LEVELOUT)RADLEVEL_F_DN(N,UM) = FLUXMULT * CUMSOURCE_DN_NEW
            CUMSOURCE_DN_OLD = CUMSOURCE_DN_NEW

!  End layer loop

         ENDDO

!  User-defined stream output, just set to the cumulative source term at BOA

         INTENSITY_F_DN(UM) = FLUXMULT * CUMSOURCE_DN_OLD

!  End recursion loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_DNUSER_INTENSITY_THERMAL

end module twostream_intensity_thermal_m
