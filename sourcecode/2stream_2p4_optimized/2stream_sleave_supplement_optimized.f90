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
! #       2STREAM_SLEAVE_MASTER (master)                        #
! #                                                             #
! ###############################################################

module twostream_sleave_supplement_m

use twostream_sleave_routines_m

!  INDWAT, MORCASIWAT Routines taken straight from Clark Weaver code
!  Compiled here by R. Spurr, 11 July 2012

!  get_fluorescence_755 Routines taken straight from Chris O'dell code
!  Compiled here by R. Spurr, 12 July 2012

public

contains

SUBROUTINE TWOSTREAM_SLEAVE_MASTER &
         ( MAXBEAMS,                                               & ! Dimension
           DO_ISOTROPIC, DO_FLUORESCENCE,                          & ! Flags
           NBEAMS,                                                 & ! Geometry
           SAL, CHL, WAV,                                          & ! Water-leaving
           FL_Wavelength, FL_Amplitude755, FL_GAUSSIANS, Fs755,    & ! Fluorescence
           SLTERM_ISOTROPIC, SLTERM_F_0 )                            ! Output

!  Prepares the Surface Leaving supplement, necessary for 2Stream.

   implicit none

!  Parameter

   integer, parameter :: fpk = KIND( 1.0D0 )

!  Input arguments
!  ===============

!  Dimension

      INTEGER, intent(in) :: MAXBEAMS

!  Isotropic flag

      LOGICAL, intent(in) :: DO_ISOTROPIC

!  Flo flag

      LOGICAL, intent(in) :: DO_FLUORESCENCE

!  Local angle control

      INTEGER, intent(in) :: NBEAMS

!  Water-leaving variables
!  -----------------------

!  Input Salinity in [ppt]

      REAL, intent(in) :: SAL

!  Input Chlorophyll concentration in [mg/M]

      REAL, intent(in) :: CHL

!  Input wavelenth in [Microns]

      REAL, intent(in) :: WAV

!  Fluorescence variables
!  ----------------------

!  Input wavelength in [nm]

      REAL(fpk), intent(in) :: FL_Wavelength

!  Input F755 Amplitude

      REAL(fpk), intent(in)  :: FL_Amplitude755

!  Input Gaussians

      REAL(fpk), intent(in) ::  FL_GAUSSIANS(3,2)

!  Fluorescence at 755 nm

      REAL(fpk), intent(in) :: Fs755(MAXBEAMS)

!  Output arguments
!  ================

!  No exact term, as 2Streams is a multiple scattering code.

!  Fourier components of Surface-leaving terms:
!    Every solar direction, SL-transmitted quadrature streams

      REAL(fpk), intent(out), dimension ( MAXBEAMS, 0:1 )   :: SLTERM_F_0

!  Isotropic Surface leaving term (if flag set)

      REAL(fpk), intent(out) :: SLTERM_ISOTROPIC (MAXBEAMS)

!  Other local variables
!  =====================

!  Water-leaving model. SINGLE PRECISION VARIABLES.

      REAL :: RW, REFR, REFI, N12, RWB

!  Solar spectral radiance model wavelength

      REAL(FPK) :: ssr_wvl

!  Fluorescence model

      INTEGER   :: IB, K
      REAL(FPK) :: FL_SunSpec, FsSum, Sleave
      REAL(FPK) :: lamda, sigma, arg, Gauss, var
      !REAL(FPK) :: solar_spec_irradiance

!  Main code
!  ---------

!  Zero the output

      SLTERM_ISOTROPIC  = 0.0_fpk
      SLTERM_F_0        = 0.0_fpk

!  Fluorescence
!  ============

      IF ( DO_FLUORESCENCE ) THEN

!  Get solar spectral irradiance, in (W m−2 μm−1), to normalize data

         !FL_SunSpec = 1.0d0  ! Temporary

         ssr_wvl = FL_Wavelength*1.0e-3_fpk !convert from nm to um
         FL_SunSpec = solar_spec_irradiance( ssr_wvl )

!  Apply Gaussians

         FsSum = 0.0_fpk
         do k = 1, 2
            lamda = FL_Gaussians(2,k)
            sigma = FL_Gaussians(3,k)
            var = 0.5_fpk/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            IF ( arg .lt. 88.0_fpk ) THEN
               Gauss = FL_Gaussians(1,k) * dexp ( - arg )
            ELSE
               Gauss = 0.0_fpk
            ENDIF
            FsSum = FsSum + Gauss
         enddo

!  For each Solar zenith angle

         DO IB = 1, NBEAMS

!  Assign output Fluorescence (Apply Amplitude)
!  Multiply by Fs755, and normalize to solar spectrum

            SLTERM_ISOTROPIC(IB) = FL_Amplitude755 * FsSum * Fs755(IB) / FL_SunSpec

!  End Beam loop

        ENDDO

      ENDIF

!  WATER-LEAVING
!  =============

      IF ( .NOT. DO_FLUORESCENCE ) THEN

!  INDWAT call. Uses single precision routine

        CALL INDWAT(WAV, SAL, refr, refi)

!  MORCASIWAT call (6S Eric Vermote)
!  Returns the ocean radiance/Irradiance ratio Rw

        CALL MORCASIWAT(WAV, CHL, RW, .false.)

!  Set the isotropic term
!  Code from Clark Weaver, assumes perfect Transmittance
!  Add change in solid angle from under to above to surface
!  that accounts for 1/(n12*n12) decrease in directional reflectivity

        IF ( DO_ISOTROPIC ) THEN
!           a   = 0.485
!           tds = 1.0
!           tdv = 1.0
!           n12 = refr*refr + refi*refi  ; n12 = sqrt(n12)
           n12 = refr*refr + refi*refi
!           Rwb=(1.0/(n12*n12))*tds*tdv*Rw/(1-a*Rw)
           Rwb=(1.0/n12)*Rw/(1-0.485*Rw)
           Sleave = DBLE(Rwb)
        ENDIF

!  Set output (same for all solar beams - this is temporary)

       SLTERM_ISOTROPIC(1:NBEAMS) = Sleave

!  PLACEHOLDERS for other Water-leaving options

      ENDIF

!  Finish routine

    RETURN
END SUBROUTINE TWOSTREAM_SLEAVE_MASTER

!  End module

END MODULE twostream_sleave_supplement_m

