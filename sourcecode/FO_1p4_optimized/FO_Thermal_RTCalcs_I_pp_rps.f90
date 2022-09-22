! ###########################################################
! #                                                         #
! #                    THE LIDORT FAMILY                    #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Author :      Robert. J. D. Spurr                      #
! #                                                         #
! #  Address :     RT Solutions, Inc.                       #
! #                9 Channing Street                        #
! #                Cambridge, MA 02138, USA                 #
! #                                                         #
! #  Tel:          (617) 492 1183                           #
! #  Email :        rtsolutions@verizon.net                 #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #                 FIRST-ORDER MODEL                       #
! #     (EXACT SINGLE-SCATTERING and DIRECT-THERMAL)        #
! #                                                         #
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3a, 29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.3b, 24 January  2013, BRDF/SL Supplements   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_Thermal_RTCalcs_I_PP_RPS_m

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities(I):

!  (1) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!  This will perform Enhanced-PS calculations (outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel LOS-path)

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!  Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!  Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!  Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!  Version 3,  29 October  2012, Extension to multiple geometries

!  For Thermal Emission sources, the subroutines are
!  DTE_Integral_I_UP_PP_RPS   (Upwelling only)
!  DTE_Integral_I_DN_PP_RPS   (Downwelling only)
!  DTE_Integral_I_UPDN_PP_RPS (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine DTE_Integral_I_UP_PP_RPS &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling,              & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,           & ! Inputs (Thermal)
     deltaus, omega, truncfac, Mu1,               & ! Inputs (Optical, Geometry)
     intensity_dta_up, intensity_dts )              ! Outputs

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  Extension to multiple geometries, 29 October 2012

!  No partials

   implicit none         

!  Parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM)

   logical, Intent(in) :: Do_Thermset

!  Flags

   logical, Intent(in) ::  DO_DELTAM_SCALING

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  deltaus, omega, truncfac

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Mu1 = cos(alpha_boa)

   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out) :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out) :: intensity_dts    ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Thermal setup

   real(fpk)  :: tcom(2,maxlayers)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, v
   real(fpk)  :: bb_inputn1, bb_inputn, help, tms, thermcoeffs
   real(fpk)  :: Mu1v, lostau, deltausn, lostrans_upn, tcom2n, t_mult_up1, t_mult_up0
   real(fpk)  :: cumsource_dste, CUMSOURCE_UP_OLD, CUMSOURCE_UP_NEW

   real(fpk), parameter  :: cutoff = 88.0d0
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1

!  Thermal setup factors
!  TMS, Initial set of thermal coefficients and TCOM variable

   if ( do_Thermset ) then
      bb_inputn1 = bb_input(0)
      do n = 1, nlayers
         tms = one - omega(n)
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
         endif
         bb_inputn = bb_input(n)
         thermcoeffs  = (bb_inputn-bb_inputn1) / deltaus(n)
         tcom(1,n) = bb_inputn1 * tms
         tcom(2,n) = thermcoeffs * tms
         bb_inputn1 = bb_inputn
      ENDDO
!      do_Thermset = .false.
   endif

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Mu1(v)

      Mu1v = Mu1(v)

!  Initialize

      nstart = nlayers

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

      do n = 1, nut-1
         lostau = deltaus(n) / Mu1v
         if ( lostau .lt. cutoff ) then
            lostrans_up(n) = exp( - lostau )
         else
            lostrans_up(n) = zero
         endif
         sources_up(n) = zero
      enddo
      do n = nut, nlayers
         deltausn = deltaus(n)
         lostau = deltausn / Mu1v
         if ( lostau .lt. cutoff ) then
            lostrans_upn = exp( - lostau )
         else
            lostrans_upn = zero
         endif
         lostrans_up(n) = lostrans_upn
         tcom2n = tcom(2,n)
         t_mult_up1 = tcom(1,n) + tcom2n * Mu1v
         t_mult_up0 = - t_mult_up1 - tcom2n * deltausn
         sources_up(n) = t_mult_up0 * lostrans_upn  + t_mult_up1
      enddo

!  Source function integration
!  ===========================

!  Start recursion ( For DSTE term, Use surface emissivity )

      CUMSOURCE_UP_OLD = zero
      CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY(v)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!  Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            LOSTRANS_UPN =  LOSTRANS_UP(N)
            CUMSOURCE_DSTE     = LOSTRANS_UPN * CUMSOURCE_DSTE
            CUMSOURCE_UP_NEW   = LOSTRANS_UPN * CUMSOURCE_UP_OLD + SOURCES_UP(N)
            CUMSOURCE_UP_OLD = CUMSOURCE_UP_NEW
         ENDDO
         INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP_OLD
         INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_I_UP_PP_RPS

!

subroutine DTE_Integral_I_DN_PP_RPS &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling,              & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, deltaus, omega, truncfac, Mu1,     & ! Inputs (Optical, Geometry)
     intensity_dta_dn )                             ! Outputs

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  Extension to multiple geometries, 29 October 2012

   implicit none         

!  Parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM)

   logical, Intent(in) ::  Do_Thermset

!  Flags

   logical, Intent(in) ::  DO_DELTAM_SCALING

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric thermal BB functions

   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Mu1 = cos(alpha_boa), required for the Regular PS only

   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Thermal setup

   real(fpk)  :: tcom(2,maxlayers)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, v
   real(fpk)  :: bb_inputn1, bb_inputn, help, tms, thermcoeffs
   real(fpk)  :: Mu1v, lostau, deltausn, lostrans_dnn, tcom2n, t_mult_dn1
   real(fpk)  :: CUMSOURCE_DN_OLD, CUMSOURCE_DN_NEW

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS

!  Thermal setup factors
!  TMS, Initial set of thermal coefficients and TCOM variable

  if ( do_Thermset ) then
      bb_inputn1 = bb_input(0)
      do n = 1, nlayers
         tms = one - omega(n)
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n)
            tms = tms / help
         endif  
         bb_inputn = bb_input(n)
         thermcoeffs  = (bb_inputn-bb_inputn1) / deltaus(n)
         tcom(1,n) = bb_inputn1 * tms
         tcom(2,n) = thermcoeffs * tms
         bb_inputn1 = bb_inputn
      ENDDO
!      do_Thermset = .false. 
   endif

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Mu1(v)

      Mu1v = Mu1(v)

!  Initialize

      nstart = nlayers

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

      do n = nlayers, nut+1, -1
         lostau = deltaus(n) / Mu1(v)
         if ( lostau .lt. cutoff ) then
            lostrans_dn(n) = exp( - lostau )
         else
            lostrans_dn(n) = zero
         endif
         sources_dn(n) = zero
      enddo
      do n = nut, 1, -1
         deltausn = deltaus(n)
         lostau = deltausn / Mu1(v)
         if ( lostau .lt. cutoff ) then
            lostrans_dnn = exp( - lostau )
         else
            lostrans_dnn = zero
         endif
         lostrans_dn(n) = lostrans_dnn
         tcom2n = tcom(2,n)
         t_mult_dn1   = tcom(1,n) - tcom2n * Mu1v
         sources_dn(n)  = - t_mult_dn1 * lostrans_dnn + t_mult_dn1 + tcom2n * deltausn
      enddo

!  Source function integration
!  ===========================

!  Start recursion

      CUMSOURCE_DN_OLD = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!  Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            CUMSOURCE_DN_NEW = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN_OLD
            CUMSOURCE_DN_OLD = CUMSOURCE_DN_NEW
         ENDDO
         INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN_OLD
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_I_DN_PP_RPS

subroutine DTE_Integral_I_UPDN_PP_RPS &
   ( maxgeoms, maxlayers, max_user_levels,                       & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_Thermset, do_deltam_scaling, & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels,                & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                          & ! Inputs (Thermal)
     deltaus, omega, truncfac, Mu1,                              & ! Inputs (Optical, Geometry)
     intensity_dta_up, intensity_dts, intensity_dta_dn )           ! Outputs

!  Stand alone routine for Upwelling and Downwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

!  This version, revised by R. Spurr, 01 June 2012
!  Extension to multiple geometries, 29 October 2012

!  No partials

   implicit none         

!  Parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

   logical, Intent(in) ::  Do_Thermset

!  flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING
   logical, Intent(in) ::  DO_DELTAM_SCALING

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Mu1 = cos(alpha_boa), required for the Regular PS only

   real(fpk), Intent(in) :: Mu1(maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts        ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels,maxgeoms )

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      call DTE_Integral_I_UP_PP_RPS &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling,              & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,           & ! Inputs (Thermal)
     deltaus, omega, truncfac, Mu1,               & ! Inputs (Optical, Geometry)
     intensity_dta_up, intensity_dts )              ! Outputs
   endif

   if ( do_dnwelling ) then
       call DTE_Integral_I_DN_PP_RPS &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling,              & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, deltaus, omega, truncfac, Mu1,     & ! Inputs (Optical, Geometry)
     intensity_dta_dn )                             ! Outputs
   endif

!  Finish

   return
end subroutine DTE_Integral_I_UPDN_PP_RPS

!  End module

end module FO_Thermal_RTCalcs_I_PP_RPS_m
