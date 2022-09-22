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

module FO_Thermal_RTCalcs_I_EPS_m

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
!  DTE_Integral_I_UP_EPS   (Upwelling only)
!  DTE_Integral_I_DN_EPS   (Downwelling only)
!  DTE_Integral_I_UPDN_EPS (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine DTE_Integral_I_UP_EPS &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,    & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, doNadir,                & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                      & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                   & ! Inputs (Optical)
     NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts )                         ! Outputs

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

!  No partials

   implicit none         

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

   logical, Intent(in) :: Do_Thermset

!  Flags

   logical, Intent(in) :: DO_DELTAM_SCALING
   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) :: NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere extinction and deltaus

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfinelayers,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

   real(fpk)  :: Solutions_fine (maxfinelayers)
   real(fpk)  :: Wtrans_fine    (maxfinelayers)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Thermal setup

   real(fpk)  :: tcom(2,maxlayers)

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, j, nj, v, nCritv
   logical    :: doNadirv
   real(fpk)  :: rayconv, bb_inputn1, bb_inputn, help, tms, thermcoeffs
   real(fpk)  :: cot_1, lostau, cot_2, kn, tcom1n, tcom2n, xjkn, ke, tran
   real(fpk)  :: cumsource_dste, CUMSOURCE_UP_OLD, CUMSOURCE_UP_NEW, LOSTRANS_UPN

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

      doNadirv = doNadir(v)

!  Raycon(v)

      if ( .not. doNadirv ) rayconv = raycon(v)

!  Criticality

      nstart = nlayers ; nCritv = Ncrit(v); if (Ncritv .ne. 0) nstart = nCritv

!  LOS-spherical Layer integrated source terms
!  ===========================================

         if (.not. doNadirv ) cot_1 = cota(nlayers,v)
         do n = nlayers, nstart+1, -1
            if (  doNadirv ) then 
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
               else
                  lostrans_up(n) = zero
               endif
            else
               cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
               lostau = rayconv * extinction(n) * ( cot_2 - cot_1 )
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
                              else
                  lostrans_up(n) = zero
               endif
               cot_1 = cot_2
            endif
            sources_up(n) = zero
         enddo
         do n = nstart, nut, -1
            kn = extinction(n) ; nj = nfinedivs(n,v)
            tcom1n = tcom(1,n) ; tcom2n = tcom(2,n)
            if (  doNadirv ) then
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
               else
                  lostrans_up(n) = zero
               endif
               do j = 1, nj
                  xjkn = xfine(j,n,v) * kn
                  solutions_fine(j) = tcom1n + xjkn * tcom2n
                  wtrans_fine(j)    = kn * exp ( - xjkn )* wfine(j,n,v)
              enddo
            else
               cot_2 = cota(n-1,v)
               ke = rayconv * kn  ; lostau = ke * ( cot_2 - cot_1 )
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
               else
                  lostrans_up(n) = zero
               endif
               do j = 1, nj
                  xjkn = xfine(j,n,v) * kn
                  tran = exp ( - ke * ( cot_2 - cotfine(j,n,v) ) )
                  solutions_fine(j) = tcom1n + xjkn * tcom2n
                  wtrans_fine(j)    = ke * tran * csqfine(j,n,v) * wfine(j,n,v)
               enddo
               cot_1 = cot_2
            endif
            sources_up(n) = dot_product(solutions_fine(1:nj),wtrans_fine(1:nj))
         enddo
         do n = nut-1, 1, -1
            if (  doNadirv ) then 
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
               else
                  lostrans_up(n) = zero
               endif
            else
               cot_2 = cota(n-1,v)
               lostau = rayconv * extinction(n) * ( cot_2 - cot_1 )
               if ( lostau .lt. cutoff ) then
                  lostrans_up(n) = exp( - lostau )
                              else
                  lostrans_up(n) = zero
               endif
               cot_1 = cot_2
            endif
            sources_up(n) = zero
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
            CUMSOURCE_DSTE   = LOSTRANS_UPN * CUMSOURCE_DSTE
            CUMSOURCE_UP_NEW = LOSTRANS_UPN * CUMSOURCE_UP_OLD + SOURCES_UP(N)
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
end subroutine DTE_Integral_I_UP_EPS

!

subroutine DTE_Integral_I_DN_EPS &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,    & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, doNadir,                & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,         & ! Inputs (Optical)
     NCrit, RadCrit, CotCrit, Raycon, radii,                 & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                   & ! Inputs (Geometry)
     intensity_dta_dn )                                        ! Outputs

!  Stand alone routine for Downwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

   implicit none         

!  Parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

   logical, Intent(in) ::  Do_Thermset

!  Flags

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric thermal BB functions

   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfinelayers,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Thermal setup

   real(fpk)  :: tcom(2,maxlayers)

!  Local solutions (enhanced_ps case)

   real(fpk)  :: Solutions_fine (maxfinelayers)
   real(fpk)  :: Wtrans_fine    (maxfinelayers)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, j, nj, v, nCritv
   logical    :: doNadirv
   real(fpk)  :: rayconv, bb_inputn1, bb_inputn, help, tms, thermcoeffs
   real(fpk)  :: RadCritv, CotCritv
   real(fpk)  :: cot_1, lostau, cot_2, kn, tcom1n, tcom2n, radiin, radiin1, rdiff
   real(fpk)  :: trand,  xfjnv, ke, tran
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

      doNadirv = doNadir(v)

!  Raycon(v)

      if ( .not. doNadir(v) ) rayconv = raycon(v)

!  RadCrit(v)

      if ( doNadir(v) ) RadCritv = RadCrit(v)

!  RadCrit(v)

      if ( .not. doNadir(v) ) CotCritv = CotCrit(v)

!  Criticality

      nstart = nlayers ; nCritv = Ncrit(v); if (Ncritv .ne. 0) nstart = nCritv

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( .not. doNadirv ) cot_1 = cota(nlayers,v)
      do n = nlayers, min(nut,nstart)+1, -1
         if (  doNadir(v) ) then
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) then
               lostrans_dn(n) = exp( - lostau )
            else
               lostrans_dn(n) = zero
            endif
         else    
            cot_2 = cota(n-1,v)
            lostau = rayconv * extinction(n) * ( cot_2 - cot_1 )
            if ( lostau .lt. cutoff ) then
               lostrans_dn(n) = exp( - lostau )
            else
               lostrans_dn(n) = zero
            endif
            cot_1 = cot_2
         endif
         sources_dn(n) = zero
      enddo

      if ( doNadirv ) radiin = radii(min(nut,nstart))
      do n = min(nut,nstart), 1, -1
         kn = extinction(n) ; nj = nfinedivs(n,v)
         tcom1n = tcom(1,n) ; tcom2n = tcom(2,n)
         if (  doNadir(v) ) then
            lostau = deltaus(n)
            if ( lostau .lt. cutoff ) then
               lostrans_dn(n) = exp( - lostau )
            else
               lostrans_dn(n) = zero
            endif
            radiin1 = radii(n-1)
            rdiff = radiin1 - radiin ; if ( n .eq. NCritv ) rdiff = radiin1 - RadCritv
            trand = one ; if ( n .eq. NCritv ) trand = exp ( - kn * ( RadCritv - radiin ) )
            do j = 1, nj
               xfjnv = xfine(j,n,v)
               solutions_fine(j) = tcom1n + xfjnv * kn * tcom2n
               wtrans_fine(j)    = kn * exp ( - kn * (rdiff - xfjnv) ) * wfine(j,n,v)
            enddo
         else
            cot_2 = cota(n-1,v)
            ke = rayconv * kn  ; lostau = ke * ( cot_2 - cot_1 )
            if ( lostau .lt. cutoff ) then
               lostrans_dn(n) = exp( - lostau )
            else
               lostrans_dn(n) = zero
            endif
            trand = one  ; if ( n .eq. NCritv ) trand = exp ( - ke * ( CotCritv - cot_1 ) )
            do j = 1, nj
               tran = exp ( - ke * ( cotfine(j,n,v) - cot_1 ) )   !  Down
               solutions_fine(j) = tcom1n + xfine(j,n,v) * kn * tcom2n
               wtrans_fine(j)    = ke * tran * csqfine(j,n,v) * wfine(j,n,v)
            enddo
            cot_1 = cot_2
         endif
         sources_dn(n) = dot_product(solutions_fine(1:nj),wtrans_fine(1:nj))
         if ( n .eq. NCritv ) sources_dn(n) = sources_dn(n) * trand         !@@ Robfix
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
end subroutine DTE_Integral_I_DN_EPS

subroutine DTE_Integral_I_UPDN_EPS &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,         & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_Thermset, do_deltam_scaling,  & ! Inputs (Flags)
     doNadir,                                                     & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                           & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
     NCrit, RadCrit, CotCrit, Raycon, radii,                      & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                        & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, intensity_dta_dn )            ! Outputs

!  Stand alone routine for Upwelling and Downwelling Direct-thermal-emission (DTE)
!  Computation of radiance. Can be derived from input Planck functions.

   implicit none         

!  Parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs
!  ======

!  Dimensions

   integer, Intent(in) :: maxgeoms
   integer, Intent(in) :: maxlayers
   integer, Intent(in) :: maxfinelayers
   integer, Intent(in) :: max_user_levels

!  Thermal setup flag (for TCOM1)

   logical, Intent(in) ::  Do_Thermset

!  Flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING
   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: OMEGA       ( MAXLAYERS )
   real(fpk), Intent(in) :: TRUNCFAC    ( MAXLAYERS )

!  Atmospheric BB functions and Surface BB and emissivity

   real(fpk), Intent(in) :: SURFBB, USER_EMISSIVITY(MAXGEOMS)
   real(fpk), Intent(in) :: BB_INPUT (0:MAXLAYERS)

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfinelayers,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfinelayers,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts        ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels,maxgeoms )

   if ( do_upwelling ) then
      call DTE_Integral_I_UP_EPS &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,    & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, doNadir,                & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                      & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                   & ! Inputs (Optical)
     NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts )                         ! Outputs
   endif

   if ( do_dnwelling ) then
       call DTE_Integral_I_DN_EPS &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,    & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, doNadir,                & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,         & ! Inputs (Optical)
     NCrit, RadCrit, CotCrit, Raycon, radii,                 & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                   & ! Inputs (Geometry)
     intensity_dta_dn )                                        ! Outputs
   endif

!  Finish

   return
end subroutine DTE_Integral_I_UPDN_EPS

!  End module

end module FO_Thermal_RTCalcs_I_EPS_m
