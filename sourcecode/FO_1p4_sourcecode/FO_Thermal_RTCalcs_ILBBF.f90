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

module FO_Thermal_RTCalcs_ILBBF_m

!  For a given wavelength, this Module will calculate First-Order
!    upwelling+downwelling Intensities(I) and Blackbody Jacobians

!     (1) For the Atmospheric and Surface Direct Thermal Emission (DTE) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!  This will perform Enhanced-PS calculations (outgoing LOS-path sphericity) 
!  This will perform Regular-PS  calculations (plane-parallel LOS-path)

!  HERE ARE THE ORIGINAL COMMENTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!    Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!    Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!    Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!    Version 3,  29 October  2012, Extension to multiple geometries

!  For Thermal Emission sources, the subroutines are
!       DTE_Integral_ILBBF_UP   (Upwelling only)
!       DTE_Integral_ILBBF_DN   (Downwelling only)
!       DTE_Integral_ILBBF_UPDN (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine DTE_Integral_ILBBF_UP &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,         & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_abbwf, do_sbbwf,          & ! Inputs (Flags)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,          & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                           & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
     Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, abbwf_up, sbbwf,            & ! Outputs (Main)
     cumsource_up, Tcom, dTcom1, dTcom2 )                           ! Outputs (Auxiliary)

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of radiance and blackbody Jacobians. 

!  This version, the first of its kind, programmed by R. Spurr 4/8/19.
!   Based on the DTE_Integral_I_UP subroutine from Version 1.4 of the FO CODE.

!  Control comes throught the two flags do_abbwf/do_sbbwf, which determine
!  calculations of the atmospheric and surface BB functions respectively.

!  Outputs are the arrays abbwf_up and sbbwf.  
  
!    No partials

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

!  Thermal setup flag (for TCOM)

   logical, Intent(inout) ::  Do_Thermset

!  4/8/19. New blackbody Jacobian flags

   logical, intent(in)  :: do_abbwf, do_sbbwf
   
!  flags

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_REGULAR_PS
   logical, Intent(in) ::  DO_ENHANCED_PS
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
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
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_up ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts    ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_up      ( 0:maxlayers,maxgeoms )

!  4/8/19. New blackbody Jacobian output

   real(fpk), Intent(Out)  :: abbwf_up (max_user_levels,maxgeoms,0:maxlayers)
   real(fpk), Intent(Out)  :: sbbwf    (max_user_levels,maxgeoms)
   
!  Thermal setup

   real(fpk), Intent(InOut)   :: Tcom  (maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom2(maxlayers,2)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

   real(fpk)  :: L_Solutions_fine (maxfinelayers,2)
   real(fpk)  :: Solutions_fine   (maxlayers,maxfinelayers)
   real(fpk)  :: Wtrans_fine      (maxlayers,maxfinelayers)

!  Source function integration results

   real(fpk)  :: L_sources_up     ( maxlayers, 2 )
   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, nj, v, jb
   logical    :: layermask_up(maxlayers)
   real(fpk)  :: help, tran, kn, ke, xjkn, cot_1, cot_2, L_source, L_cumsource
   real(fpk)  :: cumsource_dste,  tms, tms_odelt, cosmum, u1, u2, lostau

   real(fpk), parameter  :: cutoff = 88.0d0
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output

   CUMSOURCE_UP = zero ; INTENSITY_dta_up = zero ; INTENSITY_dts = zero
   abbwf_up = zero ; sbbwf = zero
   
!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1
   LAYERMASK_UP = .false.
   LAYERMASK_UP(NUT:NLAYERS) = .true.

!  Thermal setup factors
!     TMS,  TCOM and its derivatives

   if ( do_Thermset ) then
      tcom = zero ; dTcom1 = zero ; dTcom2 = zero
      do n = 1, nlayers
         tms = one - omega(n)
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n) ; tms = tms / help
         endif
         tms_odelt = tms / deltaus(n)
         dTcom1(n,1) = tms  ; dTcom1(n,2) = - tms_odelt
         dTcom2(n,1) = zero ; dTcom2(n,2) = + tms_odelt
         Tcom(n,1) = bb_input(n-1) * dTcom1(n,1) + bb_input(n) * dTcom2(n,1)
         Tcom(n,2) = bb_input(n-1) * dTcom1(n,2) + bb_input(n) * dTcom2(n,2)
      enddo
   endif

!  Start Geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_up = zero  ; sources_up = zero ; L_sources_up = zero

!  Criticality

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(v)
      if (nstart.lt.nlayers) LAYERMASK_UP(nstart+1:nlayers) = .false.

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced

      if ( do_RegPSorPP ) then
         Cosmum = Mu1(v) ; if ( doNadir(v) ) Cosmum = one
         DO n = 1, nlayers
            lostau = deltaus(n) / cosmum
            if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
            if ( layermask_up(n) ) then
               u1 = one - lostrans_up(n)  ; u2 = Cosmum - lostrans_up(n) * ( cosmum + deltaus(n) )
               sources_up(n)       = u1 * Tcom(n,1)   + u2 * TCom(n,2)
               L_sources_up(n,1) = u1 * dTcom1(n,1) + u2 * dTcom1(n,2)
               L_sources_up(n,2) = u1 * dTcom2(n,1) + u2 * dTcom2(n,2)
            endif
         enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n,v)
            if (  doNadir(v)  .and. layermask_up(n) ) then
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j,v) * kn
                  solutions_fine(n,j)   = Tcom(n,1)   + xjkn * Tcom(n,2)
                  L_solutions_fine(j,1) = dTcom1(n,1) + xjkn * dTcom1(n,2)
                  L_solutions_fine(j,2) = dTcom2(n,1) + xjkn * dTcom2(n,2)
                  wtrans_fine(n,j)    = kn * exp ( - xjkn )* wfine(n,j,v)
              enddo
            else if ( .not. doNadir(v) .and. layermask_up(n) ) then
               cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
               ke = raycon(v) * kn  ; lostau = ke * ( cot_2 - cot_1 )
               if ( lostau .lt. cutoff ) lostrans_up(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j,v) * kn
                  tran = exp ( - ke * ( cot_2 - cotfine(n,j,v) ) )
                  solutions_fine(n,j)   = Tcom(n,1)   + xjkn * Tcom(n,2)
                  L_solutions_fine(j,1) = dTcom1(n,1) + xjkn * dTcom1(n,2)
                  L_solutions_fine(j,2) = dTcom2(n,1) + xjkn * dTcom2(n,2)
                  wtrans_fine(n,j)    = ke * tran * csqfine(n,j,v) * wfine(n,j,v)
               enddo
            endif
            sources_up(n)     = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            L_sources_up(n,1) = dot_product(L_solutions_fine(1:nj,1),wtrans_fine(n,1:nj))
            L_sources_up(n,2) = dot_product(L_solutions_fine(1:nj,2),wtrans_fine(n,1:nj))
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion ( For DSTE term, Use surface emissivity )

      NC =  0
      CUMSOURCE_UP(NC,v) = zero
      CUMSOURCE_DSTE   = SURFBB * USER_EMISSIVITY(v)
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1
         
!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working upwards from NSTART to level NUT,
!     Check for updating the recursion

      DO UTA = N_USER_LEVELS, 1, -1
         NUT    = USER_LEVELS(UTA) + 1
         DO N = NSTART, NUT, -1
            NC = NLAYERS + 1 - N
            CUMSOURCE_DSTE     = LOSTRANS_UP(N) * CUMSOURCE_DSTE
            CUMSOURCE_UP(NC,V) = LOSTRANS_UP(N) * CUMSOURCE_UP(NC-1,V) + SOURCES_UP(N)
        ENDDO
         INTENSITY_DTA_UP(UTA,V) = CUMSOURCE_UP(NC,V)
         INTENSITY_DTS(UTA,V)    = CUMSOURCE_DSTE
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  SBBWF is trivial

      IF ( DO_SBBWF ) then
         DO UTA = N_USER_LEVELS, 1, -1
            SBBWF(UTA,V)    = INTENSITY_DTS(UTA,V) / SURFBB
         ENDDO
      ENDIF

!  ABBWF is more involved

      IF ( DO_ABBWF ) then
         DO JB = NLAYERS, 0, -1
            L_CUMSOURCE = zero ; nstart = nlayers ; nut_prev = nstart + 1
            DO UTA = N_USER_LEVELS, 1, -1
               NUT    = USER_LEVELS(UTA) + 1
               DO N = NSTART, NUT, -1
                  L_SOURCE = zero
                  IF ( N.EQ.JB )   L_SOURCE = L_SOURCES_UP(N,2)
                  IF ( N.EQ.JB+1 ) L_SOURCE = L_SOURCES_UP(N,1)
                  L_CUMSOURCE = LOSTRANS_UP(N) * L_CUMSOURCE + L_SOURCE
               ENDDO
               ABBWF_UP(UTA,V,JB) = L_CUMSOURCE
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
               NUT_PREV = NUT
            ENDDO
         ENDDO
      ENDIF

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_ILBBF_UP

!

subroutine DTE_Integral_ILBBF_DN &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,     & ! Inputs (dimensioning)
     Do_Thermset, do_deltam_scaling, do_abbwf,                & ! Inputs (Flags)
     do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,      & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
     BB_input, extinction, deltaus, omega, truncfac,          & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,             & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                    & ! Inputs (Geometry)
     intensity_dta_dn, abbwf_dn,                              & ! Outputs (Main)
     cumsource_dn, Tcom, dTcom1, dTcom2 )                       ! Outputs (Auxiliary)

!  Stand alone routine for Upwelling Direct-thermal-emission (DTE)
!    computation of radiance and blackbody Jacobians. 

!  This version, the first of its kind, programmed by R. Spurr 4/8/19.
!   Based on the DTE_Integral_I_DN subroutine from Version 1.4 of the FO CODE.

!  Control comes throught the flag do_abbwf, which determines
!  calculation of the atmospheric  BB functions.

!  Output is the array abbwf_dn  

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

   logical, Intent(inout) ::  Do_Thermset

!  4/8/19. New blackbody Jacobian flag

   logical, intent(in)  :: do_abbwf

!  flags

   logical, Intent(in) ::  DO_DELTAM_SCALING
   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_REGULAR_PS
   logical, Intent(in) ::  DO_ENHANCED_PS
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
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
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dta_dn ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_dn     ( 0:maxlayers,maxgeoms )
   
!  4/8/19. New blackbody Jacobian output

   real(fpk), Intent(Out)  :: abbwf_dn (max_user_levels,maxgeoms,0:maxlayers)
   
!  Thermal setup

   real(fpk), Intent(InOut)   :: Tcom  (maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom2(maxlayers,2)

!  LOCAL
!  -----

!  Local solutions (enhanced_ps case)

   real(fpk)  :: L_Solutions_fine (maxfinelayers,2)
   real(fpk)  :: Solutions_fine   (maxlayers,maxfinelayers)
   real(fpk)  :: Wtrans_fine      (maxlayers,maxfinelayers)

!  Source function integration results

   real(fpk)  :: L_sources_dn     ( maxlayers, 2 )
   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Regular_PS or plane-parallel flag

   logical    :: do_RegPSorPP

!  Help

   integer    :: n, uta, nstart, nc, nut, nut_prev, j, nj, v, jb
   logical    :: layermask_dn(maxlayers)
   real(fpk)  :: help, tran, kn, ke, xjkn, trand, argum(maxfinelayers), tms_odelt
   real(fpk)  :: cot_1, cot_2, rdiff, cot_c, tms, d1, d2, cosmum, lostau, L_SOURCE, L_CUMSOURCE

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Zero the output and the local sources

   CUMSOURCE_DN = zero ; INTENSITY_DTA_DN = zero ; abbwf_dn = zero

!  Regular_PS or plane-parallel flag

   do_RegPSorPP = (do_regular_ps .or. do_PlanPar)

!  Bookkeeping

   NUT = USER_LEVELS(N_USER_LEVELS) + 1
   IF ( NUT > NLAYERS ) NUT = NLAYERS
   LAYERMASK_DN = .false.
   LAYERMASK_DN(1:NUT) = .true.

!  Thermal setup factors
!     TMS,  TCOM and its derivatives

   if ( do_Thermset ) then
      tcom = zero ; dTcom1 = zero ; dTcom2 = zero
      do n = 1, nlayers
         tms = one - omega(n)
         if ( do_deltam_scaling ) then
            help = one - truncfac(n) * omega(n) ; tms = tms / help
         endif
         tms_odelt = tms / deltaus(n)
         dTcom1(n,1) = tms  ; dTcom1(n,2) = - tms_odelt
         dTcom2(n,1) = zero ; dTcom2(n,2) = + tms_odelt
         Tcom(n,1) = bb_input(n-1) * dTcom1(n,1) + bb_input(n) * dTcom2(n,1)
         Tcom(n,2) = bb_input(n-1) * dTcom1(n,2) + bb_input(n) * dTcom2(n,2)
      enddo
   endif

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

!  Zero the local sources

      lostrans_dn  = zero ; sources_dn   = zero

!  Criticality

      nstart = nlayers ; if (Ncrit(v).ne.0) nstart = nCrit(V)
      if (nstart.lt.nlayers) LAYERMASK_DN(nstart+1:nlayers) = .false.

!  Plane/Parallel or Regular-PS Layer integrated source terms
!  ==========================================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced

      if ( do_RegPSorPP ) then
         Cosmum = Mu1(v) ; if ( doNadir(v) ) Cosmum = one
         DO n = 1, nlayers
            lostau = deltaus(n) / cosmum
            if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
            if ( layermask_dn(n) ) then
               d1 = one - lostrans_dn(n)  ; d2 = deltaus(n) - cosmum * d1
               sources_dn(n)     = d1 * Tcom(n,1)   + d2 * TCom(n,2)
               L_sources_dn(n,1) = d1 * dTcom1(n,1) + d2 * dTcom1(n,2)
               L_sources_dn(n,2) = d1 * dTcom2(n,1) + d2 * dTcom2(n,2)
            endif
         enddo
      endif

!  LOS-spherical Layer integrated source terms
!  ===========================================

!  Bug Fixed 23 January 2013 (nadir case). Old code commented out and replaced. Defined argum(j)

      if ( do_enhanced_ps ) then
         do n = nlayers, 1, -1
            kn = extinction(n) ; nj = nfinedivs(n,v)
            if (  doNadir(v) .and. layermask_dn(n) ) then
               rdiff = radii(n-1) - radii(n) ; if ( n.eq.NCrit(v) ) rdiff = radii(n-1) - RadCrit(v)
               trand = one ; if ( n.eq.NCrit(v) ) trand = exp ( - kn * (RadCrit(v) - radii(n) ) )
               lostau = deltaus(n)
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               do j = 1, nj
                  argum(j) = rdiff - xfine(n,j,v) ; xjkn = xfine(n,j,v) * kn
                  solutions_fine(n,j)   = Tcom(n,1)   + xjkn * Tcom(n,2)
                  L_solutions_fine(j,1) = dTcom1(n,1) + xjkn * dTcom1(n,2)
                  L_solutions_fine(j,2) = dTcom2(n,1) + xjkn * dTcom2(n,2)
                  wtrans_fine(n,j)    = kn * exp ( - kn * argum(j) ) * wfine(n,j,v)
               enddo
            else if ( .not. doNadir(v) .and. layermask_dn(n)  ) then
               cot_2 = cota(n-1,v) ; cot_1 = cota(n,v)
               cot_c = cot_1     ; if ( n.eq.NCrit(v) ) cot_c = CotCrit(v)
               ke = raycon(v) * kn  ; lostau = ke * ( cot_2 - cot_1 )
               trand = one  ; if ( n.eq.NCrit(v) ) trand = exp ( - ke * ( CotCrit(v) - cot_1 ) )
               if ( lostau .lt. cutoff ) lostrans_dn(n) = exp( - lostau )
               do j = 1, nj
                  xjkn = xfine(n,j,v) * kn
                  tran = exp ( - ke * ( cotfine(n,j,v) - cot_1 ) )   !  Down
                  solutions_fine(n,j)   = Tcom(n,1)   + xjkn * Tcom(n,2)
                  L_solutions_fine(j,1) = dTcom1(n,1) + xjkn * dTcom1(n,2)
                  L_solutions_fine(j,2) = dTcom2(n,1) + xjkn * dTcom2(n,2)
                  wtrans_fine(n,j)    = ke * tran * csqfine(n,j,v) * wfine(n,j,v)
               enddo
            endif
            sources_dn(n) = dot_product(solutions_fine(n,1:nj),wtrans_fine(n,1:nj))
            L_sources_dn(n,1) = dot_product(L_solutions_fine(1:nj,1),wtrans_fine(n,1:nj))
            L_sources_dn(n,2) = dot_product(L_solutions_fine(1:nj,2),wtrans_fine(n,1:nj))
            if ( n.eq.NCrit(v) ) then
              sources_dn(n)       = sources_dn(n) * trand         !@@ Robfix
              L_sources_dn(n,1:2) = L_sources_dn(n,1:2) * trand 
            endif
         enddo
      endif

!  Source function integration
!  ===========================

!  start recursion

      NC =  0
      CUMSOURCE_DN(NC,V) = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!     NLEVEL = Layer index for given optical depth
!     Cumulative source terms : Loop over layers working Downn from NSTART to NUT
!     Check for dndating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            NC = N
            CUMSOURCE_DN(NC,V) = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN(NC-1,V)
         ENDDO
         INTENSITY_DTA_DN(UTA,V) = CUMSOURCE_DN(NC,V)
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  ABBWF is more involved

      IF ( DO_ABBWF ) then
         DO JB = 0, NLAYERS
            L_CUMSOURCE = zero ; nstart = 1 ; NUT_PREV = NSTART - 1
            DO UTA = 1, N_USER_LEVELS
               NUT    = USER_LEVELS(UTA)
               DO N = NSTART, NUT
                  L_SOURCE = zero
                  IF ( N.EQ.JB )   L_SOURCE = L_SOURCES_DN(N,2)
                  IF ( N.EQ.JB+1 ) L_SOURCE = L_SOURCES_DN(N,1)
                  L_CUMSOURCE = LOSTRANS_DN(N) * L_CUMSOURCE + L_SOURCE
               ENDDO
               ABBWF_DN(UTA,V,JB) = L_CUMSOURCE
               IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
               NUT_PREV = NUT
            ENDDO
         ENDDO
      ENDIF

!  End geometry loop

   enddo

!  Finish

   return
end subroutine DTE_Integral_ILBBF_DN



subroutine DTE_Integral_ILBBF_UPDN   &
   ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,             & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, do_Thermset, do_deltam_scaling,      & ! Inputs (Flags)
     do_abbwf, do_sbbwf, do_PlanPar, do_regular_ps, do_enhanced_ps,   & ! Inputs (Flags)
     doNadir, ngeoms, nlayers, nfinedivs, n_user_levels, user_levels, & ! Inputs (control output)
     bb_input, surfbb, user_emissivity,                               & ! Inputs (Thermal)
     extinction, deltaus, omega, truncfac,                            & ! Inputs (Optical)
     Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,                     & ! Inputs (Geometry)
     cota, xfine, wfine, csqfine, cotfine,                            & ! Inputs (Geometry)
     intensity_dta_up, intensity_dts, abbwf_up, sbbwf,                & ! Outputs (Main, up)
     intensity_dta_dn, abbwf_dn,                                      & ! Outputs (Main, down)
     cumsource_up, cumsource_dn, Tcom, dTcom1, dTcom2 )                 ! Outputs (Auxiliary)
  
!  Stand alone routine for Upwelling and Downwelling Direct-thermal-emission (DTE)
!    computation of radiance and blackbody Jacobians. 

!  This version, the first of its kind, programmed by R. Spurr 4/8/19.
!   Based on the DTE_Integral_I_UPDN subroutine from Version 1.4 of the FO CODE.

!  Control comes throught the two flags do_abbwf/do_sbbwf, which determine
!  calculations of the atmospheric and surface BB functions respectively.

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

   logical, Intent(inout) ::  Do_Thermset

!  4/8/19. New blackbody Jacobian flags

   logical, intent(in)  :: do_abbwf, do_sbbwf

!  flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING
   logical, Intent(in) ::  DO_DELTAM_SCALING

   logical, Intent(in) ::  DO_PLANPAR
   logical, Intent(in) ::  DO_REGULAR_PS
   logical, Intent(in) ::  DO_ENHANCED_PS
   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NGEOMS, NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  optical inputs
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
!    Mu1 = cos(alpha_boa), required for the Regular PS only

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu1(maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxlayers,maxfinelayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxlayers,maxfinelayers,maxgeoms)

!  outputs
!  -------

!  Original output
   
   real(fpk), Intent(Out)  :: intensity_dta_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dts        ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_up         ( 0:maxlayers,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dta_dn     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: cumsource_dn         ( 0:maxlayers,maxgeoms )

!  4/8/19. New blackbody Jacobian output

   real(fpk), Intent(Out)  :: abbwf_up (max_user_levels,maxgeoms,0:maxlayers)
   real(fpk), Intent(Out)  :: abbwf_dn (max_user_levels,maxgeoms,0:maxlayers)
   real(fpk), Intent(Out)  :: sbbwf    (max_user_levels,maxgeoms)

!  Thermal setup

   real(fpk), Intent(InOut)   :: Tcom  (maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom1(maxlayers,2)
   real(fpk), Intent(InOut)   :: dTcom2(maxlayers,2)

!  Upwelling
!  ---------

   if ( do_upwelling ) then
      call DTE_Integral_ILBBF_UP &
        ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,      & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_abbwf, do_sbbwf,          & ! Inputs (Flags)
          do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,          & ! Inputs (Flags)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,      & ! Inputs (control output)
          bb_input, surfbb, user_emissivity,                           & ! Inputs (Thermal)
          extinction, deltaus, omega, truncfac,                        & ! Inputs (Optical)
          Mu1, NCrit, Raycon, cota, xfine, wfine, csqfine, cotfine,    & ! Inputs (Geometry)
          intensity_dta_up, intensity_dts, abbwf_up, sbbwf,            & ! Outputs (Main)
          cumsource_up, Tcom, dTcom1, dTcom2 )                           ! Outputs (Auxiliary)
      do_Thermset = .false.
   endif

   if ( do_dnwelling ) then
       call DTE_Integral_ILBBF_DN &
        ( maxgeoms, maxlayers, maxfinelayers, max_user_levels,    & ! Inputs (dimensioning)
          Do_Thermset, do_deltam_scaling, do_abbwf,                & ! Inputs (Flags)
          do_PlanPar, do_regular_ps, do_enhanced_ps, doNadir,      & ! Inputs (Flags)
          ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
          BB_input, extinction, deltaus, omega, truncfac,          & ! Inputs (Optical)
          Mu1, NCrit, RadCrit, CotCrit, Raycon, radii,             & ! Inputs (Geometry)
          cota, xfine, wfine, csqfine, cotfine,                    & ! Inputs (Geometry)
          intensity_dta_dn, abbwf_dn,                              & ! Outputs (Main)
          cumsource_dn, Tcom, dTcom1, dTcom2 )                       ! Outputs (Auxiliary)
      do_Thermset = .false.
   endif

!  Finish

   return
end subroutine DTE_Integral_ILBBF_UPDN

!  End module

end module FO_Thermal_RTCalcs_ILBBF_m

