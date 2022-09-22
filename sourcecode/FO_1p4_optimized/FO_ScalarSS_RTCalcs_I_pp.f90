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
! #     FIRST-ORDER SCALAR MODEL (EXACT SINGLE SCATTERING)  #
! #                                                         #
! #  This Version :   1.3 F90                               #
! #  Release Date :   March 2013                            #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Multiple geometries   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_ScalarSS_RTCalcs_I_PP_m

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities (I):

!  (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform plane-parallel calculations

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!  Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!  Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!  Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!  Version 3,  29 October  2012, Extension to multiple geometries
!  Version 4,  31 July     2013, Lattice Multi-geometry

!  For Solar sources, the subroutines are
!       SS_Integral_I_UP_PP   (Upwelling only)
!       SS_Integral_I_DN_PP   (Downwelling only)
!       SS_Integral_I_UPDN_PP (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SS_Integral_I_UP_PP &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     reflec, deltaus, exactscat_up, flux,         & ! Inputs (Optical)
     Mu0, Mu1,                                    & ! Inputs (Geometry)
     intensity_up, intensity_db )                   ! Outputs

!  Stand-alone routine for Upwelling Solar-beam Single-scatter (SS)
!  Computation of radiance. Inputs: geometry, optical properties

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

!  Numbers

   integer, Intent(in) :: NLAYERS
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS, MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Geometrical inputs
!  ------------------

!  Mu0 = cos(theta_boa), required for PP
!  Mu1 = cos(alpha_boa), required for PP

   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)

!  Solutions

   real(fpk)  :: Solutions (maxlayers)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Average secant type arrays (Plane-Parallel)

   real(fpk)  :: factor1       ( maxlayers )
   real(fpk)  :: factor2       ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, v
   real(fpk)  :: Mu0v, Mu1v
   real(fpk)  :: cumtau, Attenuationsn1, Attenuationsn
   real(fpk)  :: sumd, multiplier, lostau
   real(fpk)  :: cumsource_db, CUMSOURCE_UP_OLD, CUMSOURCE_UP_NEW, LOSTRANS_UPN

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1

!  Start Geometry loop

   do v = 1, ngeoms

!  Mu0(v), Mu1(v)

      Mu0v = Mu0(v)
      Mu1v = Mu1(v)

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      nstart = nlayers

!  Attenuations to End points (including TOA)
!  MUST go all the way to NLAYERS (surface term required)

      Attenuations(0) = one
      cumtau = zero
      do n = 1, nlayers
         cumtau = cumtau+deltaus(n)
         sumd = cumtau / Mu0v
         if (sumd .lt. cutoff ) then
            Attenuations(n) = exp( - sumd )
         else
            Attenuations(n) = zero
         endif
      enddo

!  Plane/Parallel
!  Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      if ( Mu1v .eq. zero ) then
         factor2(1:nlayers) = zero
         lostrans_up(1:nlayers) = zero
         Attenuationsn1 = Attenuations(0)
         do n = 1, nut-1
            Solutions(n) = exactscat_up(n,v) * Attenuationsn1
            Attenuationsn1 = Attenuations(n)
            factor1(n) = zero
         enddo
         do n = nut, nlayers
            Solutions(n) = exactscat_up(n,v) * Attenuationsn1
            Attenuationsn = Attenuations(n)
            Attenuationsn1 = Attenuationsn
            if ( Attenuationsn1 .ne. zero ) then
               factor1(n) = Attenuationsn/Attenuationsn1
               nstart = n
            else
               factor1(n) = zero
            endif
         enddo
      else
         Attenuationsn1 = Attenuations(0)
         do n = 1, nut-1
            lostau = deltaus(n) / Mu1v
            Solutions(n) = exactscat_up(n,v) * Attenuationsn1
            if ( lostau .lt. cutoff ) then
               lostrans_up(n) = exp( - lostau )
            else
               lostrans_up(n) = zero
            endif
            Attenuationsn1 = Attenuations(n)
            factor1(n) = zero
            factor2(n) = zero
         enddo
         do n = nut, nlayers
            lostau = deltaus(n) / Mu1v
            Solutions(n) = exactscat_up(n,v) * Attenuationsn1
            if ( lostau .lt. cutoff ) then
               lostrans_up(n) = exp( - lostau )
            else
               lostrans_up(n) = zero
            endif
            Attenuationsn = Attenuations(n)
            if ( Attenuationsn1 .ne. zero ) then
               factor1(n) = Attenuationsn/Attenuationsn1
               factor2(n) = Mu1v/Mu0v
               nstart = n
            else
               factor1(n) = zero
               factor2(n) = zero
            endif
            Attenuationsn1 = Attenuationsn
         enddo
      endif

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel Layer integrated source terms

      do n = nlayers, nstart+1, -1
         sources_up(n) = zero
      enddo
      do n = nstart, nut, -1
         multiplier = ( one - Factor1(n)*lostrans_up(n) ) / (factor2(n) + one)
         sources_up(n) = solutions(n) * multiplier
      enddo
      do n = nut-1, 1, -1
         sources_up(n) = solutions(n)
      enddo

!  Source function integration
!  ===========================

!  Start recursion ( For Direct Beam, use PI.mu0.R.Atten )

      CUMSOURCE_UP_OLD = zero
      CUMSOURCE_DB     = 4.0_fpk * Mu0v * REFLEC(v) * Attenuationsn
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
            CUMSOURCE_DB     = LOSTRANS_UPN * CUMSOURCE_DB
            CUMSOURCE_UP_NEW = LOSTRANS_UPN * CUMSOURCE_UP_OLD + SOURCES_UP(N)
            CUMSOURCE_UP_OLD = CUMSOURCE_UP_NEW
         ENDDO
         INTENSITY_UP(UTA,V) = FLUX * CUMSOURCE_UP_OLD
         INTENSITY_DB(UTA,V) = FLUX * CUMSOURCE_DB
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
         NUT_PREV = NUT
      ENDDO

!  End Geometry Loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_UP_PP

!

subroutine SS_Integral_I_DN_PP &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     deltaus, exactscat_dn, flux, Mu0, Mu1,       & ! Inputs (Optical,Geometry)
     intensity_dn )                                 ! Outputs

!  Stand-alone routine for Downwelling Solar-beam Single-scatter (SS)
!  Computation of radiance. Inputs: geometry, optical properties

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

!  Numbers

   integer, Intent(in) ::  NLAYERS
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux

   real(fpk), Intent(in) ::  FLUX

!  Geometrical inputs
!  ------------------

!  Mu0 = cos(theta_boa), required for PP
!  Mu1 = cos(alpha_boa), required for PP

   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out) :: intensity_dn     ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuations      (0:maxlayers)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Average secant type variables (Plane-parallel)

   real(fpk)  :: factor1
   real(fpk)  :: factor2

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, v
   real(fpk)  :: Mu0v, Mu1v, cumtau, Attenuationsn, Attenuationsn1
   real(fpk)  :: sumd, multiplier, lostau, lostrans_dnn
   real(fpk)  :: CUMSOURCE_DN_OLD, CUMSOURCE_DN_NEW

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS

!  Start Geometry loop

   do v = 1, ngeoms

!  Mu0(v), Mu1(v)

      Mu0v = Mu0(v)
      Mu1v = Mu1(v)

!  Attenuations and Solar solutions
!  ================================

!  Initialize

      nstart = nlayers

!  Attenuations to End points (including TOA)
!  MUST go all the way to NLAYERS (surface term required)

      Attenuations(0) = one
      cumtau = zero
      do n = 1, nlayers  
         cumtau = cumtau+deltaus(n)
         sumd = cumtau / Mu0v   
         if (sumd .lt. cutoff ) then
            Attenuations(n) = exp( - sumd )
         else
            Attenuations(n) = zero
         endif
      enddo

!  Layer integrated Solar sources
!  ==============================

!  Plane/Parallel
!  Special treatment for the horizonal case --> Factor2 = 0, lostrans = 0

      do n = nlayers, nut+1, -1
         if ( Mu1v .gt. zero ) then
            lostau = deltaus(n) / Mu1v
            if ( lostau .lt. cutoff ) then
               lostrans_dn(n) = exp( - lostau )
            else
               lostrans_dn(n) = zero
            endif
         else
            lostrans_dn(n) = zero
         endif
         sources_dn(n) = zero
      enddo
      Attenuationsn = Attenuations(nut)
      do n = nut, 1, -1 
         Attenuationsn1 = Attenuations(n-1)
         if ( Mu1v .gt. zero ) then
            lostau = deltaus(n) / Mu1v
            if ( lostau .lt. cutoff ) then
               lostrans_dnn = exp( - lostau )
            else
               lostrans_dnn = zero
            endif
            lostrans_dn(n) = lostrans_dnn
            factor1 = Attenuationsn1*lostrans_dnn - Attenuationsn
            factor2 = Mu1v / Mu0v - one
            multiplier = factor1 / factor2
            sources_dn(n) = exactscat_dn(n,v) * multiplier
         else
            lostrans_dn(n) = zero
            factor1 = Attenuations(n)
            sources_dn(n) = exactscat_dn(n,v) * factor1
         endif
         Attenuationsn = Attenuationsn1
      enddo

!  Source function integration
!  ===========================

!  Start recursion

      CUMSOURCE_DN_OLD = zero
      NSTART = 1
      NUT_PREV = NSTART - 1

!  Main loop over all output optical depths
!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working Down from NSTART to NUT
!  Check for updating the recursion

      DO UTA = 1, N_USER_LEVELS
         NUT    = USER_LEVELS(UTA)
         DO N = NSTART, NUT
            CUMSOURCE_DN_NEW = SOURCES_DN(N) + LOSTRANS_DN(N) * CUMSOURCE_DN_OLD
            CUMSOURCE_DN_OLD = CUMSOURCE_DN_NEW
         ENDDO
         INTENSITY_DN(UTA,V) = FLUX * CUMSOURCE_DN_OLD
         IF ( NUT .NE. NUT_PREV ) NSTART = NUT + 1
         NUT_PREV = NUT
      ENDDO

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SS_Integral_I_DN_PP

!

subroutine SS_Integral_I_UPDN_PP &
   ( maxgeoms, maxlayers, max_user_levels,              & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling,                        & ! Inputs (Flags)
     ngeoms, nlayers, n_user_levels, user_levels,       & ! Inputs (control output)
     reflec, deltaus, exactscat_up, exactscat_dn, flux, & ! Inputs (Optical)
     Mu0, Mu1,                                          & ! Inputs (Geometry)
     intensity_up, intensity_db, intensity_dn )           ! Outputs

!  Stand-alone routine for Upwelling and Downwelling Solar-beam Single-scatter (SS)
!  Computation of radiance. Inputs: geometry, optical properties

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

!  Flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING

!  Numbers

   integer, Intent(in) ::  NLAYERS
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS,MAXGEOMS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Geometrical inputs
!  ------------------

!  Mu0 = cos(theta_boa), required for PP
!  Mu1 = cos(alpha_boa), required for PP

   real(fpk), Intent(in)  :: Mu0(maxgeoms), Mu1(maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )

   if ( do_upwelling ) then
       call SS_Integral_I_UP_PP &
   ( maxgeoms, maxlayers, max_user_levels,          & ! Inputs (dimensioning)
     ngeoms, nlayers, n_user_levels, user_levels,   & ! Inputs (control output)
     reflec, deltaus, exactscat_up, flux, Mu0, Mu1, & ! Inputs (Optical,Geometry)
     intensity_up, intensity_db )                     ! Outputs
   endif

   if ( do_dnwelling ) then
       call SS_Integral_I_DN_PP &
   ( maxgeoms, maxlayers, max_user_levels,        & ! Inputs (dimensioning)
     ngeoms, nlayers, n_user_levels, user_levels, & ! Inputs (control output)
     deltaus, exactscat_dn, flux, Mu0, Mu1,       & ! Inputs (Optical,Geometry)
     intensity_dn )                                 ! Outputs
   endif

!  Finish

   return
end subroutine SS_Integral_I_UPDN_PP

!  End module

end module FO_ScalarSS_RTCalcs_I_PP_m
