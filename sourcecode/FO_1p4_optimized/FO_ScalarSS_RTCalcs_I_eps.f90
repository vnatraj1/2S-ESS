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

module FO_ScalarSS_RTCalcs_I_EPS_m

!  For a given wavelength, this routine will calculate First-Order upwelling+downwelling Intensities (I):

!  (1) For the Atmospheric Solar Single-scatter and Surface Direct-Beam (SS) sources.

!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.

!  This will perform Enhanced-PS calculations (incoming solar and outgoing LOS-path sphericity) 

!  This is Versions 1-3, without Partials. Code is stand alone with no dependencies.
!  Version 1a, 01 December 2011, R. Spurr, RT Solutions Inc.
!  Version 1b, 13 February 2012, R. Spurr, RT Solutions Inc.
!  Version 2,  01 June     2012, R. Spurr, RT Solutions Inc.
!  Version 3,  29 October  2012, Extension to multiple geometries
!  Version 4,  31 July     2013, Lattice Multi-geometry

!  For Solar sources, the subroutines are
!  SS_Integral_I_UP_EPS    (Upwelling only)
!  SS_Integral_I_DN_EPS    (Downwelling only)
!  SS_Integral_I_UPDN_EPS  (Upwelling and Downwelling)

!  All subroutines public

public

contains

subroutine SS_Integral_I_UP_EPS &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, doNadir,     & ! Inputs (dimensioning, flag)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,     & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, flux,            & ! Inputs (Optical)
     Mu0, NCrit, xfine, wfine, csqfine, cotfine, Raycon,         & ! Inputs (Geometry)
     cota, sunpathsnl, ntraversenl, sunpathsfine, ntraversefine, & ! Inputs (Geometry)
     intensity_up, intensity_db )                                  ! Outputs

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
   integer, Intent(in) :: maxfine
   integer, Intent(in) :: max_user_levels

!  Flag

   logical, Intent(in) :: DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) :: NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) :: NGEOMS

   integer, Intent(in) :: N_USER_LEVELS
   integer, Intent(in) :: USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS, MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!  Mu0 = cos(theta_boa), required for surface term

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: Mu0(maxgeoms)

!  Solar paths 

   integer  , Intent(in)  :: ntraversenl   (maxgeoms)
   real(fpk), Intent(in)  :: sunpathsnl    (maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine  (maxlayers,maxfine,maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfine,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuationsnl
   real(fpk)  :: attenuationsfine

!  Solutions

   real(fpk)  :: Solutionsfine (maxfine,maxlayers)

!  Source function integration results

   real(fpk)  :: sources_up       ( maxlayers )
   real(fpk)  :: lostrans_up      ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, j, v
   integer    :: nCritv, ntraversenlv, ntraversefinejnv
   logical    :: doNadirv
   real(fpk)  :: rayconv, exactscat_upnv, deltausn
   real(fpk)  :: sumd, sum, tran_1, tran, func, kn, ke, xjkn
   real(fpk)  :: cot_1, cot_2
   real(fpk)  :: cumsource_db, CUMSOURCE_UP_OLD, CUMSOURCE_UP_NEW, LOSTRANS_UPN

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

   NUT = USER_LEVELS(1) + 1

!  Start Geometry loop

   do v = 1, ngeoms

      doNadirv = doNadir(v)

!  Raycon(v)

      if ( .not. doNadirv ) rayconv = raycon(v)

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      nstart = nlayers ; nCritv = Ncrit(v); if (nCritv .ne. 0) nstart = nCritv

!  BOA Attenuation
!  MUST go all the way to NLAYERS (surface term required)

      ntraversenlv = ntraversenl(v)
      sumd = dot_product(extinction(1:ntraversenlv),sunpathsnl(1:ntraversenlv,v))
      if (sumd .lt. cutoff ) then
         Attenuationsnl = exp( - sumd )
      else
         Attenuationsnl = zero
      endif

!  Enhanced-spherical, fine-layer attenuations

      do n = nut, nstart
         exactscat_upnv = exactscat_up(n,v)
         do j = 1, nfinedivs(n,v)
            ntraversefinejnv = ntraversefine(j,n,v)
            sumd = dot_product(extinction(1:ntraversefinejnv),sunpathsfine(1:ntraversefinejnv,j,n,v))
            if (sumd .lt. cutoff ) then
               Attenuationsfine = exp( - sumd )
            else
               Attenuationsfine = zero
            endif
            Solutionsfine(j,n) = exactscat_upnv * Attenuationsfine
         enddo
      enddo

!  Layer integrated Solar sources
!  ==============================

!  Enhanced PS: special case (nadir viewing)

      if ( doNadirv ) then
         do n = nlayers, 1, -1
            deltausn = deltaus(n)
            if (deltausn .lt. cutoff ) then
               lostrans_up(n)  = exp ( - deltausn )
            else
               lostrans_up(n)  = zero
            endif
         enddo
         do n = nlayers, nstart+1, -1
            sources_up(n) = zero
         enddo
         do n = nstart, nut, -1
            kn = extinction(n)
            sum = zero
            do j = 1, nfinedivs(n,v)
               xjkn = xfine(j,n,v) * kn
               func = solutionsfine(j,n) * exp ( - xjkn )
               sum = sum + func * wfine(j,n,v)
            enddo
            sources_up(n) = sum * kn
         enddo
         do n = nut-1, 1, -1
            sources_up(n) = zero
         enddo
      endif

!  Enhanced PS: General case

      if ( .not. doNadirv ) then
         cot_1 = cota(nlayers,v)
         do n = nlayers, nstart+1, -1
            cot_2 = cota(n-1,v)
            ke = rayconv * extinction(n)
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            sources_up(n) = zero
            cot_1 = cot_2
         enddo
         do n = nstart, nut, -1
            cot_2 = cota(n-1,v)
            ke = rayconv * extinction(n)
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            sum = zero
            do j = 1, nfinedivs(n,v)
               tran = exp ( - ke * ( cot_2 - cotfine(j,n,v) ) )
               func = solutionsfine(j,n) * csqfine(j,n,v) * tran
               sum  = sum + func * wfine(j,n,v)
            enddo
            sources_up(n) = sum * ke
            cot_1 = cot_2
         enddo
         do n = nut-1, 1, -1
            cot_2 = cota(n-1,v)
            ke = rayconv * extinction(n)
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_up(n) = tran_1
            sources_up(n) = zero
            cot_1 = cot_2
         enddo
      endif

!  Source function integration
!  ===========================

!  Start recursion ( For Direct Beam, use PI.mu0.R.Atten )

      CUMSOURCE_UP_OLD = zero
      CUMSOURCE_DB     = 4.0_fpk * Mu0(v) * REFLEC(v) * attenuationsnl
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

!  Main loop over all output optical depths
!  NLEVEL = Layer index for given optical depth
!  Cumulative source terms : Loop over layers working upwards from NSTART to level NUT
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
end subroutine SS_Integral_I_UP_EPS

!

subroutine SS_Integral_I_DN_EPS &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, doNadir,  & ! Inputs (dimensioning,flag)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
     extinction, deltaus, exactscat_dn, flux,                 & ! Inputs (Optical)
     NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     Raycon, radii, cota, sunpathsfine, ntraversefine,        & ! Inputs (Geometry)
     intensity_dn )                                             ! Outputs

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
   integer, Intent(in) :: maxfine
   integer, Intent(in) :: max_user_levels

!  Flag

   logical, Intent(in) ::  doNadir(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux

   real(fpk), Intent(in) ::  FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)

!  Solar paths 

   integer  , Intent(in)  :: ntraversefine (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine  (maxlayers,maxfine,maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfine,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )

!  LOCAL
!  -----

!  Attenuations

   real(fpk)  :: attenuationsfine

!  Solutions

   real(fpk)  :: Solutionsfine (maxfine,maxlayers)

!  Source function integration results

   real(fpk)  :: sources_dn       ( maxlayers )
   real(fpk)  :: lostrans_dn      ( maxlayers )

!  Help

   integer    :: n, uta, nstart, nut, nut_prev, j, v
   integer    :: nCritv, ntraversefinejnv
   logical    :: doNadirv
   real(fpk)  :: rayconv, exactscat_dnnv, deltausn
   real(fpk)  :: RadCritv, CotCritv
   real(fpk)  :: radiin, radiin1
   real(fpk)  :: sumd, sum, tran_1, tran, func, kn, ke, xjkn, trand
   real(fpk)  :: cot_1, cot_2, rdiff, cot_c
   real(fpk)  :: CUMSOURCE_DN_OLD, CUMSOURCE_DN_NEW

   real(fpk), parameter  :: cutoff = 88.0_fpk
   real(fpk), parameter  :: zero   = 0.0_fpk
   real(fpk), parameter  :: one    = 1.0_fpk

!  Bookkeeping

      NUT = USER_LEVELS(N_USER_LEVELS) + 1
      IF ( NUT > NLAYERS ) NUT = NLAYERS

!  Start Geometry loop

   do v = 1, ngeoms

      doNadirv = doNadir(v)

!  Raycon(v)

      if ( .not. doNadir(v) ) rayconv = raycon(v)

!  Attenuations and Solar solutions
!  ================================

!  Initialize, only to layer Ncrit if applicable

      nstart = nlayers ; nCritv = Ncrit(v); if (nCritv .ne. 0) nstart = nCritv
      if ( doNadirv ) RadCritv = RadCrit(v)
      if ( .not. doNadirv ) CotCritv = CotCrit(v)

!  Enhanced-spherical, fine-layer attenuations

      do n = 1, min(nut,nstart)
         exactscat_dnnv = exactscat_dn(n,v)
         do j = 1, nfinedivs(n,v)
            ntraversefinejnv = ntraversefine(j,n,v)
            sumd = dot_product(extinction(1:ntraversefinejnv),sunpathsfine(1:ntraversefinejnv,j,n,v))
            if (sumd .lt. cutoff ) then
               Attenuationsfine = exp( - sumd )
            else
               Attenuationsfine = zero
            endif
            Solutionsfine(j,n) = exactscat_dnnv * Attenuationsfine
         enddo
      enddo

!  Layer integrated Solar sources
!  ==============================

!  Enhanced PS: special case (nadir viewing)

      if ( doNadirv ) then
         do n = nlayers, 1, -1
            deltausn = deltaus(n)
            if (deltausn .lt. cutoff ) then 
               lostrans_dn(n)  = exp ( - deltausn )
            else
               lostrans_dn(n)  = zero
            endif
         enddo
         radiin = radii(min(nut,nstart))
         do n = min(nut,nstart), 1, -1
            kn = extinction(n)
            radiin1 = radii(n-1)
            rdiff = radiin1 - radiin ; if ( n .eq. NCritv) rdiff = radiin1 - RadCritv      
            trand = one ; if ( n .eq. NCritv) trand = exp ( -kn * ( RadCritv -radiin ) )
            radiin = radiin1
            sum = zero
            do j = 1, nfinedivs(n,v)
               xjkn = ( rdiff - xfine(j,n,v) ) * kn
               func = solutionsfine(j,n) * exp ( - xjkn )
               sum = sum + func * wfine(j,n,v)
            enddo
            sources_dn(n) = sum * kn
!  @@ Robfix, add following line
            if ( n .eq. NCritv ) sources_dn(n) = sources_dn(n) * trand
!  @@ End Robfix add line
         enddo
         do n = min(nut,nstart)+1, nlayers
            sources_dn(n) = zero
         enddo
      endif

!  Enhanced PS: General case

      if ( .not. doNadirv ) then
         cot_1 = cota(nlayers,v)
         do n = nlayers, min(nut,nstart)+1, -1
            cot_2 = cota(n-1,v)
            cot_c = cot_1  ; if ( n .eq. NCritv ) cot_c = CotCritv
            ke = rayconv * extinction(n)
            trand = one  ; if ( n .eq. NCritv ) trand = exp ( - ke * ( CotCritv - cot_1 ) )
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            cot_1 = cot_2
            sources_dn(n) = zero
         enddo
         do n = min(nut,nstart), 1, -1
            cot_2 = cota(n-1,v)
            cot_c = cot_1  ; if ( n .eq. NCritv ) cot_c = CotCritv
            ke = rayconv * extinction(n)
            trand = one  ; if ( n .eq. NCritv ) trand = exp ( - ke * ( CotCritv - cot_1 ) )
            tran_1 = exp ( - ke * ( cot_2 - cot_1 ) )
            lostrans_dn(n) = tran_1
            sum = zero
            do j = 1, nfinedivs(n,v)
               tran = exp ( - ke * ( cotfine(j,n,v) - cot_c ) )   !  Down
               func = solutionsfine(j,n) * csqfine(j,n,v) * tran
               sum  = sum + func * wfine(j,n,v)
            enddo
            sources_dn(n) = sum * ke * trand
            cot_1 = cot_2
         enddo
      endif

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
end subroutine SS_Integral_I_DN_EPS

!

subroutine SS_Integral_I_UPDN_EPS   &
   ( maxgeoms, maxlayers, maxfine, max_user_levels,                      & ! Inputs (dimensioning)
     do_upwelling, do_dnwelling, doNadir,                                & ! Inputs (Flags)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,             & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, exactscat_dn, flux,      & ! Inputs (Optical)
     Mu0, NCrit, RadCrit, CotCrit,                                       & ! Inputs (Geometry)
     xfine, wfine, csqfine, cotfine, Raycon, radii, cota,                & ! Inputs (Geometry)
     sunpaths_up_nl, ntraverse_up_nl, sunpathsfine_up, ntraversefine_up, & ! Inputs (Geometry)
     sunpathsfine_dn, ntraversefine_dn,                                  & ! Inputs (Geometry)
     intensity_up, intensity_db, intensity_dn )                            ! Outputs

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
   integer, Intent(in) :: maxfine
   integer, Intent(in) :: max_user_levels

!  Flags

   logical, Intent(in) ::  DO_UPWELLING
   logical, Intent(in) ::  DO_DNWELLING

   logical, Intent(in) ::  DONADIR(MAXGEOMS)

!  Numbers

   integer, Intent(in) ::  NLAYERS, NFINEDIVS(MAXLAYERS,MAXGEOMS)
   integer, Intent(in) ::  NGEOMS

   integer, Intent(in) ::  N_USER_LEVELS
   integer, Intent(in) ::  USER_LEVELS ( MAX_USER_LEVELS )

!  Optical inputs
!  --------------

!  Atmosphere

   real(fpk), Intent(in) :: EXTINCTION  ( MAXLAYERS )
   real(fpk), Intent(in) :: DELTAUS     ( MAXLAYERS )
   real(fpk), Intent(in) :: EXACTSCAT_UP( MAXLAYERS,MAXGEOMS )
   real(fpk), Intent(in) :: EXACTSCAT_DN( MAXLAYERS,MAXGEOMS )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(fpk), Intent(in) :: REFLEC(MAXGEOMS), FLUX

!  Geometrical inputs
!  ------------------

!  Ray constant, Cotangents, Critical layer
!  Mu0 = cos(theta_boa), required for surface term

   integer  , Intent(in)  :: NCrit(maxgeoms)
   real(fpk), Intent(in)  :: Raycon(maxgeoms), cota(0:maxlayers,maxgeoms), radii(0:maxlayers)
   real(fpk), Intent(in)  :: Mu0(maxgeoms), RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Solar paths 

   integer  , Intent(in)  :: ntraverse_up_nl  (maxgeoms)
   real(fpk), Intent(in)  :: sunpaths_up_nl   (maxlayers,maxgeoms)
   integer  , Intent(in)  :: ntraversefine_up (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine_up  (maxlayers,maxfine,maxlayers,maxgeoms)

   integer  , Intent(in)  :: ntraversefine_dn (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: sunpathsfine_dn  (maxlayers,maxfine,maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS

   real(fpk), Intent(in)  :: xfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: wfine   (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: csqfine (maxfine,maxlayers,maxgeoms)
   real(fpk), Intent(in)  :: cotfine (maxfine,maxlayers,maxgeoms)

!  Outputs
!  -------

   real(fpk), Intent(Out)  :: intensity_up     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_db     ( max_user_levels,maxgeoms )
   real(fpk), Intent(Out)  :: intensity_dn     ( max_user_levels,maxgeoms )

   if ( do_upwelling ) then
       call SS_Integral_I_UP_EPS &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, doNadir,             & ! Inputs (dimensioning, flag)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,             & ! Inputs (control output)
     reflec, extinction, deltaus, exactscat_up, flux,                    & ! Inputs (Optical)
     Mu0, NCrit, xfine, wfine, csqfine, cotfine, Raycon, cota,           & ! Inputs (Geometry)
     sunpaths_up_nl, ntraverse_up_nl, sunpathsfine_up, ntraversefine_up, & ! Inputs (Geometry)
     intensity_up, intensity_db )                                          ! Outputs
   endif

   if ( do_dnwelling ) then
       call SS_Integral_I_DN_EPS &
   ( maxgeoms, maxlayers, maxfine, max_user_levels, doNadir,  & ! Inputs (dimensioning)
     ngeoms, nlayers, nfinedivs, n_user_levels, user_levels,  & ! Inputs (control output)
     extinction, deltaus, exactscat_dn, flux,                 & ! Inputs (Optical)
     NCrit, RadCrit, CotCrit, xfine, wfine, csqfine, cotfine, & ! Inputs (Geometry)
     Raycon, radii, cota, sunpathsfine_dn, ntraversefine_dn,  & ! Inputs (Geometry)
     intensity_dn )                                             ! Outputs
   endif

!  Finish

   return
end subroutine SS_Integral_I_UPDN_EPS

!  End module

end module FO_ScalarSS_RTCalcs_I_EPS_m
