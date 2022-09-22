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
! #  This Version :   1.4 F90                               #
! #  Release Date :   August 2013                           #
! #                                                         #
! #   Version 1.1,  13 February 2012, First Code            #
! #   Version 1.2,  01 June     2012, Modularization        #
! #   Version 1.3,  29 October  2012, Observation geometry  #
! #   Version 1.4,  31 July     2013, Lattice geometries    #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_geometry_Routines_Obs_EPS_m

use FO_geometry_Generic_PS_m

!  Following routines are for Enhanced PS Obs case
!  -----------------------------------------------

!  LOS-Outgoing: Apply equally to Thermal and Solar-scatter Cases

!  subroutine LosOut_EnhancedPS_Initial
!  subroutine LosOut_EnhancedPS_Quadrature
!  subroutine LosOut_EnhancedPS_QUpgrade
!  subroutine LosOut_EnhancedPS_FindCrit

!  Solar-incoming routines (solar scatter only)

!  subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit
!  subroutine SolarIn_EnhancedPS_Obsgeom_SunPaths

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains

subroutine LosOut_EnhancedPS_Initial  &
          ( maxgeoms, maxlayers, dtr, ngeoms, nlayers, & ! Inputs
            heights, eradius, alpha_boa,               & ! Inputs
            doNadir, radii, Raycon, Lospaths,          & ! Outputs
            alpha, sina, cosa, cota )                    ! Outputs

!  Completely stand-alone geometry routine for the outgoing STD correction
!  This is applicable to Both path geometries (up and down)
!  No Partial layer stuff here

!  This routine: Initial LOS path setup

!  starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                      - height grid, earth radius

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer  , intent(in)   :: maxlayers, maxgeoms
   real(ffp), intent(in)   :: dtr

!  number of geometries

   integer  , intent(in)   :: ngeoms

!  Layer control

   integer  , intent(in)   :: nlayers
   real(ffp), intent(in)   :: eradius, heights (0:maxlayers)

!  input angles

   real(ffp), intent(in)   :: alpha_boa(maxgeoms)

!  Flag for the Nadir case

   logical  , intent(out)  :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant, Lospaths

   real(ffp), intent(out)  :: radii    (0:maxlayers)
   real(ffp), intent(out)  :: Raycon   (maxgeoms)
   real(ffp), intent(out)  :: Lospaths (maxlayers,maxgeoms)

   real(ffp), intent(out)  :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: sina     (0:maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cota     (0:maxlayers,maxgeoms)

!  Local
!  -----

   integer      :: n, n1, v
   real(ffp)    :: radiin, radiin1, Rayconv
   real(ffp)    :: salpha_boa, alpha_boa_R
   real(ffp)    :: calpha, calpha1
   real(ffp)    :: sinanv, alphanv

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Radii

   do n = 0, nlayers
     radii(n) = eradius + heights(n)
   enddo

!  START LOOP
!  ==========

   do v = 1, ngeoms

      doNadir(v) = .false.

!  Special case. Direct nadir viewing. Compute everything and Exit.

      if ( alpha_boa(v) .eq. zero ) then

         doNadir(v) = .true.
         radiin = radii(nlayers)
         do n = nlayers,1,-1
            radiin1 = radii(n-1)
            Lospaths(n,v) = radiin1 - radiin
            radiin = radiin1
         enddo

      else

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  start at BOA

         alpha_boa_R    = alpha_boa(v) * DTR

         if ( alpha_boa(v) .eq. 90.0_ffp ) then
            salpha_boa     = one
            calpha1        = zero
         else
            salpha_boa     = sin(alpha_boa_R)
            calpha1        = cos(alpha_boa_R)
         endif
         cosa(nlayers,v)  = calpha1
         sina(nlayers,v)  = salpha_boa
         cota(nlayers,v)  = calpha1 / salpha_boa

         alpha(nlayers,v) = alpha_boa_R

!  Ray constant

         Rayconv = salpha_boa * radii(nlayers)
         Raycon(v) = Rayconv

!  Whole layer values

         radiin1 = radii(nlayers)
         do n = nlayers - 1, 0, -1
            n1 = n + 1
            radiin = radii(n)
            sinanv = Rayconv / radiin
            alphanv = asin(sinanv) ;  alpha(n,v) = alphanv
            sina(n,v) = sinanv
            calpha  = cos(alphanv) ; cosa(n,v) = calpha
            cota(n,v) = calpha / sinanv
            Lospaths(n1,v) = radiin*calpha - radiin1*calpha1
            radiin1 = radiin
            calpha1 = calpha
         enddo

      endif

!  End loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Initial

!

subroutine LosOut_EnhancedPS_Quadrature  &
       ( maxgeoms, maxlayers, maxfine,     & ! Inputs
         ngeoms, nlayers, nfine,           & ! Inputs
         doNadir, radii, alpha, Raycon,    & ! Inputs
         nfinedivs, radiifine, alphafine,  & ! Outputs
         xfine, wfine, csqfine, cotfine )    ! Outputs

!  Completely stand-alone geometry routine for the outgoing STD correction
!  This is applicable to Both path geometries (up and down)
!  No Partial layer stuff here

!  Extension to Multiple Geometries, 29 October 2012
!  Extension to Lattice  Geometries, 31 July    2013 
!  Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

!  starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                      - height grid, earth radius, Layer control

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)       :: maxlayers, maxfine, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in)       :: nlayers, ngeoms, nfine

!  Flag for the Nadir case

   logical  , intent(in)     :: doNadir(maxgeoms)
  
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon      (maxgeoms)

!  Outputs
!  =======

!  Finelayer divisions is output here

   integer  , intent(out)  :: nfinedivs (maxlayers,maxgeoms)

!  Fine layering

   real(ffp), intent(out)  :: alphafine (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)  :: radiifine (maxfine,maxlayers,maxgeoms)

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxfine,maxlayers,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: csqfine  (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cotfine  (maxfine,maxlayers,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, v
   real(ffp)          :: difh, csfine, tfinej, rayconv, rf, alphanv, alphan1v
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Start geometry loop
!  ===================

   do v = 1, ngeoms

      if (.not. doNadir(v) ) rayconv = raycon(v)

!  Special case. Direct nadir viewing
!  ==================================

!  Compute everything and Exit. Quadratures are height-oriented
!  (This should be the same as the regular pseudo-spherical )

      if ( doNadir(v) ) then

         do n = nlayers, 1, -1
            difh  = radii(n-1) - radii(n)
            nfinedivs(n,v) = nfine
            call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
            do j = 1, nfine
               tfinej = tfine(j)
               radiifine(j,n,v) = radii(n-1) - tfinej
               xfine(j,n,v) = tfinej
               wfine(j,n,v) = afine(j)
            enddo
         enddo

      else

!  Outgoing sphericity geometry (General case)
!  ===========================================

!  Whole layer values

         alphanv = alpha(nlayers,v)
         do n = nlayers, 1, -1
            n1 = n - 1
            alphan1v = alpha(n1,v)
            nfinedivs(n,v) = nfine
            call gauleg_ng (alphan1v,alphanv,tfine,afine,nfine,maxfine)
            do j = 1,  nfine
               tfinej = tfine(j)
               csfine = one / sin(tfinej)
               rf = rayconv * csfine
               radiifine(j,n,v) = rf
               alphafine(j,n,v) = tfinej
               xfine(j,n,v)   = radii(n1) - rf
               wfine(j,n,v)   = afine(j)
               cotfine(j,n,v) = cos(tfinej) * csfine
               csqfine(j,n,v) = csfine * csfine
            enddo
            alphanv = alphan1v
         enddo

      endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_Quadrature

!

SUBROUTINE LosOut_EnhancedPS_FindCrit &
       ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
         extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Inputs, Input/Output
         Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs

!  Purpose: Given a list of Maximum extinctions and LOS angles
!  Then find Critical layers (NCrit) and point where LOS attenuation wipe-outs (Acrit) are achieved
!  Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points

!  Extension to Multiple Geometries, 29 October 2012
!  Equally valid for the Lattice case, if we understand ngeoms = nvzas in this case.

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer  , intent(in)  :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer  , intent(in)  :: nlayers, ngeoms

!  Attenuation and other parameters

   real(ffp), intent(in)  :: Acrit, Cutoff

!  Special case, Nadir viewing

   logical  , intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries

   real(ffp), intent(in)  :: sina (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: cosa (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)

!  Extinctions

   real(ffp), intent(in)  :: Lospaths(maxlayers,maxgeoms)
   real(ffp), intent(in)  :: extinc(maxlayers)

!  Modified inputs
!  ---------------

!  Number of Fine divisions

   integer, intent(inout) :: nfinedivs(maxlayers,maxgeoms)

!  outputs
!  -------

!  Critical layer, Number of Fine divisions for this layer

   integer  , intent(out)  :: Ncrit(maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(out)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

   real(ffp), parameter :: zero = 0.0_ffp
   real(ffp), parameter :: one  = 1.0_ffp

!  Other variables

   logical    ::  trawl
   integer    ::  N, N1, ntrans, v, nfinedivsnv
   real(ffp)  ::  opdep, cumtrans_0, cumtrans_1, trans, transcrit, dcrit, radiin1, TanCrit

!  Initialize

   Ncrit     = 0
   RadCrit   = zero

   fail    = .false.
   message = ' '

!  Start Geometry loop

   do v = 1, ngeoms

      if ( .not. doNadir(v) ) then
         AlphaCrit(v) = zero
         CotCrit(v) = zero
      endif

!  Set trawl
!  Tested March 17th

      trawl = .true. ; n = 0 ; cumtrans_0 = one
      do while (trawl .and. n .lt. nlayers) 
         n = n + 1
         opdep = Lospaths(n,v) * extinc(n)
         if ( opdep .lt. cutoff ) then
            trans = exp ( - opdep )
         else
            trans = zero
         endif
         cumtrans_1 = cumtrans_0 * trans
         nfinedivsnv = nfinedivs(n,v)
         if ( cumtrans_1 .gt. Acrit ) then
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivsnv,ntrans)
            cumtrans_0 = cumtrans_1
         else
            NCrit(v) = n ; trawl = .false.
            transcrit    = Acrit / cumtrans_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivsnv,ntrans)
            dcrit        = - log(transcrit) / extinc(n)
            if ( doNadir(v) ) then
               Radcrit(v)    = radii(n-1) - dcrit      
            else
               n1 = n-1 ; radiin1 = radii(n1)
               TanCrit = radiin1*sina(n1,v)/(radiin1*cosa(n1,v)-dcrit)
               Cotcrit(v)    = one / Tancrit
               AlphaCrit(v)  = atan( TanCrit)
               RadCrit(v)    = sina(n,v) * radii(n) / sin(alphacrit(v))    
            endif
         endif
      enddo

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_FindCrit

!

subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit &
       ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
         cutoff, alpha, radii, extinc, Raycon, theta_boa,           & ! Inputs
         doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Inputs/Outputs
         fail, message )                                              ! Outputs

!  Purpose: Given a list of Maximum extinctions and solar angles at BOA
!           Then find Critical layers (NCrit) and points where TOA attenuation wipe-outs (Acrit) are achieved
!           Then find the LOS angles and Radii (AlphaCrit,RadCrit) for these Critical Points
!           Nadir case, Alpha = 0.0, find only the radius (RadCrit)

!  Extension to Multiple Geometries, 29 October 2012

!  Find the Critical Radius (or angle) in layer Ncrit_i, Use Bisection 
!  based on the Function F(x) = opdep(x) - Crit_opdep

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxlayers, maxgeoms

!  Layer and geometry numbers control

   integer, intent(in) :: nlayers, ngeoms

!  Special case, Nadir viewing

   logical, intent(in)  :: doNadir(maxgeoms)

!  View angles and Radii at layer boundaries, Ray constant

   real(ffp), intent(in)  :: alpha(0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii(0:maxlayers)
   real(ffp), intent(in)  :: Raycon(maxgeoms)

!  Extinctions

   real(ffp), intent(in)  :: extinc(maxlayers)

!  Solar control and other parameters

   real(ffp), intent(in)  :: Acrit, theta_boa(maxgeoms), dtr, cutoff

!  Modified inputs (outputs)
!  -------------------------

!  Overall control (May be switched off if Critical test is negative for all Geometries)

   logical  , intent(inout)  :: DoCrit

!  Critical layer

   integer  , intent(inout)  :: Ncrit(maxgeoms)

!  Number of Fine divisions. This is updated according to Criticality

   integer  , intent(inout)  :: nfinedivs(maxlayers,maxgeoms)

!  Critical angle and radius and cotangent

   real(ffp), intent(inout)  :: AlphaCrit(maxgeoms)
   real(ffp), intent(inout)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Outputs
!  -------

!  Exception handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local variables
!  ---------------

!  Bisection accuracy (lower for the Nadir case, as using distances)

   real(ffp), parameter  :: BisectionAccuracy_Nadir   = 1.0d-6
   real(ffp), parameter  :: BisectionAccuracy_General = 1.0d-9
   integer  , parameter  :: jmax = 50

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Other variables

   character*2::  c2
   logical    ::  Finding, trawl, do_ZeroSunBOA, doCrit_local, doNadirv
   integer    ::  j, n, ntrans, ncrit_i, k, v, cc
   integer    ::  Ncritv, nfinedivsnv
   real(ffp)  ::  theta_boav, radiiNi1, Rayconv, radiik1
   real(ffp)  ::  s0, x1, x2, xmid, rtbis, dx, f, fmid, suncon, radii_x, sx, accuracy
   real(ffp)  ::  dist, theta_0, theta_1, stheta_0, stheta_1, ground_R, sunpaths(maxlayers)
   real(ffp)  ::  opdep, trans, atten_0, atten_1, transcrit, theta_n, stheta_n, theta_boa_R

!  Initialize

   fail    = .false.
   message = ' '

!  Crit Count, initialize

   CC = 0

!  Start Geometry loop

   do v = 1, ngeoms

!  Initial setups

      doCrit_local = .true.
      atten_0 = one ; theta_boa_R = theta_boa(v) * dtr
      doNadirv = doNadir(v)
      Ncritv = Ncrit(v)
      theta_boav = theta_boa(v)
      if ( doNadirv ) then
         s0 = sin(theta_boa_R) ; ground_R = zero ; accuracy = BisectionAccuracy_Nadir
      else
         s0 = zero ; ground_R = alpha(nlayers,v) + theta_boa_R ; accuracy = BisectionAccuracy_General
      endif

!  Trawl through layers until Critical layer is reached. Nfinedivs is updated.
!  Only go down to Initial (LOS-path) Critical layer 
!  Condition on ZerosunBOA changed 27 March 2012

      NCrit_i = 0 ; trawl = .true. ; n = 0
      do while ( trawl .and. ( (Ncritv .eq. 0 .and. n .lt. nlayers) .or. (NCritv .ne. 0 .and. n .lt. NCritv) ) )
         n = n + 1
         do_ZeroSunBOA = (n .eq. nlayers .and. theta_boav .eq. zero) .or. (doNadirv .and. theta_boav .eq. zero)
         if ( doNadirv ) then
            theta_n = theta_boa_R ; stheta_n = s0
         else
            theta_n = ground_R - alpha(n,v) ; stheta_n = sin(theta_n)
         endif
         call FindSunPaths_D(do_ZeroSunBOA,maxlayers,radii(n),radii,&
                             theta_n,stheta_n,n,sunpaths)
         opdep = dot_product(extinc(1:n),sunpaths(1:n))
         if ( opdep .lt. cutoff ) then
            atten_1 = exp ( - opdep )
         else
            atten_1 = zero
         endif
         nfinedivsnv = nfinedivs(n,v)
         if ( atten_1 .gt. Acrit ) then
            trans = atten_1 / atten_0
            ntrans = int(-log10(trans) + 1)
            nfinedivs(n,v) = max(nfinedivsnv,ntrans)
            atten_0 = atten_1
         else
            NCrit_i = n ; trawl = .false.
            transcrit    = Acrit / atten_0
            ntrans = int(-log10(transcrit) + 1)
            nfinedivs(n,v) = max(nfinedivsnv,ntrans)
         endif
      enddo

!  Nothing to do if No criticality (previous Critical values are unchanged)

      if ( trawl .and. NCrit_i .eq. 0 ) then

         if ( trawl .and. NCritv .eq. 0 ) DoCrit_local = .false.

!  If there is criticality ...

      else

         radiiNi1 = radii(NCrit_i-1)

         if ( .not. doNadirv ) then
            Rayconv = Raycon(v)
         endif

!  Bisection: set Highest/Lowest value of Function (layer bottom/top). 

         if ( doNadirv ) then
            x1 = zero               ; x2 = radiiNi1 - radii(NCrit_i)
         else
            x1 = alpha(NCrit_i-1,v) ; x2 = alpha(NCrit_i,v) 
         endif
         F     = log(Atten_0) - cutoff
         FMID  = opdep        - cutoff

!  Bisection: Check bracketing, if OK, perform bisection

         IF ( F*FMID .ge. zero ) THEN
            write(c2,'(I2)') v
            fail = .true. ; message = 'Root must be bracketed for bisection, geometry #'//c2 ; return
         ENDIF
         IF ( F .lt. zero ) THEN
            RTBIS = X1 ; DX = X2-X1
         ELSE
            RTBIS = X2 ; DX = X1-X2
         ENDIF

!  Bisection: Iterate to find the answer

         Finding = .true. ; J = 0
         DO WHILE (Finding .and. j .lt. JMAX)
            J = J + 1 ; dx = 0.5_ffp * dx ; XMID = RTBIS + DX
            if ( doNadirv ) then
               theta_0 = theta_boa_R ; stheta_0 = s0
               radii_x = radiiNi1 - xmid 
            else
               theta_0 = ground_R - xmid ; stheta_0 = sin(theta_0)
               radii_x = Rayconv / sin(xmid)
            endif
            suncon = radii_x * stheta_0
            stheta_1 = suncon / radiiNi1 ;  theta_1 = asin(stheta_1)
            dist = radiiNi1 * sin(theta_0-theta_1) / stheta_0
            opdep = dist * extinc(NCrit_i)
            theta_0 = theta_1 ; stheta_1 = stheta_0
            do k = n - 1, 1, -1
               radiik1 = radii(k-1)
               stheta_1 = suncon / radiik1 ; theta_1 = asin(stheta_1)
               dist = radiik1 * sin(theta_0-theta_1) / stheta_0
               opdep = opdep + dist * extinc(k)
               theta_0 = theta_1 ; stheta_0 = stheta_1
            enddo
            fmid = opdep - cutoff
            IF ( FMID .le. zero ) RTBIS = XMID
            IF ( ABS(DX) .lt. accuracy .or. FMID .eq. zero ) Finding = .false.
         ENDDO

!  Exception (too many bisections)

         if ( Finding ) then
            write(c2,'(I2)') v
            fail = .true. ; message = 'Too many Bisections (540); Root not found, geometry #'//c2 ; return
         endif

!  Set final output if successful

         CC = CC + 1
         if ( doNadirv ) then
            RadCrit(v)   =  radiiNi1 - RTBIS
         else
            AlphaCrit(v) = RTBIS ;  SX = sin(RTBIS)
            RadCrit(v)   = Rayconv / SX
            CotCrit(v)   = sqrt(one-SX*SX) / SX
         endif
         NCrit(v)     = NCrit_i

!  End criticality if loop

      endif

!  Zero the rest

      if ( NCrit(v) .ne. 0 ) nfinedivs(NCrit(v)+1:nlayers,v) = 0

!  End geometry loop

   enddo

!  Reset criticality

   doCrit = ( CC .gt. 0 )

!  Finish

end subroutine SolarIn_EnhancedPS_Obsgeom_FindCrit

!

subroutine LosOut_EnhancedPS_QUpgrade &
       ( maxgeoms, maxlayers, maxfine, ngeoms,     & ! Inputs
         doNadir, radii, alpha, Raycon, nfinedivs, & ! Input LOS layer quantities
         Ncrit, AlphaCrit, RadCrit,                & ! Input Criticality variables
         radiifine, alphafine, xfine,              & ! Inputs/Outputs, Adjusted fine-layering
         wfine, csqfine, cotfine )                   ! Inputs/Outputs, Adjusted fine-layering

!  Completely stand-alone geometry routine for the outgoing STD correction
!  This is applicable to Both path geometries (up and down)
!  No Partial layer stuff here

!  starting inputs are - BOA values of VZA (alpha_boa), in degrees
!                      - height grid, earth radius, Layer control
!                      - Critical Layer control

!  Regular Quadrature need not be done if LOSPATHS is set

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs
!  ------

!  Dimensions

   integer, intent(in)    :: maxlayers, maxfine, maxgeoms

!  geometry numbers control

   integer, intent(in)    :: ngeoms

!  Flag for the Nadir case

   logical  , intent(in)  :: doNadir(maxgeoms)
 
!  Alphas, Radii, Ray constant

   real(ffp), intent(in)  :: alpha      (0:maxlayers,maxgeoms)
   real(ffp), intent(in)  :: radii      (0:maxlayers)
   real(ffp), intent(in)  :: Raycon     (maxgeoms)

!  Critical stuff

   integer  , intent(in)  :: Ncrit     (maxgeoms)
   real(ffp), intent(in)  :: AlphaCrit (maxgeoms)
   real(ffp), intent(in)  :: RadCrit   (maxgeoms)

!  Fine divisions (these are already criticality-adjusted)

   integer  , intent(in)  :: nfinedivs (maxlayers,maxgeoms)

!  Outputs
!  =======

!  Fine layering, Quadratures and local geometry, Adjusted values

   real(ffp), intent(inout)  :: alphafine (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: radiifine (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: xfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: wfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: csqfine   (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(inout)  :: cotfine   (maxfine,maxlayers,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, nfine, g
   real(ffp)          :: raycong, difh, csfine, radiin1, tfinej
   real(ffp)          :: tfine(maxfine), afine(maxfine)

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Start LOS geometry loop

   do g = 1, ngeoms

      if ( .not. doNadir(g) ) raycong = raycon(g)

!  Special case. Direct nadir viewing
!  ==================================

!  Quadratures are height-oriented
!  (This should be the same as the regular pseudo-spherical )
!  -- Adjust quadrature for the Critical layer

      if ( doNadir(g) ) then
         n = NCrit(g) ; radiin1 = radii(n-1) ; difh = radiin1 - Radcrit(g) ; nfine = nfinedivs(n,g)
         call gauleg_ng (zero,difh,tfine,afine,nfine,maxfine)
         do j = 1, nfine
            tfinej = tfine(j)
            radiifine(j,n,g) = radiin1 - tfinej
            xfine(j,n,g) = tfinej
            wfine(j,n,g) = afine(j)
         enddo

!  Outgoing sphericity geometry (General case)
!  ===========================================

      else
         n = NCrit(g) ; n1 = n - 1 ; nfine = nfinedivs(n,g) ; radiin1 = radii(n1)
         call gauleg_ng (alpha(n1,g),AlphaCrit(g),tfine,afine,nfine,maxfine)
         do j = 1, nfine
            tfinej = tfine(j)
            csfine = one / sin(tfinej)
            radiifine(j,n,g) = raycong * csfine
            alphafine(j,n,g) = tfinej
            xfine(j,n,g)   = radiin1 - radiifine(j,n,g)
            wfine(j,n,g)   = afine(j)
            cotfine(j,n,g) = cos(tfinej) * csfine
            csqfine(j,n,g) = csfine * csfine
         enddo
       endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_QUpgrade

!

subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths &
       ( maxgeoms, maxlayers, maxfine,                         & ! Input dimensions
         vsign, dtr, Pie, ngeoms, nlayers,                     & ! Input Control
         obsgeom_boa, doNadir, radii, alpha,                   & ! Input layer/level quantities
         nfinedivs, radiifine, alphafine,                      & ! Input Finelayer variables
         DoCrit, NCrit, RadCrit, AlphaCrit,                    & ! Input Criticality variables
         sunpathsnl, ntraversenl, sunpathsfine, ntraversefine, & ! Outputs
         Mu0, cosscat )                                          ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!  This is for the incoming Solar Beams
!  This is applicable to Both Upwelling and Downwelling LOS-path geometries
!  No partials, this routine

!  starting inputs are the BOA values of SZA, VZA and PHI
!  need also the height grids, earth radius and control
!  need also the complete values of all VZAs along outgoing paths

!  This routine has the fine gridding treatment

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  inputs

!  Dimensions and constants and flag

   integer  , intent(In)    :: maxgeoms, maxlayers, maxfine
   real(ffp), intent(In)    :: vsign, dtr, pie

!  BOA angles

   integer  , intent(In)    :: ngeoms
   real(ffp), intent(In)    :: obsgeom_boa(maxgeoms,3)

!  Layer quantities

   integer  , intent(In)    :: nlayers
   logical  , intent(In)    :: doNadir  (maxgeoms)
   real(ffp), intent(In)    :: alpha    (0:maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radii    (0:maxlayers)

!  Finelayer quantities

   integer  , intent(In)    :: nfinedivs  (maxlayers,maxgeoms)
   real(ffp), intent(In)    :: alphafine (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(In)    :: radiifine (maxfine,maxlayers,maxgeoms)

!  Criticality quantities

   Logical  , intent(In)    :: DoCrit
   integer  , intent(In)    :: NCrit     (maxgeoms)
   real(ffp), intent(In)    :: AlphaCrit (maxgeoms)
   real(ffp), intent(In)    :: RadCrit   (maxgeoms)

!  OUTPUTS
!  =======

!  main outputs (geometry)

   integer  , intent(Out)  :: ntraversenl(maxgeoms)
   real(ffp), intent(Out)  :: sunpathsnl (maxlayers,maxgeoms)

!  Fine level output (geometry)

   integer  , intent(Out)  :: ntraversefine(maxfine,maxlayers,maxgeoms)
   real(ffp), intent(Out)  :: sunpathsfine (maxlayers,maxfine,maxlayers,maxgeoms)

!  scattering angle cosines and associated angles

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat    (maxgeoms)

!  Local

   logical       :: DirectSun, Do_OverheadSun, Do_ZeroSunBOA, Do_Normal, doNadirv
   integer       :: n, j, k, v, NCritv, nfinedivsnv
   real(ffp)     :: SolarDirection(3), Radstart, term1, term2
   real(ffp)     :: salpha_boa, calpha_boa, phi_boa_R, sphi_boa
   real(ffp)     :: theta_boa_R, stheta_boa, ctheta_boa, cphi_boa
   real(ffp)     :: AlphaCritv, RadCritv, alphanlv
   real(ffp)     :: stheta, CumAngle, diffhts(maxlayers)
   real(ffp)     :: theta_all

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Check
!   real(ffp)  :: sumd, sume, sth1

!  Local arrays associated with fine grid output

   logical         :: DirectSunf(maxfine)
   real(ffp)       :: thetaf(maxfine)
   real(ffp)       :: sthetaf(maxfine)

!  precompute

   do k = 1, nlayers
      diffhts(k) = radii(k-1) - radii(k)
   enddo

!  Start geometry loop

   do v = 1, ngeoms

      doNadirv = doNadir(v)
      NCritv = NCrit(v)

      if ( doCrit .and. n .eq. NCritv ) then
         AlphaCritv = AlphaCrit(v) ; RadCritv = RadCrit(v)
      endif

!  Special case

      Do_OverheadSun = obsgeom_boa(v,1) .eq. zero

!  BOA angles

      if ( obsgeom_boa(v,2) .eq. 90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         salpha_boa  = sin(alpha(nlayers,v))
         calpha_boa  = cos(alpha(nlayers,v))
      endif

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      if ( obsgeom_boa(v,1) .eq. 90.0_ffp ) then
         ctheta_boa     = zero
         stheta_boa     = one
      else
         stheta_boa     = sin(theta_boa_R)
         ctheta_boa     = cos(theta_boa_R)
      endif
      Mu0(v) = ctheta_boa

      phi_boa_R   = obsgeom_boa(v,3) * dtr
      cphi_boa    = cos(phi_boa_R)
      sphi_boa    = sin(phi_boa_R)

!  define Unit solar vector at BOA

      if ( Do_OverheadSun ) then
         SolarDirection = zero
      else
         SolarDirection(1) = - stheta_boa * cphi_boa * vsign
         SolarDirection(2) = - stheta_boa * sphi_boa
         SolarDirection(3) = - ctheta_boa
      endif

!  Cosine of scattering angle at boa

      if ( Do_OverheadSun ) then
         cosscat(v) = - vsign * calpha_boa ; if (calpha_boa .eq. zero) cosscat(v) = calpha_boa
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2 + term1 
      endif

      alphanlv = alpha(nlayers,v)

!  General case: LOS path in spherical geometry
!  ============================================

!  Start loop over all layers

      do n = nlayers, 1, -1

         nfinedivsnv = nfinedivs(n,v)

!  Special cases

        DO_ZeroSunBOA  = Do_OverheadSun .and. ( n .eq. nlayers .or. doNadirv )
        DO_Normal      = .not. doCrit .or. ( doCrit .and. n .le. NCritv )

!  Layer boundary Sun position
!  * Local save of angles, cosines, sines and  illumination flags
!  * Use critical ALPHA and RADIUS if N = NCrit
!  * Use Bottom-of-layer values if N < NCrit (BOA values if illuminated)

        if ( vsign .gt. zero .and. n .eq. nlayers ) then
           if ( doCrit .and. n .eq. NCritv ) then
              CumAngle = alphanlv - AlphaCritv ; Radstart = RadCritv
              call FindSun(doNadirv,Do_OverheadSun,Radstart,SolarDirection,CumAngle,theta_boa_R,&
                           theta_all,stheta,DirectSun)
           else
              Radstart = radii(n)
              theta_all = theta_boa_R ; stheta = stheta_boa ; DirectSun = .true.
           endif
        endif

!  Fine-layer sun positions

        if ( Do_Normal ) then
           do j = 1, nfinedivsnv
              CumAngle = alphanlv - alphafine(j,n,v)
              call FindSun(doNadirv,Do_OverheadSun,radiifine(j,n,v),SolarDirection,CumAngle,theta_boa_R,&
                           thetaf(j),sthetaf(j),DirectSunf(j))
           enddo
        endif

!  Sun paths in layer

        if ( vsign .gt. zero .and. n .eq. nlayers ) then
           if ( DirectSun ) then
              call FindSunPaths_D(Do_ZeroSunBOA,maxlayers,Radstart,radii,&
                                  theta_all,stheta,N,sunpathsnl(:,v))
              ntraversenl(v) = nlayers
           else
              call FindSunPaths_T(maxlayers,Pie,Radstart,radii,theta_all,&
                                  stheta,N,sunpathsnl(:,v),ntraversenl(v))
           endif
        endif
        if ( Do_Normal ) then                                !   @@RTSFix 9/5/12 (Addline)
           do j = 1, nfinedivsnv 
              if ( DirectSunf(j) ) then
                 call FindSunPaths_D &
                     (Do_ZeroSunBOA,maxlayers,radiifine(j,n,v),radii,&
                      thetaf(j),sthetaf(j),N,sunpathsfine(:,j,n,v))
                 ntraversefine(j,n,v) = n
              else
                 call FindSunPaths_T &
                     (maxlayers,Pie,radiifine(j,n,v),radii,thetaf(j),&
                      sthetaf(j),N,sunpathsfine(:,j,n,v),ntraversefine(j,n,v))
              endif
!             if ( n.eq.14 ) write(*,*)j,n,radiifine(j,n,v)-radii(n)
           enddo
        endif

!  debugging

!       if ( n.eq.14) then
!       sumd = SUM(sunpaths(1:ntraverse(n),n,v))
!       sth1 = stheta*RadCrit(v)/radii(0)
!       sume = sin(theta_all(n,v) - asin(sth1))*radii(0)/stheta
!       write(*,*)n,sumd,sume
!       do j = 1, nfinedivs(n)
!         sumd = SUM(sunpathsfine(1:ntraversefine(n,j,v),j,n,v))
!         sth1 = sthetaf(j)*radiifine(j,n,v)/radii(0)
!         sume = sin(thetaf(j) - asin(sth1))*radii(0)/sthetaf(j)
!         write(*,*)j,sumd,sume
!       enddo
!       pause
!       endif

!  End layer loop

      enddo

!  End geometry loop

   enddo

!  Finish

   return
end subroutine SolarIn_EnhancedPS_ObsGeom_SunPaths

!  Finish Module

end module FO_geometry_Routines_Obs_EPS_m
