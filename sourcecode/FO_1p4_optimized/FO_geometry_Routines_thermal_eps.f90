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

module FO_geometry_Routines_Thermal_EPS_m

use FO_geometry_Generic_PS_m, only: gauleg_ng

!  Following routines are for Enhanced PS Thermal case
!  ---------------------------------------------------

!  LOS-Outgoing: Thermal

!  subroutine LosOut_EnhancedPS_Quadrature_Thermal
!  subroutine LosOut_EnhancedPS_QUpgrade_Thermal

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains

subroutine LosOut_EnhancedPS_Quadrature_Thermal  &
       ( maxgeoms, maxlayers, maxfine,     & ! Inputs
         ngeoms, nlayers, nfine,           & ! Inputs
         doNadir, radii, alpha, Raycon,    & ! Inputs
         nfinedivs, xfine, wfine, csqfine, cotfine ) ! Outputs

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

!  Quadratures

   real(ffp), intent(out)  :: xfine    (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)  :: wfine    (maxfine,maxlayers,maxgeoms)

!  Local geoemetry arrays

   real(ffp), intent(out)  :: csqfine  (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)  :: cotfine  (maxfine,maxlayers,maxgeoms)

!  Local
!  -----

   integer            :: n, n1, j, v
   real(ffp)          :: difh, csfine, tfinej, rayconv, alphanv, alphan1v
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
               xfine(j,n,v) = tfine(j)
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
               xfine(j,n,v)   = radii(n1) - rayconv * csfine
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
end subroutine LosOut_EnhancedPS_Quadrature_Thermal

!

subroutine LosOut_EnhancedPS_QUpgrade_Thermal &
       ( maxgeoms, maxlayers, maxfine, ngeoms,     & ! Inputs
         doNadir, radii, alpha, Raycon, nfinedivs, & ! Input LOS layer quantities
         Ncrit, AlphaCrit, RadCrit,                & ! Input Criticality variables
         xfine, wfine, csqfine, cotfine )            ! Inputs/Outputs, Adjusted fine-layering

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
            xfine(j,n,g) = tfine(j)
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
            xfine(j,n,g)   = radiin1 - raycong * csfine
            wfine(j,n,g)   = afine(j)
            cotfine(j,n,g) = cos(tfinej) * csfine
            csqfine(j,n,g) = csfine * csfine
         enddo
       endif

!  End geometry loop

   enddo

!  Finish

   return
end subroutine LosOut_EnhancedPS_QUpgrade_Thermal

!  Finish Module

end module FO_geometry_Routines_Thermal_EPS_m
