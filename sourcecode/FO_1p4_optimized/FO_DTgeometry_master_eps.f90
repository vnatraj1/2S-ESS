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

module FO_DTGeometry_Master_EPS_m

!  Stand alone geometry for Direct Thermal only

   use FO_geometry_Routines_Obs_EPS_m, only: LosOut_EnhancedPS_Initial, LosOut_EnhancedPS_FindCrit
   use FO_geometry_Routines_Thermal_EPS_m

private
public FO_DTGeometry_Master_EPS

contains

subroutine FO_DTGeometry_Master_EPS &
       ( maxgeoms, maxlayers, maxfine,              & ! Input dimensions
         ngeoms, nlayers, nfine,                    & ! Inputs
         dtr, eradius, heights, alpha_boa,          & ! Inputs
         doCrit, Acrit, extinc,                     & ! Input/Output, Input
         doNadir, Raycon, radii, cota,              & ! Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine, & ! Output(Fine)
         NCrit, RadCrit, CotCrit,                   & ! Output(Critical)
         fail, message, trace )                       ! Output(Status)

   implicit none

!  Parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxlayers, maxfine, maxgeoms

!  Layer and geometry control

   integer  , intent(in)     :: ngeoms, nlayers, nfine

!  dtr = degrees-to-Radians

   real(ffp), intent(in)     :: dtr

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees)

   real(ffp), intent(in)     :: alpha_boa(maxgeoms)

!  Critical adjustment for cloud layers

   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit
   logical  , intent(inout)  :: doCrit

!  Output arguments
!  ================

!  Flag for the Nadir case.

   logical  , intent(out)    :: doNadir(maxgeoms)
  
!  Cotangents, Radii, Ray constant. Intent(inout)

   real(ffp), intent(out)    :: radii    (0:maxlayers)
   real(ffp), intent(out)    :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS, fine layering output

   integer  , intent(out)    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(out)    :: xfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: wfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: csqfine   (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: cotfine   (maxfine,maxlayers,maxgeoms)

!  Critical layer

   integer  , intent(out)    :: Ncrit(maxgeoms)
   real(ffp), intent(out)    :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths  (maxlayers,maxgeoms)

!  Alphas  

   real(ffp)  :: alpha     (0:maxlayers,maxgeoms)

!  Other angles

   real(ffp)  :: cosa      (0:maxlayers,maxgeoms)
   real(ffp)  :: sina      (0:maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables
!  --------------

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp)              :: cutoff

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   cutoff = zero; if (ACrit .gt. zero) cutoff = -log(ACrit)

!  Enhanced PS; proceed in 3 Steps
!  ===============================

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!  1a.  Given heights and BOA LOS angles, compute path angles and radii
!  1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!  1c.  Find Critical-layer adjustments (Optional)

   CALL LosOut_EnhancedPS_Initial &
      ( maxgeoms, maxlayers, dtr, ngeoms, nlayers, & ! Input
        heights, eradius, alpha_boa,               & ! Input
        doNadir, radii, Raycon, Lospaths,          & ! Output
        alpha, sina, cosa, cota )                    ! Output

   CALL LosOut_EnhancedPS_Quadrature_Thermal &
      ( maxgeoms, maxlayers, maxfine,     & ! Input 
        ngeoms, nlayers, nfine,           & ! Input
        doNadir, radii, alpha, Raycon,    & ! Input
        nfinedivs, xfine, wfine, csqfine, cotfine )    ! Output

   if ( doCrit) then
      CALL LosOut_EnhancedPS_FindCrit &
         ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
           extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
           Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
      if ( Fail ) then
         trace = 'Error from LosOut_EnhancedPS_FindCrit in FO_DTGeometry_Master' ; return
      endif
   endif

   if ( doCrit) then
      call LosOut_EnhancedPS_QUpgrade_Thermal &
         ( maxgeoms, maxlayers, maxfine, ngeoms,      & ! Input
           doNadir, radii, alpha, Raycon, nfinedivs,  & ! Input LOS layer quantities
           Ncrit, AlphaCrit, RadCrit,                 & ! Input Criticality variables
           xfine, wfine, csqfine, cotfine )             ! Output, Adjusted fine-layering
   endif
 
!  Finish

   return
end subroutine FO_DTGeometry_Master_EPS

!  Finish

end module FO_DTGeometry_Master_EPS_m
