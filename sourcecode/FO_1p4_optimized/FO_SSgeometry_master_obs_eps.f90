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
! #   Version 1.3,  29 October  2012, Obsgeom Multi-geom.   #
! #   Version 1.4,  31 July     2013, Lattice Multi-geom.   #
! #                                                         #
! ###########################################################

! ##########################################################
! #                                                        #
! #   This Version of FIRST_ORDER comes with a GNU-style   #
! #   license. Please read the license carefully.          #
! #                                                        #
! ##########################################################

module FO_SSGeometry_Master_Obs_EPS_m

!  Stand alone geometry for solar scattering only

   use FO_geometry_Generic_PS_m
   use FO_geometry_Routines_Obs_EPS_m

private
public FO_SSGeometry_Master_Obs_EPS

contains

subroutine FO_SSGeometry_Master_Obs_EPS &
       ( maxgeoms, maxlayers, maxfine,                         & ! Input dimensions
         ngeoms, nlayers, nfine, dtr, Pie, vsign,              & ! Input control and constants
         eradius, heights, obsgeom_boa,                        & ! Input
         doCrit, Acrit, extinc,                                & ! Input/Output, Input
         doNadir, Raycon, radii, cota,                         & ! Output(level)
         nfinedivs, xfine, wfine, csqfine, cotfine,            & ! Output(Fine)
         NCrit, RadCrit, CotCrit, Mu0, cosscat,                & ! Output(Crit/scat)
         sunpathsnl, ntraversenl, sunpathsfine, ntraversefine, & ! Output(Sunpaths)
         fail, message, trace )                                  ! Output(Status)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxlayers, maxfine

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers, nfine

!  dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, Pie, vsign

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees)
!  Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)

   real(ffp), intent(in)     :: Obsgeom_boa(maxgeoms,3)

!  Critical adjustment control for cloud layers

   real(ffp), intent(in)     :: extinc(maxlayers)
   real(ffp), intent(in)     :: Acrit
   logical  , intent(inout)  :: doCrit

!  Output arguments
!  ================

!  Flag for the Nadir case

   logical  , intent(out)    :: doNadir(maxgeoms)

!  Cotangents, Radii, Ray constant
!  WARNING: Adjusted geometry will require maxgeoms dimension

   real(ffp), intent(out)    :: radii    (0:maxlayers)
   real(ffp), intent(out)    :: Raycon   (maxgeoms)
   real(ffp), intent(inout)  :: cota     (0:maxlayers,maxgeoms)

!  LOS Quadratures for Enhanced PS, find layering output
!  WARNING: Adjusted geometry will require maxgeoms dimension

   integer  , intent(out)    :: nfinedivs (maxlayers,maxgeoms)
   real(ffp), intent(out)    :: xfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: wfine     (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: csqfine   (maxfine,maxlayers,maxgeoms)
   real(ffp), intent(out)    :: cotfine   (maxfine,maxlayers,maxgeoms)

!  Critical layer
!  WARNING: Adjusted geometry will require maxgeoms dimension

   integer  , intent(out)  :: Ncrit(maxgeoms)
   real(ffp), intent(out)  :: RadCrit(maxgeoms), CotCrit(maxgeoms)

!  Solar paths

   integer  , Intent(out)  :: ntraversenl     (maxgeoms)
   real(ffp), Intent(out)  :: sunpathsnl      (maxlayers,maxgeoms)

   integer  , Intent(out)  :: ntraversefine (maxfine,maxlayers,maxgeoms)
   real(ffp), Intent(out)  :: sunpathsfine  (maxlayers,maxfine,maxlayers,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp), Intent(out)  :: cosscat (maxgeoms)
   real(ffp), intent(Out)  :: Mu0     (maxgeoms)

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Local arguments
!  ===============

!  LOS path lengths

   real(ffp)  :: Lospaths (maxlayers,maxgeoms)

!  Alphas

   real(ffp)  :: alpha    (0:maxlayers,maxgeoms)

!  Other angles and cos/sin values.

   real(ffp)  :: cosa     (0:maxlayers,maxgeoms)
   real(ffp)  :: sina     (0:maxlayers,maxgeoms)

   real(ffp)  :: alphafine (maxfine,maxlayers,maxgeoms)
   real(ffp)  :: radiifine (maxfine,maxlayers,maxgeoms)

!  Critical values

   real(ffp)  :: AlphaCrit(maxgeoms)

!  Help variables
!  --------------

   real(ffp), parameter   :: zero = 0.0_ffp
   real(ffp)              :: cutoff

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '
   cutoff = zero; if (ACrit .gt. zero) cutoff = -log(ACrit)

!  Enhanced PS; proceed in 4 Steps
!  ===============================

!   doCrit = .false.

!  Step 1; Initial LOS-path quantities, OUTGOING Beam
!  --------------------------------------------------

!  1a.  Given heights and BOA LOS angles, compute path angles and radii
!  1b.  LOS fine-layer quadratures. Non-adjusted, no Criticality
!  1c.  Find Critical-layer adjustments (Optional)

   CALL LosOut_EnhancedPS_Initial &
      ( maxgeoms, maxlayers, dtr, ngeoms, nlayers,  & ! Input
        heights, eradius, obsgeom_boa(:,2),         & ! Input
        doNadir, radii, Raycon, Lospaths,           & ! Output
        alpha, sina, cosa, cota )                     ! Output

   CALL LosOut_EnhancedPS_Quadrature &
      ( maxgeoms, maxlayers, maxfine,    & ! Input 
        ngeoms, nlayers, nfine,          & ! Input
        doNadir, radii, alpha, Raycon,   & ! Input
        nfinedivs, radiifine, alphafine, & ! Output
        xfine, wfine, csqfine, cotfine )   ! Output

   if ( doCrit) then
      CALL LosOut_EnhancedPS_FindCrit &
         ( maxgeoms, maxlayers, ngeoms, nlayers, Acrit, Cutoff, doNadir, & ! Inputs
           extinc, Lospaths, sina, cosa, radii, nfinedivs,               & ! Input
           Ncrit, AlphaCrit, RadCrit, CotCrit, fail, message )             ! Outputs
      if ( Fail ) then
         trace = 'Error from LosOut_EnhancedPS_FindCrit in SSGeometry_Master_Obs_EPS' ; return
      endif
   endif

!  Step 2, INCOMING SOLAR BEAMS.
!  -----------------------------
   
!  2a.  Incoming, find Critical-layer adjustments        (Optional)
!  2b.  Critical-layer Upgrade of Quadrature done here   (Optional)
!  2c.  Incoming, solar pathlengths

!  Step 2a
!  =======

   if ( doCrit) then
      call SolarIn_EnhancedPS_Obsgeom_FindCrit &
         ( maxgeoms, maxlayers, ngeoms, nlayers, doNadir, dtr, Acrit, & ! Input
           cutoff, alpha, radii, extinc, Raycon, obsgeom_boa(:,1),    & ! Inputs
           doCrit, Ncrit, nfinedivs, AlphaCrit, RadCrit, CotCrit,     & ! Outputs
           fail, message )                                              ! Outputs
      if ( Fail ) then
         trace = 'Error from SolarIn_EnhancedPS_Obsgeom_FindCrit in SSGeometry_Master_Obs_EPS' ; return
      endif
   endif

!  Step 2b
!  =======

   if ( doCrit) then
      call LosOut_EnhancedPS_QUpgrade &
         ( maxgeoms, maxlayers, maxfine, ngeoms,     & ! Input
           doNadir, radii, alpha, Raycon, nfinedivs, & ! Input LOS layer quantities
           Ncrit, AlphaCrit, RadCrit,                & ! Input Criticality variables
           radiifine, alphafine, xfine,              & ! Output, Adjusted fine-layering
           wfine, csqfine, cotfine )                   ! Output, Adjusted fine-layering
   endif

!  Step 2c
!  =======

   if (.not. doCrit) NCrit = 0

   CALL SolarIn_EnhancedPS_ObsGeom_SunPaths &
      ( maxgeoms, maxlayers, maxfine,                         & ! Input dimensions
        vsign, dtr, Pie, ngeoms, nlayers,                     & ! Input Control
        obsgeom_boa, doNadir, radii, alpha,                   & ! Input layer/level quantities
        nfinedivs, radiifine, alphafine,                      & ! Input Finelayer variables
        DoCrit, NCrit, RadCrit, AlphaCrit,                    & ! Input Criticality variables
        sunpathsnl, ntraversenl, sunpathsfine, ntraversefine, & ! Output
        Mu0, cosscat )                                          ! Output

!  Finish

   return
end subroutine FO_SSGeometry_Master_Obs_EPS

!  Finish

end module FO_SSGeometry_Master_Obs_EPS_m


