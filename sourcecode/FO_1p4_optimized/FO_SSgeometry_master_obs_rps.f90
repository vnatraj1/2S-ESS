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

module FO_SSGeometry_Master_Obs_RPS_m

!  Stand alone geometry for solar scattering only

   use FO_geometry_Generic_PS_m
   use FO_geometry_Routines_Obs_RPS_m

private
public FO_SSGeometry_Master_Obs_RPS

contains

subroutine FO_SSGeometry_Master_Obs_RPS &
       ( maxgeoms, maxlayers,           & ! Input dimensions
         ngeoms, nlayers, dtr, vsign,   & ! Input control and constants
         eradius, heights, obsgeom_boa, & ! Input
         Mu0, Mu1, cosscat,             & ! Output(scat)
         sunpaths )                       ! Output(Sunpaths)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms, maxlayers

!  Layer and geometry control. Finelayer divisions may be changed

   integer  , intent(in)     :: ngeoms, nlayers

!  dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, vsign

!  Radius + heights

   real(ffp), intent(in)     :: eradius, heights (0:maxlayers)

!  input angles (Degrees). Enough information for Lattice or Obsgeom.
!  Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)

   real(ffp), intent(In)  :: Obsgeom_boa(maxgeoms,3)

!  Output arguments
!  ================

!  Solar paths

   real(ffp), Intent(out)  :: sunpaths      (maxlayers,maxlayers,maxgeoms)

!  Cosine scattering angle, other cosines

   real(ffp), Intent(out)  :: cosscat (maxgeoms)
   real(ffp), intent(Out)  :: Mu0     (maxgeoms)
   real(ffp), intent(Out)  :: Mu1     (maxgeoms)

!  Regular PS, One routine only
!  ----------------------------

   CALL Obsgeom_RegularPS &
      ( maxgeoms, maxlayers, dtr, vsign,                & ! Inputs
        ngeoms, nlayers, obsgeom_boa, heights, eradius, & ! Inputs
        sunpaths, Mu0, Mu1, cosscat )                     ! Outputs
!  Finish

   return
end subroutine FO_SSGeometry_Master_Obs_RPS

!  Finish

end module FO_SSGeometry_Master_Obs_RPS_m
