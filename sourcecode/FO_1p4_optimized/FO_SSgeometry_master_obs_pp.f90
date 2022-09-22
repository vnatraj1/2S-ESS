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

module FO_SSGeometry_Master_Obs_PP_m

!  Stand alone geometry for solar scattering only

!   Plane-parallel option, September 2012               (Version 1.2)

   use FO_geometry_Routines_Obs_PP_m

private
public FO_SSGeometry_Master_Obs_PP

contains

subroutine FO_SSGeometry_Master_Obs_PP &
       ( maxgeoms,           & ! Input dimensions
         ngeoms, dtr, vsign, & ! Input control and constants
         obsgeom_boa,        & ! Input
         Mu0, Mu1, cosscat )   ! Output(scat)

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms

!  Geometry control

   integer  , intent(in)     :: ngeoms

!  dtr = degrees-to-Radians. VSIGN = +1 (Up); -1(Down)

   real(ffp), intent(in)     :: dtr, vsign

!  Input angles (Degrees). Enough information for Lattice or Obsgeom.
!  Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)

   real(ffp), intent(In)  :: Obsgeom_boa(maxgeoms,3)

!  Output arguments
!  ================

!  Cosine scattering angle, other cosines

   real(ffp), Intent(out)  :: cosscat (maxgeoms)
   real(ffp), intent(Out)  :: Mu0     (maxgeoms)
   real(ffp), intent(Out)  :: Mu1     (maxgeoms)

!  Plane-parallel, One routine only
!  --------------------------------

   CALL Obsgeom_PlanPar &
      ( maxgeoms, dtr, vsign, & ! Inputs
        ngeoms, obsgeom_boa,  & ! Inputs
        Mu0, Mu1, cosscat )     ! Outputs

!  Finish

   return
end subroutine FO_SSGeometry_Master_Obs_PP

!  Finish

end module FO_SSGeometry_Master_Obs_pp_m
