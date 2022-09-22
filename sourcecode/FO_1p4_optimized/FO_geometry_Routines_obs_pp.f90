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

module FO_geometry_Routines_Obs_PP_m

!  Plane-parallel routine
!  ----------------------

!    subroutine Obsgeom_PlanPar

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  EVERYTHING PUBLIC HERE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

public

contains

subroutine ObsGeom_PlanPar &
          ( maxgeoms, dtr, vsign, & ! Inputs
            ngeoms, obsgeom_boa,  & ! Inputs
            Mu0, Mu1, cosscat )     ! Outputs

!  Completely stand-alone geometry routine for Accurate SS
!  This is for the Plane-parallel choice
!  This is applicable to the Upwelling and/or/Downwelling LOS-path geometries
!  No partials, this routine

!  Starting inputs are the BOA values of SZA, VZA and PHI

      implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs

   integer  , intent(In)   :: maxgeoms
   real(ffp), intent(In)   :: dtr, vsign

   integer  , intent(In)   :: ngeoms

   real(ffp), intent(In)   :: obsgeom_boa(maxgeoms,3)

!  Main outputs (geometry)

   real(ffp), intent(Out)  :: Mu0        (maxgeoms)
   real(ffp), intent(Out)  :: Mu1        (maxgeoms)
   real(ffp), intent(Out)  :: cosscat(maxgeoms)

!  Local

   logical        :: Do_OverheadSun
   integer        :: v
   real(ffp)      :: alpha_boa_R, theta_boa_R
   real(ffp)      :: salpha_boa, calpha_boa
   real(ffp)      :: stheta_boa, ctheta_boa, cphi_boa
   real(ffp)      :: term1, term2

   real(ffp), parameter  :: zero = 0.0_ffp
   real(ffp), parameter  :: one  = 1.0_ffp

!  Start geometry loop

   do v = 1, ngeoms

!  BOA angles

      alpha_boa_R    = obsgeom_boa(v,2) * DTR
      if ( obsgeom_boa(v,2).eq.90.0_ffp ) then
         calpha_boa     = zero
         salpha_boa     = one
      else
         calpha_boa     = cos(alpha_boa_R)
         salpha_boa     = sin(alpha_boa_R)
      endif
      Mu1(v) = calpha_boa

      theta_boa_R    = obsgeom_boa(v,1) * DTR
      stheta_boa     = sin(theta_boa_R)
      ctheta_boa     = cos(theta_boa_R)
      Mu0(v) = ctheta_boa

      cphi_boa       = cos(obsgeom_boa(v,3) * dtr)

!  Overhead Sun

      Do_OverheadSun = (obsgeom_boa(v,1) .eq. zero)

!  Set Alpha, scattering angle

      if ( Do_OverheadSun ) then
         cosscat(v) = - vsign * calpha_boa ; if (calpha_boa .eq. zero) cosscat(v) = calpha_boa
      else
         term1 = salpha_boa * stheta_boa * cphi_boa
         term2 = calpha_boa * ctheta_boa
         cosscat(v) = - vsign * term2 + term1
      endif

!  End geometry routine

   enddo

!  Finish

   return
end subroutine Obsgeom_PlanPar

!  Finish Module

end module FO_geometry_Routines_Obs_PP_m
