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

module FO_DTGeometry_Master_PP_RPS_m

!  Stand alone geometry for Direct Thermal only

private
public FO_DTGeometry_Master_PP_RPS

contains

subroutine FO_DTGeometry_Master_PP_RPS &
       ( maxgeoms, ngeoms, dtr, alpha_boa, & ! Inputs
         Mu1 )                               ! Output

   implicit none

!  parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)  :: maxgeoms

!  Geometry control

   integer  , intent(in)  :: ngeoms

!  dtr = degrees-to-Radians

   real(ffp), intent(in)  :: dtr

!  input angles (Degrees)

   real(ffp), intent(In)  :: alpha_boa(maxgeoms)

!  Output arguments
!  ================

!  Cosine

   real(ffp), intent(Out) :: Mu1     (maxgeoms)

!  Help variables
!  --------------

   integer                :: v
   real(ffp)              :: alpha_boa_R

!  Set cosines if you are here

   do v = 1, ngeoms
      alpha_boa_R  = alpha_boa(v) * DTR
      Mu1(v)       = cos(alpha_boa_R)
   enddo

!  Finish

   return
end subroutine FO_DTGeometry_Master_PP_RPS

!  Finish

end module FO_DTGeometry_Master_PP_RPS_m
