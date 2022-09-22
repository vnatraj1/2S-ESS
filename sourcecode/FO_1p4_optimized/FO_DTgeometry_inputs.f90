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

module FO_DTGeometry_inputs_m

!  Stand alone geometry for Direct Thermal only

private
public FO_DTGeometry_inputs

contains

subroutine FO_DTGeometry_inputs  &
       ( maxgeoms, ngeoms, do_planpar, do_enhanced_ps, & ! Input dimensions/flags
         alpha_boa,                                    & ! Inputs
         fail, message, trace )                          ! Output(Status)

   implicit none

!  Parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms

!  Flags (sphericity flag is mutually exclusive)

   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  Geometry control

   integer  , intent(in)     :: ngeoms

!  Input angles (Degrees)

   real(ffp), intent(In)     :: alpha_boa(maxgeoms)

!  Output arguments
!  ================

!  Exception handling

   logical      , intent(out)    :: fail
   character*(*), intent(out)    :: message
   character*(*), intent(out)    :: trace

!  Help variables
!  --------------

   real(ffp), parameter   :: zero = 0.0_ffp
   character*2            :: c2
   integer                :: v

!  Initialize output

   fail = .false. ; message = ' ' ; trace = ' '

!  Check sphericity control
!  ------------------------

!  Cannot have Plane-parallel and Enhanced PS

   if ( do_planpar .and. do_enhanced_ps ) then
      message = 'Cannot have BOTH Plane-parallel and Enhanced PS options'
      trace   = 'Initial Flag Check in FO_DTGeometry_Master'
      fail    = .true. ;  return
   endif

!  Check geometry angles
!  ---------------------

!  VZA can be 0-90 degrees inclusive, but not outside this range

   do v = 1, ngeoms
      if ( alpha_boa(v) .gt. 90.0_ffp .or. alpha_boa(v) .lt. zero ) then
         write(c2,'(I2)')v
         message = 'Boa LOS angle outside range [0,90]); Check it!'
         trace   = 'Geometry # '//c2//'; Initial Angle Check in FO_DTGeometry_inputs'
         fail    = .true. ;  return
      endif
    enddo

!  Finish

   return
end subroutine FO_DTGeometry_inputs

!  Finish

end module FO_DTGeometry_inputs_m
