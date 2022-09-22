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

module FO_SSGeometry_inputs_m

!  Stand alone geometry for solar scattering only

private
public FO_SSGeometry_inputs

contains

subroutine FO_SSGeometry_inputs &
       ( maxgeoms, ngeoms,           & ! Input dimension, control
         do_planpar, do_enhanced_ps, & ! Input flags
         obsgeom_boa,                & ! Input
         fail, message, trace )        ! Output(Status)

   implicit none

!  Parameter argument

   integer, parameter :: ffp = selected_real_kind(15)

!  Input arguments
!  ===============

!  Dimensions

   integer  , intent(in)     :: maxgeoms

!  Geometry control

   integer  , intent(in)     :: ngeoms

!  Flags (sphericity flag is mutually exclusive)

   logical  , intent(in)     :: do_enhanced_ps
   logical  , intent(in)     :: do_planpar

!  input angles (Degrees)
!  Convention for ObsGeom = same as VLIDORT/LIDORT (1=sza,2=vza,3=azm)
!  In both cases, the Phi angle may be changed.....

   real(ffp), intent(InOut)  :: Obsgeom_boa(maxgeoms,3)

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
      trace   = 'Initial Flag Check in SSGeometry_inputs'
      fail    = .true. ;  return
   endif

!  Check geometry angles
!  ---------------------

      do v = 1, ngeoms

         if ( obsgeom_boa(v,2) .gt. 90.0_ffp .or. obsgeom_boa(v,2) .lt. zero ) then
            write(c2,'(I2)') v
            message = 'Boa LOS angle outside range [0,90]); Check it!'
            trace   = 'Geometry # '//c2//'; Initial Angle Check in SSGeometry_inputs, OBSGEOM mode'
            fail    = .true. ;  return
         endif

         if ( do_planpar ) then
            if ( obsgeom_boa(v,1) .ge. 90.0_ffp .or. obsgeom_boa(v,1) .lt. zero ) then
               write(c2,'(I2)') v
               message = 'Plane-parallel: Boa SZA angle outside range [0,90)); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in SSGeometry_inputs, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         else
            if ( obsgeom_boa(v,1) .gt. 90.0_ffp .or. obsgeom_boa(v,1) .lt. zero ) then
               write(c2,'(I2)') v
               message = 'Pseudo-spherical : Boa SZA angle outside range [0,90]); Check it!'
               trace   = 'Geometry # '//c2//'; Initial Angle Check in SSGeometry_inputs, OBSGEOM mode'
               fail    = .true. ;  return
            endif
         endif

         if ( obsgeom_boa(v,3) .lt. zero )    obsgeom_boa(v,3) = - obsgeom_boa(v,3)
         if ( obsgeom_boa(v,3) .gt. 360.0_ffp ) obsgeom_boa(v,3) = obsgeom_boa(v,3) - 360.0_ffp

!  End loop over Obs geometries

      enddo

!  Finish

   return
end subroutine FO_SSGeometry_inputs

!  Finish

end module FO_SSGeometry_inputs_m
