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
! #   Version 1.3a, 29 October  2012, Multiple geometries   #
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

module FO_ScalarSS_spherfuncs_optimized_m

!  Legendre polynomials

public

contains

SUBROUTINE FO_ScalarSS_spherfuncs ( MAXMOMS, MAXGEOMS, NMOMS, NGEOMS, COSSCAT, SS_PLEG )

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  I/O

   INTEGER  , intent(in)    :: MAXMOMS,  NMOMS
   INTEGER  , intent(in)    :: MAXGEOMS, NGEOMS
   REAL(fpk), intent(in)    :: COSSCAT(MAXGEOMS)
   REAL(fpk), intent(out)   :: SS_PLEG(0:MAXMOMS,MAXGEOMS)

!  Local

   integer   :: L, V 
   real(fpk) :: MU, SSPL, SSPL1, SSPL2, LFPK
     
!  Legendre

   DO V = 1, NGEOMS
     MU = COSSCAT(V)
     SSPL2 = 1.0_fpk
     SSPL1 = MU
     SS_PLEG(0,V) = SSPL2
     SS_PLEG(1,V) = SSPL1
     DO L = 2, NMOMS
        LFPK = REAL(L,KIND=fpk)
        SSPL = REAL(2*L-1,KIND=fpk) / LFPK * SSPL1 * MU - REAL(L-1,KIND=fpk) / LFPK * SSPL2
        SS_PLEG(L,V) = SSPL
        SSPL2 = SSPL1
        SSPL1 = SSPL
     ENDDO
   ENDDO

!  Finish

   RETURN
END SUBROUTINE FO_ScalarSS_spherfuncs

!  Finish

end module FO_ScalarSS_spherfuncs_optimized_m
