! ###########################################################
! #                                                         #
! #             THE TWOSTREAM LIDORT MODEL                  #
! #                                                         #
! #      (LInearized Discrete Ordinate Radiative Transfer)  #
! #       --         -        -        -         -          #
! #                                                         #
! ###########################################################

! ###########################################################
! #                                                         #
! #  Authors :      Robert. J. D. Spurr (1)                 #
! #                 Vijay Natraj        (2)                 #
! #                                                         #
! #  Address (1) :     RT Solutions, Inc.                   #
! #                    9 Channing Street                    #
! #                    Cambridge, MA 02138, USA             #
! #  Tel:             (617) 492 1183                        #
! #  Email :           rtsolutions@verizon.net              #
! #                                                         #
! #  Address (2) :     CalTech                              #
! #                    Department of Planetary Sciences     #
! #                    1200 East California Boulevard       #
! #                    Pasadena, CA 91125                   #
! #  Tel:             (626) 395 6962                        #
! #  Email :           vijay@gps.caltech.edu                #
! #                                                         #
! #  Version 1.0-1.3 :                                      #
! #     Mark 1: October  2010                               #
! #     Mark 2: May      2011, with BRDFs                   #
! #     Mark 3: October  2011, with Thermal sources         #
! #                                                         #
! #  Version 2.0-2.1 :                                      #
! #     Mark 4: November 2012, LCS/LPS Split, Fixed Arrays  #
! #     Mark 5: December 2012, Observation Geometry option  #
! #                                                         #
! #  Version 2.2-2.3 :                                      #
! #     Mark 6: July     2013, Level outputs + control      #
! #     Mark 7: December 2013, Flux outputs  + control      #
! #     Mark 8: January  2014, Surface Leaving + control    #
! #     Mark 9: June     2014, Inverse Pentadiagonal        #
! #                                                         #
! #  Version 2.4 :                                          #
! #     Mark 10: August  2014, Green's function Regular     #
! #     Mark 11: January 2015, Green's function Linearized  #
! #                            Taylor, dethreaded, OpenMP   #
! #                                                         #
! ###########################################################

! #############################################################
! #                                                           #
! #   This Version of LIDORT-2STREAM comes with a GNU-style   #
! #   license. Please read the license carefully.             #
! #                                                           #
! #############################################################

! ###############################################################
! #                                                             #
! # Regular BVP: Subroutines in this Module                     #
! #                                                             #
! #            TWOSTREAM_BVP_MATSETUP_PENTADIAG                 #
! #            TWOSTREAM_BVP_SOLUTION_PENTADIAG                 #
! #                                                             #
! ###############################################################
 
module twostream_bvproblem_m

PUBLIC

contains

SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG &
         ( MAXLAYERS, MAXTOTAL, FF, DO_INVERSE,                     & ! Dimensions
           DO_INCLUDE_SURFACE, NLAYERS, NTOTAL,                     & ! input
           DO_BRDF_SURFACE, SURFACE_FACTORM, ALBEDO, BRDF_FM,       & ! input
           XPOS, T_DELT_EIGEN, STREAM_VALUE,                        & ! input
           MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM,               & ! output
           STATUS, MESSAGE )                                          ! output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

!  input
!  -----

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  control

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTORM
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_FM

!  Eigenvector solutions

      REAL(kind=dp), INTENT(IN)  :: XPOS(2,MAXLAYERS)

!  transmittance factors for +/- eigenvalues

      REAL(kind=dp), INTENT(IN)  :: T_DELT_EIGEN(MAXLAYERS)

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Output
!  ------

!  Pentadiagonal Matrix entries for solving BCs
!  Split into 5 vectors

      REAL(kind=dp), INTENT(OUT) :: MAT1(3:MAXTOTAL)
      REAL(kind=dp), INTENT(OUT) :: MAT22

!  Pentadiagonal elimination marix
!  Split into 4 vectors

      REAL(kind=dp), INTENT(OUT) :: ELM1 (MAXTOTAL-1)
      REAL(kind=dp), INTENT(OUT) :: ELM2 (MAXTOTAL-2)
      REAL(kind=dp), INTENT(OUT) :: ELM3 (MAXTOTAL)
      REAL(kind=dp), INTENT(OUT) :: ELM4 (3:MAXTOTAL)

!  single layer elimination matrix

      REAL(kind=dp), INTENT(OUT) :: SELM (2,2)

!  status

      INTEGER, INTENT(OUT)       :: STATUS
      CHARACTER*(*), INTENT(OUT) :: MESSAGE

!  local variables
!  ---------------

!  Pentadiagonal Matrix entries for solving BCs
!  Local vectors

      REAL(kind=dp)              :: MAT2(2:MAXTOTAL)
      REAL(kind=dp)              :: MAT3(MAXTOTAL)
      REAL(kind=dp)              :: MAT4(MAXTOTAL)
      REAL(kind=dp)              :: MAT5(MAXTOTAL)

!  square matrix for the single layer case

      REAL(kind=dp)     :: SMAT11, SMAT12, SMAT21, SMAT22

!  Help

      INTEGER           :: I, N, N1, NM, NP, INM, INP
      REAL(kind=dp)     :: XPNET, XMNET, FACTOR, BET, DEN
      REAL(kind=dp)     :: H_HOMP, H_HOMM, R2_HOMP, R2_HOMM
      CHARACTER(LEN=3)  :: CI

!  Stability check value

      REAL(kind=dp)     :: SMALLNUM=1.0D-20

!  Other variables

      REAL(kind=dp)     :: T_DELT_EIGENN, T_DELT_EIGENN1
      REAL(kind=dp)     :: XPOS1N, XPOS2N, XPOS1N1, XPOS2N1
      REAL(kind=dp)     :: ELM31, ELM1i1, ELM1i2, ELM2i1, ELM2i2, MAT1i

!  Status

      STATUS = 0
      MESSAGE = ' '

!  Additional setups for the lowest layer
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      H_HOMP = XPOS(1,NLAYERS) * STREAM_VALUE
      H_HOMM = XPOS(2,NLAYERS) * STREAM_VALUE

      R2_HOMP = ZERO
      R2_HOMM = ZERO
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTORM * BRDF_FM
        ELSE
          FACTOR = SURFACE_FACTORM * ALBEDO
        ENDIF
        R2_HOMP = FACTOR * H_HOMP
        R2_HOMM = FACTOR * H_HOMM
      ENDIF

!  Inclusion of surface contribution in BV Problem matrix

      IF ( DO_INCLUDE_SURFACE ) THEN
        XPNET = XPOS(2,NLAYERS) - R2_HOMP
        XMNET = XPOS(1,NLAYERS) - R2_HOMM
      ELSE
        XPNET = XPOS(2,NLAYERS)
        XMNET = XPOS(1,NLAYERS)
      ENDIF

!  set up BVP matrix
!  -----------------

!  If Nlayers = 1, special case, ordinary 2x2

      IF ( NLAYERS .EQ. 1 ) THEN

!  Only Top and bottom BC (with surface reflection)

        T_DELT_EIGENN = T_DELT_EIGEN(NLAYERS)
        SMAT11 = XPOS(1,NLAYERS)
        SMAT12 = XPOS(2,NLAYERS) * T_DELT_EIGENN
        SMAT21 = XPNET * T_DELT_EIGENN
        SMAT22 = XMNET

!  If NLAYERS > 1, set up Pentadiagonal matrix

      ELSE

!  Top BC for layer 1: no downward diffuse radiation
!  Intermediate layer boundaries
!  Bottom BC (including surface reflected term)

!  Regular set
!  Don't need XNEG (can derive from XPOS). Define XPOS1N, XPOS1N1, XPOS2N and XPOS2N1 and recursively read them 
!  for all layers for efficient matrix reading.
!  Use separate arrays for MAT1-MAT5

        IF ( .NOT. DO_INVERSE ) THEN
          XPOS1N = XPOS(1,1)
          XPOS2N = XPOS(2,1)
          T_DELT_EIGENN = T_DELT_EIGEN(1)
          MAT3(1)  = XPOS1N
          MAT4(1)  = XPOS2N * T_DELT_EIGENN
          MAT5(1)  = ZERO
          DO N = 2, NLAYERS
            N1 =  N - 1
            XPOS1N1 = XPOS1N
            XPOS1N = XPOS(1,N)
            XPOS2N1 = XPOS2N
            XPOS2N = XPOS(2,N)
            T_DELT_EIGENN1 = T_DELT_EIGENN
            T_DELT_EIGENN = T_DELT_EIGEN(N)
            NM = 2*N1
            NP = NM + 1
            IF (N .GT. 2) MAT1(NM) = ZERO
            MAT2(NM) =   XPOS1N1 * T_DELT_EIGENN1
            MAT3(NM) =   XPOS2N1
            MAT4(NM) = - XPOS1N
            MAT5(NM) = - XPOS2N  * T_DELT_EIGENN
            MAT1(NP) =   XPOS2N1 * T_DELT_EIGENN1
            MAT2(NP) =   XPOS1N1
            MAT3(NP) = - XPOS2N
            MAT4(NP) = - XPOS1N  * T_DELT_EIGENN
            MAT5(NP) =   ZERO
          ENDDO
          MAT1(NTOTAL) = ZERO
          MAT2(NTOTAL) = XPNET * T_DELT_EIGEN(NLAYERS)
          MAT3(NTOTAL) = XMNET
          MAT4(NTOTAL) = ZERO
          MAT5(NTOTAL) = ZERO
        endif

!  Inverted set

        IF ( DO_INVERSE ) THEN
          XPOS1N = XPOS(1,1)
          XPOS2N = XPOS(2,1)
          T_DELT_EIGENN = T_DELT_EIGEN(1)
          MAT1(NTOTAL) = ZERO
          MAT2(NTOTAL) = XPOS2N * T_DELT_EIGENN
          MAT3(NTOTAL) = XPOS1N
          MAT4(NTOTAL) = ZERO
          MAT5(NTOTAL) = ZERO
          DO N = 2, NLAYERS
            N1 =  N - 1
            XPOS1N1 = XPOS1N  
            XPOS1N = XPOS(1,N)
            XPOS2N1 = XPOS2N  
            XPOS2N = XPOS(2,N)
            T_DELT_EIGENN1 = T_DELT_EIGENN
            T_DELT_EIGENN = T_DELT_EIGEN(N)
            INM = NTOTAL - 2*N1  ; INP = INM + 1
            IF (N .LT. NLAYERS) MAT1(INM) = ZERO
            MAT2(INM) =  - XPOS1N * T_DELT_EIGENN
            MAT3(INM) =  - XPOS2N
            MAT4(INM) =    XPOS1N1
            MAT5(INM) =    XPOS2N1 * T_DELT_EIGENN1
            MAT1(INP) =  - XPOS2N  * T_DELT_EIGENN
            MAT2(INP) =  - XPOS1N
            MAT3(INP) =    XPOS2N1
            MAT4(INP) =    XPOS1N1 * T_DELT_EIGENN1
            MAT5(INP) =    ZERO
          ENDDO
          MAT3(1)  = XMNET
          MAT4(1)  = XPNET * T_DELT_EIGEN(NLAYERS)
          MAT5(1)  = ZERO
        endif

      ENDIF

!  Scaling

      MAT1 = FF * MAT1
      MAT2 = FF * MAT2
      MAT3 = FF * MAT3
      MAT4 = FF * MAT4
      MAT5 = FF * MAT5

!  Elimination of BVP pentadiagonal matrix
!  ---------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

!  Row 1

        ELM31 = ONE / MAT3(1)
        ELM3(1) = ELM31
        ELM1i2  = - MAT4(1) * ELM31
        ELM1(1) = ELM1i2
        ELM2i2 = - MAT5(1) * ELM31
        ELM2(1) = ELM2i2

!  Row 2; includes first check for singularity

        MAT22   = MAT2(2)
        BET = MAT3(2) + MAT22*ELM1i2
        IF ( DABS(BET) .LT. SMALLNUM ) THEN
          MESSAGE = 'Singularity in Pentadiagonal Matrix, Row #  2'
          STATUS = 1
          RETURN
        endif
        BET = - ONE / BET
        ELM1i1  = (MAT4(2) + MAT22*ELM2i2) * BET
        ELM1(2) = ELM1i1
        ELM2i1  = MAT5(2) * BET
        ELM2(2) = ELM2i1
        ELM3(2) = BET

!  Rows 3-NT-2: reduce to upper triangular; includes checks for singularity

        DO I = 3, NTOTAL-2
          MAT1i = MAT1(i)
          BET = MAT2(i) + MAT1i * ELM1i2
          DEN = MAT3(i) + MAT1i * ELM2i2 + BET * ELM1i1
          IF ( DABS(DEN) .LT. SMALLNUM ) THEN
            WRITE(CI, '(I3)' ) I
            MESSAGE = 'Singularity in Pentadiagonal Matrix, Row #'//CI
            STATUS = 1
            RETURN
          ENDIF
          DEN = - ONE / DEN
          ELM1i2  = ELM1i1
          ELM1i1  = (MAT4(i) + BET*ELM2i1) * DEN
          ELM1(i) = ELM1i1
          ELM2i2  = ELM2i1
          ELM2i1  = MAT5(i) * DEN
          ELM2(i) = ELM2i1
          ELM3(i) = BET
          ELM4(i) = DEN
        ENDDO

!  Row NT-1

        I = NTOTAL-1
        MAT1i = MAT1(i)
        BET = MAT2(i) + MAT1i * ELM1i2
        DEN = MAT3(i) + MAT1i * ELM2i2 + BET * ELM1i1
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          WRITE(CI, '(I3)' ) I
          MESSAGE = 'Singularity in Pentadiagonal Matrix, Row #'//CI
          STATUS = 1
          RETURN
        ENDIF
        DEN = - ONE / DEN
        ELM1i2  = ELM1i1
        ELM1i1  = (MAT4(i) + BET*ELM2i1) * DEN
        ELM1(i) = ELM1i1
        ELM2i2  = ELM2i1
        ELM3(i) = BET
        ELM4(i) = DEN

!  Row NT

        I = NTOTAL
        MAT1i = MAT1(i)
        BET = MAT2(i) + MAT1i * ELM1i2
        DEN = MAT3(i) + MAT1i * ELM2i2 + BET * ELM1i1
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          WRITE(CI, '(I3)' ) I
          MESSAGE = 'Singularity in Pentadiagonal Matrix, Row #'//CI
          STATUS = 1
          RETURN
        ENDIF
        DEN = - ONE / DEN
        ELM3(i) = BET
        ELM4(i) = DEN

!  Elimination for Single Layer only
!  ----------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

        DEN = SMAT11*SMAT22 - SMAT12*SMAT21
        IF ( DABS(DEN) .LT. SMALLNUM ) THEN
          MESSAGE = 'Singularity in 1-layer 2x2 Matrix'
          STATUS = 1
          RETURN
        ENDIF
        DEN = ONE / DEN
        SELM(1,1) =   SMAT22 * DEN
        SELM(1,2) = - SMAT12 * DEN
        SELM(2,1) = - SMAT21 * DEN
        SELM(2,2) =   SMAT11 * DEN

      ENDIF

!  Finish

      RETURN
 END SUBROUTINE TWOSTREAM_BVP_MATSETUP_PENTADIAG

!

SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG &
      ( MAXLAYERS, MAXTOTAL, FF, DO_INVERSE,             & ! Dimensions
        DO_INCLUDE_SURFACE, DO_INCLUDE_DIRECTBEAMB,      & ! inputs
        DO_INCLUDE_SURFEMISS, DO_BRDF_SURFACE,           & ! inputs
        NLAYERS, NTOTAL,                                 & ! inputs
        SURFACE_FACTORM, ALBEDO, BRDF_FM, EMISS, SURFBB, & ! inputs
        DIRECT_BEAMB, WUPPER, WLOWER,                    & ! inputs
        STREAM_VALUE, MAT1, MAT22, ELM1, ELM2, ELM3, ELM4, SELM, & ! inputs
        LCON, MCON )                                       ! Output

      implicit none

!  precision and parameters

      INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
      REAL(kind=dp), parameter :: zero = 0.0_dp

!  Input arguments
!  ---------------

!  Dimensions

      INTEGER, INTENT(IN)        :: MAXLAYERS, MAXTOTAL

!  Inverse control

      LOGICAL      , INTENT(IN)  :: DO_INVERSE
      REAL(kind=dp), INTENT(IN)  :: FF

!  Inclusion flags

      LOGICAL, INTENT(IN)        :: DO_INCLUDE_DIRECTBEAMB
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFACE
      LOGICAL, INTENT(IN)        :: DO_INCLUDE_SURFEMISS

!  Surface control

      LOGICAL      , INTENT(IN)  :: DO_BRDF_SURFACE
      REAL(kind=dp), INTENT(IN)  :: SURFACE_FACTORM
      REAL(kind=dp), INTENT(IN)  :: ALBEDO
      REAL(kind=dp), INTENT(IN)  :: BRDF_FM
      REAL(kind=dp), INTENT(IN)  :: EMISS, SURFBB

!  Numbers

      INTEGER, INTENT(IN)        :: NLAYERS, NTOTAL

!  Direct beam

      REAL(kind=dp), INTENT(IN)  :: DIRECT_BEAMB

!  Particular solutions

      REAL(kind=dp), INTENT(IN)  :: WLOWER ( 2, MAXLAYERS )
      REAL(kind=dp), INTENT(IN)  :: WUPPER ( 2, MAXLAYERS )

!  Stream

      REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

!  Pentadiagonal Matrix entries for solving BCs
!  Only need MAT1 and MAT22

      REAL(kind=dp), INTENT(IN)  :: MAT1(3:MAXTOTAL)
      REAL(kind=dp), INTENT(IN)  :: MAT22

!  Pentadiagonal elimination matrix
!  Split into 4 vectors

      REAL(kind=dp), INTENT(IN)  :: ELM1 (MAXTOTAL-1)
      REAL(kind=dp), INTENT(IN)  :: ELM2 (MAXTOTAL-2)
      REAL(kind=dp), INTENT(IN)  :: ELM3 (MAXTOTAL)
      REAL(kind=dp), INTENT(IN)  :: ELM4 (3:MAXTOTAL)

!  Single layer elimination matrix

      REAL(kind=dp), INTENT(IN)  :: SELM (2,2)

!  Output
!  ------

!  Solution constants of integration, and related quantities

      REAL(kind=dp), INTENT(OUT) :: LCON(MAXLAYERS)
      REAL(kind=dp), INTENT(OUT) :: MCON(MAXLAYERS)

!  Local variables
!  ---------------

!  Column vectors for solving BCs. Not saved.

      REAL(kind=dp)       :: COL    (MAXTOTAL)
      REAL(kind=dp)       :: SCOL   (2)

!  Other variables

      INTEGER             :: N, N1, I, NM, NP, INP, INM, NI, NLAY1
      REAL(kind=dp)       :: FACTOR, H_PARTIC, R2_PARTIC, TOA_TERM, SURFACE_TERM
      REAL(kind=dp)       :: NEW_SCOL1

      REAL(kind=dp)       :: COLi, COLi1, COLi2
      REAL(kind=dp)       :: SCOL1, SCOL2

!  Additional setups for the surface and TOA levels
!  ------------------------------------------------

!  Zero total reflected contribution (R2_PARTIC) before calculation
!  For Lambertian reflectance, all streams are the same
!  For BRDF, code added 4 May 2009 by R. Spurr

      R2_PARTIC = ZERO
      H_PARTIC = WLOWER(1,NLAYERS) * STREAM_VALUE
      IF ( DO_INCLUDE_SURFACE ) THEN
        IF ( DO_BRDF_SURFACE ) THEN
          FACTOR = SURFACE_FACTORM * BRDF_FM
        ELSE
          FACTOR = SURFACE_FACTORM * ALBEDO
        ENDIF
        R2_PARTIC = H_PARTIC * FACTOR
      ENDIF

!  Surface Term------------
!  lowest (surface) boundary with albedo (diffuse radiation terms only)
!  with non-zero albedo, include integrated downward reflectances
!  no albedo, similar code excluding integrated reflectance
!  Add direct beam solution (only to final level)
!  Add thermal emission of ground surface (only to final level)
!  Treatment of emissivity consistent with LIDORT & VLIDORT

      SURFACE_TERM = - WLOWER(2,NLAYERS)
      IF ( DO_INCLUDE_SURFACE ) THEN
        SURFACE_TERM = SURFACE_TERM + R2_PARTIC
        IF ( DO_INCLUDE_DIRECTBEAMB ) THEN
          SURFACE_TERM = SURFACE_TERM + DIRECT_BEAMB
        ENDIF
      ENDIF
      IF ( DO_INCLUDE_SURFEMISS ) THEN
        SURFACE_TERM = SURFACE_TERM + SURFBB * EMISS
      ENDIF

!  Upper boundary for layer 1: no downward diffuse radiation

      TOA_TERM  = - WUPPER(1,1)

!  --set up Column for solution vector (the "B" as in AX=B)
!  --------------------------------------------------------

!  Pentadiagonal case

      IF ( NLAYERS .GT. 1 ) THEN

!  Fill vector

        IF ( DO_INVERSE ) THEN
          COL(1)  = SURFACE_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            INM = NTOTAL - 2*N1 ; INP = INM + 1
            COL(INP) = WUPPER(1,N) - WLOWER(1,N1)
            COL(INM) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = TOA_TERM
        ELSE
          COL(1)  = TOA_TERM
          DO N = 2, NLAYERS
            N1 = N - 1
            NM = 2*N1 ; NP = NM + 1
            COL(NM) = WUPPER(1,N) - WLOWER(1,N1)
            COL(NP) = WUPPER(2,N) - WLOWER(2,N1)
          ENDDO
          COL(NTOTAL) = SURFACE_TERM
        ENDIF

!  Scaling

        COL(1:NTOTAL) = FF * COL(1:NTOTAL)

!  Debug
!        do n = 1, ntotal
!          write(44,'(2i3,1p3e24.12)')n,1,COL(N)
!        enddo
!       pause

!  If Nlayers = 1, special case

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        SCOL(1) = TOA_TERM
        SCOL(2) = SURFACE_TERM
      ENDIF

!  Solve the boundary problem for this Fourier component (back substitution)
!  -------------------------------------------------------------------------

!  Pentadiagonal back-substitution

      IF ( NLAYERS .GT. 1 ) THEN

!  Fill up back-substitution array
!  Use separate arrays for ELM1-ELM4

        COLi2  = COL(1) * ELM3(1)
        COL(1) = COLi2
        COLi1  = (MAT22*COLi2 - COL(2)) * ELM3(2)
        COL(2) = COLi1
         
        DO I = 3, NTOTAL
          COLi = (MAT1(i)*COLi2+ELM3(i)*COLi1-COL(i)) * ELM4(i)
          COLi2 = COLi1
          COLi1 = COLi
          COL(i) = COLi
        ENDDO

!  Back-substitution

        N1 = NTOTAL-1
        COLi = COLi2 + ELM1(N1) * COLi1
        COL(N1) = COLi
        COLi2 = COLi1
        COLi1 = COLi
        DO I = NTOTAL-2, 1, -1
          COLi = COL(i) + ELM1(i) * COLi1 + ELM2(i) * COLi2
          COL(i) = COLi
          COLi2 = COLi1
          COLi1 = COLi
        ENDDO

!  Set integration constants LCON and MCON for -/+ eigensolutions, all layers

        IF ( DO_INVERSE ) THEN
           NLAY1 = 1 + NLAYERS
           DO N = 1, NLAYERS
              NI = NLAY1 - N
              INP = 2*NI ; INM = INP - 1
              LCON(N) = COL(INP)
              MCON(N) = COL(INM)
           ENDDO
        ELSE
           DO N = 1, NLAYERS
              NM = 2*N-1 ; NP = NM + 1
              LCON(N) = COL(NM)
              MCON(N) = COL(NP)
           ENDDO
        ENDIF

!  Solve the boundary problem: Single Layer only
!  Defined NEW_SCOL1 so a MODIFIED version of SCOL(1) is not being used in computing SCOL(2)

      ELSE IF ( NLAYERS .EQ. 1 ) THEN
        SCOL1 = SCOL(1)
        SCOL2 = SCOL(2)
        NEW_SCOL1 = SELM(1,1) * SCOL1 + SELM(1,2) * SCOL2
        SCOL(2)   = SELM(2,1) * SCOL1 + SELM(2,2) * SCOL2
        LCON(1) = NEW_SCOL1
        MCON(1) = SCOL(2)
      ENDIF

!  Debug
!      if ( fourier.eq.0 .and.ipartic.eq.1) then
!        do n = 1, nlayers
!          write(34,'(i3,1p2e24.12)')n,LCON(N),MCON(N)
!        enddo
!      endif

!  Finish

      RETURN
END SUBROUTINE TWOSTREAM_BVP_SOLUTION_PENTADIAG

end module twostream_bvproblem_m
