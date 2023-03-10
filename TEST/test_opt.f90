program test_opt

!  Two-stream (Version 2.4) optimized modules

   USE twostream_inputs_solar_obs_m
   USE twostream_master_solar_obs_m

!  First Order (Version 1.4) optimized modules

   use FO_SSGeometry_inputs_m
   use FO_SSGeometry_Master_Obs_PP_m
   use FO_SSGeometry_Master_Obs_RPS_m
   use FO_SSGeometry_Master_Obs_EPS_m
   use FO_ScalarSS_spherfuncs_optimized_m
   use FO_ScalarSS_RTCalcs_I_PP_m , only : SS_Integral_I_UP_PP
   use FO_ScalarSS_RTCalcs_I_RPS_m, only : SS_Integral_I_UP_RPS
   use FO_ScalarSS_RTCalcs_I_EPS_m, only : SS_Integral_I_UP_EPS

   implicit none 

!  precision parameters
!  ====================

   INTEGER, PARAMETER :: sp     = selected_real_kind(6)
   INTEGER, PARAMETER :: dp     = selected_real_kind(15)

!  Dimensions
!  ==========

   INTEGER, parameter ::  E_nlayers   = 114
   INTEGER, parameter ::  E_nmoms_all = 5000

   INTEGER, parameter ::  E_ngeoms = 1
   INTEGER, parameter ::  E_nszas = 1
   INTEGER, parameter ::  E_nvzas = 1
   INTEGER, parameter ::  E_nazms = 1

!  Monitoring flag

   logical, parameter :: Monitor_CPU = .true.

!  Control integers
!  nlayers   = Number of layers
!  nmoms_all = Number of Moments in FO calculation

   integer :: nlayers
   integer :: ndat
   integer :: nmoms_all

!  Observational-geometry flag

   logical :: do_obsgeoms

!  Enhanced sphericity flag
!  Number of LIDORT discrete ordinates
!  2STREAM BVP INdex (0 = LAPACK,1 = PentaDiag) and inverse flag

   logical :: do_enhanced_ps

   logical :: DO_PDINVERSE
   INTEGER :: BVPINDEX

!  Earth radius

   real(kind=dp) :: eradius

!  ngeoms              = Number of geometries (New, Mk 6)
!  nszas, nvzas, nazms = Number of angles     (New, 8/6/13)

   integer :: ngeoms
   integer :: nszas
   integer :: nvzas
   integer :: nazms

!  Geometry

   real(kind=dp) :: sza_boa(E_nszas), vza_boa(E_nvzas), azm_boa(E_nazms)

   real(kind=dp) :: albedo

   real(kind=dp) :: taudp  (E_nlayers)
   real(kind=dp) :: omega  (E_nlayers)

   real(kind=dp) :: heights(0:E_nlayers)

!  Intensity and Status output results
!  ===================================

!  Variables for Exact intensity

   real(kind=dp) :: intensity_2S_Exact(100000,E_ngeoms)
   real(kind=dp) :: intensity_FO_Exact(100000,E_ngeoms)

!  timing variables. Updated 5/7/15, 9 entries

   real(kind=sp) :: Exacttimes(7), ExactRunTime

!  Exception handling variables

   integer, parameter :: b_max_messages   = 25
   logical            :: b_fail
   integer            :: b_nmessages
   character*100      :: b_messages (b_max_messages)

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          2STREAM ARGUMENTS
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Version 2.4 code.

!  MAXLAYERS

   integer, parameter :: MAXLAYERS = E_nlayers

!  MAXTOTAL

   integer, parameter :: MAXTOTAL = 2 * MAXLAYERS

!  MAX_MESSAGES

   integer, parameter :: MAX_MESSAGES = 25

!  New for  Version 2.1: Observational Geometry

   INTEGER, PARAMETER :: MAX_USER_OBSGEOMS = E_ngeoms

!  BRDF dimension

   INTEGER, PARAMETER :: MAXSTREAMS_BRDF = 50

!  ntotal       = 2 * nlayers

   integer   :: ntotal

!  2Stream Directional Flags

      LOGICAL           :: DO_UPWELLING, DO_DNWELLING

!  2Stream Plane parallel and deltam-2stream scaling flags

      LOGICAL           :: DO_PLANE_PARALLEL, DO_D2S_SCALING

!  2Stream BRDF surface flag

      LOGICAL           :: DO_BRDF_SURFACE

!  @@ Rob Spurr, 17 July 2013, Version 2.2, Levelout flag

      LOGICAL           :: DO_2S_LEVELOUT     ! @@ 2p2

!  @@ Rob Spurr, 05 November 2013, Version 2.3, Flux option flags

      LOGICAL           :: DO_MVOUT_ONLY       ! @@ 2p3
      LOGICAL           :: DO_ADDITIONAL_MVOUT ! @@ 2p3

!  @@ Rob Spurr, 23 January 2014, Version 2.3, SLEAVE option flags

      LOGICAL           :: DO_SURFACE_LEAVING
      LOGICAL           :: DO_SL_ISOTROPIC

!  Other flags

      LOGICAL           :: DO_INCLUDE_MVOUT(0:1)

!  Surface reflectance factors

      REAL(kind=dp)     :: SURFACE_FACTOR(0:1)
      REAL(kind=dp)     :: DELTA_FACTOR(0:1)

!  Order of Taylor series (including terms up to EPS^n). Version 2.4

      INTEGER           :: TAYLOR_ORDER
      REAL(kind=dp)     :: TAYLOR_SMALL

!  BVP control --- New 6/25/14, Version 2.3 and higher
!  * BVP Scale Factor. Debug only. Set this to 1.0 on input

      REAL(kind=dp)   :: BVPSCALEFACTOR

!  2Stream Geometry. Observational Geometry Input (New for  Version 2.1)
!  If set, this will override Lattice-geometry input

      INTEGER       :: N_USER_OBSGEOMS                      !@@
      REAL(kind=dp) :: USER_OBSGEOMS(MAX_USER_OBSGEOMS,3) !@@

!  2Stream Stream value

      REAL(kind=dp)   :: STREAM_VALUE

!  2Stream Lambertian Surface control

      REAL(kind=dp)   :: LAMBERTIAN_ALBEDO

!  2Stream BRDF Fourier components
!  0 and 1 Fourier components of BRDF, following order
!  incident solar directions,  reflected quadrature stream
!  incident quadrature stream, reflected quadrature stream
!  incident quadrature stream, reflected user streams

      REAL(kind=dp)   :: BRDF_F_0  ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp)   :: BRDF_F    ( 0:1 )
      REAL(kind=dp)   :: UBRDF_F   ( MAX_USER_OBSGEOMS, 0:1 )

!  Version 2p3. 1/23/14. Introduce SLEAVE stuff
!  Do not require any first-order inputs (exact or Fourier)
!  Isotropic Surface leaving term (if flag set)
!  Fourier components of Surface-leaving terms:

      REAL(kind=dp)   ::  SLTERM_ISOTROPIC ( MAX_USER_OBSGEOMS )
      REAL(kind=dp)   ::  SLTERM_F_0 ( MAX_USER_OBSGEOMS, 0:1 )

!  Number of Fourier terms

      INTEGER         :: N_FOURIERS

!  Constants

      REAL(kind=dp)   :: DEG_TO_RAD, PI4

!  Post-processing flag (new for Version 2p3)

      LOGICAL         :: DO_POSTPROCESSING

!  Average secant (plane-parallel case)

      REAL(kind=dp)   :: AVERAGE_SECANT_PP ( MAX_USER_OBSGEOMS )

!  Chapman factors (from pseudo-spherical geometry)

      REAL(kind=dp)   :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAX_USER_OBSGEOMS )

!  Cosines and sines

      REAL(kind=dp)   :: X0  ( MAX_USER_OBSGEOMS )
      REAL(kind=dp)   :: USER_STREAMS ( MAX_USER_OBSGEOMS )
      REAL(kind=dp)   :: USER_SECANTS ( MAX_USER_OBSGEOMS )

!  Polynomial expansions

      REAL(kind=dp)   :: PX11, PXSQ(0:1)
      REAL(kind=dp)   :: PX0X ( MAX_USER_OBSGEOMS, 0:1 )
      REAL(kind=dp)   :: ULP ( MAX_USER_OBSGEOMS )

!  Azimuth cosine factor

      REAL(kind=dp)   :: AZMFAC ( MAX_USER_OBSGEOMS )

!  2Stream Flux factor

      REAL(kind=dp)   :: FLUX_FACTOR

!  2Stream height and earth radius

      REAL(kind=dp)   :: EARTH_RADIUS
      REAL(kind=dp)   :: HEIGHT_GRID ( 0:MAXLAYERS )

!  2Stream Atmospheric Optical properties

      REAL(kind=dp)   :: DELTAU_INPUT(MAXLAYERS)
      REAL(kind=dp)   :: OMEGA_INPUT (MAXLAYERS)
      REAL(kind=dp)   :: ASYMM_INPUT (MAXLAYERS)
      REAL(kind=dp)   :: D2S_SCALING (MAXLAYERS)

!  2Stream Results

      REAL(kind=dp)   :: INTENSITY_TOA(MAX_USER_OBSGEOMS)
      REAL(kind=dp)   :: INTENSITY_BOA(MAX_USER_OBSGEOMS)

!  Flux output
!  ! @@ Rob Spurr, 05 November 2013, Version 2.3 --> Flux Output

     REAL(kind=dp)    :: FLUXES_TOA_SOLAR(2,MAX_USER_OBSGEOMS)
     REAL(kind=dp)    :: FLUXES_BOA_SOLAR(2,MAX_USER_OBSGEOMS)

!  output solutions at ALL levels
!  ! @@ Rob Spurr, 17 July 2013, Version 2.2 --> Optional Output at ALL LEVELS

      REAL(kind=dp)    :: RADLEVEL_UP (MAX_USER_OBSGEOMS,0:MAXLAYERS)
      REAL(kind=dp)    :: RADLEVEL_DN (MAX_USER_OBSGEOMS,0:MAXLAYERS)

!  2Stream Exception handling

!    1. Check Messages and actions

      INTEGER          :: STATUS_INPUTCHECK
      INTEGER          :: C_NMESSAGES
!mick fix 5/29/2015 - adjusted indexing to include 0
      CHARACTER*100    :: C_MESSAGES(0:MAX_MESSAGES)
      CHARACTER*100    :: C_ACTIONS (0:MAX_MESSAGES)

!    2. Execution message and 2 Traces

      INTEGER          :: STATUS_EXECUTION
      CHARACTER*100    :: E_MESSAGE, E_TRACE_1, E_TRACE_2

!   2Stream  ARGUMENTS (Not I/O to the model)

!      LOGICAL, parameter :: DO_FULLQUADRATURE = .false.
      LOGICAL, parameter :: DO_FULLQUADRATURE = .true.

!  BRDF variables
!  ==============

!  Number of azimuth quadrature streams for BRDF

      INTEGER          :: NSTREAMS_BRDF

!  Constants

      REAL(kind=dp)    :: PIE

!  2-Stream angle sine

      REAL(kind=dp)    :: STREAM_SINE

!  Solar zenith sine

      REAL(kind=dp)    :: SX0(MAX_USER_OBSGEOMS)

!  Viewing zenith sine

      REAL(kind=dp)    :: USER_SINES(MAX_USER_OBSGEOMS)

!  BRDF azimuth quadrature streams

      INTEGER          :: NBRDF_HALF
      REAL(kind=dp)    :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=dp)    :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)    :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=dp)    :: A_BRDF  ( MAXSTREAMS_BRDF )

!  SLEAVE variables
!  ================

!  Isotropic flag

      LOGICAL          :: DO_ISOTROPIC

!  Flo flag

      LOGICAL          :: DO_FLUORESCENCE

!  Fluorescence variables
!  ----------------------

!  Input Latitude/Longitude in [degs]

      REAL(kind=dp)    :: FL_Latitude, FL_Longitude

!  Input Epoch

      INTEGER          :: FL_Epoch(6)

!  Flag for using Data Gaussians

      LOGICAL          :: FL_DO_DataGaussian

!  Input Gaussians (alternative to Data Gaussians)

      REAL(kind=dp)    ::  FL_InputGAUSSIANS(3,2)

!  Output Gaussians

      REAL(kind=dp)    ::  FL_GAUSSIANS(3,2)

!  Fluorescence at 755 nm

      REAL(kind=dp)   :: Fs755(MAX_USER_OBSGEOMS)

!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!          FIRST ORDER VARIABLES
!  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  FIFTH REVISION: Complete relabeling
!  --- Keep distinction of variables

!  Dimensions - Use Fine Layering, One output level

   integer, parameter :: FO_maxgeoms      = MAX_USER_OBSGEOMS
   integer, parameter :: FO_maxlayers     = MAXLAYERS
   integer, parameter :: FO_maxmoments    = E_nmoms_all
   integer, parameter :: FO_maxuserlevels = 1                      ! TOA only
   integer, parameter :: FO_maxfine       = 6                      ! This will be sufficient

!  Critical attenuation

   real(dp), parameter :: FO_Acrit = 1.0d-10

!  Constants

   real(dp), parameter  :: FO_Pie = ACOS(-1.0d0)
   real(dp), parameter  :: FO_dtr = FO_Pie/180.0d0

!  2. Control Variables
!  --------------------

!  Derived flags, optical settings
!  @@ Rob 8/6/13. Modification

   logical    :: FO_do_enhanced_ps
   logical    :: FO_do_planpar
   logical    :: FO_do_deltam_scaling
   logical    :: FO_do_lambertian
   logical    :: FO_do_obsgeom

!  Layer control. Finelayer input
!  @@ Rob 8/6/13. Modification

   integer    :: FO_nlayers
   integer    :: FO_nfineinput

!  Number of user levels

   integer    :: FO_n_user_levels
   integer    :: FO_user_levels(FO_maxuserlevels)

!  3. Geometry variables
!  ---------------------

!  Number of geometries and angles (Lattice option)
!  @@ Rob 8/6/13. Modification

   integer    :: FO_ngeoms

!  Radius + heights

   real(dp)  :: FO_eradius, FO_heights (0:FO_maxlayers)

!  input angles (Degrees), VSIGN = +1 (Up); -1(Down)
!  @@ Rob 8/6/13. Modification

   real(dp)  :: FO_vsign
   real(dp)  :: FO_obsgeom_boa(FO_maxgeoms,3)

!  Critical adjustment for cloud layers

   logical   :: FO_doCrit
   real(dp)  :: FO_extincs(FO_maxlayers)

!  Flag for the Nadir case.

   logical   :: FO_doNadir(FO_maxgeoms)
  
!  Alphas,  Cotangents, Radii, Ray constant.

   real(dp)  :: FO_radii    (0:FO_maxlayers)
   real(dp)  :: FO_Raycon   (FO_maxgeoms)
   real(dp)  :: FO_alpha    (0:FO_maxlayers,FO_maxgeoms)
   real(dp)  :: FO_cota     (0:FO_maxlayers,FO_maxgeoms)

!  LOS Quadratures for Enhanced PS

   integer   :: FO_nfinedivs(FO_maxlayers,FO_maxgeoms)
   real(dp)  :: FO_xfine    (FO_maxfine,FO_maxlayers,FO_maxgeoms)
   real(dp)  :: FO_wfine    (FO_maxfine,FO_maxlayers,FO_maxgeoms)
   real(dp)  :: FO_csqfine  (FO_maxfine,FO_maxlayers,FO_maxgeoms)
   real(dp)  :: FO_cotfine  (FO_maxfine,FO_maxlayers,FO_maxgeoms)

!  Critical layer

   integer  :: FO_Ncrit(FO_maxgeoms)
   real(dp) :: FO_RadCrit(FO_maxgeoms), FO_CotCrit(FO_maxgeoms)

!  solar paths, Intent out 

   integer  :: FO_ntraversenl    (FO_maxgeoms)
   real(dp) :: FO_sunpaths       (FO_maxlayers,FO_maxlayers,FO_maxgeoms)
   real(dp) :: FO_sunpathsnl     (FO_maxlayers,FO_maxgeoms)

   integer  :: FO_ntraverse_fine (FO_maxfine,FO_maxlayers,FO_maxgeoms)
   real(dp) :: FO_sunpaths_fine  (FO_maxlayers,FO_maxfine,FO_maxlayers,FO_maxgeoms)

!  Cosine scattering angle and other cosines

   real(dp) :: FO_cosscat(FO_maxgeoms) 
   real(dp) :: FO_Mu1    (FO_maxgeoms)
   real(dp) :: FO_Mu0    (FO_maxgeoms)

!  Spherical Functions: Inputs

   integer   :: FO_NMOMENTS

!  Spherical Functions: Output (Legendre Polynomials)

   real(dp) :: FO_GENSPHER(0:FO_MAXMOMENTS,FO_maxgeoms)

!  4. Optical properties
!  ---------------------

!  optical inputs Atmosphere

   real(dp) :: FO_DELTAUS     ( FO_MAXLAYERS )
   real(dp) :: FO_EXACTSCAT   ( FO_MAXLAYERS, FO_maxgeoms )

!  Solar Flux and Surface reflectivity (Could be the albedo)

   real(dp) :: FO_REFLEC(FO_maxgeoms), FO_FLUX

!  5. Output Variables and arrays
!  ------------------------------

!  Solar SS: First-order Radiances

   real(dp) :: FO_intensity_ss     ( FO_maxuserlevels, FO_maxgeoms )
   real(dp) :: FO_intensity_db     ( FO_maxuserlevels, FO_maxgeoms )

!  Exception handling on Geometry routine

   logical         :: FO_fail
   character*100   :: FO_message
   character*100   :: FO_trace

!  HELP variables for FO

   real(kind=dp)  :: dnm1
   real(kind=dp)  :: FO_diffgrid(maxlayers)
   real(kind=dp)  :: TRUNCFAC, OMW, TMS, bigdelta, raypf(FO_maxgeoms)

!  Numbers

   REAL(kind=dp), parameter :: zero    = 0.0_dp
   real(kind=dp), parameter :: epsilon = 1.0e-08_dp

!  Other variables
!  ===============

   real(kind=dp) :: asym
   real(kind=dp) :: phasfunc(E_ngeoms)
   real(kind=dp) :: sum1
   real(kind=dp) :: extinc

!  Help variables

   integer        :: l,n,w,m,v
   character*5    :: c5

!  timing variables

   real(kind=sp)  :: e1,e2,e3 

   logical, parameter :: DO_debug_output     = .false.
!   logical, parameter :: DO_debug_output     = .true.

!  START CODE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   write(*,*)
   write(*,*) 'Doing optimized 2S-ESS calculation'
   write(*,*)

!  Initialize times (necessary here)

   ExactTimes = 0.0
   ExactRunTime = 0.0

!  Initialize output

   Intensity_FO_Exact  = zero
   Intensity_2S_Exact  = zero

!  Initial CPU call

   if ( Monitor_CPU ) call cpu_time(e1)

!  1a. 2STREAM and FO initial Setups
!      =============================

   if ( Monitor_CPU ) call cpu_time(e2)

   ngeoms = E_ngeoms
   do_obsgeoms = .true.

!  FO Settings
!  -----------

!  Set Local sphericity for FO

   FO_do_enhanced_ps = .true.
   FO_do_planpar     = .false.

!  Fine layer control. Now fixed input 8/6/13

   FO_nfineinput = 3

!  General flags

   FO_do_obsgeom    = do_obsgeoms !  Use optional input (Mark 6A, 8/6/13)

   FO_do_lambertian     = .true.         ! Lambertian only in this run
   FO_do_deltam_scaling = .true.

!  Toa output, Trig and Flux. You must set this before the "Angles section:

   FO_n_user_levels  = 1
   FO_user_levels(1) = 0

   DNM1 = 5.0_dp

!  Local angles

   sza_boa(:) = 30.0_dp
   vza_boa(:) = 0.0_dp
   azm_boa(:) = 0.0_dp
   FO_ngeoms  = ngeoms
   do n = 1, FO_ngeoms
      FO_obsgeom_boa(n,1) = sza_boa(n)
      FO_obsgeom_boa(n,2) = vza_boa(n)
      FO_obsgeom_boa(n,3) = azm_boa(n)
   enddo
   FO_vsign = +1.0_dp   ! scattering angle control, Upwelling only

!  Attenuation control

   FO_doCrit = .false.        ! Criticality not set (fast option)

!  Copy the earth radius, and height grid

   nlayers = E_nlayers
   FO_nlayers = nlayers
   eradius = 6371.0_dp
   FO_eradius = eradius
   do n = 0, nlayers
      heights(n) = real(nlayers-n,kind=dp)
   enddo
   FO_heights(0:nlayers) = heights(0:nlayers)

!  Find the maximum extinction for criticality (Clouds as layers)

   FO_extincs = 0.0_dp
   do n = 1, FO_nlayers
      taudp(n) = 0.01_dp
      FO_diffgrid(n) = heights(n-1) - heights(n)
      extinc = taudp(n)/FO_diffgrid(n)
      FO_extincs(n) = max(FO_extincs(n),extinc)
   enddo

!  Number of moments

   FO_NMOMENTS = nmoms_all

!  2stream settings
!  ----------------

!  Set BVP Scale factor

   BVPSCALEFACTOR = 1.0d0

!  Taylor control

   TAYLOR_SMALL = 1.0d-3
   TAYLOR_ORDER = 3

!  Set ntotal

   NTOTAL  = 2 * NLAYERS

!  Copy the Geometry settings

   N_USER_OBSGEOMS = ngeoms

   USER_OBSGEOMS(1:ngeoms,1) = sza_boa(1:ngeoms)
   USER_OBSGEOMS(1:ngeoms,2) = vza_boa(1:ngeoms)
   USER_OBSGEOMS(1:ngeoms,3) = azm_boa(1:ngeoms)

!  Earth radius and heights

   HEIGHT_GRID  = heights(0:nlayers)
   EARTH_RADIUS = eradius

!  1/15/16 No longer using sun-normalized automatically
!  Set FO/2S Fluxes within wavelength loop; they depend on input solar spectrum

!  Copy to 2stream control

   DO_UPWELLING = .true.
   DO_DNWELLING = .false.
   DO_BRDF_SURFACE     = .false.
   DO_PLANE_PARALLEL   = .false.

!  Deltam-scaling flag MUST BE SET
!  Otherwise, results are unphysical. 2/5/15

   DO_D2S_SCALING      = .true.

!  Flux flags

   DO_MVOUT_ONLY       = .false.
   DO_ADDITIONAL_MVOUT = .false.

!  5/6/15 upgrade - initialize BRDF inputs

!  Surface BRDF inputs

   BRDF_F_0 = 0.0_dp
   BRDF_F   = 0.0_dp
   UBRDF_F  = 0.0_dp

!  SLEAVE Option Flags

   DO_SURFACE_LEAVING  = .false.
   DO_SL_ISOTROPIC     = .true.

!  SLEAVE Terms

   SLTERM_ISOTROPIC    = 0.0_dp
   SLTERM_F_0          = 0.0_dp

!  LEVELOUT input, exclusive to 2-stream

   DO_2S_LEVELOUT      = .false.

!  Only remaining exclusive 2-stream input

   if ( DO_FULLQUADRATURE ) then
      STREAM_VALUE = 1.0_dp / sqrt(3.0_dp)
   else
      STREAM_VALUE = 0.5_dp
   endif

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(1) = e3 - e2  ! setuptime_1a
   endif

!  1b. First-Order Preliminary Geometry calculations
!  =================================================

!  Call FO geometry routine. New FO code, Version 1.4, Lattice option

   if (monitor_CPU) call cpu_time(e2)

!  Call FO input check routine

   call FO_SSGeometry_inputs &
      ( FO_maxgeoms, FO_ngeoms,           & ! Input dimension, control
        FO_do_planpar, FO_do_enhanced_ps, & ! Input flags
        FO_obsgeom_boa,                   & ! Input
        FO_fail, FO_message, FO_trace )     ! Output(Status)

!  Exception handling for FO inputs

   if ( FO_Fail ) then
      B_FAIL = .true.
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(FO_message) - '//Adjustl(TRIM(FO_message))
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(FO_Trace)   - '//Adjustl(TRIM(FO_trace))
      go to 69
   endif

!  Call appropriate FO geometry routine. New FO code, Version 1.4
      
   if ( FO_do_enhanced_ps ) then
      call FO_SSGeometry_Master_Obs_EPS &
         ( FO_maxgeoms, FO_maxlayers, FO_maxfine,                          & ! Input Dimensions
           FO_ngeoms, FO_nlayers, FO_nfineinput, FO_dtr, FO_Pie, FO_vsign, & ! Input control and constants
           FO_eradius, FO_heights, FO_obsgeom_boa,                         & ! Input geometry/heights
           FO_doCrit, FO_Acrit, FO_extincs,                                & ! Input/Output, Input
           FO_doNadir, FO_Raycon, FO_radii, FO_cota,                       & ! Output(level)
           FO_nfinedivs, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine,       & ! Output(Fine)
           FO_NCrit, FO_RadCrit, FO_CotCrit, FO_Mu0, FO_cosscat,           & ! Output(Crit/scat)
           FO_sunpathsnl, FO_ntraversenl,                                  & ! Output(Sunpaths)
           FO_sunpaths_fine, FO_ntraverse_fine,                            & ! Output(Sunpaths)
           FO_fail, FO_message, FO_trace )                                   ! Output(Status)
   else 
      if ( FO_do_planpar ) then
         call FO_SSGeometry_Master_Obs_PP &
            ( FO_maxgeoms,                 & ! Input dimensions
              FO_ngeoms, FO_dtr, FO_vsign, & ! Input control and constants
              FO_obsgeom_boa,              & ! Input
              FO_Mu0, FO_Mu1, FO_cosscat )   ! Output(scat)
      else
         call FO_SSGeometry_Master_Obs_RPS &
            ( FO_maxgeoms, FO_maxlayers,               & ! Input dimensions
              FO_ngeoms, FO_nlayers, FO_dtr, FO_vsign, & ! Input control and constants
              FO_eradius, FO_heights, FO_obsgeom_boa,  & ! Input
              FO_Mu0, FO_Mu1, FO_cosscat,              & ! Output(scat)
              FO_sunpaths )                              ! Output(Sunpaths)
      endif
   endif

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(2) = e3 - e2 ! FOgeomtime_1b
   endif

!  Exception handling for FO geometry

   if ( FO_Fail ) then
      B_FAIL = .true.
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(FO_message) - '//Adjustl(TRIM(FO_message))
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(FO_Trace)   - '//Adjustl(TRIM(FO_trace))
      go to 69
   endif

!  1c. First-Order Spherical Function calculations
!  ===============================================

!  Get the Legendre polynomials

   call cpu_time(e2)

   call FO_ScalarSS_spherfuncs &
      ( FO_MAXMOMENTS, FO_MAXGEOMS, FO_NMOMENTS, FO_NGEOMS, & ! Inputs
        FO_COSSCAT, FO_GENSPHER )                             ! Outputs

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(3) = e3 - e2 ! spherfntime_1c
   endif

!  Compute phase function for all geometries
   
   asym = 0.9_dp
   do V = 1, FO_ngeoms
      sum1 = zero
      do L = 0, FO_nmoments
         sum1 = sum1 + FO_GENSPHER(L,V) * (2.0_dp * real(l,kind=dp) + 1.0_dp) * asym**l
      enddo
      phasfunc(v) = sum1
   enddo

!  1d. 2S Preliminary Geometry/WL independent calculations
!  =======================================================

!  Call 2S geometry/input check routine. Optimized 2S code

   if (monitor_CPU) call cpu_time(e2)

   call TWOSTREAM_INPUTS_SOLAR_OBS &
      ( MAXLAYERS, MAXTOTAL, MAX_USER_OBSGEOMS, MAX_MESSAGES,           & ! Dimensions
        DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL,                  & ! Inputs
        DO_MVOUT_ONLY, DO_ADDITIONAL_MVOUT,                             & ! Inputs
        NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS,  & ! Inputs
        EARTH_RADIUS, HEIGHT_GRID,                                      & ! Inputs
        MAXSTREAMS_BRDF, NSTREAMS_BRDF, DO_BRDF_SURFACE,                & ! BRDF Inputs
        DO_SURFACE_LEAVING, DO_FLUORESCENCE, FL_DO_DataGaussian, DO_ISOTROPIC, & ! SLEAVE Inputs
        FL_Latitude, FL_Longitude, FL_Epoch, FL_InputGAUSSIANS,         & ! SLEAVE Inputs
        N_FOURIERS, DEG_TO_RAD, PI4, DO_POSTPROCESSING, DO_INCLUDE_MVOUT, & ! Outputs
        AVERAGE_SECANT_PP, CHAPMAN_FACTORS, X0,                         & ! Outputs
        USER_STREAMS, USER_SECANTS, AZMFAC,                             & ! Outputs
        PX11, PXSQ, PX0X, ULP,                                          & ! Outputs
        SURFACE_FACTOR, DELTA_FACTOR,                                   & ! Outputs
        PIE, STREAM_SINE, SX0, USER_SINES, NBRDF_HALF,                  & ! BRDF Outputs
        X_BRDF, CX_BRDF, SX_BRDF, A_BRDF,                               & ! BRDF Outputs
        FL_GAUSSIANS, Fs755,                                            & ! SLEAVE Outputs
        STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS )           ! Exception handling

   if ( Monitor_CPU ) then
      call cpu_time(e3) ; Exacttimes(4) = e3 - e2 ! 2Sgeomtime_1d
   endif

!  Exception handling for 2S geometry/input check

   if ( STATUS_INPUTCHECK .NE. 0 ) then
      B_FAIL = .true.
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(2S_message) - '//Adjustl(TRIM(C_messages(0)))
      B_NMESSAGES = B_NMESSAGES + 1
      B_MESSAGES(B_NMESSAGES) = '(2S_Trace)   - '//Adjustl(TRIM(C_actions(0)))
      go to 69
   endif

!  2. MAIN SET OF CALCULATIONS (call RT code 100,000 times)
!     ========================

   do w = 1, 100000

!  2a. Compute layer input optical properties for 2stream and First Order
!  ----------------------------------------------------------------------

      if ( Monitor_CPU ) call cpu_time(e2)

!  Preliminary zeroing

      FO_Exactscat = zero
      FO_deltaus   = zero
      ASYMM_INPUT  = zero ; D2S_SCALING = zero
      deltau_input = zero ; omega_input = zero

      if (monitor_CPU) call cpu_time(e2)

      albedo = 0.3_dp

!  First Order optical properties (New Style)

      do n = 1, FO_nlayers
         FO_deltaus(n) = taudp(n)
         FO_extincs(n) = FO_deltaus(n) / FO_diffgrid(n)
         omega(n) = 0.5_dp
      enddo
      do v = 1, ngeoms
         do n = 1, FO_nlayers
            FO_exactscat(n,v) = phasfunc(v)
            omw = omega(n)
            if ( FO_do_deltam_scaling ) then
               truncfac = asym**2
               tms = omw / (1.0_dp - truncfac * omw)
            else
               tms = omw
            endif
            FO_exactscat(n,v) = FO_exactscat(n,v) * tms
         enddo
      enddo
      FO_reflec(1:FO_ngeoms) = albedo

!  1/15/16. Solar Flux. Sun-normalized.

      FO_FLUX  = 0.25_dp / FO_Pie

!  2-stream optical properties (Use FO results where possible)

      do n = 1, nlayers
         deltau_input(n) = FO_deltaus(n)
         omega_input(n)  = omega(n)
         ASYMM_INPUT(n) = asym
         D2S_SCALING(n) = asym**2
      enddo
      LAMBERTIAN_ALBEDO = albedo

!  1/15/16. Solar Flux. Sun-normalized.

      FLUX_FACTOR  = 1.0_dp

!  time

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(5) = Exacttimes(5) + e3 - e2  !  ExactOpTime_2a_FO2S
      endif

!  2c. 2STREAM call
!  ----------------

!  5/6/15 upgrade  - modified control structure, and added timing
!  5/6/2015 Upgrade - added new 2S call, Version 2.4

      if (monitor_CPU) call cpu_time(e2)

      CALL TWOSTREAM_MASTER_SOLAR_OBS &
         ( MAXLAYERS, MAXTOTAL, MAX_MESSAGES, MAX_USER_OBSGEOMS,          & ! Dimensions
           DO_UPWELLING, DO_DNWELLING, DO_PLANE_PARALLEL, DO_2S_LEVELOUT, & ! Inputs     !@@ 2p2
           DO_MVOUT_ONLY, DO_INCLUDE_MVOUT,                               & ! Inputs     !@@ 2p3
           DO_POSTPROCESSING, DO_D2S_SCALING, DO_BRDF_SURFACE,            & ! Inputs     !@@ 2p1
           DO_SURFACE_LEAVING, DO_SL_ISOTROPIC, DO_PDINVERSE,             & ! Input !@@ 2p3 6/25/14
           BVPSCALEFACTOR, TAYLOR_ORDER, TAYLOR_SMALL,                    & ! Input !@@ 2p3 6/25/14, 8/15/14
           NLAYERS, NTOTAL, STREAM_VALUE, N_USER_OBSGEOMS, USER_OBSGEOMS, & ! Inputs     !@@ 2p1
           FLUX_FACTOR, SURFACE_FACTOR, DELTA_FACTOR,                     & ! Inputs
           DELTAU_INPUT, OMEGA_INPUT, ASYMM_INPUT, D2S_SCALING,           & ! Inputs
           LAMBERTIAN_ALBEDO, BRDF_F_0, BRDF_F, UBRDF_F,                  & ! Inputs
           SLTERM_ISOTROPIC, SLTERM_F_0,                                  & ! Inputs  !@@ 2p3 (Sleave)
           N_FOURIERS, PI4, AVERAGE_SECANT_PP,                            & ! Inputs (geometry)
           CHAPMAN_FACTORS, X0, USER_STREAMS, USER_SECANTS, AZMFAC,       & ! Inputs (geometry)
           PX11, PXSQ, PX0X, ULP,                                         & ! Inputs
           INTENSITY_TOA, INTENSITY_BOA,                                  & ! Outputs
           RADLEVEL_UP, RADLEVEL_DN,                                      & ! Outputs !@@ 2p2
           FLUXES_TOA_SOLAR, FLUXES_BOA_SOLAR,                            & ! Outputs !@@ 2p3 (Fluxes)
           STATUS_INPUTCHECK, C_NMESSAGES, C_MESSAGES, C_ACTIONS,         & ! Exception handling
           STATUS_EXECUTION,  E_MESSAGE, E_TRACE_1, E_TRACE_2 )             ! Exception handling

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(6) = Exacttimes(6) + e3 - e2 ! Exact2Stime_2c 
      endif

!  Exception handling for Input checks. 5/6/15 upgrade, added wavelength number

      if ( STATUS_INPUTCHECK .eq. 1 ) then
         DO M = 1, C_NMESSAGES
            B_MESSAGES(B_NMESSAGES+2*M-1) = '(2STREAM_2p4 message) '//Adjustl(TRIM(C_MESSAGES(M)))
            B_MESSAGES(B_NMESSAGES+2*M)   = '(2STREAM_2p4 action ) '//Adjustl(TRIM(C_ACTIONS(M)))
         ENDDO
         B_NMESSAGES = B_NMESSAGES + 2*C_NMESSAGES
         write(C5,'(I5)')w ; B_MESSAGES(B_NMESSAGES+1) = '(Optional 2S Input  Check, RT Exact # = '//C5//')'
         B_NMESSAGES = B_NMESSAGES + 1 ; B_fail = .true.;  go to 69
      endif

!  Exception handling for Calculation. 5/6/15 upgrade, added wavelength number

      IF ( STATUS_EXECUTION .eq. 1 ) THEN
         B_MESSAGES(B_NMESSAGES+1)   = '(2STREAM_2p4 message) '//Adjustl(TRIM(E_MESSAGE))
         B_MESSAGES(B_NMESSAGES+2)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(E_TRACE_1))
         B_MESSAGES(B_NMESSAGES+3)   = '(2STREAM_2p4 trace  ) '//Adjustl(TRIM(E_TRACE_2))
         B_NMESSAGES = B_NMESSAGES + 3
         write(C5,'(I5)')w ; B_MESSAGES(B_NMESSAGES+1) = '(Optional 2S Execution, RT Exact # = '//C5//')'
         B_NMESSAGES = B_NMESSAGES + 1 ; B_fail = .true.;  go to 69
      endif

!  Save exact 2stream results

      do v = 1, ngeoms
         intensity_2S_Exact(w,v) = intensity_TOA(v)
      enddo

!  2d. first order calculation
!  ---------------------------

      if (monitor_CPU) call cpu_time(e2)

!  Call

      if ( FO_do_enhanced_ps ) then
         call SS_Integral_I_UP_EPS &
           ( FO_maxgeoms, FO_maxlayers, FO_maxfine, FO_maxuserlevels, FO_doNadir,     & ! Inputs (dimension)
             FO_ngeoms, FO_nlayers, FO_nfinedivs, FO_n_user_levels, FO_user_levels,   & ! Inputs (control)  
             FO_reflec, FO_extincs, FO_deltaus, FO_exactscat, FO_flux,                & ! Inputs (Optical)
             FO_Mu0, FO_NCrit, FO_xfine, FO_wfine, FO_csqfine, FO_cotfine, FO_Raycon, & ! Inputs (Geometry)
             FO_cota, FO_sunpathsnl, FO_ntraversenl,                                  & ! Inputs (Geometry)
             FO_sunpaths_fine, FO_ntraverse_fine,                                     & ! Inputs (Geometry)
             FO_intensity_ss, FO_intensity_db )                                         ! Outputs
      else
         if ( FO_do_planpar ) then
            call SS_Integral_I_UP_PP &
               ( FO_maxgeoms, FO_maxlayers, FO_maxuserlevels,             & ! Inputs (dimensioning)
                 FO_ngeoms, FO_nlayers, FO_n_user_levels, FO_user_levels, & ! Inputs (control output)
                 FO_reflec, FO_deltaus, FO_exactscat, FO_flux,            & ! Inputs (Optical)
                 FO_Mu0, FO_Mu1,                                          & ! Inputs (Geometry) 
                 FO_intensity_ss, FO_intensity_db )                         ! Outputs
         else
            call SS_Integral_I_UP_RPS &
               ( FO_maxgeoms, FO_maxlayers, FO_maxuserlevels,              & ! Inputs (dimensioning)
                 FO_ngeoms, FO_nlayers, FO_n_user_levels, FO_user_levels,  & ! Inputs (control output)
                 FO_reflec, FO_extincs, FO_deltaus, FO_exactscat, FO_flux, & ! Inputs (Optical)
                 FO_Mu0, FO_Mu1, FO_sunpaths,                              & ! Inputs (Geometry)
                 FO_intensity_ss, FO_intensity_db )                          ! Outputs
         endif
      endif

      if (monitor_CPU) then
         call cpu_time(e3) ; Exacttimes(7) = Exacttimes(7) + e3 - e2 ! ExactFOtime_2d
      endif

!  Save FO results

      do v = 1, ngeoms
         intensity_FO_Exact(w,v) = FO_intensity_ss(1,v) + FO_intensity_db(1,v)
      enddo

!  4. Debug
!  --------

!  FD Testing: Debug First Order output. FORT.75/76/77
!
!      if ( w.lt.11 .and.do_debug_output ) then
!         write(*,202)w,intensity_FO_Exact(w,1),FO_intensity_ss(1,1), FO_intensity_db(1,1)
!      endif
!202   format(i5,1p3e20.10)

! Debug output to FORT.810

      if ( do_debug_output ) then
         v = 1
         write(810,250)w,&
             intensity_2S_Exact(w,v),intensity_FO_Exact(w,v)
250      format(i5,1p2e20.10)
      endif
!      STOP '810 first'

!  End of main wavelength loop

   enddo

!  Final CPU call

   if (monitor_CPU) then
      call cpu_time(e3) ; ExactRuntime = ExactRuntime + e3-e1
   endif

write(0,*) 'Initial setup time = ', Exacttimes(1)
write(0,*) 'ESS geometry calc time = ', Exacttimes(2)
write(0,*) 'Spherical function calc time = ', Exacttimes(3)
write(0,*) '2S geometry calc time = ', Exacttimes(4)
write(0,*) 'Optical property calc time = ', Exacttimes(5)
write(0,*) '2S calc time = ', Exacttimes(6)
write(0,*) 'ESS calc time = ', Exacttimes(7)
write(0,*) 'Total run time = ', ExactRuntime
write(0,*) '2S intensity = ', intensity_2S_Exact(1,1)
write(0,*) 'ESS intensity = ', intensity_FO_Exact(1,1)


   ! output data into a file 
open(1, file = 'test_opt_output.dat', status = 'old')   
   write(1,*) '2S intensity = ', intensity_2S_Exact(1,1)
   write(1,*) 'ESS intensity = ', intensity_FO_Exact(1,1) 
close(1) 

open(2, file = 'test_opt_timings.dat', status = 'old')     
   write(2,*) 'Initial setup time = ', Exacttimes(1)
   write(2,*) 'ESS geometry calc time = ', Exacttimes(2)
   write(2,*) 'Spherical function calc time = ', Exacttimes(3)
   write(2,*) 'Optical property calc time = ', Exacttimes(4)
   write(2,*) '2S calc time = ', Exacttimes(5)
   write(2,*) 'ESS calc time = ', Exacttimes(6)
   write(2,*) 'Total run time = ', ExactRuntime 
close(2) 

!  Normal Finish

   STOP

!  Error Finish

69 continue

   write(0,*) 'Error in running'

!  End program

end program test_opt
