module twostream_external_calc_m

    Use twostream_miscsetups_m

public

contains

subroutine twostream_external_calc ( &
    MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,    & ! Dimensions
    NLAYERS, NBEAMS, N_USER_STREAMS, USER_ANGLES, & ! Input
    DO_PLANE_PARALLEL,                        & ! Input
    BEAM_SZAS, EARTH_RADIUS, HEIGHT_GRID,     & ! Input
    STREAM_VALUE,                             & ! INPUT
    CHAPMAN_FACTORS, LOCAL_SZA,               & ! OUTPUT
    SinStream, X0,                            & ! OUTPUT  
    USER_STREAMS, USER_SECANTS                & ! OUTPUT
    )
    !  precision and parameters

    implicit none

    INTEGER      , PARAMETER :: dp   = KIND( 1.0D0 )
    REAL(kind=dp), parameter :: zero = 0.0_dp, one = 1.0_dp

    
    INTEGER, INTENT(IN)        :: MAXLAYERS
    INTEGER, INTENT(IN)        :: MAXBEAMS
    INTEGER, INTENT(IN)        :: MAX_USER_STREAMS
    INTEGER, INTENT(IN)        :: NLAYERS
    INTEGER, INTENT(IN)        :: NBEAMS
    INTEGER, INTENT(IN)        :: N_USER_STREAMS
    REAL(kind=dp), INTENT(IN)  :: USER_ANGLES  ( MAX_USER_STREAMS )
    
    LOGICAL, INTENT(IN)        :: DO_PLANE_PARALLEL
    REAL(kind=dp), INTENT(IN)  :: BEAM_SZAS ( MAXBEAMS )
    REAL(kind=dp), INTENT(IN)  :: EARTH_RADIUS
    REAL(kind=dp), INTENT(IN)  :: HEIGHT_GRID ( 0:MAXLAYERS )
    REAL(kind=dp), INTENT(IN)  :: STREAM_VALUE

    REAL(kind=dp), intent(out) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
    REAL(kind=dp), intent(out) :: LOCAL_SZA       ( 0:MAXLAYERS, MAXBEAMS )
    REAL(kind=dp), intent(out) :: SINSTREAM
    REAL(kind=dp), intent(out)    :: X0  ( MAXBEAMS )
    REAL(kind=dp), intent(out)    :: USER_STREAMS ( MAX_USER_STREAMS )
    REAL(kind=dp), intent(out)    :: USER_SECANTS ( MAX_USER_STREAMS )


    !Local
    Integer             :: IB, I
    real(kind=dp)       :: deg_to_rad, PI4
    DEG_TO_RAD = ACOS( - one ) / 180.0_dp
    PI4 = DEG_TO_RAD * 720.0_dp    

    DO IB = 1, NBEAMS
       CALL TWOSTREAM_BEAM_GEOMETRY_PREPARE &
          ( MAXLAYERS, MAXBEAMS,                      & ! Dimensions
            NLAYERS, DO_PLANE_PARALLEL, IB,           & ! Input
            BEAM_SZAS(IB), EARTH_RADIUS, HEIGHT_GRID, & ! Input
            CHAPMAN_FACTORS, LOCAL_SZA )                ! In/Out
    ENDDO

!  Get derived inputs
!  ==================

!  Quadrature

!     MUSTREAM  = STREAM_VALUE
    SINSTREAM = SQRT ( ONE - STREAM_VALUE * STREAM_VALUE )

!  Solar zenith angle cosine

    DO IB = 1, NBEAMS
       X0(IB) = COS ( BEAM_SZAS(IB) * DEG_TO_RAD )
    ENDDO

!  User stream cosines. 11/5/13 2p3 Post-processing control

    DO I = 1, N_USER_STREAMS
          USER_STREAMS(I) = COS(DEG_TO_RAD * USER_ANGLES(I))
          USER_SECANTS(I) = ONE / USER_STREAMS(I)
    ENDDO

    end subroutine twostream_external_calc

end module twostream_external_calc_m