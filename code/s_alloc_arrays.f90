SUBROUTINE alloc_arrays

USE constants
USE quantitys

IMPLICIT NONE



!----------------------------------------------------------------------------------------------------



ALLOCATE(     &
              a_g             (  0:N+1 )        , &

              eta             ( 0:N+1 )         , &
              eta_tilde       ( 0:N   )         , &
              zeta            ( 0:N+1 )         , &
              zeta_tilde      ( 0:N   )         , &

              d_radius        (  0:N+1 )        , &
              d_radius_tilde  ( -1:N   )        , &
              radius          ( -1:N+1 )        , &
              radius_tilde    ( -2:N+1 )        , &
              ring_area       (  0:N+1 )        , &

              d_u             ( 1:3, 0:N+1 )    , &
              u               ( 1:3, 0:N+1 )    , &
              u_dev           ( 1:3, 0:N+1 )    , &
              u_old           ( 1:3, 0:N+1 )    , &

              u_tilde         ( 1:3, 0:N   )    , &
              u_tilde_m       ( 1:3, 0:N   )    , &
              u_tilde_p       ( 1:3, 0:N   )    , &

              A               ( 1:k, 1:N )      , &
              D               ( 1:k, 1:N )      , &
              D_old           ( 1:k, 1:N )      , &

              F_num           ( 1:3, 0:N )      , &
              V               ( 1:3, 0:N )      , &
              W               ( 1:3, 1:N )      , &
              S               ( 1:3, 1:N )      , &

              G               ( 0:N+1, 0:N+1 )  , &
              H               ( 1:k*N, 1:k*N )  , &

              d2              ( 1:3, 0:N )     ,  &
              f               ( 1:N )                )



!----------------------------------------------------------------------------------------------------

END SUBROUTINE alloc_arrays
