SUBROUTINE dealloc_arrays

USE constants
USE quantitys

IMPLICIT NONE

!----------------------------------------------------------------------------------------------------



DEALLOCATE(   &
              a_g             , &

              eta             , &
              eta_tilde       , &
              zeta            , &
              zeta_tilde      , &

              d_radius        , &
              d_radius_tilde  , &
              radius          , &
              radius_tilde    , &
              ring_area       , &

              d_u             , &
              u               , &
              u_dev           , &
              u_old           , &

              u_tilde         , &
              u_tilde_m       , &
              u_tilde_p       , &

              A               , &
              D               , &
              D_old           , &

              F_num           , &
              V               , &
              W               , &
              S               , &

              G               , &
              H               , &

              d2              ,  &
              f                    )



!----------------------------------------------------------------------------------------------------

END SUBROUTINE dealloc_arrays
