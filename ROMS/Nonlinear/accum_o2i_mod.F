#include "cppdefs.h"
      MODULE accum_o2i_mod
#if defined CICE_OCEAN
      use mod_kinds, only: r8

! Accumulate the ROMS-to-CICE coupling fields over the given coupling timestep.
! The frazil array qfraz_accum is handled separately, in frazil_ice_prod.

      implicit none

      PRIVATE
      PUBLIC accum_o2i

      CONTAINS

      SUBROUTINE accum_o2i (ng, tile)

      USE mod_param
      USE mod_ocean
      USE mod_ice

      integer, intent(in) :: ng, tile
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 44)
# endif
!
      CALL accum_o2i_tile (ng, tile, LBi, UBi, LBj, UBj,                &
     &                     OCEAN(ng) % t, OCEAN(ng) % u, OCEAN(ng) % v, &
     &                     OCEAN(ng) % zeta, ICE(ng) % sst_accum,       &
     &                     ICE(ng) % sss_accum, ICE(ng) % u_accum,      &
     &                     ICE(ng) % v_accum, ICE(ng) % ssh_accum)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 44)
# endif
      RETURN
      END SUBROUTINE accum_o2i

      subroutine accum_o2i_tile (ng, tile, LBi, UBi, LBj, UBj,          &
     &                           t, u, v, zeta, sst_accum, sss_accum,   &
     &                           u_accum, v_accum, ssh_accum)

      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE exchange_2d_mod
      USE exchange_3d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: sst_accum(LBi:,LBj:)
      real(r8), intent(inout) :: sss_accum(LBi:,LBj:)
      real(r8), intent(inout) :: u_accum(LBi:,LBj:,:)
      real(r8), intent(inout) :: v_accum(LBi:,LBj:,:)
      real(r8), intent(inout) :: ssh_accum(LBi:,LBj:)
# else
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: sst_accum(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: sss_accum(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: u_accum(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: v_accum(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ssh_accum(LBi:UBi,LBj:UBj)
# endif
!
! Local variable definitions
!
      integer :: i, j, k

# include "set_bounds.h"

      DO j=JstrR,JendR
        DO i=IstrR,IendR
          sst_accum(i,j)=sst_accum(i,j)+t(i,j,N(ng),NOUT,itemp)
          sss_accum(i,j)=sss_accum(i,j)+t(i,j,N(ng),NOUT,isalt)
          ssh_accum(i,j)=ssh_accum(i,j)+zeta(i,j,NOUT)
          DO k=1,N(ng)
            u_accum(i,j,k)=u_accum(i,j,k)+0.5_r8*                       &
     &                     (u(i+1,j,k,NOUT)+u(i,j,k,NOUT))
            v_accum(i,j,k)=v_accum(i,j,k)+0.5_r8*                       &
     &                     (v(i,j+1,k,NOUT)+v(i,j+1,k,NOUT))
          END DO
        END DO
      END DO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          sst_accum)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 1, LBi, UBi, LBj, UBj,      &
     &                      NghostPoints, EWperiodic(ng),               &
     &                      NSperiodic(ng), sst_accum)
# endif
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          sss_accum)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 1, LBi, UBi, LBj, UBj,      &
     &                      NghostPoints, EWperiodic(ng),               &
     &                      NSperiodic(ng), sss_accum)
# endif
        CALL exchange_r2d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          ssh_accum)
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, iNLM, 1, LBi, UBi, LBj, UBj,      &
     &                      NghostPoints, EWperiodic(ng),               &
     &                      NSperiodic(ng), ssh_accum)
# endif
        CALL exchange_r3d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          1, N(ng), u_accum)
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 1, LBi, UBi, LBj, UBj,      &
     &                      1, N(ng), NghostPoints, EWperiodic(ng),     &
     &                      NSperiodic(ng), u_accum)
# endif
        CALL exchange_r3d_tile (ng, tile, LBi, UBi, LBj, UBj,           &
     &                          1, N(ng), v_accum)
# ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 1, LBi, UBi, LBj, UBj,      &
     &                      1, N(ng), NghostPoints, EWperiodic(ng),     &
     &                      NSperiodic(ng), v_accum)
# endif
      END IF

      RETURN
      END SUBROUTINE accum_o2i_tile

#endif
      END MODULE accum_o2i_mod

      
      
      
      

