#include "cppdefs.h"
      MODULE set_vbc_mod
#ifdef NONLINEAR
!
!svn $Id: set_vbc.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: set_vbc
!
      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_ice
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % Hz,                                 &
#  if defined UV_LOGDRAG
     &                   GRID(ng) % ZoBot,                              &
#  elif defined UV_LDRAG
     &                   GRID(ng) % rdrag,                              &
#  elif defined UV_QDRAG
     &                   GRID(ng) % rdrag2,                             &
#  endif
#  if !defined BBL_MODEL
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
#  endif
     &                   OCEAN(ng) % t,                                 &
#  if !defined BBL_MODEL
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
#  endif
#  ifdef QCORRECTION
     &                   FORCES(ng) % dqdt,                             &
     &                   FORCES(ng) % sst,                              &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                   FORCES(ng) % sss,                              &
#   if defined SAVE_SSFREST
     &                   FORCES(ng) % ssflx_rest,                       &
#   endif
#  endif
#  ifdef SCORRECTION_MASK
     &                   GRID(ng)%latr, GRID(ng)%h, GRID(ng)%zice,      &
#  endif
#  ifdef SSFLUX_EXTRA
     &                   FORCES(ng) % ssflux_extra,                     &
#  endif
#  ifdef ICEBERGS
     &                   FORCES(ng) % icebergs,                         &
#  endif
#  if defined CICE_COUPLING
#   if defined SHORTWAVE 
     &                   FORCES(ng) % srflx,                            &
#   endif
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
     &                     ICE(ng) % aice,                              &
     &                     ICE(ng) % aice_u,                            &
     &                     ICE(ng) % aice_v,                            &
     &                     ICE(ng) % fhocnAI,                           &
     &                     ICE(ng) % fswthruAI,                         &
     &                     ICE(ng) % freshAI,                           &
     &                     ICE(ng) % fsaltAI,                           &
     &                     ICE(ng) % stru,                              &
     &                     ICE(ng) % strv,                              &
# endif
#  ifndef BBL_MODEL
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
#  endif
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % btflx)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif

      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         Hz,                                      &
#  if defined UV_LOGDRAG
     &                         ZoBot,                                   &
#  elif defined UV_LDRAG
     &                         rdrag,                                   &
#  elif defined UV_QDRAG
     &                         rdrag2,                                  &
#  endif
#  if !defined BBL_MODEL
     &                         z_r, z_w,                                &
#  endif
     &                         t,                                       &
#  if !defined BBL_MODEL
     &                         u, v,                                    &
#  endif
#  ifdef QCORRECTION
     &                         dqdt, sst,                               &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                         sss,                                     &
#   if defined SAVE_SSFREST
     &                         ssflx_rest,                              &
#   endif
#  endif
#  ifdef SCORRECTION_MASK
     &                         latr, h,  zice,                          &
#  endif
#  ifdef SSFLUX_EXTRA
     &                         ssflux_extra,                            &
#  endif
#  ifdef ICEBERGS
     &                         icebergs,                                &
#  endif
#  if defined CICE_COUPLING
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
     &                         sustr, svstr,                            &
     &                         aice,                                    &
     &                         aice_u,                                  &
     &                         aice_v,                                  &
     &                         fhocnAI,                                 &
     &                         fswthruAI,                               &
     &                         freshAI,                                 &
     &                         fsaltAI,                                 &
     &                         stru,                                    &
     &                         strv,                                    &
#  endif
#  ifndef BBL_MODEL
     &                         bustr, bvstr,                            &
#  endif
     &                         stflx, btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
#ifdef CICE_COUPLING
      USE exchange_2d_mod
#endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
#   if !defined BBL_MODEL
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   if !defined BBL_MODEL
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#    ifdef SAVE_SSFREST
      real(r8), intent(inout) :: ssflx_rest(LBi:,LBj:)
#    endif
#   endif
#   ifdef SCORRECTION_MASK
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
#   ifdef SSFLUX_EXTRA
      real(r8), intent(in) :: ssflux_extra(LBi:,LBj:)
#   endif
#   ifdef ICEBERGS
      real(r8), intent(in) :: icebergs(LBi:,LBj:)
#   endif
#   if defined CICE_COUPLING
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: aice(LBi:,LBj:)
      real(r8), intent(in) :: aice_u(LBi:,LBj:)
      real(r8), intent(in) :: aice_v(LBi:,LBj:)
      real(r8), intent(in) :: fhocnAI(LBi:,LBj:)
      real(r8), intent(in) :: fswthruAI(LBi:,LBj:)
      real(r8), intent(in) :: freshAI(LBi:,LBj:)
      real(r8), intent(in) :: fsaltAI(LBi:,LBj:)
      real(r8), intent(in) :: stru(LBi:,LBj:)
      real(r8), intent(in) :: strv(LBi:,LBj:)
#  endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   if defined UV_LOGDRAG
      real(r8), intent(in) :: ZoBot(LBi:UBi,LBj:UBj)
#   elif defined UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   elif defined UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
#   if !defined BBL_MODEL
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if !defined BBL_MODEL
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#    ifdef SAVE_SSFREST
      real(r8), intent(inout) :: ssflx_rest(LBi:UBi,LBj:UBj)
#    endif
#   endif
#   ifdef SCORRECTION_MASK
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SSFLUX_EXTRA
      real(r8), intent(in) :: ssflux_extra(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ICEBERGS
      real(r8), intent(in) :: icebergs(LBi:UBi,LBj:UBj)
#   endif
#   if defined CICE_COUPLING
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: aice(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: aice_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: aice_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: fhocnAI(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: fswthruAI(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: freshAI(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: fsaltAI(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: stru(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: strv(LBi:UBi,LBj:UBj)
#  endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
!
!  Local variable declarations.
!
      integer :: i, j, itrc

#  if !defined BBL_MODEL || defined CICE_COUPLING || defined LIMIT_STFLX_COOLING
      real(r8) :: cff, cff1, cff2, cff3
#  endif

#  if !defined BBL_MODEL && defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
#  endif

#  if defined SCORRECTION || defined SRELAXATION
      real(r8) :: flag
#  endif

#  ifdef SCORRECTION_MASK
      real(r8), parameter :: lat0 = -60.0_r8
      real(r8), parameter :: h0 = 1500.0_r8
#  endif

#  include "set_bounds.h"

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Add in flux correction to surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
! Add in net heat flux correction.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,itemp)=stflx(i,j,itemp)+                            &
     &                     dqdt(i,j)*(t(i,j,N(ng),nrhs,itemp)-sst(i,j))
        END DO
      END DO
#  endif

#  ifdef LIMIT_STFLX_COOLING
!
!-----------------------------------------------------------------------
!  If net heat flux is cooling and SST is at freezing point or below
!  then suppress further cooling. Note: stflx sign convention is that
!  positive means heating the ocean (J Wilkin).
!-----------------------------------------------------------------------
!
!  Below the surface heat flux stflx(:,:,itemp) is ZERO if cooling AND
!  the SST is cooler that the threshold.  The value is retained if
!  warming.
!
!    cff3 = 0      if SST warmer than threshold (cff1) - change nothing
!    cff3 = 1      if SST colder than threshold (cff1)
!
!    0.5*(cff2-ABS(cff2)) = 0                        if flux is warming
!                         = stflx(:,:,itemp)         if flux is cooling
!
      cff1=-2.0_r8              ! nominal SST threshold to cease cooling
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          cff2=stflx(i,j,itemp)
          cff3=0.5_r8*(1.0_r8+SIGN(1.0_r8,cff1-t(i,j,N(ng),nrhs,itemp)))
          stflx(i,j,itemp)=cff2-cff3*0.5_r8*(cff2-ABS(cff2))
        END DO
      END DO
#  endif

#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Multiply fresh water flux with surface salinity.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        END DO
      END DO
#  endif

#  ifdef CICE_COUPLING
!-----------------------------------------------------------------------
!  Modify heat and salt fluxes fluxes due to sea ice 
!-----------------------------------------------------------------------
      cff1=1.0_r8/(rho0*Cp)
      cff=1.0_r8/rhow

      DO j=JstrR,JendR
        DO i=IstrR,IendR
! merge with fluxes from ice model
          srflx(i,j)=srflx(i,j)*(1.0_r8-aice(i,j))+fswthruAI(i,j)*cff1
          stflx(i,j,itemp)=stflx(i,j,itemp)*(1.0_r8-aice(i,j))+         &
     &         (fswthruAI(i,j)+fhocnAI(i,j))*cff1

! Calculation of virtual salt flux 
          stflx(i,j,isalt) = stflx(i,j,isalt)*(1.0_r8-aice(i,j))        &
     &         -cff*( freshAI(i,j)*t(i,j,N(ng),nrhs,isalt)              &
! Fsalt is given as kg/s/m^2. Should be converted to salinity (factor 1000).
     &         -fsaltAI(i,j)*1000.0_r8)  

          ! merge with stress from ice model
          sustr(i,j)=sustr(i,j)*(1.0_r8-aice_u(i,j))+stru(i,j)*cff
          svstr(i,j)=svstr(i,j)*(1.0_r8-aice_v(i,j))+strv(i,j)*cff

        END DO
      END DO

        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)

      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx, stflx(:,:,itemp),stflx(:,:,isalt))
#   endif
#  endif

#  ifdef ICEBERGS
! Add extra freshwater flux from iceberg melt.
      cff=1.0_r8/rhow
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)-                            &
     &                     cff*t(i,j,N(ng),nrhs,isalt)*icebergs(i,j)
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,isalt))
#   endif
#  endif

#  ifdef SSFLUX_EXTRA
! Add extra surface salt flux.
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)+ssflux_extra(i,j)
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,isalt))
#   endif
#  endif

#  if defined SALINITY && (defined SCORRECTION || defined SRELAXATION)
!
!-----------------------------------------------------------------------
!  Apply correction.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#   ifdef SCORRECTION_MASK
          ! Don't apply correction if lat<lat0 and h<h0, or in ice shelf
          ! cavities. In these cases, flag=0; otherwise, flag=1
          flag=(1.0_r8-MIN(SIGN(1.0_r8,latr(i,j)-lat0),0.0_r8)*         &
     &                 MIN(SIGN(1.0_r8,h(i,j)-h0),0.0_r8))*             &
     &         (1.0_r8+MIN(SIGN(1.0_r8,zice(i,j)),0.0_r8))
#   else
          flag=1.0_r8
#   endif
#   if defined SCORRECTION
          stflx(i,j,isalt)=stflx(i,j,isalt)-                            &
     &                     flag*SSnudg(ng)*Hz(i,j,N(ng))*               &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#   elif defined SRELAXATION
          stflx(i,j,isalt)=-flag*SSnudg(ng)*Hz(i,j,N(ng))*              &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#   endif
#   ifdef SAVE_SSFREST
          ssflx_rest(i,j)=-flag*SSnudg(ng)*Hz(i,j,N(ng))*               &
     &                    (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#   endif
        END DO
      END DO
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,isalt))
#   endif
#  endif

#  ifndef BBL_MODEL
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------

#   ifdef LIMIT_BSTRESS
!
!  Set limiting factor for bottom stress. The bottom stress is adjusted
!  to not change the direction of momentum.  It only should slow down
!  to zero.  The value of 0.75 is arbitrary limitation assigment.
!
      cff=0.75_r8/dt(ng)
#   endif

#   if defined UV_LOGDRAG
!
!  Set logarithmic bottom stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_r(i,j,1)-z_w(i,j,0))/ZoBot(i,j))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*                      &
     &               u(i,j,1,nrhs)*cff2
#    ifdef LIMIT_BSTRESS
          cff3=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
             
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          bvstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*                      &
     &               v(i,j,1,nrhs)*cff2
#    ifdef LIMIT_BSTRESS
          cff3=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
     &               u(i,j,1,nrhs)*cff2
#    ifdef LIMIT_BSTRESS
          cff3=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
     &               v(i,j,1,nrhs)*cff2
#    ifdef LIMIT_BSTRESS
          cff3=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff3)
#    endif
        END DO
      END DO
#    elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
     &               u(i,j,1,nrhs)
#    ifdef LIMIT_BSTRESS
          cff1=cff*0.5_r8*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.0_r8, bustr(i,j))*                          &
     &               MIN(ABS(bustr(i,j)),                               &
     &                   ABS(u(i,j,1,nrhs))*cff1)
#    endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
     &               v(i,j,1,nrhs)
#    ifdef LIMIT_BSTRESS
          cff1=cff*0.5_r8*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.0_r8, bvstr(i,j))*                          &
     &               MIN(ABS(bvstr(i,j)),                               &
     &                   ABS(v(i,j,1,nrhs))*cff1)
#    endif
        END DO
      END DO
#   endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr)
#   endif
#  endif

      RETURN
      END SUBROUTINE set_vbc_tile

# else

!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
#  if defined UV_LDRAG
     &                   GRID(ng) % rdrag,                              &
#  elif defined UV_QDRAG
     &                   GRID(ng) % rdrag2,                             &
#  endif
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif

      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs, kstp, knew,                        &
#  if defined UV_LDRAG
     &                         rdrag,                                   &
#  elif defined UV_QDRAG
     &                         rdrag2,                                  &
#  endif
     &                         ubar, vbar, bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, knew
!
#  ifdef ASSUMED_SHAPE
#   ifdef UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:,LBj:)
#   endif
#   ifdef UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#  else
#   ifdef UV_LDRAG
      real(r8), intent(in) :: rdrag(LBi:UBi,LBj:UBj)
#   endif
#   ifdef UV_QDRAG
      real(r8), intent(in) :: rdrag2(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff1, cff2

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic barotropic bottom momentum stress (m2/s2).
!-----------------------------------------------------------------------

#  if defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          bustr(i,j)=0.5_r8*(rdrag(i-1,j)+rdrag(i,j))*                  &
     &               ubar(i,j,krhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          bvstr(i,j)=0.5_r8*(rdrag(i,j-1)+rdrag(i,j))*                  &
     &               vbar(i,j,krhs)
        END DO
      END DO
#  elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(vbar(i  ,j  ,krhs)+                             &
     &                  vbar(i  ,j+1,krhs)+                             &
     &                  vbar(i-1,j  ,krhs)+                             &
     &                  vbar(i-1,j+1,krhs))
          cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(rdrag2(i-1,j)+rdrag2(i,j))*                &
     &               ubar(i,j,krhs)*cff2
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
     &                  ubar(i+1,j  ,krhs)+                             &
     &                  ubar(i  ,j-1,krhs)+                             &
     &                  ubar(i+1,j-1,krhs))
          cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
          bvstr(i,j)=0.5_r8*(rdrag2(i,j-1)+rdrag2(i,j))*                &
     &               vbar(i,j,krhs)*cff2
        END DO
      END DO
#  endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bustr, bvstr)
#  endif

      RETURN
      END SUBROUTINE set_vbc_tile
# endif
#endif
      END MODULE set_vbc_mod
