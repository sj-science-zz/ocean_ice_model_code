        MODULE mod_ice
!
!   temporary storage for ice data
!
!   ice_strx - ice stress in i-direction on B-grid T-point
!   ice_stry - ice stress in j-direction, B-grid T-point
!   aice - ice fraction on rho-point

!  CICE B-grid
!
!   U(i-1,j) --------------U(i,j)
!        |                    |
!        |                    |
!        |      T(i,j)        |
!        |                    |
!        |                    |
!   U(i-1,j-1)-------------U(i,j-1)
!
!
! ROMS C-grid
!
!         ------V(i,j+1)------
!        |                    |
!        |                    |
!     U(i,j)    T(i,j)   U(i+1,j)
!        |                    |
!        |                    |
!         ------V(i,j)--------
! 
!  Notice the difference in boundary condition at i=0 or j=0 

        USE mod_kinds

        implicit none

        TYPE T_ICE 
          real(r8), pointer :: aice(:,:)      ! ice concentration
          real(r8), pointer :: aice_u(:,:)    ! ice concentration on u-points
          real(r8), pointer :: aice_v(:,:)    ! ice concentration on v-points
          real(r8), pointer :: fhocnAI(:,:)   ! nonradiative heat flux under ice
          real(r8), pointer :: fswthruAI(:,:) ! radiative heat flux under ice
          real(r8), pointer :: freshAI(:,:)   ! freshwater flux ice to ocean (kg/s/m**2)
          real(r8), pointer :: fsaltAI(:,:)   ! salt flux ice to ocean (kg/s/m**2)
          real(r8), pointer :: strx(:,:)      ! x-stress on CICE psi-points
          real(r8), pointer :: stry(:,:)      ! y-stress on CICE psi-points
          real(r8), pointer :: stru(:,:)      ! x-stress on ROMS u-points
          real(r8), pointer :: strv(:,:)      ! y-stress on ROMS v-points
#if defined CICE_COUPLING
          real(r8), pointer :: qfraz(:,:)     ! Frazil formed in K*m
! ROMS-to-CICE coupling fields accumulated over coupling interval
          real(r8), pointer :: qfraz_accum(:,:) ! Frazil
          real(r8), pointer :: melt_accum(:,:) ! Melt potential
          real(r8), pointer :: sst_accum(:,:)   ! Surface temperature
          real(r8), pointer :: sss_accum(:,:)   ! Surface salinity
          real(r8), pointer :: u_accum(:,:,:)   ! u (3D for now)
          real(r8), pointer :: v_accum(:,:,:)   ! v (3D for now)
          real(r8), pointer :: ssh_accum(:,:)   ! Sea surface height
#endif
        END TYPE T_ICE



! Note that all fluxes come multiplied by the (CICE) ice concentration, so e.g. strx = a*strx

        TYPE (T_ICE), allocatable :: ICE(:)

        CONTAINS

        SUBROUTINE allocate_ice (ng, LBi, UBi, LBj, UBj)

        USE mod_param
      
        integer, intent(in) :: ng, LBi, UBi, LBj, UBj

        real(r8), parameter :: IniVal = 0.0_r8


        if (ng.eq.1) allocate ( ICE(Ngrids) )
        allocate ( ICE(ng)%aice(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%aice_u(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%aice_v(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%fhocnAI(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%fswthruAI(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%freshAI(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%fsaltAI(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%strx(LBi:UBi,LBj:UBj) ) ! Could be local variable in ocean_coupler? 
        allocate ( ICE(ng)%stry(LBi:UBi,LBj:UBj) ) ! Could be local variable in ocean_coupler?
        allocate ( ICE(ng)%stru(LBi:UBi,LBj:UBj) )
        allocate ( ICE(ng)%strv(LBi:UBi,LBj:UBj) )
#if defined CICE_COUPLING
        allocate ( ICE(ng) % qfraz(LBi:UBi,LBj:UBj))
        allocate ( ICE(ng) % qfraz_accum(LBi:UBi,LBj:UBj))
        allocate ( ICE(ng) % melt_accum(LBi:UBi,LBj:UBj))
        allocate ( ICE(ng) % sst_accum(LBi:UBi,LBj:UBj))
        allocate ( ICE(ng) % sss_accum(LBi:UBi,LBj:UBj))
        allocate ( ICE(ng) % u_accum(LBi:UBi,LBj:UBj,N(ng)))
        allocate ( ICE(ng) % v_accum(LBi:UBi,LBj:UBj,N(ng)))
        allocate ( ICE(ng) % ssh_accum(LBi:UBi,LBj:UBj))
#endif
        
        ICE(ng)%aice = IniVal
        ICE(ng)%aice_u = IniVal
        ICE(ng)%aice_v = IniVal
        ICE(ng)%fhocnAI = IniVal
        ICE(ng)%freshAI = IniVal
        ICE(ng)%fsaltAI = IniVal
        ICE(ng)%fswthruAI = IniVal
        ICE(ng)%strx = IniVal
        ICE(ng)%stry = IniVal
        ICE(ng)%stru = IniVal
        ICE(ng)%strv = IniVal
# if defined CICE_COUPLING
        ICE(ng) % qfraz = IniVal
        ICE(ng) % qfraz_accum = IniVal
        ICE(ng) % melt_accum = IniVal
        ICE (ng) % sst_accum = IniVal
        ICE (ng) % sss_accum = IniVal
        ICE (ng) % u_accum = IniVal
        ICE (ng) % v_accum = IniVal
        ICE (ng) % ssh_accum = IniVal
#endif
        
        RETURN
        END SUBROUTINE allocate_ice
        END MODULE mod_ice
  
  
