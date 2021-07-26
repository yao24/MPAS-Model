! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts3
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This module contains the LTS3 time integration and LTS init routines
!
!-----------------------------------------------------------------------

module ocn_time_integration_lts3 

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled
   use ocn_wetting_drying

   use ocn_effective_density_in_land_ice

   use ocn_time_integration_ssprk3

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_lts3, ocn_time_integration_lts3_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_lts3
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This routine integrates one timestep (dt) using an LTS3 time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_lts3(domain,dt)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    !   Local Time Stepping of order 3 scheme
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    type (domain_type), intent(inout) :: domain
    real (kind=RKIND), intent(in) :: dt

    integer :: iCell, iEdge, iRegion, k, ic, ie, im, M, nRegions, err

    type (block_type), pointer :: block

    type (mpas_pool_type), pointer :: tendPool
    type (mpas_pool_type), pointer :: tracersTendPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: scratchPool
    type (mpas_pool_type), pointer :: swForcingPool

    ! LTS Pools
    type (mpas_pool_type), pointer :: LTSPool
    type (mpas_pool_type), pointer :: firstStagePool
    type (mpas_pool_type), pointer :: tracersPoolFirstStage
    type (mpas_pool_type), pointer :: secondStagePool
    type (mpas_pool_type), pointer :: tracersPoolSecondStage
    type (mpas_pool_type), pointer :: tendSum1stPool
    type (mpas_pool_type), pointer :: tracersTendSum1stPool
    type (mpas_pool_type), pointer :: tendSum2ndPool
    type (mpas_pool_type), pointer :: tracersTendSum2ndPool
    type (mpas_pool_type), pointer :: tendSum3rdPool
    type (mpas_pool_type), pointer :: tracersTendSum3rdPool
    type (mpas_pool_type), pointer :: prevFirstStagePool, nextFirstStagePool, prevTendSum1stPool, nextTendSum1stPool, prevTendSum2ndPool, nextTendSum2ndPool
    type (mpas_pool_type), pointer :: prevSecondStagePool, nextSecondStagePool, prevTendSum3rdPool, nextTendSum3rdPool

    ! Dimensions
    integer, pointer :: nCells, nEdges, nVertLevels, haveFine, haveCoarse, haveInterface

    ! Config options
    integer, pointer :: config_dt_scaling_LTS
    logical, pointer :: config_prescribe_velocity, config_prescribe_thickness
    logical, pointer :: config_filter_btr_mode, config_use_freq_filtered_thickness
    logical, pointer :: config_use_GM
    logical, pointer :: config_use_cvmix_kpp
    logical, pointer :: config_use_tracerGroup
    logical, pointer :: config_disable_thick_all_tend
    logical, pointer :: config_disable_vel_all_tend
    logical, pointer :: config_disable_tr_all_tend
    real (kind=RKIND), pointer :: config_mom_del4
    real (kind=RKIND), pointer :: config_drying_min_cell_height
    logical, pointer :: config_use_wetting_drying
    logical, pointer :: config_verify_not_dry
    logical, pointer :: config_prevent_drying
    logical, pointer :: config_zero_drying_velocity
    character (len=StrKIND), pointer :: config_land_ice_flux_mode

    ! State indices
    integer, pointer :: indexTemperature
    integer, pointer :: indexSalinity

    ! Diagnostics Indices
    integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
    integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional

    ! Mesh array pointers
    integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop
    real (kind=RKIND), dimension(:), pointer :: bottomDepth

    ! Tend Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend, lowFreqDivergenceTend, normalVelocityTend, layerThicknessTend
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum1st, lowFreqDivergenceTendSum1st, normalVelocityTendSum1st, layerThicknessTendSum1st
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum2nd, lowFreqDivergenceTendSum2nd, normalVelocityTendSum2nd, layerThicknessTendSum2nd
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTendSum3rd, lowFreqDivergenceTendSum3rd, normalVelocityTendSum3rd, layerThicknessTendSum3rd
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend, tracersGroupTendSum1st, tracersGroupTendSum2nd, tracersGroupTendSum3rd

    ! Diagnostics Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge
    real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
    real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity, normalGMBolusVelocity
    real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
    real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
    real (kind=RKIND), dimension(:), pointer :: gradSSH
    real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
    real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
    real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, sshGradient

    ! State Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocityMPIHaloUpdate, layerThicknessMPIHaloUpdate, tracersGroupMPIHaloUpdate
    real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew, normalVelocityFirstStage, normalVelocitySecondStage
    real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew, layerThicknessFirstStage, layerThicknessSecondStage
    real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew, highFreqThicknessFirstStage, highFreqThicknessSecondStage
    real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew, lowFreqDivergenceFirstStage, lowFreqDivergenceSecondStage 
    real (kind=RKIND), dimension(:,:), pointer :: wettingVelocity
    real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew, sshFirstStage, sshSecondStage
    real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupCur, tracersGroupNew, tracersGroupFirstStage, tracersGroupSecondStage

    ! Diagnostics Field Pointers
    type (field1DReal), pointer :: boundaryLayerDepthField, effectiveDensityField
    type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField

    ! Tracer Group Iteartion
    type (mpas_pool_iterator_type) :: groupItr
    character (len=StrKIND) :: modifiedGroupName
    character (len=StrKIND) :: configName

    ! Tidal boundary condition
    logical, pointer :: config_use_tidal_forcing
    character (len=StrKIND), pointer :: config_tidal_forcing_type
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue
    real (kind=RKIND), dimension(:,:), pointer :: restingThickness
    real (kind=RKIND) :: totalDepth

    ! LTS objects
    real (kind=RKIND) :: dtFine, alpha, alphaHat, beta, betaHat, gam, gamHat
    integer, dimension(:,:), pointer :: nCellsInLTSHalo, nEdgesInLTSHalo
    integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
    real (kind=RKIND) :: weightOld, weightNew, weightTend
    real (kind=RKIND) :: weightTendSum1st, weightTendSum2nd, weightTendSum3rd

    ! Get config options
    call mpas_pool_get_config(domain % configs, 'config_dt_scaling_LTS', config_dt_scaling_LTS)
    call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
    call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
    call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
    call mpas_pool_get_config(domain % configs, 'config_use_GM', config_use_GM)
    call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
    call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
    call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
    call mpas_pool_get_config(domain % configs, 'config_use_wetting_drying', config_use_wetting_drying)
    call mpas_pool_get_config(domain % configs, 'config_prevent_drying', config_prevent_drying)
    call mpas_pool_get_config(domain % configs, 'config_verify_not_dry', config_verify_not_dry)
    call mpas_pool_get_config(domain % configs, 'config_drying_min_cell_height', config_drying_min_cell_height)
    call mpas_pool_get_config(domain % configs, 'config_zero_drying_velocity', config_zero_drying_velocity)
    call mpas_pool_get_config(domain % configs, 'config_use_tidal_forcing', config_use_tidal_forcing)
    call mpas_pool_get_config(domain % configs, 'config_tidal_forcing_type', config_tidal_forcing_type)

    M = config_dt_scaling_LTS
    nRegions = 2

    err = 0

    dtFine = dt / M

    weightOld = 0.75
    weightNew = 0.25
    weightTend = 0.25

    weightTendSum1st = 1.0 / 6.0
    weightTendSum2nd = 1.0 / 6.0
    weightTendSum3rd = 2.0 / 3.0

    call mpas_timer_start("create extra pools LTS")
    block => domain % blocklist
    do while (associated(block))
    call mpas_pool_get_subpool(block % structs, 'state', statePool)
    call mpas_pool_get_subpool(block % structs, 'tracers', tracersPool)
    call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

    call mpas_pool_create_pool(firstStagePool)
    call mpas_pool_clone_pool(statePool, firstStagePool, 1)
    call mpas_pool_create_pool(secondStagePool)
    call mpas_pool_clone_pool(statePool, secondStagePool, 1)

    call mpas_pool_create_pool(tendSum1stPool)
    call mpas_pool_clone_pool(tendPool, tendSum1stPool, 1)
    call mpas_pool_create_pool(tendSum2ndPool)
    call mpas_pool_clone_pool(tendPool, tendSum2ndPool, 1)
    call mpas_pool_create_pool(tendSum3rdPool)
    call mpas_pool_clone_pool(tendPool, tendSum3rdPool, 1)

    call mpas_pool_add_subpool(block % structs, 'first_stage_state', firstStagePool)
    call mpas_pool_add_subpool(block % structs, 'second_stage_state', secondStagePool)
    call mpas_pool_add_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
    call mpas_pool_add_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
    call mpas_pool_add_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)

    call mpas_pool_get_array(tendSum1stPool, 'normalVelocity', normalVelocityTendSum1st)
    call mpas_pool_get_array(tendSum1stPool, 'layerThickness', layerThicknessTendSum1st)
    call mpas_pool_get_subpool(tendSum1stPool, 'tracersTend', tracersTendSum1stPool)

    call mpas_pool_get_array(tendSum2ndPool, 'normalVelocity', normalVelocityTendSum2nd)
    call mpas_pool_get_array(tendSum2ndPool, 'layerThickness', layerThicknessTendSum2nd)
    call mpas_pool_get_subpool(tendSum2ndPool, 'tracersTend', tracersTendSum2ndPool)

    call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', normalVelocityTendSum3rd)
    call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', layerThicknessTendSum3rd)
    call mpas_pool_get_subpool(tendSum3rdPool, 'tracersTend', tracersTendSum3rdPool)

    normalVelocityTendSum1st(:,:) = 0.0
    layerThicknessTendSum1st(:,:) = 0.0

    normalVelocityTendSum2nd(:,:) = 0.0
    layerThicknessTendSum2nd(:,:) = 0.0

    normalVelocityTendSum3rd(:,:) = 0.0
    layerThicknessTendSum3rd(:,:) = 0.0

    call mpas_pool_begin_iteration(tracersPool)
       do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
          if ( groupItr % memberType == MPAS_POOL_FIELD ) then
             configName = 'config_use_' // trim(groupItr % memberName)
             call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

             if ( config_use_tracerGroup ) then
                modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                call mpas_pool_get_array(tracersTendSum1stPool, modifiedGroupName, tracersGroupTendSum1st)
                call mpas_pool_get_array(tracersTendSum2ndPool, modifiedGroupName, tracersGroupTendSum2nd)
                call mpas_pool_get_array(tracersTendSum3rdPool, modifiedGroupName, tracersGroupTendSum3rd)
                if ( associated(tracersGroupTendSum1st) ) then
                   tracersGroupTendSum1st(:,:,:) = 0.0
                else if (associated(tracersGroupTendSum2nd) ) then
                   tracersGroupTendSum2nd(:,:,:) = 0.0
                else if (associated(tracersGroupTendSum3rd) ) then
                   tracersGroupTendSum3rd(:,:,:) = 0.0
                end if
             end if
          end if
       end do

       block => block % next
    end do

    block => domain % blocklist
    do while(associated(block))
       if (associated(block % prev)) then
          call mpas_pool_get_subpool(block % prev % structs, 'first_stage_state', prevFirstStagePool)
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_1st', tendSum1stPool)
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_2nd', tendSum2ndPool)         
       else
          nullify(prevFirstStagePool)
          nullify(prevTendSum1stPool)
          nullify(prevTendSum2ndPool)
       end if

       if (associated(block % next)) then
          call mpas_pool_get_subpool(block % next % structs, 'first_stage_state', nextFirstStagePool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_1st', nextTendSum1stPool)
            call mpas_pool_get_subpool(block % next % structs, 'tend_sum_2nd', nextTendSum2ndPool)
       else
          nullify(nextFirstStagePool)
          nullify(nextTendSum1stPool)
          nullify(nextTendSum2ndPool)
       end if

       call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)
       call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
       call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)

       if (associated(prevFirstStagePool) .and. associated(nextFirstStagePool)) then
          call mpas_pool_link_pools(firstStagePool, prevFirstStagePool, nextFirstStagePool)
       else if (associated(prevFirstStagePool)) then
          call mpas_pool_link_pools(firstStagePool, prevFirstStagePool)
       else if (associated(nextFirstStagePool)) then
          call mpas_pool_link_pools(firstStagePool, nextPool=nextFirstStagePool)
       else
          call mpas_pool_link_pools(firstStagePool)
       end if

       if (associated(prevTendSum1stPool) .and. associated(nextTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool, nextTendSum1stPool)
       else if (associated(prevTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool)
       else if (associated(nextTendSum1stPool)) then
          call mpas_pool_link_pools(tendSum1stPool,nextPool=nextTendSum1stPool)
       else
          call mpas_pool_link_pools(tendSum1stPool)
       end if

       if (associated(prevTendSum2ndPool) .and. associated(nextTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool, nextTendSum2ndPool)
       else if (associated(prevTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool)
       else if (associated(nextTendSum2ndPool)) then
          call mpas_pool_link_pools(tendSum2ndPool,nextPool=nextTendSum2ndPool)
       else
          call mpas_pool_link_pools(tendSum2ndPool)
       end if

       call mpas_pool_link_parinfo(block, firstStagePool)
       call mpas_pool_link_parinfo(block, tendSum1stPool)
       call mpas_pool_link_parinfo(block, tendSum2ndPool)

       if (associated(block % prev)) then
          call mpas_pool_get_subpool(block % prev % structs, 'second_stage_state', prevSecondStagePool)
          call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_3rd', prevTendSum3rdPool)
       else
          nullify(prevSecondStagePool)
          nullify(prevTendSum3rdPool)
       end if

       if (associated(block % next)) then 
          call mpas_pool_get_subpool(block % next % structs, 'second_stage_state', nextSecondStagePool)
          call mpas_pool_get_subpool(block % next % structs, 'tend_sum_3rd', nextTendSum3rdPool)
       else
          nullify(nextSecondStagePool)
          nullify(nextTendSum3rdPool)
       end if

       call mpas_pool_get_subpool(block % structs, 'second_stage_state', secondStagePool)
       call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)

       if (associated(prevSecondStagePool) .and. associated(nextSecondStagePool)) then
          call mpas_pool_link_pools(secondStagePool, prevSecondStagePool, nextSecondStagePool)
       else if (associated(prevSecondStagePool)) then
          call mpas_pool_link_pools(secondStagePool, prevSecondStagePool)
       else if (associated(nextSecondStagePool)) then
          call mpas_pool_link_pools(secondStagePool, nextPool=nextSecondStagePool)
       else
          call mpas_pool_link_pools(secondStagePool)
       end if

       if (associated(prevTendSum3rdPool) .and. associated(nextTendSum3rdPool)) then
          call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool, nextTendSum3rdPool)
       else if (associated(prevTendSum3rdPool)) then
          call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool)
       else if (associated(nextTendSum3rdPool)) then 
          call mpas_pool_link_pools(tendSum3rdPool,nextPool=nextTendSum3rdPool)
       else
          call mpas_pool_link_pools(tendSum3rdPool)
       end if

       call mpas_pool_link_parinfo(block, secondStagePool)
       call mpas_pool_link_parinfo(block, tendSum3rdPool)

       block => block % next
    end do

call mpas_timer_stop("create extra pools LTS")

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! BEGIN LTS SCHEME
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! --------------------------------------  STEP  ------------------------------------------
     ! --- compute the first stage of SSPRK for interface layer 1, interface layer 2 and coarse

      ! --- update halos for diagnostic variables
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("LTS3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("LTS3-boundary layer depth halo update")
      end if

      call mpas_timer_start("LTS3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("LTS3-diagnostic halo update")
      call mpas_threading_barrier()

      ! --- compute tendencies for high frequency thickness
      if (config_use_freq_filtered_thickness) then
         !TODO: MISSING FOR NOW (it would add too much complexity to something already complex)
      end if

      ! --- require that cells don't dry out
      if (config_use_wetting_drying) then
         !TODO: MISSING FOR NOW
      end if

      ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
      call mpas_timer_start("LTS3-update diagnostic variables and compute tendencies")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool) 


         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

         ! DIAGNOSTICS UPDATE ---

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
         call mpas_threading_barrier()

         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalVelocityCur(:, iEdge)
         end do
         !$omp end do

         ! --- compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
         if (config_use_GM) then
            call ocn_gm_compute_Bolus_velocity(statePool, diagnosticsPool, meshPool, scratchPool, timeLevelIn=1)
         end if
         call mpas_threading_barrier()

         if (config_use_GM) then
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalTransportVelocity(:,iEdge) + normalGMBolusVelocity(:,iEdge)
         end do
         !$omp end do
         end if

         ! TENDENCIES COMPUTATION ---

         call ocn_time_integrator_compute_vel_tends(block, dt, err, statePool, 1)
         call ocn_time_integrator_compute_thick_tends(block, dt, err, statePool, 1)
         call ocn_time_integrator_compute_tracer_tends( block, dt, err, statePool, 1)

         block => block % next
      end do
      call mpas_timer_stop("LTS3-update diagnostic variables and compute tendencies")
      call mpas_threading_barrier()

      call mpas_timer_start("LTS advance soln")
      ! --- advance solution for interface layer 1, interface layer 2 and coarse with first stage of SSPRK
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(firstStagePool, 'tracers', tracersPoolFirstStage)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)

         call mpas_pool_get_array(firstStagePool, 'normalVelocity', normalVelocityFirstStage, 1)
         call mpas_pool_get_array(firstStagePool, 'layerThickness', layerThicknessFirstStage, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', wettingVelocity)

         ! --- NORMAL VELOCITY

         ! --- interface layers
         do iRegion =1,nRegions
            do ie = 1, nEdgesInLTSHalo(iRegion,2)
               iEdge = edgesInLTSHalo(iRegion,2,ie)
               normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
               normalVelocityNew(:, iEdge) = normalVelocityNew(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
            end do
         end do
         ! --- coarse
         do ie = 1, nEdgesInLTSHalo(2,1)
            iEdge = edgesInLTSHalo(2,1,ie)
            normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
            normalVelocityNew(:, iEdge) = normalVelocityNew(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
         end do
         ! --- fine layers close to interface layers
         do ie = 1, nEdgesInLTSHalo(1,3)
            iEdge = edgesInLTSHalo(1,3,ie)
            normalVelocityFirstStage(:,iEdge) =  normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge)
             normalVelocityNew(:, iEdge) = normalVelocityNew(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
         end do

     
         ! --- LAYER THICKNESS

         ! --- interface layers
         do iRegion =1,nRegions
            do ic = 1, nCellsInLTSHalo(iRegion,2)
               iCell = cellsInLTSHalo(iRegion,2,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
               end do
            end do
         end do
         ! --- coarse
         do ic = 1, nCellsInLTSHalo(2,1)
            iCell = cellsInLTSHalo(2,1,ic)
            do k = 1, maxLevelCell(iCell)
               layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
            end do
         end do
         ! --- fine layers close to interface layers
         do ic = 1, nCellsInLTSHalo(1,3)
            iCell = cellsInLTSHalo(2,1,ic)
            do k = 1, maxLevelCell(iCell)
               layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
            end do
         end do

         ! --- TRACERS

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               configName = 'config_use_' // trim(groupItr % memberName)
               call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

               if ( config_use_tracerGroup ) then
                  call mpas_pool_get_array(tracersPoolFirstStage, groupItr % memberName, tracersGroupFirstStage, 1)
                  call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                  modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                  call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                  if ( associated(tracersGroupFirstStage) .and. associated(tracersGroupTend) ) then
                     !$omp do schedule(runtime) private(k)
                     ! --- interface layers
                     do iRegion =1,nRegions
                        do ic = 1, nCellsInLTSHalo(iRegion,2)
                           iCell = cellsInLTSHalo(iRegion,2,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) & 
                                                                +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                           end do
                        end do
                     end do
                     ! --- coarse
                     do ic = 1, nCellsInLTSHalo(2,1)
                        iCell = cellsInLTSHalo(2,1,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                        end do
                     end do
                     ! --- fine layers close to interface layers
                     do ic = 1, nCellsInLTSHalo(1,3)
                        iCell = cellsInLTSHalo(1,3,ic)
                        do k = 1, maxLevelCell(iCell)
                           tracersGroupFirstStage(:,k,iCell) = ( tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                +   dt * tracersGroupTend(:,k,iCell) ) / layerThicknessFirstStage(k,iCell)
                        end do
                     end do
                     !$omp end do
                  end if
               end if
            end if
         end do

         block => block % next
      end do
      call mpas_timer_stop("LTS3 advance soln")
      call mpas_threading_barrier()

      call mpas_timer_start("copy for halo update  LTS")
      ! --- copy for MPI halo update of the solution
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)

         call mpas_pool_get_array(LTSPool, 'normalVelocityMPIHaloUpdate', normalVelocityMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'layerThicknessMPIHaloUpdate', layerThicknessMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'tracersGroupMPIHaloUpdate', tracersGroupMPIHaloUpdate) !  TODO this will likely have to be modified!!!!!!!!!!!

         call mpas_pool_get_array(firstStagePool, 'normalVelocity', normalVelocityFirstStage, 1)
         call mpas_pool_get_array(firstStagePool, 'layerThickness', layerThicknessFirstStage, 1)

         ! TODO: take care of the tracers!!!!!!!!!!!!!!!!!!

         normalVelocityMPIHaloUpdate(:,:) = normalVelocityFirstStage(:,:)
         layerThicknessMPIHaloUpdate(:,:) = layerThicknessFirstStage(:,:)

         block => block % next
      end do
      ! ---  copy for MPI halo update of the solution
      call mpas_timer_stop("copy for halo update  LTS")

      call mpas_timer_start("halo updates  LTS")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocityMPIHaloUpdate')
      call mpas_dmpar_field_halo_exch(domain, 'layerThicknessMPIHaloUpdate')
      !call mpas_dmpar_field_halo_exch(domain, 'tracersMPIHaloUpdate') TODO
      call mpas_timer_stop("halo updates  LTS")

      call mpas_timer_start("copy for halo update  LTS")
      ! ---  copy back for MPI halo update of the solution
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool) 
         call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)

         call mpas_pool_get_array(LTSPool, 'normalVelocityMPIHaloUpdate', normalVelocityMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'layerThicknessMPIHaloUpdate', layerThicknessMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'tracersGroupMPIHaloUpdate', tracersGroupMPIHaloUpdate) !  TODO this will likely have to be modified!!!!!!!!!!!

         call mpas_pool_get_array(firstStagePool, 'normalVelocity', normalVelocityFirstStage, 1)
         call mpas_pool_get_array(firstStagePool, 'layerThickness', layerThicknessFirstStage, 1)

         ! TODO: take care of the tracers!!!!!!!!!!!!!!!!!!

         normalVelocityFirstStage(:,:) = normalVelocityMPIHaloUpdate(:,:)
         layerThicknessFirstStage(:,:) = layerThicknessMPIHaloUpdate(:,:)

         block => block % next
      end do
      ! ---  copy back for MPI halo update of the solution
      call mpas_timer_stop("copy for halo update  LTS")


      ! --------------------------------------  STEP  ------------------------------------------
      ! --- compute the second stage of SSPRK for interface layer 1, interface layer 2 and coarse

      ! --- update halos for diagnostic variables.
      if (config_use_cvmix_kpp) then
         call mpas_timer_start("LTS3-boundary layer depth halo update")
         call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
         call mpas_timer_stop("LTS3-boundary layer depth halo update")
      end if

      call mpas_timer_start("LTS3-diagnostic halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
      if (config_mom_del4 > 0.0_RKIND) then
         call mpas_dmpar_field_halo_exch(domain, 'divergence')
         call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
      end if
      call mpas_timer_stop("LTS3-diagnostic halo update")
      call mpas_threading_barrier()

      ! --- compute tendencies for high frequency thickness
      if (config_use_freq_filtered_thickness) then
         !TODO: MISSING FOR NOW (it would add too much complexity to something
         !already complex)
      end if

      ! --- require that cells don't dry out
      if (config_use_wetting_drying) then
         !TODO: MISSING FOR NOW
      end if

      ! --- compute diagnostics and tendencies for velocity, thickness, and tracers
      call mpas_timer_start("LTS3-update diagnostic variables and compute tendencies")
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
         call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

         call mpas_pool_get_array(firstStagePool, 'normalVelocity', normalVelocityFirstStage, 1)
         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

         if(haveCoarse == 1 .or. haveInterface == 1) then

            ! DIAGNOSTICS UPDATE ---

            call ocn_diagnostic_solve(dt, firstStagePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 1)
            call mpas_threading_barrier()

            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalVelocityFirstStage(:, iEdge)
            end do
            !$omp end do

            ! --- compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
            if (config_use_GM) then
               call ocn_gm_compute_Bolus_velocity(firstStagePool, diagnosticsPool, meshPool, scratchPool, timeLevelIn=1)
            end if
            call mpas_threading_barrier()

            if (config_use_GM) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalTransportVelocity(:, iEdge) = normalTransportVelocity(:,iEdge) + normalGMBolusVelocity(:,iEdge)
            end do
            !$omp end do
            end if

            ! TENDENCIES COMPUTATION ---

            call ocn_time_integrator_compute_vel_tends(block, dt, err, firstStagePool, 1)
            call ocn_time_integrator_compute_thick_tends(block, dt, err, firstStagePool, 1)
            call ocn_time_integrator_compute_tracer_tends( block, dt, err, firstStagePool, 1)
           
         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3-update diagnostic variables and compute tendencies")
      call mpas_threading_barrier()

      call mpas_timer_start("LTS advance soln")
      ! --- advance solution (2nd stage of SSPRK for the coarse using dt)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'first_stage_state', firstStagePool)
         call mpas_pool_get_subpool(block % structs, 'second_stage_state', secondStagePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(firstStagePool, 'tracers', tracersPoolFirstStage)
         call mpas_pool_get_subpool(secondStagePool, 'tracers', tracersPoolSecondStage)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo',cellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo',nCellsInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo',nEdgesInLTSHalo)
         call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)

         call mpas_pool_get_array(firstStagePool, 'normalVelocity', normalVelocityFirstStage, 1)
         call mpas_pool_get_array(firstStagePool, 'layerThickness', layerThicknessFirstStage, 1)

         call mpas_pool_get_array(secondStagePool, 'normalVelocity', normalVelocitySecondStage, 1)
         call mpas_pool_get_array(secondStagePool, 'layerThickness', layerThicknessSecondStage, 1)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)

         call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
         call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

         call mpas_pool_get_array(diagnosticsPool, 'wettingVelocity', wettingVelocity)

         if(haveCoarse == 1 .or. haveInterface == 1) then

            ! --- NORMAL VELOCITY

            ! --- coarse
            do ie = 1, nEdgesInLTSHalo(2,1)
               iEdge = edgesInLTSHalo(2,1,ie)
               normalVelocitySecondStage(:,iEdge) = weightOld * normalVelocityCur(:,iEdge) &
                                                  + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                  + weightTend * dt * normalVelocityTend(:,iEdge)
               normalVelocitySecondStage(:, iEdge) = normalVelocitySecondStage(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
            end do

            ! --- interface layers
            do iRegion = 1, nRegions
               do ie = 1, nEdgesInLTSHalo(iRegion,2)
                  iEdge = edgesInLTSHalo(iRegion,2,ie)
                  normalVelocitySecondStage(:,iEdge) = weightOld * normalVelocityCur(:,iEdge) &
                                                  + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                  + weightTend * dt * normalVelocityTend(:,iEdge)
                  normalVelocitySecondStage(:, iEdge) = normalVelocitySecondStage(:, iEdge) * (1.0_RKIND - wettingVelocity(:, iEdge))
               end do
            end do

            ! LAYER THICKNESS

            ! --- coarse
            do ic = 1, nCellsInLTSHalo(2,1)
               iCell = cellsInLTSHalo(2,1,ic)
               do k = 1, maxLevelCell(iCell)
                  layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) & 
                                                     + weightNew * layerThicknessFirstStage(k,iCell) &
                                                     + weightTend * dt * layerThicknessTend(k,iCell)
               end do
            end do

            ! --- interface layers
            do iRegion = 1, nRegions
               do ic = 1, nCellsInLTSHalo(iRegion,2)
                  iCell = cellsInLTSHalo(iRegion,2,ic)
                  do k = 1, maxLevelCell(iCell)
                     layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) &
                                                        + weightNew * layerThicknessFirstStage(k,iCell) &
                                                        + weightTend * dt * layerThicknessTend(k,iCell)
                  end do
               end do
            end do

            ! TRACERS

            call mpas_pool_begin_iteration(tracersPool)
            do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
               if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                  configName = 'config_use_' // trim(groupItr % memberName)
                  call mpas_pool_get_config(block % configs, configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPoolFirstStage, groupItr % memberName, tracersGroupFirstStage, 1)
                     call mpas_pool_get_array(tracersPoolSecondStage, groupItr % memberName, tracersGroupSecondStage, 1)
                     call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)

                     modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)
                     if ( associated(tracersGroupSecondStage) .and. associated(tracersGroupTend) ) then
                        !$omp do schedule(runtime) private(k)
                        ! --- interface layers
                        do iRegion =1,nRegions
                           do ic = 1, nCellsInLTSHalo(iRegion,2)
                              iCell = cellsInLTSHalo(iRegion,2,ic)
                              do k = 1, maxLevelCell(iCell)
                                 tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                        + weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell)  &
                                                                        + weightTend * dt * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                              end do
                           end do
                        end do
                        ! --- coarse
                        do ic = 1, nCellsInLTSHalo(2,1)
                           iCell = cellsInLTSHalo(2,1,ic)
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupSecondStage(:,k,iCell) = ( weightOld * tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) &
                                                                     + weightNew * tracersGroupFirstStage(:,k,iCell) * layerThicknessFirstStage(k,iCell)  &
                                                                     + weightTend * dt * tracersGroupTend(:,k,iCell) ) / layerThicknessSecondStage(k,iCell)
                           end do
                        end do
                        !$omp end do
                     end if
                  end if
               end if
            end do

         end if

         block => block % next
      end do
      call mpas_timer_stop("LTS3 advance soln")
      call mpas_threading_barrier()

      call mpas_timer_start("copy for halo update  LTS")
      ! --- copy for MPI halo update of the solution
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'second_stage_state', secondStagePool)

         call mpas_pool_get_array(LTSPool, 'normalVelocityMPIHaloUpdate', normalVelocityMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'layerThicknessMPIHaloUpdate', layerThicknessMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'tracersGroupMPIHaloUpdate', tracersGroupMPIHaloUpdate) !  TODO this will likely have to be modified!!!!!!!!!!!

         call mpas_pool_get_array(secondStagePool, 'normalVelocity', normalVelocitySecondStage, 1)
         call mpas_pool_get_array(secondStagePool, 'layerThickness', layerThicknessSecondStage, 1)

         ! TODO: take care of the tracers!!!!!!!!!!!!!!!!!!

         normalVelocityMPIHaloUpdate(:,:) = normalVelocitySecondStage(:,:)
         layerThicknessMPIHaloUpdate(:,:) = layerThicknessSecondStage(:,:)

         block => block % next
      end do
      ! ---  copy for MPI halo update of the solution
      call mpas_timer_stop("copy for halo update  LTS")

      call mpas_timer_start("halo updates  LTS")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocityMPIHaloUpdate')
      call mpas_dmpar_field_halo_exch(domain, 'layerThicknessMPIHaloUpdate')
      !call mpas_dmpar_field_halo_exch(domain, 'tracersMPIHaloUpdate') TODO
      call mpas_timer_stop("halo updates  LTS")

      call mpas_timer_start("copy for halo update  LTS")
      ! ---  copy back for MPI halo update of the solution
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
         call mpas_pool_get_subpool(block % structs, 'second_stage_state', secondStagePool)

         call mpas_pool_get_array(LTSPool, 'normalVelocityMPIHaloUpdate', normalVelocityMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'layerThicknessMPIHaloUpdate', layerThicknessMPIHaloUpdate)
         call mpas_pool_get_array(LTSPool, 'tracersGroupMPIHaloUpdate', tracersGroupMPIHaloUpdate) !  TODO this will likely have to be modified!!!!!!!!!!!

         call mpas_pool_get_array(secondStagePool, 'normalVelocity', normalVelocitySecondStage, 1)
         call mpas_pool_get_array(secondStagePool, 'layerThickness', layerThicknessSecondStage, 1)

         ! TODO: take care of the tracers!!!!!!!!!!!!!!!!!!

         normalVelocitySecondStage(:,:) = normalVelocityMPIHaloUpdate(:,:)
         layerThicknessSecondStage(:,:) = layerThicknessMPIHaloUpdate(:,:)

         block => block % next
      end do
      ! ---  copy back for MPI halo update of the solution

    end subroutine ocn_time_integrator_lts3

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts3_init
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   July 2021
!> \details
!>  This routine integrates one timestep (dt) using an LTS3 time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integration_lts3_init(domain)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set local time stepping (LTS) regions and lists
    !
    ! Output: LTS instances are written
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    type (domain_type), intent(inout) :: domain

    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: statePool, meshPool, tendPool, LTSPool
    integer, dimension(:), allocatable :: isLTSRegionEdgeAssigned
    integer :: nLTSHalos
    integer :: i, iCell, iEdge, iVertex, iRegion, iHalo, iHaloP1, iCellHalo, cell1, cell2, nLTSHalosCopy, haveFineCounter, haveCoarseCounter, moreCellsOnInterface
    integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve, haveFine, haveCoarse, haveInterface, doExtraLayers
    real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, areaCell
    integer, dimension(:), pointer :: LTSRegion, LTSRegionLocal, nEdgesOnCell, minMaxLTSRegion, minMaxAfterHaloUpdate, indexToCellID
    integer, dimension(:,:), pointer :: cellsOnEdge,  edgesOnCell, verticesOnCell, nCellsInLTSHalo, nEdgesInLTSHalo
    integer, dimension(:,:,:), pointer :: cellsInLTSHalo, edgesInLTSHalo
    integer, dimension(:), pointer :: mpiRankCell, blockPartitionCell
    real (kind=RKIND) :: latPoint, lonPoint, distLTS

    nLTSHalos = 4

    ! set to 1 if need more cells in interface layers (for load balancing)
    moreCellsOnInterface = 1

    nLTSHalosCopy = nLTSHalos
    if (moreCellsOnInterface == 1) then
       nLTSHalosCopy = 9 !the extra layers you get are actually nLTSHalosCopy - 3 
    end if

    !NOTE:
    ! the word "halo" in LTS instances is not related to the MPI halo.

    !NOTE:
    ! for performance, LTS should to be run with either one of the following setups:
    ! 1) N MPI partitions and 3N blocks, each blocks owns an LTS region (fine, coarse or interface)
    ! 2) N MPI partitions and N blocks,  each blocks owns an LTS region (fine, coarse or interface)
    ! the first case is better theoretically because it optimizes the parallel efficiency by minimizing the wait time among processors.
    ! There are python scripts to modify the graph.info and the graph.info.part.X files to achieve 1) or 2).
    ! Nevertheless, this code is written so that LTS can run with virtually any graph.info file.

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'mpiRankCell', mpiRankCell)
       call mpas_pool_get_array(meshPool, 'blockPartitionCell', blockPartitionCell)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

       call mpas_pool_get_array(LTSPool, 'LTSRegion', LTSRegion)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
       call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)
       call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       ! the loop below initializes the LTS regions for the current processor (or block)
       ! LTS Regions code:
       ! 1 = fine
       ! 2 = coarse
       ! 3 = interface layer 1
       ! 4 = interface layer 2
       ! 5 = fine (to advance when doing 1st stage on interface)

       LTSRegion(:) = 1
       do iCell = 1, nCells
          !if (2.0 * sqrt(areaCell(iCell) / pii) .ge. 300000) then
          !   LTSRegion(iCell) = 2
          !end if

          !latPoint = pii / 6.0
          !lonPoint = 3.0 * pii * 0.5
          latPoint = pii * 0.5
          lonPoint = 135 * pii / 180
          distLTS = sphere_distance(latCell(iCell), lonCell(iCell), latPoint, lonPoint, 1.0)
          ! NOTE: the number on the right hand side of the inequality dictates
          ! where is the interface between coarse and fine on the mesh
          !if (distLTS > 1.2) then  !this was for the convergence test for test 5
          !if (distLTS > 0.35) then   !this is for the CPU ratio test for test 5
          if (distLTS > 0.7) then   !this is for the CPU ratio test with multiple interface layers for test 5
          !if (distLTS > 0.2) then   !this is for the speed up test for test 5
          !if (distLTS > 0.55) then   !this is for the speed up test with multiple interface layers for test 5 (this goes with 38 or 40 layers)
          !if (distLTS > 0.3) then   !this is for the speed up test with multiple interface layers for test 5 (this goes with 18 or 20 layers)
             LTSRegion(iCell) = 2
          end if
          !BEGIN: test initializaton
          !if (latCell(iCell) < 0.) then
          !   LTSRegion(iCell) = 1
          !else
          !   LTSRegion(iCell) = 2
          !end if
          !END: test initialization
       end do

       nCellsInLTSHalo(:,:) = 0
       nEdgesInLTSHalo(:,:) = 0
       LTSRegionLocal(:) = LTSRegion(:)

       haveFine = 1
       haveCoarse = 1
       haveInterface = 1
       doExtraLayers = 0

       ! this is to visualize MPI partitions on Paraview
       mpiRankCell(1:nCells) = domain % dminfo % my_proc_id

       ! this is to visualize block partitions on Paraview
       do iCell = 1, nCells
          blockPartitionCell(iCell) = block % blockID
       end do

       minMaxLTSRegion(1) = minval(LTSRegionLocal(1:nCells))
       minMaxLTSRegion(2) = maxval(LTSRegionLocal(1:nCells))

       if (minMaxLTSRegion(2) - minMaxLTSRegion(1) > 1) then
          call mpas_log_write('minMaxLTSRegion(2) - minMaxLTSRegion(1) > 1 not allowed. Change LTSRegion initialization in ocn_time_integrator_lts3_init subroutine.')
          stop
       end if

       if (minMaxLTSRegion(2) - minMaxLTSRegion(1) == 1) then  !this means there are two regions, a fine and a coarse

          doExtraLayers = 1

          do iEdge = 1, nEdges
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
             if ((indexToCellID(cell1) .ne. 0) .and. (indexToCellID(cell2) .ne. 0)) then !this is to take out of consideration land cells
                ! here we write the first layers of the LTS halo (interface layer 1 and interface layer 2) 
                if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2)) then
                   if (cell1 .le. nCells) then
                      LTSRegionLocal(cell1) = minMaxLTSRegion(1) + 2
                   end if
                   if (cell2 .le. nCells) then
                      LTSRegionLocal(cell2) = minMaxLTSRegion(2) + 2
                   end if
                else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2)) then
                   if (cell2 .le. nCells) then
                      LTSRegionLocal(cell2) = minMaxLTSRegion(1) + 2
                   end if
                   if (cell1 .le. nCells) then
                      LTSRegionLocal(cell1) = minMaxLTSRegion(2) + 2
                   end if
                ! the next if statements are there because above we are possibly changing two LTSRegionLocal at the time, 
                ! so we have to take into consideration the fact that an edge might be shared by two elements, one of which had its LTSRegionLocal already changed
                else if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) + 2 .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2)) then
                   if (cell2 .le. nCells) then
                      LTSRegionLocal(cell2) = minMaxLTSRegion(2) + 2
                   end if
                else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) + 2 .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2)) then
                   if (cell1 .le. nCells) then
                      LTSRegionLocal(cell1) = minMaxLTSRegion(2) + 2
                   end if
                else if (LTSRegionLocal(cell1) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell2) == minMaxLTSRegion(2) + 2) then
                   if ( cell1 .le. nCells) then
                      LTSRegionLocal(cell1) = minMaxLTSRegion(1) + 2
                   end if
                else if (LTSRegionLocal(cell2) == minMaxLTSRegion(1) .and. LTSRegionLocal(cell1) == minMaxLTSRegion(2) + 2) then
                   if (cell2 .le. nCells) then
                      LTSRegionLocal(cell2) = minMaxLTSRegion(1) + 2
                   end if
                end if
             end if
          end do

       end if

      block => block % next
    end do

    call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')

    ! this below is to see if a processor that started with only fine or only coarse needs to fill out some extra LTS halo layers
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

       call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
       call mpas_pool_get_array(LTSPool, 'minMaxAfterHaloUpdate', minMaxAfterHaloUpdate)
       call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)

       minMaxAfterHaloUpdate(1) = minval(LTSRegionLocal(1:nCells))
       minMaxAfterHaloUpdate(2) = maxval(LTSRegionLocal(1:nCells))

       if (minMaxLTSRegion(2) - minMaxLTSRegion(1) == 0) then

          if (minMaxAfterHaloUpdate(1) .ne. minMaxAfterHaloUpdate(2)) then
             doExtraLayers = 1
          end if

          if (minMaxAfterHaloUpdate(1) .ne. minMaxLTSRegion(1)) then
             call mpas_log_write('minMaxAfterHaloUpdate(1) not equal to minMaxLTSRegion(1), something is wrong, see ocn_time_integrator_lts3_init subroutine')
          end if
       end if

       block => block % next
    end do

    ! this is a loop to change LTSRegionLocal for the extra LTS regions
    do iRegion = 1,2
       do iHalo = 1, nLTSHalosCopy - 1

          block => domain % blocklist
          do while (associated(block))
             call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

             call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

             call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
             call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
             call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
             call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
             call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
             call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)

             if (doExtraLayers == 1) then

                do iCell = 1, nCells
                   if (LTSRegionLocal(iCell) == (minMaxLTSRegion(iRegion) + 2 * iHalo)) then
                      do i = 1, nEdgesOnCell(iCell)
                         iEdge = edgesOnCell(i, iCell)
                         cell1 = cellsOnEdge(1,iEdge)
                         cell2 = cellsOnEdge(2,iEdge)
                         if (LTSRegionLocal(cell1) == minMaxLTSRegion(iRegion) ) then
                            if (cell1 .le. nCells) then
                               LTSRegionLocal(cell1) = LTSRegionLocal(cell1) +  2 * (iHalo+1)
                            end if
                         else if (LTSRegionLocal(cell2) == minMaxLTSRegion(iRegion)) then
                            if (cell2 .le. nCells) then
                               LTSRegionLocal(cell2) = LTSRegionLocal(cell2) +  2 * (iHalo+1)
                            end if
                         end if
                      end do
                   end if
                end do

             end if

             block => block % next
          end do


          call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')


          ! this below is to see if a processor that started with only fine or only coarse needs to fill out some extra LTS halo layers
          block => domain % blocklist
          do while (associated(block))
             call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
             call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
             call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

             call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
             call mpas_pool_get_array(LTSPool, 'minMaxAfterHaloUpdate', minMaxAfterHaloUpdate)
             call mpas_pool_get_array(LTSPool, 'doExtraLayers', doExtraLayers)
             call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)

             minMaxAfterHaloUpdate(1) = minval(LTSRegionLocal(1:nCells))
             minMaxAfterHaloUpdate(2) = maxval(LTSRegionLocal(1:nCells))

             if (doExtraLayers == 0) then

                if (minMaxAfterHaloUpdate(1) .ne. minMaxAfterHaloUpdate(2)) then
                   doExtraLayers = 1
                end if

                if (minMaxAfterHaloUpdate(1) .ne. minMaxLTSRegion(1)) then
                   call mpas_log_write('minMaxAfterHaloUpdate(1) not equal to minMaxLTSRegion(1), something is wrong, see ocn_time_integrator_lts3_init subroutine')
                  stop
                end if
             end if

             block => block % next
          end do

       end do
    end do


    call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)

       do iCell = 1, nCells
          if( (mod(LTSRegionLocal(iCell),2) == 1) .and. (LTSRegionLocal(iCell) .ne. 1) ) then
             ! if we are here it could be either interface 1 or those two layers of fine we need for the third order
             if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy + 1) then
                LTSRegionLocal(iCell) = 5 !this goes in cellsInLTSHalo(1,3)
             else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy - 1) then
                LTSRegionLocal(iCell) = 5 !this goes in cellsInLTSHalo(1,3)
             else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy - 3) then
                LTSRegionLocal(iCell) = 5 !this goes in cellsInLTSHalo(1,3)
             else
                LTSRegionLocal(iCell) = 3 !this goes in cellsInLTSHalo(1,2)
             end if
          else if ( (mod(LTSRegionLocal(iCell),2) == 0) .and. (LTSRegionLocal(iCell) .ne. 2) ) then
             ! if we are here it could be either interface 2 or those two layers of coarse
             if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy + 2) then
                LTSRegionLocal(iCell) = 2 !this goes in cellsInLTSHalo(2,1)
             else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy) then
                LTSRegionLocal(iCell) = 2 !this goes in cellsInLTSHalo(2,1)
             else if (LTSRegionLocal(iCell) == 2 * nLTSHalosCopy - 2) then
                LTSRegionLocal(iCell) = 2 !this goes in cellsInLTSHalo(2,1)
             else
                LTSRegionLocal(iCell) = 4 !this goes in cellsInLTSHalo(2,2)
             end if
          end if
       end do

       block => block % next
    end do

    call mpas_dmpar_field_halo_exch(domain, 'LTSRegionLocal')

    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
       call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
       call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)

       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(LTSPool, 'LTSRegionLocal', LTSRegionLocal)
       call mpas_pool_get_array(LTSPool, 'minMaxLTSRegion', minMaxLTSRegion)
       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)
       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       ! this is a loop to build the lists of elements in the fine, coarse, and interface regions
       do iCell = 1, nCellsSolve !this loops up to nCellsSolve because in the time stepping we only want to advance the cells owned by the processor (or block)
          do iRegion = 1,2
             if (iRegion == minMaxLTSRegion(iRegion)) then
                if(LTSRegionLocal(iCell) == minMaxLTSRegion(iRegion)) then
                   nCellsInLTSHalo(iRegion,1) = nCellsInLTSHalo(iRegion,1) + 1
                   cellsInLTSHalo(iRegion,1,nCellsInLTSHalo(iRegion,1)) = iCell
                end if
                if(LTSRegionLocal(iCell) == (minMaxLTSRegion(iRegion) + 2) ) then
                   nCellsInLTSHalo(iRegion,2) = nCellsInLTSHalo(iRegion,2) + 1
                   cellsInLTSHalo(iRegion,2,nCellsInLTSHalo(iRegion,2)) = iCell
                end if
             end if
          end do
          if (LTSRegionLocal(iCell) == 5) then
             nCellsInLTSHalo(1,3) = nCellsInLTSHalo(1,3) + 1
             cellsInLTSHalo(1,3,nCellsInLTSHalo(1,3)) = iCell
          end if
       end do

       if ((nCellsInLTSHalo(1,1) == 0) .and. (nCellsInLTSHalo(1,3) == 0)) then
             haveFine = 0
       end if

       if (nCellsInLTSHalo(2,1) == 0) then
             haveCoarse = 0
       end if

       if ((nCellsInLTSHalo(1,2) == 0) .and. (nCellsInLTSHalo(2,2) == 0)) then
          haveInterface = 0
       end if

       ! below we fill out the lists for the edges, according to the
       ! LTSRegionLocal that have been assigned to the cells.
       ! we move from the fine to the coarse (i.e. from the fine to the nearest
       ! LTS region in the direction of the coarse).
       ! Note that edges shared between cells of different LTS regions are owned
       ! by the cell in the LTS region closest to the fine region,
       ! see Figure 3 in "Conservative explicit local time-stepping schemes for
       ! the shallow water equations" by Hoang et al. (halo edges however
       ! are owned by whatever processor they are initially assigned to)

       allocate(isLTSRegionEdgeAssigned(nEdgesSolve))
       isLTSRegionEdgeAssigned(:) = 0

       do iCell = 1, nCellsInLTSHalo(1,1)
          do i = 1, nEdgesOnCell(cellsInLTSHalo(1,1,iCell))
             iEdge = edgesOnCell(i,cellsInLTSHalo(1,1,iCell))
             if (iEdge .le. nEdgesSolve) then
                if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                   nEdgesInLTSHalo(1,1) = nEdgesInLTSHalo(1,1) + 1
                   edgesInLTSHalo(1,1, nEdgesInLTSHalo(1,1)) = iEdge
                   isLTSRegionEdgeAssigned(iEdge) = 1
                end if
             end if
          end do
       end do

       do iHalo = 1, nLTSHalos
          do iCell = 1, nCellsInLTSHalo(1, nLTSHalos - iHalo + 2)
             do i = 1, nEdgesOnCell(cellsInLTSHalo(1, nLTSHalos- iHalo + 2, iCell))
                iEdge = edgesOnCell(i,cellsInLTSHalo(1, nLTSHalos - iHalo + 2, iCell))
                if (iEdge .le. nEdgesSolve) then
                   if (isLTSRegionEdgeAssigned(iEdge) == 0) then 
                      nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2) = nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2) + 1
                      edgesInLTSHalo(1, nLTSHalos - iHalo + 2, nEdgesInLTSHalo(1, nLTSHalos - iHalo + 2)) = iEdge
                      isLTSRegionEdgeAssigned(iEdge) = 1
                   end if
                end if
             end do
          end do
       end do

       do iHalo = 1, nLTSHalos
          do iCell = 1, nCellsInLTSHalo(2, iHalo+1)
             do i = 1, nEdgesOnCell(cellsInLTSHalo(2,iHalo+1,iCell))
                iEdge = edgesOnCell(i,cellsInLTSHalo(2,iHalo+1,iCell))
                if (iEdge .le. nEdgesSolve) then
                   if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                      nEdgesInLTSHalo(2, iHalo+1) = nEdgesInLTSHalo(2, iHalo+1) + 1
                      edgesInLTSHalo(2, iHalo+1,  nEdgesInLTSHalo(2, iHalo+1)) = iEdge
                      isLTSRegionEdgeAssigned(iEdge) = 1
                   end if
                end if
             end do
          end do
       end do

       do iCell = 1, nCellsInLTSHalo(2,1)
          do i = 1, nEdgesOnCell(cellsInLTSHalo(2,1,iCell))
             iEdge = edgesOnCell(i,cellsInLTSHalo(2,1,iCell))
             if (iEdge .le. nEdgesSolve) then
                if(isLTSRegionEdgeAssigned(iEdge) == 0) then
                   nEdgesInLTSHalo(2,1) = nEdgesInLTSHalo(2,1) + 1
                   edgesInLTSHalo(2,1,nEdgesInLTSHalo(2,1)) = iEdge
                   isLTSRegionEdgeAssigned(iEdge) = 1
                end if
             end if
          end do
       end do

       deallocate(isLTSRegionEdgeAssigned)

       block => block % next
    end do

    ! BEGIN: CHECK TO SEE IF THE INTERFACE BLOCKS HAVE ENOUGH CELLS
    block => domain % blocklist
    do while (associated(block))

       call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

       call mpas_pool_get_array(LTSPool, 'haveFine', haveFine)
       call mpas_pool_get_array(LTSPool, 'haveCoarse', haveCoarse)
       call mpas_pool_get_array(LTSPool, 'haveInterface', haveInterface)

       call mpas_pool_get_array(LTSPool, 'cellsInLTSHalo', cellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nCellsInLTSHalo', nCellsInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'edgesInLTSHalo', edgesInLTSHalo)
       call mpas_pool_get_array(LTSPool, 'nEdgesInLTSHalo', nEdgesInLTSHalo)

       !if(block % blockID == 5) then
       !if(haveInterface  == 1) then
       !if (domain % dminfo % my_proc_id == 0) then
       !   print*, block % blockID
       !   print*, 'haveFine=', haveFine, 'haveCoarse=', haveCoarse,
       !   'haveInterface=', haveInterface
       !   print*, 'nCellsInLTSHalo(1,1)=', nCellsInLTSHalo(1,1),
       !   'nCellsInLTSHalo(1,3)=', nCellsInLTSHalo(1,3),
       !   'nCellsInLTSHalo(1,4)=', nCellsInLTSHalo(1,4)
       !   print*, 'nCellsInLTSHalo(2,1)=', nCellsInLTSHalo(2,1),
       !   'nCellsInLTSHalo(2,3)=', nCellsInLTSHalo(2,3),
       !   'nCellsInLTSHalo(2,4)=', nCellsInLTSHalo(2,4)
       !   print*, 'nCellsInLTSHalo(1,2)=', nCellsInLTSHalo(1,2),
       !   'nCellsInLTSHalo(2,2)=', nCellsInLTSHalo(2,2)
       !   print*, 'nEdgessInLTSHalo(1,1)=', nEdgesInLTSHalo(1,1),
       !   'nEdgesInLTSHalo(1,3)=', nEdgesInLTSHalo(1,3),
       !   'nEdgesInLTSHalo(1,4)=', nEdgesInLTSHalo(1,4)
       !   print*, 'nEdgessInLTSHalo(2,1)=', nEdgesInLTSHalo(2,1),
       !   'nEdgesInLTSHalo(2,3)=', nEdgesInLTSHalo(2,3),
       !   'nEdgesInLTSHalo(2,4)=', nEdgesInLTSHalo(2,4)
       !   print*, 'nEdgessInLTSHalo(1,2)=', nEdgesInLTSHalo(1,2),
       !   'nEdgesInLTSHalo(2,2)=', nEdgesInLTSHalo(2,2)
       !end if

       !TO ERASE
       !if(haveFine == 1) then
       !   if (nCellsInLTSHalo(1,1) + nCellsInLTSHalo(1,3) + nCellsInLTSHalo(1,4)
       !    < 100) then
       !   !if (nCellsInLTSHalo(1,1) < 100) then
       !      print*, 'WARNING: fine block', block % blockID, 'has',
       !      nCellsInLTSHalo(1,1) + nCellsInLTSHalo(1,3) + nCellsInLTSHalo(1,4),
       !      'cells'
       !      !print*, 'WARNING: fine block', block % blockID, 'has',
       !      nCellsInLTSHalo(1,1), 'cells'
       !   end if
       !end if
       !END

       if(haveInterface == 1) then
          if (nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2) < 100) then
             print*, 'WARNING: interface block', block % blockID, 'has', nCellsInLTSHalo(1,2) + nCellsInLTSHalo(2,2), 'cells'
          end if
       end if

       block => block % next
    end do

    ! END: CHECK TO SEE IF THE INTERFACE BLOCKS HAVE ENOUGH CELLS

    end subroutine ocn_time_integration_lts3_init

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on
   ! a sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

end module ocn_time_integration_lts3
