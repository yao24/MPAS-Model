










!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver





  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar
  use mpas_timer
  use mpas_log, only: mpas_log_write, mpas_log_info

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver, &
       seaice_run_velocity_solver

  ! strain scheme type
  integer :: &
       strainSchemeType

  integer, parameter :: &
       WEAK_STRAIN_SCHEME = 1, &
       VARIATIONAL_STRAIN_SCHEME = 2

  logical :: &
       averageVariationalStrains

  ! stress divergence scheme type
  integer :: &
       stressDivergenceSchemeType

  integer, parameter :: &
       WEAK_STRESS_DIVERGENCE_SCHEME = 1, &
       VARIATIONAL_STRESS_DIVERGENCE_SCHEME = 2

  ! ocean stress type
  integer :: &
       oceanStressType

  integer, parameter :: &
       QUADRATIC_OCEAN_STRESS = 1, &
       LINEAR_OCEAN_STRESS = 2

  ! velocity solver constants
  real(kind=RKIND), parameter, private :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       seaiceAreaMinimum = 0.001_RKIND, &
       seaiceMassMinimum = 0.01_RKIND

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver(&
       domain)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         seaice_init_evp

    use seaice_velocity_solver_weak, only: &
         seaice_init_velocity_solver_weak

    use seaice_velocity_solver_variational, only: &
         seaice_init_velocity_solver_variational

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_strain_scheme, &
         config_stress_divergence_scheme, &
         config_variational_basis, &
         config_variational_denominator_type, &
         config_wachspress_integration_type, &
         config_ocean_stress_type

    logical, pointer :: &
         config_use_velocity_solver, &
         config_average_variational_strain, &
         config_rotate_cartesian_grid, &
         config_include_metric_terms, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_use_c_grid, &
         config_use_projected_b_grid

    type (MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         velocitySolver, &
         velocity_weak, &
         velocity_variational, &
         velocity_pwl

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep, &
         elasticTimeStep, &
         config_dt

    integer, pointer :: &
         config_dynamics_subcycle_number, &
         config_elastic_subcycle_number, &
         config_wachspress_integration_order

    integer :: &
         ierr

    ! set up the dynamically locked cell mask
    call dynamically_locked_cell_mask(domain)

    ! set timesteps even with velocity turned off
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
       call MPAS_pool_get_config(block % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)
       call MPAS_pool_get_array(velocitySolver, "elasticTimeStep", elasticTimeStep)

       dynamicsTimeStep = config_dt / real(config_dynamics_subcycle_number,RKIND)

       elasticTimeStep = dynamicsTimeStep / real(config_elastic_subcycle_number,RKIND)

       block => block % next
    enddo

    ! check if we initialize velocity solver
    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then

       ! options
       call MPAS_pool_get_config(domain % configs, "config_strain_scheme", config_strain_scheme)
       if (trim(config_strain_scheme) == "weak") then
          strainSchemeType = WEAK_STRAIN_SCHEME
       else if (trim(config_strain_scheme) == "variational") then
          strainSchemeType = VARIATIONAL_STRAIN_SCHEME
       else
          call MPAS_log_write("config_strain_scheme unknown: "//trim(config_strain_scheme), MPAS_LOG_CRIT)
       endif

       call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
       if (trim(config_stress_divergence_scheme) == "weak") then
          stressDivergenceSchemeType = WEAK_STRESS_DIVERGENCE_SCHEME
       else if (trim(config_stress_divergence_scheme) == "variational") then
          stressDivergenceSchemeType = VARIATIONAL_STRESS_DIVERGENCE_SCHEME
       else
          call MPAS_log_write("config_stress_divergence_scheme unknown: "//trim(config_stress_divergence_scheme), MPAS_LOG_CRIT)
       endif

       call MPAS_pool_get_config(domain % configs, "config_ocean_stress_type", config_ocean_stress_type)
       if (trim(config_ocean_stress_type) == "quadratic") then
          oceanStressType = QUADRATIC_OCEAN_STRESS
       else if (trim(config_ocean_stress_type) == "linear") then
          oceanStressType = LINEAR_OCEAN_STRESS
       else
          call MPAS_log_write("config_ocean_stress_type unknown: "//trim(config_ocean_stress_type), MPAS_LOG_CRIT)
       endif

       if (strainSchemeType           == VARIATIONAL_STRAIN_SCHEME .and. &
           stressDivergenceSchemeType == WEAK_STRESS_DIVERGENCE_SCHEME) then
          call MPAS_log_write("Cannot have variational strain scheme with variational stress divergence scheme", MPAS_LOG_CRIT)
       endif

       call MPAS_pool_get_config(domain % configs, "config_average_variational_strain", config_average_variational_strain)
       averageVariationalStrains = config_average_variational_strain

       ! initialize the evp solver
       call seaice_init_evp(domain)

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
          call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocity_weak)
          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocity_variational)

          call MPAS_pool_get_config(block % configs, "config_variational_basis", config_variational_basis)
          call MPAS_pool_get_config(block % configs, "config_variational_denominator_type", config_variational_denominator_type)
          call MPAS_pool_get_config(block % configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
          call MPAS_pool_get_config(block % configs, "config_include_metric_terms", config_include_metric_terms)
          call MPAS_pool_get_config(block % configs, "config_wachspress_integration_type", config_wachspress_integration_type)
          call MPAS_pool_get_config(block % configs, "config_wachspress_integration_order", config_wachspress_integration_order)

          ! init solvers
          if (strainSchemeType           == WEAK_STRAIN_SCHEME .or. &
              stressDivergenceSchemeType == WEAK_STRESS_DIVERGENCE_SCHEME) then

             call seaice_init_velocity_solver_weak(&
                  mesh, &
                  boundary, &
                  velocity_weak, &
                  config_rotate_cartesian_grid)

          endif

          if (strainSchemeType           == VARIATIONAL_STRAIN_SCHEME .or. &
              stressDivergenceSchemeType == VARIATIONAL_STRESS_DIVERGENCE_SCHEME) then

             call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)
             call MPAS_pool_get_config(block % configs, "config_use_projected_b_grid", config_use_projected_b_grid)

             call seaice_init_velocity_solver_variational(&
                  mesh, &
                  velocity_variational, &
                  boundary, &
                  config_rotate_cartesian_grid, &
                  config_include_metric_terms, &
                  config_variational_basis, &
                  config_variational_denominator_type, &
                  config_wachspress_integration_type, &
                  config_wachspress_integration_order, &
                  config_use_c_grid, &
                  config_use_projected_b_grid)

          endif

          block => block % next
       enddo

    endif

    ! initialize the land ice shelve mask
    call init_ice_shelve_vertex_mask(domain)

    ! prep for aggregated halo exchanges
    call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
    if (config_aggregate_halo_exch) then

       ! create the velocity aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'velocityHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'velocityHaloExchangeGroup', 'uVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add uVelocity to velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'velocityHaloExchangeGroup', 'vVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add vVelocity to velocityHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the iceAreaTotalMassExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'iceAreaTotalMassExchangeGroup', 'iceAreaCellInitial', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add iceAreaCellInitial to iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'iceAreaTotalMassExchangeGroup', 'totalMassCell', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add totalMassCell to iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the solveVelocityExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'solveVelocityExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create solveVelocityExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'solveVelocityExchangeGroup', 'solveVelocity', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add solveVelocity to solveVelocityExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the icePressureExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'icePressureExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create icePressureExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'icePressureExchangeGroup', 'icePressure', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add icePressure to icePressureExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the airStressHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'airStressHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'airStressHaloExchangeGroup', 'airStressCellU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add airStressCellU to airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'airStressHaloExchangeGroup', 'airStressCellV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add airStressCellV to airStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the seaSurfaceTiltHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'seaSurfaceTiltHaloExchangeGroup', 'seaSurfaceTiltU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add seaSurfaceTiltU to seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'seaSurfaceTiltHaloExchangeGroup', 'seaSurfaceTiltV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add seaSurfaceTiltV to seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! create the oceanStressHaloExchangeGroup aggregated halo exchange
       call mpas_dmpar_exch_group_create(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to create oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'oceanStressHaloExchangeGroup', 'oceanStressU', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add oceanStressU to oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif
       call mpas_dmpar_exch_group_add_field(domain, 'oceanStressHaloExchangeGroup', 'oceanStressV', iErr=ierr)
       if (ierr /= MPAS_DMPAR_NOERR) then
          call MPAS_log_write("failure to add oceanStressV to oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
       endif

       ! build reusable buffers
       call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
       if (config_reuse_halo_exch) then

          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "velocityHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "iceAreaTotalMassExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "solveVelocityExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "icePressureExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for icePressureExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "airStressHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "seaSurfaceTiltHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
          endif
          call mpas_dmpar_exch_group_build_reusable_buffers(domain, "oceanStressHaloExchangeGroup", iErr=ierr)
          if (ierr /= MPAS_DMPAR_NOERR) then
             call MPAS_log_write("failure to build reusable buffers for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
          endif

       endif

    endif ! config_aggregate_halo_exch

  end subroutine seaice_init_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  dynamically_locked_cell_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 6th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine dynamically_locked_cell_mask(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool, &
         boundaryPool

    integer, dimension(:), pointer :: &
         dynamicallyLockedCellsMask, &
         nEdgesOnCell, &
         interiorVertex

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    integer, pointer :: &
         nCells

    integer :: &
         iCell, &
         iVertexOnCell, &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)

       call MPAS_pool_get_array(velocitySolverPool, "dynamicallyLockedCellsMask", dynamicallyLockedCellsMask)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "verticesOnCell", verticesOnCell)

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       do iCell = 1, nCells

          dynamicallyLockedCellsMask(iCell) = 1

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             iVertex = verticesOnCell(iVertexOnCell,iCell)

             if (interiorVertex(iVertex) == 1) then
                dynamicallyLockedCellsMask(iCell) = 0
                exit
             endif

          enddo ! iVertexOnCell

       enddo ! iCell

       block => block % next
    enddo

  end subroutine dynamically_locked_cell_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_shelve_vertex_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_shelve_vertex_mask(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         oceanCouplingPool, &
         meshPool

    integer, dimension(:), pointer :: &
         landIceMaskVertex, &
         landIceMaskEdge, &
         landIceMask

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         cellsOnEdge

    integer, pointer :: &
         vertexDegree, &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iCell, &
         iVertex, &
         iEdge, &
         iCellOnVertex

    logical, pointer :: &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)
       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskVertex", landIceMaskVertex)
       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskEdge", landIceMaskEdge)

       call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)
       call MPAS_pool_get_array(meshPool, "cellsOnEdge", cellsOnEdge)

       call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)
       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          landIceMaskVertex(:) = 0

          do iVertex = 1, nVerticesSolve

             do iCellOnVertex = 1, vertexDegree

                iCell = cellsOnVertex(iCellOnVertex,iVertex)

                if (landIceMask(iCell) == 1) then

                   landIceMaskVertex(iVertex) = 1

                endif

             enddo ! iCellOnVertex

          enddo ! iVertex

       else 

          landIceMaskEdge(:) = 0

          do iEdge = 1, nEdgesSolve

             if(landIceMask(cellsOnEdge(1, iEdge)) == 1 .or. landIceMask(cellsOnEdge(2, iEdge)) == 1) then

                landIceMaskEdge(iEdge) = 1  

             end if

          end do ! iEdge

       end if

       block => block % next
    enddo

  end subroutine init_ice_shelve_vertex_mask

!-----------------------------------------------------------------------
! Time stepping
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_run_velocity_solver(domain, clock)!{{{

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_velocity_solver

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_velocity_solver", config_use_velocity_solver)

    if (config_use_velocity_solver) then

       ! pre subcycle
       call mpas_timer_start("Velocity solver pre-cycle")
       call velocity_solver_pre_subcycle(domain)
       call mpas_timer_stop("Velocity solver pre-cycle")

       ! subcycle the dynamics
       call mpas_timer_start("Velocity solver sub-cycle")
       call subcycle_velocity_solver(domain, clock)
       call mpas_timer_stop("Velocity solver sub-cycle")

       ! post subcycle
       call mpas_timer_start("Velocity solver post-cycle")
       call velocity_solver_post_subcycle(domain)
       call mpas_timer_stop("Velocity solver post-cycle")

    endif ! config_use_velocity_solver

  end subroutine seaice_run_velocity_solver!}}}

!-----------------------------------------------------------------------
! Pre sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_solver_pre_subcycle
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date January 13th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine velocity_solver_pre_subcycle(domain)

    use seaice_mesh_pool, only: &
         seaice_mesh_pool_update

    type(domain_type), intent(inout) :: &
         domain

    ! aggregate categories for area and volume into total mass
    call mpas_timer_start("agg mass and area")
    call aggregate_mass_and_area(domain)
    call mpas_timer_stop("agg mass and area")

    ! calculate computational masks
    call mpas_timer_start("calc masks")
    call calculation_masks(domain)
    call mpas_timer_stop("calc masks")

    ! set new ice velocities to ocean velocity
    call mpas_timer_start("new ice vel")
    call new_ice_velocities(domain)
    call mpas_timer_stop("new ice vel")

    ! calculate the ice strength
    call mpas_timer_start("ice strength")
    call ice_strength(domain)
    call mpas_timer_stop("ice strength")

    ! calculate air stress
    call mpas_timer_start("air stress")
    call air_stress(domain)
    call mpas_timer_stop("air stress")

    ! calculate the coriolis force coefficient
    call mpas_timer_start("coliolis force coef")
    call coriolis_force_coefficient(domain)
    call mpas_timer_stop("coliolis force coef")

    ! calculate the ocean stress
    call mpas_timer_start("ocean stress")
    call ocean_stress(domain)
    call mpas_timer_stop("ocean stress")

    ! calculate the surface tilt force
    call mpas_timer_start("surface tilt")
    call surface_tilt(domain)
    call mpas_timer_stop("surface tilt")

    ! initialize subcycle variables
    call mpas_timer_start("init subcycle var")
    call init_subcycle_variables(domain)
    call mpas_timer_stop("init subcycle var")

    ! update mesh pool variables
    call mpas_timer_start("update mesh pool")
    call seaice_mesh_pool_update(domain)
    call mpas_timer_stop("update mesh pool")

  end subroutine velocity_solver_pre_subcycle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  aggregate_mass_and_area
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine aggregate_mass_and_area(domain)!{{{

    use seaice_constants, only: &
         seaiceDensityIce, &
         seaiceDensitySnow

    type(domain_type), intent(in) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         tracersPool, &
         tracersAggregatePool, &
         icestatePool

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), dimension(:), pointer :: &
         totalMassCell, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
       call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracersAggregatePool, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)

       do iCell = 1, nCellsSolve

          iceAreaCell(iCell)    = sum(iceAreaCategory(1,:,iCell))
          iceVolumeCell(iCell)  = sum(iceVolumeCategory(1,:,iCell))
          snowVolumeCell(iCell) = sum(snowVolumeCategory(1,:,iCell))

          totalMassCell(iCell)  = iceVolumeCell(iCell)  * seaiceDensityIce + &
                                  snowVolumeCell(iCell) * seaiceDensitySnow

       enddo ! iCell

       block => block % next
    enddo

  end subroutine aggregate_mass_and_area!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculation_masks
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine calculation_masks(domain)

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         tracersAggregatePool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaEdge, &
         iceAreaCellInitial, &
         iceAreaCell, &
         totalMassVertex, &
         totalMassEdge, &
         totalMassCell

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_calc_velocity_masks, &
         config_use_c_grid

    integer, pointer :: &
         nEdgesSolve

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    integer :: &
         cell1, &
         cell2, &
         iEdge, &
         ierr

    ! set initial ice area if not have column physics
    call MPAS_pool_get_config(domain % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", &
                                                 config_use_column_vertical_thermodynamics)

    if (.not. config_use_column_package .or. &
         (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
          call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
          call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
          call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)

          iceAreaCellInitial = iceAreaCell

          block => block % next
       end do

    endif

    ! halo exchange of initial ice area of cell
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("ice area halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'iceAreaTotalMassExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for iceAreaTotalMassExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'iceAreaCellInitial')
          call MPAS_dmpar_field_halo_exch(domain, 'totalMassCell')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("ice area halo")

    call seaice_load_balance_timers(domain, "vel prep after")

    ! interpolate area and mass from cells to vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "iceAreaEdge", iceAreaEdge)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)
       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)

       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)
       call MPAS_pool_get_array(meshPool, "cellsOnEdge", cellsOnEdge)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid ) then

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               iceAreaVertex, &
               iceAreaCellInitial)

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               totalMassVertex, &
               totalMassCell)

       else

          do iEdge = 1, nEdgesSolve

             cell1 = cellsOnEdge(1, iEdge)
             cell2 = cellsOnEdge(2, iEdge)

             iceAreaEdge(iEdge) = ( iceAreaCellInitial(cell1) + iceAreaCellInitial(cell2) ) * 0.5_RKIND
             totalMassEdge(iEdge) = ( totalMassCell(cell1) + totalMassCell(cell2) ) * 0.5_RKIND

          end do

       end if

       block => block % next
    end do

    ! calculate computational masks
    call MPAS_pool_get_config(domain % configs, "config_calc_velocity_masks", config_calc_velocity_masks)
    if (config_calc_velocity_masks) then
       call stress_calculation_mask(domain)
       call velocity_calculation_mask(domain)
    endif

    ! halo exchange velocity mask
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("velocity mask halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'solveVelocityExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'solveVelocityExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for solveVelocityExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'solveVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'solveVelocityCGrid')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("velocity mask halo")

    call seaice_load_balance_timers(domain, "vel prep after")

  end subroutine calculation_masks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  stress_calculation_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine stress_calculation_mask(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool, &
         icestatePool, &
         oceanCouplingPool

    integer, dimension(:), pointer :: &
         solveStress, &
         solveStressTri

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCellInitial, &
         totalMassCell, &
         iceAreaVertex, &
         totalMassVertex

    integer :: &
         iCell, &
         iCellOnCell, &
         iCellNeighbour, &
         iVertex, &
         iCellOnVertex

    integer, pointer :: &
         nCells, &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell, &
         landIceMask, &
         landIceMaskVertex

    integer, dimension(:,:), pointer :: &
         cellsOnCell, &
         cellsOnVertex

    logical, pointer :: &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

       call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
       call MPAS_pool_get_array(meshPool, "cellsOnCell", cellsOnCell)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestatePool, "totalMassCell", totalMassCell)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)

       do iCell = 1, nCells

          solveStress(iCell) = 0

          if (iceAreaCellInitial(iCell) > seaiceAreaMinimum .and. &
              totalMassCell(iCell) > seaiceMassMinimum .and. &
              landIceMask(iCell) == 0) then

             ! this cell has sufficient ice
             solveStress(iCell) = 1

          else

             ! test neighbouring cells to see if have sufficient ice
             do iCellOnCell = 1, nEdgesOnCell(iCell)

                iCellNeighbour = cellsOnCell(iCellOnCell,iCell)

                if (iceAreaCellInitial(iCellNeighbour) > seaiceAreaMinimum .and. &
                    totalMassCell(iCellNeighbour) > seaiceMassMinimum .and. &
                    landIceMask(iCellNeighbour) == 0) then

                   solveStress(iCell) = 1
                   exit

                endif

             enddo ! iCellOnCell

          endif

       enddo ! iCell

       do iCell = nCells+1, nCells

          solveStress(iCell) = 0

       enddo ! iCell

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (config_use_c_grid) then

          call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)
          call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)
          call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)

          call MPAS_pool_get_array(velocitySolverPool, "solveStressTri", solveStressTri)

          solveStressTri(:) = 0

          do iVertex = 1, nVertices

             ! if a vertex belongs to a cell that has enough ice, then solveStressTri is set to true
             do iCellOnVertex = 1, vertexDegree
                if (iceAreaCellInitial(cellsOnVertex(iCellOnVertex,iVertex)) > seaiceAreaMinimum .and. &
                   totalMassCell(cellsOnVertex(iCellOnVertex,iVertex)) > seaiceMassMinimum .and. &
                   landIceMask(cellsOnVertex(iCellOnVertex,iVertex)) == 0) then

                   solveStressTri(iVertex) = 1
                   exit

                end if
             end do

          end do          

       end if
 
       block => block % next
    enddo

  end subroutine stress_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_calculation_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine velocity_calculation_mask(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool, &
         boundaryPool, &
         oceanCouplingPool, &
         meshPool

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaEdge, &
         totalMassVertex, &
         totalMassEdge

    integer, dimension(:), pointer :: &
         interiorVertex, &
         interiorEdge, &
         landIceMaskVertex, &
         landIceMaskEdge

    integer :: &
         iVertex, &
         iEdge, &
         iVertOnEdge

    integer, pointer :: &
         nVerticesSolve, &
         nVertices, &
         nEdgesSolve, &
         nEdges

    logical, pointer :: &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdges", nEdges)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)

       call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(icestatePool, "iceAreaEdge", iceAreaEdge)
       call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)
       call MPAS_pool_get_array(boundaryPool, "interiorEdge", interiorEdge)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskVertex", landIceMaskVertex)
       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskEdge", landIceMaskEdge)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          do iVertex = 1, nVerticesSolve

             solveVelocity(iVertex) = 0

             if (interiorVertex(iVertex) == 1 .and. &
                 landIceMaskVertex(iVertex) == 0 .and. &
                 iceAreaVertex(iVertex) > seaiceAreaMinimum .and. &
                 totalMassVertex(iVertex) > seaiceMassMinimum) then

                ! this vertex has sufficient ice
                solveVelocity(iVertex) = 1

             endif

          end do ! iVertex

          do iVertex = nVerticesSolve+1, nVertices

             solveVelocity(iVertex) = 0

          end do ! iVertex

       else

          solveVelocityCGrid(:) = 0

          do iEdge = 1, nEdgesSolve 

             if (interiorEdge(iEdge) == 1 .and. &
                 landIceMaskEdge(iEdge) == 0 .and. &
                 iceAreaEdge(iEdge) > seaiceAreaMinimum .and. &
                 totalMassEdge(iEdge) > seaiceMassMinimum) then

             ! this edge has sufficient ice
             solveVelocityCGrid(iEdge) = 1

             end if

          end do

       end if

       block => block % next
    end do


  end subroutine velocity_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  new_ice_velocities
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 29th June 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine new_ice_velocities(domain)!{{{

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool, &
         oceanCouplingPool

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocity, &
         vOceanVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uOceanVelocityEdge, &
         vOceanVelocityEdge, &
         uVelocity, &
         vVelocity, &
         uVelocityCGrid, &
         vVelocityCGrid, &
         uVelocityInitial, &
         vVelocityInitial, &
         uVelocityInitialCGrid, &
         vVelocityInitialCGrid, &
         stressDivergenceU, &
         stressDivergenceV, &
         stressDivergenceUCGrid, &
         stressDivergenceVCGrid, &
         oceanStressU, &
         oceanStressV, &
         oceanStressEdgeU, &
         oceanStressEdgeV

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityPrevious, &
         solveVelocityCGrid, &
         solveVelocityPreviousCGrid

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    integer, pointer :: &
         nVerticesSolve, & 
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge, &
         cell1, &
         cell2, &
         ierr

    logical, pointer :: &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)

       call MPAS_pool_get_array(meshPool, "cellsOnEdge", cellsOnEdge)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid", uVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid", vVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitialCGrid", uVelocityInitialCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitialCGrid", vVelocityInitialCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityPrevious", solveVelocityPrevious)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityPreviousCGrid", solveVelocityPreviousCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityEdge", uOceanVelocityEdge)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityEdge", vOceanVelocityEdge)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeU", oceanStressEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeV", oceanStressEdgeV)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGrid", stressDivergenceUCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGrid", stressDivergenceVCGrid)

       call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
       call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          ! interpolate cell ocean velocity to vertex
          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               uOceanVelocityVertex, &
               uOceanVelocity)

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               vOceanVelocityVertex, &
               vOceanVelocity)

          ! set new ice to ocean velocity
          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                if (solveVelocityPrevious(iVertex) == 0) then

                   uVelocity(iVertex) = uOceanVelocityVertex(iVertex)
                   vVelocity(iVertex) = vOceanVelocityVertex(iVertex)

                endif

             else

                uVelocity(iVertex) = 0.0_RKIND
                vVelocity(iVertex) = 0.0_RKIND
                stressDivergenceU(iVertex) = 0.0_RKIND
                stressDivergenceV(iVertex) = 0.0_RKIND
                oceanStressU(iVertex) = 0.0_RKIND
                oceanStressV(iVertex) = 0.0_RKIND

             endif

          enddo ! iVertex

          solveVelocityPrevious = solveVelocity

          uVelocityInitial = uVelocity
          vVelocityInitial = vVelocity

       else 

          do iEdge = 1, nEdgesSolve

             cell1 = cellsOnEdge(1, iEdge)
             cell2 = cellsOnEdge(2, iEdge)
             uOceanVelocityEdge(iEdge) = ( uOceanVelocity(cell1) + uOceanVelocity(cell2) ) * 0.5_RKIND
             vOceanVelocityEdge(iEdge) = ( vOceanVelocity(cell1) + vOceanVelocity(cell2) ) * 0.5_RKIND

             if(solveVelocityCGrid(iEdge) == 1) then

                if(solveVelocityPreviousCGrid(iEdge) == 0 ) then

                   uVelocityCGrid(iEdge) = uOceanVelocityEdge(iEdge)
                   vVelocityCGrid(iEdge) = vOceanVelocityEdge(iEdge)

                end if
         
             else

                uVelocityCGrid(iEdge) = 0.0_RKIND !not sure this line is necessary see init_subcycle_variables routine
                vVelocityCGrid(iEdge) = 0.0_RKIND !not sure this line is necessary see init_subcycle_variables routine
                stressDivergenceUCGrid(iEdge) = 0.0_RKIND !not sure this line is necessary see init_subcycle_variables routine
                stressDivergenceVCGrid(iEdge) = 0.0_RKIND !not sure this line is necessary see init_subcycle_variables routine
                oceanStressEdgeU(iEdge) = 0.0_RKIND
                oceanStressEdgeV(iEdge) = 0.0_RKIND

             end if

          end do ! iEdge

          solveVelocityPreviousCGrid = solveVelocityCGrid

          uVelocityInitialCGrid = uVelocityCGrid
          vVelocityInitialCGrid = vVelocityCGrid

       end if
       block => block % next
    enddo

    ! halo exchange velocities
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("velocity halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocityCGrid')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocityCGrid')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("velocity halo")

    call seaice_load_balance_timers(domain, "vel prep after")

  end subroutine new_ice_velocities!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ice_strength
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ice_strength(domain)

    use ice_colpkg, only: &
         colpkg_ice_strength

    use seaice_constants, only: &
         seaiceIceStrengthConstantHiblerP, &
         seaiceIceStrengthConstantHiblerC

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         tracersAggregatePool, &
         icestatePool, &
         tracersPool

    real(kind=RKIND), dimension(:), pointer:: &
         icePressure, &
         iceAreaCell, &
         iceVolumeCell, &
         openWaterArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory

    integer, dimension(:), pointer :: &
         solveStress

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell, &
         ierr

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)
       call MPAS_pool_get_config(block % configs, "config_use_column_vertical_thermodynamics", &
                                                   config_use_column_vertical_thermodynamics)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)

       call MPAS_pool_get_array(icestatePool, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)

       if (.not. config_use_column_package .or. &
            (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

          do iCell = 1, nCellsSolve

             if (solveStress(iCell) == 1) then

                icePressure(iCell) = seaiceIceStrengthConstantHiblerP * iceVolumeCell(iCell) * &
                     exp(-seaiceIceStrengthConstantHiblerC*(1.0_RKIND-iceAreaCell(iCell)))

             else

                icePressure(iCell) = 0.0_RKIND

             endif ! solveStress

          enddo ! iCell

       else

          do iCell = 1, nCellsSolve

             icePressure(iCell) = 0.0_RKIND

             if (solveStress(iCell) == 1) then

                ! this routine doesnt reset icePressure
                call colpkg_ice_strength(&
                     nCategories, &
                     iceAreaCell(iCell), &
                     iceVolumeCell(iCell), &
                     openWaterArea(iCell), &
                     iceAreaCategory(1,:,iCell), &
                     iceVolumeCategory(1,:,iCell), &
                     icePressure(iCell))

             endif ! solveStress

          enddo ! iCell

       endif

       block => block % next
    enddo

    ! halo exchange ice strength
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("ice strength halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'icePressureExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for icePressureExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'icePressureExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for icePressureExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'icePressure')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("ice strength halo")

    call seaice_load_balance_timers(domain, "vel prep after")

  end subroutine ice_strength

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  air_stress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine air_stress(domain)

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_air_stress, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_use_c_grid

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressVertexU, &
         airStressVertexV, &
         airStressEdgeU, &
         airStressEdgeV, &
         airStressCellU, &
         airStressCellV

    integer :: &
         iEdge, &
         cell1, &
         cell2, &
         ierr

    integer, pointer :: &
         nEdgesSolve

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    ! calculate the air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_vertical_thermodynamics", &
                                                             config_use_column_vertical_thermodynamics)

    if (.not. config_use_column_package .or. &
         (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then
       call constant_air_stress(domain)
    endif

    ! check for no air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_air_stress", config_use_air_stress)
    if (.not. config_use_air_stress) then
       block => domain % blocklist
       do while (associated(block))
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
          airStressCellU = 0.0_RKIND
          airStressCellV = 0.0_RKIND
          block => block % next
       end do
    endif ! .not. config_use_air_stress

    ! halo exchange air stress
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("air stress halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'airStressHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'airStressHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for airStressHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'airStressCellU')
          call MPAS_dmpar_field_halo_exch(domain, 'airStressCellV')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("air stress halo")

    call seaice_load_balance_timers(domain, "vel prep after")

    ! interpolate air stress
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressEdgeU", airStressEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressEdgeV", airStressEdgeV)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if ( .not. config_use_c_grid ) then

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               airStressVertexU, &
               airStressCellU)

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               airStressVertexV, &
               airStressCellV)

       else

          call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve) 
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)   
          call MPAS_pool_get_array(meshPool, "cellsOnEdge", cellsOnEdge)

          do iEdge = 1, nEdgesSolve

             cell1 = cellsOnEdge(1, iEdge)
             cell2 = cellsOnEdge(2, iEdge)
             airStressEdgeU(iEdge) = ( airStressCellU(cell1) + airStressCellU(cell2) ) * 0.5_RKIND
             airStressEdgeV(iEdge) = ( airStressCellV(cell1) + airStressCellV(cell2) ) * 0.5_RKIND

          end do

       end if

       block => block % next
    end do

  end subroutine air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  constant_air_stress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine constant_air_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         atmosCouplingPool, &
         tracersAggregatePool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressCellU, &
         airStressCellV, &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         iceAreaCell

    real(kind=RKIND) :: &
         windSpeed

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "airDensity", airDensity)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

       do iCell = 1, nCellsSolve

          windSpeed = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

          airStressCellU(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * uAirVelocity(iCell) * iceAreaCell(iCell)
          airStressCellV(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * vAirVelocity(iCell) * iceAreaCell(iCell)

       enddo ! iCell

       block => block % next
    end do

  end subroutine constant_air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  coriolis_force_coefficient
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine coriolis_force_coefficient(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         totalMassVertexfVertex, &
         totalMassEdgefEdge, &
         totalMassVertex, &
         totalMassEdge, &
         fVertex, &
         fEdge

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge

    logical, pointer :: &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassEdgefEdge", totalMassEdgefEdge)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)
       call MPAS_pool_get_array(meshPool, "fEdge", fEdge)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)
 
       if ( .not. config_use_c_grid) then

          do iVertex = 1, nVerticesSolve

             totalMassVertexfVertex(iVertex) = totalMassVertex(iVertex) * fVertex(iVertex)

          enddo ! iVertex

       else

          do iEdge = 1, nEdgesSolve

             totalMassEdgefEdge(iEdge) = totalMassEdge(iEdge) * fEdge(iEdge)

          end do

       end if

       block => block % next
    end do

  end subroutine coriolis_force_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocean_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressU, &
         oceanStressV, &
         oceanStressEdgeU, &
         oceanStressEdgeV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uOceanVelocityEdge, &
         vOceanVelocityEdge, &
         fVertex, &
         fEdge

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    logical, pointer :: &
         configUseOceanStress, &
         config_use_c_grid

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeU", oceanStressEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeV", oceanStressEdgeV)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (configUseOceanStress) then

          if (.not. config_use_c_grid) then

             call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

             call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

             do iVertex = 1, nVerticesSolve

                if (solveVelocity(iVertex) == 1) then

                   oceanStressU(iVertex) = uOceanVelocityVertex(iVertex) * cosOceanTurningAngle - &
                                           vOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex))
                   oceanStressV(iVertex) = uOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)) + &
                                           vOceanVelocityVertex(iVertex) * cosOceanTurningAngle

                 else

                   oceanStressU(iVertex) = 0.0_RKIND
                   oceanStressV(iVertex) = 0.0_RKIND

                endif ! solvePoints

             enddo ! iVertex

          else

             call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityEdge", uOceanVelocityEdge)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityEdge", vOceanVelocityEdge)

             call MPAS_pool_get_array(meshPool, "fEdge", fEdge)

             do iEdge = 1, nEdgesSolve

                if (solveVelocityCGrid(iEdge) == 1) then

                   oceanStressEdgeU(iEdge) = uOceanVelocityEdge(iEdge) * cosOceanTurningAngle - &
                                              vOceanVelocityEdge(iEdge) * sinOceanTurningAngle * sign(1.0_RKIND,fEdge(iEdge))
                   oceanStressEdgeV(iEdge) = uOceanVelocityEdge(iEdge) * sinOceanTurningAngle * sign(1.0_RKIND,fEdge(iEdge)) + &
                                              vOceanVelocityEdge(iEdge) * cosOceanTurningAngle

                 else

                   oceanStressEdgeU(iEdge) = 0.0_RKIND
                   oceanStressEdgeV(iEdge) = 0.0_RKIND

                endif ! solvePoints

             enddo ! iEdge

          end if

       else

          if (.not. config_use_c_grid) then

             ! no ocean stress
             oceanStressU = 0.0_RKIND
             oceanStressV = 0.0_RKIND
 
          else

             ! no ocean stress
             oceanStressEdgeU = 0.0_RKIND
             oceanStressEdgeV = 0.0_RKIND

          end if    

       endif

       block => block % next
    end do

  end subroutine ocean_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine surface_tilt(domain)

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         configUseSurfaceTilt, &
         configGeostrophicSurfaceTilt

    call MPAS_pool_get_config(domain % configs, "config_use_surface_tilt", configUseSurfaceTilt)
    call MPAS_pool_get_config(domain % configs, "config_geostrophic_surface_tilt", configGeostrophicSurfaceTilt)

    if (configUseSurfaceTilt) then

       if (configGeostrophicSurfaceTilt) then

          call surface_tilt_geostrophic(domain)

       else

          call surface_tilt_ssh_gradient(domain)

       endif

    else

       call no_surface_tilt(domain)

    endif

  end subroutine surface_tilt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt_geostrophic
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine surface_tilt_geostrophic(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         surfaceTiltForceEdgeU, &
         surfaceTiltForceEdgeV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uOceanVelocityEdge, &
         vOceanVelocityEdge, &
         totalMassVertex, &
         totalMassEdge, &
         fVertex, &
         fEdge

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    logical, pointer :: &
         config_use_c_grid

    integer :: &
         iVertex, &
         iEdge

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeU", surfaceTiltForceEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeV", surfaceTiltForceEdgeV)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)

       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityEdge", uOceanVelocityEdge)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityEdge", vOceanVelocityEdge)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)
       call MPAS_pool_get_array(meshPool, "fEdge", fEdge)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          ! calculate surface tilt force from geostrophic currents
          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then
  
                surfaceTiltForceU(iVertex) = -fVertex(iVertex) * totalMassVertex(iVertex) * vOceanVelocityVertex(iVertex)
                surfaceTiltForceV(iVertex) =  fVertex(iVertex) * totalMassVertex(iVertex) * uOceanVelocityVertex(iVertex)

             else

                surfaceTiltForceU(iVertex) = 0.0_RKIND
                surfaceTiltForceV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

       else 

          ! calculate surface tilt force from geostrophic currents
          do iEdge = 1, nEdgesSolve

             if (solveVelocityCGrid(iEdge) == 1) then

                surfaceTiltForceEdgeU(iEdge) = - fEdge(iEdge) * totalMassEdge(iEdge) * vOceanVelocityEdge(iEdge)
                surfaceTiltForceEdgeV(iEdge) =  fEdge(iEdge) * totalMassEdge(iEdge) * uOceanVelocityEdge(iEdge)

             else

                surfaceTiltForceEdgeU(iEdge) = 0.0_RKIND
                surfaceTiltForceEdgeV(iEdge) = 0.0_RKIND

             endif ! solveVelocityCGrid

          enddo ! iEdge

       end if

       block => block % next
    end do

  end subroutine surface_tilt_geostrophic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  surface_tilt_ssh_gradient
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine surface_tilt_ssh_gradient(domain)

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    use seaice_constants, only: &
         seaiceGravity

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool, &
         oceanCouplingPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         surfaceTiltForceEdgeU, &
         surfaceTiltForceEdgeV, &
         totalMassVertex, &
         totalMassEdge, &
         seaSurfaceTiltU, &
         seaSurfaceTiltV, &
         seaSurfaceTiltVertexU, &
         seaSurfaceTiltVertexV, &
         seaSurfaceTiltEdgeU, &
         seaSurfaceTiltEdgeV

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    integer, dimension(:,:), pointer :: &
         cellsOnEdge

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge, &
         cell1, &
         cell2, &
         ierr

    logical, pointer :: &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_use_c_grid

    ! halo exchange surface tilt
    call seaice_load_balance_timers(domain, "vel prep before")

    call mpas_timer_start("sea surface tilt halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'seaSurfaceTiltHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for seaSurfaceTiltHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltU')
          call MPAS_dmpar_field_halo_exch(domain, 'seaSurfaceTiltV')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("sea surface tilt halo")

    call seaice_load_balance_timers(domain, "vel prep after")

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(meshPool, "cellsOnEdge", cellsOnEdge)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeU", surfaceTiltForceEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeV", surfaceTiltForceEdgeV)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)

       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltU", seaSurfaceTiltU)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltV", seaSurfaceTiltV)

       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexU", seaSurfaceTiltVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexV", seaSurfaceTiltVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltEdgeU", seaSurfaceTiltEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltEdgeV", seaSurfaceTiltEdgeV)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          ! interpolate sea surface tilt from cells to vertices
          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               seaSurfaceTiltVertexU, &
               seaSurfaceTiltU)

          call seaice_interpolate_cell_to_vertex(&
               meshPool, &
               seaSurfaceTiltVertexV, &
               seaSurfaceTiltV)

          ! calculate surface tilt from sea surface gradient
          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                surfaceTiltForceU(iVertex) = -seaiceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexU(iVertex)
                surfaceTiltForceV(iVertex) = -seaiceGravity * totalMassVertex(iVertex) * seaSurfaceTiltVertexV(iVertex)

             else

                surfaceTiltForceU(iVertex) = 0.0_RKIND
                surfaceTiltForceV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

          else 

             do iEdge = 1, nEdgesSolve
                cell1 = cellsOnEdge(1, iEdge)
                cell2 = cellsOnEdge(2, iEdge)
                seaSurfaceTiltEdgeU(iEdge) = ( seaSurfaceTiltU(cell1) + seaSurfaceTiltU(cell2) ) * 0.5_RKIND
                seaSurfaceTiltEdgeV(iEdge) = ( seaSurfaceTiltV(cell1) + seaSurfaceTiltV(cell2) ) * 0.5_RKIND

                if(solveVelocityCGrid(iEdge) == 1) then

                   surfaceTiltForceEdgeU(iEdge) = -seaiceGravity * totalMassEdge(iEdge) * seaSurfaceTiltEdgeU(iEdge)
                   surfaceTiltForceEdgeV(iEdge) = -seaiceGravity * totalMassEdge(iEdge) * seaSurfaceTiltEdgeV(iEdge)

                else 

                   surfaceTiltForceEdgeU(iEdge) = 0.0_RKIND
                   surfaceTiltForceEdgeV(iEdge) = 0.0_RKIND 

                end if ! solveVelocityCGrid
 
             end do ! iEdge
      
          end if

       block => block % next
    end do

  end subroutine surface_tilt_ssh_gradient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  no_surface_tilt
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine no_surface_tilt(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         surfaceTiltForceEdgeU, &
         surfaceTiltForceEdgeV

    logical, pointer :: &
         config_use_c_grid

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeU", surfaceTiltForceEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeV", surfaceTiltForceEdgeV)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          ! no surface tilt
          surfaceTiltForceU = 0.0_RKIND
          surfaceTiltForceV = 0.0_RKIND

       else

          ! no surface tilt
          surfaceTiltForceEdgeU = 0.0_RKIND
          surfaceTiltForceEdgeV = 0.0_RKIND

       end if

       block => block % next
    end do

  end subroutine no_surface_tilt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_subcycle_variables
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 18th October 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_subcycle_variables(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         velocityVariationalPool, &
         velocityWeakPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11Var,       &
         strain22Var,       &
         strain12Var,       &
         strain11VarTri,    &
         strain22VarTri,    &
         strain12VarTri,    &
         stress11Var,       &
         stress22Var,       &
         stress12Var,       &
         stress11VarTri,    &
         stress22VarTri,    &
         stress12VarTri

    real(kind=RKIND), dimension(:), pointer :: &
         strain11Weak, &
         strain22Weak, &
         strain12Weak, &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:), pointer :: &
         stressDivergenceU, &
         stressDivergenceV, &
         stressDivergenceUCGrid, &
         stressDivergenceVCGrid, &
         uVelocity, &
         vVelocity, &
         uVelocityCGrid, &
         vVelocityCGrid, &
         oceanStressCoeff, &
         oceanStressCoeffEdge

    integer, dimension(:), pointer :: &
         solveStress,    &
         solveStressTri, &
         solveVelocity, &
         solveVelocityCGrid

    integer, pointer :: &
         nCells, &
         nVertices, &
         nVerticesSolve, &
         nEdgesSolve

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    logical, pointer :: &
         config_use_c_grid

    integer :: &
         iCell, &
         iVertex, &
         iEdge

    call MPAS_pool_get_config(domain % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       if ( .not. config_use_c_grid ) then

          ! divergence of stress
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

          stressDivergenceU = 0.0_RKIND
          stressDivergenceV = 0.0_RKIND

       else

          ! divergence of stress
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGrid", stressDivergenceUCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGrid", stressDivergenceVCGrid)

          stressDivergenceUCGrid = 0.0_RKIND
          stressDivergenceVCGrid = 0.0_RKIND
       
       end if

       ! sea ice velocity
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid", uVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid", vVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeffEdge", oceanStressCoeffEdge)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       if ( .not. config_use_c_grid ) then

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) /= 1) then

                uVelocity(iVertex) = 0.0_RKIND
                vVelocity(iVertex) = 0.0_RKIND

                oceanStressCoeff(iVertex) = 0.0_RKIND

             endif

          enddo ! iVertex

       else

          do iEdge = 1, nEdgesSolve

             if (solveVelocityCGrid(iEdge) /= 1) then

                uVelocityCGrid(iEdge) = 0.0_RKIND
                vVelocityCGrid(iEdge) = 0.0_RKIND

                oceanStressCoeffEdge(iEdge) = 0.0_RKIND
 
             end if

          end do

       end if

       if (trim(config_stress_divergence_scheme) == "variational") then
          ! variational

          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

          ! strains
          call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12Var)

          strain11Var = 0.0_RKIND
          strain22Var = 0.0_RKIND
          strain12Var = 0.0_RKIND

          if (config_use_c_grid) then

             ! strains
             call MPAS_pool_get_array(velocityVariationalPool, "strain11Tri", strain11VarTri)
             call MPAS_pool_get_array(velocityVariationalPool, "strain22Tri", strain22VarTri)
             call MPAS_pool_get_array(velocityVariationalPool, "strain12Tri", strain12VarTri)

             strain11VarTri = 0.0_RKIND
             strain22VarTri = 0.0_RKIND
             strain12VarTri = 0.0_RKIND

          end if

          ! stresses
          call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

          call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12Var)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (solveStress(iCell) /= 1) then

                stress11Var(:,iCell) = 0.0_RKIND
                stress22Var(:,iCell) = 0.0_RKIND
                stress12Var(:,iCell) = 0.0_RKIND

             endif

          enddo ! iCell

          if (config_use_c_grid) then

             ! stresses
             call MPAS_pool_get_array(velocitySolverPool, "solveStressTri", solveStressTri)

             call MPAS_pool_get_array(velocityVariationalPool, "stress11Tri", stress11VarTri)
             call MPAS_pool_get_array(velocityVariationalPool, "stress22Tri", stress22VarTri)
             call MPAS_pool_get_array(velocityVariationalPool, "stress12Tri", stress12VarTri)

             call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)

             do iVertex = 1, nVertices

                if (solveStressTri(iVertex) /= 1) then
  
                   stress11VarTri(:,iVertex) = 0.0_RKIND
                   stress22VarTri(:,iVertex) = 0.0_RKIND
                   stress12VarTri(:,iVertex) = 0.0_RKIND

                endif

             end do ! iVertex

          end if

       else if (trim(config_stress_divergence_scheme) == "weak") then

          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)

          ! strains
          call MPAS_pool_get_array(velocityWeakPool, "strain11", strain11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain22", strain22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain12", strain12Weak)

          strain11Weak = 0.0_RKIND
          strain22Weak = 0.0_RKIND
          strain12Weak = 0.0_RKIND

          ! stresses
          call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12Weak)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (solveStress(iCell) /= 1) then

                stress11Weak(iCell) = 0.0_RKIND
                stress22Weak(iCell) = 0.0_RKIND
                stress12Weak(iCell) = 0.0_RKIND

             endif

          enddo ! iCell

       endif ! config_stress_divergence_scheme

       block => block % next
    enddo

  end subroutine init_subcycle_variables

!-----------------------------------------------------------------------
! Sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  subcycle_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine subcycle_velocity_solver(&
       domain, &
       clock)!{{{

    use seaice_special_boundaries, only: &
         seaice_set_special_boundaries_velocity, &
         seaice_set_special_boundaries_velocity_masks

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, pointer :: &
         config_elastic_subcycle_number

    integer :: &
         iElasticSubcycle


    ! special boundaries
    call seaice_set_special_boundaries_velocity(domain)
    call seaice_set_special_boundaries_velocity_masks(domain)

    call MPAS_pool_get_config(domain % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)

    do iElasticSubcycle = 1, config_elastic_subcycle_number

       call single_subcycle_velocity_solver(&
            domain, &
            clock, &
            iElasticSubcycle)

       ! special boundaries
       call seaice_set_special_boundaries_velocity(domain)
       call seaice_set_special_boundaries_velocity_masks(domain)

    enddo


  end subroutine subcycle_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  single_subcycle_velocity_solver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine single_subcycle_velocity_solver(&
       domain, &
       clock, &
       iElasticSubcycle)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         constitutiveRelationType, &
         EVP_CONSTITUTIVE_RELATION, &
         REVISED_EVP_CONSTITUTIVE_RELATION

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    use seaice_mesh_pool, only: &
         stressDivergenceU, &
         stressDivergenceV, &
         uVelocity, &
         vVelocity

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         iElasticSubcycle !< Input: !! testing

    logical, pointer :: &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch

    integer :: &
         ierr

    ! calculate internal stresses
    call mpas_timer_start("Velocity solver internal stress")
    call seaice_internal_stress(domain)
    call mpas_timer_stop("Velocity solver internal stress")


    ! ocean stress coefficient
    call mpas_timer_start("ocn stress coef")
    call ocean_stress_coefficient(domain)
    call mpas_timer_stop("ocn stress coef")

    ! solve for velocity
    if (constitutiveRelationType == EVP_CONSTITUTIVE_RELATION) then

       call mpas_timer_start("Velocity solver compute")
       call solve_velocity(domain)
       call mpas_timer_stop("Velocity solver compute")

    else if (constitutiveRelationType == REVISED_EVP_CONSTITUTIVE_RELATION) then

       call mpas_timer_start("Velocity solver compute")
       call solve_velocity_revised(domain)
       call mpas_timer_stop("Velocity solver compute")

    endif

    ! halo exchange
    call seaice_load_balance_timers(domain, "vel before")

    call mpas_timer_start("Velocity solver halo")

    call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
    if (config_use_halo_exch) then

       call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
       if (config_aggregate_halo_exch) then

          ! aggregated halo exchange
          call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
          if (.not. config_reuse_halo_exch) then

             ! without reuse
             call mpas_dmpar_exch_group_full_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          else

             ! with reuse
             call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'velocityHaloExchangeGroup', iErr=ierr)
             if (ierr /= MPAS_DMPAR_NOERR) then
                call MPAS_log_write("failure to perform reuse halo exchange for velocityHaloExchangeGroup", MPAS_LOG_CRIT)
             endif

          endif ! config_reuse_halo_exch

       else

          ! no aggregated halo exchange
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')
          call MPAS_dmpar_field_halo_exch(domain, 'uVelocityCGrid')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocityCGrid')

       endif ! config_aggregate_halo_exch

    endif ! config_use_halo_exch

    call mpas_timer_stop("Velocity solver halo")

    call seaice_load_balance_timers(domain, "vel after")


  end subroutine single_subcycle_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_internal_stress
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 28th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_internal_stress(domain)

    use seaice_mesh_pool, only: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisGradientUNew, &
         basisGradientVNew, &
         basisIntegralsMetricNew, &
         basisGradientUTriNew, &
         basisGradientVTriNew, &
         basisIntegralsMetricTriNew, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsUNew, &
         basisIntegralsVNew, &
         basisIntegralsUTriNew, &
         basisIntegralsVTriNew, &
         cellVerticesAtVertex, &
         cellEdgesAtEdge, &
         triangleEdgesAtEdge, &
         icePressure, &
         solveStress, &
         solveStressTri, &
         solveVelocity, &
         solveVelocityCGrid, &
         stress11var => stress11, &
         stress22var => stress22, &
         stress12var => stress12, &
         stress11varTri => stress11Tri, &
         stress22varTri => stress22Tri, &
         stress12varTri => stress12Tri, &
         tanLatVertexRotatedOverRadius, &
         tanLatEdgeRotatedOverRadius, &
         uVelocity, &
         vVelocity, &
         uVelocityCGrid, &
         vVelocityCGrid

    use seaice_velocity_solver_weak, only: &
         seaice_strain_tensor_weak, &
         seaice_stress_tensor_weak, &
         seaice_stress_divergence_weak

    use seaice_velocity_solver_variational, only: &
         seaice_strain_tensor_variational, &
         seaice_strain_tensor_variational_c_grid, &
         seaice_stress_tensor_variational, &
         seaice_stress_tensor_variational_c_grid, &
         seaice_stress_divergence_variational, &
         seaice_stress_divergence_variational_c_grid, &
         seaice_average_strains_on_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type (MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         velocityWeakPool, &
         velocityVariationalPool, &
         velocityWeakVariationalPool, &
         velocitySolverPool

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    real(kind=RKIND), dimension(:), pointer :: &
         stressDivergenceU, &
         stressDivergenceV, &
         stressDivergenceUNoMetric, &
         stressDivergenceVNoMetric, &
         stressDivergenceUMetric, &
         stressDivergenceVMetric, &
         stressDivergenceUCGrid, &
         stressDivergenceVCGrid, &
         stressDivergenceUCGridNoMetric, &
         stressDivergenceVCGridNoMetric, &
         stressDivergenceUCGridMetric, &
         stressDivergenceVCGridMetric

    real(kind=RKIND), dimension(:), pointer :: &
         replacementPressureWeak, &
         strain11weak, &
         strain22weak, &
         strain12weak, &
         stress11weak, &
         stress22weak, &
         stress12weak, &
         latCellRotated, &
         latVertexRotated, &
         areaCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         normalVectorPolygon, &
         normalVectorTriangle

    real(kind=RKIND), dimension(:,:), pointer :: &
         replacementPressureVar, &
         replacementPressureVarTri, &
         strain11var, &
         strain22var, &
         strain12var, &
         strain11varTri, &
         strain22varTri, &
         strain12varTri

    real(kind=RKIND), dimension(:), pointer :: &
         variationalDenominator, &
         variationalDenominatorCGrid 

    logical, pointer :: &
         config_use_c_grid

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "boundary", boundaryPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_weak", velocityWeakPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_weak_variational", velocityWeakVariationalPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGrid", stressDivergenceUCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGrid", stressDivergenceVCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGridNoMetric", stressDivergenceUCGridNoMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGridNoMetric", stressDivergenceVCGridNoMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGridMetric", stressDivergenceUCGridMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGridMetric", stressDivergenceVCGridMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUNoMetric", stressDivergenceUNoMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVNoMetric", stressDivergenceVNoMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUMetric", stressDivergenceUMetric)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVMetric", stressDivergenceVMetric)

       if (strainSchemeType           == WEAK_STRAIN_SCHEME .or. &
           stressDivergenceSchemeType == WEAK_STRESS_DIVERGENCE_SCHEME) then

          call MPAS_pool_get_array(velocityWeakPool, "normalVectorPolygon", normalVectorPolygon)
          call MPAS_pool_get_array(velocityWeakPool, "normalVectorTriangle", normalVectorTriangle)
          call MPAS_pool_get_array(velocityWeakPool, "latCellRotated", latCellRotated)
          call MPAS_pool_get_array(velocityWeakPool, "latVertexRotated", latVertexRotated)
          call MPAS_pool_get_array(velocityWeakPool, "strain11", strain11weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain22", strain22weak)
          call MPAS_pool_get_array(velocityWeakPool, "strain12", strain12weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12weak)
          call MPAS_pool_get_array(velocityWeakPool, "replacementPressure", replacementPressureWeak)

       endif

       if (strainSchemeType           == VARIATIONAL_STRAIN_SCHEME .or. &
           stressDivergenceSchemeType == VARIATIONAL_STRESS_DIVERGENCE_SCHEME) then

          call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22var)
          call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12var)

          call MPAS_pool_get_array(velocityVariationalPool, "replacementPressure", replacementPressureVar)
          call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominator", variationalDenominator)

          call MPAS_pool_get_config(blockPtr % configs, "config_use_c_grid", config_use_c_grid)

          if (config_use_c_grid) then
             call MPAS_pool_get_array(velocityVariationalPool, "strain11Tri", strain11varTri)
             call MPAS_pool_get_array(velocityVariationalPool, "strain22Tri", strain22varTri)
             call MPAS_pool_get_array(velocityVariationalPool, "strain12Tri", strain12varTri)
             call MPAS_pool_get_array(velocityVariationalPool, "replacementPressureTri", replacementPressureVarTri)
             call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominatorCGrid", variationalDenominatorCGrid)
          end if

       endif

       ! strain
       if (strainSchemeType == WEAK_STRAIN_SCHEME) then

          call mpas_timer_start("Velocity solver strain tensor")
          call seaice_strain_tensor_weak(&
               meshPool, &
               strain11weak, &
               strain22weak, &
               strain12weak, &
               uVelocity, &
               vVelocity, &
               normalVectorPolygon, &
               latCellRotated, &
               solveStress)
          call mpas_timer_stop("Velocity solver strain tensor")

       else if (strainSchemeType == VARIATIONAL_STRAIN_SCHEME) then

          call mpas_timer_start("Velocity solver strain tensor")
          if ( .not. config_use_c_grid ) then
             call seaice_strain_tensor_variational(&
                  meshPool, &
                  strain11var, &
                  strain22var, &
                  strain12var, &
                  uVelocity, &
                  vVelocity, &
                  basisGradientU, &
                  basisGradientV, &
                  tanLatVertexRotatedOverRadius, &
                  solveStress)
          else 

             call seaice_strain_tensor_variational_c_grid(&
                  meshPool, &
                  strain11var, &
                  strain22var, &
                  strain12var, &
                  strain11varTri, &
                  strain22varTri, &
                  strain12varTri, &
                  uVelocityCGrid, &
                  vVelocityCGrid, &
                  basisGradientUNew, &
                  basisGradientVNew, &
                  basisGradientUTriNew, &
                  basisGradientVTriNew, &
                  tanLatEdgeRotatedOverRadius, &
                  solveStress, &
                  solveStressTri, &
                  solveVelocityCGrid)
          end if
          call mpas_timer_stop("Velocity solver strain tensor")

       endif

       ! average variational strains around vertex
       if (strainSchemeType == VARIATIONAL_STRAIN_SCHEME .and. &
           averageVariationalStrains) then

          call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
          call seaice_average_strains_on_vertex(&
                areaCell, &
                strain11var, &
                strain22var, &
                strain12var)
       endif

       ! weak strain / variational stress divergence
       if (strainSchemeType           == WEAK_STRAIN_SCHEME .and. &
           stressDivergenceSchemeType == VARIATIONAL_STRESS_DIVERGENCE_SCHEME) then

          call mpas_timer_start("Velocity solver interpolate strain")
          call interpolate_strains_weak_to_variational(&
               meshPool, &
               velocityWeakVariationalPool, &
               strain11weak, &
               strain22weak, &
               strain12weak, &
               strain11var, &
               strain22var, &
               strain12var)
          call mpas_timer_stop("Velocity solver interpolate strain")

       endif

       ! consitutive relation
       if (stressDivergenceSchemeType == WEAK_STRESS_DIVERGENCE_SCHEME) then

          call mpas_timer_start("Velocity solver stress tensor")
          call seaice_stress_tensor_weak(&
               meshPool, &
               stress11weak, &
               stress22weak, &
               stress12weak, &
               strain11weak, &
               strain22weak, &
               strain12weak, &
               icePressure, &
               replacementPressureWeak, &
               solveStress, &
               elasticTimeStep)
          call mpas_timer_stop("Velocity solver stress tensor")

       else if (stressDivergenceSchemeType == VARIATIONAL_STRESS_DIVERGENCE_SCHEME) then

          call mpas_timer_start("Velocity solver stress tensor")
          if ( .not. config_use_c_grid ) then
             call seaice_stress_tensor_variational(&
                  meshPool, &
                  stress11var, &
                  stress22var, &
                  stress12var, &
                  strain11var, &
                  strain22var, &
                  strain12var, &
                  icePressure, &
                  replacementPressureVar, &
                  solveStress, &
                  elasticTimeStep)
          else 
             call seaice_stress_tensor_variational_c_grid(&
                  meshPool, &
                  velocityVariationalPool, &
                  stress11var, &
                  stress22var, &
                  stress12var, &
                  stress11varTri, &
                  stress22varTri, &
                  stress12varTri, &
                  strain11var, &
                  strain22var, &
                  strain12var, &
                  strain11varTri, &
                  strain22varTri, &
                  strain12varTri, &
                  icePressure, &
                  replacementPressureVar, &
                  replacementPressureVarTri, &
                  solveStress, &
                  solveStressTri, &
                  elasticTimeStep)
          end if
          call mpas_timer_stop("Velocity solver stress tensor")

       endif

       ! stress divergence
       if (stressDivergenceSchemeType == WEAK_STRESS_DIVERGENCE_SCHEME) then

          call mpas_timer_start("Velocity solver stress divergence")
          call seaice_stress_divergence_weak(&
               meshPool, &
               stressDivergenceU, &
               stressDivergenceV, &
               stress11weak, &
               stress22weak, &
               stress12weak, &
               normalVectorTriangle, &
               latVertexRotated, &
               solveVelocity)
          call mpas_timer_stop("Velocity solver stress divergence")

       else if (stressDivergenceSchemeType == VARIATIONAL_STRESS_DIVERGENCE_SCHEME) then

          call mpas_timer_start("Velocity solver stress divergence")
          if ( .not. config_use_c_grid ) then
             call seaice_stress_divergence_variational(&
                  meshPool, &
                  stressDivergenceU, &
                  stressDivergenceV, &
                  stressDivergenceUNoMetric, &
                  stressDivergenceVNoMetric, &
                  stressDivergenceUMetric, &
                  stressDivergenceVMetric, &
                  stressDivergenceUCGrid, &
                  stressDivergenceVCGrid, &
                  stress11var, &
                  stress22var, &
                  stress12var, &
                  basisIntegralsU, &
                  basisIntegralsV, &
                  basisIntegralsMetric, &
                  variationalDenominator, &
                  tanLatVertexRotatedOverRadius, &
                  cellVerticesAtVertex, &
                  solveVelocity)
          else
             call seaice_stress_divergence_variational_c_grid(&
                  meshPool, &
                  boundaryPool, &
                  stressDivergenceU, &
                  stressDivergenceV, &
                  stressDivergenceUCGrid, &
                  stressDivergenceVCGrid, &
                  stressDivergenceUCGridNoMetric, &
                  stressDivergenceVCGridNoMetric, &
                  stressDivergenceUCGridMetric, &
                  stressDivergenceVCGridMetric, &
                  stress11var, &
                  stress22var, &
                  stress12var, &
                  stress11varTri, &
                  stress22varTri, &
                  stress12varTri, &
                  basisIntegralsUNew, &
                  basisIntegralsVNew, &
                  basisIntegralsMetricNew, &
                  basisIntegralsUTriNew, &
                  basisIntegralsVTriNew, &
                  basisIntegralsMetricTriNew, &
                  variationalDenominatorCGrid, &
                  tanLatEdgeRotatedOverRadius, &
                  cellEdgesAtEdge, &
                  triangleEdgesAtEdge, &
                  solveVelocityCGrid)
          end if 
          call mpas_timer_stop("Velocity solver stress divergence")

       endif

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_internal_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interpolate_strains_weak_to_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 29th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine interpolate_strains_weak_to_variational(&
       meshPool, &
       velocityWeakVariationalPool, &
       strain11weak, &
       strain22weak, &
       strain12weak, &
       strain11var, &
       strain22var, &
       strain12var)

    use seaice_mesh_pool, only: &
         nCells, &
         nVerticesSolve, &
         vertexDegree, &
         cellsOnVertex, &
         nEdgesOnCell, &
         verticesOnCell

    type (MPAS_pool_type), pointer :: &
         meshPool, &
         velocityWeakVariationalPool

    real(kind=RKIND), dimension(:), intent(in) :: &
         strain11weak, & !< Input/Output:
         strain22weak, & !< Input/Output:
         strain12weak    !< Input/Output:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11var, & !< Input/Output:
         strain22var, & !< Input/Output:
         strain12var    !< Input/Output:

    real(kind=RKIND), dimension(:), pointer :: &
         strain11Vertex, &
         strain22Vertex, &
         strain12Vertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    real(kind=RKIND) :: &
         denom

    integer :: &
         iVertex, &
         iCellOnVertex, &
         iCell, &
         iVertexOnCell

    call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

    call MPAS_pool_get_array(velocityWeakVariationalPool, "strain11Vertex", strain11Vertex)
    call MPAS_pool_get_array(velocityWeakVariationalPool, "strain22Vertex", strain22Vertex)
    call MPAS_pool_get_array(velocityWeakVariationalPool, "strain12Vertex", strain12Vertex)

    do iVertex = 1, nVerticesSolve

       strain11Vertex(iVertex) = 0.0_RKIND
       strain22Vertex(iVertex) = 0.0_RKIND
       strain12Vertex(iVertex) = 0.0_RKIND
       denom = 0.0_RKIND

       do iCellOnVertex = 1, vertexDegree

          iCell = cellsOnVertex(iCellOnVertex,iVertex)

          if (iCell >= 1 .and. iCell <= nCells) then
             strain11Vertex(iVertex) = strain11Vertex(iVertex) + areaCell(iCell) * strain11weak(iCell)
             strain22Vertex(iVertex) = strain22Vertex(iVertex) + areaCell(iCell) * strain22weak(iCell)
             strain12Vertex(iVertex) = strain12Vertex(iVertex) + areaCell(iCell) * strain12weak(iCell)
             denom = denom + areaCell(iCell)
          endif

       enddo ! iCellOnVertex

       strain11Vertex(iVertex) = strain11Vertex(iVertex) / denom
       strain22Vertex(iVertex) = strain22Vertex(iVertex) / denom
       strain12Vertex(iVertex) = strain12Vertex(iVertex) / denom

    enddo ! iVertex

    do iCell = 1, nCells

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell,iCell)

          strain11var(iVertexOnCell,iCell) = strain11Vertex(iVertex)
          strain22var(iVertexOnCell,iCell) = strain22Vertex(iVertex)
          strain12var(iVertexOnCell,iCell) = strain12Vertex(iVertex)

       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine interpolate_strains_weak_to_variational

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress_coefficient
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocean_stress_coefficient(domain)

    use seaice_constants, only: &
         seaiceDensitySeaWater, &
         seaiceIceOceanDragCoefficient

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressCoeff, &
         oceanStressCoeffEdge, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uOceanVelocityEdge, &
         vOceanVelocityEdge, &
         uVelocity, &
         vVelocity, &
         iceAreaVertex, &
         uVelocityCGrid, &
         vVelocityCGrid, &
         iceAreaEdge

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    logical, pointer :: &
         configUseOceanStress, &
         config_use_c_grid

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeffEdge", oceanStressCoeffEdge)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

          if (configUseOceanStress) then

             call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
             call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

             call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

             if (oceanStressType == QUADRATIC_OCEAN_STRESS) then

                do iVertex = 1, nVerticesSolve
 
                   if (solveVelocity(iVertex) == 1) then

                      oceanStressCoeff(iVertex) = seaiceIceOceanDragCoefficient * seaiceDensitySeaWater * iceAreaVertex(iVertex) * &
                           sqrt((uOceanVelocityVertex(iVertex) - uVelocity(iVertex))**2 + &
                                (vOceanVelocityVertex(iVertex) - vVelocity(iVertex))**2)

                   endif

                enddo ! iVertex

             else if (oceanStressType == LINEAR_OCEAN_STRESS) then

                do iVertex = 1, nVerticesSolve

                   if (solveVelocity(iVertex) == 1) then

                      oceanStressCoeff(iVertex) = seaiceIceOceanDragCoefficient * seaiceDensitySeaWater * iceAreaVertex(iVertex)

                   endif

                enddo ! iVertex

             end if ! oceanStressType

          else 
   
             ! no ocean stress
             oceanStressCoeff = 0.0_RKIND

          endif ! configUseOceanStress

       else

          if (configUseOceanStress) then


             call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityEdge", uOceanVelocityEdge)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityEdge", vOceanVelocityEdge)
             call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid", uVelocityCGrid)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid", vVelocityCGrid)

             call MPAS_pool_get_array(icestatePool, "iceAreaEdge", iceAreaEdge)

             if (oceanStressType == QUADRATIC_OCEAN_STRESS) then

                do iEdge = 1, nEdgesSolve

                   if (solveVelocityCGrid(iEdge) == 1) then

                      oceanStressCoeffEdge(iEdge) = seaiceIceOceanDragCoefficient * seaiceDensitySeaWater * iceAreaEdge(iEdge) * &
                           sqrt((uOceanVelocityEdge(iEdge) - uVelocityCGrid(iEdge))**2 + &
                                (vOceanVelocityEdge(iEdge) - vVelocityCGrid(iEdge))**2)

                   endif

                enddo ! iEdge

             else if (oceanStressType == LINEAR_OCEAN_STRESS) then

                do iEdge = 1, nEdgesSolve

                   if (solveVelocityCGrid(iEdge) == 1) then

                      oceanStressCoeffEdge(iEdge) = seaiceIceOceanDragCoefficient * seaiceDensitySeaWater * iceAreaEdge(iEdge)

                   endif

                enddo ! iEdge

             end if ! oceanStressType

          else 

             !no ocean stress
             oceanStressCoeffEdge = 0.0_RKIND

          end if ! configUseOceanStress
 
       end if

       block => block % next
    end do

  end subroutine ocean_stress_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  solve_velocity
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine solve_velocity(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool, &
         velocityVariationalPool, &
         meshPool

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    integer, dimension(:,:), pointer :: &
         verticesOnEdge, &
         edgesOnVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityCGrid, &
         vVelocityCGrid, &
         totalMassVertex, &
         totalMassEdge, &
         totalMassVertexfVertex, &
         totalMassEdgefEdge, &
         stressDivergenceU, &
         stressDivergenceV, &
         stressDivergenceUCgrid, &
         stressDivergenceVCGrid, &
         airStressVertexU, &
         airStressVertexV, &
         airStressEdgeU, &
         airStressEdgeV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         surfaceTiltForceEdgeU, &
         surfaceTiltForceEdgeV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressEdgeU, &
         oceanStressEdgeV, &
         oceanStressCoeff, &
         oceanStressCoeffEdge, &
         variationalDenominatorCGrid

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator, &
         sumOfDiamonds, &
         coeff

    integer, pointer :: &
         nVerticesSolve, &
         nVertices, &
         nEdgesSolve, &
         nEdges, &
         vertexDegree

    integer :: &
         iVertex, &
         iEdge, &
         vertexOnEdge1, &
         vertexOnEdge2, &
         iEdgeOnVertex

    logical, pointer :: &
         config_use_c_grid

    logical :: &
         computeVertex !to erase, here only to visualize edge-to-vertex projection

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if (.not. config_use_c_grid) then

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  totalMassVertex(iVertex) / elasticTimeStep + oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) - &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) + &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  totalMassVertex(iVertex) / elasticTimeStep  + oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + &
                  (totalMassVertex(iVertex) * uVelocity(iVertex)) / elasticTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * vVelocity(iVertex)) / elasticTimeStep

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          endif ! solveVelocity

       enddo ! iVertex

       else

          call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)
          call MPAS_pool_get_dimension(block % dimensions, "nEdges", nEdges)
          call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

          call MPAS_pool_get_array(icestatePool, "totalMassEdge", totalMassEdge)
          call MPAS_pool_get_array(meshPool, "verticesOnEdge", verticesOnEdge)
          call MPAS_pool_get_array(meshPool, "edgesOnVertex", edgesOnVertex)
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceUCGrid", stressDivergenceUCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceVCGrid", stressDivergenceVCGrid)
          call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominatorCGrid", variationalDenominatorCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid", uVelocityCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid", vVelocityCGrid)
          call MPAS_pool_get_array(velocitySolverPool, "airStressEdgeU", airStressEdgeU)
          call MPAS_pool_get_array(velocitySolverPool, "airStressEdgeV", airStressEdgeV)
          call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeU", surfaceTiltForceEdgeU)
          call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceEdgeV", surfaceTiltForceEdgeV)
          call MPAS_pool_get_array(velocitySolverPool, "totalMassEdgefEdge", totalMassEdgefEdge)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeU", oceanStressEdgeU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeV", oceanStressEdgeV)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeffEdge", oceanStressCoeffEdge)

          do iEdge = 1, nEdgesSolve

             if (solveVelocityCGrid(iEdge) == 1) then

                ! U equation
                leftMatrix(1,1) =  totalMassEdge(iEdge) / elasticTimeStep + oceanStressCoeffEdge(iEdge) * cosOceanTurningAngle
                leftMatrix(1,2) =  - totalMassEdgefEdge(iEdge) - & 
                                   oceanStressCoeffEdge(iEdge) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassEdgefEdge(iEdge))

                ! V equation
                leftMatrix(2,1) = totalMassEdgefEdge(iEdge) + &
                                  oceanStressCoeffEdge(iEdge) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassEdgefEdge(iEdge))
                leftMatrix(2,2) = totalMassEdge(iEdge) / elasticTimeStep  + oceanStressCoeffEdge(iEdge) * cosOceanTurningAngle

                ! right hand side of matrix solve
                rightHandSide(1) = stressDivergenceUCGrid(iEdge) + airStressEdgeU(iEdge) + surfaceTiltForceEdgeU(iEdge) + &
                                   oceanStressCoeffEdge(iEdge) * oceanStressEdgeU(iEdge) + & 
                                   (totalMassEdge(iEdge) * uVelocityCGrid(iEdge)) / elasticTimeStep

                rightHandSide(2) = stressDivergenceVCGrid(iEdge) + airStressEdgeV(iEdge) + surfaceTiltForceEdgeV(iEdge) + &
                                   oceanStressCoeffEdge(iEdge) * oceanStressEdgeV(iEdge) + &
                                   (totalMassEdge(iEdge) * vVelocityCGrid(iEdge)) / elasticTimeStep

                ! solve the equation
                solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

                uVelocityCGrid(iEdge) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
                vVelocityCGrid(iEdge) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

             end if

          end do

          ! -BEGIN: to be removed
          ! INTERPOLATION OF VELOCITY FROM VERTICES BACK TO EDGES 

          call MPAS_dmpar_field_halo_exch(domain, 'uVelocityCGrid')
          call MPAS_dmpar_field_halo_exch(domain, 'vVelocityCGrid')
          call MPAS_dmpar_field_halo_exch(domain, 'stressDivergenceUCGrid')
          call MPAS_dmpar_field_halo_exch(domain, 'stressDivergenceVCGrid')

          do iVertex = 1, nVerticesSolve

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND
             stressDivergenceU(iVertex) = 0.0_RKIND
             stressDivergenceV(iVertex) = 0.0_RKIND
             computeVertex = .true.
             sumOfDiamonds = 0.0_RKIND
             do iEdgeOnVertex = 1, vertexDegree
                sumOfDiamonds = sumOfDiamonds + variationalDenominatorCGrid(edgesOnVertex(iEdgeOnVertex,iVertex))              
                if (solveVelocityCGrid(edgesOnVertex(iEdgeOnVertex,iVertex)) == 0) then
                   computeVertex = .false.
                end if
             end do
             if (computeVertex) then
                do iEdgeOnVertex = 1, vertexDegree
                   iEdge = edgesOnVertex(iEdgeOnVertex,iVertex)
                   coeff = variationalDenominatorCGrid(iEdge) / sumOfDiamonds
                   uVelocity(iVertex) = uVelocity(iVertex) + coeff * uVelocityCGrid(iEdge) 
                   vVelocity(iVertex) = vVelocity(iVertex) + coeff * vVelocityCGrid(iEdge) 
                   stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + coeff * stressDivergenceUCGrid(iEdge) 
                   stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + coeff * stressDivergenceVCGrid(iEdge) 
                end do    
             end if

          end do 

          ! -END: to be removed 

       end if

       block => block % next
    end do

  end subroutine solve_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  solve_velocity_revised
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine solve_velocity_revised(domain)

    use seaice_velocity_solver_constitutive_relation, only: &
         numericalInertiaCoefficient

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         totalMassVertex, &
         totalMassVertexfVertex, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
       call MPAS_pool_get_array(velocitySolverPool, "dynamicsTimeStep", dynamicsTimeStep)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dynamicsTimeStep) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) &
                             -  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) &
                             +  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dynamicsTimeStep) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * uVelocity(iVertex) + &
                  uVelocityInitial(iVertex))) / dynamicsTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * vVelocity(iVertex) + &
                  vVelocityInitial(iVertex))) / dynamicsTimeStep

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine solve_velocity_revised

!-----------------------------------------------------------------------
! Post sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  velocity_solver_post_subcycle
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date January 13th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine velocity_solver_post_subcycle(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! calculate the final divergence and shear for ridging
    call mpas_timer_start("final div shear")
    call final_divergence_shear(domain)
    call mpas_timer_stop("final div shear")

    ! calculate principal stresses
    call mpas_timer_start("principal stress")
    call principal_stresses_driver(domain)
    call mpas_timer_stop("principal stress")

    ! calculate final stress after subcycling
    call mpas_timer_start("ocn stress final")
    call ocean_stress_final(domain)
    call mpas_timer_stop("ocn stress final")

  end subroutine velocity_solver_post_subcycle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_divergence_shear
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date July 9th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine final_divergence_shear(domain)

    use seaice_velocity_solver_variational, only: &
         seaice_final_divergence_shear_variational

    use seaice_velocity_solver_weak, only: &
         seaice_final_divergence_shear_weak

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

       if (trim(config_stress_divergence_scheme) == "weak") then

          call seaice_final_divergence_shear_weak(block)

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call seaice_final_divergence_shear_variational(block)

       endif

       block => block % next
    enddo

  end subroutine final_divergence_shear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  principal_stresses_driver
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine principal_stresses_driver(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocityWeakPool, &
         velocityVariationalPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         replacementPressureWeak, &
         principalStress1Weak, &
         principalStress2Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         replacementPressureVar, &
         principalStress1Var, &
         principalStress2Var

    real(kind=RKIND), dimension(:), pointer :: &
         stress11Weak, &
         stress22Weak, &
         stress12Weak

    real(kind=RKIND), dimension(:,:), pointer :: &
         stress11Var, &
         stress22Var, &
         stress12Var

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell, &
         iVertexOnCell

    character(len=strKIND), pointer :: &
         config_stress_divergence_scheme

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       ! calculate the principal stresses
       if (trim(config_stress_divergence_scheme) == "weak") then

          call MPAS_pool_get_subpool(block % structs, "velocity_weak", velocityWeakPool)

          call MPAS_pool_get_array(velocityWeakPool, "stress11", stress11Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress22", stress22Weak)
          call MPAS_pool_get_array(velocityWeakPool, "stress12", stress12Weak)
          call MPAS_pool_get_array(velocityWeakPool, "replacementPressure", replacementPressureWeak)
          call MPAS_pool_get_array(velocityWeakPool, "principalStress1", principalStress1Weak)
          call MPAS_pool_get_array(velocityWeakPool, "principalStress2", principalStress2Weak)

          do iCell = 1, nCellsSolve

             call principal_stresses(&
                  principalStress1Weak(iCell), &
                  principalStress2Weak(iCell), &
                  stress11Weak(iCell), &
                  stress22Weak(iCell), &
                  stress12Weak(iCell), &
                  replacementPressureWeak(iCell))

          enddo ! iCell

       else if (trim(config_stress_divergence_scheme) == "variational") then

          call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

          call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)

          call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22Var)
          call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12Var)
          call MPAS_pool_get_array(velocityVariationalPool, "replacementPressure", replacementPressureVar)
          call MPAS_pool_get_array(velocityVariationalPool, "principalStress1", principalStress1Var)
          call MPAS_pool_get_array(velocityVariationalPool, "principalStress2", principalStress2Var)

          do iCell = 1, nCellsSolve
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call principal_stresses(&
                     principalStress1Var(iVertexOnCell,iCell), &
                     principalStress2Var(iVertexOnCell,iCell), &
                     stress11Var(iVertexOnCell,iCell), &
                     stress22Var(iVertexOnCell,iCell), &
                     stress12Var(iVertexOnCell,iCell), &
                     replacementPressureVar(iVertexOnCell,iCell))

             enddo ! iVertexOnCell
          enddo ! iCell

       endif

       block => block % next
    enddo

  end subroutine principal_stresses_driver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  principal_stresses
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine principal_stresses(&
       principalStress1, &
       principalStress2, &
       stress11, &
       stress22, &
       stress12, &
       replacementPressure)

    use seaice_constants, only: &
         seaicePuny

    real(kind=RKIND), intent(out) :: &
         principalStress1, & !< Output:
         principalStress2    !< Output:

    real(kind=RKIND), intent(in) :: &
         stress11, &         !< Input:
         stress22, &         !< Input:
         stress12, &         !< Input:
         replacementPressure !< Input:

    real(kind=RKIND) :: &
         sqrtContents

    if (replacementPressure > seaicePuny) then

       sqrtContents = (stress11 + stress22)**2 - &
                      4.0_RKIND * stress11 * stress22 + &
                      4.0_RKIND * stress12**2

       principalStress1 = 0.5_RKIND * (stress11 + stress22) + 0.5_RKIND * sqrt(sqrtContents)
       principalStress2 = 0.5_RKIND * (stress11 + stress22) - 0.5_RKIND * sqrt(sqrtContents)

       principalStress1 = principalStress1 / replacementPressure
       principalStress2 = principalStress2 / replacementPressure

    else

       principalStress1 = 1.0e30_RKIND
       principalStress2 = 1.0e30_RKIND

    endif

  end subroutine principal_stresses

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocean_stress_final
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine ocean_stress_final(domain)

    use seaice_mesh, only: &
         seaice_interpolate_vertex_to_cell

    use seaice_diagnostics, only: &
         seaice_load_balance_timers

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         icestatePool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         iceAreaEdge, &
         fVertex, &
         fEdge, &
         uVelocity, &
         vVelocity, &
         uVelocityCGrid, &
         vVelocityCGrid, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uOceanVelocityEdge, &
         vOceanVelocityEdge, &
         oceanStressU, &
         oceanStressV, &
         oceanStressEdgeU, &
         oceanStressEdgeV, &
         oceanStressCellU, &
         oceanStressCellV, &
         oceanStressCoeff, &
         oceanStressCoeffEdge

    integer, dimension(:), pointer :: &
         solveVelocity, &
         solveVelocityCGrid

    logical, pointer :: &
         configUseOceanStress, &
         config_use_halo_exch, &
         config_aggregate_halo_exch, &
         config_reuse_halo_exch, &
         config_use_c_grid

    integer, pointer :: &
         nVerticesSolve, &
         nEdgesSolve

    integer :: &
         iVertex, &
         iEdge, &
         ierr

    ! get ocean stress coefficient
    call ocean_stress_coefficient(domain)

    ! use stress config
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_ocean_stress", configUseOceanStress)

    ! get ocean stress on vertices or edges
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeU", oceanStressEdgeU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeV", oceanStressEdgeV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

       call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

       if ( .not. config_use_c_grid) then

          if (configUseOceanStress) then

             call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolve)

             call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

             call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

             call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
             call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
             call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

             do iVertex = 1, nVerticesSolve

                if (solveVelocity(iVertex) == 1) then

                   oceanStressU(iVertex) = oceanStressCoeff(iVertex) * &
                        ((uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * cosOceanTurningAngle &
                        - (vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * sinOceanTurningAngle * &
                        sign(1.0_RKIND,fVertex(iVertex)))
                   oceanStressV(iVertex) = oceanStressCoeff(iVertex) * &
                        ((vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * cosOceanTurningAngle &
                        + (uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * sinOceanTurningAngle * &
                        sign(1.0_RKIND,fVertex(iVertex)))

                   oceanStressU(iVertex) = oceanStressU(iVertex) / iceAreaVertex(iVertex)
                   oceanStressV(iVertex) = oceanStressV(iVertex) / iceAreaVertex(iVertex)

                else

                   oceanStressU(iVertex) = 0.0_RKIND
                   oceanStressV(iVertex) = 0.0_RKIND

                endif ! solveVelocity

             enddo ! iVertex

          else

             ! no ocean stress
             oceanStressU = 0.0_RKIND
             oceanStressV = 0.0_RKIND

          endif

       else

          if (configUseOceanStress) then

             call MPAS_pool_get_dimension(meshPool, "nEdgesSolve", nEdgesSolve)

             call MPAS_pool_get_array(meshPool, "fEdge", fEdge)

             call MPAS_pool_get_array(icestatePool, "iceAreaEdge", iceAreaEdge)

             call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid", uVelocityCGrid)
             call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid", vVelocityCGrid)
             call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityEdge", uOceanVelocityEdge)
             call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityEdge", vOceanVelocityEdge)
             call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeffEdge", oceanStressCoeffEdge)
             call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)

             do iEdge = 1, nEdgesSolve

                if (solveVelocityCGrid(iEdge) == 1) then

                   oceanStressEdgeU(iEdge) = oceanStressCoeffEdge(iEdge) * &
                        ((uOceanVelocityEdge(iEdge) - uVelocityCGrid(iEdge)) * cosOceanTurningAngle &
                        - (vOceanVelocityEdge(iEdge) - vVelocityCGrid(iEdge)) * sinOceanTurningAngle * &
                        sign(1.0_RKIND,fEdge(iEdge)))
                   oceanStressEdgeV(iEdge) = oceanStressCoeffEdge(iEdge) * &
                        ((vOceanVelocityEdge(iEdge) - vVelocityCGrid(iEdge)) * cosOceanTurningAngle &
                        + (uOceanVelocityEdge(iEdge) - uVelocityCGrid(iEdge)) * sinOceanTurningAngle * &
                        sign(1.0_RKIND,fEdge(iEdge)))

                   oceanStressEdgeU(iEdge) = oceanStressEdgeU(iEdge) / iceAreaEdge(iEdge)
                   oceanStressEdgeV(iEdge) = oceanStressEdgeV(iEdge) / iceAreaEdge(iEdge)

                else

                   oceanStressEdgeU(iEdge) = 0.0_RKIND
                   oceanStressEdgeV(iEdge) = 0.0_RKIND

                endif ! solveVelocityCGrid

             enddo ! iEdge

          else

             ! no ocean stress
             oceanStressEdgeU = 0.0_RKIND
             oceanStressEdgeV = 0.0_RKIND

          endif

       end if
 
       block => block % next
    end do

    ! get ocean stress on cells
    if (configUseOceanStress) then

       ! halo exchange ocean stress
       call seaice_load_balance_timers(domain, "vel prep before")

       call mpas_timer_start("ocean stress halo")

       call MPAS_pool_get_config(domain % configs, "config_use_halo_exch", config_use_halo_exch)
       if (config_use_halo_exch) then

          call MPAS_pool_get_config(domain % configs, "config_aggregate_halo_exch", config_aggregate_halo_exch)
          if (config_aggregate_halo_exch) then

             ! aggregated halo exchange
             call MPAS_pool_get_config(domain % configs, "config_reuse_halo_exch", config_reuse_halo_exch)
             if (.not. config_reuse_halo_exch) then

                ! without reuse
                call mpas_dmpar_exch_group_full_halo_exch(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
                if (ierr /= MPAS_DMPAR_NOERR) then
                   call MPAS_log_write("failure to perform halo exchange for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
                endif

             else

                ! with reuse
                call mpas_dmpar_exch_group_reuse_halo_exch(domain, 'oceanStressHaloExchangeGroup', iErr=ierr)
                if (ierr /= MPAS_DMPAR_NOERR) then
                   call MPAS_log_write("failure to perform reuse halo exchange for oceanStressHaloExchangeGroup", MPAS_LOG_CRIT)
                endif

             endif ! config_reuse_halo_exch

          else

             ! no aggregated halo exchange
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressU')
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressV')
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressEdgeU')
             call MPAS_dmpar_field_halo_exch(domain, 'oceanStressEdgeV')

          endif ! config_aggregate_halo_exch

       endif ! config_use_halo_exch

       call mpas_timer_stop("ocean stress halo")

       call seaice_load_balance_timers(domain, "vel prep after")

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
          call MPAS_pool_get_dimension(block % dimensions, "nEdgesSolve", nEdgesSolve)

          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)
          call MPAS_pool_get_array(icestatePool, "iceAreaEdge", iceAreaEdge)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeU", oceanStressEdgeU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressEdgeV", oceanStressEdgeV)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid", solveVelocityCGrid)

          call MPAS_pool_get_config(block % configs, "config_use_c_grid", config_use_c_grid)

          if (.not. config_use_c_grid) then

             ! interpolate ocean stress to cells
             call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellU, oceanStressU)
             call seaice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellV, oceanStressV)

             ! multiply ocean stress back by ice area
             do iVertex = 1, nVerticesSolve

                if (solveVelocity(iVertex) == 1) then

                   oceanStressU(iVertex) = oceanStressU(iVertex) * iceAreaVertex(iVertex)
                   oceanStressV(iVertex) = oceanStressV(iVertex) * iceAreaVertex(iVertex)

                endif ! solveVelocity

             enddo ! iVertex

          else

             ! interpolate ocean stress to cell from edges
             ! TODO: I don't think it is necessary right now              

             do iEdge = 1, nEdgesSolve

                if (solveVelocityCGrid(iEdge) == 1) then

                   oceanStressEdgeU(iEdge) = oceanStressEdgeU(iEdge) * iceAreaEdge(iEdge)
                   oceanStressEdgeV(iEdge) = oceanStressEdgeV(iEdge) * iceAreaEdge(iEdge)

                end if

             end do

          end if

          block => block % next
       end do

    else

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          oceanStressCellU = 0.0_RKIND
          oceanStressCellV = 0.0_RKIND

          block => block % next
       end do

    endif ! configUseOceanStress

  end subroutine ocean_stress_final

  !-----------------------------------------------------------------------

end module seaice_velocity_solver
