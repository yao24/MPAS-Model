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
