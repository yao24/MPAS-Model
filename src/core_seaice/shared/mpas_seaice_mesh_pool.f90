










!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool
!
!> \brief
!> \date 2020
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_mesh_pool





   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   implicit none
   private

   integer, public :: &
      nCells,         &
      nEdges,         &
      nEdgesSolve,    &
      nVertices,      &
      nVerticesSolve, &
      vertexDegree,   &
      TWO

   integer, public, dimension(:), pointer :: & 
      nEdgesOnCell,       &
      solveStress,        &
      solveStressTri,     &
      solveVelocity,      &
      solveVelocityCGrid

   integer, public, dimension(:,:), pointer :: &
      verticesOnCell, &
      cellsOnVertex,  &
      edgesOnCell,    &
      edgesOnVertex,  &
      verticesOnEdge, &
      cellsOnEdge, &
      cellVerticesAtVertex, &
      cellEdgesAtEdge, &
      triangleEdgesAtEdge

   real(kind=RKIND), public, dimension(:), pointer :: &
      areaTriangle,   &
      tanLatVertexRotatedOverRadius, &
      tanLatEdgeRotatedOverRadius, &
      icePressure,    &
      uVelocity,      &
      vVelocity,      &
      uVelocityCGrid,      &
      vVelocityCGrid,      &
      stressDivergenceU, &
      stressDivergenceV

   real(kind=RKIND), public, dimension(:,:), pointer :: &
      stress11,         &
      stress12,         &
      stress22,         &
      stress11Tri,      &
      stress12Tri,      &
      stress22Tri,      &
      kiteAreasOnVertex   

   real(kind=RKIND), public, dimension(:,:,:), pointer :: &
      basisGradientU, &
      basisGradientV, &
      basisIntegralsU,&
      basisIntegralsV,&
      basisIntegralsMetric

   real(kind=RKIND), public, dimension(:,:,:,:), pointer :: &
      basisGradientUNew, &
      basisGradientVNew, &
      basisIntegralsUNew,&
      basisIntegralsVNew,&
      basisIntegralsMetricNew, &
      basisGradientUTriNew, &
      basisGradientVTriNew, &
      basisIntegralsUTriNew,&
      basisIntegralsVTriNew,&
      basisIntegralsMetricTriNew

   public ::                   &
      seaice_mesh_pool_create, &
      seaice_mesh_pool_update, &
      seaice_mesh_pool_destroy

!-----------------------------------------------------------------------

contains



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_create
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_create(&
       domain)!{{{

    type(domain_type) :: &
         domain

    integer :: &
         blockCount

    type(block_type), pointer :: &
         block

    type (mpas_pool_type), pointer :: &
         meshPool,               &
         velocitySolverPool,     &
         velocityVariationalPool

    integer, pointer ::   &
         nCellsTmp,         &
         nEdgesTmp,         &
         nEdgesSolveTmp,    &
         nVerticesTmp,      &
         nVerticesSolveTmp, &
         vertexDegreeTmp,   &
         TWOTmp

    blockCount = 0
    block => domain % blocklist
    do while ( associated(block) )

       blockCount = blockCount + 1
       if (blockCount > 1) then
          call mpas_log_write('seaice_mesh_pool_create: more than one block is no longer supported', MPAS_LOG_CRIT)
       endif

       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'velocity_solver', velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

       ! convert mesh dimensions from pointers to scalars
       call mpas_pool_get_dimension(meshPool, 'nCells', nCellsTmp)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdgesTmp)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolveTmp)
       call MPAS_pool_get_dimension(meshPool, "nVertices", nVerticesTmp)
       call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolveTmp)
       call MPAS_pool_get_dimension(meshPool, "vertexDegree", vertexDegreeTmp)
       call MPAS_pool_get_dimension(meshPool, "TWO", TWOTmp)

       nCells         = nCellsTmp
       nEdges         = nEdgesTmp
       nEdgesSolve    = nEdgesSolveTmp
       nVertices      = nVerticesTmp
       nVerticesSolve = nVerticesSolveTmp
       vertexDegree   = vertexDegreeTmp
       TWO = TWOTmp

       ! point to existing arrays
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell',   nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'cellsOnVertex',  cellsOnVertex)
       call mpas_pool_get_array(meshPool, 'edgesOnCell',  edgesOnCell)
       call mpas_pool_get_array(meshPool, 'edgesOnVertex',  edgesOnVertex)
       call mpas_pool_get_array(meshPool, 'verticesOnEdge',  verticesOnEdge)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge',  cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'areaTriangle',   areaTriangle)
       call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex',   kiteAreasOnVertex)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress",       solveStress)
       call MPAS_pool_get_array(velocitySolverPool, "solveStressTri",       solveStressTri)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity",     solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocityCGrid",     solveVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "icePressure",       icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity",         uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity",         vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityCGrid",         uVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityCGrid",         vVelocityCGrid)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientU",       basisGradientU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientV",       basisGradientV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsU",      basisIntegralsU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsV",      basisIntegralsV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetric", basisIntegralsMetric)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientUNew",       basisGradientUNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientVNew",       basisGradientVNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsUNew",      basisIntegralsUNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsVNew",      basisIntegralsVNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetricNew", basisIntegralsMetricNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientUTriNew",       basisGradientUTriNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientVTriNew",       basisGradientVTriNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsUTriNew",      basisIntegralsUTriNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsVTriNew",      basisIntegralsVTriNew)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetricTriNew", basisIntegralsMetricTriNew)
       call MPAS_pool_get_array(velocityVariationalPool, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
       call MPAS_pool_get_array(velocityVariationalPool, "tanLatEdgeRotatedOverRadius", tanLatEdgeRotatedOverRadius)
       call MPAS_pool_get_array(velocityVariationalPool, "cellVerticesAtVertex", cellVerticesAtVertex)
       call MPAS_pool_get_array(velocityVariationalPool, "cellEdgesAtEdge", cellEdgesAtEdge)
       call MPAS_pool_get_array(velocityVariationalPool, "triangleEdgesAtEdge", triangleEdgesAtEdge)
       call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11)
       call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22)
       call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12)
       call MPAS_pool_get_array(velocityVariationalPool, "stress11Tri", stress11Tri)
       call MPAS_pool_get_array(velocityVariationalPool, "stress22Tri", stress22Tri)
       call MPAS_pool_get_array(velocityVariationalPool, "stress12Tri", stress12Tri)


       block => block % next
    end do

  end subroutine seaice_mesh_pool_create!}}}
!-----------------------------------------------------------------------



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_destroy
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_destroy(&
       err)!{{{

    integer, intent(out) :: &
      err                   ! returned error flag

    err = 0


    ! then nullify on host
    nullify(nEdgesOnCell,         &
         verticesOnCell,          &
         cellsOnVertex,           &
         edgesOnCell,             &
         edgesOnVertex,           &
         verticesOnEdge,          &
         cellsOnEdge,             &
         areaTriangle,            &
         solveStress,             &
         solveStressTri,          &
         solveVelocity,           &
         solveVelocityCGrid,      &
         icePressure,             &
         uVelocity,               &
         vVelocity,               &
         uVelocityCGrid,          &
         vVelocityCGrid,          &
         stressDivergenceU,       &
         stressDivergenceV,       &
         basisGradientU,          &
         basisGradientV,          &
         basisIntegralsU,         &
         basisIntegralsV,         &
         basisIntegralsMetric,    &
         basisGradientUNew,       &
         basisGradientVNew,       &
         basisIntegralsUNew,      &
         basisIntegralsVNew,      &
         basisIntegralsMetricNew, &
         basisGradientUTriNew,       &
         basisGradientVTriNew,       &
         basisIntegralsUTriNew,      &
         basisIntegralsVTriNew,      &
         basisIntegralsMetricTriNew, &
         tanLatVertexRotatedOverRadius, &
         tanLatEdgeRotatedOverRadius,   &
         cellVerticesAtVertex, &
         cellEdgesAtEdge,      &
         triangleEdgesAtEdge,  &
         stress11,             &
         stress12,             &
         stress22,             &
         stress11Tri,          &
         stress12Tri,          &
         stress22Tri,          &
         kiteAreasOnVertex     &   
    )

  end subroutine seaice_mesh_pool_destroy!}}}
!-----------------------------------------------------------------------



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_update
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_update(&
       domain)!{{{

    type(domain_type) :: &
         domain


  end subroutine seaice_mesh_pool_update!}}}
!-----------------------------------------------------------------------



end module seaice_mesh_pool
