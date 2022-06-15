










!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_variational





  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver_variational, &
       seaice_strain_tensor_variational, &
       seaice_strain_tensor_variational_c_grid, &
       seaice_average_strains_on_vertex, &
       seaice_stress_tensor_variational, &
       seaice_stress_tensor_variational_c_grid, &
       seaice_stress_divergence_variational, &
       seaice_stress_divergence_variational_c_grid, &
       seaice_final_divergence_shear_variational

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver_variational(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms, &
       variationalBasisType, &
       variationalDenominatorType, &
       integrationType, &
       integrationOrder, &
       useCGrid,useProjectedBGrid)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output:
         boundary                !< Input/Output:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms, &  !< Input:
         useCGrid, &               !< Input:
         useProjectedBGrid

    character(len=*), intent(in) :: &
         variationalBasisType, &       !< Input:
         variationalDenominatorType, & !< Input:
         integrationType               !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

    call init_velocity_solver_variational_primary_mesh(&
         mesh, &
         velocity_variational, &
         boundary, &
         rotateCartesianGrid, &
         includeMetricTerms, &
         variationalBasisType, &
         variationalDenominatorType, &
         integrationType, &
         integrationOrder, &
         useCGrid,useProjectedBGrid)

  end subroutine seaice_init_velocity_solver_variational

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_velocity_solver_variational_primary_mesh
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_velocity_solver_variational_primary_mesh(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms, &
       variationalBasisType, &
       variationalDenominatorType, &
       integrationType, &
       integrationOrder, &
       useCGrid, &
       useProjectedBGrid)!{{{

    use seaice_mesh, only: &
         seaice_cell_vertices_at_vertex, &
         seaice_cell_edges_at_edge, &
         seaice_triangle_edges_at_edge

    use seaice_velocity_solver_variational_shared, only: &
         seaice_calc_local_coords, &
         seaice_calc_local_coords_projected_b_grid, &
         seaice_calc_local_coords_c_grid, &
         seaice_calc_variational_metric_terms

    use seaice_velocity_solver_wachspress, only: &
         seaice_init_velocity_solver_wachspress

    use seaice_velocity_solver_pwl, only: &
         seaice_init_velocity_solver_pwl, &
         seaice_init_velocity_solver_pwl_c_grid

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output:
         boundary                !< Input/Output:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms, &  !< Input:
         useCGrid, &               !< Input:
         useProjectedBGrid

    character(len=*), intent(in) :: &
         variationalBasisType, &       !< Input:
         variationalDenominatorType, & !< Input:
         integrationType               !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex, &
         cellEdgesAtEdge, &
         triangleEdgesAtEdge

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius, &
         tanLatEdgeRotatedOverRadius, &
         variationalDenominator, &
         variationalDenominatorCGrid, &
         areaCellCGrid, &
         areaTriCGrid

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV

    real(kind=RKIND), dimension(:,:,:,:), pointer :: &
         basisGradientUNew, &
         basisGradientVNew, &
         basisIntegralsMetricNew, &
         basisIntegralsUNew, &
         basisIntegralsVNew, &
         basisGradientUTriNew, &
         basisGradientVTriNew, &
         basisIntegralsMetricTriNew, &
         basisIntegralsUTriNew, &
         basisIntegralsVTriNew

    integer, pointer :: &
         nCells, &
         nVertices, &
         nEdges, &
         nEdgesSolve, &
         vertexDegree, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell, &
         interiorVertex

    integer, dimension(:,:), pointer :: &
         verticesOnCell, &
         cellsOnVertex, &
         edgesOnCell, &
         cellsOnEdge, &
         edgesOnVertex, &
         verticesOnEdge

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, &
         yVertex, &
         zVertex, &
         xEdge, &
         yEdge, &
         zEdge, &
         xCell, &
         yCell, &
         zCell, &
         latCell, &
         lonCell, &
         latVertex, &
         lonVertex, &
         latEdge, &
         lonEdge, &
         areaCell, &
         areaTriangle, &
         dvEdge

    logical, pointer :: &
         on_a_sphere

    real(kind=RKIND), pointer :: &
         sphere_radius

    integer, dimension(:), allocatable :: & !needed to use the existing refactored functions
         nEdgesOnVertex

    real(kind=RKIND), dimension(:,:), allocatable :: &
         xLocal, &
         yLocal

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         xLocalNew, &
         yLocalNew, &
         xLocalNewNew, &
         yLocalNewNew, &
         xLocalTriNew, &
         yLocalTriNew

    integer :: iEdge, iCell, iVertex, iCellOnEdge, iVertexOnEdge

    call MPAS_pool_get_config(mesh, "on_a_sphere", on_a_sphere)
    call MPAS_pool_get_config(mesh, "sphere_radius", sphere_radius)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "nEdges", nEdges)
    call MPAS_pool_get_dimension(mesh, "nEdgesSolve", nEdgesSolve)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnEdge", verticesOnEdge)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnVertex", edgesOnVertex)
    call MPAS_pool_get_array(mesh, "cellsOnEdge", cellsOnEdge)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)
    call MPAS_pool_get_array(mesh, "yVertex", yVertex)
    call MPAS_pool_get_array(mesh, "zVertex", zVertex)
    call MPAS_pool_get_array(mesh, "xEdge", xEdge)
    call MPAS_pool_get_array(mesh, "yEdge", yEdge)
    call MPAS_pool_get_array(mesh, "zEdge", zEdge)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)
    call MPAS_pool_get_array(mesh, "latCell", latCell)
    call MPAS_pool_get_array(mesh, "lonCell", lonCell)
    call MPAS_pool_get_array(mesh, "latVertex", latVertex)
    call MPAS_pool_get_array(mesh, "lonVertex", lonVertex)
    call MPAS_pool_get_array(mesh, "latEdge", latEdge)
    call MPAS_pool_get_array(mesh, "lonEdge", lonEdge)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)
    call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(velocity_variational, "cellVerticesAtVertex", cellVerticesAtVertex)
    call MPAS_pool_get_array(velocity_variational, "cellEdgesAtEdge", cellEdgesAtEdge)
    call MPAS_pool_get_array(velocity_variational, "triangleEdgesAtEdge", triangleEdgesAtEdge)
    call MPAS_pool_get_array(velocity_variational, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "tanLatEdgeRotatedOverRadius", tanLatEdgeRotatedOverRadius)
    call MPAS_pool_get_array(velocity_variational, "basisGradientU", basisGradientU)
    call MPAS_pool_get_array(velocity_variational, "basisGradientV", basisGradientV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsU", basisIntegralsU)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsV", basisIntegralsV)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetric", basisIntegralsMetric)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsUNew", basisIntegralsUNew)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsVNew", basisIntegralsVNew)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetricNew", basisIntegralsMetricNew)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsUTriNew", basisIntegralsUTriNew)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsVTriNew", basisIntegralsVTriNew)
    call MPAS_pool_get_array(velocity_variational, "basisIntegralsMetricTriNew", basisIntegralsMetricTriNew)
    call MPAS_pool_get_array(velocity_variational, "basisGradientUNew", basisGradientUNew)
    call MPAS_pool_get_array(velocity_variational, "basisGradientVNew", basisGradientVNew)
    call MPAS_pool_get_array(velocity_variational, "basisGradientUTriNew", basisGradientUTriNew)
    call MPAS_pool_get_array(velocity_variational, "basisGradientVTriNew", basisGradientVTriNew)
    call MPAS_pool_get_array(velocity_variational, "variationalDenominator", variationalDenominator)
    call MPAS_pool_get_array(velocity_variational, "variationalDenominatorCGrid", variationalDenominatorCGrid)
    call MPAS_pool_get_array(velocity_variational, "areaCellCGrid", areaCellCGrid)
    call MPAS_pool_get_array(velocity_variational, "areaTriCGrid", areaTriCGrid)

    call mpas_timer_start("variational calc_metric_terms")
    if ( .not. useCGrid ) then
       call seaice_calc_variational_metric_terms(&
            tanLatVertexRotatedOverRadius, &
            nVertices, &
            xVertex, &
            yVertex, &
            zVertex, &
            sphere_radius, &
            rotateCartesianGrid, &
            includeMetricTerms)
    else
           call seaice_calc_variational_metric_terms(&
            tanLatEdgeRotatedOverRadius, &
            nEdges, &
            xEdge, &
            yEdge, &
            zEdge, &
            sphere_radius, &
            rotateCartesianGrid, &
            includeMetricTerms)
    end if
    call mpas_timer_stop("variational calc_metric_terms")

    call mpas_timer_start("variational vertices_at_vertex")
    if ( .not. useCGrid ) then
       call seaice_cell_vertices_at_vertex(&
            cellVerticesAtVertex, &
            nVertices, &
            vertexDegree, &
            nEdgesOnCell, &
            verticesOnCell, &
            cellsOnVertex)
    else
       call seaice_cell_edges_at_edge(&
            cellEdgesAtEdge, &
            nEdges, &
            nEdgesOnCell, &
            edgesOnCell, &
            cellsOnEdge)
       call seaice_triangle_edges_at_edge(&
            triangleEdgesAtEdge, &
            nEdges, &
            vertexDegree, &
            edgesOnVertex, &
            verticesOnEdge)
    end if
    call mpas_timer_stop("variational vertices_at_vertex")

    call mpas_timer_start("variational calc_local_coords")

    if (.not. useCGrid) then
       if (.not. useProjectedBGrid) then 
            allocate(xLocal(maxEdges,nCells))
            allocate(yLocal(maxEdges,nCells))
        else
            allocate(xLocalNewNew(nVertices,maxEdges,3))
            allocate(yLocalNewNew(nVertices,maxEdges,3))
        endif
       
    else 

       allocate(xLocalNew(nEdges,maxEdges,2))
       allocate(yLocalNew(nEdges,maxEdges,2))

       allocate(xLocalTriNew(nEdges,vertexDegree,2))
       allocate(yLocalTriNew(nEdges,vertexDegree,2))

       allocate(nEdgesOnVertex(nVertices))
       nEdgesOnVertex(:) = vertexDegree

    end if

    if ( .not. useCGrid ) then
       if (.not. useProjectedBGrid) then 
            call seaice_calc_local_coords(&
                 xLocal, &
                 yLocal, &
                 nCells, &
                 nEdgesOnCell, &
                 verticesOnCell, &
                 xVertex, &
                 yVertex, &
                 zVertex, &
                 xCell, &
                 yCell, &
                 zCell, &
                 rotateCartesianGrid, &
                 on_a_sphere)
        else
            call seaice_calc_local_coords_projected_b_grid(&
                xLocalNewNew, &
                yLocalNewNew, &
                nCells, &
                nEdgesOnCell, &
                verticesOnCell, &
                nVertices, & 
                cellsOnVertex, & 
                vertexDegree, &
                xVertex, &
                yVertex, &
                zVertex, &
                xCell, &
                yCell, &
                zCell, &
                rotateCartesianGrid, &
                on_a_sphere)
        endif
    else

          call seaice_calc_local_coords_c_grid(&
             xLocalNew, &
             yLocalNew, &
             xLocalTriNew, &
             yLocalTriNew, &
             nEdges, &
             cellsOnEdge, &
             verticesOnEdge, &
             edgesOnCell, &
             edgesOnVertex, &
             nEdgesOnCell, &
             vertexDegree, &
             xCell, &
             yCell, &
             zCell, &
             xVertex, &
             yVertex, &
             zVertex, &
             xEdge, &
             yEdge, &
             zEdge, &
             rotateCartesianGrid, &
             on_a_sphere, &
             dvEdge, &
             interiorVertex)

    end if

    call mpas_timer_stop("variational calc_local_coords")

    if (useCGrid) then
       call mpas_timer_start("compute_c_grid_polygon_and_tri_areas")
       call compute_c_grid_polygon_and_tri_areas(&
           on_a_sphere, &
           sphere_radius, &
           edgesOnCell, &
           edgesOnVertex, &
           nEdges, &
           nCells, &
           nVertices, &
           nEdgesOnCell, &
           vertexDegree, &
           xCell, &
           yCell, &
           zCell, &
           xVertex, &
           yVertex, &
           zVertex, &
           xEdge, &
           yEdge, &
           zEdge, &
           latCell, &
           lonCell, &
           latVertex, &
           lonVertex, &
           latEdge, &
           lonEdge, &
           areaCellCGrid, &
           areaTriCGrid)
       call mpas_timer_stop("compute_c_grid_polygon_and_tri_areas")
    end if

    if (trim(variationalBasisType) == "wachspress") then

       if( .not. useCGrid) then

          do iCell = 1, nCells
             call seaice_init_velocity_solver_wachspress(&
                iCell, &
                maxEdges, &
                nEdgesOnCell, &
                xLocal(:,iCell), &
                yLocal(:,iCell), &
                rotateCartesianGrid, &
                includeMetricTerms, &
                on_a_sphere, &
                integrationType, &
                integrationOrder, &
                sphere_radius, &
                basisGradientU(:,:,iCell), &
                basisGradientV(:,:,iCell), &
                basisIntegralsU(:,:,iCell), &
                basisIntegralsV(:,:,iCell), &
                basisIntegralsMetric(:,:,iCell))
          end do

       else

          do iEdge = 1, nEdges

             do iCellOnEdge = 1, 2

                iCell = cellsOnEdge(iCellOnEdge,iEdge)

                if (iCell <= nCells) then

                   call seaice_init_velocity_solver_wachspress(&
                      iCell, &
                      maxEdges, &
                      nEdgesOnCell, &
                      xLocalNew(iEdge,:,iCellOnEdge), &
                      yLocalNew(iEdge,:,iCellOnEdge), &
                      rotateCartesianGrid, &
                      includeMetricTerms, &
                      on_a_sphere, &
                      integrationType, &
                      integrationOrder, &
                      sphere_radius, &
                      basisGradientUNew(iEdge,:,:,iCellOnEdge), &
                      basisGradientVNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsUNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsVNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsMetricNew(iEdge,:,:,iCellOnEdge))


                 end if

              end do

              do iVertexOnEdge = 1, 2

                 iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                 call seaice_init_velocity_solver_wachspress(&
                    iVertex, &
                    vertexDegree, &
                    nEdgesOnVertex, &
                    xLocalTriNew(iEdge,:,iVertexOnEdge), &
                    yLocalTriNew(iEdge,:,iVertexOnEdge), &
                    rotateCartesianGrid, &
                    includeMetricTerms, &
                    on_a_sphere, &
                    integrationType, &
                    integrationOrder, &
                    sphere_radius, &
                    basisGradientUTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisGradientVTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsUTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsVTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsMetricTriNew(iEdge,:,:,iVertexOnEdge))

               end do

          end do

       end if

    else if (trim(variationalBasisType) == "pwl") then

       if ( .not. useCGrid ) then
          call seaice_init_velocity_solver_pwl(&
               nCells, &
               maxEdges, &
               nEdgesOnCell, &
               verticesOnCell, &  !NOTE: this is actually not used
               edgesOnCell, &     !NOTE: this can be avoided if c in the function is computed
               dvEdge, &          !NOTE: this can be avoided if c in the function is computed 
               areaCell, &
               xLocal, &
               yLocal, &
               basisGradientU, &
               basisGradientV, &
               basisIntegralsMetric, &
               basisIntegralsU, &
               basisIntegralsV)
       else

          do iEdge = 1, nEdges

             do iCellOnEdge = 1, 2

                iCell = cellsOnEdge(iCellOnEdge,iEdge)

                if (iCell <= nCells) then

                   call seaice_init_velocity_solver_pwl_c_grid(&
                      iCell, &
                      maxEdges, &
                      nEdgesOnCell, &
                      areaCellCGrid, &
                      xLocalNew(iEdge,:,iCellOnEdge), &
                      yLocalNew(iEdge,:,iCellOnEdge), &
                      basisGradientUNew(iEdge,:,:,iCellOnEdge), &
                      basisGradientVNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsMetricNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsUNew(iEdge,:,:,iCellOnEdge), &
                      basisIntegralsVNew(iEdge,:,:,iCellOnEdge))

                end if

             end do

             do iVertexOnEdge = 1, 2

                 iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

                 call seaice_init_velocity_solver_pwl_c_grid(&
                    iVertex, &
                    vertexDegree, &
                    nEdgesOnVertex, &
                    areaTriCGrid, &
                    xLocalTriNew(iEdge,:,iVertexOnEdge), &
                    yLocalTriNew(iEdge,:,iVertexOnEdge), &
                    basisGradientUTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisGradientVTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsMetricTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsUTriNew(iEdge,:,:,iVertexOnEdge), &
                    basisIntegralsVTriNew(iEdge,:,:,iVertexOnEdge))

              end do

          end do

       end if

    else if (trim(variationalBasisType) == "none") then

       continue

    else

       call MPAS_log_write("Unknown variational basis type: "//trim(variationalBasisType), MPAS_LOG_CRIT)

    endif

    if ( .not. useCGrid ) then
       call mpas_timer_start("variational denominator")
       call variational_denominator(&
            nVertices, &
            vertexDegree, &
            nEdgesOnCell, &
            interiorVertex, &
            areaTriangle, &
            cellsOnVertex, &
            cellVerticesAtVertex, &
            basisIntegralsMetric, &
            variationalDenominatorType, &
            variationalDenominator)
       call mpas_timer_stop("variational denominator")
    else
       call mpas_timer_start("variational denominator_c_grid")
       call variational_denominator_c_grid(&
          on_a_sphere, &
          sphere_radius, &
          cellsOnEdge, &
          edgesOnCell, &
          verticesOnEdge, &
          edgesOnVertex, &
          nEdges, &
          nCells, &
          nVertices, &
          nEdgesOnCell, &
          vertexDegree, &
          xCell, &
          yCell, &
          zCell, &
          xVertex, &
          yVertex, &
          zVertex, &
          xEdge, &
          yEdge, &
          zEdge, &
          latCell, &
          lonCell, &
          latVertex, &
          lonVertex, &
          latEdge, &
          lonEdge, &
          variationalDenominatorCGrid, &
          variationalDenominatorType, &
          basisIntegralsMetricNew, &
          basisIntegralsMetricTriNew, &
          cellEdgesAtEdge, &
          triangleEdgesAtEdge)
       call mpas_timer_stop("variational denominator_c_grid")
    end if

    ! clean up
    if (.not. useCGrid) then
       if (.not. useProjectedBGrid) then 
           deallocate(xLocal)
           deallocate(yLocal)
       else
           deallocate(xLocalNewNew)
           deallocate(yLocalNewNew)
       endif
           
    else 
       deallocate(xLocalNew)
       deallocate(yLocalNew)
       deallocate(xLocalTriNew)
       deallocate(yLocalTriNew)
       deallocate(nEdgesOnVertex)
    end if

    !call homogenize_variational_basis_field()

  end subroutine init_velocity_solver_variational_primary_mesh

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  variational_denominator
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine variational_denominator(&
       nVertices, &
       vertexDegree, &
       nEdgesOnCell, &
       interiorVertex, &
       areaTriangle, &
       cellsOnVertex, &
       cellVerticesAtVertex, &
       basisIntegralsMetric, &
       variationalDenominatorType, &
       variationalDenominator)

    integer, intent(in) :: &
         nVertices, & !< Input:
         vertexDegree !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell, & !< Input:
         interiorVertex  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaTriangle !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellsOnVertex !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisIntegralsMetric !< Input:

    character(len=*), intent(in) :: &
         variationalDenominatorType !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         variationalDenominator

    integer :: &
         iVertex, &
         iSurroundingCell, &
         iStressVertex, &
         iCell, &
         iVelocityVertex

    if (trim(variationalDenominatorType) == "alternate") then

       do iVertex = 1, nVertices

          variationalDenominator(iVertex) = 0.0_RKIND

          ! loop over surrounding cells
          do iSurroundingCell = 1, vertexDegree

             ! get the cell number of this cell
             iCell = cellsOnVertex(iSurroundingCell, iVertex)

             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             iVelocityVertex = cellVerticesAtVertex(iSurroundingCell,iVertex)

             ! loop over the vertices of the surrounding cell
             do iStressVertex = 1, nEdgesOnCell(iCell)

                variationalDenominator(iVertex) = variationalDenominator(iVertex) + &
                     basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell)

             enddo ! iStressVertex

          enddo ! iSurroundingCell

          ! inverse
          variationalDenominator(iVertex) = variationalDenominator(iVertex)

       enddo ! iVertex

    else if (trim(variationalDenominatorType) == "original") then

       do iVertex = 1, nVertices
          variationalDenominator(iVertex) = areaTriangle(iVertex)
       enddo ! iVertex

    else

       call MPAS_log_write("Unknown variational denominator type: "//trim(variationalDenominatorType), MPAS_LOG_CRIT)

    endif

  end subroutine variational_denominator

!-----------------------------------------------------------------------

  subroutine compute_c_grid_polygon_and_tri_areas(&
       on_a_sphere, &
       sphere_radius, &
       edgesOnCell, &
       edgesOnVertex, &
       nEdges, &
       nCells, &
       nVertices, &
       nEdgesOnCell, &
       vertexDegree, &
       xCell, &
       yCell, &
       zCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       xEdge, &
       yEdge, &
       zEdge, &
       latCell, &
       lonCell, &
       latVertex, &
       lonVertex, &
       latEdge, &
       lonEdge, &
       areaCellCGrid, &
       areaTriCGrid)

    use seaice_mesh, only: &
         spherical_triangle_area, &
         planar_triangle_area

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

    integer, intent(in) :: &
         nEdges, &     !< Input:
         nCells, &     !< Input:
         nVertices, &  !< Input:
         vertexDegree  !< Input:

    real(kind=RKIND), intent(in) :: &
         sphere_radius !<Input:

    logical, intent(in) :: &
         on_a_sphere !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell     !<Input:

    integer, dimension(:,:), intent(in) :: &
         edgesOnCell, &     !<Input:
         edgesOnVertex      !<Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell,         & !< Input:
         yCell,         & !< Input:
         zCell,         & !< Input:
         xVertex,       & !< Input:
         yVertex,       & !< Input:
         zVertex,       & !< Input:
         xEdge,         & !< Input:
         yEdge,         & !< Input:
         zEdge,         & !< Input:
         latCell,       & !< Input:
         lonCell,       & !< Input:
         latVertex,     & !< Input:
         lonVertex,     & !< Input:
         latEdge,       & !< Input:
         lonEdge          !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         areaCellCGrid, & !< Output:
         areaTriCGrid     !< Output:

    integer :: &
         iEdge, &
         jEdge, &
         iCell, &
         iVertex, &
         iEdgeOnVertex, &
         jEdgeOnVertex, &
         iCellOnEdge, &
         iTriangleOnEdge, &
         cellOnEdge, &
         vertOnEdge1, &
         vertOnEdge2, &
         iStressEdge, &
         iVelocityEdge

    real(kind=RKIND) :: &
         laCell, &
         loCell, &
         laVertex, &
         loVertex, &
         laVert1, &
         loVert1, &
         laVert2, &
         loVert2, &
         laEdge1, &
         loEdge1, &
         laEdge2, &
         loEdge2, &
         xc, &
         yc, &
         zc, &
         xv1, &
         yv1, &
         zv1, &
         xv2, &
         yv2, &
         zv2, &
         xe1, &
         ye1, &
         ze1, &
         xe2, &
         ye2, &
         ze2, &
         xe3, &
         ye3, &
         ze3
     
    real(kind=RKIND) :: &
         triangleArea 

    areaCellCGrid(:) = 0.0_RKIND
    areaTriCGrid(:) = 0.0_RKIND

    if(on_a_sphere) then

       do iCell = 1, nCells

          do iEdge = 1, nEdgesOnCell(iCell)

             jEdge = seaice_wrapped_index(iEdge + 1, nEdgesOnCell(iCell))

             laEdge1 = latEdge(edgesOnCell(iEdge,iCell))
             loEdge1 = lonEdge(edgesOnCell(iEdge,iCell))
             laEdge2 = latEdge(edgesOnCell(jEdge,iCell))
             loEdge2 = lonEdge(edgesOnCell(jEdge,iCell))
             laCell = latCell(iCell)
             loCell = lonCell(iCell)

             areaCellCGrid(iCell) = areaCellCGrid(iCell) + spherical_triangle_area(laCell, loCell, laEdge1, loEdge1, laEdge2, loEdge2, sphere_radius)

          end do

       end do 

       do iVertex = 1, nVertices

          do iEdgeOnVertex = 1, vertexDegree

             jEdgeOnVertex = seaice_wrapped_index(iEdgeOnVertex + 1, vertexDegree)

             laEdge1 = latEdge(edgesOnVertex(iEdgeOnVertex,iVertex))
             loEdge1 = lonEdge(edgesOnVertex(iEdgeOnVertex,iVertex))
             laEdge2 = latEdge(edgesOnVertex(jEdgeOnVertex,iVertex))
             loEdge2 = lonEdge(edgesOnVertex(jEdgeOnVertex,iVertex))
             laVertex = latVertex(iVertex)
             loVertex = lonVertex(iVertex)

             areaTriCGrid(iVertex) = areaTriCGrid(iVertex) + spherical_triangle_area(laVertex, loVertex, laEdge1, loEdge1, laEdge2, loEdge2, sphere_radius)       

          end do

       end do

    else

       do iCell = 1, nCells

          do iEdge = 1, nEdgesOnCell(iCell)

             jEdge = seaice_wrapped_index(iEdge + 1, nEdgesOnCell(iCell))

             xe1 = xEdge(edgesOnCell(iEdge,iCell))
             ye1 = yEdge(edgesOnCell(iEdge,iCell)) 
             ze1 = zEdge(edgesOnCell(iEdge,iCell))
             xe2 = xEdge(edgesOnCell(jEdge,iCell))
             ye2 = yEdge(edgesOnCell(jEdge,iCell))
             ze2 = zEdge(edgesOnCell(jEdge,iCell))
             xc = xCell(iCell)
             yc = yCell(iCell)
             zc = zCell(iCell)

             areaCellCGrid(iCell) = areaCellCGrid(iCell) + planar_triangle_area(xc, yc, zc, xe1, ye1, ze1, xe2, ye2, ze2) 

          end do

       end do

       do iVertex = 1, nVertices

          do iEdgeOnVertex = 1, vertexDegree

             jEdgeOnVertex = seaice_wrapped_index(iEdgeOnVertex + 1, vertexDegree)

             xe1 = xEdge(edgesOnVertex(iEdgeOnVertex,iVertex))
             ye1 = yEdge(edgesOnVertex(iEdgeOnVertex,iVertex))
             ze1 = zEdge(edgesOnVertex(iEdgeOnVertex,iVertex))
             xe2 = xEdge(edgesOnVertex(jEdgeOnVertex,iVertex))
             ye2 = yEdge(edgesOnVertex(jEdgeOnVertex,iVertex))
             ze2 = zEdge(edgesOnVertex(jEdgeOnVertex,iVertex))
             xc = xVertex(iVertex)
             yc = yVertex(iVertex)
             zc = zVertex(iVertex)

             areaTriCGrid(iVertex) = areaTriCGrid(iVertex) + planar_triangle_area(xc, yc, zc, xe1, ye1, ze1, xe2, ye2, ze2)

          end do

       end do

    end if

end subroutine compute_c_grid_polygon_and_tri_areas

!-----------------------------------------------------------------------

  subroutine variational_denominator_c_grid(&
       on_a_sphere, &
       sphere_radius, &
       cellsOnEdge, &
       edgesOnCell, &
       verticesOnEdge, &
       edgesOnVertex, &
       nEdges, &
       nCells, &
       nVertices, &
       nEdgesOnCell, &
       vertexDegree, &
       xCell, &
       yCell, &
       zCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       xEdge, &
       yEdge, &
       zEdge, &
       latCell, &
       lonCell, &
       latVertex, &
       lonVertex, &
       latEdge, &
       lonEdge, &
       variationalDenominatorCGrid, &
       variationalDenominatorType, &
       basisIntegralsMetricNew, &
       basisIntegralsMetricTriNew, &
       cellEdgesAtEdge, &
       triangleEdgesAtEdge)

    use seaice_mesh, only: &
         spherical_triangle_area, &
         planar_triangle_area

    use seaice_velocity_solver_variational_shared, only: &
         seaice_wrapped_index

    integer, intent(in) :: &
         nEdges, &     !< Input:
         nCells, &     !< Input:
         nVertices, &  !< Input:
         vertexDegree  !< Input:

    real(kind=RKIND), intent(in) :: &
         sphere_radius !<Input:

    integer, dimension(:,:), intent(in) :: &
         cellEdgesAtEdge,     &   !< Input:
         triangleEdgesAtEdge      !< Input:

    character(len=*), intent(in) :: &
         variationalDenominatorType !< Input:

    real(kind=RKIND), dimension(:,:,:,:), intent(in) :: &
         basisIntegralsMetricNew, & !< Input:
         basisIntegralsMetricTriNew !< Input:

    logical, intent(in) :: &
         on_a_sphere !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell     !<Input:

    integer, dimension(:,:), intent(in) :: &
         cellsOnEdge, &     !<Input:
         edgesOnCell, &     !<Input:
         verticesOnEdge, &  !<Input:
         edgesOnVertex      !<Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell,         & !< Input:
         yCell,         & !< Input:
         zCell,         & !< Input:
         xVertex,       & !< Input:
         yVertex,       & !< Input:
         zVertex,       & !< Input:
         xEdge,         & !< Input:
         yEdge,         & !< Input:
         zEdge,         & !< Input:
         latCell,       & !< Input:
         lonCell,       & !< Input:
         latVertex,     & !< Input:
         lonVertex,     & !< Input:
         latEdge,       & !< Input:
         lonEdge          !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         variationalDenominatorCGrid !< Output:

    integer :: &
         iEdge, &
         jEdge, &
         iCell, &
         iVertex, &
         iEdgeOnVertex, &
         jEdgeOnVertex, &
         iCellOnEdge, &
         iTriangleOnEdge, &
         cellOnEdge, &
         vertOnEdge1, &
         vertOnEdge2, &
         iStressEdge, &
         iVelocityEdge

    real(kind=RKIND) :: &
         laCell, &
         loCell, &
         laVertex, &
         loVertex, &
         laVert1, &
         loVert1, &
         laVert2, &
         loVert2, &
         laEdge1, &
         loEdge1, &
         laEdge2, &
         loEdge2, &
         xc, &
         yc, &
         zc, &
         xv1, &
         yv1, &
         zv1, &
         xv2, &
         yv2, &
         zv2, &
         xe1, &
         ye1, &
         ze1, &
         xe2, &
         ye2, &
         ze2, &
         xe3, &
         ye3, &
         ze3
     
    real(kind=RKIND) :: &
         triangleArea 

    variationalDenominatorCGrid(:) = 0.0_RKIND

    if(on_a_sphere) then

       if (trim(variationalDenominatorType) == "original") then

          do iEdge = 1, nEdges

             vertOnEdge1 = verticesOnEdge(1, iEdge)
             vertOnEdge2 = verticesOnEdge(2, iEdge)
             laVert1 = latVertex(vertOnEdge1)
             loVert1 = lonVertex(vertOnEdge1)
             laVert2 = latVertex(vertOnEdge2)
             loVert2 = lonVertex(vertOnEdge2)

             do iCellOnEdge = 1, 2

                cellOnEdge = cellsOnEdge(iCellOnEdge, iEdge)
                laCell = latCell(cellOnEdge)             
                loCell = lonCell(cellOnEdge)

                triangleArea = spherical_triangle_area(laCell, loCell, laVert1, loVert1, laVert2, loVert2, sphere_radius)                 
                variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + triangleArea 

             end do
          
          end do

       else if (trim(variationalDenominatorType) == "alternate") then

          do iEdge = 1, nEdges

             do iCellOnEdge = 1, 2

                iCell = cellsOnEdge(iCellOnEdge, iEdge)

                iVelocityEdge = cellEdgesAtEdge(iCellOnEdge,iEdge)

                ! loop over the vertices of the surrounding cell
                do iStressEdge = 1, nEdgesOnCell(iCell)
 
                   variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge)

                end do ! iStressEdge

             end do ! iSurroundingCell

             do iTriangleOnEdge = 1, 2

                iVelocityEdge = triangleEdgesAtEdge(iTriangleOnEdge,iEdge)

                do iStressEdge = 1, vertexDegree

                   variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iTriangleOnEdge)

                end do

             end do

          end do ! iEdge

       else

          call MPAS_log_write("Unknown variational denominator type: "//trim(variationalDenominatorType), MPAS_LOG_CRIT)

       end if

    else

       if (trim(variationalDenominatorType) == "original") then

          do iEdge = 1, nEdges

             vertOnEdge1 = verticesOnEdge(1, iEdge)
             vertOnEdge2 = verticesOnEdge(2, iEdge)
             xv1 = xVertex(vertOnEdge1)
             yv1 = yVertex(vertOnEdge1)
             zv1 = zVertex(vertOnEdge1)
             xv2 = xVertex(vertOnEdge2)
             yv2 = yVertex(vertOnEdge2)
             zv2 = zVertex(vertOnEdge2)

             do iCellOnEdge = 1, 2

                cellOnEdge = cellsOnEdge(iCellOnEdge, iEdge)
                xc = xCell(cellOnEdge)
                yc = yCell(cellOnEdge)
                zc = zCell(cellOnEdge)

               triangleArea = planar_triangle_area(xc, yc, zc, xv1, yv1, zv1, xv2, yv2, zv2)
               variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + triangleArea

             end do

         end do

       else if (trim(variationalDenominatorType) == "alternate") then

          do iEdge = 1, nEdges

             do iCellOnEdge = 1, 2

                iCell = cellsOnEdge(iCellOnEdge, iEdge)

                iVelocityEdge = cellEdgesAtEdge(iCellOnEdge,iEdge)

                ! loop over the vertices of the surrounding cell
                do iStressEdge = 1, nEdgesOnCell(iCell)

                   variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge)

                end do ! iStressEdge

             end do ! iSurroundingCell

             do iTriangleOnEdge = 1, 2

                iVelocityEdge = triangleEdgesAtEdge(iTriangleOnEdge,iEdge)

                do iStressEdge = 1, vertexDegree

                   variationalDenominatorCGrid(iEdge) = variationalDenominatorCGrid(iEdge) + basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iTriangleOnEdge)

                end do

             end do

          end do ! iEdge

       else

          call MPAS_log_write("Unknown variational denominator type: "//trim(variationalDenominatorType), MPAS_LOG_CRIT)

       end if

    end if

end subroutine variational_denominator_c_grid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  homogenize_variational_basis_field
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine homogenize_variational_basis_field()

    use seaice_mesh_pool, only: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV, &
         nCells

    integer :: &
         iCell

    integer, parameter :: &
         iCellHomogenize = 1111

    !call homogenize_cell(basisGradientU, iCellHomogenize)
    !call homogenize_cell(basisGradientV, iCellHomogenize)
    !call homogenize_cell(basisIntegralsMetric, iCellHomogenize)
    !call homogenize_cell(basisIntegralsU, iCellHomogenize)
    !call homogenize_cell(basisIntegralsV, iCellHomogenize)

    do iCell = 1, nCells

       basisGradientU(:,:,iCell)       = basisGradientU(:,:,iCellHomogenize)
       basisGradientV(:,:,iCell)       = basisGradientV(:,:,iCellHomogenize)
       basisIntegralsMetric(:,:,iCell) = basisIntegralsMetric(:,:,iCellHomogenize)
       basisIntegralsU(:,:,iCell)      = basisIntegralsU(:,:,iCellHomogenize)
       basisIntegralsV(:,:,iCell)      = basisIntegralsV(:,:,iCellHomogenize)

    enddo ! iCell

  end subroutine homogenize_variational_basis_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  homogenize_variational_basis_field
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine homogenize_cell(&
       field, &
       iCell)

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         field

    integer, intent(in) :: &
         iCell

    integer :: &
         iVertexOnCell1, &
         iVertexOnCell2, &
         nCanonicalValues, &
         iCV

    real(kind=RKIND), dimension(100) :: &
         canonicalValues

    logical :: &
         lFound

    nCanonicalValues = 0

    do iVertexOnCell1 = 1, 6
       do iVertexOnCell2 = 1, 6

          lFound = .false.
          do iCV = 1, nCanonicalValues

             if (abs(abs(canonicalValues(iCV) - abs(field(iVertexOnCell1,iVertexOnCell2,iCell)))) < 1e-12) then
                lFound = .true.
                field(iVertexOnCell1,iVertexOnCell2,iCell) = canonicalValues(iCV) * sign(1.0_RKIND, field(iVertexOnCell1,iVertexOnCell2,iCell))
                exit
             endif

          enddo ! iCV

          if (.not. lFound) then
             nCanonicalValues = nCanonicalValues + 1
             canonicalValues(nCanonicalValues) = abs(field(iVertexOnCell1,iVertexOnCell2,iCell))
          endif

       enddo ! iVertexOnCell2
    enddo ! iVertexOnCell1

    !do iVertexOnCell1 = 1, 6
    !   do iVertexOnCell2 = 1, 6
    !      write(*,*) iVertexOnCell1, iVertexOnCell2, field(iVertexOnCell1,iVertexOnCell2,iCell)
    !   enddo ! iVertexOnCell2
    !enddo ! iVertexOnCell1

  end subroutine homogenize_cell

!-----------------------------------------------------------------------
! time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_strain_tensor_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_strain_tensor_variational(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       uVelocity, &
       vVelocity, &
       basisGradientU, &
       basisGradientV, &
       tanLatVertexRotatedOverRadius, &
       solveStress)!{{{

    use seaice_mesh_pool, only: &
         nCells, &
         verticesOnCell, &
         nEdgesOnCell

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, & !< Output:
         strain22, & !< Output:
         strain12    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input:
         vVelocity, & !< Input:
         tanLatVertexRotatedOverRadius !< Input:

    real(kind=RKIND), dimension(:,:,:), contiguous, intent(in) :: &
         basisGradientU, & !< Input:
         basisGradientV    !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    integer :: &
         iCell, &
         iGradientVertex, &
         iBasisVertex, &
         iVertex, &
         jVertex

    real(kind=RKIND) :: &
         strain11Tmp, &
         strain22Tmp, &
         strain12Tmp

    ! loop over cells
!$omp parallel do default(shared) private(iGradientVertex, iBasisVertex, iVertex, jVertex, &
!$omp&                                    strain11Tmp, strain22Tmp, strain12Tmp)

    do iCell = 1, nCells
   
       if (solveStress(iCell) == 1) then

          ! loop over velocity points surrounding cell - location of stress and derivative
          do iGradientVertex = 1, nEdgesOnCell(iCell)

             strain11Tmp = 0.0_RKIND
             strain22Tmp = 0.0_RKIND
             strain12Tmp = 0.0_RKIND

             ! loop over basis functions
             do iBasisVertex = 1, nEdgesOnCell(iCell)

                iVertex = verticesOnCell(iBasisVertex,iCell)

                strain11Tmp = strain11Tmp + uVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell)
                strain22Tmp = strain22Tmp + vVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell)
                strain12Tmp = strain12Tmp + 0.5_RKIND * (&
                            uVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell) + &
                            vVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell))

             enddo ! iVertexOnCell

             ! metric terms
             jVertex = verticesOnCell(iGradientVertex,iCell)

             strain11(iGradientVertex,iCell) = strain11Tmp - vVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex)
             strain12(iGradientVertex,iCell) = strain12Tmp + uVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex) * 0.5_RKIND
             strain22(iGradientVertex,iCell) = strain22Tmp

          enddo ! jVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine seaice_strain_tensor_variational!}}}

!-----------------------------------------------------------------------

  subroutine seaice_strain_tensor_variational_c_grid(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       strain11Tri, &
       strain22Tri, &
       strain12Tri, &
       uVelocityCGrid, &
       vVelocityCGrid, &
       basisGradientUNew, &
       basisGradientVNew, &
       basisGradientUTriNew, &
       basisGradientVTriNew, &
       tanLatEdgeRotatedOverRadius, &
       solveStress, &
       solveStressTri, &
       solveVelocityCGrid)!{{{

    use seaice_mesh_pool, only: &
         nCells, &
         nEdges, &
         nVertices, & 
         vertexDegree, &
         cellsOnEdge, &
         edgesOnCell, &
         edgesOnVertex, &
         nEdgesOnCell, &
         verticesOnEdge

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11,    & !< Output:
         strain22,    & !< Output:
         strain12,    &!< Output:
         strain11Tri, & !< Output:
         strain22Tri, & !< Output:
         strain12Tri    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         tanLatEdgeRotatedOverRadius !< Input:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         uVelocityCGrid, & !< Input:
         vVelocityCGrid    !< Input:

    real(kind=RKIND), dimension(:,:,:,:), contiguous, intent(in) :: &
         basisGradientUNew, & !< Input:
         basisGradientVNew, & !< Input:
         basisGradientUTriNew, & !< Input:
         basisGradientVTriNew    !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress, &  !< Input:
         solveStressTri, &  !< Input:
         solveVelocityCGrid !< Input:

    integer :: &
         iCellOnEdge, &
         iVertexOnEdge, &
         iCell, &
         iGradientEdge, &
         iBasisEdge, &
         iEdge, &
         jEdge, &
         kEdge, &
         vertexOnEdge1, & 
         vertexOnEdge2, &
         iVertex

    real(kind=RKIND) :: &
         strain11Tmp, &
         strain22Tmp, &
         strain12Tmp


    do iEdge = 1, nEdges

       do iCellOnEdge = 1, 2
          iCell = cellsOnEdge(iCellOnEdge,iEdge)
   
          if (solveStress(iCell) == 1) then

             ! loop over velocity points surrounding cell - location of stress and derivative
             do iGradientEdge = 1, nEdgesOnCell(iCell)

                strain11Tmp = 0.0_RKIND
                strain22Tmp = 0.0_RKIND
                strain12Tmp = 0.0_RKIND

                ! loop over basis functions
                do iBasisEdge = 1, nEdgesOnCell(iCell)

                   jEdge = edgesOnCell(iBasisEdge,iCell)

                   strain11Tmp = strain11Tmp + uVelocityCGrid(jEdge) * basisGradientUNew(iEdge,iBasisEdge,iGradientEdge,iCellOnEdge)
                   strain22Tmp = strain22Tmp + vVelocityCGrid(jEdge) * basisGradientVNew(iEdge,iBasisEdge,iGradientEdge,iCellOnEdge)
                   strain12Tmp = strain12Tmp + 0.5_RKIND * (&
                              uVelocityCGrid(jEdge) * basisGradientVNew(iEdge,iBasisEdge,iGradientEdge,iCellOnEdge) + &
                              vVelocityCGrid(jEdge) * basisGradientUNew(iEdge,iBasisEdge,iGradientEdge,iCellOnEdge))

                 end do ! jEdgeOnCell

                 ! metric terms
                 kEdge = edgesOnCell(iGradientEdge,iCell)

                  strain11(iGradientEdge,iCell) = strain11Tmp - vVelocityCGrid(kEdge) * tanLatEdgeRotatedOverRadius(kEdge)
                  strain12(iGradientEdge,iCell) = strain12Tmp + uVelocityCGrid(kEdge) * tanLatEdgeRotatedOverRadius(kEdge) * 0.5_RKIND
                  strain22(iGradientEdge,iCell) = strain22Tmp

               end do ! kEdgeOnCell

           end if !solveStress

       end do

       do iVertexOnEdge = 1, 2

          iVertex = verticesOnEdge(iVertexOnEdge,iEdge)
          if (solveStressTri(iVertex) == 1) then

             ! loop over velocity points surrounding triangle - location of stress and derivative
             do iGradientEdge = 1, vertexDegree

                strain11Tmp = 0.0_RKIND
                strain22Tmp = 0.0_RKIND
                strain12Tmp = 0.0_RKIND

                ! loop over basis functions
                do iBasisEdge = 1, vertexDegree

                   jEdge = edgesOnVertex(iBasisEdge,iVertex)
   
                   strain11Tmp = strain11Tmp + uVelocityCGrid(jEdge) * basisGradientUTriNew(iEdge,iBasisEdge,iGradientEdge,iVertexOnEdge)
                   strain22Tmp = strain22Tmp + vVelocityCGrid(jEdge) * basisGradientVTriNew(iEdge,iBasisEdge,iGradientEdge,iVertexOnEdge)
                   strain12Tmp = strain12Tmp + 0.5_RKIND * (&
                              uVelocityCGrid(jEdge) * basisGradientVTriNew(iEdge,iBasisEdge,iGradientEdge,iVertexOnEdge) + &
                              vVelocityCGrid(jEdge) * basisGradientUTriNew(iEdge,iBasisEdge,iGradientEdge,iVertexOnEdge))

                 end do ! jEdgeOnVertex

                 ! metric terms
                 kEdge = edgesOnVertex(iGradientEdge,iVertex)

                 strain11Tri(iGradientEdge,iVertex) = strain11Tmp - vVelocityCGrid(kEdge) * tanLatEdgeRotatedOverRadius(kEdge)
                 strain12Tri(iGradientEdge,iVertex) = strain12Tmp + uVelocityCGrid(kEdge) * tanLatEdgeRotatedOverRadius(kEdge) * 0.5_RKIND
                 strain22Tri(iGradientEdge,iVertex) = strain22Tmp

              end do ! kEdgeOnVertex

          end if !solveStressTri

       end do ! iVertex

    end do !iEdge

  end subroutine seaice_strain_tensor_variational_c_grid !}}}

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_average_strains_on_vertex
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_average_strains_on_vertex(&
       areaCell, &
       strain11, &
       strain22, &
       strain12)

    use seaice_mesh_pool, only: &
         nCells, &
         nVerticesSolve, &
         cellVerticesAtVertex, &
         cellsOnVertex, &
         vertexDegree

    real(kind=RKIND), dimension(:), intent(in) :: &
         areaCell !< Input

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         strain11, & !< Input/Output:
         strain22, & !< Input/Output:
         strain12    !< Input/Output:

    real(kind=RKIND) :: &
         strain11avg, &
         strain22avg, &
         strain12avg, &
         denominator

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         iVertexOnCell

    do iVertex = 1, nVerticesSolve

       strain11avg = 0.0_RKIND
       strain22avg = 0.0_RKIND
       strain12avg = 0.0_RKIND
       denominator = 0.0_RKIND

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          if (iCell <= nCells) then

             iVertexOnCell = cellVerticesAtVertex(iVertexDegree,iVertex)

             strain11avg = strain11avg + strain11(iVertexOnCell,iCell) * areaCell(iCell)
             strain22avg = strain22avg + strain22(iVertexOnCell,iCell) * areaCell(iCell)
             strain12avg = strain12avg + strain12(iVertexOnCell,iCell) * areaCell(iCell)
             denominator = denominator + areaCell(iCell)

          endif

       enddo ! iVertexDegree

       strain11avg = strain11avg / denominator
       strain22avg = strain22avg / denominator
       strain12avg = strain12avg / denominator

       do iVertexDegree = 1, vertexDegree

          iCell = cellsOnVertex(iVertexDegree,iVertex)

          if (iCell <= nCells) then

             iVertexOnCell = cellVerticesAtVertex(iVertexDegree,iVertex)

             strain11(iVertexOnCell,iCell) = strain11avg
             strain22(iVertexOnCell,iCell) = strain22avg
             strain12(iVertexOnCell,iCell) = strain12avg

          endif

       enddo ! iCellOnVertex

    enddo ! iVertex

  end subroutine seaice_average_strains_on_vertex

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_tensor_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_tensor_variational(&
       mesh, &
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       solveStress, &
       dtElastic)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         constitutiveRelationType, &
         EVP_CONSTITUTIVE_RELATION, &
         REVISED_EVP_CONSTITUTIVE_RELATION, &
         LINEAR_CONSTITUTIVE_RELATION, &
         seaice_evp_constitutive_relation, &
         seaice_evp_constitutive_relation_revised, &
         seaice_linear_constitutive_relation, &
         eccentricitySquared, puny, dampingTimescale

    use seaice_mesh_pool, only: &
         nCells, &
         nEdgesOnCell

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), contiguous, intent(inout) :: &
         stress11, & !< Input/Output:
         stress22, & !< Input/Output:
         stress12    !< Input/Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(out) :: &
         replacementPressure !< Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input:

    integer :: &
         iCell, &
         iVertexOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    real(kind=RKIND) :: &
         strainDivergence,    &
         strainTension,       &
         strainShearing,      &
         stress1,             &
         stress2,             &
         Delta,               &
         pressureCoefficient, &
         denominator

    ! init variables
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    if (constitutiveRelationType == EVP_CONSTITUTIVE_RELATION) then

       denominator = 1.0_RKIND + (0.5_RKIND * dtElastic) / dampingTimescale

!$omp parallel do default(shared) private(iVertexOnCell)
       do iCell = 1, nCells

          replacementPressure(:,iCell) = 0.0_RKIND

          if (solveStress(iCell) == 1) then

             !$omp simd
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call seaice_evp_constitutive_relation(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell), &
                     icePressure(iCell), &
                     replacementPressure(iVertexOnCell,iCell), &
                     areaCell(iCell), &
                     dtElastic)
             enddo ! iVertexOnCell

          endif ! solveStress

       enddo ! iCell

    else if (constitutiveRelationType == REVISED_EVP_CONSTITUTIVE_RELATION) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             !$omp simd
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call seaice_evp_constitutive_relation_revised(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell), &
                     icePressure(iCell), &
                     replacementPressure(iVertexOnCell,iCell), &
                     areaCell(iCell))

             enddo ! iVertexOnCell

          endif ! solveStress

       enddo ! iCell

    else if (constitutiveRelationType == LINEAR_CONSTITUTIVE_RELATION) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             !$omp simd
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call seaice_linear_constitutive_relation(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell))

             enddo ! iVertexOnCell

          endif ! solveStress

       enddo ! iCell

    endif ! constitutiveRelationType

  end subroutine seaice_stress_tensor_variational!}}}

!-----------------------------------------------------------------------

  subroutine seaice_stress_tensor_variational_c_grid(&
       mesh, &
       velocityVariational, &
       stress11, &
       stress22, &
       stress12, &
       stress11Tri, &
       stress22Tri, &
       stress12Tri, &
       strain11, &
       strain22, &
       strain12, &
       strain11Tri, &
       strain22Tri, &
       strain12Tri, &
       icePressure, &
       replacementPressure, &
       replacementPressureTri, &
       solveStress, &
       solveStressTri, &
       dtElastic)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         constitutiveRelationType, &
         EVP_CONSTITUTIVE_RELATION, &
         REVISED_EVP_CONSTITUTIVE_RELATION, &
         LINEAR_CONSTITUTIVE_RELATION, &
         seaice_evp_constitutive_relation, &
         seaice_evp_constitutive_relation_revised, &
         seaice_linear_constitutive_relation, &
         eccentricitySquared, puny, dampingTimescale

    use seaice_mesh_pool, only: &
         nCells, &
         nVertices, &
         nEdgesOnCell, &
         cellsOnVertex, &
         vertexDegree, &
         kiteAreasOnVertex, &
         areaTriangle

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh, & !< Input:
         velocityVariational !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress, & !< Input:
         solveStressTri !< Input:

    real(kind=RKIND), dimension(:,:), contiguous, intent(inout) :: &
         stress11,    & !< Input/Output:
         stress22,    & !< Input/Output:
         stress12,    & !< Input/Output:
         stress11Tri, & !< Input/Output:
         stress22Tri, & !< Input/Output:
         stress12Tri    !< Input/Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(out) :: &
         replacementPressure, &    !< Output:
         replacementPressureTri    !< Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(in) :: &
         strain11,    & !< Input:
         strain22,    & !< Input:
         strain12,    & !< Input:
         strain11Tri, & !< Input:
         strain22Tri, & !< Input:
         strain12Tri    !< Input:

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input:

    real(kind=RKIND), dimension(:), pointer :: &
         areaCellCGrid, &
         areaTriCGrid

    integer :: &
         iCell, &
         iVertex, &
         iEdgesOnCell, &
         iEdgesOnVertex, &
         iCellsOnVertex

    real(kind=RKIND), dimension(:), allocatable :: &
         icePressureCellCGrid, &  !TODO: give as input!!!!
         icePressureTriCGrid

    real(kind=RKIND) :: &
         strainDivergence,    &
         strainTension,       &
         strainShearing,      &
         stress1,             &
         stress2,             &
         Delta,               &
         pressureCoefficient, &
         denominator

    ! TODO: I removed the OpenMP parts!

    allocate(icePressureCellCGrid(nCells))
    allocate(icePressureTriCGrid(nVertices))

    ! TODO: this has to be handled properly in the future
    do iCell = 1, nCells
       icePressureCellCGrid(iCell) = icePressure(iCell)
    end do

    icePressureTriCGrid(:) = 0.0_RKIND
    do iVertex = 1, nVertices
       do iCellsOnVertex = 1, vertexDegree
          iCell = cellsOnVertex(iCellsOnVertex,iVertex)
          icePressureTriCGrid(iVertex) = icePressureTriCGrid(iVertex) + kiteAreasOnVertex(iCellsOnVertex, iVertex) * icePressure(iCell) / areaTriangle(iVertex)
       end do
    end do

    ! init variables
    call MPAS_pool_get_array(velocityVariational, "areaCellCGrid", areaCellCGrid)
    call MPAS_pool_get_array(velocityVariational, "areaTriCGrid", areaTriCGrid)

    if (constitutiveRelationType == EVP_CONSTITUTIVE_RELATION) then

       denominator = 1.0_RKIND + (0.5_RKIND * dtElastic) / dampingTimescale

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             do iEdgesOnCell = 1, nEdgesOnCell(iCell)

                call seaice_evp_constitutive_relation(&
                     stress11(iEdgesOnCell,iCell), &
                     stress22(iEdgesOnCell,iCell), &
                     stress12(iEdgesOnCell,iCell), &
                     strain11(iEdgesOnCell,iCell), &
                     strain22(iEdgesOnCell,iCell), &
                     strain12(iEdgesOnCell,iCell), &
                     icePressureCellCGrid(iCell), &
                     replacementPressure(iEdgesOnCell,iCell), &
                     areaCellCGrid(iCell), &
                     dtElastic)

              end do ! iEdgesOnCell
 
          end if 
     
       end do ! iCell

       do iVertex = 1, nVertices

          if (solveStressTri(iVertex) == 1) then

             do iEdgesOnVertex = 1, vertexDegree

                call seaice_evp_constitutive_relation(&
                     stress11Tri(iEdgesOnVertex,iVertex), &
                     stress22Tri(iEdgesOnVertex,iVertex), &
                     stress12Tri(iEdgesOnVertex,iVertex), &
                     strain11Tri(iEdgesOnVertex,iVertex), &
                     strain22Tri(iEdgesOnVertex,iVertex), &
                     strain12Tri(iEdgesOnVertex,iVertex), &
                     icePressureTriCGrid(iVertex), &
                     replacementPressureTri(iEdgesOnVertex,iVertex), &
                     areaTriCGrid(iVertex), &
                     dtElastic)

              end do ! iEdgesOnVertex

           end if

       end do ! iVertex

    else if (constitutiveRelationType == REVISED_EVP_CONSTITUTIVE_RELATION) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             do iEdgesOnCell = 1, nEdgesOnCell(iCell)

                call seaice_evp_constitutive_relation_revised(&
                     stress11(iEdgesOnCell,iCell), &
                     stress22(iEdgesOnCell,iCell), &
                     stress12(iEdgesOnCell,iCell), &
                     strain11(iEdgesOnCell,iCell), &
                     strain22(iEdgesOnCell,iCell), &
                     strain12(iEdgesOnCell,iCell), &
                     icePressureCellCGrid(iCell), &
                     replacementPressure(iEdgesOnCell,iCell), &
                     areaCellCGrid(iCell))

              end do ! iEdgesOnCell

          end if

       end do ! iCell

       do iVertex = 1, nVertices

          if (solveStressTri(iVertex) == 1) then

             do iEdgesOnVertex = 1, vertexDegree

                call seaice_evp_constitutive_relation_revised(&
                     stress11Tri(iEdgesOnVertex,iVertex), &
                     stress22Tri(iEdgesOnVertex,iVertex), &
                     stress12Tri(iEdgesOnVertex,iVertex), &
                     strain11Tri(iEdgesOnVertex,iVertex), &
                     strain22Tri(iEdgesOnVertex,iVertex), &
                     strain12Tri(iEdgesOnVertex,iVertex), &
                     icePressureTriCGrid(iVertex), &
                     replacementPressureTri(iEdgesOnVertex,iVertex), &
                     areaTriCGrid(iVertex))

              end do ! iEdgesOnVertex

          end if  

       end do ! iVertex

    else if (constitutiveRelationType == LINEAR_CONSTITUTIVE_RELATION) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             do iEdgesOnCell = 1, nEdgesOnCell(iCell)

                call seaice_linear_constitutive_relation(&
                     stress11(iEdgesOnCell,iCell), &
                     stress22(iEdgesOnCell,iCell), &
                     stress12(iEdgesOnCell,iCell), &
                     strain11(iEdgesOnCell,iCell), &
                     strain22(iEdgesOnCell,iCell), &
                     strain12(iEdgesOnCell,iCell))

             end do ! iEdgesOnCell

          end if  

       end do ! iCell

       do iVertex = 1, nVertices

          if (solveStressTri(iVertex) == 1) then

             do iEdgesOnVertex = 1, vertexDegree

                call seaice_linear_constitutive_relation(&
                     stress11Tri(iEdgesOnVertex,iVertex), &
                     stress22Tri(iEdgesOnVertex,iVertex), &
                     stress12Tri(iEdgesOnVertex,iVertex), &
                     strain11Tri(iEdgesOnVertex,iVertex), &
                     strain22Tri(iEdgesOnVertex,iVertex), &
                     strain12Tri(iEdgesOnVertex,iVertex))

             end do ! iEdgesOnVertex

           end if

       end do ! iVertex

    end if ! constitutiveRelationType

  deallocate(icePressureCellCGrid)
  deallocate(icePressureTriCGrid)

  end subroutine seaice_stress_tensor_variational_c_grid!}}}

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_tensor_variational_linear
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 3rd October 2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_tensor_variational_linear(&
       mesh, &
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       solveStress)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         seaice_linear_constitutive_relation

    use seaice_mesh_pool, only: &
         nCells, &
         nEdgesOnCell

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), contiguous, intent(out) :: &
         stress11, & !< Input/Output:
         stress22, & !< Input/Output:
         stress12    !< Input/Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    integer :: &
         iCell, &
         iVertexOnCell

    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          !$omp simd
          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             call seaice_linear_constitutive_relation(&
                  stress11(iVertexOnCell,iCell), &
                  stress22(iVertexOnCell,iCell), &
                  stress12(iVertexOnCell,iCell), &
                  strain11(iVertexOnCell,iCell), &
                  strain22(iVertexOnCell,iCell), &
                  strain12(iVertexOnCell,iCell))

          enddo ! iVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine seaice_stress_tensor_variational_linear!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_divergence_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_divergence_variational(&
       mesh, &
       stressDivergenceU, &
       stressDivergenceV, &
       stressDivergenceUNoMetric, &
       stressDivergenceVNoMetric, &
       stressDivergenceUMetric, &
       stressDivergenceVMetric, &
       stressDivergenceUCGrid, &
       stressDivergenceVCGrid, &
       stress11, &
       stress22, &
       stress12, &
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       variationalDenominator, &
       tanLatVertexRotatedOverRadius, &
       cellVerticesAtVertex, &
       solveVelocity)!{{{

    use seaice_mesh_pool, only: &
         nVerticesSolve, &
         nEdges, &
         cellsOnVertex, &
         nEdgesOnCell, &
         areaTriangle, &
         vertexDegree, &
         verticesOnEdge

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output:
         stressDivergenceV, &    !< Output:
         stressDivergenceUNoMetric, & !< Output:
         stressDivergenceVNoMetric, &    !< Output:
         stressDivergenceUMetric, & !< Output:
         stressDivergenceVMetric, &    !< Output:
         stressDivergenceUCGrid, & !< Output:
         stressDivergenceVCGrid    !< Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(in) :: &
         stress11, & !< Input:
         stress22, & !< Input:
         stress12    !< Input:

    real(kind=RKIND), dimension(:,:,:), contiguous, intent(in) :: &
         basisIntegralsU, &   !< Input:
         basisIntegralsV, &   !< Input:
         basisIntegralsMetric !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         tanLatVertexRotatedOverRadius, & !< Input:
         variationalDenominator           !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex !< Input:

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input:

    real(kind=RKIND) :: &
         stressDivergenceUVertex, &
         stressDivergenceVVertex, &
         stressDivergenceUCell, &
         stressDivergenceVCell, &
         stressDivergenceUVertexNoMetric, &
         stressDivergenceVVertexNoMetric, &
         stressDivergenceUCellNoMetric, &
         stressDivergenceVCellNoMetric, &
         stressDivergenceUVertexMetric, &
         stressDivergenceVVertexMetric, &
         stressDivergenceUCellMetric, &
         stressDivergenceVCellMetric, &
         areaTriangleSum, coeff1, coeff2

    integer :: &
         iVertex, &
         iSurroundingCell, &
         iCell, &
         iEdge, &
         iStressVertex, &
         iVelocityVertex, &
         vert1, &
         vert2

    ! loop over velocity positions
!$omp parallel do default(shared) private(stressDivergenceUVertex, stressDivergenceVVertex, &
!$omp&   iSurroundingCell, iCell, iVelocityVertex, stressDivergenceUCell, stressDivergenceVCell, iStressVertex)
    do iVertex = 1, nVerticesSolve

       if (solveVelocity(iVertex) == 1) then

          stressDivergenceUVertex = 0.0_RKIND
          stressDivergenceVVertex = 0.0_RKIND
          stressDivergenceUVertexNoMetric = 0.0_RKIND
          stressDivergenceVVertexNoMetric = 0.0_RKIND
          stressDivergenceUVertexMetric = 0.0_RKIND
          stressDivergenceVVertexMetric = 0.0_RKIND

          ! loop over surrounding cells
          do iSurroundingCell = 1, vertexDegree

             ! get the cell number of this cell
             iCell = cellsOnVertex(iSurroundingCell, iVertex)

             ! get the vertexOnCell number of the iVertex velocity point from cell iCell
             iVelocityVertex = cellVerticesAtVertex(iSurroundingCell,iVertex)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND
             stressDivergenceUCellNoMetric = 0.0_RKIND
             stressDivergenceVCellNoMetric = 0.0_RKIND
             stressDivergenceUCellMetric = 0.0_RKIND
             stressDivergenceVCellMetric = 0.0_RKIND

             ! loop over the vertices of the surrounding cell
             do iStressVertex = 1, nEdgesOnCell(iCell)

                ! normal & metric terms
                stressDivergenceUCell = stressDivergenceUCell - &
                     stress11(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                     tanLatVertexRotatedOverRadius(iVertex)

                stressDivergenceVCell = stressDivergenceVCell - &
                     stress22(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) + &
                     stress11(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                     tanLatVertexRotatedOverRadius(iVertex)

                stressDivergenceUCellNoMetric = stressDivergenceUCellNoMetric - &
                     stress11(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell)  

                stressDivergenceVCellNoMetric = stressDivergenceVCellNoMetric - &
                     stress22(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell)  

                stressDivergenceUCellMetric = stressDivergenceUCellMetric - &
                     stress12(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                     tanLatVertexRotatedOverRadius(iVertex)

                stressDivergenceVCellMetric = stressDivergenceVCellMetric + &
                     stress11(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                     tanLatVertexRotatedOverRadius(iVertex)



             enddo ! iStressVertex

             stressDivergenceUVertex = stressDivergenceUVertex + stressDivergenceUCell
             stressDivergenceVVertex = stressDivergenceVVertex + stressDivergenceVCell
             stressDivergenceUVertexNoMetric = stressDivergenceUVertexNoMetric + stressDivergenceUCellNoMetric
             stressDivergenceVVertexNoMetric = stressDivergenceVVertexNoMetric + stressDivergenceVCellNoMetric
             stressDivergenceUVertexMetric = stressDivergenceUVertexMetric + stressDivergenceUCellMetric
             stressDivergenceVVertexMetric = stressDivergenceVVertexMetric + stressDivergenceVCellMetric

          enddo ! iSurroundingCell

          stressDivergenceU(iVertex) = stressDivergenceUVertex / variationalDenominator(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceVVertex / variationalDenominator(iVertex)
          stressDivergenceUNoMetric(iVertex) = stressDivergenceUVertexNoMetric / variationalDenominator(iVertex)
          stressDivergenceVNoMetric(iVertex) = stressDivergenceVVertexNoMetric / variationalDenominator(iVertex)
          stressDivergenceUMetric(iVertex) = stressDivergenceUVertexMetric / variationalDenominator(iVertex)
          stressDivergenceVMetric(iVertex) = stressDivergenceVVertexMetric / variationalDenominator(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

    do iEdge = 1, nEdges
       vert1 = verticesOnEdge(1,iEdge)
       vert2 = verticesOnEdge(2,iEdge)
       areaTriangleSum = areaTriangle(vert1) + areaTriangle(vert2)
       coeff1 = areaTriangle(vert1) / areaTriangleSum   
       coeff2 = areaTriangle(vert2) / areaTriangleSum   
       !coeff1 = 0.5_RKIND
       !coeff2 = 0.5_RKIND
       stressDivergenceUCGrid(iEdge) = coeff1 * stressDivergenceU(vert1) + coeff2 * stressDivergenceU(vert2)
       stressDivergenceVCGrid(iEdge) = coeff1 * stressDivergenceV(vert1) + coeff2 * stressDivergenceV(vert2)
    end do


  end subroutine seaice_stress_divergence_variational!}}}

!-----------------------------------------------------------------------

  subroutine seaice_stress_divergence_variational_c_grid(&
       mesh, &
       boundary, &
       stressDivergenceU, &
       stressDivergenceV, &
       stressDivergenceUCGrid, &
       stressDivergenceVCGrid, &
       stressDivergenceUCGridNoMetric, &
       stressDivergenceVCGridNoMetric, &
       stressDivergenceUCGridMetric, &
       stressDivergenceVCGridMetric, &
       stress11, &
       stress22, &
       stress12, &
       stress11Tri, &
       stress22Tri, &
       stress12Tri, &
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
       solveVelocityCGrid)!{{{

    use seaice_mesh_pool, only: &
         nEdgesSolve, &
         nEdgesOnCell, &
         nVerticesSolve, &
         edgesOnVertex, &
         cellsOnEdge,  &
         edgesOnCell, &
         verticesOnEdge, &
         vertexDegree 

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh, &   !< Input:
         boundary  !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceUCGrid,    & !< Output:
         stressDivergenceVCGrid,    &    !< Output:
         stressDivergenceUCGridNoMetric,    & !< Output:
         stressDivergenceVCGridNoMetric,    &    !< Output:
         stressDivergenceUCGridMetric,    & !< Output:
         stressDivergenceVCGridMetric,    &    !< Output:
         stressDivergenceU,    & !< Output:
         stressDivergenceV       !< Output:

    real(kind=RKIND), dimension(:,:), contiguous, intent(in) :: &
         stress11,        & !< Input:
         stress22,        & !< Input:
         stress12,        & !< Input:
         stress11Tri,     & !< Input:
         stress22Tri,     & !< Input:
         stress12Tri        !< Input:

    real(kind=RKIND), dimension(:,:,:,:), contiguous, intent(in) :: &
         basisIntegralsUNew,        &   !< Input:
         basisIntegralsVNew,        &   !< Input:
         basisIntegralsMetricNew,   &   !< Input:
         basisIntegralsUTriNew,     &   !< Input:
         basisIntegralsVTriNew,     &   !< Input:
         basisIntegralsMetricTriNew     !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         tanLatEdgeRotatedOverRadius, & !< Input:
         variationalDenominatorCGrid     !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellEdgesAtEdge,     &   !< Input:
         triangleEdgesAtEdge      !< Input:

    integer, dimension(:), intent(in) :: &
         solveVelocityCGrid !< Input:

    real(kind=RKIND), dimension(:), pointer :: xCell, xVertex

    real(kind=RKIND) :: &
         stressDivergenceUEdge, &
         stressDivergenceVEdge, &
         stressDivergenceUCell, &
         stressDivergenceVCell, &
         stressDivergenceUEdgeNoMetric, &
         stressDivergenceVEdgeNoMetric, &
         stressDivergenceUCellNoMetric, &
         stressDivergenceVCellNoMetric, &
         stressDivergenceUEdgeMetric, &
         stressDivergenceVEdgeMetric, &
         stressDivergenceUCellMetric, &
         stressDivergenceVCellMetric

    integer :: &
         iEdge, &
         jEdge, &
         iEdgeOnVertex, &
         iCell, &
         iVertex, &
         iCellOnEdge, &
         iVertexOnEdge, &
         iVelocityEdge, &
         iStressEdge

    logical :: computeVertex
    real(kind=RKIND) :: coeff, sumOfDiamonds

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "xVertex", xVertex)

!!!!!!!!!! TODO: I REMOVED THE STUFF ON OpenMP

    do iEdge = 1, nEdgesSolve

       if (solveVelocityCGrid(iEdge)==1) then

          stressDivergenceUCGrid(iEdge) = 0.0_RKIND
          stressDivergenceVCGrid(iEdge) = 0.0_RKIND
          stressDivergenceUCGridNoMetric(iEdge) = 0.0_RKIND
          stressDivergenceVCGridNoMetric(iEdge) = 0.0_RKIND
          stressDivergenceUCGridMetric(iEdge) = 0.0_RKIND
          stressDivergenceVCGridMetric(iEdge) = 0.0_RKIND

          stressDivergenceUEdge = 0.0_RKIND
          stressDivergenceVEdge = 0.0_RKIND
          stressDivergenceUEdgeNoMetric = 0.0_RKIND
          stressDivergenceVEdgeNoMetric = 0.0_RKIND
          stressDivergenceUEdgeMetric = 0.0_RKIND
          stressDivergenceVEdgeMetric = 0.0_RKIND

          ! loop over surrounding cells
          do iCellOnEdge = 1, 2

             ! get the cell number of this cell
             iCell = cellsOnEdge(iCellOnEdge, iEdge)

             ! get the local index of iEdge on iCellOnEdge 
             iVelocityEdge = cellEdgesAtEdge(iCellOnEdge,iEdge)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND
             stressDivergenceUCellNoMetric = 0.0_RKIND
             stressDivergenceVCellNoMetric = 0.0_RKIND
             stressDivergenceUCellMetric = 0.0_RKIND
             stressDivergenceVCellMetric = 0.0_RKIND

             ! loop over the edges of the surrounding cell
             do iStressEdge = 1, nEdgesOnCell(iCell)

                jEdge = edgesOnCell(iStressEdge,iCell)                 
                ! normal & metric terms
                stressDivergenceUCell = stressDivergenceUCell - &
                     stress11(iStressEdge,iCell) * basisIntegralsUNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) - &
                     stress12(iStressEdge,iCell) * basisIntegralsVNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) - &
                     2.0_RKIND * stress12(iStressEdge,iCell) * basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

                stressDivergenceVCell = stressDivergenceVCell - &
                     stress22(iStressEdge,iCell) * basisIntegralsVNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) - &
                     stress12(iStressEdge,iCell) * basisIntegralsUNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) + &
                     stress11(iStressEdge,iCell) * basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge) - &
                     stress22(iStressEdge,iCell) * basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge) 

                stressDivergenceUCellNoMetric = stressDivergenceUCellNoMetric - &
                     stress11(iStressEdge,iCell) * basisIntegralsUNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) - &
                     stress12(iStressEdge,iCell) * basisIntegralsVNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) 

                stressDivergenceVCellNoMetric = stressDivergenceVCellNoMetric - &
                     stress22(iStressEdge,iCell) * basisIntegralsVNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) - &
                     stress12(iStressEdge,iCell) * basisIntegralsUNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge)

                stressDivergenceUCellMetric = stressDivergenceUCellMetric - &
                     stress12(iStressEdge,iCell) * basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

                stressDivergenceVCellMetric = stressDivergenceVCellMetric + &
                     stress11(iStressEdge,iCell) * basisIntegralsMetricNew(iEdge,iStressEdge,iVelocityEdge,iCellOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

             enddo ! iStressEdge

             stressDivergenceUEdge = stressDivergenceUEdge + stressDivergenceUCell
             stressDivergenceVEdge = stressDivergenceVEdge + stressDivergenceVCell

             stressDivergenceUEdgeNoMetric = stressDivergenceUEdgeNoMetric + stressDivergenceUCellNoMetric
             stressDivergenceVEdgeNoMetric = stressDivergenceVEdgeNoMetric + stressDivergenceVCellNoMetric

             stressDivergenceUEdgeMetric = stressDivergenceUEdgeMetric + stressDivergenceUCellMetric
             stressDivergenceVEdgeMetric = stressDivergenceVEdgeMetric + stressDivergenceVCellMetric

          enddo ! iCellOnEdge

          stressDivergenceUCGrid(iEdge) = stressDivergenceUCGrid(iEdge) + stressDivergenceUEdge / variationalDenominatorCGrid(iEdge)
          stressDivergenceVCGrid(iEdge) = stressDivergenceVCGrid(iEdge) + stressDivergenceVEdge / variationalDenominatorCGrid(iEdge)
          stressDivergenceUCGridNoMetric(iEdge) = stressDivergenceUCGridNoMetric(iEdge) + stressDivergenceUEdgeNoMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceVCGridNoMetric(iEdge) = stressDivergenceVCGridNoMetric(iEdge) + stressDivergenceVEdgeNoMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceUCGridMetric(iEdge) = stressDivergenceUCGridMetric(iEdge) + stressDivergenceUEdgeMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceVCGridMetric(iEdge) = stressDivergenceVCGridMetric(iEdge) + stressDivergenceVEdgeMetric / variationalDenominatorCGrid(iEdge)

          stressDivergenceUEdge = 0.0_RKIND
          stressDivergenceVEdge = 0.0_RKIND
          stressDivergenceUEdgeNoMetric = 0.0_RKIND
          stressDivergenceVEdgeNoMetric = 0.0_RKIND
          stressDivergenceUEdgeMetric = 0.0_RKIND
          stressDivergenceVEdgeMetric = 0.0_RKIND

          ! loop over surrounding triangles
          do iVertexOnEdge = 1, 2

             ! get the vertex number of this triangle
             iVertex = verticesOnEdge(iVertexOnEdge, iEdge)

             ! get the local index of iEdge on iTriangleOnEdge
             iVelocityEdge = triangleEdgesAtEdge(iVertexOnEdge,iEdge)

             stressDivergenceUCell = 0.0_RKIND
             stressDivergenceVCell = 0.0_RKIND
             stressDivergenceUCellNoMetric = 0.0_RKIND
             stressDivergenceVCellNoMetric = 0.0_RKIND
             stressDivergenceUCellMetric = 0.0_RKIND
             stressDivergenceVCellMetric = 0.0_RKIND

             ! loop over the edges of the surrounding triangle
             do iStressEdge = 1, vertexDegree

                jEdge = edgesOnVertex(iStressEdge,iVertex)
                ! normal & metric terms
                stressDivergenceUCell = stressDivergenceUCell - & 
                     stress11Tri(iStressEdge,iVertex) * basisIntegralsUTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) - &
                     stress12Tri(iStressEdge,iVertex) * basisIntegralsVTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) - &
                     2.0_RKIND * stress12Tri(iStressEdge,iVertex) * basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

                stressDivergenceVCell = stressDivergenceVCell - &
                     stress22Tri(iStressEdge,iVertex) * basisIntegralsVTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) - &
                     stress12Tri(iStressEdge,iVertex) * basisIntegralsUTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) + &
                     stress11Tri(iStressEdge,iVertex) * basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge) - &
                     stress22Tri(iStressEdge,iVertex) * basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

                stressDivergenceUCellNoMetric = stressDivergenceUCellNoMetric - &
                     stress11Tri(iStressEdge,iVertex) * basisIntegralsUTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) - &
                     stress12Tri(iStressEdge,iVertex) * basisIntegralsVTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge)

                stressDivergenceVCellNoMetric = stressDivergenceVCellNoMetric - &
                     stress22Tri(iStressEdge,iVertex) * basisIntegralsVTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) - &
                     stress12Tri(iStressEdge,iVertex) * basisIntegralsUTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge)

                stressDivergenceUCellMetric = stressDivergenceUCellMetric - &
                     stress12Tri(iStressEdge,iVertex) * basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

                stressDivergenceVCellMetric = stressDivergenceVCellMetric + &
                     stress11Tri(iStressEdge,iVertex) * basisIntegralsMetricTriNew(iEdge,iStressEdge,iVelocityEdge,iVertexOnEdge) * &
                     tanLatEdgeRotatedOverRadius(iEdge)

             enddo ! iStressEdge

             stressDivergenceUEdge = stressDivergenceUEdge + stressDivergenceUCell
             stressDivergenceVEdge = stressDivergenceVEdge + stressDivergenceVCell

             stressDivergenceUEdgeNoMetric = stressDivergenceUEdgeNoMetric + stressDivergenceUCellNoMetric
             stressDivergenceVEdgeNoMetric = stressDivergenceVEdgeNoMetric + stressDivergenceVCellNoMetric

             stressDivergenceUEdgeMetric = stressDivergenceUEdgeMetric + stressDivergenceUCellMetric
             stressDivergenceVEdgeMetric = stressDivergenceVEdgeMetric + stressDivergenceVCellMetric

          enddo ! iTriangleOnEdge

          stressDivergenceUCGrid(iEdge) = stressDivergenceUCGrid(iEdge) + stressDivergenceUEdge / variationalDenominatorCGrid(iEdge) 
          stressDivergenceVCGrid(iEdge) = stressDivergenceVCGrid(iEdge) + stressDivergenceVEdge / variationalDenominatorCGrid(iEdge) 
          stressDivergenceUCGridNoMetric(iEdge) = stressDivergenceUCGridNoMetric(iEdge) + stressDivergenceUEdgeNoMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceVCGridNoMetric(iEdge) = stressDivergenceVCGridNoMetric(iEdge) + stressDivergenceVEdgeNoMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceUCGridMetric(iEdge) = stressDivergenceUCGridMetric(iEdge) + stressDivergenceUEdgeMetric / variationalDenominatorCGrid(iEdge)
          stressDivergenceVCGridMetric(iEdge) = stressDivergenceVCGridMetric(iEdge) + stressDivergenceVEdgeMetric / variationalDenominatorCGrid(iEdge)

       endif ! solveVelocityCGrid

    enddo ! iEdge



    ! -BEGIN: to be removed
    ! INTERPOLATION OF DIV OF STRESS FROM EDGES BACK TO VERTICES
    ! - note: to use this need to do the necessary modification to subroutine def and variable declarations

    !do iVertex = 1, nVerticesSolve
    !   stressDivergenceU(iVertex) = 0.0_RKIND
    !   stressDivergenceV(iVertex) = 0.0_RKIND
    !   sumOfDiamonds = 0.0_RKIND
    !   computeVertex = .true.
    !   do iEdgeOnVertex = 1, vertexDegree
    !      sumOfDiamonds = sumOfDiamonds + variationalDenominatorCGrid(edgesOnVertex(iEdgeOnVertex,iVertex))
    !      if (solveVelocityCGrid(edgesOnVertex(iEdgeOnVertex,iVertex)) == 0) then
    !         computeVertex = .false.
    !      end if
    !   end do
    !   if (computeVertex) then
    !      do iEdgeOnVertex = 1, vertexDegree
    !         iEdge = edgesOnVertex(iEdgeOnVertex,iVertex)
    !         coeff = variationalDenominatorCGrid(iEdge) / sumOfDiamonds
    !         stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + coeff * stressDivergenceUCGrid(iEdge)
    !         stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + coeff * stressDivergenceVCGrid(iEdge)
    !      end do
    !   end if
    !end do

    ! -END: to be removed

  end subroutine seaice_stress_divergence_variational_c_grid!}}}

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_divergence_shear_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date July 9th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_final_divergence_shear_variational(block)

    use seaice_velocity_solver_constitutive_relation, only: &
         eccentricitySquared

    use seaice_mesh_pool, only: &
         nCells, &
         solveStress, &
         nEdgesOnCell

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocityVariationalPool, &
         velocitySolverPool, &
         ridgingPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), pointer :: &
         divergence, &
         shear, &
         ridgeConvergence, &
         ridgeShear

    real(kind=RKIND), dimension(:), allocatable :: &
         DeltaAverage

    real(kind=RKIND) :: &
         strainDivergenceSum, &
         strainTensionSum, &
         strainShearingSum, &
         strainDivergence, &
         strainTension, &
         strainShearing, &
         Delta

    logical, pointer :: &
         config_use_column_package

    integer :: &
         iCell, &
         iVertexOnCell

    call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

    call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
    call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
    call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

    call MPAS_pool_get_array(velocitySolverPool, "divergence", divergence)
    call MPAS_pool_get_array(velocitySolverPool, "shear", shear)

    allocate(DeltaAverage(nCells))

    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          strainDivergenceSum = 0.0_RKIND
          strainTensionSum    = 0.0_RKIND
          strainShearingSum   = 0.0_RKIND
          DeltaAverage(iCell) = 0.0_RKIND

          do iVertexOnCell = 1, nEdgesOnCell(iCell)

             strainDivergence = strain11(iVertexOnCell,iCell) + strain22(iVertexOnCell,iCell)
             strainTension    = strain11(iVertexOnCell,iCell) - strain22(iVertexOnCell,iCell)
             strainShearing   = strain12(iVertexOnCell,iCell) * 2.0_RKIND

             Delta = sqrt(strainDivergence**2 + (strainTension**2 + strainShearing**2) / eccentricitySquared)

             strainDivergenceSum = strainDivergenceSum + strainDivergence
             strainTensionSum    = strainTensionSum    + strainTension
             strainShearingSum   = strainShearingSum   + strainShearing
             DeltaAverage(iCell) = DeltaAverage(iCell) + Delta

          enddo ! iVertexOnCell

          divergence(iCell)   = strainDivergenceSum                              / real(nEdgesOnCell(iCell),RKIND)
          shear(iCell)        = sqrt(strainTensionSum**2 + strainShearingSum**2) / real(nEdgesOnCell(iCell),RKIND)
          DeltaAverage(iCell) = DeltaAverage(iCell)                              / real(nEdgesOnCell(iCell),RKIND)

       else

          divergence(iCell)   = 0.0_RKIND
          shear(iCell)        = 0.0_RKIND

       endif

    enddo ! iCell

    ! ridging parameters
    call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       call MPAS_pool_get_subpool(block % structs, "ridging", ridgingPool)

       call MPAS_pool_get_array(ridgingPool, "ridgeConvergence", ridgeConvergence)
       call MPAS_pool_get_array(ridgingPool, "ridgeShear", ridgeShear)

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             ridgeConvergence(iCell) = -min(divergence(iCell),0.0_RKIND)
             ridgeShear(iCell)       = 0.5_RKIND * (DeltaAverage(iCell) - abs(divergence(iCell)))

          else

             ridgeConvergence(iCell) = 0.0_RKIND
             ridgeShear(iCell)       = 0.0_RKIND

          endif

       enddo ! iCell

    endif ! config_use_column_package

    ! units - for comparison to CICE
    divergence = divergence * 100.0_RKIND * 86400.0_RKIND
    shear      = shear      * 100.0_RKIND * 86400.0_RKIND

    ! cleanup
    deallocate(DeltaAverage)

  end subroutine seaice_final_divergence_shear_variational

!-----------------------------------------------------------------------

end module seaice_velocity_solver_variational
