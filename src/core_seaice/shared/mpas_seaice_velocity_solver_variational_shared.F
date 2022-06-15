!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_variational_shared
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_variational_shared

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       seaice_calc_local_coords, &
       seaice_calc_local_coords_projected_b_grid, &
       seaice_calc_local_coords_c_grid, &
       seaice_calc_variational_metric_terms, &
       seaice_wrapped_index

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_calc_local_coords
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_calc_local_coords(&
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
       onASphere)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:
         
    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell   !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         onASphere           !< Input:

    if (onASphere) then
       call calc_local_coords_spherical(&
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
            rotateCartesianGrid)
    else
       call calc_local_coords_planar(&
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
            zCell)
    endif

  end subroutine seaice_calc_local_coords!}}}


  subroutine seaice_calc_local_coords_projected_b_grid(&
       xLocalNewBGrid, &
       yLocalNewBGrid, &
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
       onASphere)!{{{

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xLocalNewBGrid, & !< Output:
         yLocalNewBGrid    !< Output:
         
    integer, intent(in) :: &
         nCells, & !< Input:
         nVertices, & !< Input:
         vertexDegree !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell   !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell, & !< Input:
         cellsOnVertex

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         onASphere           !< Input:

    if (onASphere) then

        call calc_local_coords_spherical_projected_b_grid(& 
             xLocalNewBGrid, &
             yLocalNewBGrid, &
             nVertices, &
             cellsOnVertex, &
             verticesOnCell, &
             nEdgesOnCell, &
             vertexDegree, &
             xCell, &
             yCell, &
             zCell, &
             xVertex, &
             yVertex, &
             zVertex, &
             rotateCartesianGrid)

     else

        call mpas_log_write('the projected B-grid method is defined only on a sphere')
        call abort

     endif
  end subroutine seaice_calc_local_coords_projected_b_grid!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_local_coords_planar
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine calc_local_coords_planar(&
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
       zCell)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell    !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell

    do iCell = 1, nCells
       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          xLocal(iVertexOnCell,iCell) = xVertex(iVertex) - xCell(iCell)
          yLocal(iVertexOnCell,iCell) = yVertex(iVertex) - yCell(iCell)

        end do
    end do

  end subroutine calc_local_coords_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calc_local_coords_spherical
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine calc_local_coords_spherical(&
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
       rotateCartesianGrid)!{{{

    use seaice_mesh, only: &
         seaice_project_3D_vector_onto_local_2D, &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, dimension(:,:), intent(in) :: &
         verticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex, & !< Input:
         xCell, &   !< Input:
         yCell, &   !< Input:
         zCell      !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    real(kind=RKIND), dimension(3) :: &
         normalVector3D

    real(kind=RKIND), dimension(2) :: &
         normalVector2D

    integer :: &
         iCell, &
         iVertex, &
         iVertexOnCell, &
         jVertexOnCell

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated

    do iCell = 1, nCells

       call seaice_grid_rotation_forward(&
            xCellRotated, yCellRotated, zCellRotated, &
            xCell(iCell), yCell(iCell), zCell(iCell), &
            rotateCartesianGrid)

       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          call seaice_grid_rotation_forward(&
               normalVector3D(1), normalVector3D(2), normalVector3D(3), &
               xVertex(iVertex),  yVertex(iVertex),  zVertex(iVertex), &
               rotateCartesianGrid)

          call seaice_project_3D_vector_onto_local_2D(&
               normalVector2D, &
               normalVector3D, &
               xCellRotated, &
               yCellRotated, &
               zCellRotated)

          xLocal(iVertexOnCell,iCell) = normalVector2D(1) 
          yLocal(iVertexOnCell,iCell) = normalVector2D(2)


       enddo ! iVertexOnCell

    enddo ! iCell

  end subroutine calc_local_coords_spherical!}}}

!-----------------------------------------------------------------------

  subroutine seaice_calc_local_coords_c_grid(&
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
       onASphere, &
       dvEdge, &
       interiorVertex)!{{{

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xLocalNew, & !< Output:
         yLocalNew, & !< Output:
         xLocalTriNew, & !< Output:
         yLocalTriNew  !< Output:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, intent(in) :: &
         nEdges, & !< Input:
         vertexDegree

    integer, dimension(:,:), intent(in) :: &
         verticesOnEdge, & !< Input:
         cellsOnEdge, &    !< Input:
         edgesOnCell, & !< Input:
         edgesOnVertex       !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell, & !< Input:
         yCell, & !< Input:
         zCell, & !< Input:
         xVertex, &
         yVertex, &
         zVertex, &
         xEdge, &   !< Input:
         yEdge, &   !< Input:
         zEdge      !< Input:

    integer, dimension(:), intent(in) :: &
         interiorVertex  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         dvEdge  !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         onASphere              !< Input:

    if (onASphere) then
       call calc_local_coords_spherical_c_grid(&
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
            rotateCartesianGrid)
    else
       call calc_local_coords_planar_c_grid(&
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
            dvEdge, &
            interiorVertex )
    endif


  end subroutine seaice_calc_local_coords_c_grid

!-----------------------------------------------------------------------

  subroutine calc_local_coords_planar_c_grid(&
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
       dvEdge, &
       interiorVertex )!{{{

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xLocalNew, & !< Output:
         yLocalNew, & !< Output:
         xLocalTriNew, & !< Output:
         yLocalTriNew    !< Output:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, intent(in) :: &
         nEdges, & !< Input:
         vertexDegree

    integer, dimension(:,:), intent(in) :: &
         verticesOnEdge, & !< Input:
         cellsOnEdge, &    !< Input:
         edgesOnCell, & !< Input:
         edgesOnVertex       !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell, & !< Input:
         yCell, & !< Input:
         zCell, & !< Input:
         xVertex, &
         yVertex, &
         zVertex, &
         xEdge, &   !< Input:
         yEdge, &   !< Input:
         zEdge      !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         dvEdge  !< Input:

    integer, dimension(:), intent(in) :: &
         interiorVertex    !< Input:

    integer :: &
         iCell, &
         iCellOnEdge, &
         iVertex, &
         iVertexOnEdge, &
         iEdgeOnCell, &
         iEdgeOnVertex, &
         iEdge, &
         jEdge

    real(kind=RKIND) :: &
         halfMeshSize

    do iEdge = 1, nEdges

      do iCellOnEdge = 1, 2

         iCell = cellsOnEdge(iCellOnEdge,iEdge)

         do iEdgeOnCell = 1, nEdgesOnCell(iCell)

            jEdge = edgesOnCell(iEdgeOnCell, iCell)

            xLocalNew(iEdge, iEdgeOnCell, iCellOnEdge) = xEdge(jEdge) - xCell(iCell)
            yLocalNew(iEdge, iEdgeOnCell, iCellOnEdge) = yEdge(jEdge) - yCell(iCell)

         end do

      end do

      do iVertexOnEdge = 1, 2

         iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

         do iEdgeOnVertex = 1, vertexDegree

            jEdge = edgesOnVertex(iEdgeOnVertex, iVertex)

            xLocalTriNew(iEdge, iEdgeOnVertex, iVertexOnEdge) = xEdge(jEdge) - xVertex(iVertex)
            yLocalTriNew(iEdge, iEdgeOnVertex, iVertexOnEdge) = yEdge(jEdge) - yVertex(iVertex)

         end do

         if(size(xLocalTriNew(iEdge,:,iVertexOnEdge)) == 4) then ! the mesh is planar with square cells

            halfMeshSize = abs(maxval(dvEdge(:))) * 0.5_RKIND

            if(interiorVertex(iVertex) == 0) then

               xLocalTriNew(iEdge,1,iVertexOnEdge) = 0.0_RKIND
               yLocalTriNew(iEdge,1,iVertexOnEdge) = halfMeshSize
               xLocalTriNew(iEdge,2,iVertexOnEdge) = - halfMeshSize
               yLocalTriNew(iEdge,2,iVertexOnEdge) = 0.0_RKIND
               xLocalTriNew(iEdge,3,iVertexOnEdge) = 0.0_RKIND
               yLocalTriNew(iEdge,3,iVertexOnEdge) = - halfMeshSize
               xLocalTriNew(iEdge,4,iVertexOnEdge) = halfMeshSize
               yLocalTriNew(iEdge,4,iVertexOnEdge) = 0.0_RKIND

            end if

         end if

      end do

    end do

  end subroutine calc_local_coords_planar_c_grid

  subroutine calc_local_coords_spherical_c_grid(&
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
       rotateCartesianGrid)!{{{

    use seaice_mesh, only: &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xLocalNew, & !< Output:
         yLocalNew, & !< Output:
         xLocalTriNew, & !< Output:
         yLocalTriNew    !< Output:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, intent(in) :: &
         nEdges, & !< Input:
         vertexDegree

    integer, dimension(:,:), intent(in) :: &
         verticesOnEdge, & !< Input:
         cellsOnEdge, &    !< Input:
         edgesOnCell, & !< Input:
         edgesOnVertex       !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell, & !< Input:
         yCell, & !< Input:
         zCell, & !< Input:
         xVertex, &
         yVertex, &
         zVertex, &
         xEdge, &   !< Input:
         yEdge, &   !< Input:
         zEdge      !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    integer :: &
         iCell, &
         iCellOnEdge, &
         iVertex, &
         iVertexOnEdge, &
         iEdgeOnCell, &
         iEdgeOnVertex, &
         iEdge, &
         jEdge

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated, &
         xEdgeRotated, &
         yEdgeRotated, &
         zEdgeRotated, &
         xCellNew, yCellNew, zCellNew, &
         xCenterRotated, yCenterRotated, zCenterRotated, latCenterRotated, lonCenterRotated, xCenter, yCenter, zCenter, &
         e1, e2, e3, n1, n2, n3, norm, latCellRotated, lonCellRotated 

   do iEdge = 1, nEdges

      xCenter = xEdge(iEdge)
      yCenter = yEdge(iEdge)
      zCenter = zEdge(iEdge)

      call seaice_grid_rotation_forward(&
        xCenterRotated, yCenterRotated, zCenterRotated, &
        xCenter, yCenter, zCenter, &
        rotateCartesianGrid)

      norm = sqrt(xCenterRotated * xCenterRotated + yCenterRotated * yCenterRotated + zCenterRotated * zCenterRotated)
      latCenterRotated = asin(zCenterRotated / norm)
      lonCenterRotated = atan2(yCenterRotated, xCenterRotated)
      if (lonCenterRotated < 0.0_RKIND) then
         lonCenterRotated = 2.0_RKIND * 3.14159265359_RKIND + lonCenterRotated
      end if

      e1 = - sin(lonCenterRotated)
      e2 = cos(lonCenterRotated)
      e3 = 0.0_RKIND

      n1 = - cos(lonCenterRotated) * sin(latCenterRotated)
      n2 = - sin(latCenterRotated) * sin(lonCenterRotated)
      n3 = cos(latCenterRotated) 

      do iCellOnEdge = 1, 2

         iCell = cellsOnEdge(iCellOnEdge,iEdge)

          call seaice_grid_rotation_forward(&
               xCellRotated, yCellRotated, zCellRotated, &
               xCell(iCell), yCell(iCell), zCell(iCell), &
               rotateCartesianGrid)

          xCellNew = (xCellRotated - xCenterRotated) * e1 + &
                     (yCellRotated - yCenterRotated) * e2 + &
                     (zCellRotated - zCenterRotated) * e3

          yCellNew = (xCellRotated - xCenterRotated) * n1 + &
                     (yCellRotated - yCenterRotated) * n2 + &
                     (zCellRotated - zCenterRotated) * n3 

         do iEdgeOnCell = 1, nEdgesOnCell(iCell)

            jEdge = edgesOnCell(iEdgeOnCell, iCell)

            call seaice_grid_rotation_forward(&
                 xEdgeRotated, yEdgeRotated, zEdgeRotated, &
                 xEdge(jEdge), yEdge(jEdge), zEdge(jEdge), &
                 rotateCartesianGrid)

             xLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) = (xEdgeRotated - xCenterRotated) * e1 + &
                                           (yEdgeRotated - yCenterRotated) * e2 + &
                                           (zEdgeRotated - zCenterRotated) * e3 

             yLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) = (xEdgeRotated - xCenterRotated) * n1 + &
                                           (yEdgeRotated - yCenterRotated) * n2 + &
                                           (zEdgeRotated - zCenterRotated) * n3

             xLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) = xLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) - xCellNew
             yLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) = yLocalNew(iEdge,iEdgeOnCell,iCellOnEdge) - yCellNew

          end do

      end do

      do iVertexOnEdge = 1,2

         iVertex = verticesOnEdge(iVertexOnEdge,iEdge)

         call seaice_grid_rotation_forward(&
               xCellRotated, yCellRotated, zCellRotated, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)

          xCellNew = (xCellRotated - xCenterRotated) * e1 + &
                     (yCellRotated - yCenterRotated) * e2 + &
                     (zCellRotated - zCenterRotated) * e3

          yCellNew = (xCellRotated - xCenterRotated) * n1 + &
                     (yCellRotated - yCenterRotated) * n2 + &
                     (zCellRotated - zCenterRotated) * n3

         do iEdgeOnVertex = 1, vertexDegree

            jEdge = edgesOnVertex(iEdgeOnVertex, iVertex)

            call seaice_grid_rotation_forward(&
                 xEdgeRotated, yEdgeRotated, zEdgeRotated, &
                 xEdge(jEdge), yEdge(jEdge), zEdge(jEdge), &
                 rotateCartesianGrid)

             xLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) = (xEdgeRotated - xCenterRotated) * e1 + &
                                           (yEdgeRotated - yCenterRotated) * e2 + &
                                           (zEdgeRotated - zCenterRotated) * e3

             yLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) = (xEdgeRotated - xCenterRotated) * n1 + &
                                           (yEdgeRotated - yCenterRotated) * n2 + &
                                           (zEdgeRotated - zCenterRotated) * n3

             xLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) = xLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) - xCellNew
             yLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) = yLocalTriNew(iEdge,iEdgeOnVertex,iVertexOnEdge) - yCellNew

          end do

      end do

   end do

  end subroutine calc_local_coords_spherical_c_grid
  
  subroutine calc_local_coords_spherical_projected_b_grid(&
       xLocalNew, &
       yLocalNew, &
       nVertices, &
       cellsOnVertex, &
       verticesOnCell, &
       nEdgesOnCell, &
       vertexDegree, &
       xCell, &
       yCell, &
       zCell, &
       xVertex, &
       yVertex, &
       zVertex, &
       rotateCartesianGrid)!{{{

    use seaice_mesh, only: &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:,:,:), intent(out) :: &
         xLocalNew, & !< Output:
         yLocalNew !< Output:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell !< Input:

    integer, intent(in) :: &
         nVertices, & !< Input:
         vertexDegree

    integer, dimension(:,:), intent(in) :: &
         cellsOnVertex, &    !< Input:
         VerticesOnCell !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         xCell, & !< Input:
         yCell, & !< Input:
         zCell, & !< Input:
         xVertex, &
         yVertex, &
         zVertex

    logical, intent(in) :: &
         rotateCartesianGrid !< Input:

    integer :: &
         iCell, &
         iCellOnVertex, &
         iVertex, &
         iVertexOnCell, &
         jVertex

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated, &
         xVertexRotated, &
         yVertexRotated, &
         zVertexRotated, &
         xCellNew, yCellNew, zCellNew, &
         xCenterRotated, yCenterRotated, zCenterRotated, latCenterRotated, lonCenterRotated, xCenter, yCenter, zCenter, &
         e1, e2, e3, n1, n2, n3, norm, latCellRotated, lonCellRotated 

   do iVertex = 1, nVertices

      xCenter = xVertex(iVertex)
      yCenter = yVertex(iVertex)
      zCenter = zVertex(iVertex)

      call seaice_grid_rotation_forward(&
        xCenterRotated, yCenterRotated, zCenterRotated, &
        xCenter, yCenter, zCenter, &
        rotateCartesianGrid)

      norm = sqrt(xCenterRotated * xCenterRotated + yCenterRotated * yCenterRotated + zCenterRotated * zCenterRotated)
      latCenterRotated = asin(zCenterRotated / norm)
      lonCenterRotated = atan2(yCenterRotated, xCenterRotated)
      if (lonCenterRotated < 0.0_RKIND) then
         lonCenterRotated = 2.0_RKIND * 3.14159265359_RKIND + lonCenterRotated
      end if

      e1 = - sin(lonCenterRotated)
      e2 = cos(lonCenterRotated)
      e3 = 0.0_RKIND

      n1 = - cos(lonCenterRotated) * sin(latCenterRotated)
      n2 = - sin(latCenterRotated) * sin(lonCenterRotated)
      n3 = cos(latCenterRotated) 

      do iCellOnVertex = 1, vertexDegree

         iCell = cellsOnVertex(iCellOnVertex,iVertex)

          call seaice_grid_rotation_forward(&
               xCellRotated, yCellRotated, zCellRotated, &
               xCell(iCell), yCell(iCell), zCell(iCell), &
               rotateCartesianGrid)

          xCellNew = (xCellRotated - xCenterRotated) * e1 + &
                     (yCellRotated - yCenterRotated) * e2 + &
                     (zCellRotated - zCenterRotated) * e3

          yCellNew = (xCellRotated - xCenterRotated) * n1 + &
                     (yCellRotated - yCenterRotated) * n2 + &
                     (zCellRotated - zCenterRotated) * n3 

         do iVertexOnCell = 1, nEdgesOnCell(iCell)  !double check

            jVertex = verticesOnCell(iVertexOnCell, iCell)  !double check

            call seaice_grid_rotation_forward(&
                 xVertexRotated, yVertexRotated, zVertexRotated, &
                 xVertex(jVertex), yVertex(jVertex), zVertex(jVertex), &
                 rotateCartesianGrid)

             xLocalNew(iVertex,iVertexOnCell,iCellOnVertex) = (xVertexRotated - xCenterRotated) * e1 + &
                                           (yVertexRotated - yCenterRotated) * e2 + &
                                           (zVertexRotated - zCenterRotated) * e3 

             yLocalNew(iVertex,iVertexOnCell,iCellOnVertex) = (xVertexRotated - xCenterRotated) * n1 + &
                                           (yVertexRotated - yCenterRotated) * n2 + &
                                           (zVertexRotated - zCenterRotated) * n3

             xLocalNew(iVertex,iVertexOnCell,iCellOnVertex) = xLocalNew(iVertex,iVertexOnCell,iCellOnVertex) - xCellNew
             yLocalNew(iVertex,iVertexOnCell,iCellOnVertex) = yLocalNew(iVertex,iVertexOnCell,iCellOnVertex) - yCellNew

          end do

      end do

   end do

  end subroutine calc_local_coords_spherical_projected_b_grid

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_calc_variational_metric_terms
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 22 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_calc_variational_metric_terms(&
       tanLatVertexRotatedOverRadius, &
       nVertices, &
       xVertex, &
       yVertex, &
       zVertex, &
       sphereRadius, &
       rotateCartesianGrid, &
       includeMetricTerms)

    use seaice_mesh, only: &
         seaice_grid_rotation_forward

    real(kind=RKIND), dimension(:), intent(out) :: &
         tanLatVertexRotatedOverRadius !< Output:

    integer, intent(in) :: &
         nVertices !< Input:

    real(kind=RKIND), dimension(:), pointer :: &
         xVertex, & !< Input:
         yVertex, & !< Input:
         zVertex    !< Input:

    real(kind=RKIND), pointer :: &
         sphereRadius !< Input:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms     !< Input:

    integer :: &
         iVertex

    real(kind=RKIND) :: &
         xVertexRotated, &
         yVertexRotated, &
         zVertexRotated, &
         latVertexRotated

    if (includeMetricTerms) then

       do iVertex = 1, nVertices

          call seaice_grid_rotation_forward(&
               xVertexRotated,   yVertexRotated,   zVertexRotated, &
               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
               rotateCartesianGrid)

          latVertexRotated = asin(zVertexRotated / sphereRadius)

          tanLatVertexRotatedOverRadius(iVertex) = tan(latVertexRotated) / sphereRadius

       enddo ! iVertex

    else

       do iVertex = 1, nVertices

          tanLatVertexRotatedOverRadius(iVertex) = 0.0_RKIND

       enddo ! iVertex

    endif

  end subroutine seaice_calc_variational_metric_terms

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_wrapped_index
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  function seaice_wrapped_index(&
       input, &
       nelements) &
       result(output)!{{{

    integer, intent(in) :: &
         input, &  !< Input:
         nelements !< Input:

    integer :: output

    output = modulo(input - 1, nelements) + 1

  end function seaice_wrapped_index!}}}

!-----------------------------------------------------------------------

end module seaice_velocity_solver_variational_shared
