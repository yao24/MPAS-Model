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
       onASphere, &
       dvEdge, &
       interiorVertex)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell  !< Input:

    integer, dimension(:), optional, intent(in) :: &
         interiorVertex  !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         dvEdge  !< Input:

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
         onASphere              !< Input:

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
            zCell, &
            dvEdge, &
            interiorVertex)
    endif

  end subroutine seaice_calc_local_coords!}}}

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
       zCell, &
       dvEdge, &
       interiorVertex )!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell    !< Input:

    integer, dimension(:), optional, intent(in) :: &
         interiorVertex    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         dvEdge  !< Input:

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

    real(kind=RKIND) :: &
         halfMeshSize

    do iCell = 1, nCells
       do iVertexOnCell = 1, nEdgesOnCell(iCell)

          iVertex = verticesOnCell(iVertexOnCell, iCell)

          xLocal(iVertexOnCell,iCell) = xVertex(iVertex) - xCell(iCell)
          yLocal(iVertexOnCell,iCell) = yVertex(iVertex) - yCell(iCell)

        end do
    end do


    if(size(xLocal(:,1)) == 4 .and. present(interiorVertex)) then ! the mesh is planar with square cells

       halfMeshSize = abs(maxval(dvEdge(:))) * 0.5_RKIND    

       do iCell = 1, nCells

          if(interiorVertex(iCell) == 0) then

             !put in some junk so Wachspress can run (but the edges with at
             !least one vertex on the boundary should not be considered because
             !they won't produce a correct value with the trick below)
             xLocal(1,iCell) = 0.0_RKIND
             yLocal(1,iCell) = halfMeshSize
             xLocal(2,iCell) = - halfMeshSize
             yLocal(2,iCell) = 0.0_RKIND
             xLocal(3,iCell) = 0.0_RKIND
             yLocal(3,iCell) = - halfMeshSize
             xLocal(4,iCell) = halfMeshSize
             yLocal(4,iCell) = 0.0_RKIND
  
          end if

       end do

    end if

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
         iVertexOnCell

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated, &
         e1, e2, e3, n1, n2, n3, norm, latCellRotated, lonCellRotated, &
         xVertexRotated, yVertexRotated, zVertexRotated

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

    !do iCell = 1, nCells

       !if (size(xLocal(:,iCell)) == 3 .and. iCell == 1116) then

    !   do iVertexOnCell = 1, nEdgesOnCell(iCell)
    !      print*, 'local coords old=', iVertexOnCell, xLocal(iVertexOnCell,iCell), yLocal(iVertexOnCell,iCell)
    !   end do

    !   call seaice_grid_rotation_forward(&
    !      xCellRotated, yCellRotated, zCellRotated, &
    !      xCell(iCell), yCell(iCell), zCell(iCell), &
    !      rotateCartesianGrid)

    !   norm = sqrt(xCellRotated * xCellRotated + yCellRotated * yCellRotated + zCellRotated * zCellRotated)
    !   latCellRotated = asin(zCellRotated / norm) 
    !   lonCellRotated = atan2(yCellRotated, xCellRotated)
    !   if (lonCellRotated < 0.0_RKIND) then
    !      lonCellRotated = 2.0_RKIND * 3.14159265359_RKIND + lonCellRotated
    !   end if

    !   !components of the eastwards unit vector on tg plane at (xCellRoated, yCellRotated, zCellRotated)
    !   !obtained after a coordinate change on (0,1,0)
    !   e1 = xCellRotated * sin(lonCellRotated) + (1.0_RKIND - yCellRotated) * cos(lonCellRotated)
    !   e2 = xCellRotated * cos(lonCellRotated) * sin(latCellRotated) &
    !        - sin(lonCellRotated) * sin(latCellRotated) * (1.0_RKIND - yCellRotated) &
    !        - zCellRotated * cos(latCellRotated)
    !   e3 = - xCellRotated * cos(lonCellRotated) * cos(latCellRotated) &
    !        + sin(lonCellRotated) * cos(latCellRotated) * (1.0_RKIND - yCellRotated) &
    !        - zCellRotated * sin(latCellRotated)

    !   norm = sqrt(e1*e1 + e2*e2 + e3*e3)

    !   e1 = e1 / norm
    !   e2 = e2 / norm
    !   e3 = e3 / norm

    !   !components of the northwards unit vector on tg plane at (xCellRoated, yCellRotated, zCellRotated)
    !   !obtained after a coordinate change on (0,0,1)
    !   n1 = xCellRotated * sin(lonCellRotated) - yCellRotated * cos(lonCellRotated)
    !   n2 = xCellRotated * cos(lonCellRotated) * sin(latCellRotated) &
    !        + sin(lonCellRotated) * sin(latCellRotated) * yCellRotated &
    !        + (1.0_RKIND - zCellRotated) * cos(latCellRotated)
    !   n3 = - xCellRotated * cos(lonCellRotated) * cos(latCellRotated) &
    !        - sin(lonCellRotated) * cos(latCellRotated) * yCellRotated &
    !        + (1.0_RKIND - zCellRotated) * sin(latCellRotated)

    !   norm = sqrt(n1*n1 + n2*n2 + n3*n3)

    !   n1 = n1 / norm
    !   n2 = n2 / norm
    !   n3 = n3 / norm
         
    !   do iVertexOnCell = 1, nEdgesOnCell(iCell)

    !      iVertex = verticesOnCell(iVertexOnCell, iCell)

    !      call seaice_grid_rotation_forward(&
    !             xVertexRotated, yVertexRotated, zVertexRotated, &
    !             xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
    !             rotateCartesianGrid)

    !      xLocal(iVertexOnCell,iCell) = (xVertexRotated - xCellRotated) * e1 + &
    !                                    (yVertexRotated - yCellRotated) * e2 + &
    !                                    (zVertexRotated - zCellRotated) * e3
    !      yLocal(iVertexOnCell,iCell) = (xVertexRotated - xCellRotated) * n1 + &
    !                                    (yVertexRotated - yCellRotated) * n2 + &
    !                                    (zVertexRotated - zCellRotated) * n3
            
    !   end do

    !   do iVertexOnCell = 1, nEdgesOnCell(iCell)
    !      xLocal(iVertexOnCell,iCell) = -   xLocal(iVertexOnCell,iCell)
    !      print*, 'local coords new=', iVertexOnCell, xLocal(iVertexOnCell,iCell), yLocal(iVertexOnCell,iCell)
    !   end do

       !end if

    !enddo ! iCell

  end subroutine calc_local_coords_spherical!}}}

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
