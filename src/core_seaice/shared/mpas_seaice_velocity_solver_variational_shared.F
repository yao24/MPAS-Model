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
       interiorVertex)!{{{

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         xLocal, & !< Output:
         yLocal    !< Output:

    integer, intent(in) :: &
         nCells !< Input:

    integer, dimension(:), intent(in) :: &
         nEdgesOnCell   !< Input:

    integer, dimension(:), optional, intent(in) :: &
         interiorVertex  !< Input:

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

       halfMeshSize = abs(xVertex(verticesOnCell(1,iCell)) - xCell(iCell))

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
         iVertexOnCell, &
         jVertexOnCell

    real(kind=RKIND) :: &
         xCellRotated, &
         yCellRotated, &
         zCellRotated, &
         e1, e2, e3, n1, n2, n3, norm, latCellRotated, lonCellRotated, &
         xVertexRotated, yVertexRotated, zVertexRotated, difference

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

    !   e1 = - sin(lonCellRotated)
    !   e2 = cos(lonCellRotated)
    !   e3 = 0.0_RKIND

    !   !norm = sqrt(e1*e1 + e2*e2 + e3*e3)
    !   !print*, 'norm e=', norm

    !   n1 = - cos(lonCellRotated) * sin(latCellRotated)
    !   n2 = - sin(latCellRotated) * sin(lonCellRotated)
    !   n3 = cos(latCellRotated)

    !   !norm = sqrt(n1*n1 + n2*n2 + n3*n3)
    !   !print*, 'norm n=', norm

    !   !check to se if e and u are perpendicular to normal vector of plane:
    !   !print*, "e:", e1 * xCellRotated + e2 * yCellRotated + e3 * zCellRotated
    !   !print*, "n:", n1 * xCellRotated + n2 * yCellRotated + n3 * zCellRotated
    !   !print*, "e dot n:", e1*n1 + e2*n2 + e3*n3

    !   do iVertexOnCell = 1, nEdgesOnCell(iCell)

    !      iVertex = verticesOnCell(iVertexOnCell, iCell)

    !      call seaice_grid_rotation_forward(&
    !             xVertexRotated, yVertexRotated, zVertexRotated, &
    !             xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
    !             rotateCartesianGrid)

    !      ! ------------BEGIN: REMOVE THIS IN DUE TIME
    !      difference = xLocal(iVertexOnCell,iCell) - ((xVertexRotated - xCellRotated)* e1 + &
    !                                    (yVertexRotated - yCellRotated) * e2 + &
    !                                    (zVertexRotated - zCellRotated) * e3)

    !      if(abs(difference) > 1.0e-10) then
    !        print*, iCell,'is no good for xLocal', ' difference:', difference
    !      end if

    !      !xLocal(iVertexOnCell,iCell) = (xVertexRotated * e1 + &
    !      !                               yVertexRotated * e2 + &
    !      !                               zVertexRotated * e3)
          
    !      difference = yLocal(iVertexOnCell,iCell) - ((xVertexRotated - xCellRotated) * n1 + &
    !                                    (yVertexRotated - yCellRotated) * n2 + &
    !                                    (zVertexRotated - zCellRotated) * n3)
    !      if(abs(difference) > 1.0e-10 ) then
    !        print*, iCell,'is no good for yLocal', 'difference:', difference 
    !      end if

    !      !yLocal(iVertexOnCell,iCell) = (xVertexRotated * n1 + &
    !      !                               yVertexRotated * n2 + &
    !      !                               zVertexRotated * n3)
    !      ! ------------END: REMOVE THIS IN DUE TIME

    !      xLocal(iVertexOnCell,iCell) = (xVertexRotated - xCellRotated) * e1 + &
    !                                    (yVertexRotated - yCellRotated) * e2 + &
    !                                    (zVertexRotated - zCellRotated) * e3
    !      yLocal(iVertexOnCell,iCell) = (xVertexRotated - xCellRotated) * n1 + &
    !                                    (yVertexRotated - yCellRotated) * n2 + &
    !                                    (zVertexRotated - zCellRotated) * n3

    !   end do

       !difference = xLocal(1,iCell) * yLocal(2,iCell) - yLocal(1,iCell) * xLocal(2,iCell)

       !if (difference < 0.0_RKIND) then

          !print*, "HERE BUT SHOULD NOT BE ------------>"

          !iVertexOnCell = verticesOnCell(2,iCell) 
          !jVertexOnCell = verticesOnCell(3,iCell)
          !verticesOnCell(2,iCell) = jVertexOnCell
          !verticesOnCell(3,iCell) = iVertexOnCell

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

       !end if

       !difference = xLocal(1,iCell) * yLocal(2,iCell) - yLocal(1,iCell) * xLocal(2,iCell)
       !if( difference < 0.0_RKIND) then
       !!   print*, 'fix did not work:',  iCell, difference
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
