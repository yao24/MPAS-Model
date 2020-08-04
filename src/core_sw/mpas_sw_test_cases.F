! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module sw_test_cases

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_log


   contains


   subroutine setup_sw_test_case(domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Configure grid metadata and model state for the shallow water test case 
   !   specified in the namelist
   !
   ! Output: block - a subset (not necessarily proper) of the model domain to be
   !                 initialized
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: statePool, meshPool

      integer, pointer :: config_test_case      

      call mpas_pool_get_config(domain % configs, 'config_test_case', config_test_case)

      if (config_test_case == 0) then
         call mpas_log_write('Using initial conditions supplied in input file')

      else if (config_test_case == 1) then
         call mpas_log_write('Setting up shallow water test case 1')
         call mpas_log_write(' -- Advection of Cosine Bell over the Pole')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_1(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 2) then
         call mpas_log_write('Setting up shallow water test case 2')
         call mpas_log_write(' -- Setup shallow water test case 2: Global Steady State Nonlinear Zonal Geostrophic Flow')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_2(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 5) then
         call mpas_log_write('Setting up shallow water test case 5')
         call mpas_log_write(' -- Setup shallow water test case 5: Zonal Flow over an Isolated Mountain')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_5(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else if (config_test_case == 6) then
         call mpas_log_write('Setting up shallow water test case 6')
         call mpas_log_write(' -- Rossby-Haurwitz Wave')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

            call sw_test_case_6(meshPool, statePool)

            call mpas_pool_initialize_time_levels(statePool)

            block_ptr => block_ptr % next
         end do

      else
         call mpas_log_write('Only test case 1, 2, 5, and 6 are currently supported.')
         stop
      end if

      if (config_use_local_time_stepping) then
         call mpas_log_write('Setting up local time step regions and halos')

         block_ptr => domain % blocklist
         do while (associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

            call mpas_init_LTSRegions(meshPool)

            block_ptr => block_ptr % next
         end do
      end if

   end subroutine setup_sw_test_case


   subroutine sw_test_case_1(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 1: Advection of Cosine Bell over the Pole
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: h0 = 1000.0
      real (kind=RKIND), parameter :: theta_c = 0.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: alpha = pii/4.0

      integer :: iCell, iEdge, iVtx
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize cosine bell at (theta_c, lambda_c)
      !
      do iCell = 1, nCells
         r = sphere_distance(theta_c, lambda_c, latCell(iCell), lonCell(iCell), a) 
         if (r < a/3.0) then
            h(1,iCell) = (h0 / 2.0) * (1.0 + cos(pii*r*3.0/a))
         else
            h(1,iCell) = h0 / 2.0
         end if
      end do

   end subroutine sw_test_case_1


   subroutine sw_test_case_2(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 2: Global Steady State Nonlinear Zonal 
   !                                  Geostrophic Flow
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 2.0 * pii * a / (12.0 * 86400.0)
      real (kind=RKIND), parameter :: gh0 = 29400.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h


      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0
      

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                 sin(latVertex(iVtx)) * cos(alpha) - &
                                 cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                               )
      end do
      do iEdge = 1,nEdges
         u(1,iEdge) = -1.0 * ( &
                                 psiVertex(verticesOnEdge(2,iEdge)) - &
                                 psiVertex(verticesOnEdge(1,iEdge)) &
                               ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge = 1, nEdges
         fEdge(iEdge) = 2.0 * omega * &
                            ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                              sin(latEdge(iEdge)) * cos(alpha) &
                            )
      end do
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                          sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                             (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                              sin(latCell(iCell)) * cos(alpha) &
                                             )**2.0 &
                                      ) / &
                                      gravity
      end do

   end subroutine sw_test_case_2


   subroutine sw_test_case_5(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 5: Zonal Flow over an Isolated Mountain
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: u0 = 20.
      real (kind=RKIND), parameter :: gh0 = 5960.0*gravity
      real (kind=RKIND), parameter :: hs0 = 2000.
      real (kind=RKIND), parameter :: theta_c = pii/6.0
      real (kind=RKIND), parameter :: lambda_c = 3.0*pii/2.0
      real (kind=RKIND), parameter :: rr = pii/9.0
      real (kind=RKIND), parameter :: alpha = 0.0

      integer :: iCell, iEdge, iVtx, nTracers
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: r, v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: latEdge, lonEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s, fEdge, fVertex
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
      call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)
      call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
      call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)

      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * u0 * ( &
                                       sin(latVertex(iVtx)) * cos(alpha) - &
                                       cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) &
                                     )
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Generate rotated Coriolis field
      !
      do iEdge = 1, nEdges
         fEdge(iEdge) = 2.0 * omega * &
                                        (-cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha) + &
                                          sin(latEdge(iEdge)) * cos(alpha) &
                                        )
      end do
      do iVtx = 1, nVertices
         fVertex(iVtx) = 2.0 * omega * &
                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha) + &
                                           sin(latVertex(iVtx)) * cos(alpha) &
                                         )
      end do

      !
      ! Initialize mountain
      !
      do iCell = 1, nCells
         if (lonCell(iCell) < 0.0) lonCell(iCell) = lonCell(iCell) + 2.0 * pii
         r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
         h_s(iCell) = hs0 * (1.0 - r/rr)
      end do

      !
      ! Initialize tracer fields
      !
      do iCell = 1, nCells
         r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + (latCell(iCell) - theta_c)**2.0))
         tracers(1,1,iCell) = 1.0 - r/rr
      end do
      if (nTracers > 1) then
         do iCell = 1, nCells
            r = sqrt(min(rr**2.0, (lonCell(iCell) - lambda_c)**2.0 + &
                         (latCell(iCell) - theta_c - pii/6.0)**2.0 &
                        ) &
                    )
            tracers(2,1,iCell) = 1.0 - r/rr
         end do
      end if

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gh0 - (a * omega * u0 + 0.5 * u0**2.0) * &
                                         (-cos(lonCell(iCell)) * cos(latCell(iCell)) * sin(alpha) + &
                                           sin(latCell(iCell)) * cos(alpha) &
                                         )**2.0 &
                                      ) / &
                                      gravity
         h(1,iCell) = h(1,iCell) - h_s(iCell)
      end do

   end subroutine sw_test_case_5


   subroutine sw_test_case_6(meshPool, statePool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Setup shallow water test case 6: Rossby-Haurwitz Wave
   !
   ! Reference: Williamson, D.L., et al., "A Standard Test Set for Numerical 
   !            Approximations to the Shallow Water Equations in Spherical 
   !            Geometry" J. of Comp. Phys., 102, pp. 211--224
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (mpas_pool_type), intent(inout) :: meshPool
      type (mpas_pool_type), intent(inout) :: statePool

      real (kind=RKIND), parameter :: h0 = 8000.0
      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      integer :: iCell, iEdge, iVtx, nTracers
      integer, pointer :: nVertices, nEdges, nCells
      real (kind=RKIND) :: v
      real (kind=RKIND), allocatable, dimension(:) :: psiVertex

      integer, dimension(:,:), pointer :: verticesOnEdge

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:), pointer :: latVertex, lonVertex
      real (kind=RKIND), dimension(:), pointer :: areaCell, areaTriangle, h_s
      real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex

      real (kind=RKIND), dimension(:,:), pointer :: u, h
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers



      !
      ! Scale all distances and areas from a unit sphere to one with radius a
      !

      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
      call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
      call mpas_pool_get_array(meshPool, 'h_s', h_s)

      call mpas_pool_get_array(statePool, 'u', u, 1)
      call mpas_pool_get_array(statePool, 'h', h, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
      nTracers = size(tracers, dim=1)


      xCell = xCell * a
      yCell = yCell * a
      zCell = zCell * a
      xVertex = xVertex * a
      yVertex = yVertex * a
      zVertex = zVertex * a
      xEdge = xEdge * a
      yEdge = yEdge * a
      zEdge = zEdge * a
      dvEdge = dvEdge * a
      dcEdge = dcEdge * a
      areaCell = areaCell * a**2.0
      areaTriangle = areaTriangle * a**2.0
      kiteAreasOnVertex = kiteAreasOnVertex * a**2.0

      !
      ! Initialize wind field
      !
      allocate(psiVertex(nVertices))
      do iVtx = 1, nVertices
         psiVertex(iVtx) = -a * a * w * sin(latVertex(iVtx)) + &
                            a *a * K * (cos(latVertex(iVtx))**R) * &
                            sin(latVertex(iVtx)) * cos(R * lonVertex(iVtx))
      end do
      do iEdge = 1, nEdges
         u(1,iEdge) = -1.0 * ( &
                                               psiVertex(verticesOnEdge(2,iEdge)) - &
                                               psiVertex(verticesOnEdge(1,iEdge)) &
                                             ) / dvEdge(iEdge)
      end do
      deallocate(psiVertex)

      !
      ! Initialize height field (actually, fluid thickness field)
      !
      do iCell = 1, nCells
         h(1,iCell) = (gravity * h0 + a*a*aa(latCell(iCell)) + &
                                                      a*a*bb(latCell(iCell)) * cos(R*lonCell(iCell)) + &
                                                      a*a*cc(latCell(iCell)) * cos(2.0*R*lonCell(iCell)) &
                                      ) / gravity
      end do

   end subroutine sw_test_case_6

   subroutine mpas_init_LTSRegions(meshPool)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Set local time stepping (LTS) regions on the mesh
   !
   ! Output: LTSRegions is written
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none   

   type (mpas_pool_type), intent(inout) :: meshPool

   integer :: iCell, iEdge, cell1, cell2  
   integer, pointer :: nCells, nEdges
   real (kind=RKIND), dimension(:), pointer :: latCell
   integer, dimension(:), pointer :: LTSRegion
   integer, dimension(:,:), pointer :: cellsOnEdge

   call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
   call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

   call mpas_pool_get_array(meshPool, 'latCell', latCell)
   call mpas_pool_get_array(meshPool, 'LTSRegion', LTSRegion)
   call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

   do iCell = 1, nCells
      if (latCell(iCell) < 0.) then
         LTSRegion(iCell) = 1 ! fine region
      else 
         LTSRegion(iCell) = 2 ! coarse region
      end if
   end do

   do iEdge = 1, nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)
      if (LTSRegion(cell1) == 1 .and.  LTSRegion(cell2) == 2) then
         LTSRegion(cell1) = 3 ! interface layer 1
         LTSRegion(cell2) = 4 ! interface layer 2
      else if (LTSRegion(cell1) == 2 .and.  LTSRegion(cell2) == 1) then
            LTSRegion(cell1) = 4
            LTSRegion(cell2) = 3
      else if (LTSRegion(cell1) == 1 .and.  LTSRegion(cell2) == 4) then
            LTSRegion(cell1) = 3
      else if (LTSRegion(cell1) == 4 .and.  LTSRegion(cell2) == 1) then
            LTSRegion(cell2) = 3
      else if (LTSRegion(cell1) == 2 .and.  LTSRegion(cell2) == 3) then
            LTSRegion(cell1) = 4
      else if (LTSRegion(cell1) == 3 .and.  LTSRegion(cell2) == 2) then
            LTSRegion(cell2) = 4
      end if
   end do

   ! the flags set there are needed for the third order LTS scheme
   do iEdge = 1, nEdges
      cell1 = cellsOnEdge(1,iEdge)
      cell2 = cellsOnEdge(2,iEdge)
      if (LTSRegion(cell1) == 1 .and.  LTSRegion(cell2) == 3) then
         LTSRegion(cell1) = 5 ! fine cells that share an edge with interface layer 1
      else if (LTSRegion(cell1) == 3 .and.  LTSRegion(cell2) == 1) then
         LTSRegion(cell2) = 5
      else if (LTSRegion(cell1) == 2 .and.  LTSRegion(cell2) == 4) then
         LTSRegion(cell1) = 6 ! coarse cells that share and edge with interface layer 2
      else if (LTSRegion(cell1) == 4 .and.  LTSRegion(cell2) == 2) then
         LTSRegion(cell2) = 6
      end if
   end do

   end subroutine mpas_init_LTSRegions


   real function sphere_distance(lat1, lon1, lat2, lon2, radius)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on a
   !   sphere with given radius.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

      real (kind=RKIND) :: arg1

      arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                   cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
      sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance


   real function aa(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! A, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      aa = 0.5 * w * (2.0 * omega + w) * cos(theta)**2.0 + &
          0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 + 2.0*R**2.0 - R - 2.0 - 2.0*R**2.0 * cos(theta)**(-2.0))

   end function aa

   
   real function bb(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! B, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      bb = (2.0*(omega + w)*K / ((R+1.0)*(R+2.0))) * cos(theta)**R * ((R**2.0 + 2.0*R + 2.0) - ((R+1.0)*cos(theta))**2.0)

   end function bb


   real function cc(theta)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! C, used in height field computation for Rossby-Haurwitz wave
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real (kind=RKIND), parameter :: w = 7.848e-6
      real (kind=RKIND), parameter :: K = 7.848e-6
      real (kind=RKIND), parameter :: R = 4.0

      real (kind=RKIND), intent(in) :: theta

      cc = 0.25 * K**2.0 * cos(theta)**(2.0*R) * ((R+1.0)*cos(theta)**2.0 - R - 2.0)

   end function cc

end module sw_test_cases 
