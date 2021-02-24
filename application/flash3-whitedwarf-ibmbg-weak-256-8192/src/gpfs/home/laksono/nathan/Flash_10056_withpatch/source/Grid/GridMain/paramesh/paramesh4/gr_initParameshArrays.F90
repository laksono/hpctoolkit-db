!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_initParameshArrays
!!
!! NAME
!!
!!  gr_initParameshArrays
!!
!! SYNOPSIS
!!
!!  call gr_initParameshArrays(logical(IN) :: restart,
!!                             integer(IN) :: nprocs,
!!                             integer(IN) :: mype,
!!                             integer(IN) :: xlboundary,
!!                             integer(IN) :: xrboundary,
!!                             integer(IN) :: ylboundary,
!!                             integer(IN) :: yrboundary,
!!                             integer(IN) :: zlboundary,
!!                             integer(IN) :: zrboundary
!!                             )
!!
!! DESCRIPTION
!!
!!  Perform early initialization of some Grid data structures.
!!
!!  This routine prepares the Grid for being filled with
!!  meaningful data.
!!
!! ARGUMENTS
!!
!!   restart -   Is the grid being prepared for initialization with
!!               data from a checkpoint file?
!!   nprocs  -   Number of PEs.
!!   mype    -   PE executing this code.
!!   xlboundary - coordinate of outer domain boundary in lower X direction.
!!   xrboundary - coordinate of outer domain boundary in upper X direction.
!!   ylboundary - coordinate of outer domain boundary in lower Y direction.
!!   yrboundary - coordinate of outer domain boundary in upper Y direction.
!!   zlboundary - coordinate of outer domain boundary in lower Z direction.
!!   zrboundary - coordinate of outer domain boundary in upper Z direction.
!!
!!
!!
!!***
subroutine gr_initParameshArrays(restart,nprocs,     mype, &
                                     &  xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use paramesh_dimensions
   use physicaldata
   use workspace
   use tree
   use paramesh_mpi_interfaces, ONLY : mpi_amr_global_domain_limits,&
                                       mpi_amr_boundary_block_info
   use paramesh_interfaces, ONLY : amr_refine_derefine,amr_guardcell

   implicit none
#include "constants.h"
#include "Flash.h"
   integer,intent(IN) :: nprocs, mype
   logical,intent(IN) :: restart
   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary
   integer :: lnblocks_old 
   integer :: five,six

! Make sure that blocks produced by divide_domain are in strict
! morton order



   if(restart) then
      call mpi_amr_global_domain_limits()
   else
      call amr_refine_derefine()

      call mpi_amr_global_domain_limits
      call amr_reorder_grid()
   end if
   ! x boundaries
   boundary_box(1,2:3,1:2) = -1.e30
   boundary_box(2,2:3,1:2) =  1.e30
   boundary_box(1,1,1) = -1.e30
   boundary_box(2,1,1) = grid_xmin
   boundary_box(1,1,2) = grid_xmax
   boundary_box(2,1,2) = 1.e30
   boundary_index(1) = xlboundary
   boundary_index(2) = xrboundary
   if (xlboundary.eq.PERIODIC) then ! periodic
      boundary_box(1,2:3,1:2) = 0.
      boundary_box(2,2:3,1:2) =  0.
      boundary_box(1,1,1) = 0.
      boundary_box(2,1,1) = 0.
      boundary_box(1,1,2) = 0.
      boundary_box(2,1,2) = 0.
      boundary_index(1) = 0
    end if
    if (xrboundary.eq.PERIODIC) then ! periodic
      boundary_box(1,2:3,1:2) = 0.
      boundary_box(2,2:3,1:2) =  0.
      boundary_box(1,1,1) = 0.
      boundary_box(2,1,1) = 0.
      boundary_box(1,1,2) = 0.
      boundary_box(2,1,2) = 0.
      boundary_index(2) = 0
    end if
! y boundaries
    if(ndim.ge.2) then

      boundary_box(1,1,3:4) = -1.e30
      boundary_box(2,1,3:4) =  1.e30
      boundary_box(1,3,3:4) = -1.e30
      boundary_box(2,3,3:4) =  1.e30
      boundary_box(1,2,3) = -1.e30
      boundary_box(2,2,3) = grid_ymin
      boundary_box(1,2,4) = grid_ymax
      boundary_box(2,2,4) = 1.e30
      boundary_index(3) = ylboundary
      boundary_index(4) = yrboundary
      if (ylboundary.eq.PERIODIC) then ! periodic
         boundary_box(1,1,3:4) = 0.
         boundary_box(2,1,3:4) = 0.
         boundary_box(1,3,3:4) = 0.
         boundary_box(2,3,3:4) = 0.
         boundary_box(1,2,3) = 0.
         boundary_box(2,2,3) = 0.
         boundary_box(1,2,4) = 0.
         boundary_box(2,2,4) = 0.
         boundary_index(3) = 0
      end if
      if (yrboundary.eq.PERIODIC) then ! periodic
         boundary_box(1,1,3:4) = 0.
         boundary_box(2,1,3:4) = 0.
         boundary_box(1,3,3:4) = 0.
         boundary_box(2,3,3:4) = 0.
         boundary_box(1,2,3) = 0.
         boundary_box(2,2,3) = 0.
         boundary_box(1,2,4) = 0.
         boundary_box(2,2,4) = 0.
         boundary_index(4) = 0
       end if
    endif
! z boundaries
    if(ndim.ge.3) then
       five = 1+(4*k3d)
       six  = five + k3d
       boundary_box(1,1:2,five:six) = -1.e30
       boundary_box(2,1:2,five:six) =  1.e30
       boundary_box(1,3,five) = -1.e30
       boundary_box(2,3,five) = grid_zmin
       boundary_box(1,3,six) = grid_zmax
       boundary_box(2,3,six) = 1.e30
       boundary_index(five) = zlboundary
       boundary_index(six) = zrboundary
       if (zlboundary.eq.PERIODIC) then ! periodic
          boundary_box(1,1:2,five:six) = 0.
          boundary_box(2,1:2,five:six) = 0.
          boundary_box(1,3,five) = 0.
          boundary_box(2,3,five) = 0.
          boundary_box(1,3,six) = 0.
          boundary_box(2,3,six) = 0.
          boundary_index(five) = 0
       end if
       if (zrboundary.eq.PERIODIC) then ! periodic
          boundary_box(1,1:2,five:six) = 0.
          boundary_box(2,1:2,five:six) =  0.
          boundary_box(1,3,five) = 0.
          boundary_box(2,3,five) = 0.
          boundary_box(1,3,six) = 0.
          boundary_box(2,3,six) = 0.
          boundary_index(six) = 0.
       end if
    end if
  
    call amr_morton_process()
    if(restart) then
       grid_analysed_mpi=1
       call mpi_amr_boundary_block_info(mype,nprocs)
    end if
  ! reset for quadratic interpolation
  
  interp_mask_unk(:) = 1
  interp_mask_work(:) = 1
  
  return
end subroutine gr_initParameshArrays

