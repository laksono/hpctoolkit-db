!!****if* source/Grid/GridMain/paramesh/paramesh4/Grid_fillGuardCells
!!
!! NAME
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  Grid_fillGuardCells(integer(IN) :: myPE, 
!!                      integer(IN) :: gridDataStruct,
!!                      integer(IN) :: idir,
!!             optional,integer(IN) :: minLayers,
!!             optional,integer(IN) :: eosMode,
!!             optional,logical(IN) :: doEos,
!!             optional,integer(IN) :: maskSize,
!!             optional,logical(IN) :: mask(maskSize),
!!             optional,logical(IN) :: makeMaskConsistent)
!!
!!  
!! DESCRIPTION 
!!
!!  This routine fills the guardcells of the specified data structure of 
!!  each block present in the mesh. If the adjacent blocks are at the same
!!  level of refinement, then the values are copied from the appropriate
!!  internal cells of the neighboring block. If the blocks are at different
!!  levels of refinement, then in addition to copying, there will be 
!!  restriction or prolongation. 
!!  
!!  For performance reason, users may choose to use masking in filling the 
!!  guardcells, which is a feature supported only with Paramesh 4. They can 
!!  do so by using the optional arguments mask, maskSize and makeMaskConsistent. However, 
!!  users should exercise extreme caution in using masking. Please see the 
!!  NOTES section for potential side effects. If masking is present, and 
!!  both makeMaskConsistent, and doEos are true, the Eos calculation is done only
!!  if one of the variables selected in the mask is affected by Eos calculation.
!!  A local routine gr_makeMaskConsistent handles this processing.
!!
!!  The argument "gridDataStruct" can take on one of many valid 
!!  values to determine a specific grid data structure on which to apply
!!  the guardcell fill operation. The currently available options are listed with
!!  the arguments. Most users will use CENTER as the option,
!!  since applications typically use the cell centered grid data, and they want
!!  guardcells to be filled for all the variables.
!!  More specialized applications, such as the unsplit methods, may want to use
!!  other options. 
!!  The user can also choose to fill guardcells either in a single direction,
!!  or all of them. For most of Flash solvers supplied with the release
!!  guardcell are filled in all directions. 
!!
!!  The optional arguments related to the Eos calls are provided because if
!!  a solver relies on a call to Eos to establish thermodynamic equilibrium
!!  in various variables, then either the solver can call Eos explicitly itself,
!!  or more conveniently, it can instruct the gcfill process to do it. This 
!!  feature is especially useful when mesh quantities are used to derive 
!!  for example tracer particle attributes. The the user opts to enable doEos,
!!  they will not have to determine whether Eos call is needed, the code will
!!  do it for them
!!  
!!
!!
!! ARGUMENTS 
!!  
!!  myPE - local processor number
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   Paramesh has 5 data structures for grid variables, the first
!!                   four include all physical variables defined on the mesh. The 
!!                   fifth one includes a single variable.
!!
!!                   unk                cell centered, 
!!                   facex,facey,facez  face centered along i,j,k 
!!                                      direction respectively
!!                   work               cell centered, single variable.
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               work 
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   FACES              facex,facey, and facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if for example the algorithm only does one directional
!!         sweep at a time then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants defined
!!         in constants.h IAXIS, JAXIS or KAXIS to fill guardcells in only one 
!!         direction.        
!!         All layers of guardcells in the given direction(s) are filled if one
!!         of IAXIS, JAXIS, or KAXIS is specified, or if ALLDIR is specified
!!         and minLayers (see below) is not present.
!!
!!  minLayers - minimum number of guardcell layers requested for all directions.
!!              The caller requests at least this many layers of
!!              guardcells to be filled.  If idir is given as one of IAXIS,
!!              JAXIS, or KAXIS, this applies to any directions NOT selected by
!!              idir. On the other hand, if idir is given as ALLDIR, then
!!              minLayers appliers to all directions equally.
!!              If not specified, the default is 0 if idir is given as IAXIS,
!!              JAXIS, or KAXIS, meaning no guardcells need to be filled in for
!!              the directions perpendicular to what idir selects; the default
!!              is NGUARD, the full number of available guardcells, if idir
!!              is ALLDIR.
!!              Note that the caller can specify, using minLayers, how many
!!              layers it needs filled, but the implementation may do
!!              more and actually fill all or some additional layers.
!!              The caller must therefore not rely on some guardcells
!!              remaining unchanged.
!!
!!   eosMode -  the EOS mode being used by the solver that is calling the 
!!              routine. The valid values are :
!!              MODE_DEFAULT     the default eos mode being used by the grid unit
!!              MODE_DENS_EI     density and energy input, pressure and temp output
!!              MODE_DENS_PRES   density/pressure input, temperature/energy output
!!              MODE_DENS_TEMP   density/temperature input, pressure/energy output
!!
!!  doEos    - a logical variable indicating if the calling routine wants the
!!             gcfill process to also make sure that Eos is applied to achieve
!!             thermodynamically consistent values of all variables.
!! 
!!  maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!!  
!! makeMaskConsistent - If true when mask is applied, it is made sure that for
!!          all the selected variables in the mask, the ones they are dependent
!!          on are true too. It is also determined whether there is a need to 
!!          apply Eos if doEos argument is true.
!!                  
!! EXAMPLE 
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells(myPE, CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Flash.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells(myPE, WORK, ALLDIR)
!!     
!!     This call fills guardcells along all directions. The operation is applied
!!     to the WORK data structure available in paramesh only.
!!
!! SIDE EFFECTS
!!
!!  After this function returns, all parents of leaf blocks will have current and 
!!  valid solution data (at least for the variables determined by the gridDataStruct
!!  and mask dummy arguments). This is because amr_guardcell calls amr_restrict
!!  internally.
!!
!! NOTES
!!
!!  In the default mode, this routine fills guardcells of all the variables
!!  in the specified data structure. However, if masking is used
!!  it is possible to fill guardcells of only 
!!  selected variables. The users must exercise a great deal of caution in 
!!  using masking in the guardcell filling, since masking out some variables may
!!  have unintended consequences.  For example if any conserved variable
!!  is set to true in the mask, density must be true, otherwise the interpolated
!!  values of the the conserved variable at fine-coarse boundaries will be 
!!  wrong. It is highly recommended that the argument makeMaskConsistent be set to 
!!  true when using masking, since that provides error checking for
!!  the variable sets included in the released solvers. However,  
!!  these checks may not be enough for new solvers introduced by the users.
!!
!!  This function, or one of the lower-level functions invoked by it, MUST
!!  be called (with updated solution data) before child blocks may be removed
!!  in the course of derefinement! See side effects, above, for the reason.
!!
!! 
!!***



subroutine Grid_fillGuardCells(myPE, gridDataStruct, idir,&
     minLayers, eosMode, doEos, maskSize, mask, makeMaskConsistent)

#include "Flash.h"

  use Grid_data, ONLY : gr_blkList, gr_justExchangedGC, &
       gr_convertToConsvdInMeshInterp, gr_numProcs, &
       gr_enableMaskedGCFill, gr_eosModeNow

  use Logfile_interface, ONLY : Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use gr_interface, ONLY : gr_setGcFillNLayers
  use tree, ONLY : lnblocks
  use physicaldata, ONLY : gcell_on_cc,gcell_on_fc, no_permanent_guardcells
  use paramesh_interfaces, ONLY : amr_guardcell, amr_restrict
  use paramesh_mpi_interfaces, ONLY: mpi_amr_comm_setup
  use Eos_interface, ONLY : Eos_guardCells
  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(in) :: myPE
  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional, intent(in) :: eosMode
  logical, optional,intent(IN) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:),intent(IN) :: mask
  logical, optional,intent(IN) :: makeMaskConsistent
  
  integer :: iopt,guard, numLeafBlocks,i,offset,gcEosMode
  integer,dimension(MDIM) :: layers
  integer :: nlayers_transverse

  logical :: lcc, lfc, lec, lnc, lguard, lprolong, lflux, ledge, lrestrict, lfulltree
  integer :: ierr

  logical :: needEos

  logical, save :: maskWarningDone = .FALSE.


#ifdef DEBUG_GRID
  logical:: validDataStructure,validMaskSize
  validDataStructure = (gridDataStruct==CENTER).or.&
                       (gridDataStruct==FACES).or.&
                       (gridDataStruct==FACEX).or.&
                       (gridDataStruct==FACEY).or.&
                       (gridDataStruct==FACEZ).or.&
                       (gridDataStruct==WORK).or.&
                       (gridDataStruct==CENTER_FACES)
  if(.not.validDataStructure)then
     call Driver_abortFlash("GCfill: invalid data structure")
  end if

#endif

  if(present(eosMode)) then
     gcEosMode=eosMode
  else
     gcEosMode=gr_eosModeNow
  end if

  !! If masking is not done then Eos should be applied since it is not known
  !! which variables are of interest
  needEos=.true.

  if (.NOT. gr_enableMaskedGCFill) then

     !! If masking is disabled then a warning is issued and all masking related
     !! processing is skipped

     if (.NOT. maskWarningDone) then
        call Logfile_stampMessage(myPE, 'INFO: Grid_fillGuardCells is ignoring masking.')
        if (myPE==MASTER_PE) print*,    'INFO: Grid_fillGuardCells is ignoring masking.'
        maskWarningDone = .TRUE.
     end if

  else
     
     !! if masking is not explicitly disabled then the presence of a mask allows 
     !! masking to proceed

     if(present(mask))then
        if(present(maskSize)) then

           !! If both mask and masksize are present, apply the mask
           call gr_setMasks(gridDataStruct,maskSize,mask)
           if(present(makeMaskConsistent))then
              if(makeMaskConsistent) then
                 !! if the caller routine is asking for a consistency check
                 !! then mask may be modified, and also determine if eos needs
                 !! to be applied based upon the mask consistency
                 call gr_makeMaskConsistent(gridDataStruct,gcEosMode,needEos)
              end if
           end if
        else  !! if mask is present without the maskSize, abort
           call Driver_abortFlash("gcfill :: maskSize must be present with mask")
        end if
     end if
  end if
  
  !! if guardcells were just exchanged, don't do it again!
  !!if (gr_justExchangedGC) return
  call Timers_start("guardcell Barrier")
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call Timers_stop("guardcell Barrier")
  call Timers_start("guardcell internal")
  !! appropriately mask the data structures to ensure that only the correct data
  !! structure is filled.
  if((gridDataStruct/=CENTER_FACES).and.(gridDataStruct/=CENTER))gcell_on_cc = .false.
  if(gridDataStruct==CENTER)gcell_on_fc=.false.
  if(gridDataStruct==FACEX)gcell_on_fc(JAXIS:KAXIS,:)=.false.
  if(gridDataStruct==FACEY)then
     gcell_on_fc(IAXIS,:)=.false.
     gcell_on_fc(KAXIS,:)=.false.
  end if
  if(gridDataStruct==FACEZ)gcell_on_fc(IAXIS:JAXIS,:)=.false.

  
  guard = NGUARD
  
  iopt = 1
  if(gridDataStruct==WORK)iopt = 2
  
  ! If called for WORK, at this point
  ! - mask will be completely ignored (including its presence or absence)
  ! - the gcell_on_{cc,fc} arrays will not be modified. - KW
  
  
  if (no_permanent_guardcells) then

     ! We don't really need this according to K.O., since it is included in amr_restrict
     ! as long as we don't use the mode in which all tree levels are always propagated in
     ! time and kept always valid. - KW
     call amr_1blk_copy_soln(-1)
     
     iopt = 1
     call amr_restrict(mype, iopt, 1)
     call gr_commSetUp(gridDataStruct)
  else
     if((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
        call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
        call gr_primitiveToConserve(gr_blkList,numLeafBlocks)
     end if
     
     !----------------------------------------------------------------
     ! Figure out nlayers arguments to amr_guardcell based on our arguments
     call gr_setGcFillNLayers(layers, idir, guard, minLayers)
     

#ifdef DEBUG_GRID
     print*, 'amr_guardcell(PE', iopt, guard,layers')'
#endif
     call amr_guardcell(myPE, iopt, guard,layers(IAXIS),layers(JAXIS),layers(KAXIS))  
     !---------------------------------------------------------------------------------
     
     if ((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
        call gr_conserveToPrimitive(gr_blkList,numLeafBlocks, .TRUE.)
        
        if (gr_convertToConsvdInMeshInterp) then
           call gr_sanitizeDataAfterInterp(gr_blkList, numLeafBlocks, 'after gc filling', layers)
        end if
     end if
     
     gr_justExchangedGC = .true.
     
  end if

  if (iopt == 1) then
     ! Reset the global gcell_on_cc and gcell_on_fc to their defaults.
     ! This seems to be needed, otherwise the behavior of PARAMESH routines
     ! called elsewhere is changed. - KW
     gcell_on_cc = .true.
     gcell_on_fc = .true.
  end if

  if(present(doEos)) then
     if(doEos.and.needEos) then
        call Timers_start("eos gc")
        if ((gridDataStruct.NE.CENTER).AND.(gridDataStruct.NE.CENTER_FACES)) &
             call Grid_getListOfBlocks(ACTIVE_BLKS, gr_blkList, numLeafBlocks)
        do i = 1,numLeafBlocks
           call Eos_guardCells(gcEosMode,gr_blkList(i),corners=.true.,layers=layers)
        end do
        call Timers_stop("eos gc")
     end if
  end if

  call Timers_stop("guardcell internal")

end subroutine Grid_fillGuardCells
