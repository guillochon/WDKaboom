!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter, gr_globalComm, &
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_max, lrefine_min, &
                   lnblocks
  use Grid_interface, ONLY : Grid_fillGuardCells
  use Driver_interface, ONLY: Driver_getSimTime
  use Grid_interface, ONLY : Grid_getMinCellSize
  use Simulation_data

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  real :: ref_cut,deref_cut,ref_filter,t,blockwidth
  integer       :: l,i,iref,jLo,ierr,max_blocks
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask

  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,&
       selectBlockType=ACTIVE_BLKS)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.


  call Driver_getSimTime(t)

  if (t .eq. sim_tInitial) then
      call gr_markInRadius(sim_xCenter, sim_yCenter, sim_zCenter, wd_radius, lrefine_max, 0)                                
  else
      call MPI_ALLREDUCE(lnblocks,max_blocks,1,MPI_INTEGER,MPI_SUM,gr_globalComm,ierr)
      if (max_blocks .gt. sim_maxBlocks) then
          sim_dynRefineMax = sim_dynRefineMax - 1
      endif
      do l = 1,gr_numRefineVars
         iref = gr_refine_var(l)
         ref_cut = gr_refine_cutoff(l)
         deref_cut = gr_derefine_cutoff(l)
         ref_filter = min(int(gr_refine_filter(l)), sim_dynRefineMax)
         call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
      end do

!      call gr_markVarThreshold(DENS_VAR, sim_rhoCutoff, 1, lrefine_max, .false.)
!      call gr_markVarThreshold(DENS_VAR, 1.0e-2*sim_rhoCutoff, -1, lrefine_min, .true.)
#ifdef FLASH_GRID_PARAMESH2
      ! Make sure lrefine_min and lrefine_max are obeyed - KW
      if (gr_numRefineVars .LE. 0) then
         call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
      end if
#endif

      if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

      if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

      if(gr_lrefineMaxRedDoByLogR) &
           call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
           gr_lrefineCenterJ,gr_lrefineCenterK)
  endif
  
  return
end subroutine Grid_markRefineDerefine

