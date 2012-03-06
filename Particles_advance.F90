!!****if* source/Particles/ParticlesMain/Particles_advance
!!
!! NAME
!!
!!  Particles_advance
!!
!! SYNOPSIS
!!
!!  Particles_advance(real(in) :: dtOld,
!!                    real(in) :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  Calls passive and active versions
!!  
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles (indirectly).
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!***

!===============================================================================

subroutine Particles_advance (dtOld,dtNew)
    
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, useParticles, &
       pt_gcMaskForAdvance, pt_gcMaskSizeForAdvance, pt_myPe, pt_typeInfo

  use pt_interface, ONLY: pt_advancePassive, pt_advanceActive,  pt_updateTypeDS
  use Grid_interface, ONLY : Grid_moveParticles, Grid_fillGuardCells, &
                             Grid_mapMeshToParticles, Grid_sortParticles
  use Driver_data, ONLY : dr_simTime
  use Simulation_data, ONLY : sim_tExplode
  
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Particles.h"
#include "GridParticles.h"
  
  real, INTENT(in)  :: dtOld, dtNew

  integer       :: i,nstep, myPE
  integer       :: totalPassive, totalActive
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  real          :: jumpx,jumpy,jumpz
  integer       ::p_begin,p_end, p_count
  logical       :: allTypesNotDone
  logical       :: regrid=.false.
  integer, dimension(GRPT_NOSLOPE) :: index_list

!!------------------------------------------------------------------------------
  ! Don't do anything if runtime parameter isn't set
  if (dr_simTime .lt. sim_tExplode) return
  if (.not.useParticles ) return

  
  !We need a work around for pure active particle simulations.
  !NONEXISTENT is known to be -1, which does not conflict with 
  !real particle types.
  !Note the undefine at the end of this subroutine.
#ifdef PASSIVE_PART_TYPE
#define PASSIVE_PART_TYPE_LOCAL PASSIVE_PART_TYPE
#else
#define PASSIVE_PART_TYPE_LOCAL NONEXISTENT
#endif


  ! Prepare guardcell data needed for particle interpolation.
  !
  ! Experimentation with passive particles (with the old way of advancing particles)
  ! has shown that at least 2 layers of guardcells need to be filled
  ! with updated data for vel[xyz] and density, in order to get the
  ! same results as for a full guardcell fill, when using native grid interpolation. - KW
  ! With "monotonic" interpolation, even more layers are needed. - KW

  call Grid_fillGuardCells(pt_myPe, CENTER, ALLDIR,&
       maskSize=pt_gcMaskSizeForAdvance,mask=pt_gcMaskForAdvance)

  ! Sort particles there so we only have to move the minimum of them.
  !  Sort by type, and then within each sort by block

#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif

  ! Now update the pt_typeInfo data structure
  call pt_updateTypeDS(particlesPerBlk)
  

  !! Now do actual movement, advance particles in time
  !! Here we assume that all passive particles have identical
  !! integration. The active particles may also chose to have
  !! identical integration, but they also have to option of picking
  !! different integration methods for different types.
  i=1
  allTypesNotDone=.true.
  do while(allTypesNotDone)
     if(pt_typeInfo(PART_TYPE,i)==PASSIVE_PART_TYPE_LOCAL) then

        !! find the first and the last particle of type passive 
        p_count=pt_typeInfo(PART_LOCAL,i)
        p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
        p_end=p_begin + p_count -1

        !! call the passive advance routine. Only one type of
        !! integration is included for passive type. 
        allTypesNotDone = (i<NPART_TYPES)
        call pt_advancePassive(dtOld,dtNew,particles(:,p_begin:p_end),&
             p_count)
     else
        p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
        p_count=sum(pt_typeInfo(PART_LOCAL,i:NPART_TYPES))
        allTypesNotDone=.false.
        p_end=p_begin+p_count-1
        call pt_advanceActive(dtOld,dtNew,particles(:,p_begin:p_end),&
             p_count)
     end if
     i=i+1
  end do

#ifdef DEBUG_PARTICLES
  print*,' ready to move Particles'
#endif

  index_list(GRPT_POSX_IND)=POSX_PART_PROP
  index_list(GRPT_POSY_IND)=POSY_PART_PROP
  index_list(GRPT_POSZ_IND)=POSZ_PART_PROP
  index_list(GRPT_BLK_IND) = BLK_PART_PROP
  index_list(GRPT_PROC_IND)= BLK_PART_PROP
  index_list(GRPT_TAG_IND) = TAG_PART_PROP

  ! Put the particles in the appropriate blocks if they've moved off
  call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
       index_list, GRPT_NOSLOPE,&
       regrid) 

#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
  ! Now update the pt_typeInfo data structure
  call pt_updateTypeDS(particlesPerBlk)

  
#ifdef DEBUG_PARTICLES
  print*,' back from Grid_moveParticles'
#endif

  ! If predictive routines are used, they will need to sort and prepare for the
  !  next time step.  Since sorting is so expensive, we suffer code duplication
  !  and do it in the pt_preparePassive routines.
  ! Many algorithms use the stub routines.

  allTypesNotDone=.true.
  i=1
  do while(allTypesNotDone)
     if(pt_typeInfo(PART_TYPE,i)==PASSIVE_PART_TYPE_LOCAL) then
        p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
        p_count=pt_typeInfo(PART_LOCAL,i)
        p_end=p_begin+p_count-1
        call pt_preparePassive(dtOld,dtNew,particles(:,p_begin:p_end),p_count)
        allTypesNotDone=.false.
     else
        allTypesNotDone=(i<NPART_TYPES)
     end if
     i=i+1
  end do


#undef PASSIVE_PART_TYPE_LOCAL

  return

  !!-----------------------------------------------------------------------
end subroutine Particles_advance


