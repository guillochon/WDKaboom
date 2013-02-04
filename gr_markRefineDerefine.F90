!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_markRefineDerefine
!!
!! NAME
!!  gr_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  gr_markRefineDerefine(integer(IN) :: myPE,
!!                        integer(IN) :: iref,
!!                        real(IN) :: refine_cutoff,
!!                        real(IN) :: derefine_cutoff,
!!                        real(IN) :: refine_filter)
!!  
!!  DESCRIPTION
!!  
!!    Blocks are marked for refining or derefining.
!!    This version uses the second derivative calculations on the specified variable to 
!!    determine if the block needs more resoultion (refine) or less resolution (derefine)
!!    de/refine_cutoff are the thresholds for triggering the corresponding action.
!!    Once the blocks have been marked, the control is passed to Paramesh to update refinement.
!!
!!  ARGUMENTS 
!!
!!    myPE - local processor nnumber
!!
!!    iref - index of the refinement variable in data structure "unk"
!!
!!    refine_cutoff - the threshold value for triggering refinement 
!!
!!    derefine_cutoff - the threshold for triggereing derefinement
!!
!!    refine_filter - makes sure that error calculations to determine refinement
!!                    don't diverge numerically 
!! 
!!  NOTES
!!  
!!    See Grid_markRefineDerefine
!!
!!  SEE ALSO
!!  
!!    Grid_markRefineDerefine
!!
!!***

!!REORDER(5): unk, unk1


subroutine gr_markRefineDerefine(myPE,&
                              iref,refine_cutoff,derefine_cutoff,refine_filter)

  use paramesh_dimensions
  use physicaldata, ONLY : gcell_on_cc, unk, unk1, no_permanent_guardcells
  use tree
  use Grid_data, ONLY: gr_geometry, gr_oneBlock, gr_maxRefine, gr_lrefineMaxRedDoByTime

  use paramesh_interfaces, ONLY: amr_1blk_guardcell

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"  

  integer, intent(IN) :: myPE,iref
  real, intent(IN) :: refine_cutoff, derefine_cutoff, refine_filter
  integer, parameter :: SQNDIM = NDIM*NDIM
  
  real delx,dely,delz
  real dely_f, delz_f
  real delu(mdim,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
  real delua(mdim,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

  real delu2(SQNDIM), delu3(SQNDIM), delu4(SQNDIM)

  real num,denom,error(maxblocks),error_par(maxblocks)
  real J_sq,xx,dd
  real error_max
  real vx,vy,vz,bv,bx,by,bz,bx1,by1,bz1,rhbi,rhi,bk,et
  real px,py,pz,rho,b2,v2,gf,vv,gr,dx,vt,press
  real times,time_exe

  integer lb,i,j,k,nref
  integer ineigh,ineigh_proc,ierr
  integer kstart,kend,jstart,jend,istart,iend
  integer nsend,nrecv
  integer reqr(maxblocks),reqs(maxblocks*nchild)
!
  integer :: ii, jj, kk, ip, jp, kp
  integer :: ndim2
  integer :: statr(MPI_STATUS_SIZE,maxblocks)
  integer :: stats(MPI_STATUS_SIZE,maxblocks*nchild)

  real, pointer :: solnData(:,:,:)
  logical :: gcell_on_cc_backup(NUNK_VARS)
  integer :: idest, iopt, nlayers, icoord
  logical :: lcc, lfc, lec, lnc, l_srl_only, ldiag

  integer :: max_refine

!==============================================================================

  max_refine = nint(refine_filter) !Hack to set maximum refinement level for a given criteria

  if (no_permanent_guardcells) gcell_on_cc_backup = gcell_on_cc

  if (.not. no_permanent_guardcells) then
! A non-directional guardcell fill for CENTER (and also EOS calls for
! all block cells, including guardcells, if any refinement variables
! refine_var_# require this to be current) must have been performed
! when this routine is invoked. Moreover, there must not be any
! intervening calls that would modify the solution data in unk (at
! least, for the variables to be used for refinement criteria).
! Finally, this should be true for both LEAF and PARENT blocks
! (node types 1 and 2).
! Normally the caller of this routine, Grid_markRefineDerefine, takes care
! of all that.
!
! If this routine must be used in a situation where the conditions above
! are not true, the simplest (but probably inefficient) way of adapting
! this code to that situation would be uncommenting the following line:
!!$  call Grid_fillGuardCells(myPE,CENTER_FACES,ALLDIR)

! The following is unused code - it copies only the inner cells to work.
!!$  do lb = 1,lnblocks
!!$     do k = NGUARD*K3D+1,NGUARD*K3D+NZB
!!$        do j = NGUARD*K2D+1,NGUARD*K2D+NYB
!!$           do i = NGUARD+1,NGUARD+NXB
!!$              work(i,j,k,lb,1) = unk(iref,i,j,k,lb)
!!$           end do
!!$        end do
!!$     end do
!!$  end do

! We are using more cell layers, including guardcells, from unk.

     
!     work(:,:,:,:,1)=unk(iref,:,:,:,:)

  end if
  !==============================================================================

  ndim2 = ndim*ndim


#define XCOORD(I) (gr_oneBlock(lb)%firstAxisCoords(CENTER,I))
#define YCOORD(I) (gr_oneBlock(lb)%secondAxisCoords(CENTER,I))

  do lb = 1,lnblocks
      error(lb) = 0.d0

      if (nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
          error(lb) = maxval(unk(iref,:,:,:,lb))
      end if
  end do
     
  
! MARK FOR REFINEMENT OR DEREFINEMENT

! first communicate error of parent to its children
! Children collect messages from parents.

  error_par(1:lnblocks) = 0.
  nrecv = 0
  do lb = 1,lnblocks
     if(parent(1,lb).gt.-1) then
        if (parent(2,lb).ne.myPE) then
           nrecv = nrecv + 1
           call MPI_IRecv(error_par(lb),1, &
                MPI_DOUBLE_PRECISION, &
                parent(2,lb), &
                lb, &
                MPI_COMM_WORLD, &
                reqr(nrecv), &
                ierr)
        else
           error_par(lb) = error(parent(1,lb))
        end if
     end if
  end do
 
  ! parents send error to children

  nsend = 0
  do lb = 1,lnblocks
     do j = 1,nchild
        if(child(1,j,lb).gt.-1) then
           if (child(2,j,lb).ne.myPE) then
              nsend = nsend + 1
              call MPI_ISend(error(lb), &
                   1, &
                   MPI_DOUBLE_PRECISION, &
                   child(2,j,lb), &  ! PE TO SEND TO
                   child(1,j,lb), &  ! THIS IS THE TAG
                   MPI_COMM_WORLD, &
                   reqs(nsend), &
                   ierr)
           end if
        end if
     end do
  end do

  if (nsend.gt.0) then
     call MPI_Waitall (nsend, reqs, stats, ierr)
  end if
  if (nrecv.gt.0) then
     call MPI_Waitall (nrecv, reqr, statr, ierr)
  end if

  do lb = 1,lnblocks

     if (nodetype(lb).eq.1) then
        
        ! test for derefinement
        
        if (.not.refine(lb).and..not.stay(lb) &
             &          .and.error(lb).le.derefine_cutoff &
             &          .and.error_par(lb).le.derefine_cutoff) then
           derefine(lb) = .TRUE.
        else
           derefine(lb) = .FALSE.
        end if
        
        ! test for refinement
        if (error(lb).gt.refine_cutoff) then
           derefine(lb) = .FALSE.
           refine(lb) = .TRUE.
        end if

        if (error(lb).gt.derefine_cutoff.or.error_par(lb).gt.derefine_cutoff)  &
             &           stay(lb) = .TRUE.

        if (lrefine(lb).ge.max_refine) refine(lb) = .FALSE.
        if (lrefine(lb).gt.max_refine) derefine(lb) = .TRUE.

     end if
     
  end do

  !restore to the state when we came in
  if (no_permanent_guardcells) gcell_on_cc = gcell_on_cc_backup
  !=========================================================================
  return
end subroutine gr_markRefineDerefine














