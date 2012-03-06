!!****if* source/Simulation/SimulationMain/Sedov/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       integer :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sedov spherical
!!  explosion problem.
!!
!!  References:  Sedov, L. I., 1959, Similarity and Dimensional Methods
!!                 in Mechanics (New York:  Academic)
!!
!!               Landau, L. D., & Lifshitz, E. M., 1987, Fluid Mechanics,
!!                 2d ed. (Oxford:  Pergamon)
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  myPE -          current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_expEnergy      Explosion energy (distributed over 2^dimen central zones)
!!  sim_rInit          Radial position of inner edge of grid (for 1D )
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones` in cells for applying 1d profile
!!
!!
!!***

subroutine Simulation_initBlock (blockId, myPE)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Multispecies_interface, ONLY:  Multispecies_getSumFrac, Multispecies_getSumInv, Multispecies_getAvg
  use Eos_interface, ONLY: Eos
  
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
  
  integer,intent(IN) ::  blockId
  integer,intent(IN) ::  myPE
  
  integer  ::  i, j, k, n, jLo, jHi, tLo, tHi
  integer  ::  ii, jj, kk, put, cnt, nrtemp
  double precision     ::  xDist, yDist, zDist
  double precision     ::  sumRho, sumE, sumVx, sumVy, phi
  double precision     ::  vel, diagonal
  double precision     ::  xx, dxx, yy, dyy, zz, dzz, frac, tfrac
  double precision     ::  vx, vy, vz, p, rho, ei, t, mp, kb, G, vtot
  double precision     ::  dist, th, cdist, cth, cphi, bg_omega, mcell, tcell, crho, comega, cvx, cvy, cei
  double precision     ::  lrho, lomega
  logical  ::  validGeom
  integer  ::  istat

  double precision,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  double precision, dimension(:,:,:,:),pointer :: solnData
  double precision, save    :: xn(SPECIES_BEGIN:SPECIES_END)
  double precision, dimension(EOS_NUM) :: eosData

  logical :: gcell = .true., is_core, is_bg

  double precision  mtot,mu
  integer mode
  !
  !  Construct the radial samples needed for the initialization.
  !

  call PhysicalConstants_get("proton mass", mp)
  call PhysicalConstants_get("Boltzmann", kb)
  call PhysicalConstants_get("Newton", G)
  
  ! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     ! Find a double precision difference between z's if problem is >= 3D
     if (NDIM > 2) then
        if (k .eq. 1) then
           dzz = zCoord(2) - zCoord(1) 
        else
           dzz = zCoord(k) - zCoord(k-1) 
        endif
     ! Otherwise this problem is <= 2D, so dzz is meaningless
     else
       dzz = 0.0
     endif
     zz = zCoord(k)
     
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        ! Find a double precision difference between y's if problem is >= 2D
        if (NDIM > 1) then
           if (j .eq. 1) then
              dyy = yCoord(2) - yCoord(1) 
           else
              dyy = yCoord(j) - yCoord(j-1) 
           endif
        ! Otherwise this problem is <= 1D, so dyy is meaningless
        else
          dyy = 0.0
        endif
        yy = yCoord(j)
        

        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
           if (i .eq. 1) then
              dxx = xCoord(2) - xCoord(1) 
           else
              dxx = xCoord(i) - xCoord(i-1) 
           endif

           xx = xCoord(i) - sim_xCenter
           yy = yCoord(j) - sim_yCenter
           zz = zCoord(k) - sim_zCenter
           cdist = sqrt(xx**2 + yy**2 + zz**2)
           cth = atan2(abs(zz), sqrt(xx**2.d0 + yy**2.d0))
           cphi = atan2(yy, xx)
           call sim_find (theta(1:ns), ns, cth, tLo, tHi, tfrac)
           nrtemp = max(ipos(tLo), ipos(tHi))
           call sim_find (radius(1:nrtemp), nrtemp, cdist, jLo, jHi, frac)

           cnt = 0

           is_bg = .false.
           if (jLo .eq. jHi) then
               is_bg = .true.
               bg_omega = (omega(ipos(tLo), tLo) + tfrac*(omega(ipos(tHi), tHi) - omega(ipos(tLo), tLo)))*(cdist/radius(jLo))**(-3.0d0/2.0d0)
           else
               mcell = m(jLo) + frac*(m(jHi) - m(jLo))
               tcell = (temper(jLo, tLo) + frac*(temper(jHi, tLo) - temper(jLo, tLo)) + &
                   tfrac*(temper(jLo, tHi) - temper(jLo, tLo)) + &
                   frac*tfrac*(temper(jLo, tLo) - temper(jHi, tLo) - temper(jLo, tHi) + temper(jHi, tHi)))
               crho = (rhop(jLo, tLo) + frac*(rhop(jHi, tLo) - rhop(jLo, tLo)) + &
                   tfrac*(rhop(jLo, tHi) - rhop(jLo, tLo)) + &
                   frac*tfrac*(rhop(jLo, tLo) - rhop(jHi, tLo) - rhop(jLo, tHi) + rhop(jHi, tHi)))
               cei = ((eint(jLo, tLo) + frac*(eint(jHi, tLo) - eint(jLo, tLo)) + &
                   tfrac*(eint(jLo, tHi) - eint(jLo, tLo)) + &
                   frac*tfrac*(eint(jLo, tLo) - eint(jHi, tLo) - eint(jLo, tHi) + eint(jHi, tHi))))
               comega = (omega(jLo, tLo) + frac*(omega(jHi, tLo) - omega(jLo, tLo)) + &
                   tfrac*(omega(jLo, tHi) - omega(jLo, tLo)) + &
                   frac*tfrac*(omega(jLo, tLo) - omega(jHi, tLo) - omega(jLo, tHi) + omega(jHi, tHi)))
               cvx = -lrho*comega*cdist*cos(cth)*sin(cphi)
               cvy =  lrho*comega*cdist*cos(cth)*cos(cphi)

               if (mcell .lt. sim_accMass) then
                   is_core = .true.
               else
                   is_core = .false.
               endif
               
               sumRho = 0.d0
               sumE   = 0.d0
               sumVx  = 0.d0
               sumVy  = 0.d0
               
               !
               !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
               !       appropriate quantities along the 1d profile for that subzone.  
               !
               !       Have the final values for the zone be equal to the average of
               !       the subzone values.
               ! 

               do kk = 0, (sim_nSubZones-1)*K3D
                  zz    = zCoord(k) + (kk*sim_inSubzm1-.5)*dzz 
                  zDist = (zz - sim_zCenter) * K3D
                  
                  do jj = 0, (sim_nSubZones-1)*K2D
                     yy    = yCoord(j) + (jj*sim_inSubzm1-.5)*dyy
                     yDist = (yy - sim_yCenter) * K2D
                     
                     do ii = 0, (sim_nSubZones-1)
                        xx    = xCoord(i) + (ii*sim_inSubzm1-.5)*dxx
                        xDist = xx - sim_xCenter
                        
                        dist = sqrt( xDist**2 + yDist**2 + zDist**2 )
                        th = atan2(abs(zDist), sqrt(xDist**2.d0 + yDist**2.d0))
                        phi = atan2(yDist, xDist)
                        call sim_find (theta(1:ns), ns, th, tLo, tHi, tfrac)
                        nrtemp = min(ipos(tLo), ipos(tHi))
                        call sim_find (radius(1:nrtemp), nrtemp, dist, jLo, jHi, frac)

                        if (jLo .eq. jHi) cycle

                        if (is_core .and. m(jHi) .lt. sim_accMass) then
                            lrho = (rhop(jLo, tLo) + frac*(rhop(jHi, tLo) - rhop(jLo, tLo)) + &
                                tfrac*(rhop(jLo, tHi) - rhop(jLo, tLo)) + &
                                frac*tfrac*(rhop(jLo, tLo) - rhop(jHi, tLo) - rhop(jLo, tHi) + rhop(jHi, tHi)))
                            sumRho = sumRho + lrho
                            sumE = sumE + lrho* &
                                ((eint(jLo, tLo) + frac*(eint(jHi, tLo) - eint(jLo, tLo)) + &
                                tfrac*(eint(jLo, tHi) - eint(jLo, tLo)) + &
                                frac*tfrac*(eint(jLo, tLo) - eint(jHi, tLo) - eint(jLo, tHi) + eint(jHi, tHi))))
                            cnt = cnt + 1
                        endif
                        if (.not. is_core .and. m(jLo) .ge. sim_accMass) then
                            lrho = (rhop(jLo, tLo) + frac*(rhop(jHi, tLo) - rhop(jLo, tLo)) + &
                                tfrac*(rhop(jLo, tHi) - rhop(jLo, tLo)) + &
                                frac*tfrac*(rhop(jLo, tLo) - rhop(jHi, tLo) - rhop(jLo, tHi) + rhop(jHi, tHi)))
                            sumRho = sumRho + lrho
                            sumE = sumE + lrho* &
                                ((eint(jLo, tLo) + frac*(eint(jHi, tLo) - eint(jLo, tLo)) + &
                                tfrac*(eint(jLo, tHi) - eint(jLo, tLo)) + &
                                frac*tfrac*(eint(jLo, tLo) - eint(jHi, tLo) - eint(jLo, tHi) + eint(jHi, tHi))))
                            lomega = (omega(jLo, tLo) + frac*(omega(jHi, tLo) - omega(jLo, tLo)) + &
                                tfrac*(omega(jLo, tHi) - omega(jLo, tLo)) + &
                                frac*tfrac*(omega(jLo, tLo) - omega(jHi, tLo) - omega(jLo, tHi) + omega(jHi, tHi)))
                            sumVx = sumVx - lrho*lomega*dist*cos(th)*sin(phi)
                            sumVy = sumVy + lrho*lomega*dist*cos(th)*cos(phi)
                            cnt = cnt + 1
                        endif
                     enddo
                  enddo
               enddo
           endif
           
           if (.not. is_bg .and. cnt .gt. 0.0) then
              t = tcell
              if (cnt .gt. 0) then
                  rho = max (sumRho / cnt, sim_rhoAmbient)
                  if (is_core) then
                      xn = core_xn
                      eosData(EOS_DENS) = rho
                      eosData(EOS_TEMP) = tcell
                      call Eos(MODE_DENS_TEMP,1,eosData,core_xn)
                      ei = eosData(EOS_EINT)
                      vx = 0.d0
                      vy = 0.d0
                  else 
                      xn = torus_xn
                      vx = sumVx / sumRho
                      vy = sumVy / sumRho 
                      ei = sumE / sumRho
                  endif
              else
                  rho = crho
                  if (is_core) then
                      xn = core_xn
                      eosData(EOS_DENS) = rho
                      eosData(EOS_TEMP) = tcell
                      call Eos(MODE_DENS_TEMP,1,eosData,core_xn)
                      ei = eosData(EOS_EINT)
                      vx = 0.d0
                      vy = 0.d0
                  else 
                      xn = torus_xn
                      vx = cvx
                      vy = cvy
                      ei = cei
                  endif
              endif
           else
              xn = torus_xn
              rho = sim_rhoAmbient
              t = sim_tAmbient
              eosData(EOS_DENS) = rho
              eosData(EOS_TEMP) = t
              call Eos(MODE_DENS_TEMP,1,eosData,torus_xn)
              ei = eosData(EOS_EINT)
              vx = 0.0
              vy = 0.0
              !vx = -bg_omega*cdist*sin(cphi)*cos(cth)
              !vy =  bg_omega*cdist*cos(cphi)*cos(cth)
           endif

           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
           solnData(DENS_VAR,i,j,k)=rho
           solnData(TEMP_VAR,i,j,k)=t
           solnData(EINT_VAR,i,j,k)=ei
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=0.0
           do put=1,NSPECIES
              solnData(SPECIES_BEGIN+put-1,i,j,k)=xn(SPECIES_BEGIN+put-1)
           enddo
#else
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, t)
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, ei)    
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, 0.0)
           do put=1,NSPECIES
              call Grid_putPointData(blockID,CENTER,SPECIES_BEGIN+put-1,&
                   EXTERIOR,axis,xn(SPECIES_BEGIN+put-1))
           enddo
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock


!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(nn) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, nn, x0, lo, hi, f)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: nn
  integer, intent(OUT):: lo, hi
  double precision, intent(IN)    :: x(nn), x0
  double precision, intent(OUT)   :: f

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     lo = 0

  elseif (x0 .gt. x(nn)) then

     lo = nn

  else

     il = 1
     ir = nn
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   lo = il

  endif
  if (lo .eq. 0) then
     lo = 1
     hi = 1
     f = 0.d0
  else if (lo .eq. nn) then
     hi = nn
     f = 0.d0
  else
     hi = lo + 1
     f = (x0 - x(lo)) / (x(hi)-x(lo))
  endif

  return
end subroutine sim_find
