!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt) 
    use Driver_data, ONLY: dr_simTime, dr_initialSimTime
    !use Flame_interface, ONLY : Flame_step
    use Stir_interface, ONLY : Stir
    use Heat_interface, ONLY : Heat
    use Burn_interface, ONLY : Burn
    use Cool_interface, ONLY : Cool
    use Simulation_data, ONLY: sim_smallX, sim_rhoAmbient, sim_tRelax, wd_radius, sim_accTemp, rhop, &
        sim_tInitial, ipos, sim_accMass, sim_relaxRate, radius, m, exploded, sim_tExplode, core_pos, &
        sim_rotFac, sim_tSpinup, core_xn, torus_xn, sim_detDens, sim_detTemp, sim_detRadius, &
        sim_critDens, sim_critKine, sim_explodeCore, sim_dbleDetTemp, sim_detHeight
    use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
        Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize
    use Eos_interface, ONLY : Eos_wrapped, Eos
    use gr_mpoleData, ONLY: gr_mpoleXcenterOfMass, gr_mpoleYcenterOfMass, gr_mpoleZcenterOfMass, gr_mpoleTotalMass
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, RuntimeParameters_get
    use Logfile_interface, ONLY : Logfile_stampMessage
    use Grid_data, ONLY : gr_meshMe, gr_globalComm
    implicit none
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  
    real, intent(IN)    :: dt
    integer, intent(IN) :: blockCount
    integer, dimension(blockCount), intent(IN):: blockList
    
    integer  ::  i, j, k, l, put, lb, ierr
    real     ::  x, y, z, vx, vy, vz, G, dist, mcoor
    real     ::  relax_rate, core_dist, mcs, distxy, vperp, vpara, vspin
    integer  ::  istat
  
    real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
    integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
    integer :: sizeX,sizeY,sizeZ
    real, dimension(:,:,:,:),pointer :: solnData
    real, dimension(EOS_NUM) :: eosData
    real, dimension(MDIM) :: avg_vel, new_vel, tot_avg_vel
    integer,dimension(MDIM) :: axis
    double precision, dimension(2) :: expl_coord, mpi_input, mpi_output
    double precision :: max_expl_kine
  
    logical :: gcell = .true.
    logical :: changedGrid

    character(len=128) :: str

    call PhysicalConstants_get("Newton", G)
    call Grid_getMinCellSize(mcs)
  
    if (dr_simTime .lt. sim_tInitial + sim_tRelax) then
        !adj_cfactor = sqrt(1.0-(dt*(1.0-sim_relaxRate)))
        relax_rate = max(0.0d0, dr_simTime - sim_tSpinup)/(sim_tRelax - sim_tSpinup)*(1.0 - sim_relaxRate) + sim_relaxRate 
        do lb = 1, blockCount
            call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
            sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
            allocate(xCoord(sizeX),stat=istat)
            sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
            allocate(yCoord(sizeY),stat=istat)
            sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
            allocate(zCoord(sizeZ),stat=istat)
  
            if (NDIM == 3) call Grid_getCellCoords&
                                (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
            if (NDIM >= 2) call Grid_getCellCoords&
                                (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
            call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
            call Grid_getBlkPtr(blockList(lb),solnData)
            do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                    do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                        x = xCoord(i) - gr_mpoleXcenterOfMass
                        y = yCoord(j) - gr_mpoleYcenterOfMass
                        z = zCoord(k) - gr_mpoleZcenterOfMass
                        vx = solnData(VELX_VAR,i,j,k)
                        vy = solnData(VELY_VAR,i,j,k)
                        vz = solnData(VELZ_VAR,i,j,k)
                        dist = sqrt(x**2.d0 + y**2.d0 + z**2.d0)
                        if (solnData(HE4_SPEC,i,j,k) .lt. 0.5d0) then
                            solnData(TEMP_VAR,i,j,k) = sim_accTemp
                            solnData(VELX_VAR,i,j,k) = vx*relax_rate
                            solnData(VELY_VAR,i,j,k) = vy*relax_rate
                            solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = core_xn
                        else
                            distxy = sqrt(x**2 + y**2)
                            vpara = (x*vx + y*vy)/distxy
                            if (dr_simTime .lt. sim_tSpinup) then
                                vspin = dr_simTime/sim_tSpinup*&
                                    min(sqrt(G*sim_accMass/distxy)*sim_rotFac, &
                                        sqrt(G*sim_accMass/radius(core_pos)**3.d0)*distxy)
                                !vperp = sqrt((vx - x/distxy*vpara)**2.d0 + (vy - y/distxy*vpara)**2.d0)
                                !if (vperp .gt. vspin) then
                                !    vperp = vperp - (sim_tSpinup - dr_simTime)/sim_tSpinup*(vperp - vspin)
                                !else
                                !    vperp = vperp + (sim_tSpinup - dr_simTime)/sim_tSpinup*(vspin - vperp)
                                !endif
                                solnData(VELX_VAR,i,j,k) = -vspin*y/distxy + x/distxy*vpara*relax_rate
                                solnData(VELY_VAR,i,j,k) =  vspin*x/distxy + y/distxy*vpara*relax_rate
                            else
                                solnData(VELX_VAR,i,j,k) = vx - x/distxy*vpara*(1.d0 - relax_rate)
                                solnData(VELY_VAR,i,j,k) = vy - y/distxy*vpara*(1.d0 - relax_rate)
                            endif
                            solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = torus_xn
                        endif
                        solnData(VELZ_VAR,i,j,k) = vz*relax_rate
                        
                        solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                            0.5*(solnData(VELX_VAR,i,j,k)**2. + solnData(VELY_VAR,i,j,k)**2. + solnData(VELZ_VAR,i,j,k)**2.)
                    enddo
                enddo
            enddo
  
            call Grid_releaseBlkPtr(blockList(lb), solnData)
            deallocate(xCoord)
            deallocate(yCoord)
            deallocate(zCoord)
            call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
        enddo
    else
        if (dr_simTime .ge. sim_tExplode .and. .not. exploded) then
            exploded = .true.
            core_dist = radius(core_pos)
            do lb = 1, blockCount
                call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
                sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
                allocate(xCoord(sizeX),stat=istat)
                sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
                allocate(yCoord(sizeY),stat=istat)
                sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
                allocate(zCoord(sizeZ),stat=istat)
      
                if (NDIM == 3) call Grid_getCellCoords&
                                    (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
                if (NDIM >= 2) call Grid_getCellCoords&
                                    (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
                call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
                call Grid_getBlkPtr(blockList(lb),solnData)
                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            if (solnData(HE4_SPEC,i,j,k) .lt. 0.5d0) then
                                solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = core_xn
                            else
                                solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = torus_xn
                            endif
                            if ((solnData(HE4_SPEC,i,j,k) .ge. 0.5d0) .and. &
                                (minval(solnData(HE4_SPEC,i,j+sim_detHeight:j+sim_detHeight+1,k)) .lt. 0.5d0) .and. &
                                (sqrt((xCoord(i) - gr_mpoleXcenterOfMass)**2.d0 + (zCoord(k) - gr_mpoleZcenterOfMass)**2.d0) .le. sim_detRadius*mcs)) then
                                solnData(DENS_VAR,i,j,k) = sim_detDens
                                solnData(TEMP_VAR,i,j,k) = sim_detTemp
                            endif
                        enddo
                    enddo
                enddo
      
                call Grid_releaseBlkPtr(blockList(lb), solnData)
                deallocate(xCoord)
                deallocate(yCoord)
                deallocate(zCoord)
                call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
            enddo
        elseif (sim_explodeCore .eq. .true.) then
            max_expl_kine = 0.0d0
            do lb = 1, blockCount
                call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
                sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
                allocate(xCoord(sizeX),stat=istat)
                sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
                allocate(yCoord(sizeY),stat=istat)
                sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
                allocate(zCoord(sizeZ),stat=istat)
      
                if (NDIM == 3) call Grid_getCellCoords&
                                    (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
                if (NDIM >= 2) call Grid_getCellCoords&
                                    (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
                call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
                call Grid_getBlkPtr(blockList(lb),solnData)
                do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                    do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                            if (solnData(HE4_SPEC,i,j,k) .lt. 0.1d0 .and. dabs(zCoord(k) - gr_mpoleZcenterOfMass) .le. mcs .and. &
                                solnData(DENS_VAR,i,j,k) .gt. sim_critDens .and. &
                                0.5d0*solnData(DENS_VAR,i,j,k)*sum(solnData(VELX_VAR:VELZ_VAR,i,j,k)**2.d0) .gt. max_expl_kine) then
                                max_expl_kine = 0.5d0*solnData(DENS_VAR,i,j,k)*sum(solnData(VELX_VAR:VELZ_VAR,i,j,k)**2.d0)
                                expl_coord = (/ xCoord(i), yCoord(j) /)
                            endif
                        enddo
                    enddo
                enddo
      
                call Grid_releaseBlkPtr(blockList(lb), solnData)
                deallocate(xCoord)
                deallocate(yCoord)
                deallocate(zCoord)
            enddo
            mpi_input = (/ max_expl_kine, dble(gr_meshMe) /)
            call MPI_ALLREDUCE (mpi_input, mpi_output, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, gr_globalComm, ierr)
            if (mpi_output(1) .gt. sim_critKine) then
                sim_explodeCore = .false.
                call MPI_BCAST(expl_coord, 2, FLASH_REAL, int(mpi_output(2)), gr_globalComm, ierr)
                if (gr_meshMe .eq. int(mpi_output(2))) then
                    call Logfile_stampMessage('Core exploded!')
                    write(str, *) 'Explosion coordinates: ', expl_coord
                    call Logfile_stampMessage(str)
                endif
                do lb = 1, blockCount
                    call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
                    sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
                    allocate(xCoord(sizeX),stat=istat)
                    sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
                    allocate(yCoord(sizeY),stat=istat)
                    sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
                    allocate(zCoord(sizeZ),stat=istat)
          
                    if (NDIM == 3) call Grid_getCellCoords&
                                        (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
                    if (NDIM >= 2) call Grid_getCellCoords&
                                        (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
                    call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
                    call Grid_getBlkPtr(blockList(lb),solnData)
                    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
                        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
                            do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                                if (sqrt((xCoord(i) - expl_coord(1))**2.d0 + (yCoord(j) - expl_coord(2))**2.d0 + &
                                         (zCoord(k) - gr_mpoleZcenterOfMass)**2.d0) .le. sim_detRadius*mcs) then
                                    solnData(TEMP_VAR,i,j,k) = sim_dbleDetTemp
                                    solnData(VELX_VAR:VELZ_VAR,i,j,k) = 0.d0
                                endif
                            enddo
                        enddo
                    enddo
          
                    call Grid_releaseBlkPtr(blockList(lb), solnData)
                    deallocate(xCoord)
                    deallocate(yCoord)
                    deallocate(zCoord)
                    call Eos_wrapped(MODE_DENS_TEMP,blkLimitsGC,blockList(lb))
                enddo
            endif
        endif
        call Stir(blockCount, blockList, dt) 
        call Burn(blockCount, blockList, dt) 
        call Heat(blockCount, blockList, dt, dr_simTime) 
        call Cool(blockCount, blockList, dt, dr_simTime)
    endif
  
  
    return
end subroutine Driver_sourceTerms
