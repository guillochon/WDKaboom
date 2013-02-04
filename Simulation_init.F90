!!****if* source/Simulation/SimulationMain/WDAccrection/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!! ARGUMENTS
!!
!!   myPE -   current processor number
!!
!! PARAMETERS
!!
!!  sim_pAmbient       Initial ambient pressure
!!  sim_rhoAmbient     Initial ambient density
!!  sim_xctr           Explosion center coordinates
!!  sim_yctr           Explosion center coordinates
!!  sim_zctr           Explosion center coordinates
!!  sim_nsubzones      Number of `sub-zones' in cells for applying 1d profile
!!
!!***

subroutine Simulation_init(myPE)

    use Simulation_data 
    use Driver_data, ONLY : dr_restart, dr_simTime
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Logfile_interface
    use Grid_interface, ONLY : Grid_getMinCellSize
    use newt_wrappers
    use sim_newt_functions

    implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

    integer, intent(in) :: myPE
    integer             :: i, j, k, ierr, istat, ii, jj, kk, gdim, jLo, jHi, n, cnt, tLo, tHi
    double precision    :: start_t, rho0, rho1, wd_mass_tot, old_mass_tot, min_grid, &
                           sumRho, mcell, xx, yy, zz, cdist, cth, &
                           xDist, yDist, zDist, xCoord, yCoord, zCoord, frac, tfrac, grid_center, dist, th, dmass_drho
    character(len=32), dimension(1,2) :: block_buff
    character(len=32) :: int_to_str
    double precision,allocatable,dimension(:,:,:) :: grid_3d
    logical             :: calc_3d = .false., is_core
    double precision, dimension(1) :: input, output

    sim_pi = PI
    call RuntimeParameters_get('sim_tAmbient', sim_tAmbient)
    call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
    call RuntimeParameters_get('sim_nsubzones',sim_nSubZones)
    call RuntimeParameters_get('sim_xctr',sim_xCenter)
    call RuntimeParameters_get('sim_yctr',sim_yCenter)
    call RuntimeParameters_get('sim_zctr',sim_zCenter)
    call RuntimeParameters_get('smallx', sim_smallX)
    call RuntimeParameters_get('smlrho', sim_smallRho)
    call RuntimeParameters_get('smallp', sim_smallP)
    call RuntimeParameters_get('xmin',sim_xMin)
    call RuntimeParameters_get('ymin',sim_yMin)
    call RuntimeParameters_get('zmin',sim_zMin)
    call RuntimeParameters_get('xmax',sim_xMax)
    call RuntimeParameters_get('ymax',sim_yMax)
    call RuntimeParameters_get('zmax',sim_zMax)
    call RuntimeParameters_get('tinitial',sim_tInitial)
    call RuntimeParameters_get('sim_tRelax',sim_tRelax)
    call RuntimeParameters_get('sim_tExplode',sim_tExplode)
    call RuntimeParameters_get('sim_tSpinup',sim_tSpinup)
    call RuntimeParameters_get('sim_relaxRate',sim_relaxRate)
    call RuntimeParameters_get('sim_profileSubdiv',sim_profileSubdiv)
    call RuntimeParameters_get('sim_minProfDelta',sim_minProfDelta)
    call RuntimeParameters_get('sim_maxProfDelta',sim_maxProfDelta)
    call RuntimeParameters_get('sim_rhoGuess',sim_rhoGuess)
    call RuntimeParameters_get('sim_accMass',sim_accMass)
    call RuntimeParameters_get('sim_accTemp',sim_accTemp)
    call RuntimeParameters_get('sim_torusMass',sim_torusMass)
    call RuntimeParameters_get('sim_virTemp',sim_virTemp)
    call RuntimeParameters_get('sim_rotFac',sim_rotFac)
    call RuntimeParameters_get('sim_massAcc',sim_massAcc)
    call RuntimeParameters_get('sim_axisRatio',sim_axisRatio)
    call RuntimeParameters_get('sim_maxBlocks',sim_maxBlocks)
    call RuntimeParameters_get('sim_detDens',sim_detDens)
    call RuntimeParameters_get('sim_detTemp',sim_detTemp)
    call RuntimeParameters_get('sim_detRadius',sim_detRadius)
    call RuntimeParameters_get('sim_critDens',sim_critDens)
    call RuntimeParameters_get('sim_critKine',sim_critKine)
    call RuntimeParameters_get('sim_dbleDetTemp',sim_dbleDetTemp)
    call RuntimeParameters_get('sim_explodeCore',sim_explodeCore)
    call RuntimeParameters_get('sim_detHeight',sim_detHeight)

    if (sim_nSubZones .le. 1) sim_nSubZones = 2

    sim_inSubZones = 1./real(sim_nSubZones)
    sim_inSubzm1   = 1./real(sim_nSubZones-1)
    sim_inszd      = sim_inSubZones**NDIM
    
    core_xn = sim_smallX
    core_xn(C12_SPEC) = 0.5d0
    core_xn(O16_SPEC) = 0.5d0

    torus_xn = sim_smallX
    torus_xn(HE4_SPEC) = 1.d0

    if (sim_rhoGuess .eq. 0.0d0) then
        if (myPE .eq. MASTER_PE) then
            call Logfile_stampMessage(myPE,'[Simulation Init] Guessing Central Density')
        endif
        call wd_guess_rho_0(sim_accMass+sim_torusMass,sim_rhoGuess)
    endif

    if (myPE .eq. MASTER_PE) then
        input = sim_rhoGuess
        call run_newt(newt_wd_mass_1d, 1, sim_massAcc, input, output)

        !massacc = 5.d-3
        !wd_mass_tot = 0.0
        !n = 0
        !do while (abs(wd_mass_tot-sim_accMass-sim_torusMass)/(sim_accMass+sim_torusMass).gt.massacc)
        !    rho1 = rho0 - (wd_mass_tot - sim_accMass - sim_torusMass)/dmass_drho
        !    if(rho1.gt.1.1e0*rho0) rho1 = (1.d0+1.d-1*(1.d3 - n)/1.d3)*rho0
        !    if(rho1.lt.0.9e0*rho0) rho1 = (1.d0-1.d-1*(1.d3 - n)/1.d3)*rho0
        !    write (int_to_str, '(i11)') n
        !    call Logfile_stampMessage(myPE,'[Simulation Init] Generating Accretor in 3D, Iteration #' // &
        !        trim(adjustl(int_to_str)))

        !    call wd_intout(rho1,sim_accTemp,wd_mass_tot,ipos,radius,rhop,m,eint,temper,omega,theta)

        !    gdim = 2*ceiling(radius(maxval(ipos))/min_grid)
        !    grid_center = gdim/2*min_grid
        !    allocate(grid_3d(0:gdim-1,0:gdim-1,0:gdim-1),stat=istat)
        !    grid_3d = 0.d0

        !    do k = 0, gdim-1
        !        zCoord = min_grid*(k + 0.5d0)
        !        do j = 0, gdim-1
        !            yCoord = min_grid*(j + 0.5d0)
        !            do i = 0, gdim-1
        !                xCoord = min_grid*(i + 0.5d0)
        !                sumRho = 0.0
        !                xx = xCoord - grid_center
        !                yy = yCoord - grid_center
        !                zz = zCoord - grid_center
        !                cdist = sqrt(xx**2 + yy**2 + zz**2)
        !                cth = atan2(abs(zz), sqrt(xx**2.d0 + yy**2.d0))
        !                call sim_find (theta, ns, cth, tLo, tHi, tfrac)
        !                call sim_find (radius, max(ipos(tLo), ipos(tHi)), cdist, jLo, jHi, frac)
        !                mcell = m(jLo) + frac*(m(jHi) - m(jLo))
        !                if (mcell .lt. sim_accMass) then
        !                    is_core = .true.
        !                else
        !                    is_core = .false.
        !                endif

        !                cnt = 0

        !                do kk = 0, sim_nSubZones-1
        !                   zDist = zCoord + (kk*sim_inSubzm1 - 0.5d0)*min_grid - grid_center
        !                   do jj = 0, sim_nSubZones-1
        !                      yDist = yCoord + (jj*sim_inSubzm1 - 0.5d0)*min_grid - grid_center
        !                      do ii = 0, sim_nSubZones-1
        !                         xDist = xCoord + (ii*sim_inSubzm1 - 0.5d0)*min_grid - grid_center
        !                         dist = sqrt(xDist**2 + yDist**2 + zDist**2)
        !                         th = atan2(abs(zDist), sqrt(xDist**2.d0 + yDist**2.d0))
        !                         call sim_find(theta, ns, th, tLo, tHi, tfrac)
        !                         call sim_find(radius, max(ipos(tLo), ipos(tHi)), dist, jLo, jHi, frac)

        !                         if (jLo .eq. jHi) cycle 
        !                        
        !                         mcell = m(jLo) + frac*(m(jHi) - m(jLo))
        !                         if ((is_core .and. mcell .lt. sim_accMass) .or. &
        !                             (.not. is_core .and. mcell .ge. sim_accMass)) then
        !                             sumRho = sumRho + (rhop(jLo, tLo) + frac*(rhop(jHi, tLo) - rhop(jLo, tLo)) + &
        !                                 tfrac*(rhop(jLo, tHi) - rhop(jLo, tLo)) + &
        !                                 frac*tfrac*(rhop(jLo, tLo) - rhop(jHi, tLo) - rhop(jLo, tHi) + rhop(jHi, tHi)))

        !                             cnt = cnt + 1
        !                         endif
        !                      enddo
        !                   enddo
        !                enddo
        !                if (cnt .gt. 0) grid_3d(i,j,k) = sumRho/cnt
        !            enddo
        !        enddo
        !    enddo
        !    wd_mass_tot = sum(grid_3d)*min_grid**3.0
        !    deallocate(grid_3d)

        !    write (int_to_str, '(e10.4)') radius(ipos(1))
        !    call Logfile_stampMessage(myPE,'Radius: (equator): ' // trim(adjustl(int_to_str)))
        !    write (int_to_str, '(e10.4)') radius(ipos(ns))
        !    call Logfile_stampMessage(myPE,'Radius: (pole): ' // trim(adjustl(int_to_str)))
        !    write (int_to_str, '(e10.4)') wd_mass_tot
        !    call Logfile_stampMessage(myPE,'Mass: ' // trim(adjustl(int_to_str)))
        !    write (int_to_str, '(i10)') ipos(1)
        !    call Logfile_stampMessage(myPE,'Profile length (equator): ' // trim(adjustl(int_to_str)))
        !    write (int_to_str, '(e10.4)') abs(wd_mass_tot-sim_accMass-sim_torusMass)/(sim_accMass+sim_torusMass)
        !    call Logfile_stampMessage(myPE,'Mass Error: ' // trim(adjustl(int_to_str)))
        !    write (int_to_str, '(e10.4)') rhop(1,1)
        !    call Logfile_stampMessage(myPE,'Central Density: ' // trim(adjustl(int_to_str)))

        !    if (n .gt. 1000) then
        !        print *, 'wd mass did not converge!'
        !        stop
        !    endif
        !    dmass_drho = (wd_mass_tot - old_mass_tot)/(rho1 - rho0)
        !    rho0 = rho1
        !    old_mass_tot = wd_mass_tot
        !    n = n + 1
        !enddo

    endif

    call MPI_BCAST(radius, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rhop, np*ns, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eint, np*ns, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(temper, np*ns, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(omega, np*ns, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(theta, ns, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(m, np, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ipos, ns, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(core_pos, ns, FLASH_INTEGER, MASTER_PE, MPI_COMM_WORLD, ierr)

    if (dr_restart .eq. .true. .and. dr_simTime .gt. sim_tExplode) then
        exploded = .true.
    else
        exploded = .false.
    endif
    wd_radius = radius(ipos(1))

end subroutine Simulation_init
