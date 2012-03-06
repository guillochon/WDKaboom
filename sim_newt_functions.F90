MODULE sim_newt_functions

implicit none

CONTAINS

    FUNCTION newt_wd_mass_1d(x)
        use nrtype
        use Simulation_data
        use Logfile_interface
        use Grid_interface, ONLY: Grid_getMyPE
        implicit none
        real, dimension(:), intent(in) :: x
        real, dimension(size(x)) :: newt_wd_mass_1d
        double precision    :: wd_mass_tot
        integer             :: myPE
        character(len=32) :: int_to_str
    
        call Grid_getMyPE(myPE)
    
        !loop until mass is correct
        wd_mass_tot = 0.0
    
        radius = 0.0e0
        ipos = 0
        rhop = 0.0e0
        call wd_intout(x(1),sim_accTemp,wd_mass_tot,ipos,radius,rhop,m,eint,temper,omega,theta)
    
        call Logfile_stampMessage(myPE,'WD Mass 1D')
        write (int_to_str, '(e15.9)') radius(ipos(1))
        call Logfile_stampMessage(myPE,'Radius (equator): ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') wd_mass_tot
        call Logfile_stampMessage(myPE,'Mass: ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(i10)') ipos(1)
        call Logfile_stampMessage(myPE,'Profile length (equator): ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') abs(wd_mass_tot-sim_accMass-sim_torusMass)/(sim_accMass + sim_torusMass)
        call Logfile_stampMessage(myPE,'Mass Error: ' // trim(adjustl(int_to_str)))
        write (int_to_str, '(e15.9)') rhop(1,1)
        call Logfile_stampMessage(myPE,'Central Density: ' // trim(adjustl(int_to_str)))
        call Logfile_stampMessage(myPE,'')
    
        newt_wd_mass_1d(1) = (wd_mass_tot - sim_accMass - sim_torusMass)/(sim_accMass + sim_torusMass)
    END FUNCTION newt_wd_mass_1d
    
END MODULE sim_newt_functions
