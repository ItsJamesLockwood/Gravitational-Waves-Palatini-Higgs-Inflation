!!The model here is V=lambda/4 phi^4 + 1/2 g^2 phi^2 chi^2 (default model in LatticeEasy and DEFROST, see hep-ph 9705347)

module model
  use define_fields
  implicit none

#include "configure.h"
!!*******************define the couplings etc. for your model *************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!the predefined constants that you can use: GeV, MPl (the reduced Planck Mass), PlanckMass (the Planck Mass), Mplsq (square of the reduced Planck Mass)
  real(dl):: lambda = 7.18d-14
  real(dl):: kappa = MPl * 1.94d-12
  real(dl):: m = MPl * 6.82d-6
  real(dl):: j = 47998.
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************define macros here;************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#define PHI f(1)
#define CHI f(2)
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!!***************** define initial conditions ******************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!HLattice always treat the first field as the inflaton \phi. 
!!To be self consistent you may want to get the background \phi and \dot\phi by solving the background equations from inflation. HLattice will do that for you.
!!Just define where you want to start the background evolution. Note this is NOT where you start the lattice simulation, which is defined in a subroutine "start_box_simulation" in this file.
 !!initial field values
  real(dl),dimension(ns)::init_fields=(/ &
       1.358d0 * MPl, 0._dl &
       /)

!!Initial field momenta
 !! if set to be greater than or equal to Mplsq (square of reduced Planck Mass), the initial field momenta will be determined by slow-roll inflationary attractor
  !! Note again these are NOT the initial field momenta where you start the lattice simulation, the are the initial values that HLattice take to evolve the inflaton. 
  real(dl),dimension(ns):: init_momenta = Mplsq


!!put initial random Gaussian perturbations in the fields when you starts lattice simulation;
!!the amplitude of fluctuations are defined in subroutine model_Power
  logical,dimension(ns)::do_init = (/ .true., .true. /)

!!Important note: init_fields and init_momenta will be changed after the initialization. After the subroutine init() is called, they will equal to the the fields and field momenta AT THE BEGINNING OF LATTICE SIMULATION. In addition, the Hubble parameter at the beginning of lattice simulation will be saved to a global variable "init_Hubble".
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

contains
  
!! the potential of the scalar fields
  function potential(f)
    real(dl) f(ns)
    real(dl) potential
    potential = (lambda/4.d0) * PHI**4 - (2.d0 * kappa/3.d0) * PHI**3 + ((m**2)/2.d0) * PHI**2 + ((lambda/4.d0) * (j + 2.d0) * (j + 3.d0) * PHI**2 - kappa * (j + 2.d0) * PHI + (m**2)/2.d0) * CHI**2
  end function potential

!! the derivative of pential w.r.t. to the fields
!! return an array (dV/df1, dV/df2, ...)
  function dVdf(f)
    real(dl) f(ns)
    real(dl),dimension(ns):: dVdf
    dVdf = (/ &
         lambda * PHI**3 - (2.d0 * kappa) * PHI**2 + (m**2) * PHI + ((lambda/2.d0) * (j + 2.d0) * (j + 3.d0) * PHI - kappa * (j + 2.d0)) * CHI**2 ,   &  
         2.d0 * ((lambda/4.d0) * (j + 2.d0) * (j + 3.d0) * PHI**2 - kappa * (j + 2.d0) * PHI + (m**2)/2.d0) * CHI &
         /)
  end function dVdf

!!d^2 V / d f_{fld}^2,, here fld can be 1,2,..., ns
  function mass_sq(f,fld)
    real(dl) f(ns)
    real(dl):: mass_sq
    integer fld
    select case(fld)
    case(1)
       mass_sq = 3.d0 * lambda * PHI**2 - (4.d0 * kappa) * PHI + m**2 + ((lambda/2.d0) * (j + 2.d0) * (j + 3.d0)) * CHI**2
    case(2)
       mass_sq = 2.d0 * ((lambda/4.d0) * (j + 2.d0) * (j + 3.d0) * PHI**2 - kappa * (j + 2.d0) * PHI + (m**2)/2.d0)
    case default
       stop "wrong argument fld in mass_sq"
    end select
  end function mass_sq

!!define the condition to start box simulation
  function start_box_simulation(f,p)
    logical start_box_simulation
    real(dl),dimension(ns)::f,p
    start_box_simulation = mass_sq(f, 2) .lt. 1.e-8*Mplsq
    !!start_box_simulation = (p(1) + sqrt((sum(p**2)/2.+potential(f))/3.)/Mpl*f(1) .lt. 0.)
    !!LATTICEEASY initial condition

    !!start_box_simulation = (sum(p**2) .ge. potential(f) ) 
    !!i.e. epsilon = 1 (the end of inflation) 

    !! or starting much later is actually ok, you can still catch the main dynamics
    !!    start_box_simulation = (sum(p**2) / 2. + potential(f) .lt. (0.35e-9_dl*Mplsq)**2)
  end function start_box_simulation


!!define the time step; you will get significant error if you use physical time and 2nd order integrator (the integrator is chosen in configure.h)
!!I strongly recommend you use conformal time for this model, for which a second order integrator will do ok (if metric perturbs are off)
  function model_dt_per_step()
    real(dl) model_dt_per_step
!!to capture the rapid oscillations of the chi field
    model_dt_per_step = metric%dx/2048._dl
  end function model_dt_per_step

!!return an array (|f_k|^2, |f'_k|^2) for fields initialization. 
!!remember always set them to be zero for k=0.
  function model_Power(fld,k)
    real(dl),dimension(2):: model_Power
    real(dl) k,omega
    integer(IB) fld
    omega=k**2+mass_sq(init_fields,fld)
    if(n*k*metric%dx .lt. const_pi .or. omega.lt.0.) then
       model_Power =0._dl
       return
    endif
    omega=sqrt(omega)
    if(omega.gt.Init_Hubble)then 
       !!put in 1/2 particle per phase space volume
       !!note this is just an approximation. If you care about time scale of one or two oscillations. The quantum->classical transition should be calculated analytically.
       model_Power(1) = 0.5_dl/omega
       model_Power(2) = 0.5_dl*omega
       return
    endif
    model_Power = 0._dl
    return
  end function model_Power


!!model outputs
  subroutine model_output(i)
    !! i = 1; called before initialization
    !! i = 2; called after initialization
    !! i = 3; called in output_to_screen()
    !! i = 4; called at the end of simulation
    integer(IB) i 
    type(file_pointer) fp
    select case(i)
    case(1)
    case(2)
       fp = open_file("data/"//trim(run_name)//"_model.info", "w")
       write(fp%unit,*) "V = lambda/4 phi^4 + g^2/2 chi^2 phi^2"
       write(fp%unit,*) "reduced Planck Mass M_p = ", Mpl
       !!write(fp%unit,*) "lambda = ",lambda
       !!write(fp%unit,*) "g^2 = ", g2l*lambda
       write(fp%unit,*) "n = ",n
       write(fp%unit,*) "dx=",metric%dx
       write(fp%unit,*) "Initial fields values:", Init_fields
       write(fp%unit,*) "Initial momenta:", Init_momenta
       write(fp%unit,*) "Initial mass (/M_p):", sqrt(mass_sq(init_fields, 1))/Mpl, sqrt(mass_sq(init_fields, 2))/Mpl
       write(fp%unit,*) "effective k_min / Hubble = ", k_unit/metric%physdx/Init_Hubble
       write(fp%unit,*) "effective k_max / Hubble = ", fft_numk*k_unit/metric%physdx/Init_Hubble
       call close_file(fp)
       write(*,*) "effective k_min / Hubble = ", k_unit/metric%physdx/Init_Hubble
       write(*,*) "effective k_max / Hubble = ", fft_numk*k_unit/metric%physdx/Init_Hubble
    case(3)
    case(4)
    case default
       write(*,*) "Unknow option in model_output"
       stop
    end select
  end subroutine model_output

end module model
