!!model in arXiv: 0902.2197

module model
  use define_fields
  implicit none

#include "configure.h"


!!*******************define the couplings etc. for your model *************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!the predefined constants that you can use: GeV, MPl (the reduced Planck Mass), PlanckMass (the Planck Mass), Mplsq (square of the reduced Planck Mass)
  real(dl),parameter:: m_s = 100._dl * GeV
  real(dl),parameter:: M_cutoff = Mpl
  real(dl),parameter:: K_run = -0.1
  real(dl),parameter:: phi0 = 1.e16_dl*GeV
  real(dl),parameter:: EpsilonPlus = 1.e-8_dl !!this number is used to avoid log(0) glitches.
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************define macros here;************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#define REPHI f(1)
#define IMPHI f(2)
#define PHI2 (REPHI **2 + IMPHI **2)
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************** define initial conditions ******************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!HLattice always treat the first field as the inflaton \phi. 
!!To be self consistent you may want to get the background \phi and \dot\phi by solving the background equations from inflation. HLattice will do that for you.
!!Just define where you want to start the background evolution. Note this is NOT where you start the lattice simulation, which is defined in a subroutine "start_box_simulation" in this file.
 !!initial field values
 real(dl),dimension(ns)::init_fields=(/ &
       phi0, 0._dl &
       /)

!!Initial field momenta
 !! if set to be greater than or equal to Mplsq (square of reduced Planck Mass), the initial field momenta will be determined by slow-roll inflationary attractor
  !! Note again these are NOT the initial field momenta where you start the lattice simulation, the are the initial values that HLattice take to evolve the inflaton. 
  real(dl),dimension(ns):: init_momenta = 0.


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
    potential = (m_s**2/2._dl) * PHI2 * (1._dl + K_run * log(PHI2/(2._dl*M_cutoff**2)+EpsilonPlus)) !!I use a epsilon+ to avoid log(zero) numerical instability.
  end function potential

!! the derivative of pential w.r.t. to the fields
!! return an array (dV/df1, dV/df2, ...)
  function dVdf(f)
    real(dl) f(ns)
    real(dl),dimension(ns):: dVdf
    dVdf = (m_s**2)*(1._dl +  K_run * (1.-EpsilonPlus/(PHI2/(2._dl*M_cutoff**2)+EpsilonPlus)+log(PHI2/(2._dl*M_cutoff**2)+EpsilonPlus))) * f
  end function dVdf

!!d^2 V / d f_{fld}^2,, here fld can be 1,2,..., ns
  function mass_sq(f,fld)
    real(dl) f(ns)
    real(dl):: mass_sq
    integer fld
    select case(fld)
    case(1)
       mass_sq = (m_s**2) * (1._dl +  K_run * (1.+ 2.*REPHI**2/PHI2 + log(PHI2/(2._dl*M_cutoff**2))))
    case(2)
       mass_sq = (m_s**2) * (1._dl +  K_run * (1.+ 2.*IMPHI**2/PHI2 + log(PHI2/(2._dl*M_cutoff**2))))
    case default
       stop "wrong argument fld in mass_sq"
    end select
  end function mass_sq

!!define the condition to start the lattice simulation
  function start_box_simulation(f,p)
    logical start_box_simulation
    real(dl),dimension(ns)::f,p
    start_box_simulation = .true.
    !!immediately start the lattice simulation
  end function start_box_simulation

!!define the time step;
  function model_dt_per_step()
    real(dl) model_dt_per_step
    model_dt_per_step = metric%dx/50._dl

  end function model_dt_per_step

!!return an array (|f_k|^2, |f'_k|^2) for fields initialization. 
!!remember always set them to be zero for k=0.
  function model_Power(fld,k)
    real(dl),dimension(2):: model_Power
    real(dl) k,omega
    integer(IB) fld
    if(n*k*metric%dx .lt. const_pi) then
       model_Power =0._dl
       return
    endif
    model_Power(1) = (const_pi**2 * 2._dl)/k*((5.e-5*phi0)/k)**2
    model_Power(2) = 0._dl
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
#if USE_CONFORMAL_TIME
       stop "you may want to turn off USE_CONFORMAL_TIME in configure.h"
#endif
    case(2)
       fp = open_file("data/"//trim(run_name)//"_model.info", "w")
       write(fp%unit,*) "flat direction arXiv: 0902.2197"
       write(fp%unit,*) "reduced Planck Mass M_p = ", Mpl
       write(fp%unit,*) "m_{3/2} = ",m_s
       write(fp%unit,*) "M = ", M_cutoff
       write(fp%unit,*) "K = ", K_run
       write(fp%unit,*) "n = ",n
       write(fp%unit,*) "dx=",metric%dx
       write(fp%unit,*) "m_{3/2} / Hubble = ", m_s/Init_Hubble
       write(fp%unit,*) "effective k_min / Hubble = ", k_unit/metric%physdx/Init_Hubble
       write(fp%unit,*) "effective k_max / Hubble = ", fft_numk*k_unit/metric%physdx/Init_Hubble
       call close_file(fp)
       write(*,*) "m_{3/2} / Hubble = ", m_s/Init_Hubble
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
