!!The model here is V=lambda/4 (h^2 - v^2)^2 (default model in LatticeEasy and DEFROST, see hep-ph 9705347)

module model
  use define_fields
  implicit none

#include "configure.h"
!!*******************define the couplings etc. for your model *************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!the predefined constants that you can use: GeV, MPl (the reduced Planck Mass), PlanckMass (the Planck Mass), Mplsq (square of the reduced Planck Mass)
  real(dl):: lambda =1.d-14
  real(dl),parameter:: viv = 1.e-3_dl * PlanckMass
  real(dl),parameter:: v2 = viv**2
  
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************define macros here;************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#define PHI f(1) 
!! Note: PHI was used to denote the Higgs field h
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!!***************** define initial conditions ******************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!HLattice always treat the first field as the inflaton \phi. 
!!To be self consistent you may want to get the background \phi and \dot\phi by solving the background equations from inflation. HLattice will do that for you.
!!Just define where you want to start the background evolution. Note this is NOT where you start the lattice simulation, which is defined in a subroutine "start_box_simulation" in this file.
 !!initial field values
  real(dl),dimension(ns)::init_fields=(/ &
       3.5_dl *PlanckMass &
       /)

!!Initial field momenta
 !! if set to be greater than or equal to Mplsq (square of reduced Planck Mass), the initial field momenta will be determined by slow-roll inflationary attractor
  !! Note again these are NOT the initial field momenta where you start the lattice simulation, the are the initial values that HLattice take to evolve the inflaton. 
  real(dl),dimension(ns):: init_momenta = Mplsq


!!put initial random Gaussian perturbations in the fields when you starts lattice simulation;
!!the amplitude of fluctuations are defined in subroutine model_Power
  logical,dimension(ns)::do_init = (/ .true. /)

!!Important note: init_fields and init_momenta will be changed after the initialization. After the subroutine init() is called, they will equal to the the fields and field momenta AT THE BEGINNING OF LATTICE SIMULATION. In addition, the Hubble parameter at the beginning of lattice simulation will be saved to a global variable "init_Hubble".
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

contains
  
!! the potential of the scalar fields
  function potential(f)
    real(dl) f(ns)
    real(dl) potential
    potential = (lambda/4.d0) * ( PHI**2 - v2)**2 
  end function potential

!! the derivative of pential w.r.t. to the fields
!! return an array (dV/df1, dV/df2, ...)
  function dVdf(f)
    real(dl) f(ns)
    real(dl),dimension(ns):: dVdf
    dVdf = (/ &
         lambda*PHI*( PHI**2 - v2) &
         /)
  end function dVdf

!!d^2 V / d f_{fld}^2,, here fld can be 1,2,..., ns
  function mass_sq(f,fld)
    real(dl) f(ns)
    real(dl):: mass_sq
    integer fld
    select case(fld)
    case(1)
       mass_sq = lambda * (3.d0* PHI**2 -v2)
    case default
       stop "wrong argument fld in mass_sq"
    end select
  end function mass_sq

!!define the condition to start box simulation
!!here array f -> the background field values (\phi_1, \phi_2, ...)
!! array p-> the background field derivative, (d\phi_1/dt, d\phi_1/dt, ... )
  function start_box_simulation(f,p)
    logical start_box_simulation
    real(dl),dimension(ns)::f,p
!!LatticeEasy default, i.e. d \phi/dt + H phi = 0 (or equivalently d(a\phi)/dt =0
    start_box_simulation = (p(1) + sqrt((sum(p**2)/2.+potential(f))/3.)/Mpl*f(1) .lt. 0.)
!! If you want to start the lattice simualtion immediately
    !!start_box_simulation = .true.
!! If you want to use DEFROST starting point, the end of inflation, i.e. \ddot a =0
    !!start_box_simulation = (potential(f) .le. sum(p**2) )
  end function start_box_simulation

!!define the time step; you will get significant error if you use physical time and 2nd order integrator (the integrator is chosen in parameters.f90)
!!I strongly recommend you use conformal time for this model, for which a second order integrator will do ok (if metric perturbs are off)
  function model_dt_per_step()
    real(dl) model_dt_per_step
#if USE_CONFORMAL_TIME
    !!when metric backreaction is off, just use conformal time; Easy and very stable;
    model_dt_per_step = metric%dx/32._dl
#else
    !!when metric backreaction is on, a bit tricky for this model cz HLattice use synchronous gauge (and physical time t).
    !!Better not continuously change dt; Use step functions.
    model_dt_per_step = metric%dx/64._dl*ceiling(metric%a/5._dl)
#endif
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
!!A check for lambda phi^4 model, for which conformal time is a better choice.
!!#if ! USE_CONFORMAL_TIME && ! METRIC_PERTURB
!!       stop "you may want to turn on USE_CONFORMAL_TIME in configure.h"
!!#endif
    case(2)
       fp = open_file("data/"//trim(run_name)//"_model.info", "w")
       write(fp%unit,*) "V = lambda/4 (h^2 - v^2)"
       write(fp%unit,*) "reduced Planck Mass M_p = ", Mpl
       write(fp%unit,*) "lambda = ",lambda
       write(fp%unit,*) "n = ",n
       write(fp%unit,*) "dx=",metric%dx
       write(fp%unit,*) "Initial fields values:", Init_fields
       write(fp%unit,*) "Initial momenta:", Init_momenta
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




