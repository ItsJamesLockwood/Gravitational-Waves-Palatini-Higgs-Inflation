!!The model here is V=lambda * Mpl**2/(4 xi**2) * tanh(xi * chi/ Mpl) (default model in LatticeEasy and DEFROST, see hep-ph 9705347)

module model
  use define_fields
  implicit none

#include "configure.h"
!!*******************define the couplings etc. for your model *************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!the predefined constants that you can use: GeV, MPl (the reduced Planck Mass), PlanckMass (the Planck Mass), Mplsq (square of the reduced Planck Mass)
  real(dl),parameter:: lambda =1.d-4
  real(dl),parameter:: Nstar = 50
  real(dl),parameter:: xi = 3.8d6 * Nstar**2 * lambda
  real(dl),parameter:: xi2 = xi**2
  real(dl),parameter:: xisqrt = xi**0.5 
  
  real(dl),parameter:: a = lambda * Mplsq**2 / 4.d0 / xi2
  real(dl),parameter:: b = xisqrt/Mpl
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
    potential = a * TANH( b * PHI )**4 
  end function potential

!! the derivative of pential w.r.t. to the fields
!! return an array (dV/df1, dV/df2, ...)
  function dVdf(f)
    real(dl) f(ns)
    real(dl),dimension(ns):: dVdf
    dVdf = (/ &
         4.d0 * a * b * TANH( b * PHI )**3 / COSH(b * PHI )**2 &
         /)
  end function dVdf

!!d^2 V / d f_{fld}^2,, here fld can be 1,2,..., ns
  function mass_sq(f,fld)
    real(dl) f(ns)
    real(dl):: mass_sq
    integer fld
    select case(fld)
    case(1)
       mass_sq = 4.d0*a* b**2 * (3.d0*TANH(b*PHI)**2 / COSH(b*PHI)**4 - 2.d0*TANH(b*PHI)**4 / COSH(b*PHI)**2 )
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
    logical,save::warning = .true.
    integer(IB) fld
    if(n*k*metric%dx .lt. const_pi) then
      model_Power =0._dl
      return
    endif
    !check that k > k_min=a * H_init (note: k in code is defined as comoving momentum):
    if(k.gt.Init_Hubble)then 
      !define omega^2:
      omega=k**2+mass_sq(init_fields,fld)
      !check if omega / effective mass squared is positive
      if(omega.gt.0.)then
         omega=sqrt(abs(omega))
         !TODO: why is this statement necessary?
         if(omega*metric%dx .le. const_2pi .and. omega*metric%dx*n .ge. const_2pi)then
            model_Power(1) = 0.5_dl/omega
            model_Power(2) = 0.5_dl*omega
            return
         endif
      !We now consider the tachyonic region (i.e. k_min < k < k_max)
      else
         !Set the initial conditions from arxiv:1902.10148. 
         model_Power(1) = Mpl/k /metric%a**3
         model_Power(2) = xisqrt/SQRT(lambda)* k/Mpl /metric%a**3
         if(warning)then
            write(*,*) "Tachyonic region initialization may be not correct"
            warning = .false.
         endif
         return
      endif
   !When k < k_min: no longer in tachyonic region
   else
      model_Power(1) = 0.5_dl/k*(init_Hubble/k)**2
      model_Power(2) = 0._dl
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
       write(fp%unit,*) "V = lambda * Mpl^2 /(4 xi^2) tanh(xi^.5 chi/Mpl)"
       write(fp%unit,*) "reduced Planck Mass M_p = ", Mpl
       write(fp%unit,*) "lambda = ",lambda
       write(fp%unit,*) "xi = ",xi
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




