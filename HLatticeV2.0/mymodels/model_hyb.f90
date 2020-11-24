!!Hybrid Inflation model. See arxiv 0812.2917
!!change NUM_SCALAR_FIELDS to 3 in configure.h before you use this model
module model
  use define_fields
  implicit none

#include "configure.h"

!!*******************define the couplings etc. for your model *************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!the predefined constants that you can use: GeV, MPl (the reduced Planck Mass), PlanckMass (the Planck Mass), Mplsq (square of the reduced Planck Mass)
  real(dl),parameter:: sqrt_lambda = 1.e-7_dl
  real(dl),parameter:: lambda = sqrt_lambda**2
  real(dl),parameter:: cpl_g = const_sqrt2 * sqrt_lambda
  real(dl),parameter:: g2 = cpl_g**2
  real(dl),parameter:: viv = 1.e-3_dl * PlanckMass
  real(dl),parameter:: v2 = viv**2
  real(dl),parameter:: V_c = 1.e-5_dl
  real(dl),parameter:: phi_c = sqrt_lambda / cpl_g * viv
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************define macros here;************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#define PHI f(1)
#define PHI2 (PHI**2)
#define RECHI f(2)
#define IMCHI f(3)
#define CHI2 (RECHI**2+IMCHI**2)
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!***************** define initial conditions ******************************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!HLattice always treat the first field as the inflaton \phi. 
!!To be self consistent you may want to get the background \phi and \dot\phi by solving the background equations from inflation. HLattice will do that for you.
!!Just define where you want to start the background evolution. Note this is NOT where you start the lattice simulation, which is defined in a subroutine "start_box_simulation" in this file.
 !!initial field values
  real(dl),dimension(ns)::init_fields=(/ &
       phi_c*0.9975, 0._dl, 0._dl &
       /)

!!Initial field momenta
 !! if set to be greater than or equal to Mplsq (square of reduced Planck Mass), the initial field momenta will be determined by slow-roll inflationary attractor
  !! Note again these are NOT the initial field momenta where you start the lattice simulation, the are the initial values that HLattice take to evolve the inflaton. 
  real(dl),dimension(ns):: init_momenta = (/ &
        -V_c * lambda*v2/cpl_g, 0._dl, 0._dl  & 
       /)

!!put initial random Gaussian perturbations in the fields when you starts lattice simulation;
!!the amplitude of fluctuations are defined in subroutine model_Power
  logical,dimension(ns)::do_init = (/ .false., .true., .true. /)

!!Important note: init_fields and init_momenta will be changed after the initialization. After the subroutine init() is called, they will equal to the the fields and field momenta AT THE BEGINNING OF LATTICE SIMULATION. In addition, the Hubble parameter at the beginning of lattice simulation will be saved to a global variable "init_Hubble".
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


contains

  !!the potential  
  function potential(f)
    real(dl) f(ns)
    real(dl) potential
    potential = (lambda/4._dl) * (CHI2 - v2)**2 + (g2/2._dl) * PHI**2 * CHI2
  end function potential

  !!the derivative of potential w.r.t. the fields
  function dVdf(f)
    real(dl) f(ns)
    real(dl),dimension(ns):: dVdf
    dVdf = (/ &
         g2*CHI2*PHI, &
         lambda *RECHI*(CHI2 - v2) + g2 * RECHI * PHI**2,     &
         lambda *IMCHI*(CHI2 - v2) + g2 * IMCHI * PHI**2     &
         /)
  end function dVdf

!!d^2V/df^2, fld can be 1, 2,.., NUM_SCALAR_FIELDS
  function mass_sq(f,fld)
    real(dl) f(ns)
    real(dl):: mass_sq
    integer fld
    select case(fld)
    case(1)
       mass_sq = g2*CHI2
    case(2)
       mass_sq = lambda*(3.*RECHI**2+IMCHI**2 - v2) + g2 * PHI2 
    case(3)
       mass_sq = lambda*(3.*IMCHI**2+RECHI**2 - v2) + g2 * PHI2
    end select
  end function mass_sq

!!define the condition to start the lattice simulation
  function start_box_simulation(f,p)
    logical start_box_simulation
    real(dl),dimension(ns)::f,p
    start_box_simulation = .true. !!here for Hybrid inflation we do not specify the inflationary potential, instead we use V_c as a free parameter, see arXiv: 0812.2917
  end function start_box_simulation

  function model_dt_per_step()
    real(dl) model_dt_per_step
    model_dt_per_step = metric%dx/64.d0
  end function model_dt_per_step

  function model_Power(fld,k)
    real(dl),dimension(2):: model_Power
    real(dl) k,omega
    logical,save::warning = .true.
    integer(IB) fld
    if(n*k*metric%dx .lt. const_pi) then
       model_Power =0._dl
       return
    endif
    if(k.gt.Init_Hubble)then
       omega=k**2+mass_sq(init_fields,fld)-2.25*Init_Hubble**2
       if(omega.gt.0.)then
          omega=sqrt(abs(omega))
          if(omega*metric%dx .le. const_2pi .and. omega*metric%dx*n .ge. const_2pi)then
             model_Power(1) = 0.5_dl/omega
             model_Power(2) = 0.5_dl*omega
             return
          endif
       else
          model_Power(1) = 0.5_dl/k
          model_Power(2) = 0.5_dl*k
          if(warning)then
             write(*,*) "Tachyonic region initialization may be not correct"
             warning = .false.
          endif
          return
       endif
    else
       model_Power(1) = 0.5_dl/k*(init_Hubble/k)**2
       model_Power(2) = 0._dl
       return
    endif
    model_Power = 0._dl
    return
  end function model_Power

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
       stop "you probably want to use physical time for this model"
#endif
    case(2)
       fp = open_file("data/"//trim(run_name)//"_model.info", "w")
       write(fp%unit,'(a)') "V = lambda/4 (chi^2-v^2)^2 + g^2/2 phi^2 chi^2"
       write(fp%unit,'(a,G16.7)') "reduced Planck Mass M_p = ", Mpl
       write(fp%unit,'(a,G16.7)') "lambda = ",lambda
       write(fp%unit,'(a,G16.7)') "g^2 = ", g2
       write(fp%unit,'(a,G16.7)') "V_c = ", V_c
       write(fp%unit,'(a,G16.7)') "v^2 = ", v2
       write(fp%unit,'(a,I5)') "n = ",n
       write(fp%unit,'(a,G16.7)') "dx=",metric%dx
       write(fp%unit,'(a,'//trim(int2str(ns))//'G16.7)') "Initial fields values:", Init_fields
       write(fp%unit,'(a,'//trim(int2str(ns))//'G16.7)') "Initial momenta:", Init_momenta
       write(fp%unit,'(a,G16.7)') "Initial Hubble = ",Init_Hubble
       write(fp%unit,'(a,G16.7)') "\dot\phi/H = ",Init_momenta(1)/Init_Hubble
       write(fp%unit,*) "k_min = ", const_2pi/n/metric%physdx
       write(fp%unit,*) "k_max = ", const_2pi/n/metric%physdx*fft_numk
       write(fp%unit,*) "k* = ", (2.*g2*phi_c*abs(Init_momenta(1)))**(1./3._dl)
       call close_file(fp)
       write(*,*) "k_min = ", const_2pi/n/metric%physdx
       write(*,*) "k_max = ", const_2pi/n/metric%physdx*fft_numk
       write(*,*) "k* = ", (2.*g2*phi_c*abs(Init_momenta(1)))**(1./3._dl)
    case(3)
    case(4)
    case default
       write(*,*) "Unknow option in model_output"
       stop
    end select
  end subroutine model_output

end module model
