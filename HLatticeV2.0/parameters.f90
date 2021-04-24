module parameters
  implicit none
!!since HLattice V2.0, all parameters are set in "configure.h", do not change this file anymore.

!!============================================================================
#include "configure.h"

  integer,parameter :: dl=kind(1.d0) !!use double precision
  integer,parameter :: dlc=kind((1.d0,1.d0)) !!double precision complex numbers
  integer,parameter :: IB = 4 !! use 4-bytes integers

  real(dl),parameter :: const_pi =   3.141592653589793238462643383_dl
  real(dl),parameter :: const_2pi =6.283185307179586476925286766_dl
  real(dl),parameter :: const_pi2 = const_pi**2
  real(dl),parameter :: const_sqrt2 = 1.4142135623730950488_dl
  real(dl),parameter :: const_sqrt3 = 1.7320508075688772935_dl
  real(dl),parameter :: const_sqrtpi = 1.7724538509055160273_dl
!!==============================================================================
#if SIMU_RESOLUTION == 16
#include "headfiles/n_params_4.h"
#elif SIMU_RESOLUTION == 32
#include "headfiles/n_params_5.h"
#elif SIMU_RESOLUTION == 64
#include "headfiles/n_params_6.h"
#elif SIMU_RESOLUTION == 128
#include "headfiles/n_params_7.h"
#elif SIMU_RESOLUTION == 256
#include "headfiles/n_params_8.h"
#elif SIMU_RESOLUTION == 512
#include "headfiles/n_params_9.h"
#elif SIMU_RESOLUTION == 1024
#include "headfiles/n_params_10.h"
#elif SIMU_RESOLUTION == 2048
#include "headfiles/n_params_11.h"
#elif SIMU_RESOLUTION == 4096
#include "headfiles/n_params_12.h"
#endif


!!======================= # of scalar fields ==================================
  integer(IB),parameter :: ns = NUM_SCALAR_FIELDS
!!======================= box size ============================================
  real(dl),parameter::boxsize_H = INIT_BOXSIZE_TIMES_H 
!!======================= feedback ==============================================
  integer(IB),parameter::n_feedback = FEEDBACK_INTERVAL 
!!======================== stop simulation =======================================
  integer(IB),parameter:: stop_at_step = MAX_STEPS 
  real(dl),parameter:: stop_at_a = MAX_SCALE_FACTOR  
!!======================== check point ============================================
  logical,parameter:: use_checkpoint = .true.
  logical,parameter::write_check_at_step0 = .true. 
  integer,parameter:: checkpoint_steps = CHECKPOINT_AND_GW_INTERVAL


!!======================== refined steps for Runge-Kutta sub-integrator ======
  integer(IB),parameter::noncanonical_Runge_Kutta_steps = SUB_INTEGRATOR_STEPS

!!============================= reduced Planck Mass M_p=======================
  real(dl),parameter :: Mpl = REDUCED_PLANCK_MASS !!1024._dl 

!!============================ default prefix of output files ================
!!the output files will be in ./data/ directory. All named with a prefix.
!!the prefix can be specified by running HLattice in command line:
!! ./HLattice prefix
!!when you run HLattice without passing a prefix to it, a default prefix will be used. It's defined here:
  character(LEN=256)::run_name="test" !! *

!!============================ output of field/momenta slices =====================================
  !! Probably a bit redundant, but section available for future use 
  integer(IB),parameter :: save_slice_interval = CHECKPOINTS_PER_SLICE
  integer(IB),parameter :: which_lattice_cut = WHICHSLICE
  
  integer(IB),parameter :: yz_lattice_slice = YZ_SLICE 
  integer(IB),parameter :: xy_lattice_slice = XY_SLICE 
  integer(IB),parameter :: xz_lattice_slice = XZ_SLICE 

  integer(IB),parameter :: save_field_interval = CHECKPOINTS_PER_FIELD
  integer(IB),parameter :: field_number_cutoff = N_FIELD_SAVES
  

!!============================ do not change things below =========================================
!! n/2
  integer,parameter::Nby2=n/2
!!n^2
  real(dl),parameter:: nsquare = (1._dl*n)**2
!!n^3
  real(dl),parameter:: ncube = (1._dl * n)**3
!!wave number unit * dx
  real(dl),parameter::k_unit = const_2pi/n

  real(dl),parameter :: PlanckMass = 2._dl *const_sqrt2 * const_sqrtpi * Mpl !!the Planck Mass (not reduced)
  real(dl),parameter :: Mplsq = Mpl**2 !!reduced Planck Mass squared
  real(dl),parameter :: Newton_G = 1./Mplsq/8./const_pi
  real(dl),parameter :: GeV = Mpl/2.43e18_dl
  
  integer(IB),parameter::n_Hamiltonian_terms = 3 ! number of non-comutable terms in Hamiltonian of the system

!!organize file units
  logical::unit_opened(7:99) = .false.


!!a global variable, will be calculated automatically. No need to change it here.
  real(dl):: init_Hubble 
!!========================================

end module parameters
