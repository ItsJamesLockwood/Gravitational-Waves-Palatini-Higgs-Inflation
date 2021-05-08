
!!******  VERSION HISTORY ************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Version 2.0 Mar 19, 2011 @ IPhT, CEA Saclay, Zhiqi Huang
!!easy-to-use version; 
!!all configuration parameters are now all defined in this file ("configure.h")
!!add HLATTICE2 discritization option (\nabla^2 accurate to (dx^3))
!!add an option to remove the gauge mode in synchronous gauge by using adaptive spatial coordinates (beta version, to be further tested).

!!Version 1.1 Feb 3, 2011 @ IPhT, CEA Saclay, Zhiqi Huang
!!add "LatticeEasy" discretization option for the case when metric perturbations are off; add option "write_check_at_step0"; add option "match_configuration_space_fluc"; add option "USE_CONFORMAL_TIME"; fix a potential bug writing null address fft_vol(0); 

!!Version 1.0 Dec 2010 @ IPhT, CEA Saclay, Zhiqi Huang
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!************* LICENCE *****************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! HLATTICE is a free software
!! You can download, use or modify the code for any non-commercial use,
!! as long as you keep this licence block here.
!! If you use the code in your research. Please cite the article
!! arXiv: 1102.0227, "The Art of Lattice and Gravity Waves from Preheating", Z. Huang (2011)
!! To report bugs or if you have questions about HLattice, contact me at
!! zhiqi.huang@cea.fr
!! or
!! zqhuang@cita.utoronto.ca
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!************ ACCKNOWLEDGEMENT ****************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!Thank Andrei Frolov who guide me to the fantastic symplectic integrator.
!!Thank my PhD supervisors J.Richard Bond and Lev Kofman who initiated this project
!!Thank IPhT@Saclay for the financial support and other facilities
!!Thank my wife Xu, who has always been my great source of motivations
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!************* import predefined constants, do not change this block***
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!#define DEBUG_MODE
#include "headfiles/predef_constants.h"
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!********** HLattice configuration section *************
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!number of scalar fields, the global variable "ns" in the code.
#define NUM_SCALAR_FIELDS 1

!!the length of each edge times Hubble: In HLattice the speed of light is defined to be 1, so this quantity is dimensionless. It should be chosen such that the dominating growing mode is captured. 
#define INIT_BOXSIZE_TIMES_H 0.100000000

!! define the discretization scheme
#define DIS_SCHEME HLATTICE2
  !!here you can use:
  !! LATTICEEASY : \partial^2_x -->  -2 ( 1 - cos(k_x) ) This option is recommended for simulations without metric perturbations. It CANNOT be used for simulations with metric perturbations, because the first-order derivatives are not defined here.
  !! HLATTICE1 : \partial_x --> i sin(k_x), faster but less accurate spatial derivatives
  !! HLATTICE2 : \partial_x --> i/3 sin(k_x) (4 - cos(k_x)), slower but more accurate spatial derivatives. This option is recommended for simulations with metric perturbations

!! define the metric
#define METRIC_OPTION FRW_BACKGROUND
  !! you can use:
  !! MINKOWSKI_BACKGROUND : minkowski spacetime
  !! FRW_BACKGROUND : FRW, evolve a(t) together with the scalar fields
  !! FRW_PERTURB_FIXED : FRW with perturbations in synchronous gauge and fixed spatial coordinates; It is not guaranteed that h_{ij} will be small for the chosen frame. The gauge modes (corresponding to the freedom of choosing arbitrary spatial coordinates in synchronous gauge) can be large enough to spoil the h_{ij}<<1 approximation. DO NOT TRUST THE RESULTS when you see h_{ij} ~ 0.01: in such case the code might still be stable, but the approximations in HLattice is no longer valid.
 !! FRW_PERTURB_ADAPTIVE: use adaptive spatial coordinates to remove gauge modes. This option can slightly spoil energy conservation (since you are now summing up the energy on different physical positions). This is a new feature in HLattice. The stability need to be further studied.


!! define the integrator: SYMPLECTIC_2ND, SYMPLECTIC_4TH, SYMPLECTIC_6TH
#define INTEGRATOR SYMPLECTIC_6TH

!! define the # of grid points along each edge of the cubical box (global variable "n" in the code):it must be integer power of 2 (and between 4 and 4096); the total number of grid points is n^3
#define SIMU_RESOLUTION 512

!! want gravitational waves? You can use it for any METRIC_OPTION. In the case METRIC_OPTION = MINKOWSKI_BACKGROUND or FRW_BACKGROUND, the metric perturbations will be integrated without giving feedback to the scalar fields.
#define WANTGW YES

#if METRIC_OPTION ==  FRW_BACKGROUND
!! use conformal time? This only works for METRIC_OPTION == FRW_BACKGROUND
#define USE_CONFORMAL_TIME YES

#else
#define USE_CONFORMAL_TIME NO
#endif

!!stop the program if the scale factor exceeds the following value
!!Note that the scale factor is normalized to be 1 at the beginning of lattice simulation
#define MAX_SCALE_FACTOR 10.

!!Determine which y-z slice of the lattice should be output for the fields and momenta
!!By default, set to the first slice
#define WHICHSLICE XZ 
!! Choose between 0 and n.
#define YZ_SLICE 32
#define XY_SLICE 1
#define XZ_SLICE 1

#define WANTFIELDS YES
#define WANTSLICES NO
!! As slices and fields are resource intensive, there is an option to save them less regularly than regular checkpoints.
#define CHECKPOINTS_PER_SLICE 5
#define CHECKPOINTS_PER_FIELD 10
#define N_FIELD_SAVES 50 
!!*** below are some settings that you can, but usually do not need to change *****

!!how often to write feedback to the screen (# of steps, recommended value is between 5 and 50)
!!IMPORTANT: If you are calculating gravity waves in the mode METRIC_OPTION = FRW_BACKGROUND or MINKOWSKI_BACKGROUND, h_{ij} will be evolved only when the feedback is written. In such case you should not make this number too big (5-10 is recommended). 
#define FEEDBACK_INTERVAL 4
#define LINES_BETWEEN_SAVES 10

!!how often to write the checkpoint files and gravity wave spectrum. It must be an interger times FEEDBACK_INTERVAL
!! you can make it a huge number (e.g. 10000000*FEEDBACK_INTERVAL) to disable checkpoints
!!IMPORTANT: for METRIC_OPTION = FRW_PERTURB_ADAPTIVE, the redefinition of spatial coordinates (i.e. gauge transformation) is done when the checkpoint files are written. So do not make it too small in this case.
#define CHECKPOINT_AND_GW_INTERVAL (FEEDBACK_INTERVAL * LINES_BETWEEN_SAVES)

!!stop the program if # of steps exceeds
#define MAX_STEPS 10000000

!!In HLattice if metric feedbacks are turned on, I split the noncanonical term e^{Adt} to (e^{A/m dt})^m and integrate e^{A/m dt} using a 4-th order Runge-Kutta integrator. Typically m=10 suffices. But feel free to change it here.
#define SUB_INTEGRATOR_STEPS 10

!!Reduced Planck Mass, it can be any number. To achieve better performance you may want to define it ~ (Planck Mass / Energy Scale of your problem), so the numbers in the simulation can have order ~ 1. (Since HLattice uses double precision numbers, this usually is not a big issue.)
#define REDUCED_PLANCK_MASS 1.00000000000

!!use k^{std} instead of k^{eff} in the TT projector;
#define USE_STANDARD_WAVENUMBER NO

!!only print the scale factor and Hubble parameter?
#define FEEDBACK_ONLYAH NO
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!**************** end of configuration section *********************!!

!!******************* include the definitions of spatial derivatives, do not change this block *****************
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "headfiles/def_derivs.h"
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
