This is HLattice V2.0 by Zhiqi Huang (zhiqi.huang@cea.fr or zqhuang@cita.utoronto.ca). This code simulates one or a few canonical scalar fields on a cubical lattice. You can choose to include metric perturbations or not. You can choose between 2nd, 4th, or 6th symplectic integrators.

Version 2.0 is a user-friendly version. I assume you are in a linux environment. And you have the latest free fortran compiler gfortran with openmp enabled (if openmp is not enabled, you can switch off the openmp option in Makefile) and the free gcc C preprocessor cpp installed. If that is true, here are four simple steps to use HLattice:

1. Read and modify the self-explanary file "configure.h".

2. Read and modify the self-explanary file "model.f90"

3. Run "make" in a terminal

4. Run the code "./HLattice run_name" in a terminal. 
Here run_name will be any string that can be used as a prefix of file names. The default run_name is "test". So if you just run "./HLattice", the output files will be named as "./data/test_...".


The outputs: I assume you use a run_name "prefix".
./data/prefix_model.info
The information about this model. You can define how to write this file in the subroutien model_output in "model.f90".

./data/prefix_screen.log
The information about background evolution and the rms fluctuations of the fields (and h_{ij} if metric perturbations are included).

./data/prefix_pw_i.log
here i can be 1, 2, ..., NUM_SCALAR_FIELDS
let's assume f is the i-th scalar field. The outputs are:
line#1: scale factor a
line#2: The power spectrum of the scalar fields k^3 [(k/a)^2|f_k|^2] /(4pi^2 \rho_{background}), for k=2pi/L, 4pi/L, .... (L is the length of each edge) 
line#3: The power spectra of the scalar fields k^3 [|\dot f_k|^2] /(4pi^2 \rho_{background}), for k=2pi/L, 4pi/L, .... 
Such 3 lines are repeatly written as the fields evolves.

./data/prefix_GW.log
the gravitational wave spectrum. The outputs are converted to today's observables using the formulas given in  arxiv 0812.2917.
line#1: scale factor a
line#2: list of frenquencies in unit of Hz
line#3: today's d (\Omega_{gw}h^2)/d(\ln frequency)

The other files are checkpoint binary files. When you run HLattice with the same run_name, these files will be loaded and the simulation is resumed.



**************************************************************************
More details about what the code actually does:

HLattice evolve the following quantities as functions of time:

fields_f(:,:,:,:), fields_p(:,:,:,:)

these variables are defined in "define_fields.f90".

The physical meaning of them are the physical field values (\phi) and physical field momenta (d\phi/dt a^3) (when metric perturbation back-reaction is switched on, the definition of field momenta is a bit more complicated, I can explain to you if you are interested).

For example, fields_f(2,3,4,5) is the value of the 2nd field at grid point (x=2, y=4, z=5). Of course all of these variables evolves as functions of time.

When you launch the code, it does NOT evolve fields_f and fields_h immediately. The code first start a background evolution. The initial BACKGROUND field values (phi_1, phi_2, ...) are defined by you in model.f90 (the init_fields variable at the beginning of the file).
The initial field momenta are defined either by "init_momenta" variable also given by you in "model.f90", right below the "init_fields". However, if you set init_momenta = Mplsq (the reduced Planck Mass^2, which of course cannot be a correct value), the code will calculate the initial BACKGROUND field momenta using the slow-roll equation:
  3H\dot\phi \approx - \partial V / \partial \phi

So far we are just talking about a background evolution, essentially just solving the dynamics of inflation (not preheating yet!).

So when does the code start a lattice simulation? This is defined by the function "start_box_simulation(f,p)", here array f is the background field values, array p is the time derivative of the background field values. The code evaluate this function every 0.005 Hubble time. When this function returns .true., the code stops the background evolution, and start to launch the lattice simulation. The current field values and time derivatives are put into lattice as the mean field values and mean field momenta. On top of it you also need initial field fluctuations and field-momenta fluctuations, this is realized by Fourier transforming a random Gaussian fields with power spectra |\phi_k|^2 and |\dot\phi_k|^2 (user-defined, function model_Power in model.f90; the default is use n_k=1/2).
