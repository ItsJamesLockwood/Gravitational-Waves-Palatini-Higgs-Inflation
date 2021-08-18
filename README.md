# Master's Project: Graviational Waves produced during the non-linear regime of Palatini Higgs inflation.

We study the model of Higgs inflation with a non-minimal coupling to gravity in the Palatini formalism and obtain gravitational wave (GW) spectra using lattice simulations. Perturbations in the Higgs field are rapidly amplified due to the tachyonic instability of the potential. As a result, non-linear dynamics quickly become dominant and a full lattice simulation is required to study the evolution. We used _HLattice_ to obtain spectra of GW produced during preheating. The spectra saturate at <img src="https://render.githubusercontent.com/render/math?math=\Omega_\textrm{gw}H^2_{100}\sim10^{-11}">, which is similar to the other inflationary models, namely the hybrid model with a Higgs-like potential. A study of lattice cross-sections revealed that large-scale 3D-structure appears due the preferred amplification of highly tachyonic modes. The formation of highly perturbed regions of the scalar field were linked to the classical production of GW.

This repository contains two distinct, but related, programs:
- the lattice simulation software _HLattice_ (Fortran 90),
- the data analysis tools to process the outputs from the simulations.

The former has been extensively modified from the original simulation software, last update in 2012 (see the [original HLattice by Z.Huang](https://www.cita.utoronto.ca/~zqhuang/hlat/)).
New features include detailed simulation output of gravitational energy density throughout the lattice, tracking the formation of pseudo-solitonic objects (i.e. _bubbles_) in the lattice, and fixing existing issues caused by improper multithreading. 

The Python modules provide tools for formatting, analysing and saving the results from the Fortran 90 _HLattice_ software. This includes, among others:
- dashboard for visualising the results from an entire simualtion: field values, energy densities, energy contained in field perturbations and/or in gravitational waves.
- animation tools for playing back the evolution of pseudo-solitonic objects inside the lattice, showing how the initial uniform Higgs field quickly fragments.
- analytics functions for quantifying the behaviour of the Higgs-Palatini model and comparing the simulation results to previously studied models.

Some examples of the dashboards and output from the Python analysis scripts are shown below.

Output from the dashboard, giving an overview of a typical simulation:
![rtanh_std_metric_h_panel](https://user-images.githubusercontent.com/33159939/129893982-10fe5b5c-48aa-4663-893a-a4a783df81be.png)

Overlay of the Floquet exponent for the field's evolution at different scale factors _a(t)_:
![image](https://user-images.githubusercontent.com/33159939/129894952-90ee34e8-3d7d-4e45-995f-4a76f504a95e.png)

Comparison of the evolution of a perturbation mode in the approximated linear regime and the simulated fully non-linear regime:  
![image](https://user-images.githubusercontent.com/33159939/129895057-9d883d32-ba6c-4b71-a79e-0cc364ac62a7.png)

Slice of the lattice showing the gravitational wave energy density late in the simulation:
![image](https://user-images.githubusercontent.com/33159939/129895208-75888da8-2c81-489b-87c2-520d457342a3.png)


Please note, that due to the ever-evolving nature of the requirements of this project, the Python scripts are not heavily documented at this time and require a reorganisation before they can be used effectively in the public domain. For any questions regarding this project, including its current status, please feel free to get in touch. 
