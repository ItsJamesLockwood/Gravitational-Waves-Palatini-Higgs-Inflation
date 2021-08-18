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

Please note, that due to the ever-evolving nature of the requirements of this project, the Python scripts are not heavily documented at this time and require a reorganisation before they can be used effectively in the public domain. For any questions regarding this project, including its current status, please feel free to get in touch. 
