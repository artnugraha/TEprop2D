# TEprop2D
Simple Fortran program to calculate thermoelectric properties of 2D materials by taking the outputs of Quantum ESPRESSO and EPW packages as the input of the program.

The code is "TEprop2D.f90", which can be compiled by standard Fortran compiler (e.g., ifort or gfortran):

`gfortran TEprop2D.f90 -o TEprop2D.out`

Execution:

`./TEprop2D.out`

## Input format

The program will read "TEprop.inp" as input. It then calls two other important input files: energy dispersion and self energy. The format of "TEprop.inp" should strictly follow this order:
```
energy dispersion filename
self energy filename
temperature
alatt (lattice constant in angstrom)
thick (thickness in angstrom)
nkx (# of kx point)
nky (# of ky point
nbnd (# of bands)
nbndvb (index of top valence band)
nbndcb (index of bottom conduction band)
nbndc (# of bands to be calculated)
Efermi0, Efermif (E_F plot range)
kappa (# lattice thermal conductivity)
```

See examples of "TEprop.inp", "band.eig-300K-200x200-0.00" (energy dispersion file), and "linewidth.elself-300K-200x200-0.00" (self energy file, please unzip it first) to better understand the format.

## Outputs

After a successful execution, we will get several output files, mainly:
- Carrier density ($10^21$ cm$^{-3}$) vs. Fermi energy (eV)
- Seebeck coefficient (mV/K) vs. Fermi energy (eV)
- Electrical conductivity ($10^8$ / ohm.m) vs. Fermi energy (eV)
- Power factor (W/K$^2$ m ) vs. Fermi energy (eV)
- Electronic thermal conductivity (W/(m K)) vs. Fermi energy (eV)
- ZT vs. Fermi energy (eV)
All thermoelectric properties contain longitudinal x- and y-directions. There are also corresponding outputs of thermoelectric properties as a function of carrier density (instead of a function of Fermi energy).

## Citations
If readers benefit from this repository, please kindly cite one (or all) of the following papers, which used the program with some adjustment for different 2D materials:
- <strong>Designing high-performance thermoelectrics in two-dimensional tetradymites</strong>, <a href="https://doi.org/10.1016/j.nanoen.2019.02.015">Nano Energy 58, 743-749 (2019)</a>.
- <strong>Thermoelectric performance of monolayer InSe improved by convergence of multivalley bands</strong>, <a href="https://aip.scitation.org/doi/10.1063/1.5040752">Journal of Applied Physics 125, 082502 (2019)</a>.
- <strong>Two-dimensional InSe as a potential thermoelectric material</strong>, <a href="https://aip.scitation.org/doi/10.1063/1.5001184">Applied Physics Letters 111, 092107 (2017)</a>.

## Contributors:
- Ahmad Ridwan Tresna Nugraha, <a href="mailto:art.nugraha@gmail.com">art.nugraha [at] gmail [dot] com</a>
- Nguyen Tuan Hung, <a href="mailto:hung@live.jp">hung [at] live [dot] jp</a>
