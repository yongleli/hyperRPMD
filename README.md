# hyperRPMD
Incarnations of RPMD-related methods.

* PIMD_1D_example.ipynb: An example of 1-D PIMD with harmonic oscillator potential, 
<img src="https://render.githubusercontent.com/render/math?math=V(x) = \frac{1}{2}k(x-x_{0})^2">
<img src="https://render.githubusercontent.com/render/math?math=k = 0.5 \, \, {\rm{a.u.}}">
<img src="https://render.githubusercontent.com/render/math?math=x_{0} = 1.0 \, \, {\rm{a.u.}}">


* cayley.f90: Purely the original Cayley propagator compatible with RPMDrate.
* BCOCB.f90: Purely the BCOCB in the two papers: 


|  Papers:  |  
|---|
| arXiv:1911.00931v3 [physics.chem-ph] 15 Mar 2020  |
|  arXiv:2011.01601v1 [physics.chem-ph] 3 Nov 2020 | 

* Fortran version 1-D PIMD example: (Intel Fortran compiler needed, either FFTW or analytical Fourier transform used)
|  1-D example |  
|---|
| main.f90  |
|  system.f90 |
|  math.f90 |
|  makefile |
