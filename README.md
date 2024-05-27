# 1-D-Supernova-shock-radius-evolution
The code I built is based upon the !["ZEUS-2D"](https://ui.adsabs.harvard.edu/abs/1992ApJS...80..753S/abstract) paper which solves the Hydrodynamics equations in 2 dimensions, I used a one dimensional version of it to reproduce the evolution of a Supernova explosion which is assumed to be spherically symmetric. The numerical approach taken is discussed in the following chapter.

# Hydro Equations

The goal of the project is, in its essence, to solve the following set of equations for a specific astrophysical problem: 

$$\frac{D\rho}{Dt}     =-\rho \nabla \cdot v $$

$$\rho \frac{Dv}{Dt}   =-\nabla p  $$

$$\frac{D\epsilon}{Dt} =-p\nabla \cdot v $$


