# 1-D-Supernova-shock-radius-evolution
The code I built is based upon the !["ZEUS-2D"](https://ui.adsabs.harvard.edu/abs/1992ApJS...80..753S/abstract) paper which solves the Hydrodynamics equations in 2 dimensions, I used a one dimensional version of it to reproduce the evolution of a Supernova explosion which is assumed to be spherically symmetric. The numerical approach taken is discussed in the following chapter.

## Hydro Equations

The goal of the project is, in its essence, to solve the following set of equations for a specific astrophysical problem: 

$$\frac{D\rho}{Dt}     =-\rho \nabla \cdot v $$

$$\rho \frac{Dv}{Dt}   =-\nabla p  $$

$$\frac{D\epsilon}{Dt} =-p\nabla \cdot v $$

When integrating the respective lagrangian derivatives,

$$ \frac{D}{Dt} \equiv \frac{\partial}{\partial t} + v\cdot \nabla$$

these three equations define the conservation laws for $\rho$ (density), v (velocity) and $\epsilon$ (energy) respectively. p (pressure) is determined by an equation of state, which in this case is set to be

$$p=(\gamma- 1)\epsilon$$



