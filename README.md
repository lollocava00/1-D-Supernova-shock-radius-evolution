# Introduction to Zeus-2D codee and relevant equations
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

where $\gamma$ is set to the value $\gamma=\frac{5}{3}$.

Equations 1-3 constitute a set of hyperbolic partial differential equations, which are tackled in the code through a finite-differences approach (upwind I); in particular, the code breaks the solving procedure in two separate steps:
- Source step, where only the soure terms of the equations are considered:

$$\rho \frac{Dv}{Dt} =-\nabla p -\nabla Q $$

 $$       \frac{D\epsilon}{Dt}=-p\nabla \cdot v -Q \nabla \cdot v$$

 where Q is an \textbf{artificial viscosity term} that will be discussed later; 
    
  - Transport step, which accounts for the effects of fluid advection:

$$ \frac{\partial }{\partial t} \int_V \rho dV =-\int_{\partial V} \rho v \;dS $$
    
$$ \frac{\partial}{\partial t} \int_V \rho v dV =-\int_{\partial V} \rho v^2 \; dS $$
    
$$\frac{\partial}{\partial t} \int_V \epsilon dV =-\int_{\partial V} \epsilon v \;dS $$
   
    (the grid velocity is assumed to be null for our purposes).
The aforementioned equations are solved on two staggered grids,  named $x_a$ and $x_b$, where $x_b(i)=x_a(i+\frac{1}{2})$. A choice is also available between two differente coordinate systems (cartesian or spherical).

![Alt text](plots/shocktube.png?raw=true)



