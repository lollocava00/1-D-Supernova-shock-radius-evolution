# Introduction to Zeus-2D code and relevant equations
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

# Introduction to the Astrophysical problem
On average, a total energy of about $10^{51} erg$ is injected into the interstellar medium (ISM) by supernovae.
The resulting shock wave propagates outwards at very high speed, but quickly loses large amounts of energy due to compression and heating of the ISM and also because of radiative losses. To allow the modelling of this rather complex process, some simplifying assumptions were made, which include radial symmetry of the entire system, a uniform ISM, the absence of a stellar wind phase before the explosion and thermal conduction is neglected. 

The code is set up to allow for the analysis of an expanding SuperNova Remnant (SNR), more specifically in the following terms:
- evolution of $\rho$,$p$,$T$ and $v$ over various time ranges;
- evolution of the shock radius generated by the supernova event in radial expansion;
- evolution of the X-RAY luminosity of the source as a function of time;
- the kinetic and thermal energy fractions injected into the surrounding ISM as a function of time.

The aforementioned characteristics are evaluated both in the presence and lack of radiative cooling as to highlight possible differences. 
The initial conditions chosen for this setting are typical values in the ISM:

$$d_0=2\cdot 10^{-24} g/cm^3 $$

$$ T_0=10^4 K $$

$$e_0=c_v d_0 T_0 \quad where \quad c_v=2\cdot 10^8 \, cgs $$

$$ p_0=e_0(\gamma -1) \quad where \quad  \gamma=5/3$$

The central energy values are set to $E_0 = 10^{51} \, erg$ in order to simulate a realistic supernova output; this allows us to define the central energy density as

   $$ e(2)=e(3)=\frac{E_0}{\frac{4}{3}\pi x_a^3} $$

The time step over which the DO loop runs is set to:

 $$   \Delta t=cfl *\frac{\Delta x}{|v|+c_s}$$

where the $cfl$ value is increased from $0.01$ to $0.5$ along with the time integration, in the code:

  $$  clf=min(0.5\, , \, clf*1.1).$$

Such values for $\Delta t$ always satisfy the stability condition of Upwind I order method; 

  $$  \Delta t \leq \frac{\Delta x}{v}$$

## Results without cooling
The DO loop is performed over a selection of time spans:

  $$ t= 2\cdot 10^4,4\cdot 10^4,6\cdot 10^4,8\cdot 10^4, 10^5, 2\cdot 10^5,3\cdot 10^5, 4\cdot 10^5,5\cdot 10^5 \, yrs $$

at first, the contribution of radiative cooling is neglected.
The first set of plots shows the evolution of the different quantities between $2\cdot 10^4$ and $10^5$years, while figure 2 only
focuses on the time span that goes from $10^5$ to $5\cdot 10^5$ years. The scale of the variables on both axis except for the velocity is logarithmic.

In this case a smooth expansion of the shock wave is observed.
As expected, the points of maximum shock wave pressure are gradually shifted to the outer regions of the cluster. These peaks are also characterised by a decreasing intensity as the shock loses energy by heating the ISM, causing it to slow down. The same applies to the maximum density; the temperature also decreases over time as the shock expands and the injected energy is distributed over a greater volume.

![Alt text](plots/hydro_new.png?raw=true)

![Alt text](plots/hydro_new2.png?raw=true)

Each quantity is shown to vary drastically when crossing the shock radius. Density and velocity increase steadily
up until they reach the point of discontinuity in the fluid, after which their values decrease; pressure remains constant up until very close to the shock radius, while temperature decreases steadily
and drops right after $R_{shock}$.

The same is true for figure 3, where all the quantities keep evolving following the same trends.

## Results with cooling
The contribution of a cooling function $\Lambda(T)$ is added; $\Lambda(T)$ defines the bolometric power emitted by the surrounding
medium in the form of thermal radiation and takes into account different properties of the plasma. The function I used was first introduced by Sutherland and Dupita:

$$\Lambda(T)=10^{-22}(8.6\cdot 10^{-3}T_{KeV}^{-1.7}+0.058T_{KeV}^{0.5}+0.063 \quad for \quad T>0.02 \,KeV $$

$$\Lambda(T)=6.72\cdot 10^{-22}(T_{KeV}/0.02)^{0.6} \quad for \quad 0.0017235 \leq T \leq  0.02 \;KeV$$

$$\Lambda(T)=1.544\cdot 10^{-22}(T_{KeV}/0.0017235)^6 \quad for \quad T<0.0017235 \, KeV$$

Density, pressure, velocity and temperature evolution is computed the same way as for the non-cooling case and
shown in figures 4 and 5. The behaviour of these quantities is visibly affected by the presence of a cooling function.
Densities do not remain constant in time, but rather decrease right behind the shock radius and increase at the discontinuity; velocities appear to show a slight increase where the previous case featured a constant decrease; pressure also shows an uneven, faster decrease before the shock radius due to the energy dissipation caused by the cooling function; temperature instead features a small peak right after the discontinuity starting from $ 6\cdot10^4 \,yrs$, which was absent in the no-cooling case.

![Alt text](plots/hydro_new3.png?raw=true)

![Alt text](plots/hydro_new4.png?raw=true)

Another observation that can be made with regards to the density is that, after the increase at the shock surface featured in Figure 4, the values stop increasing after around $10^5 \,yrs$, but the densities behind the shock keep decreasing.
%another observation that can be made with regards to the temperature is that after about $2\cdot 10^4 yrs$ the shell temperature (look at shock radius) is not decreasing, this signifies the end of the adiabatic phase and the beginning of the radiative phase in which the energy that is gained by impact with the ISM is instantly radiated away.

## Shock radius
The shock radius describes the distance from the centre of the explosion to the ”front” of the shock wave. An analytical expression exists, called themSedov solution. It allows the calculation of the shock radius via

 $$ R_s=\left( \frac{2E_0}{\rho_0} \right)^{1/5}t^{2/5}.$$

However, this solution is only relevant in the adiabatic phase of the SNR expansion and should stop working after radiative losses become important.
Numerically the front of the shock wave can not be defined in a rigorous way. In order to obtain a fixed value for the shock radius it was chosen as the point in space at which the density is maximal.

The code computes the value of the shock radius every $10^3$ years by introducing an IF statement near the
end of the DO loop; the right time step is selected by using a counter that only contains multiples of $10^3$.
Inside the IF statement the value of the Shock radius is taken as the grid point that contains the maximum value of the density. This is done by using the function \textit{maxval(d)}.

As the literature on Supernovae tells us, radiative losses start to beacome important after a few $10^4 yrs$ strongly depending on the ISM density. Figure 7 shows that in our case this happens at about $log(t) \sim 4.5$ or $3\cdot 10^4 yrs$ .


![Alt text](plots/R_shock1.png?raw=true)

The shock radius can be computed following the same procedure as the no-cooling case and compared to
the previous linear fit together with the analytical solution.
The shock radius stops following the Sedov law due to energy dissipation
via cooling and a subsequent loss of driving momentum to the shock wave. The results of the shock radius computation
are plotted in Figure 7.

![Alt text](plots/R_shock2.png?raw=true)

## X-Ray luminosity
The supernova remnant is a source of X-RAY radiation, mainly caused by bremsstrahlung emission, recombination continuum and two-photon emission; the values initially increase, but as the energy of the charged particles is irradiated through high energy photons over time the luminosity decreases. The following plot shows the evolution of the source’s X-RAY luminosity both in the presence of a cooling function and in the absence of it. One can notice how the more efficient energy dissipation caused by $\Lambda$ causes a steeper decrease over time.


![Alt text](plots/LX1.png?raw=true)

The time steps over which the luminosity values are computed are found by following the same procedure used to compute the shock radius, thereby generating a value of $L_x$ every $10^3 yrs$.
The X-ray luminosity is analitically obtained by computing the following integral:

 $$L_x=\int_{V_{SNR}}n^2\Lambda(T)4\pi r^2 dr$$

which i calculate numerically as a sum over the volume shells $dvl_{1a}(i) $:

$$L_x=L_x+4\pi dvl_{1a}(i) \cdot \Lambda(T(i))$$

This is computed only for the grid point which have $T>10^6 \, K$ as the gas only emits via bremsstrahlung in the X-ray for such temperatures.

## Energy conservation
As per the last request I found the evolution of the kinetic, thermal and total energy contents over time as

$$   E_{th}=E_{th}+ 4\pi e(i) dvl_{1a}(i) - 4\pi dvl_{1a}(i) cv \cdot v \rho_0 T_0 $$
$$    E_{kin}=E_{kin}+\frac{1}{2}\cdot 4\pi d(i) dvl_{1a}(i) v(i)^2 $$
$$    E_{tot}=E_{th}+E_{kin}$$

The objective of this section is to show whether the energy is conserved in this system. From the following plot one can deduce that, when no cooling function is in action (adiabatic case), the energy content of the system is somewhat conserved, despite a slight decrease, while the changes are more noticeable and abrupt after a certain time step close to $log(t) \sim 4.5$ when a dissipative $\Lambda$ is present; This is consistent with the previous shock radius and x-ray luminosity results in which radiative losses become important around that time. Figure 9 plots the different energy contents in units of $10^{51}\, erg$.
%they also had the initial energy of the supernova ($e_{in}$) subtracted from them.

![Alt text](plots/Energy_1.png?raw=true)

![Alt text](plots/Energy_cooling.png?raw=true)



