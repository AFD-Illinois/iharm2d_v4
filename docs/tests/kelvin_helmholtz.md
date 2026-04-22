# Kelvin-Helmholtz instability

## Overview

A dense inner layer moves in the $+x$ direction while a lighter outer fluid moves in $-x$, creating two horizontal shear interfaces. A small sinusoidal perturbation in $v^y$ seeds the Kelvin-Helmholtz instability, which rolls up into the characteristic vortices at each interface. Unlike the linear wave tests, there is no analytic solution to compare against — this is a qualitative test of the code's ability to develop and sustain fluid instabilities.

## Setup

The domain is $[0,1]\times[0,2]$ in Minkowski coordinates, periodic on all four sides and resolved with a $256\times512$ grid by default. The initial conditions follow [Lecoanet et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4274L/abstract). The density and $x$-velocity profiles are smooth step functions built from hyperbolic tangents centered at $z_1=0.5$ and $z_2=1.5$,

$$
\rho(y) = \rho_0 + \frac{\Delta\rho}{2}\left[\tanh\frac{y-z_1}{a} - \tanh\frac{y-z_2}{a}\right],
$$

$$
\tilde{u}^x(y) = u_{\rm flow}\left[\tanh\frac{y-z_1}{a} - \tanh\frac{y-z_2}{a} - 1\right]\times\texttt{tscale},
$$

where $a =$ `a_KH` sets the shear layer width. A sinusoidal perturbation in $v^y$, localized to the two interfaces, seeds the instability,

$$
\tilde{u}^y(x,y) = A\sin(2\pi x)\left[\exp\!\left(-\frac{(y-z_1)^2}{\sigma^2}\right) + \exp\!\left(-\frac{(y-z_2)^2}{\sigma^2}\right)\right]\times\texttt{tscale}.
$$

The pressure is uniform, $\mathbf{B}=0$, and `tscale` $= 0.01$ makes the flow non-relativistic.

## Parameters

Problem-specific runtime parameters are:

| Parameter | Meaning |
|---|---|
| `tscale`  | Velocity scale; controls how relativistic the flow is |
| `rho0`    | Background density |
| `Drho`    | Density contrast between inner and outer layers |
| `P0`      | Background pressure |
| `u_flow`  | Shear velocity amplitude |
| `a_KH`    | Shear layer thickness |
| `sigma`   | Width of the Gaussian localizing the $v^y$ perturbation |
| `amp`     | Amplitude of the $v^y$ perturbation |
| `z1`, `z2`| $y$-coordinates of the two shear interfaces |

Relevant compile-time parameters are:

| Parameter | Default | Notes |
|---|---|---|
| `N1TOT`             | `256`       | |
| `N2TOT`             | `512`       | 2:1 aspect ratio matches the $[0,1]\times[0,2]$ domain |
| `METRIC`            | `MINKOWSKI` | |
| `RECONSTRUCTION`    | `WENO`      | |
| `X{1,2}{L,R}_BOUND` | `PERIODIC`  | |

## Output

The video below shows $\rho$ over the simulation at the default $256\times512$ resolution. The two interfaces roll up independently into KH vortices.

<video controls width="100%">
  <source src="../../assets/kelvin_helmholtz.mp4" type="video/mp4">
</video>

A plotting script for individual dumps is provided at `prob/kelvin_helmholtz/plot_kelvin_helmholtz.py`.

## References

- [Lecoanet et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4274L/abstract) — source of the initial condition setup used here.
