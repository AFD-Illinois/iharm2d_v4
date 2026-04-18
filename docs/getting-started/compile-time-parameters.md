# Compile-time parameters

As mentioned in the Quickstart guide, compile-time parameters are baked into the executable at build time — changing any of them requires a fresh `make`. They live in `parameters.h` inside your `build_archive/` directory (copied there from `prob/PROB/parameters.h` when you first compile), so you can edit them freely without touching the original repository.

The torus problem (`prob/torus/parameters.h`) is a good representative example, and we will use it throughout this page.

## Grid resolution

```c
#define N1TOT 256
#define N2TOT 256
```

`N1TOT` and `N2TOT` set the total number of grid zones in the $x^1$ and $x^2$ coordinate directions, respectively. In Kerr-Schild coordinates, as used in the black hole accretion setup, these directions map to functions of the radial coordinate $r$ and polar angle $\theta$. Unlike [`iharm3d`](git@github.com:AFD-Illinois/iharm3d.git), there is no companion `NiCPU` parameter here — `iharm2d_v4` does not use MPI, so the `NiTOT` values fully describe the grid with no further decomposition.

## Metric and coordinate system

```c
#define METRIC MKS
#define DEREFINE_POLES 1
```

`METRIC` defines the spacetime metric used in the simulation. The two supported values are:

- `MINKOWSKI` — flat spacetime in Cartesian coordinates, appropriate for many test problems.
- `MKS` — Modified Kerr-Schild coordinates, employing an exponential radial coordinate and a modified polar coordinate that increases the resolution close to the midplane. The exact transformation is given by Equation (F1) in [the PATOKA paper](https://ui.adsabs.harvard.edu/abs/2022ApJS..259...64W/abstract).
- 
`DEREFINE_POLES` takes a value of `0` or `1`. Setting it to `1` activates FMKS (Funky MKS) coordinates, which enlarge zones in the $\theta$ direction near the polar axis and close to the event horizon. This relaxes the severe timestep restriction caused by very small cells, whose light-crossing times become extremely short. For torus problems, leaving this at `1` is strongly recommended.

## Floors

```c
#define WIND_TERM    0
#define BSQORHOMAX   100.
#define UORHOMAX     100.
```

Floors prevent the fluid state from reaching unphysical values in low-density regions. `BSQORHOMAX` and `UORHOMAX` cap the maximum allowed ratios of magnetic energy density to rest-mass density (plasma magnetization) and internal energy density to rest-mass density (proportional to temperature), respectively. The defaults of `100.` are appropriate for most torus runs; you may want to raise them if your problem has strongly magnetised jets or lower them to keep the floors from activating too aggressively. Additional floors include `BSQOUMAX` (which caps the inverse plasma beta), geometric floors (`RHOMIN` and `UUMIN`) that inject rest-mass density and internal energy as functions of radius, and a velocity ceiling (`GAMMAMAX`). If these are not set in `parameters.h`, their default values are taken from `decs.h`; see that file for the complete list. `WIND_TERM` enables a small mass-injection source term (`1`) that can help stabilise the evacuated funnel region in torus problems; it is off (`0`) by default.

## Electron thermodynamics

```c
#define ELECTRONS           0
#define ALLMODELS           0
#define SUPPRESS_HIGHB_HEAT 1
```

Setting `ELECTRONS` to `1` activates a separate electron entropy equation, enabling two-temperature plasma modelling. This is based on the entropy-tracking procedure of [Ressler et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.454.1848R/abstract). `ALLMODELS` (`0` or `1`) tracks separate electron entropies for all built-in heating prescriptions simultaneously, which avoids re-running the GRMHD simulation for each heating model. `SUPPRESS_HIGHB_HEAT` (`0` or `1`) disables electron heating in strongly magnetised regions (magnetisation $\sigma > 1$) when set to `1`. NOTE: The electron heating prescription is applicable to the black hole accretion problem for which it was devised.

## Reconstruction algorithm

```c
#define RECONSTRUCTION WENO
```

Controls the spatial interpolation scheme used to interpolate the primitive variables to zone faces for flux calculation. The two options are:

- `LINEAR` — piecewise-linear reconstruction with monotonized central (MC) slope limiter; faster but more diffusive.
- `WENO` — fifth-order weighted essentially non-oscillatory reconstruction; more accurate and the recommended default.

## Boundary conditions

```c
#define X1L_BOUND OUTFLOW
#define X1R_BOUND OUTFLOW
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR

#define X1L_INFLOW 0
#define X1R_INFLOW 0
```

The four `_BOUND` parameters set the boundary condition at each face of the domain. Supported values are:

- `OUTFLOW` — outflow boundary condition; copy the outermost active zone into all ghost zones. Material flows out freely of the domain. Standard for the radial (`X1`) boundaries in a black hole accretion simulation.
- `POLAR` — reflective polar boundary; reflect the fluid state across the polar axis. The sign of `U2` and `B2` is flipped across the polar boundary. Required for the `X2` boundaries when using MKS coordinates.
- `PERIODIC` — periodic boundary, used in most flat-space test problems.

`X1L_INFLOW` and `X1R_INFLOW` (`0` or `1`) control whether inflow through the corresponding `OUTFLOW` boundary is permitted. Setting either to `0` prevents matter from being pulled back into the domain if the velocity at that face reverses — a common safeguard at the inner radial boundary near the black hole horizon.
