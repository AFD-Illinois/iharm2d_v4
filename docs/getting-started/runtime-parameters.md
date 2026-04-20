# Runtime parameters

Runtime parameters are read from `param.dat` each time the executable is launched — you can change them between runs without recompiling. Each problem ships with a sample `param.dat` under `prob/PROB/param.dat`; copy it into your output directory before running.

The parameter file format is straightforward: each line begins with a type tag (`[dbl]`, `[int]`, or `[str]`) that tells the parser which data type to read, followed by the parameter name and its value. Lines beginning with `#` are ignored.

The torus problem (`prob/torus/param.dat`) is used as the running example throughout this page.

## Time control

```
[dbl] tf     = 1000.0
[dbl] dt     = 1.0e-06
[dbl] cour   = 0.9
```

`tf` is the simulation end time in code units (gravitational units $GM/c^3$ for black hole problems). `dt` sets the initial timestep; the code will adapt it quickly, so this value mainly prevents an enormous first step — `1e-6` is a safe conservative starting point. `cour` is the Courant–Friedrichs–Lewy (CFL) number: the timestep is set to `cour` times the minimum light-crossing time across all zones. Values between `0.5` and `0.9` are typical; smaller values are more stable but take more steps to reach `tf`.

## Output cadences

```
[dbl] DTd  = 5.0
[dbl] DTf  = 10.0
[dbl] DTl  = 0.5
[int] DTr  = 10000
[int] DTp  = 100
```

These control how often different types of output are written. `DTd`, `DTf`, and `DTl` are measured in code time units; `DTr` and `DTp` are measured in timesteps.

- `DTd` — interval between fluid-state dumps written to `dumps/`.
- `DTf` — legacy "full" dump interval, retained for historical compatibility and unused in current runs.
- `DTl` — interval between entries in `log.out` (for example, radial fluxes and maximum divergence error).
- `DTr` — interval between restart files, written every N timesteps.
- `DTp` — interval between performance reports printed to stdout, written every N timesteps.

For black hole accretion problems in EHT production runs, a typical value is `DTd = 5.0`.

## Problem-specific

In addition to the core parameters above, each problem registers its own parameters in `set_problem_params()` inside `prob/PROB/problem.c`. These are read from the same `param.dat` — there is no separate file — so what appears in `param.dat` depends entirely on the problem. For example, `mhdmodes2d` adds only a single parameter `nmode` (an integer selecting which linear MHD eigenmode to initialise), whereas the torus problem adds `rin`, `rmax`, `u_jitter`, `beta`, `mad_type`, explained below.

```
[dbl] rin      = 6.0
[dbl] rmax     = 12.0
[dbl] u_jitter = 4.0e-02
[dbl] beta     = 100.0
[int] mad_type = 0
```

`rin` and `rmax` define the [Fishbone-Moncrief torus](https://ui.adsabs.harvard.edu/abs/1976ApJ...207..962F/abstract): `rin` is the inner edge of the torus (in $r_g\equiv GM/c^2$) and `rmax` is the radius of the pressure maximum (also in $r_g). It is a solution to the relativistic Euler equations in a stationary, axisymmetric background for an isentropic, axisymmetric, purely azimuthal flow. The solutions are parameterized by a constant angular momentum mass density, $\ell\equiv u_{\phi}u^{t}$. The angular momentum of the torus is set analytically from `rmax`. Typically we set `rin` and `rmax` to 10 and 20 for SANE simulations, and 20 and 41 for MAD simulations.

`u_jitter` adds a small random perturbation to the initial internal energy field. This perturbs the equilibrium torus, seeds instabilities such as the magnetorotational instability (MRI), and helps trigger accretion.

`beta` is the target plasma $\beta = P_{\rm gas} / P_{\rm mag}$ used to normalise the initial magnetic field. `beta = 100` is the standard choice and it implies that the initial field is dynamically weak, which is standard for most GRMHD torus setups.


`mad_type` specifies the initial magnetic field configuration. `0` selects the standard and normal evolution (SANE), while `1` selects a magnetically arrested disk (MAD) setup. Although the steady-state magnetic flux is not trivially determined by the initial conditions, these configurations reflect setups identified in previous work as producing either SANE or MAD flows.

## Fluid parameters

```
[dbl] gam = 1.333333
```

`gam` is the adiabatic index $\gamma$ of the gas. For a relativistic gas, $\gamma = 4/3 \approx 1.333$; for a non-relativistic monatomic gas, $\gamma = 5/3 \approx 1.666667$. The latter is a good choice for single-temperature GRMHD.

## Coordinates and domain

```
[dbl] a           = 0.9375
[dbl] Rout        = 50.0
[dbl] hslope = 0.3
[dbl] mks_smooth  = 0.5
[dbl] poly_xt     = 0.82
[dbl] poly_alpha  = 14.0
```
`a` is the dimensionless black hole spin parameter, $a = J c / G M^2$, ranging from `0` (Schwarzschild) to `1` (maximally spinning Kerr).

`Rout` sets the outer radial boundary of the grid in gravitational radii ($r_g = GM/c^2$). The inner boundary is calculated based on the outer boundary and the resolution in the radial direction so that there are five zones inside the event horizon. 

`hslope` controls the degree of coordinate compression toward the midplane in MKS coordinates — smaller values concentrate more zones near $\theta = \pi/2$; `0.3` is a common choice for torus problems. `mks_smooth`, `poly_xt`, and `poly_alpha` are FMKS coordinate shape parameters (only used when `DEREFINE_POLES = 1`).

## Electron thermodynamics

```
[dbl] game    = 1.333333
[dbl] gamp    = 1.666667
[dbl] fel0    = 0.01
[dbl] tptemin = 0.001
[dbl] tptemax = 1000.0
```

These parameters are only read when `ELECTRONS = 1` is set in `parameters.h`. `game` and `gamp` are the adiabatic indices for electrons and ions respectively; the standard choices are $\gamma_e = 4/3$ (relativistic electrons) and $\gamma_p = 5/3$ (non-relativistic ions). `fel0` is the initial fraction of dissipated energy assigned to electrons; it also sets the electron entropy in the initial condition. `tptemin` and `tptemax` clamp the ion-to-electron temperature ratio $T_p/T_e$ to a physically reasonable range.