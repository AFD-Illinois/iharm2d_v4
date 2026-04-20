# Output data format

`iharm2d_v4` writes all output as plain-text ASCII files. There are two file types: a **grid file** written once at the start of the run, and **dump files** written periodically throughout the simulation. Both are placed inside the `dumps/` subdirectory of your output directory, which the code creates automatically on the first dump.

## Grid file

The grid file is written to `dumps/grid` at the very start of the run (before the first timestep) and never updated again. It contains one line per physical zone, in `ZLOOP_OUT` order (X2 index varies fastest, X1 slowest), with the following columns:

| Column(s) | Description |
|-----------|-------------|
| `x`, `z` | Cartesian coordinates ($r\sin\theta$, $r\cos\theta$) for MKS; code coordinates $x$, $y$ for flat space |
| `r`, `th` | Spherical Kerr-Schild radius and polar angle (zero if `METRIC` was set to `MINKOWSKI`) |
| `X1`, `X2` | Native code coordinates |
| `gdet` | $\sqrt{-g}$ at zone centre |
| `lapse` | Lapse function $\alpha$ at zone centre |
| `gcon` | All 16 components of the contravariant metric $g^{\mu\nu}$ at zone centre |
| `gcov` | All 16 components of the covariant metric $g_{\mu\nu}$ at zone centre |

Since the grid never changes, post-processing scripts read this file once and reuse it for every dump.

## Dump files

Fluid state dumps are written to `dumps/dump_NNNNNNNN` (zero-padded eight-digit counter) at the cadence set by `DTd` in `param.dat`. 
<!-- If the code aborts unexpectedly, an emergency snapshot is written to `dumps/dump_abort` instead. -->

Each dump file has two parts: a single header line followed by one data line per physical zone.

### Header

The header is a single line of space-separated fields in the following order:

1. **Problem-specific fields** — written first by `save_problem_data()` in `prob/PROB/problem.c`. These vary by problem. For the torus, this is `mad_type`, `"torus"`, `rin`, `rmax`, `beta`, `u_jitter`; for `mhdmodes2d` it is just `nmode`. When reading a reader, you need to know which problem generated the file in order to parse this section correctly.

2. **Core fields** (always present, in this order):

| Field | Type | Description |
|-------|------|-------------|
| `version` | string | Code version string |
| `has_electrons` | int | `1` if electron thermodynamics is active |
| `gridfile` | string | Name of the companion grid file (`"grid"`) |
| `metric` | string | Coordinate system: `MINKOWSKI`, `MKS`, or `FMKS` |
| `reconstruction` | string | Reconstruction scheme: `LINEAR` or `WENO` |
| `N1`, `N2` | int | Grid dimensions |
| `n_prims` | int | Number of primitive variables (NVAR) |
| `n_prims_passive` | int | Number of passive scalars (currently always `0`) |
| `game`, `gamp`, `fel0`, `tptemin`, `tptemax` | double | Electron parameters (only present if `has_electrons = 1`) |
| `gam` | double | Adiabatic index $\gamma$ |
| `cour` | double | CFL number |
| `tf` | double | Simulation end time |
| `startx[1]`, `startx[2]` | double | Starting coordinate values |
| `dx[1]`, `dx[2]` | double | Zone spacing in each direction |
| `n_dim` | int | Number of dimensions (`4`) |
| `poly_xt`, `poly_alpha`, `mks_smooth` | double | FMKS shape parameters (`MKS` metric + `DEREFINE_POLES` only) |
| `Rin`, `Rout`, `Rhor`, `Risco` | double | Radial domain boundaries and derived radii (`MKS` metric only) |
| `hslope` | double | Midplane concentration parameter (`MKS` metric only) |
| `a` | double | Black hole spin (MKS only) |
| `t` | double | Current simulation time |
| `dt` | double | Current timestep |
| `nstep` | int | Total number of timesteps taken till dump output |
| `dump_cnt` | int | Dump file index |
| `DTd`, `DTf` | double | Dump cadences |

### Zone data

After the header, there is one line per physical zone in the same `ZLOOP_OUT` order as the grid file (X2 fastest, X1 slowest). Each line contains:

| Column(s) | Description |
|-----------|-------------|
| `P[0..NVAR-1]` | All primitive variables (see table below) |
| `jcon[0..3]` | Four-current $j^\mu$ |
| `gamma` | Lorentz factor $\Gamma$ |
| `divB` | Divergence-B estimate (zero for ghost-zone-adjacent cells) |
| `fail` | `U_to_P` failure flag (`0` = success) |
| `fflag` | Floor/fixup flag |

The primitive variables in order are:

| Index | Name | Description |
|-------|------|-------------|
| 0 | `RHO` | Rest-mass density $\rho$ |
| 1 | `UU` | Internal energy density $u$ |
| 2 | `U1` | Fluid velocity along X1 $\tilde{u}^1$ |
| 3 | `U2` | Fluid velocity along X2 $\tilde{u}^2$ |
| 4 | `U3` | Fluid velocity along X3 $\tilde{u}^3$ |
| 5 | `B1` | Magnetic field component $B^1$ |
| 6 | `B2` | Magnetic field component $B^2$ |
| 7 | `B3` | Magnetic field component $B^3$ |
| 8 | `KTOT` | Total entropy (electrons only) |
| 9+ | `KEL0`, ... | Per-model electron entropies (electrons only) |

The velocity primitives are $\tilde{u}^i\equiv u^i + u^t g^{t1}\alpha^2$. These variables vary from $-\infty$ to $\infty$ and thus are more stable than fluid 3-velocity $v^i=u^i/u^t$. The magnetic field primitives are $B^i\equiv\hspace{0.1em}^*F^{it}$, where $^*F^{\mu\nu}$ is the dual of the Faraday tensor.

## Reading dumps in Python

A minimal reader that loads density from a torus dump might look like:

```python
import numpy as np

# Read grid
grid = np.loadtxt("dumps/grid")
n1, n2 = 256, 256
r  = grid[:, 2].reshape(n1, n2)
th = grid[:, 3].reshape(n1, n2)

# Read dump — skip 1-line header
# Torus header has 6 problem-specific fields before the core fields;
# adjust the skip count for other problems.
data = np.loadtxt("dumps/dump_00000000", skiprows=1)
n_prims = 8  # NVAR without electrons
rho = data[:, 0].reshape(n1, n2)  # RHO is the first primitive
```

The `analysis/` directory in the repository contains more complete Python utilities for reading and plotting dump files.
