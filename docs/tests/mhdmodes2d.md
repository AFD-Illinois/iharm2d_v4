# Linear MHD modes (2D)

## Overview

The same four MHD eigenmodes as in the 1D test are initialized here as oblique waves propagating at 45° in the $x_1$–$x_2$ plane. Because the wave is no longer aligned with the background magnetic field, the slow, Alfvén, and fast modes all have distinct phase speeds and involve non-trivial couplings between all eight primitives. This makes the 2D mode test a more stringent check of the multi-dimensional MHD transport: errors in the treatment of off-axis wave propagation, the divergence-free constraint on $\mathbf{B}$, and the flux-CT update of the magnetic field all show up as deviations from the analytic eigenmode solution.

## Setup

The domain is the unit square $[0,1]\times[0,1]$ in Minkowski coordinates with periodic boundaries. The background state is

$$
\rho_0 = 1,\quad u_0 = 1,\quad B^1_0 = 1,\quad B^2_0 = B^3_0 = 0,\quad \tilde{u}^i_0 = 0.
$$

The initial perturbation is a plane wave with wave-vector $(k_1, k_2) = (2\pi, 2\pi)$, giving one wavelength along each axis and a diagonal wavelength $\lambda = 1/\sqrt{2}$,

$$
q(x,y,0) = q_0 + A\,\delta q\,\cos(k_1 x + k_2 y),\qquad A = 10^{-4}.
$$

The perturbation eigenvector $\delta q$ for each mode is:

| `nmode` | Mode | Primary perturbations | $|\omega|$ |
|---|---|---|---|
| 0 | Entropy | $\delta\rho$ | $2\pi/5$ |
| 1 | Slow    | $\delta\rho,\,\delta u,\,\delta\tilde{u}^{1,2},\,\delta B^{1,2}$ | $2.410$ |
| 2 | Alfvén  | $\delta\tilde{u}^3,\,\delta B^3$ | $3.441$ |
| 3 | Fast    | $\delta\rho,\,\delta u,\,\delta\tilde{u}^{1,2},\,\delta B^{1,2}$ | $5.537$ |

The final time is set automatically to one wave period $t_f = 2\pi/|\omega|$ for the selected mode.

## Parameters

Problem-specific runtime parameters are:

| Parameter | Meaning |
|---|---|
| `nmode` | Eigenmode: `0`=entropy, `1`=slow, `2`=Alfvén, `3`=fast |

Relevant compile-time parameters are:

| Parameter | Default | Notes |
|---|---|---|
| `N1TOT`, `N2TOT`    | `64`        | Grid resolution; change for convergence study |
| `METRIC`            | `MINKOWSKI` | |
| `RECONSTRUCTION`    | `LINEAR`    | |
| `X{1,2}{L,R}_BOUND` | `PERIODIC`  | |

## Convergence

The L1 error is computed against the analytic eigenmode solution for all eight primitives at the final dump. Run at a sequence of resolutions placing each in a directory named after its resolution, then

```
python /path/to/iharm2d_v4/prob/mhdmodes2d/convergence_mhdmodes2d.py \
    -r 64,128,256,512 -m fast
```

The `-m` flag selects the mode (`entropy`, `slow`, `alfven`, `fast`). With `LINEAR` reconstruction the expected slope is $L_1 \propto N^{-2}$.

<!-- TODO: embed convergence plot -->

## References

- [Gammie, McKinney & Tóth (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..444G/abstract).
