# Linear MHD modes (1D)

## Overview

Four linear MHD eigenmodes — entropy, slow, Alfvén, and fast — are initialized as small-amplitude sinusoidal perturbations propagating along one coordinate axis through a magnetized background. Selecting `nmode` chooses which eigenmode is excited; selecting `idim` sets the propagation direction. Because the wave propagates strictly along a coordinate axis aligned with the background magnetic field, this is a genuinely 1D problem: the grid collapses to a single strip of cells ($N_2 = 1$ for `idim=1`, $N_1 = 1$ for `idim=2`). All eight primitives are tracked, and the analytic solution is known at all times, making this a precise test of the MHD wave speeds and the accuracy with which the code propagates each mode individually.

## Setup

The domain is the unit square $[0,1]\times[0,1]$ in Minkowski coordinates with periodic boundaries. The background state is

$$
\rho_0 = 1,\quad u_0 = 1,\quad \mathbf{B}_0 = \hat{e}_{\rm idim},\quad \tilde{u}^i_0 = 0,
$$

i.e. the field is aligned with the propagation direction. The initial state is

$$
q(x,0) = q_0 + A\,\delta q\,\cos(k\,x_{\rm idim}),
$$

with amplitude $A = 10^{-4}$ and $k = 2\pi$ (one wavelength in the box). The perturbation eigenvector $\delta q$ for each mode is:

| `nmode` | Mode | Perturbed variables |
|---|---|---|
| 0 | Entropy | $\delta\rho$ only |
| 1 | Slow | $\delta\rho,\,\delta u,\,\delta\tilde{u}_{\parallel}$ |
| 2 | Alfvén | $\delta\tilde{u}_{\perp}^{(1)},\,\delta B_{\perp}^{(1)}$ |
| 3 | Fast | $\delta\tilde{u}_{\perp}^{(2)},\,\delta B_{\perp}^{(2)}$ |

where $\parallel$ denotes the propagation direction and $\perp^{(1,2)}$ denote the two transverse directions. The final time is set automatically to one full wave period $t_f = 2\pi/|\omega|$ for each mode. With $\Gamma = 4/3$ the angular frequencies are: $|\omega_{\rm slow}| = 2.742$, $|\omega_{\rm Alfv\acute{e}n}| = |\omega_{\rm fast}| = 3.441$ (in units where $k=2\pi$).

## Parameters

Problem-specific runtime parameters are:

| Parameter | Meaning |
|---|---|
| `nmode` | Eigenmode: `0`=entropy, `1`=slow, `2`=Alfvén, `3`=fast |
| `idim`  | Propagation direction: `1`=x1, `2`=x2 |

Relevant compile-time parameters are:

| Parameter | Default | Notes |
|---|---|---|
| `N1TOT`             | `64`        | Set to $N$, keep `N2TOT=1` when `idim=1`; swap for `idim=2` |
| `N2TOT`             | `1`         | |
| `METRIC`            | `MINKOWSKI` | |
| `RECONSTRUCTION`    | `LINEAR`    | |
| `X{1,2}{L,R}_BOUND` | `PERIODIC`  | |

## Convergence

Because the analytic solution is known at all times, the L1 error is computed for all eight primitives ($\rho$, $u$, $\tilde{u}^1$–$\tilde{u}^3$, $B^1$–$B^3$) at the final dump. Run at a sequence of resolutions placing each in a directory named after its resolution, then

```
python /path/to/iharm2d_v4/prob/mhdmodes1d/convergence_mhdmodes1d.py \
    -r 64,128,256,512 -m fast -d 1
```

The `-m` flag selects the mode (`entropy`, `slow`, `alfven`, `fast`) and `-d` selects `idim`. With `LINEAR` reconstruction the expected slope is $L_1 \propto N^{-2}$.

<!-- TODO: embed convergence plot -->

## References

- [Gammie, McKinney & Tóth (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..444G/abstract).
