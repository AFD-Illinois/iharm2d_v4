# Sod shock tube

## Overview

The Sod shock tube is the standard test of a code's ability to capture discontinuous wave structures. The initial state consists of two constant states separated by a discontinuity: high-density, high-pressure left state and a low-density, low-pressure right state. The resulting Riemann problem immediately produces three families of waves: a left-going rarefaction fan, a right-going contact discontinuity, and a right-going shock. Because the problem is 1D and involves no magnetic field, it tests only the hydrodynamic Riemann solver and reconstruction scheme.

## Setup

The domain is $[0,1]$ in $x$ with outflow boundaries on both sides and $N_2=1$ (a single cell in $y$, making the problem truly 1D). The initial condition is a step function at $x=0.5$,

$$
(\rho,\, p,\, v^x) =
\begin{cases}
(1.0,\; 1.0,\; 0) & x < 0.5 \\
(0.125,\; 0.1,\; 0) & x \geq 0.5
\end{cases}
$$

with $\Gamma = 1.4$ and $\mathbf{B} = 0$. The `tscale` parameter rescales the velocities and internal energies so the flow is non-relativistic ($\texttt{tscale} = 0.01 \ll 1$), recovering the classical Sod solution to high accuracy.

## Parameters

Relevant compile-time parameters are:

| Parameter | Default | Notes |
|---|---|---|
| `N1TOT`             | `256`    | Resolution; increase for sharper profiles |
| `N2TOT`             | `1`      | Keep at 1 (problem is 1D) |
| `METRIC`            | `MINKOWSKI` | |
| `RECONSTRUCTION`    | `WENO`   | |
| `X{1,2}{L,R}_BOUND` | `OUTFLOW` | |

## Output and resolution dependence

The video below shows $\rho$, $p$, $v^x$, and specific internal energy $u/\rho$ at resolutions $N = 128$, $256$, $512$, and $4096$ as a function of time.

<video controls width="100%">
  <source src="../../assets/sod_resolution.mp4" type="video/mp4">
</video>

A plotting script for individual dumps at a single resolution is provided at `prob/sod/plot_sod.py`; the multi-resolution comparison script is `prob/sod/plot_sod_resolution.py`.

At lower resolutions ($N \lesssim 256$), small-amplitude oscillations are visible in $u/\rho$ near the contact discontinuity. These arise because at low resolution the numerical dissipation smears the density jump over only a few cells (see $\rho$), and any errors in the reconstructed density profile get amplified when computing the ratio $u/\rho$ as $u$ is continuous at the contact discontinuity. The oscillations diminish with increasing resolution and are not prominent in the $N=4096$ reference run, which is included as a high-resolution benchmark.

## References

- [Sod (1978)](https://doi.org/10.1016/0021-9991(78)90023-2) — original shock tube problem.