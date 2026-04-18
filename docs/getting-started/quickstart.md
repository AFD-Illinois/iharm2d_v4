# Quickstart

This page walks through obtaining the code, compiling it, and running your first problem. By the end you should have a working `harm` executable and a set of output dumps from one of the bundled test problems.

`iharm2d_v4` has no external library dependencies beyond the C standard library, so getting up and running is straightforward, all you need is `gcc` with OpenMP support and GNU `make`.

Before diving in, here are a few terms that will come up throughout this guide:

- **Repository (repo)**: The folder containing all of the source code, downloaded from GitHub. While you can, it is recommended you don't run simulations from inside this folder, it just holds the original code.
- **Executable (harm)**: The program that actually runs your simulation. It is produced by compiling the source code with make.
- **Output directory**: A folder you create to hold everything related to a particular simulation run — the executable, parameter files, and all of the data the simulation produces. Keeping this separate from the repository means you can run many different simulations without cluttering or modifying the original source.
- **Compile-time parameters**:Settings that get baked into the executable when you build it, so changing them means recompiling. Each problem has its own `parameters.h` file in its subdirectory under `prob/`, and these are the ones you will likely find yourself editing most often. Examples of compile-time parameters is the grid resolution: `N1TOT`, `N2TOT`.
- **Run-time parameters (param.dat)**: Options that the executable reads from a file called `param.dat` each time it is executed. You can change these between runs wihtout recompiling the code. Each problem comes with a sample `param.dat` in its subdirectory `prob/` that you can use as a starting point.
- **Dumps**: Snapshot files that the code writes out periodically during a simulation, containing the full fluid state (density, velocity, magnetic field, etc.) at that moment in time. These are what you analyze and plot after the run finishes.

## Cloning the repository

```
git clone https://github.com/AFD-Illinois/iharm2d_v4.git
```

## Compiling

It is recommended to compile the code from a separate output directory (`my_run` in the example below) rather than from inside the repository itself. This keeps your build artifacts, parameter files, and simulation output cleanly separated from the source. Just create a directory wherever you want your run to live, then call `make` from there:
```
mkdir my_run && cd my_run
make -f /path/to/iharm2d_v4/makefile PROB=orszag_tang
```
`PROB` must match the name of a subdirectory under `prob/`—for example `torus`, `bondi`, `orszag_tang`. You can also write your own problem from scratch (see [TODO]). The makefile assumes `gcc` is in your `PATH`.

This produces a `build_archive/` directory inside your output directory containing copies of all the source and header files, along with the compiled executable `harm`. A copy of `harm` is also placed in the output directory for convenience. If `build_archive/` already exists from a previous build, make will prefer any newer or modified files in that directory over their equivalents in the original repository. This is a handy feature — it means you can edit source files in `build_archive/` (for instance to tweak the grid resolution) without ever altering the original repo.

## Parameters

`iharm2d_v4` splits its configuration into compile-time and run-time parameters, following the same convention as `iharm3d`.

**Compile-time parameters** live in `build_archive/parameters.h`. These include, for example, the grid resolution (`N1TOT`, `N2TOT`), the coordinate system (`METRIC`, `DEREFINE_POLES`), the boundary conditions (`X1L_BOUND`, `X1R_BOUND`, `X2L_BOUND`, `X2R_BOUND`), and whether electron thermodynamics is enabled (`ELECTRONS`). Since these are baked into the code at compile time, changing any of them means you will need to rerun `make`.

**Run-time parameters** are read from `param.dat` when the executable is run. Each problem comes with a sample `param.dat` under `prob/PROBLEM/param.dat`—just copy it into your output directory before running:
```
cp /path/to/iharm2d_v4/prob/torus/param.dat .
```

## Running

Once everything is in place, run the code by calling the `harm` binary with the `-p` flag pointing to the parameter file:
```
./harm -p param.dat
```
Runtime logs are printed to stdout. If you would like to save them to a file:
```
./harm -p param.dat > run.log
```

You can also specify a separate output directory with the `-o` flag if you prefer not to run from within the output directory itself. The code will create `dumps/` and `restarts/` subdirectories automatically. Fluid state dumps are written as ASCII files at the cadence set by the `DTd` parameter in `param.dat`, and grid data is written once at the start of the run.