# iharm2d_v4
`iharm2d_v4` is the 2D version of [iharm3d](https://github.com/AFD-Illinois/iharm3d) with minimal dependencies. `iharm2d_v4` parallelizes computation across all cores of a single node using `OpenMP`. This means you must have C compiler with OpenMP. Data is written out to ASCII files and a sample analysis script is provided at `scripts/analysis/simple_plot.py` to get a feel for the output format (One can always browse through `core/io.c` to get a better sense for the dump output).
TODO: Documentation page for dump format.

## Compiling and running
One can `make` the program from the output directory,
```bash
$ make -f IHARM3D_DIRECTORY/makefile PROB=PROBLEM
```
where IHARM3D_DIRECTORY is the path to your local `iharm3d` respository that contains the makefile. `PROBLEM` is the name of the problem you wish to run. NOTE: This must match the problem name in the `prob` directory. The makefile by default assumes the GCC as your compiler and that it is in your `PATH`.

The `makefile` will produce a directory named `build_archive` in your output directory. This contains all the source files (source code and header files) and a copy of the executable, `harm`. The executable is copied over to the output directory. Much like `iharm3d`, if `build_archive` already exists, `make` will prefer an newer/modified files in that directory over their equivalents in the original source. This allows you to modify the source code without disrupting the original repository.

When it comes to modifying parameters, `iharm2d_v4` follows the ideology of `iharm3d`. Compile-time parameters are present in `parameters.h` located in `build_archive` while run-time parameters are provided via an additional parameter file `param.dat`. If you're modifying compile-time parameters (eg: grid size, floors, reconstruction scheme, boundary conditions), you'll have to rerun the `make` command. A sample `param.dat` file is provided at `prob/PROBLEM/param.dat` which must be copied over to your output directory. Run-time parameters include duration of the run, domain size, fluid properties, and any other problem specific parameters.

The problem can be executed by calling the harm executable,
```bash
$ ./harm -p param.dat >OUTPUT_LOG_FILE
```
where the `STDOUT` will be redirected to `OUTPUT_LOG_FILE`. If `STDOUT` is not redirected, the run-time log will be printed on the terminal.

## Simple plotting script
TODO: Upload the file and provide instructions on executing it. Mention file output.
