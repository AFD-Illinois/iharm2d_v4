# iharm2d v4

`iharm2d_v4` is a 2D GRMHD black hole accretion code that implements the HARM algorithm in C. It has minimal dependencies, requiring only a C compiler with OpenMP support and the standard math library (`libm`). It is based on the 3D production code `iharm3d`, with the intention of providing a lightweight, more approachable codebase for development, testing, and teaching. Parallelism is handled entirely through OpenMP, targeting shared-memory execution on a single node. Data is written to ASCII dump files rather than HDF5.

This documentation covers the structure of the code, the numerical algorithm and its implementation, and the format of the output data.

## References

- The original HARM paper: [Gammie et al. (2003)](https://ui.adsabs.harvard.edu/abs/2003ApJ...589..444G/abstract)
- `iharm3d` code paper: [Prather et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021JOSS....6.3336P/abstract)