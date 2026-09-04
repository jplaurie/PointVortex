# Initial-condition generator

This directory contains the geometry-aware input generator used by PointVortex:

- `generate_initial.cpp` provides the `point_vortex_initial` command-line program.
- `initial_condition.cpp/.h` provide generation, separation measurement, geometry validation,
  and file writing for the program and automated tests.

Build from the repository root. With CMake:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --target point_vortex_initial -j
```

Or with Make:

```bash
make initial
```

Generate a neutral periodic random state:

```bash
./build/release/point_vortex_initial \
    --geometry periodic \
    --case random \
    --count 400 \
    --seed 20261376 \
    --box-length 2 \
    --min-separation 0.01 \
    --output data/initial_n400.dat
```

For the Make build, use `./build/make/point_vortex_initial` instead. Run with `--help` to list
every option. The output file contains `x y circulation` rows and metadata comments understood
by the simulation driver.

In a simulation parameter file, select the same geometry and domain size:

```text
boundaryCondition periodic
boxLengthX 2.0
boxLengthY 2.0
initialConditionFile data/initial_n400.dat
```

The generated file can then be used unchanged by CPU, MPI, or CUDA. The driver rejects geometry
and domain-size mismatches before starting a simulation.
