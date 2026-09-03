# PointVortex

PointVortex is a dependency-light C++20 solver for two-dimensional point vortices in an
infinite plane, square periodic box, or circular disk. It uses direct `O(N^2)` velocity kernels and stores vortex data in
contiguous `std::vector<double>` arrays. OpenMP is optional and accelerates sufficiently large
systems by distributing target vortices between threads.

The periodic kernel follows the rapidly convergent square-torus construction of Weiss and
McWilliams (1991). The disk kernel uses the circle-theorem image vortex at
`R^2/conj(z_j)`, evaluated with equivalent Cartesian arithmetic.

## Mathematical model

For vortex `i` at `(x_i,y_i)`, the infinite-plane equations are

```text
dx_i/dt = -1/(2*pi) sum_{j != i} Gamma_j (y_i-y_j)/(r_ij^2 + epsilon^2)
dy_i/dt =  1/(2*pi) sum_{j != i} Gamma_j (x_i-x_j)/(r_ij^2 + epsilon^2)
```

where `Gamma_j` is the circulation and `epsilon` is `coreRadius`. Setting
`coreRadius 0` gives the singular point-vortex equations. Coincident vortices are rejected in
that case. A positive core radius regularizes close encounters.

The solver can use fixed-step classical RK4 or adaptive Dormand–Prince 5(4). The latter controls
the largest component-wise normalized local error using absolute and relative tolerances.

## Code organization

- `vortex.h`: `VortexSystem` and `VelocityField` structure-of-arrays containers.
- `compute.cpp/.h`: infinite-plane, Weiss–McWilliams periodic, and disk image-vortex kernels,
  plus geometry-aware Hamiltonians and conserved-quantity calculations.
- `timestep.cpp/.h`: RK4 and adaptive Dormand–Prince integrators with reusable work arrays.
- `read.cpp/.h`: parameter and initial-condition parsing.
- `print.cpp/.h`: trajectory and diagnostics CSV writers.
- `checkpoint.cpp/.h`: versioned checkpoint serialization and restart loading.
- `dipole.cpp/.h`: closest-first dipole removal and optional random reinjection.
- `main.cpp`: simulation driver.
- `tests.cpp`: analytic one- and two-vortex tests plus checkpoint round-trip tests.

The velocity kernel is separate from the integrator so additional geometries can be introduced
without changing the timestep algorithms.

## Requirements

- A C++20 compiler, such as GCC or Clang
- GNU Make or CMake 3.20 or newer
- Optional: an OpenMP-capable compiler/runtime

No external numerical library is required.

## Compiling with Make

Build the optimized OpenMP executable:

```bash
make
```

This creates `point_vortex`. To build without OpenMP:

```bash
make clean
make OPENMP=
```

Run the tests or remove generated build files with:

```bash
make test
make clean
```

## Compiling with CMake

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Disable OpenMP explicitly with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPOINT_VORTEX_OPENMP=OFF
```

The CMake executable is `build/point_vortex`.
The benchmark executable is `build/point_vortex_benchmark`.

## Running

The default parameter file is `params.txt`:

```bash
./point_vortex
```

Pass a different parameter file as the first argument:

```bash
./point_vortex my_run.params
```

Paths inside the parameter file are interpreted relative to the directory from which the
program is run, not relative to the parameter file.

## Parameter-file format

Each non-empty line contains a key and a value separated by whitespace. A `#` begins a
comment. Unknown keys and invalid values cause the program to stop with an error.

```text
N 100
timeStep 0.001
endTime 1.0
outputTime 0.1
coreRadius 0.0
integrator dopri5
absoluteTolerance 1e-10
relativeTolerance 1e-8
minimumTimeStep 1e-12
maximumTimeStep 0.05
numThreads 0
boundaryCondition infinite
boxLengthX 2.0
boxLengthY 2.0
periodicImageLayers 8
randomSeed 1234567
diskRadius 1.0
dipoleRemoval false
dipoleRemovalDistance 0.01
dipoleReinjection none
outputFile vortices.csv
diagnosticsFile diagnostics.csv
checkpointDirectory checkpoints
overwriteOutput false
overwriteCheckpoints false
```

### Physical and runtime parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `N` | `100` | Number of vortices used by the built-in initial condition. Ignored when loading an initial-condition file or checkpoint. |
| `timeStep` | `0.001` | Nominal fixed RK4 step, or initial proposed DOPRI5 step. A step is shortened when needed to land on an output time or `endTime`. |
| `endTime` | `1.0` | Absolute time at which the run stops. For a restart, it must not precede the checkpoint time. |
| `outputTime` | `0.1` | Interval used for trajectory, diagnostics, and checkpoint output. Timesteps are clipped to land on it; `endTime` is also output when it is off this cadence. |
| `coreRadius` | `0.0` | Regularization radius `epsilon`. Zero selects the singular model. |
| `boundaryCondition` | `infinite` | `infinite`, `periodic`, or `disk`. |
| `boxLengthX` | `2.0` | Periodic-box side length in x. |
| `boxLengthY` | `2.0` | Periodic-box side length in y; it must equal `boxLengthX`. |
| `periodicImageLayers` | `8` | Symmetric truncation count for the Weiss–McWilliams sums; allowed range 0–64. |
| `randomSeed` | `1234567` | Seed for the built-in random periodic initial condition and dipole reinjection. Change it to generate a different realization. |
| `diskRadius` | `1.0` | Radius of the circular impermeable boundary. |
| `numThreads` | `0` | OpenMP thread count. Zero leaves selection to the OpenMP runtime. |

### Dipole removal and reinjection

| Parameter | Default | Meaning |
|---|---:|---|
| `dipoleRemoval` | `false` | Enable removal of close, opposite-sign vortex pairs after initialization and every accepted timestep. |
| `dipoleRemovalDistance` | `0.01` | Strict distance threshold below which an opposite-sign pair is removed. It must be positive, even when removal is disabled. |
| `dipoleReinjection` | `none` | `none` permanently removes each pair; `independent` or `paired` reinjects it in a periodic box or disk. |

At each check, all eligible pairs are sorted by distance. The closest available pair is selected
first and each vortex can participate at most once in that check. Infinite-plane and disk
distances are Euclidean; periodic distances use the minimum-image convention. The initial state
is checked before its first output, and subsequent checks occur after every accepted RK4 or
DOPRI5 timestep. Rejected adaptive attempts do not trigger events.

Reinjection preserves the two removed circulation values and therefore preserves vortex count.
In `independent` mode both positions are sampled independently and uniformly by area. In
`paired` mode the first position is uniform, while the second is placed in a uniformly random
direction at the mean-spacing distance

```text
sqrt(domain area / vortex population before removal).
```

Periodic positions are wrapped into the box. Disk samples use `r = diskRadius*sqrt(U)`, and a
paired disk placement is resampled until both vortices are strictly inside the boundary.
Reinjection is rejected for the infinite plane. `randomSeed` initializes the reinjection random
generator as well as the built-in periodic distribution; checkpoints preserve its complete state
so restarted reinjection sequences are reproducible.

The periodic kernel requires zero total circulation. Removing a pair without reinjection will
preserve this requirement when the two circulations are equal and opposite, as in the built-in
periodic initial condition. Removing unequal opposite-sign strengths can violate periodic
neutrality and will then be rejected by the periodic velocity kernel.

For compatibility, `OutputTime` is an alias for `outputTime`, `coreSize` is interpreted as
`coreRadius^2`, and the legacy `numSteps` key sets `endTime = timeStep * numSteps` when
`endTime` is not also present. Prefer the names shown in `params.txt` for new runs.

### Geometry details

For `boundaryCondition periodic`, the box must currently be square and the total circulation
must be zero. The velocity uses the complementary, rapidly convergent image sums derived by
[Weiss and McWilliams (1991)](https://atoc.colorado.edu/~jweiss/website/publications/WeissMcWilliams1991.pdf).
Increasing `periodicImageLayers` reduces truncation error at additional cost. Eight layers are
normally already close to double-precision convergence for a square box.

For `boundaryCondition disk`, every vortex must lie strictly inside `diskRadius`. The kernel is
the usual complex-plane circle-theorem construction: a vortex at `z_j` has an opposite-sign
image at `diskRadius^2/conj(z_j)`. The implementation evaluates this rational expression using
real coordinates; no cotangent is involved. `cot()` is associated with periodic rows of images,
not inversion in a circular boundary.

### Integrator parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `integrator` | `dopri5` | `rk4` for fixed-step RK4 or `dopri5` for adaptive Dormand–Prince 5(4). |
| `absoluteTolerance` | `1e-10` | Absolute local-error scale used by DOPRI5. |
| `relativeTolerance` | `1e-8` | Relative local-error scale used by DOPRI5. |
| `minimumTimeStep` | `1e-12` | Smallest allowed adaptive step. |
| `maximumTimeStep` | `0.1` | Largest allowed adaptive step. The initial step is clamped to these bounds. |

The tolerance and timestep-bound parameters are parsed for both integrators, but error control is
only used by `dopri5`.

### Input, output, and restart parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `initialConditionFile` | unset | Load arbitrary vortices from a text file instead of creating the built-in initial condition. |
| `restartFile` | unset | Restore a simulation checkpoint. This takes precedence over `N` and `initialConditionFile`. |
| `outputFile` | `vortices.csv` | Trajectory CSV filename. |
| `diagnosticsFile` | `diagnostics.csv` | Conserved-quantity CSV filename. |
| `checkpointDirectory` | `checkpoints` | Directory for indexed checkpoint files. |
| `overwriteOutput` | `false` | Permit replacement of trajectory and diagnostics CSV files. |
| `overwriteCheckpoints` | `false` | Permit replacement of indexed checkpoint files in `checkpointDirectory`. |

## Initial conditions

Without `initialConditionFile`, the infinite-plane and disk geometries place `N` vortices on a
circle with alternating circulations `+1,-1,+1,...`. Its radius is 1 in the infinite plane and
one half of the disk radius. The periodic geometry instead distributes positions uniformly at
random throughout the box, using `randomSeed`, while retaining alternating circulations so the
total circulation is zero. The built-in periodic condition therefore requires an even `N`.

For a general system, provide one vortex per line as `x y circulation`. Whitespace and commas
are both accepted, and comments begin with `#`:

```text
# x, y, circulation
-1.0, 0.0, 1.0
 1.0, 0.0, 1.0
 0.0, 2.0, -0.5
```

Then add:

```text
initialConditionFile initial.csv
```

The number of rows determines the vortex count.

## Output files

The trajectory file contains one row per vortex at every output event:

```text
time,index,x,y,circulation,u,v
```

Velocities are evaluated at the saved positions. The diagnostics file contains:

```text
time,circulation,linear_impulse_x,linear_impulse_y,angular_impulse,hamiltonian,delta_circulation,delta_linear_impulse_x,delta_linear_impulse_y,delta_angular_impulse,delta_hamiltonian,segment_delta_circulation,segment_delta_linear_impulse_x,segment_delta_linear_impulse_y,segment_delta_angular_impulse,segment_delta_hamiltonian,removed_pairs,reinjected_pairs
```

The `delta_*` columns are differences from the original state before any time-zero dipole event.
The `segment_delta_*` columns instead use the state immediately after the most recent removal or
reinjection event as their reference. They reset to zero at an event and then isolate numerical
conservation drift until the next event. The final two columns contain cumulative event counts.
Removal and reinjection are discrete model operations, so they can cause physical jumps in the
original-reference `delta_*` values. Values are written with 17
significant digits and flushed after every output frame, so diagnostics are visible while a run
is still in progress. Existing CSV files are rejected by default to prevent accidental data loss.

The Hamiltonian is geometry-aware. Periodic momentum is conserved only when positions remain
unwrapped, as they do here; angular impulse is not a periodic invariant. In a disk, angular
impulse is conserved but linear impulse generally is not. Regularization is currently supported
only in the infinite plane.

At every output time, the same relevant conserved quantities, original-reference drift, and
`segmentD*` drift are printed to the terminal:

- infinite plane: circulation, Hamiltonian, both linear impulses, and angular impulse;
- periodic box: circulation, Hamiltonian, and both unwrapped linear impulses;
- disk: circulation, Hamiltonian, and angular impulse.

The diagnostics CSV keeps one fixed schema for all geometries. Columns that are not invariants
of the selected geometry are still recorded for analysis, but are omitted from terminal output.

### Creating a movie

The geometry-aware Python tool in `../scripts/movie/make_vortex_movie.py` reads the trajectory
CSV directly. It supports fixed automatic limits on the infinite plane, wrapped display in a
square periodic box, and a circular disk boundary. Frames may contain different vortex counts
after permanent dipole removals:

```bash
python3 ../scripts/movie/make_vortex_movie.py vortices.csv --geometry infinite -o infinite.mp4
python3 ../scripts/movie/make_vortex_movie.py vortices.csv --geometry periodic --box-length 2 -o periodic.mp4
python3 ../scripts/movie/make_vortex_movie.py vortices.csv --geometry disk --radius 1 -o disk.mp4
```

Use `--start`, `--stop`, and `--stride` to select output-frame indices. The tool requires
Python 3.10 or newer, Matplotlib, FFmpeg, and LaTeX; its adjacent README documents all options.

## Checkpoints and restarting

The solver writes a checkpoint for the initial state and at every subsequent output event:

```text
checkpoints/checkpoint_00000000.dat
checkpoints/checkpoint_00000001.dat
...
```

Each versioned text checkpoint stores:

- vortex positions and circulations;
- simulation time and suggested next timestep;
- accepted-step count, output index, and next output time;
- core radius, integrator type, and geometry metadata;
- the original invariant values used to measure conservation drift.
- dipole-removal configuration, cumulative event counts, random-generator state, and the
  post-event invariant reference.

Checkpoint writes use a temporary file followed by a rename. Existing indexed checkpoints are
rejected unless `overwriteCheckpoints true` is set. This option is useful when deliberately
rerunning the same configuration, but a new checkpoint directory is safer when branching a run.
Velocities are not stored because they are derived and are recomputed after loading. The
checkpoint interval is currently the same as `outputTime`.

To branch from a checkpoint, use new output filenames and a new checkpoint directory:

```text
restartFile checkpoints/checkpoint_00000005.dat
endTime 2.0
outputFile branch_vortices.csv
diagnosticsFile branch_diagnostics.csv
checkpointDirectory branch_checkpoints
```

Then run normally:

```bash
./point_vortex branch.params
```

The configured `coreRadius`, `integrator`, boundary type, and active geometry parameters must
match the checkpoint. Conservation drift remains relative to the beginning of the original run.
New version-4 checkpoints contain geometry and dipole-event metadata, including the post-event
invariant reference. Older version-1 through version-3 checkpoints remain readable when their
configuration is compatible; a restart from an older checkpoint begins a new local-drift segment.

## Performance notes

- The direct calculation scales as `O(N^2)` per velocity evaluation.
- Positions, circulations, and velocities use separate contiguous vectors for cache-friendly
  access and vectorization.
- Runge–Kutta stage storage is allocated once and reused without stage-circulation copies.
- DOPRI5 uses FSAL reuse after accepted steps and retains its first stage across rejected attempts.
- OpenMP parallelism is enabled for systems of at least 256 vortices; smaller systems avoid its
  thread-launch overhead.
- DOPRI5 needs multiple velocity evaluations per attempted step, so tolerance choice can have a
  large effect on runtime.
- `numThreads 0` uses the OpenMP runtime default. You can also control placement and affinity
  with standard environment variables such as `OMP_NUM_THREADS` and `OMP_PROC_BIND`.

The direct kernel is intended to remain the correctness reference when faster tree, FMM, or
particle-mesh methods are added later.

Build and run the standalone infinite-plane kernel benchmark with:

```bash
make benchmark
OMP_NUM_THREADS=4 ./benchmark 5000 10
```

The arguments are vortex count and repetition count. The program uses a fixed random seed and
reports elapsed time and pair interactions per second, making thread-count and compiler-option
comparisons reproducible.

## Tests

`make test` checks:

- a single stationary vortex;
- the analytic velocity of an equal-circulation pair;
- the two-vortex orbit with RK4;
- the two-vortex orbit with adaptive DOPRI5;
- checkpoint write/load round trips;
- periodic translation invariance and the disk single-vortex image solution;
- reproducible, neutral random periodic initialization;
- closest-first dipole removal, periodic minimum-image matching, both reinjection modes, and
  restartable random-generator state;
- DOPRI5 FSAL evaluation counts.

The tests should be run after changes to the kernel, integrators, or checkpoint format.

## Current limitations

- All velocity kernels are direct `O(N^2)` methods; tree and FMM acceleration are not present.
- Weiss--McWilliams periodic dynamics currently supports square boxes and zero net circulation.
- Core regularization is implemented only for the infinite-plane kernel.
- Checkpoint frequency is tied to trajectory output frequency.
- CSV output is intentionally simple and portable, but can become large for long simulations.

Source formatting is defined by `.clang-format`; run `clang-format -i *.h *.cpp` after editing.
