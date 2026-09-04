# PointVortex

PointVortex is a dependency-light C++20 solver for two-dimensional point vortices in an
infinite plane, square periodic box, or circular disk. It uses direct `O(N^2)` velocity kernels and stores vortex data in
contiguous `std::vector<double>` arrays. Three executables share the same model, integrators,
input, output, and checkpoint code: serial/OpenMP, distributed-memory MPI, and NVIDIA CUDA.

The periodic kernel follows the rapidly convergent square-torus construction of Weiss and
McWilliams (1991). The disk kernel uses the circle-theorem image vortex at
`R^2/conj(z_j)`, evaluated with equivalent Cartesian arithmetic.

## Quick start

From the repository root, configure and build all toolchains that are available:

```bash
cd /path/to/PointVortex
source /etc/profile        # makes nvcc visible on Arch/Omarchy
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
```

Run the supplied `params.txt` with one backend:

```bash
./build/release/point_vortex_cpu params.txt
mpirun -n 4 ./build/release/point_vortex_mpi params.txt
./build/release/point_vortex_cuda params.txt
```

The parameter file controls the physical problem; selecting an executable controls only how
velocities are calculated. Consequently, a parameter file and checkpoint can be moved between
the three backends. Each run prints its backend name before the first diagnostics line.

Existing output and checkpoint files are protected by default. For a new experiment, change
`outputFile`, `diagnosticsFile`, and `checkpointDirectory`. Use the overwrite options only when
discarding the old results is intentional.

## Choosing a backend

| Executable | Best use | Parallelism |
|---|---|---|
| `point_vortex_cpu` | Small runs, debugging, and the numerical reference | Serial with one thread; shared-memory OpenMP with multiple threads |
| `point_vortex_mpi` | Large runs across CPU sockets or compute nodes | Target vortices divided among MPI ranks; every rank retains all source vortices |
| `point_vortex_cuda` | Large direct-sum runs on an NVIDIA GPU | One CUDA thread per target vortex |

Start with CPU for correctness checks. Measure before choosing MPI or CUDA: the direct kernel is
quadratic, but communication and GPU-transfer overhead can dominate when `N` is small.

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

```text
PointVortex/
├── src/                  active solver and backend sources
├── initial_conditions/   generator source and focused guide
├── tests/                automated tests
├── build/                generated CMake and Make products (Git-ignored)
├── data/                 run output and checkpoints (Git-ignored)
├── scripts/              plotting and movie utilities
├── archive/              preserved historical implementations and backup
├── CMakeLists.txt
├── Makefile
├── params.txt
└── README.md
```

Within `src/`, `compute.cpp/.h` contains the three geometry kernels, `timestep.cpp/.h` contains
RK4 and DOPRI5, and the `backend_cpu.cpp`, `backend_mpi.cpp`, and `backend_cuda.cu` files provide
the three execution strategies. `tests/tests.cpp` exercises the model, restart format, dipole
events, and initial-condition generator.

The velocity kernel is separate from the integrator so additional geometries can be introduced
without changing the timestep algorithms.

## Requirements

- A C++20 compiler, such as GCC or Clang
- GNU Make or CMake 3.20 or newer
- Optional: an OpenMP-capable compiler/runtime
- Optional MPI version: an MPI C++ compiler and runtime
- Optional CUDA version: NVIDIA CUDA Toolkit and a CUDA-capable GPU

The CPU version has no external numerical-library dependency. CMake detects MPI and CUDA and
skips those optional executables when their toolchains are unavailable.

## Compiling with Make

Build the optimized CPU/OpenMP executables:

```bash
make
```

This creates `build/make/point_vortex_cpu`; `build/make/point_vortex` is retained as a compatible
name for the same CPU backend. To build without OpenMP:

```bash
make clean
make OPENMP=
```

Build either accelerator explicitly:

```bash
make mpi
source /etc/profile        # if nvcc is not yet on PATH in this shell
make cuda
```

Run the tests or remove generated build files with:

```bash
make test
make clean
```

## Compiling with CMake

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j
ctest --test-dir build/release --output-on-failure
```

Disable OpenMP explicitly with:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DPOINT_VORTEX_OPENMP=OFF
```

MPI and CUDA builds can be disabled independently with
`-DPOINT_VORTEX_MPI=OFF` and `-DPOINT_VORTEX_CUDA=OFF`. The default CUDA architecture is `120`
(`sm_120`, suitable for this system's RTX 5070). Override it for another GPU, for example:

```bash
cmake -S . -B build/release -DPOINT_VORTEX_CUDA_ARCHITECTURES=89
make cuda CUDA_ARCH=sm_89
```

The CMake executables are under `build/release/`; Make executables are under `build/make/`.
Both directories contain generated products from the same source code and may be deleted and
rebuilt at any time. They are alternative build workflows, not separate solver implementations.

## Running

The default parameter file is `params.txt`:

```bash
./build/release/point_vortex_cpu params.txt
```

Pass a different parameter file as the first argument:

```bash
./build/release/point_vortex_cpu my_run.params
```

All backends accept the same parameter file:

```bash
./build/release/point_vortex_cpu my_run.params
mpirun -n 4 ./build/release/point_vortex_mpi my_run.params
./build/release/point_vortex_cuda my_run.params
```

Only MPI rank zero writes trajectories, diagnostics, checkpoints, and console diagnostics.
Every rank retains the full vortex state while disjoint target ranges are evaluated in parallel.
The CUDA backend evaluates velocities on the GPU and currently keeps integration and I/O on the
host; this is a reliable baseline for validation before moving complete Runge--Kutta stages to
device memory.

### CPU and OpenMP

Use one thread for a strictly serial run:

```bash
OMP_NUM_THREADS=1 ./build/release/point_vortex_cpu run.params
```

For OpenMP, either set `numThreads` in the parameter file or use the runtime environment. The
parameter takes precedence when it is greater than zero:

```bash
OMP_NUM_THREADS=8 OMP_PROC_BIND=close ./build/release/point_vortex_cpu run.params
```

OpenMP is activated only for at least 256 target vortices, avoiding thread overhead for smaller
systems.

### MPI

Launch MPI through `mpirun`; do not start `point_vortex_mpi` independently once per terminal:

```bash
mpirun -n 4 ./build/release/point_vortex_mpi run.params
```

The current implementation is also OpenMP-capable inside each MPI rank. For pure MPI, use
`numThreads 1` or `OMP_NUM_THREADS=1`. For a deliberate hybrid run, choose ranks times threads
to match the available CPU cores, for example two ranks with eight threads each:

```bash
OMP_NUM_THREADS=8 mpirun -n 2 ./build/release/point_vortex_mpi run.params
```

All ranks must see the parameter, initial-condition, and restart files at the same paths. Output
paths need be writable only from rank zero. A fatal error aborts the MPI job so that other ranks
do not remain stuck in a collective operation.

### CUDA

Check the driver and toolkit before running:

```bash
nvidia-smi
nvcc --version
./build/release/point_vortex_cuda run.params
```

CUDA errors include the failed runtime operation. This backend retains reusable device buffers
between velocity evaluations, while positions and resulting velocities are transferred for each
Runge--Kutta stage. `numThreads` does not control CUDA thread blocks.

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
outputFile data/vortices.csv
diagnosticsFile data/diagnostics.csv
checkpointDirectory data/checkpoints
overwriteOutput false
overwriteCheckpoints false
```

### Physical and runtime parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `N` | `100` | Number of vortices used by the built-in initial condition. Ignored when loading an initial-condition file or checkpoint. |
| `timeStep` | `0.001` | Fixed RK4 step, or initial proposed DOPRI5 step. A step is shortened when needed to land on an output time or `endTime`; only DOPRI5 clamps it to the adaptive timestep bounds. |
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

The tolerance and adaptive timestep-bound parameters are parsed for both integrators, but they
affect only `dopri5`. Fixed-step RK4 uses `timeStep`, shortened only to land exactly on output
times and `endTime`.

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

The generator source and its focused guide live in `initial_conditions/`. Build it from the
repository root with either CMake:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release --target point_vortex_initial -j
```

or Make:

```bash
make initial
```

The resulting command is `./build/release/point_vortex_initial` for CMake or
`./build/make/point_vortex_initial` for Make.

### Generating a file

`point_vortex_initial` creates files accepted by every simulation backend. Geometry is required
so positions, periodic neutrality, and separation distances can be validated correctly:

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

The generator prints the actual minimum pair separation. In periodic geometry this is the
minimum-image distance across the boundaries. Random points are uniform in the configured
square for `infinite` and `periodic`, and uniform by area in a disk. Placement uses rejection
sampling when `--min-separation` is positive; an impossible density terminates with a useful
error rather than silently weakening the constraint.

Available cases are:

| Case | Vortices | Purpose |
|---|---:|---|
| `single` | 1 | Stationary infinite-plane vortex or disk image-vortex orbit |
| `pair` | 2, same sign | Co-rotating analytic pair; not neutral and therefore rejected for periodic geometry |
| `dipole` | 2, opposite signs | Translating infinite-plane pair and a neutral periodic smoke test |
| `ring` | `--count` | Equally spaced alternating-sign ring |
| `random` | `--count` | Seeded alternating-sign random population |

Useful generator options are:

```text
--geometry infinite|periodic|disk
--case random|single|pair|dipole|ring
--count N
--seed N
--min-separation D
--circulation G
--half-width L
--box-length L
--disk-radius R
--ring-radius R
--output FILE
--overwrite
```

Run `./build/release/point_vortex_initial --help` for the complete command summary. Existing
files are not replaced unless `--overwrite` is supplied.

Generated files include geometry and domain-size metadata in comments. At startup, the solver
checks this metadata against `boundaryCondition`, `boxLengthX`, `boxLengthY`, or `diskRadius`.
This prevents accidentally feeding a disk condition to a periodic run, for example. Legacy
hand-written files without metadata remain valid.

### Running a generated condition

Set the matching geometry and file in a simulation parameter file:

```text
boundaryCondition periodic
boxLengthX 2.0
boxLengthY 2.0
initialConditionFile data/initial_n400.dat
```

Then select any backend without regenerating the state:

```bash
./build/release/point_vortex_cpu run.params
mpirun -n 4 ./build/release/point_vortex_mpi run.params
./build/release/point_vortex_cuda run.params
```

The number of data rows determines the vortex count, so `N` in the simulation parameter file is
ignored when `initialConditionFile` is set. A `restartFile` takes precedence over both.

### Hand-written files

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

Without `initialConditionFile`, the legacy built-in initializer remains available. Infinite and
disk geometries create an alternating-sign ring; periodic geometry creates a seeded random
population. The dedicated generator is preferable for reproducible production runs because it
records the chosen geometry and can enforce separation.

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

The geometry-aware Python tool in `scripts/movie/make_vortex_movie.py` reads the trajectory
CSV directly. It supports fixed automatic limits on the infinite plane, wrapped display in a
square periodic box, and a circular disk boundary. Frames may contain different vortex counts
after permanent dipole removals:

```bash
python3 scripts/movie/make_vortex_movie.py data/vortices.csv --geometry infinite -o data/infinite.mp4
python3 scripts/movie/make_vortex_movie.py data/vortices.csv --geometry periodic --box-length 2 -o data/periodic.mp4
python3 scripts/movie/make_vortex_movie.py data/vortices.csv --geometry disk --radius 1 -o data/disk.mp4
```

Use `--start`, `--stop`, and `--stride` to select output-frame indices. The tool requires
Python 3.10 or newer, Matplotlib, FFmpeg, and LaTeX; its adjacent README documents all options.

### Plotting diagnostics and vortex configurations

The Jupyter notebook `scripts/analysis/point_vortex_analysis.ipynb` reads the trajectory and
diagnostics CSV files from any backend and creates vector PDF figures. Its first code cell selects
the geometry, periodic box length or disk radius, input files, number of physical-space snapshots,
and output directory. By default it reads `data/vortices.csv` and `data/diagnostics.csv` and writes:

```text
data/figures/diagnostics_evolution.pdf
data/figures/vortex_configurations.pdf
```

Launch it from the repository root:

```bash
jupyter lab scripts/analysis/point_vortex_analysis.ipynb
```

The diagnostics figure shows the Hamiltonian, geometry-relevant impulses, absolute conservation
drift, and dipole-event counts. The configuration figure selects evenly spaced saved times. It
wraps periodic coordinates into the fundamental square for display, draws the circular wall for
disk runs, and uses limits spanning the trajectory for infinite-plane runs. Positive and negative
circulations are colored red and blue, and marker area reflects circulation magnitude.

The notebook requires Python 3, NumPy, Matplotlib, and Jupyter; pandas and LaTeX are not required.
See `scripts/analysis/README.md` for a concise usage guide.

## Checkpoints and restarting

The solver writes a checkpoint for the initial state and at every subsequent output event:

```text
data/checkpoints/checkpoint_00000000.dat
data/checkpoints/checkpoint_00000001.dat
...
```

Each versioned text checkpoint stores:

- vortex positions and circulations;
- simulation time and suggested next timestep;
- accepted-step count, output index, and next output time;
- core radius, integrator type, and geometry metadata;
- the original invariant values used to measure conservation drift;
- dipole-removal configuration, cumulative event counts, random-generator state, and the
  post-event invariant reference.

Checkpoint writes use a temporary file followed by a rename. Existing indexed checkpoints are
rejected unless `overwriteCheckpoints true` is set. This option is useful when deliberately
rerunning the same configuration, but a new checkpoint directory is safer when branching a run.
Velocities are not stored because they are derived and are recomputed after loading. The
checkpoint interval is currently the same as `outputTime`.

To branch from a checkpoint, use new output filenames and a new checkpoint directory:

```text
restartFile data/checkpoints/checkpoint_00000005.dat
endTime 2.0
outputFile data/branch_vortices.csv
diagnosticsFile data/branch_diagnostics.csv
checkpointDirectory data/branch_checkpoints
```

Then run normally:

```bash
./build/release/point_vortex_cpu branch.params
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
OMP_NUM_THREADS=4 ./build/make/point_vortex_benchmark 5000 10
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
- DOPRI5 FSAL reuse and explicit cache invalidation after discrete events;
- reproducible geometry-aware initial-condition generation, separation enforcement, metadata,
  neutrality, and exact write/load round trips.

The tests should be run after changes to the kernel, integrators, or checkpoint format.
For backend validation, run the same short fixed-step case with all three executables and compare
their trajectory CSV files. Parallel floating-point operation ordering means CUDA results can
differ from CPU in the last few bits; compare with a numerical tolerance rather than requiring
byte-identical files. MPI target decomposition preserves each target's source summation order in
the present implementation and normally matches the CPU result exactly.

## Current limitations

- All velocity kernels are direct `O(N^2)` methods; tree and FMM acceleration are not present.
- MPI uses replicated source arrays, so its per-rank memory requirement remains `O(N)`.
- CUDA offloads velocity kernels but not the complete time integrator or invariant calculation.
- Weiss--McWilliams periodic dynamics currently supports square boxes and zero net circulation.
- Core regularization is implemented only for the infinite-plane kernel.
- Checkpoint frequency is tied to trajectory output frequency.
- CSV output is intentionally simple and portable, but can become large for long simulations.

Source formatting is defined by `.clang-format`; run
`clang-format -i src/*.h src/*.cpp src/*.cu initial_conditions/*.h initial_conditions/*.cpp`
after editing.
