# Code and documentation review — 2026-09-05

The active solver, its CPU/OpenMP, MPI and CUDA backends, initial-condition generator, build
files, checkpoint/output code, supported analysis notebook, and movie tool were reviewed against
the README. The review found and fixed the issues below. No failures remain in the checks listed
here. This is a source review and regression-test pass, not a proof for every possible input or
platform.

## Findings and fixes

| Area | Finding | Resolution |
|---|---|---|
| Periodic kernel | Subtracting `cosh(a) - cos(b)` rounded to zero for a distinct pair separated by `1e-10`. | Use the equivalent sum of squared half-angle functions in CPU/CUDA velocities and the periodic Hamiltonian. Actual coincidences still fail. |
| Disk kernel | Explicit inversion of a source near the center (`x = 1e-155`) overflowed and produced NaN velocities. | Evaluate the image contribution with normalized Cartesian coordinates without constructing the inverse point. |
| Disk dipole removal | Removal considered only real-real pairs, so a real vortex approaching the wall was never treated as a dipole with its opposite-sign image. | Add wall-image candidates using distance `(R²-r²)/r`; remove the real vortex, and reinject one real vortex when reinjection is enabled. |
| Integrators | NaN errors could be ignored by `std::max`, allowing an invalid adaptive step to be accepted. RK4 could also return a non-finite final state. | Validate numerical inputs/results and reject non-finite states. RK4 validates its proposed state before committing it. |
| DOPRI5 | The accepted fifth-order state was recomputed in a different floating-point order from the state used for its cached final derivative. Subtracting two complete solutions also lost precision in the error estimate. | Accept the stage-7 positions directly and form the embedded error from the difference of weights. |
| Output scheduling | A unit-time floor in rounding tolerance was inappropriate for very small simulation times; an end time very close to an output boundary could produce two final frames. | Scale rounding tolerance with simulation time and record a single final frame at the configured end time. |
| Input and checkpoints | Malformed initial-condition rows could be silently skipped. Corrupt checkpoint counts could request a huge allocation before rows were checked. | Reject malformed nonempty rows and invalid counters; validate checkpoint rows as they are read, schedules, and invariant values. |
| File protection | Aliased CSV/checkpoint/input paths could replace source files when overwrite flags were enabled. Checkpoint replacement removed the old destination before rename. | Reject conflicting destinations and rename the completed checkpoint directly without first unlinking its destination. |
| Make and MPI | Some Make targets missed header changes; the CPU backend name claimed OpenMP even in serial builds. Hybrid MPI used initialization without requesting thread support. | Track headers for all affected targets, report the compiled CPU mode, request/check MPI thread support, and guard MPI count conversions. |
| Generator and benchmark | Generator validation happened after placement, output parents were not created, and final stream errors could be missed. Invalid benchmark arguments could terminate through an uncaught exception. | Validate generator options before placement, create parents, check file close, and reject invalid benchmark arguments cleanly. |
| Analysis and movies | The notebook failed on empty trajectories or populations containing only zero-strength vortices. The movie reader did not reject decreasing times, and empty selections could overwrite an existing movie. | Handle those notebook cases, validate movie time order, and check input/output paths and frame availability before opening movie output. |

The README and focused guides now describe these behaviors, the generator's stationary
single vortex at the origin, impulse-column definitions, independent output schedules,
random-sampling reproducibility limits, and the conditional distribution used by paired disk
reinjection. Stored notebook outputs were cleared so they do not describe an earlier run.

## Validation performed

- Release CMake build and all three CTest suites: existing numerical tests, new numerical
  audit tests, and output/restart integration tests.
- Debug CPU build with AddressSanitizer, UndefinedBehaviorSanitizer, leak detection, and C++
  library assertions: all three suites passed. Leak checks ran outside the execution sandbox
  because its tracing restrictions prevent LeakSanitizer from operating.
- Analytic close-pair and near-center checks; finite-difference Hamiltonian derivatives and
  short Hamiltonian-conservation runs in all three geometries.
- CPU-versus-CUDA and CPU-versus-two-rank-MPI comparisons: eight cases each, covering all
  geometries, short trajectories, 258-vortex OpenMP/hybrid populations, and numerical edge cases.
  Comparisons use relative tolerance `2e-10` and absolute tolerance `2e-11`.
- Expanded output/restart tests on CPU, CUDA, and two MPI ranks, including unchanged/changed
  schedules, version-4 checkpoints, malformed input, corrupt counters, empty populations,
  path collisions, and final-time rounding. Existing readers for versions 1–3 were reviewed;
  this pass did not add separate version-1/2/3 integration fixtures.
- Make builds for CPU, MPI, CUDA, and the generator; Make numerical and output tests. A dry
  run after a simulated parameter-header change confirmed affected targets are rebuilt.
- Real generator-to-solver-to-plot/movie smoke tests for infinite, periodic, and disk geometry:
  nine notebook cases (ordinary, zero-strength, and empty trajectories) and three MP4 renders.
- Python/notebook syntax checks and `git diff --check`.

The reusable tests are `tests/audit_tests.cpp`, `tests/output_integration.py`,
`tests/backend_consistency.py`, and `tests/tooling_tests.py`. The README lists their commands.

## Scope and remaining limitations

Historical implementations in `archive/` and the old numbered-snapshot plotting scripts are
not part of the active build or certified by these tests. Their status and incompatible input
format are now explicit in the README.

Singular encounters or disk-boundary violations during intermediate integration stages stop
the run; adaptive integration does not automatically retry those exceptions at a smaller step.
The direct kernels remain quadratic in population, and extreme scales can exceed double
precision. Changing output intervals changes integration boundaries and can change dipole-event
timing. These limitations are documented rather than hidden by the tests.

Only the available local toolchains and hardware were exercised. This pass does not establish
cross-platform bitwise reproducibility or validate long scientific production runs. The existing
absence of a project license remains a separate publication decision for the authors.
