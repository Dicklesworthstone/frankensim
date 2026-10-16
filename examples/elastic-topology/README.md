# A three-dimensional, independently loaded elasticity topology study

This executable connects the existing P1 tetrahedral elasticity, Helmholtz
filter, Heaviside projection, SIMP sensitivities, and independent-load solves
to a feasible density optimizer. It produces a changed design and exports the
actual accepted displacement fields and iteration history.

It is a **fixed-mesh density study**, not the free-boundary CutFEM elasticity
successor tracked by `q61wp.16`, and does not enable 3-D elasticity in the native
`frankensim study` command. It does not provide ledger packages or a resumable
checkpoint. No claim of stationarity, optimality, mesh convergence, stress
safety, experimental validation, or a manufacturable binary solid is made.

## Run a study

From the repository root, with the workspace's usual toolchain and sibling
dependencies available:

```sh
cargo run -p fs-topopt --example elastic_topology -- \
  --output topology-run \
  --cells 2 --iterations 30 --volume 0.5 \
  --length 0.2 --width 0.1 --height 0.1 \
  --youngs 7e10 --poisson 0.3 --filter-radius 0.015 \
  --load-y 1 --load-z -1 --weight-y 0.5 \
  --seconds 30 --linear-iterations 5000 --total-linear-iterations 100000
```

`topology-run` must not exist; its parent must exist. Input admission failures
and stops before a complete baseline return an error without a solved-design
export. After a baseline, cancellation, work exhaustion or a failed numerical
trial preserve the last accepted fields for export. A numerical failure still
returns a nonzero exit after exporting that prefix. I/O errors propagate.
A failed write can leave a partial directory; the final summary is written
only after the mesh and iteration files have been flushed. A partial summary
after an I/O interruption is not a completion receipt.

The rectangular cantilever is clamped on `x=0`. Each total end load is shared
equally across the nodes at `x=length`; this is a specified nodal-force model,
not exact integration of a continuous traction. The y and z loads act in
**separate** equilibrium solves. The objective is
`weight_y * compliance_y + (1 - weight_y) * compliance_z`, not the compliance
of their vector sum. Force, coordinates, displacement, modulus and compliance
use N, m, m, Pa and N*m respectively.

Run `--help` for controls. Mesh allocation is bounded to 1..8 subdivisions per
axis (6..3072 tetrahedra), and the accepted-update budget is capped at 200.
Dimensions are admitted in 1e-4..1e3 m, Young's modulus in 1e3..1e13 Pa, and
signed loads within +/-1e9 N. These input bounds are not a solver-convergence
or physical-validity guarantee. Extreme aspect ratios or conditioning can
still be refused by the existing component solvers.

## Inspect the result

`design.vtk` is an ASCII VTK unstructured tetrahedral grid. It contains the
raw and projected density on every cell and the accepted displacement vector
for **each independent load case** at every vertex. A visualization program
can color by `projected_density` and warp by either displacement field. The
file does not threshold density into a certified solid or claim that its
boundary is watertight or mesh-converged.

`iterations.csv` includes the baseline and every accepted update. Each row
reports weighted and individual compliance, physical projected-volume
fraction, and maximum raw-density change for the **same solved design**.
`summary.txt` records physical inputs, numerical controls, terminal reason,
initial/final objective, detailed evaluation stop, actual linear solves and
Krylov iteration consumption. Work in rejected or unfinished trials is counted,
but their densities, fields and gradients never replace the accepted result.

The starting uniform density is found by inverting the actual
filter/projection at the requested volume cap. Consequently, improvement is
measured against a feasible baseline under the **same loads and material
budget**, never against an inadmissible full-material reference.

## Stops, budgets and API boundary

The optimizer enforces actual projected volume and accepts only non-increasing
compliance. Its OC denominator uses the full projected-volume derivative
`F^T(H' * V / sum(V))`. Log-space multiplier bracketing avoids an arbitrary
fixed multiplier range tied to load units. The older fixed-step
`optimality_criteria` and `robust_optimality_criteria` algorithms are unchanged.

The example uses `controlled_multi_load_optimality_criteria` with one
`SolveControl` spanning baseline preparation and optimization. The convenience
`multi_load_optimality_criteria` also polls its callback inside the same filter,
elasticity and transpose solves, with the default component iteration limits.
Public `try_*` filter/pipeline APIs allow the same controlled evaluations outside
an optimizer. Returned evaluation stops restore the prior operator and carry
no partial gradient. Invalid model inputs retain the component panic contract.

`--linear-iterations` caps each CG solve (default 50000, with the filter's
existing tighter 20000 cap). `--total-linear-iterations` caps all Krylov work
in the run (default 2000000). Both may be zero; zero is not silently replaced
with a default. The total includes baseline inversion, every load, volume
search, rejected backtracking trial and transposed filter. Exactly exhausted
work does not invalidate an already completed solve or accepted state.

`--seconds` now covers baseline inversion as well as optimization. Matrix
assembly and output I/O remain outside the wall-time budget. The callback is
polled every **at most 32 Krylov iterations**, during filter scatter tiles,
between algebraic stages, and before evaluation publication. A single matrix
application is not preemptible, so this is not a millisecond cancellation
latency guarantee. No process signal handler or on-disk restart is implied.

Terminal reasons distinguish update-budget exhaustion, the design-change
threshold, wall-time stop, linear-work exhaustion, numerical failure, and
failure to find an acceptable bounded step. **Design change is not a KKT or
optimality test.** CG convergence retains its recursively estimated residual
semantics; this work does not promote it to a certified Euclidean error bound.
A stopped run exports the last accepted solved fields without another solve.

No random generator is used. Fixed update/work budgets are intended to replay
on the same admitted toolchain/profile; wall-time stopping can select different
prefixes. Focused Rust checks exercise real filter/P1 elasticity, independent
loads, derivative consistency, operator/field rollback inside nested solves,
work limits, cancellation and repeated runs:

```sh
cargo test -p fs-topopt --lib
cargo test -p fs-topopt --test controlled_evaluation
cargo test -p fs-topopt --test multi_load_oc
cargo test -p fs-topopt --example elastic_topology
```

The editing environment did not have Rust/DSR/RCH, so these Rust tests and the
executable were **not run there**. Independent Python calculations checked the
P1 sensitivities and feasible improvement in the preceding implementation;
168 additional recurrence checks compared batched/uninterrupted arithmetic
and iteration ceilings. These are numerical cross-checks, not execution of the
Rust path or a claim that its build/tests pass.
