# A three-dimensional, independently loaded elasticity topology study

This executable connects the existing P1 tetrahedral elasticity, Helmholtz
filter, Heaviside projection, SIMP sensitivities, and independent-load solves
to a feasible density optimizer. It produces a changed design and exports the
actual accepted displacement fields and iteration history.

It is a **fixed-mesh density study**, not the free-boundary CutFEM elasticity
successor tracked by `q61wp.16`, and does not enable elasticity in the native
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
  --load-y 1 --load-z -1 --weight-y 0.5
```

`topology-run` must not exist; its parent must exist. Input admission and
numerical failures return an error rather than producing an apparent success.
I/O errors propagate. A failed write can leave a partial directory; the final
summary is written only after the mesh and iteration files have been flushed.
A partial summary after an I/O interruption is not a completion receipt.

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
`summary.txt` records the physical inputs, numerical controls, terminal reason,
and initial/final objective. No measured benchmark values are frozen into it.

The starting uniform density is found by inverting the actual
filter/projection at the requested volume cap. Consequently, improvement is
measured against a feasible baseline under the **same loads and material
budget**, never against an inadmissible full-material reference.

## Stops, evidence and API boundary

The optimizer enforces actual projected volume and accepts only non-increasing
compliance. Its OC denominator uses the full projected-volume derivative
`F^T(H' * V / sum(V))`. Log-space multiplier bracketing avoids an arbitrary
fixed multiplier range tied to load units. The older fixed-step
`optimality_criteria` and `robust_optimality_criteria` APIs are unchanged; this
example uses the new `multi_load_optimality_criteria` API.

Terminal reasons distinguish update-budget exhaustion, the design-change
threshold, a wall-time stop, and failure to find an acceptable bounded step.
**Design change is not a KKT or optimality test.** `--seconds 10` adds a
wall-time budget over optimization, excluding setup and output I/O. Checkpoints
are between component evaluations and multiplier trials; a running linear
solve cannot yet be interrupted through this API. A time stop exports the
last accepted solved fields without starting another physics solve. A stop
before any equilibrium produces an error and no solved-design export.

No random generator is used. Fixed update budgets are intended to replay on
the same admitted toolchain/profile; wall-time stopping can select different
prefixes. The added Rust regressions cover independent opposite loads, aligned
final fields, volume/move limits, replay and cancellation. They require actual
execution in the project toolchain:

```sh
cargo test -p fs-topopt --test multi_load_oc
cargo test -p fs-topopt --lib multi_load::tests
cargo test -p fs-topopt --example elastic_topology
```

The editing environment did not have Rust/DSR/RCH, so these Rust tests and the
executable were **not run there**. An independent dense-P1 Python calculation
checked the volume/compliance gradients, feasible improvement, and load-scale
invariance; that is a numerical cross-check, not execution of this Rust path.
