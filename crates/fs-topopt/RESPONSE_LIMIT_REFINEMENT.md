# Minimum-material designs: check limits, then selectively refit

`minimum_material_sdf3` can assess its complete physical limit family rather
than applying the squared-response fitting estimator to a different problem.

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 80 250000 0.0012 0.003 --estimate
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 80 250000 0.0012 0.003 --refit
```

Use the repository's remote execution lane where available. The update cap
applies per optimization phase; the original Krylov and quadrature controls
are shared by all fits, preparations, probes and rejected trials. The estimate
is diagnostic at inherited physical stiffness, not a new optimized design.

## Per-limit assessment

`limits::refinement` exposes reference displacement and reaction observations,
upper/lower bounds, scales and independent experiment laws. Call
`estimate_response_limit_enrichment` with a complete retained
`ResponseLimitEvaluation3`, an enriched backend and these same reference laws.
Every case/row and the original smooth aggregate are revalidated before
preparation. Loads, observations and prescribed motion are reintegrated on each
grid; no nodal vector is prolonged or filled by nearest-neighbor matching.

Each backend prepares once. One new primal per case and a coarse/fine adjoint
pair per original limit use the existing solver control. Even inactive coarse
limits are assessed: the worst constraint may change with refinement. The
reaction's affine boundary offset is retained alongside the existing DWR,
coarse-space, observation-transfer and finite-solve contributions. Each signed
normalized residual has its own reconstruction. The nonlinear smooth-maximum
change is reported directly, never mislabeled as a sum of linear DWR terms.

Reports keep every original value/bound/residual, the enriched smooth aggregate,
its worst original residual, normalized changes and numerical volume. Feasibility
checks include every fine original row AND the conservative smooth aggregate.
Original-grid marks use the sum of normalized absolute local contributions, so
upper/lower bounds and opposite experiments cannot cancel marking. Offset and
measure changes remain separate: empty marks are not an accuracy pass.

## Selective reoptimization

`fit_reference_minimum_material` binds reference laws to one grid and calls the
existing `MinimumMaterialStudy3`. `refit_reference_minimum_material` consumes an
explicit candidate study while borrowing its source immutably. It transfers
raw densities only, inherits protected region labels through the geometric
parent map, rebuilds the candidate's actual filter/projection and solves a new
baseline. The new optimizer does not inherit coarse fields, adjoints, gradients,
multiplier or spectral-step state.

A force/displacement window is not monotone in material density. There is no
uniform scaling that generally restores all those constraints, and minimizing
material does not introduce a separate volume cap. The candidate may therefore
start infeasible. Its report retains that fact, the accepted prefix, the least-
volume feasible incumbent (possibly absent), and the actual optimizer outcome.
A late failure does not mutate the source. `into_best_feasible` explicitly selects
matching candidate fields/scales; callers must still handle its optimizer error.

The example's `--refit` marks at most two original cells, selectively rebuilds
that tree (plus balance closure), and prints new-baseline and fitting evidence
before propagating late errors. It installs a candidate only after a successful
optimizer return and selection of a feasible accepted endpoint. Iteration limits
are not stationarity; coarse-to-new-baseline volume is not claimed to descend.
The final refitted design would need its own further error probe.

These are numerical two-grid diagnostics, not continuum constraints, saturation
bounds, follower loads or shape derivatives. Families and per-limit reports are
bounded by their existing case/limit/geometry allowances; more limits require
more adjoint solves and localized evidence. Native Rust validation is separate
from the independent algebraic development check. Focused tests are
`cargo test -p fs-topopt --features cutfem-marquee --test response_limit_refinement`.
