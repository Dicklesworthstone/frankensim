# Adaptive 3-D topology studies

Run the actual raw-implicit CutFEM density optimizer through `frankensim study`:

```bash
cargo run -p fs-cli --features sdf3-study -- --json study \
  examples/marquee/bracket-3d-adaptive.fsim /tmp/cantilever-3d.db
```

This frontier feature uses the existing `fs-topopt` adaptive continuation driver.
The separate `fsim-sdf3-study :version 1` input preserves the existing 2-D
free-boundary study format. Every section and field shown in the example is
required, in the shown order; unknown and repeated fields refuse before work.
All quantities use SI units. The example's modulus is **1 Pa**, deliberately a
numerical mechanics example rather than a physical bracket material card.

The domain is the part of the explicit one-metre cube below
`z = height-m + curvature-per-m * x * (1 m - x)`. This scalar field is an
implicit domain function, not an exact signed distance. Its boundary stays
fixed. The clamp is `x = 0`; each declared vector is an independent constant
reference body-force density in N/m³. Weights multiply independently solved
compliances, without normalization or summing the load vectors together.
Change the material, load vectors/weights, volume cap, physical filter radius,
raw-density start and continuation schedule in the input.

Each later `(penal beta)` stage first performs an enriched solve with the
accepted physical stiffness. Actual compliance goal residuals select octree
cells for refinement. The proposed background inherits raw densities, rebuilds
the same-radius filter, restores projected-volume feasibility and passes
compliance **and** volume directional derivative checks before optimization.
Failed or interrupted refined stages preserve the prior accepted background,
model, design and displacement fields. Their spent work remains charged.

The returned `run_id` is a `study-<receipt-hash>`. Export retained results with:

```bash
cargo run -p fs-cli --features sdf3-study -- --json report \
  study-RECEIPT_HASH /tmp/cantilever-3d.db
cargo run -p fs-cli --features sdf3-study -- --json package \
  study-RECEIPT_HASH /tmp/cantilever-3d.db
```

`report` writes HTML, JSON, a `.design.json` with cell identities, cut volumes,
raw/projected densities, physical node positions/connectivity and one nodal
displacement field per independent load, and `.stages.json` with stage-local
histories and gradient checks. Export does not repeat any physics and works
even in a binary built without the `sdf3-study` feature. The JSON report retains
the goal correction's residual, coarse-space, transfer and algebraic terms,
actual marking fractions, installed refinements and detailed stop reasons.
The package is the shared structural provenance envelope; its checker pass
does not certify the mechanics or turn estimates into error bounds.

`updates-per-stage` limits accepted updates in each stage. Every filter,
equilibrium and derivative solve shares the declared cumulative Krylov budget.
Initial, enriched and candidate quadrature share one box/point allowance.
Wall checks run at cooperative geometry and solver boundaries. The memory
setting controls a bounded admission envelope, not measured peak RSS; one
kernel or ledger write is indivisible. The finite model also caps the initial
level at 2, background leaves at 2048, load cases at 4 and stages at 4.

Exit 0 means the bounded schedule completed, **not** that an optimum was found.
Exit 6 retains an honest partial result when a work/tree budget stops the study.
A result stopped before any solved baseline has empty fields. No invalid or
unassessed displacement is published. Gradient or numerical failures exit 4.
The CLI has no OS signal handler; its internal cancellation seam is cooperative.
Results become durable at invocation completion. `--resume` and `--budget`
overrides explicitly refuse for this producer; use its declared work budgets.

These results are discrete numerical evidence. They do not establish a
continuum error bound, moving-boundary optimization, a mesh-independent optimum,
manufacturing suitability or physical validation. Descent applies within a
stage, with a fresh baseline for each changed mesh/material/projection model.
