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
model, design and displacement fields. Their spent work remains charged when
an invocation reaches a retained terminal or a later checkpoint.

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

## Durable stages and resume

Each **completed continuation stage** is now committed in one ledger transaction
before any geometry or solver work for the next stage. Its receipt binds the
full octree background, exact raw/projected density bits, native displacement
bits for every load, physical exported fields, complete stage history and
installed refinement evidence. The new receipt consumes the preceding receipt
as actual input lineage. A failed storage operation stops execution immediately;
there is no success record or subsequent refinement after a failed checkpoint.

After each commit, stderr receives a `frankensim.cli.sdf3-progress.v1` JSON line
with its usable `run_id` and `stages_completed`. stdout remains the final command
record. A checkpoint made while more stages are planned says `checkpointed`,
not `completed`. Its report/package can be exported even when later work fails.
Preserve the most recent progress ID to recover after process interruption.

For this producer, `--budget N` limits **new stages in this invocation**, not
individual density updates or Krylov work. `updates-per-stage` in the source
still governs the numerical updates in each stage; neither target is silently
changed. For example, build once and advance the existing three-stage example:

```bash
cargo build -p fs-cli --features sdf3-study
./target/debug/frankensim --json study \
  examples/marquee/bracket-3d-adaptive.fsim /tmp/cantilever-3d.db --budget 1
# Use the run_id printed by that command; exit 6 denotes the retained prefix.
./target/debug/frankensim --json study --resume \
  study-RECEIPT_HASH /tmp/cantilever-3d.db --budget 1
# Use the NEW run_id to continue the remaining stage.
./target/debug/frankensim --json study --resume \
  study-NEW_RECEIPT_HASH /tmp/cantilever-3d.db
```

Resume requires the **same executable bytes** and unchanged original source.
It uses exact, verified prefix replay, not direct optimizer-state restoration:
it reconstructs the original stage sequence, compares the entire accepted-state
manifest hash, then advances. A changed native field bit, different mesh,
density, history or refinement evidence refuses before any new checkpoint.
Previously completed stages are not republished. A completed receipt returns
unchanged without further physics or a new ledger operation.

**Replay is real work and consumes the remaining original allowance.** Recorded
Krylov iterations, geometry boxes/points and wall consumption are subtracted
before replay. Every new checkpoint includes the previous consumption plus the
current invocation's replay, rejected proposals and new work. Setup operator
applications, Galerkin products, solve starts and field evaluations are retained
as separate counters rather than mislabeled as outer iterations. `--budget`
does not raise any original allowance. Resume refuses when replay cannot fit;
the prior receipt remains usable for export and is never replaced by an
unverified reconstruction. Plan enough total work for replay when splitting
an expensive study into invocations.

A first stage stopped before completion still exports its honest numerical
partial, but has no resumable whole-stage checkpoint. Older receipts without
the new manifest remain exportable and explicitly refuse resume. Rebuilding
the executable can also make resume incompatible; retained-only export is
still available. Starting again from an older receipt creates a separate
continuation branch, not an exactly-once shared resource account.

Durability is stage-granular, not per density update or Krylov iteration.
A crash or a failed replay cannot durably charge work performed after the last
checkpoint. Ledger I/O and individual numerical kernels are indivisible; the
final ledger write cannot include its own elapsed time in the bytes it writes.
The CLI has no OS signal handler; its internal cancellation seam is cooperative.

## Budgets and evidence limits

`updates-per-stage` limits accepted updates in each stage. Every filter,
equilibrium and derivative solve shares the declared cumulative Krylov budget.
Initial, enriched and candidate quadrature share one box/point allowance.
Wall checks run at cooperative geometry and solver boundaries. The memory
setting controls a bounded admission envelope, not measured peak RSS. The
finite model also caps the initial level at 2, background leaves at 2048,
load cases at 4 and stages at 4.

Exit 0 means the bounded schedule completed, **not** that an optimum was found.
Exit 6 retains an honest partial result when a stage/work/tree budget stops the
study, or refuses a resume whose original allowance cannot support replay.
A result stopped before any solved baseline has empty fields. No invalid or
unassessed displacement is published. Gradient or numerical failures exit 4;
observed cancellation exits 130.

These results are discrete numerical evidence. They do not establish a
continuum error bound, moving-boundary optimization, a mesh-independent optimum,
manufacturing suitability or physical validation. Descent applies within a
stage, with a fresh baseline for each changed mesh/material/projection model.
