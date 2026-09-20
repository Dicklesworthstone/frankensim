# Refine a projected multi-load design

`fs-marquee-elasticity-robust --projected --refine` takes a retained checkpoint
and starts a NEW study on the next finer uniform grid. It uses the existing
CutFEM elasticity, projected compliance descent, optional sampled-stress gate,
exports and checkpoints. `--resume` remains exact same-grid continuation;
`--refine` is deliberately not an alias for it.

For a saved two-case level-3 study:

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected --refine /tmp/coarse/checkpoint.fscp /tmp/fine \
  --updates 30 --max-solves 362 --recovery-solves 4
```

This creates a level-4 study. Repeating the operation on its final checkpoint
provides a coarse-to-fine sequence. Add `--pause-after 0` to retain only the
fine-grid baseline, or `--pause-after N` to pause after N new accepted updates.
Fine-grid checkpoints use the existing `--resume` command without special flags.
The executable permits coarse levels 2–6 (fine levels 3–7), 1–200 new updates,
and the inherited 1–16 candidate limit. New solve and update budgets are required;
missing, repeated, unknown, or physical-policy-replacement options refuse.

## What transfers, and what restarts

The nodal SDF is prolonged bilinearly, not thresholded or re-meshed. Every old
node is injected bitwise. New edge and cell values use fixed-order, overflow-safe
midpoints. Fixed coarse nodes persist. A new node is fixed only if ALL coarse
nodes with nonzero interpolation weights were fixed: fixed edges and interior
regions survive; an isolated prescribed node does not freeze adjacent cells.

Load cases and order, objective aggregation, material parameters, area/stress
limits, and candidate controls are inherited. The load-support policy and other
controls measured in cells keep their declared cell values; doubling resolution
halves their physical cell length. Finer quadrature may measure a different
material area. The original area projection restores feasibility, and EVERY load
is solved again on the resulting fine field before it becomes the new baseline.
An overstressed finer baseline refuses; no allowable is enlarged to admit it.

The new study resets its accepted-update and nucleation schedule and initializes
the AL search state from the inherited `ell0`. Its fine-grid baseline is distinct
from the coarse endpoint. All claimed decreases are fine-baseline/fine-candidate
comparisons, never coarse-versus-fine improvements. No error or convergence bound
is inferred from the two-grid difference.

## Work, files and failures

The source checkpoint is hash-checked and numerically recovered using the existing
recovery implementation (two solves per load case). `--recovery-solves` funds
that work separately. `--max-solves` funds the NEW study, including one new
baseline solve per load. Old spent work remains in the source checkpoint and is
also recorded in the handoff. An exhausted coarse study can be refined only by
explicitly declaring new work; exact `--resume` still cannot replenish its budget.
Failed fine initialization can spend solve work without returning an optimizer.
Initialization and individual physics calls remain synchronous, not preemptible.

`refinement.json` and the summary's `refinement` object record the source payload
hash, source work/ordinal, coarse/fine levels, coarse endpoint, fine-quadrature
area before projection, maximum nodal field correction, and fine baseline.
The correction is not a geometric distance bound. `input-level-set.csv` is the
actual prolonged field; `baseline-level-set.csv` is the area-restored, re-solved
fine baseline. Other output and checkpoint files retain their existing meanings.
Source files and existing destinations are never overwritten. A later write or
solver failure leaves earlier files but prints no false successful completion;
keep `refinement.json` for ancestry when subsequently resuming a fine checkpoint.

## Verification boundary

Six core tests and four CLI tests were added for exact-node and bilinear transfer,
fixed support, refusal/extrema, both aggregates, fine re-solves, stress inheritance,
actual fine-grid accepted updates, pause/resume and new-versus-old work budgets.
They were NOT executed in the implementation environment: Cargo, rustc, rustfmt,
DSR and RCH were unavailable. Independent Python interpolation checks are not
native Rust or elasticity proof. Focused commands in a configured checkout:

```sh
cargo test -p fs-topols --lib refinement
cargo test -p fs-marquee --bin fs-marquee-elasticity-robust --test projected_resume
```

This is uniform refinement, not adaptive DWR marking, a stress adjoint, a
continuum-volume/stress certificate, physical validation or full Journey B closure.
