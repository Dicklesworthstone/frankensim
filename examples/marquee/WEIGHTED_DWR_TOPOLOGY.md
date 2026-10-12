# Weighted-load DWR checks in native structural studies

`bracket-weighted-dwr-2d.fsim` adds signed dual-weighted-residual (DWR) gain
screening to the existing stress-constrained, same-material study. It uses the
existing independent-load optimizer, CutFEM solves, stress sampler, checkpoint,
ledger, report and package paths. No new optimizer or physics implementation is
introduced. This extends the volume-only route in `DWR_CHECKED_TOPOLOGY.md` to
weighted independent loads with the original sampled-stress constraint.

```sh
cargo run --locked -p fs-cli -- --json study \
  examples/marquee/bracket-weighted-dwr-2d.fsim weighted-dwr.db --budget 1

# Use the actual run_id returned by study and the identical executable.
cargo run --locked -p fs-cli -- --json study --resume study-<receipt-hash> weighted-dwr.db
cargo run --locked -p fs-cli -- --json report study-<receipt-hash> weighted-dwr.db
cargo run --locked -p fs-cli -- --json package study-<receipt-hash> weighted-dwr.db
```

The declaration belongs INSIDE the existing load family, before `:additional`:

```lisp
:load-family (independent
  :aggregate weighted-sum
  :max-solves 512
  :max-recovery-solves 64
  :goal-check (elasticity-dwr :max-error-fraction 0.5)
  :additional (
    (case :band (0.375 0.625) :traction-pa (0.0 1.0) :weight 1.0)
  ))
```

The primary scenario retains weight one. Each additional case is a separate
operating condition: opposite forces do not cancel before equilibrium. Weights
are used exactly as declared, not normalized and not interpreted as
probabilities. Zero-weight cases still have real equilibrium and DWR solves;
they remain subject to the unweighted sampled-stress constraint.

An empty `:additional ()` is permitted with `:goal-check`, providing the same
functionality for a single-load stress-constrained study without inventing an
extra load. Without a goal check, existing declarations and canonical identities
are unchanged. The example's large stress allowance deliberately isolates the
DWR behavior; it is not evidence of a binding stress-limited optimum.

## What controls acceptance

The original numerical-area, coarse objective, and all-load sampled-stress gates
run first. For each otherwise admissible candidate, the DWR producer then solves
every declared case on the coarse and enriched grid using the same unchanged
bilinear geometry. Every coarse response must reproduce that exact candidate's
geometry, area and per-case compliance, not just its aggregate. Coarse and
enriched recomputed residuals must be finite and below `1e-12`.

Let `g0_i` and `g1_i` denote each case's coarse and enriched compliance reductions,
and `G0 = sum(w_i*g0_i)`. The weighted objective must improve by the requested
relative amount on both grids. In addition, each of the following must fit
`max-error-fraction * G0`:

```text
abs(aggregate signed-DWR correction change)
abs(aggregate enriched gain - aggregate coarse gain)
sum_i w_i * abs(signed-DWR correction change for case i)
sum_i w_i * abs(g1_i - g0_i)
```

The last two checks prevent opposite error changes in independent cases from
cancelling into an apparent pass. Individual case objectives may trade off when
the declared weighted objective improves; unweighted stress constraints may not.
A rejected candidate contracts the original bounded search. No candidate
replaces accepted geometry, stress, iteration or multiplier until all gates pass.

## Work, recovery and reports

All in-loop DWR work consumes the ORIGINAL lifetime `:max-solves` allowance. For
N cases, the initial baseline uses N original solves; a DWR baseline needs 2*N
slots, and each proposed candidate requires 3*N slots including its original
N solves. A complete family must fit before it starts.

Each entered DWR case conservatively charges two solve slots before its
indivisible estimator. Failed cases can finish fewer solves; those slots are not
refunded. With this policy, `solves_started` is a conservative work debit, not a
measured completed-solve count. The separately declared recovery allowance is
unchanged. Cooperative cancellation and returned numerical failures retain the
accepted geometry and the updated lifetime counter through the existing native
checkpoint. A retry cannot reset that counter. This does not add crash-exact
accounting for work performed after the last durable checkpoint.

JSON reports and receipts expose `constraints.load_family.goal_resolution`:
per-case load identities, coarse/enriched compliance, signed estimates, residuals,
last complete search, accepted candidate and recalculable gain arithmetic. HTML
shows the matching per-case table. Resume reconstructs arithmetic, verifies
individual coarse-state bindings and checks the source level and work counters
before recovery solves. Finished or exhausted receipts are returned unchanged.
A partial or cancelled DWR family is not published as a complete measurement.

## Scope and validation

This is Estimated numerical screening, not a continuum error bound, stress
certificate, KKT proof, optimum or physical validation. A shared discretization
bias can survive. Material area and stress remain the original coarse-grid
quadrature and deterministic sample set. Absolute indicator mass is not used as
an error radius. There is no spatial adaptation in this increment.

Native admission requires level 2..=5, at most 16 cases, at least 256 MiB memory,
and a finite fraction strictly between zero and one. Worst-weighted objectives,
stress restoration and simultaneous mesh-check policies are refused rather than
interpreted as a differentiable weighted-compliance sum. The separate final DWR
assessment remains available under its own declaration and work scope.

Original CG polling remains available during ordinary candidate solves. A DWR
case's two solves and residual integration remain indivisible between callbacks;
no intra-DWR-CG latency or hard wall-time bound is claimed.

```sh
cargo test --locked -p fs-topols weighted_dwr_ -- --nocapture
cargo test --locked -p fs-cli weighted_dwr_native_ -- --nocapture
```

Six library and five native regression sources cover non-cancelling gain
arithmetic, objective tradeoffs, case binding, zero weights, real independent
solves, lifetime budgets, cancellation, deterministic retry, report projection,
source admission and changed retained fields. Native compilation/execution is
pending: the editing environment has no Rust toolchain. The broader
`frankensim-rc-root-q61wp.16` acceptance remains open.
