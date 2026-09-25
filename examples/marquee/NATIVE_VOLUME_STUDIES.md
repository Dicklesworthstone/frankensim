# Native hard-area elasticity studies

`frankensim study` accepts an explicit `:constraint-mode projected-volume`
for the existing 2-D plane-strain level-set study. It uses the existing
`fs-topols::ProjectedOptimizer`, canonical `.fsim` source, design ledger,
`study --resume`, `report`, and `package` commands. It does not require an
artificially large stress limit and does not produce stress measurements.
The separate `projected-stress` and independent-load modes are unchanged.

## Run the supplied numerical example

```sh
cargo run -p fs-cli --bin frankensim -- --json study \
  examples/marquee/bracket-projected-volume-2d.fsim \
  /tmp/projected-volume.db --budget 1
```

Use a fresh ledger path for a fresh example. The one-update invocation returns
`budget-exhausted` with a nonzero exit and a `run_id` identifying its durable
receipt. Use that actual identifier, not a fabricated hash:

```sh
cargo run -p fs-cli --bin frankensim -- --json study \
  --resume "$RUN_ID" /tmp/projected-volume.db
cargo run -p fs-cli --bin frankensim -- --json report \
  "$FINAL_RUN_ID" /tmp/projected-volume.db
cargo run -p fs-cli --bin frankensim -- --json package \
  "$FINAL_RUN_ID" /tmp/projected-volume.db
```

The example is a small software/numerical fixture with declared material and
load scales, not a calibrated bracket or a material allowable. Native execution
is required before claiming that a particular build passes its regressions.

## What is constrained and measured

The driver first projects the declared geometry onto the requested numerical
material-area equality and independently solves that feasible baseline.
Removing excess material to reach that baseline is not counted as an
optimization improvement. Subsequent candidates use the original level-set
proposal engine, area projection, and a separate solve of the projected field.
Only a strict measured compliance decrease at the same area is accepted.
Prescribed material/void regions remain supported; their exact fixed-node
values are reconstructed from the original source during recovery.

The optimizer's additional fields are explicit and canonical, in this order:

```lisp
    :constraint-mode projected-volume
    :area-tolerance-m2 0.0001
    :max-projection-shift 2.0
    :max-area-evaluations 64
    :max-candidates 16
    :contraction 0.5
    :min-relative-improvement 0.00000001
    :cg-poll-iters 1
```

An optional `:design-regions` declaration follows these fields. Stress or
independent-load fields cannot be silently supplied to volume-only mode:
unsupported, missing, duplicate, or reordered fields refuse admission.

The retained `constraints` object uses `mode: projected-volume-v1`, keeps the
same-material/same-load baseline and every accepted measurement, and reports
actual candidate counts and terminal refusal reasons. `stress_evaluation` is
`not-requested`; there is no invented stress value, sample count, or stress
allowance. The HTML and JSON reports project these same retained results.
The package remains the existing structural evidence envelope, not a newly
certified mechanics result.

## Stopping and recovery

`completed` means the declared accepted-update count was reached while
preserving the numerical area gate. It is not a convergence, KKT, or optimum
claim. `no-feasible-descent` is a bounded candidate-search refusal, not proof
that no improving geometry exists. `cancelled` and `budget-exhausted` retain
the last accepted feasible field; interrupted trial geometry is discarded.
Stopping during initial projection/solve publishes no unmeasured feasible
baseline. CG and phase boundaries are cooperative, not a latency guarantee
inside assembly, quadrature, or ledger I/O.

A resumed segment re-admits the exact endpoint through the existing owner;
it does not replay every previous geometry update or re-project the accepted
field. Recovery time is added to the retained lifetime wall charge. The
original baseline stays in the history. Direct continuation requires the
identical executable, original source identity, field bits, and measured
history. Completed and no-descent receipts return unchanged without solving.

The older augmented-Lagrangian mode is still available when no constraint
mode is declared. Its completion now requires the final evaluated material
area to satisfy the same one-percent-of-box tolerance shown in its report.
Ending the requested updates outside that tolerance returns `constraint-unmet`
with a nonzero exit; the actual design, reports, and terminal receipt remain
available. This does not change its trajectory, target, or tolerance.

## Verification boundary

The focused tests are selected by the existing native-study lane:

```sh
cargo test -p fs-cli --lib study::elasticity::continuation
cargo test -p fs-cli --test study_checkpoint_cli
```

They exercise actual projection and elasticity, exact split/resume histories,
cancellation before publication, schema and retained-state mutations, and the
public study/resume/report/package entry points. Area is numerical cut
quadrature. No continuous-volume/stress certificate, material calibration,
physical validation, 3-D topology result, or optimality is asserted here.
