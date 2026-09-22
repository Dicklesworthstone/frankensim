# One native design, multiple independent operating conditions

`bracket-multi-load-2d.fsim` uses the existing `frankensim study` command to
optimize one geometry for separate downward and horizontal traction scenarios.
It uses the existing simultaneous multi-load CutFEM/level-set optimizer, hard
numerical area projection and per-case sampled-stress admission. It does not
sum forces before solving or optimize a different geometry for each case.

```sh
cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-multi-load-2d.fsim ./multi-load.db --budget 1
# An accepted first update returns exit 6 with a budget-exhausted receipt.
# Copy its run_id (study- followed by 64 hex digits) into RUN_ID.
cargo run --release -p fs-cli -- --json study --resume "$RUN_ID" ./multi-load.db
# Use the newly returned run_id for retained-result export.
cargo run --release -p fs-cli -- --json report "$RUN_ID" ./multi-load.db
cargo run --release -p fs-cli -- --json package "$RUN_ID" ./multi-load.db
```

A bounded search may instead return `no-feasible-descent` (exit 4) with its
accepted design and refusal reasons retained. Neither that terminal nor
`completed` (exit 0, requested update count reached with stress feasible) asserts
convergence. A solve allowance too small for another complete family gives
`budget-exhausted` (exit 6), never a partial-family objective.

## Declare the actual operating conditions

In a `projected-stress` optimizer, place this optional field immediately after
`:stress-tolerance-pa` and before optional `:design-regions`:

```lisp
    :load-family (independent
      :aggregate weighted-sum
      :max-solves 512
      :max-recovery-solves 64
      :additional (
        (case :band (0.375 0.625) :traction-pa (0.5 0.0) :weight 1.0)
      ))
```

The primary `scenario` is case zero with weight **one**, unchanged from its
existing downward load declaration. Add 1..=15 right-edge cases, each with an
explicit support interval, nonzero signed traction vector in Pa and finite
nonnegative dimensionless weight. These are separate operating conditions,
not simultaneous loads or probability samples. Opposite vectors cannot cancel
their demands. The complete declared band must remain supported by material.
Explicit restoration also permits an empty `:additional` list, using only the
primary scenario through this same numerical owner; no duplicate load is needed.

`weighted-sum` minimizes the unnormalized sum of weight times compliance;
`worst-weighted-case` minimizes the largest weighted compliance. A weight of
zero removes that case from the objective, **not** from equilibrium or stress
admission. Every case must satisfy the common unweighted sampled-stress limit
before compliance descent. There is no implicit weight normalization or inferred
allowable. Without restoration, an overstressed initial design is refused.

Protected material/void regions remain available, with their original exact
fixed-node prescriptions. Region preparation precedes area projection and
baseline solves; neither regions nor extra loads relax the feasibility gates.
See `PROJECTED_STRESS.md` for those declarations and their discrete scope.
Omitting `:load-family` keeps the original single-load path unchanged.

## Restore stress feasibility before optimizing compliance

`bracket-stress-restoration-2d.fsim` declares a bounded repair phase for an
initial design that may violate the stress limit. It uses only the primary
load; additional cases can be declared normally. Run it through the same native
study/resume/report commands:

```sh
cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-stress-restoration-2d.fsim ./restoration.db --budget 1
```

The explicit opt-in is `:stress-restoration-reduction 0.01` immediately after
`:max-recovery-solves` in the load family, before `:additional`. Its finite value
must lie in `[0,1)`. Let `e = max(0, worst_sampled_stress - admitted_limit)`.
While infeasible, an accepted candidate must reach `e = 0` or reduce excess
strictly below `(1 - reduction) * previous_e`. This example requests more than
one percent excess reduction per still-infeasible accepted update. Zero still
requires strict reduction; it does not admit unchanged stress.

The existing numerical owner proposes changes using the governing stress case's
**unweighted compliance direction**, including a zero-objective-weight case.
Acceptance uses the independently solved complete family's measured stress,
not that direction's prediction. This is a heuristic restoration search, not
a stress adjoint or a guarantee that a feasible design exists or will be found.
Compliance may increase during repair. Once stress is feasible, the ordinary
strict compliance-decrease and stress-preservation rules apply; returning to an
infeasible state is not permitted. Fixed nodes and numerical area remain enforced
through both phases, which share the original total update and solve allowances.

Reports label the input baseline `area_feasible_study_start` and expose
`constraints.stress_restoration`: current feasibility/phase, repair and compliance
update counts, first feasible update, its measured baseline, and remaining stress
excess. `relative_reduction` is **null before stress feasibility**, then measures
compliance improvement only from the first stress-feasible design. The crossing
update establishes that baseline; it is not a compliance-improvement win.

A repair can be accepted yet remain stress-infeasible. Exhausting the total
update count in that state returns `budget-exhausted`, not `completed`; exhausting
the candidate family returns `no-feasible-descent`. Both retain the actual field,
phase and measurements for reporting. Resume preserves the original first-feasible
baseline and phase through the owner's versioned checkpoint. It cannot reset
an exhausted allowance or reclassify an incomplete repair as success. The example's
1 Pa stress limit and artificial material/load scales are not engineering data,
and its bounded search may stall or finish without a feasible design.

## Retained results and recovery budgets

The existing receipt/report contains `constraints.load_family`: ordered case
loads, each case's baseline and accepted compliance/stress/location/sample count,
solve charges and the numerical owner's exact checkpoint. Aggregate compliance
columns refer to the declared objective. The aggregate sampled maximum is the
worst **unweighted** case; it need not be the objective's governing case.
Without restoration, improvement uses the original feasible baseline, not a reset
resume baseline. With restoration, it uses the first stress-feasible design above.

`:max-solves` counts actual baseline and candidate case-solve starts, including
refused/interrupted candidates, and cannot exceed 100000. Resume never refunds
this work. A separate `:max-recovery-solves` allowance (0..=100000) covers the
original-baseline and current-state replay families: 2 times the number of
cases per successful recovery. Zero disables resumed computation. Charges
persist along the returned receipt chain, including failed replay attempts;
follow the newly retained pointer in an error rather than an earlier receipt.
An explicitly selected old receipt remains an immutable branch point, not a
global account-wide compute quota. Successful uninterrupted and split runs can
therefore match design/iteration/checkpoint bytes while their recovery charges
differ. Reports and repeated completed/stalled/solve-exhausted terminals do not
run fresh physics. Resume requires the same executable and only the ledger;
the original input file is not required.

Candidate updates poll inside CG and before stress-cell sampling. **Multi-load
baseline construction and checkpoint recovery currently use synchronous owner
APIs**; they are not interruptible inside those solves. Wall checks bracket
that work: expired initialization publishes no admitted study, and expired
recovery retains only the previously accepted design with its new work charge.
The already available single-load controlled startup/recovery path is unchanged.
Assembly, individual area evaluations and ledger I/O are also non-preemptible.
No hard deadline or measured peak-memory guarantee is claimed.

The strict example's 2 Pa modulus, 1 Pa primary traction and 1e12 Pa stress limit
are a software exercise, not calibrated engineering data. Scope remains the unit
plate with the existing unit-thickness, 2-D plane-strain convention and a left
clamp. No arbitrary boundary conditions, 3-D physics, physical validation,
continuous stress bound or KKT certificate is added. The package remains
structural evidence, including when its retained design is stress-infeasible.

Focused native tests are in the existing `study::elasticity` library filter and
`study_checkpoint_cli` executable target. The new tests require native Rust
execution; source-level checks do not establish acceptance or replay success.
