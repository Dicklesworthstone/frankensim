# Adaptive final-design elasticity assessment

The native `study` command can now run the existing goal-driven local CutFEM
refiner after retaining a projected optimizer endpoint. The original bilinear
material geometry is unchanged; only the background analysis grid is refined.

```sh
cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-adaptive-dwr-2d.fsim adaptive.db

cargo run --release -p fs-cli -- --json study \
  examples/marquee/bracket-weighted-adaptive-dwr-2d.fsim weighted-adaptive.db

# Use the returned run_id and the identical executable.
cargo run --release -p fs-cli -- --json report study-<receipt-hash> adaptive.db
cargo run --release -p fs-cli -- --json study --resume study-<receipt-hash> adaptive.db
```

Append the explicit assessment after the existing optimizer block:

```lisp
(assessment :type elasticity-adaptive-dwr
  :max-solves-per-attempt 4
  :max-level 4
  :max-rounds 2
  :max-cells 256
  :absolute-tolerance-j 0.00000000000001
  :relative-tolerance 0.0
  :mark-fraction 0.3)
```

`max-solves-per-attempt` must equal `2 * case_count * max-rounds`. The weighted
example has two independent cases, so it declares eight rather than four.
All cases are solved, including zero-weight cases. Objective weights are neither
normalized nor interpreted as probabilities. Original coarse responses are
retained individually through refinement and checked against the accepted design
before publishing the assessment; an unchanged weighted sum is insufficient.

Native admission supports `projected-volume` and `projected-stress` with an
explicit weighted-sum load family. The existing stress-mode finish path requires
an area/stress-feasible endpoint. Worst-weighted objectives and unprojected
adaptive endpoints are not admitted. Existing fixed two-grid `elasticity-dwr`
assessment declarations keep their previous canonical form and behavior.

## Results and termination

The existing JSON report contains `goal_error_assessment` with method
`adaptive-elasticity-compliance-dwr`. It includes every completed round, the final
mixed-level `grid_cells`, per-case signed cell indicators and residual terms,
initial and final per-case compliances, solve counts, stopping reason and an
explicit `tolerance_satisfied` boolean. HTML contains the corresponding round
table. Original and assessed material areas are reported separately.

The top-level study status describes the optimizer, not mesh accuracy. In
particular, `completed` can coexist with `tolerance_satisfied: false`. These
examples deliberately request a tight tolerance with only two rounds, exposing
this distinction. A complete budget-limited assessment is retained rather than
reported as tolerance success. Increasing its allowance requires a new study
with a new explicit policy, not a silent budget increase on terminal resume.

The optimizer geometry, compliance history and original material-area constraint
are never replaced by finer-grid results. Stress constraints retain their
original coarse-grid sample scope; adaptive compliance assessment does not
certify finer-grid stress. The estimates are not continuum-error bounds,
physical validation, convergence proofs or optimality certificates.

## Work and interruption

The initial analysis level is the study's mesh level (2..=5). `max-level` caps
analysis leaves through level 6; each DWR estimate internally enriches once more.
`max-rounds` is 1..=16, and `max-cells` is at most 4096 and must cover the initial
uniform grid. Cell admission includes all balancing splits before mutation.
At least 256 MiB memory must be admitted; that is not a measured peak-memory
claim. Nonnegative finite tolerances require at least one positive component,
and the marking fraction must be in (0,1].

The existing final-assessment lifecycle retains the finished optimizer before
analysis. Cancellation at a case, round, refinement or publication boundary
leaves that design intact and analysis pending. Resume retries assessment only;
it does not repeat accepted optimizer updates. A cancelled analysis publishes no
partial grid/load family. Completed receipts return unchanged without new solves.

Each case's two solves and residual integration remain indivisible. Assessment
work has an explicit per-attempt allowance and shares the existing cumulative
wall budget; retries can repeat completed analysis work. This is not a mid-solve
checkpoint, crash-exact work ledger, or hard wall-time guarantee.

## Focused checks

```sh
cargo test --release -p fs-cli --lib native_adaptive_ -- --nocapture
cargo test --release -p fs-topols --lib projected_adaptive_dwr_ -- --nocapture
```

Native regression sources cover canonical policy admission, actual local
refinement, normal report/ledger output, read-only terminal resume, cancellation
and deterministic retry, unchanged geometry, full independent-load solve grants,
and an altered zero-weight response that leaves the aggregate unchanged.
