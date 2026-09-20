# Same-material multi-load elasticity studies

The `--projected` mode of the existing `fs-marquee-elasticity-robust` binary
combines the multi-load evolution kernel with numerical material-area equality.
It does not replace the unflagged whole-trajectory mode or introduce a new
elasticity, sensitivity, quadrature or advection implementation.

```sh
printf 'right,0.375,0.625,0,-1,0.7\nright,0.375,0.625,0.5,0,0.3\n' > loads.csv
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected /tmp/multiload-study loads.csv 4 30 0.45 6 worst 362
```

Positional arguments after `--projected`: new output directory, load CSV, grid
level (default 4), accepted-update target (30), material area (0.45), candidates
per update (6), aggregate (`sum` or `worst`, default `worst`), total case-solve
allowance (default `cases * (1 + updates * candidates)`), and optional initial
level-set CSV. The last example admits two baseline solves plus up to 360
candidate case solves. Each independent solve retains the existing 60,000-CG-
iteration cap. This allowance is not a bound on wall time, allocations, or all
floating-point work; projections separately admit at most 64 area evaluations
per candidate and sensitivities reuse the existing bounded kernel.

## Geometry and loads

The CSV is the existing `edge,start,end,fx,fy,weight` format. Each row is an
independent physical loading experiment; forces are not added before solving.
Weights are not normalized or interpreted as probabilities. Zero-weight rows
still require a successful solve and valid support. `sum` may trade off
individual cases; `worst` minimizes the largest WEIGHTED compliance. Neither
mode adds a stress limit by default or a per-case compliance constraint. An
explicit sample-scoped stress constraint is available as described below.

Coordinates, traction, Young's modulus and compliance use the existing
normalized unit-square plane-strain model (`E=1`, `nu=0.3`), not implicit SI units.
The left edge is clamped. All four boundary nodal traces are frozen from the
initial field, so optimization cannot shorten its authored load-support trace.
The default strip supports mid-height right-edge loads; other edge loads need
an initial field containing their full supported segments. Unsupported loads
refuse rather than being dropped or silently moved.

The optional initial field uses the existing exporter header
`x_normalized,y_normalized,phi_normalized`, followed by exactly `(2^level+1)^2`
rows in x-fastest order. Coordinates must match the declared lattice and values
must be finite. Missing, extra, reordered or non-finite samples refuse. Input
reads are bounded to 8 MiB for fields and 1 MiB for load cases. An exported field
can warm-start a NEW study; it is not an exact checkpoint of the previous
multiplier, budget or global nucleation schedule.

## Acceptance and retained results

First, numerical area is projected to the target (absolute tolerance `1e-4`),
then every baseline load is solved. Improvement is measured against this
feasible baseline, not against the initial overfilled strip. For each update,
the existing kernel constructs simultaneous-load sensitivities. Each proposal
is projected, then solved under EVERY load before the area and strict aggregate-
decrease gates are evaluated. The worst case is recomputed on the candidate;
it is not inherited from the previous design. Rejections contract interface
travel and retry from the same accepted state. Numerical refusals remain visible.

`input-level-set.csv`, `baseline-level-set.csv`, `level-set.csv` and
`load-cases.csv` retain actual geometry and loads. `trajectory.jsonl` records
accepted before/after states; `attempts.jsonl` records candidate projections,
complete-family results and refusals. Proposal drift and nucleation diagnostics
are explicitly proposal-only: projection can change that geometry afterward.
The last accepted field is exported even when the solver budget is exhausted,
the bounded candidate search stalls, or a later sensitivity calculation refuses.
`summary.json` is written last; write failure returns an error, never a printed
success. This is not a crash-atomic multi-file transaction or a design-ledger
package. Existing output directories are refused without overwrite.

Exit 0 / `iteration_limit` means only that the requested update count was reached.
Exit 11 / `no_descent` means the candidate allowance found no acceptable update.
Exit 13 / `solve_budget` means another complete scenario family cannot be started.
Exit 12 / `refused` preserves the last valid field after a later numerical error.
None is a KKT or global-optimality claim.

## Library and verification boundary

`fs_topols::robust_descent::MultiLoadProjectedOptimizer` retains the accepted
field, complete independent solutions, ordinal and search multiplier between
updates. `advance_one_controlled` checks before sensitivities and proposals,
at each projection boundary, before/after each case solve, and before publication.
Cancellation leaves scientific state unchanged, but already attempted solves
remain charged. A failed or interrupted final case cannot publish a favorable
partial-family objective. Baseline construction and each individual solve,
smoothing pass and area evaluation remain synchronous/non-preemptible.

Focused commands:

```sh
cargo test -p fs-topols --lib robust_descent
cargo test -p fs-marquee --bin fs-marquee-elasticity-robust --test projected_multiload
```

The patch adds nine kernel/controller regressions, two field-input unit tests,
and five real-binary integration cases, including non-vacuous accepted-update
checks. They were authored but NOT executed in the implementation environment:
Cargo, rustc, rustfmt, DSR and RCH were unavailable. Patch/source identity,
lexical, delimiter, whitespace and format-string checks are not a Rust build or
PDE result. All metrics remain numerical estimates; no continuum-volume
certificate, experimental validation, 3-D result or full Journey B closure is
claimed.

## Optional sampled-stress constraint

Pass `--stress-limit MAX` and optionally `--stress-tolerance ABS` anywhere after
`--projected`. The limit must be positive and finite; tolerance must be finite
and nonnegative, with a finite sum. Tolerance without a limit, repeated options
and unknown options refuse before input-file reads or physics. For example:

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected /tmp/stress-limited-study loads.csv 4 30 0.45 6 worst 362 \
  --stress-limit 100 --stress-tolerance 0.01
```

The illustrative bound is in the SAME normalized stress units as the existing
`E=1`, `nu=0.3` model. It is not a material allowable in pascals, a calibration,
or a physical safety recommendation. The caller supplies a meaningful bound.

The area-projected baseline must satisfy the stress constraint before the output
directory is created. Every accepted candidate must preserve both area and
sampled-stress feasibility while strictly reducing the selected compliance
aggregate. All cases participate in stress admission, including zero-weight
cases and cases that do not govern worst-weighted compliance. An overstressed
baseline is refused: no stress-feasibility restoration or stress adjoint is
implemented. The existing compliance-generated search may stall at a stress
boundary; this is `no_descent`, never an optimality claim.

`projected-multiload-stress-v1` summaries add `baseline_stress` and `final_stress`;
accepted and attempted rows add `sampled_stress`. These carry the exact limit
and tolerance, every case's sampled maximum, first maximum location and sample
count, the governing case, and the matching geometry snapshot. Missing samples
are `unavailable`, never feasible. Runs without these options keep the original
`projected-multiload-v1` output and do not perform or claim stress assessment.

Stress uses the existing robust sampler over four Gauss probes plus the centre
of full material cells and retained positive-weight bulk/interface points of
positive-volume cut cells. Gradients use each quadrature cell's own material
trace; interface-only exterior neighbours are not material sampling owners.
Missing owning-cell displacement refuses instead of silently omitting a probe.
This is a numerical sample maximum, not a bound on stress between probes or at
singularities. The report explicitly denies a continuum-maximum certificate and
validation of the user-supplied physical allowable.

The library builder `with_sampled_stress_limit` uses already solved displacement
fields and installs the immutable constraint before candidate solves/updates.
Sampling adds no hidden PDE solves. `advance_one_controlled` additionally polls
at every stress cell; cancellation preserves prior accepted geometry and stress
evidence while leaving attempted solves charged. Baseline sampling is synchronous.

Eight added core regressions and six added CLI regressions cover these paths,
including actual acceptance and independent final re-solves. They remain
unexecuted in the implementation environment: Cargo, rustc, rustfmt, DSR and RCH
are absent. Independent affine-Q1/principal-stress reference controls are not
native Rust or CutFEM execution, and do not establish continuum validity.
