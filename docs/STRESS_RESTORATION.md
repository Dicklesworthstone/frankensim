# Restore sampled stress feasibility before compliance optimization

`fs-marquee-elasticity-robust --projected` now accepts an explicitly opted-in
stress-restoration phase. The ordinary `--stress-limit` mode still refuses an
initially overstressed design. Add `--restore-stress` to attempt repairing it
at the SAME numerical material area and with the SAME loads and stress limit.

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected /tmp/restoration loads.csv 4 30 0.45 6 worst 362 \
  --stress-limit 100 --restore-stress --restoration-reduction 0.01 \
  --checkpoint
```

The illustrative bound uses the existing normalized stress units (`E=1`,
`nu=0.3`), not a validated physical allowable in pascals. Supply a limit suitable
for the declared problem; this example does not promise that restoration is
needed or that a feasible design will be found.

## Two different acceptance phases

While the retained field exceeds the sampled limit, the currently governing
stress case supplies an UNWEIGHTED compliance-based design direction. This
reuses the existing energy, Sobolev smoothing, advection, nucleation and area
projection. It is a heuristic proposal, NOT a stress adjoint or stress gradient.
Zero-objective-weight cases can govern this phase. No physical load or objective
weight is changed in the actual independent case solves.

Every projected candidate is solved under EVERY load and sampled for stress.
Only a strict decrease in the complete family's worst sampled stress excess
can be accepted while infeasible. With reduction `r`, the remaining excess must
be strictly less than `(1-r)` times the previous excess, or reach the original
admitted bound. The default `r` is `0.01`; the valid range is `[0,1)`. Reaching
the bound ends restoration. Thereafter each accepted update must preserve
sampled-stress feasibility AND strictly improve the original compliance
aggregate. A proposal that helps one case but increases the family maximum
cannot pass restoration.

Compliance may increase during restoration. Accepted rows explicitly state
`acceptance_phase: "stress_restoration"` or `"compliance"`. The original
area-feasible baseline remains unchanged and its stress is reported honestly;
it is not renamed a stress-feasible baseline or credited as an objective win.
Both phases share the original update goal, candidate limits and total solve
allowance. No hidden elasticity solves or budget refunds are introduced.

## Stops, retained results and continuation

Restoration-enabled summaries use `projected-multiload-restoration-v1`. Their
`stress_restoration` object records the immutable reduction policy, current
phase, restoration-update count and compliance-update count. The original
per-case baseline/final stresses retain measured maxima, locations, counts,
geometry identities and over-limit status.

An update limit reached while still infeasible returns `stress_infeasible`
(exit 15); an exhausted candidate search while infeasible returns
`stress_restoration_stalled` (exit 16). Solve-budget exhaustion remains exit 13,
with the actual phase and remaining violation visible. Deliberate pause remains
exit 14. Reaching an update limit after feasibility returns the existing exit 0;
this still means an iteration limit, not convergence or optimality.

Accepted fields and optional checkpoints remain available after these stops.
Restoration checkpoints use explicit version 2, preserving the original problem,
area-feasible baseline, accepted field, phase-related count, search state and
spent work. The existing resume path independently re-solves baseline/current
and rechecks stress before restoring the exact recorded state. Its separate
recovery allowance does not replenish study work. Resume cannot replace policy.
Refinement inherits the same restoration policy and allowable, but retains its
existing explicitly funded NEW-study semantics. Strict/no-stress studies keep
their version-1 bytes and previous numerical path.

## Verification boundary

Eight core tests and six CLI tests were added: strict policy boundaries,
zero-weight governing loads, partial infeasibility, real accepted restoration
and independent endpoint checks, cancellation, versioned restart, refined policy,
option validation, honest terminal codes, original-path parity and process
restart without a work refund. Native execution is pending: Cargo, rustc,
rustfmt, DSR and RCH were unavailable in the implementation environment.
Independent scalar acceptance-policy checks are not an executed elasticity run.

```sh
cargo test -p fs-topols --lib robust_descent::engine::projected::restoration
cargo test -p fs-marquee --bin fs-marquee-elasticity-robust --test projected_resume --test projected_multiload
```

The compliance-generated direction may stall even when a feasible design exists.
This feature does not establish stress-gradient correctness, guaranteed feasible
recovery, a bound on stress between probes, a physical safety assessment, KKT
convergence, global optimality or full Journey B completion.
