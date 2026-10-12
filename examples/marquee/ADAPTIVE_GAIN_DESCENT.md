# Local adaptive analysis inside projected topology descent

`MultiLoadProjectedOptimizer::advance_one_adaptive_dwr_polling` checks a
proposed design **before acceptance**, using actual locally refined elasticity
analysis. It keeps the existing proposal, numerical-area projection, independent
load solves, sampled-stress requirements, and atomic accepted-state update.
It does not introduce another optimizer or a finite-difference physics path.

```sh
cargo run --release -p fs-topols --example adaptive_descent
cargo run --release -p fs-topols --example adaptive_descent -- --strict
cargo test --release -p fs-topols --lib projected_adaptive_ -- --nocapture
```

The example starts from a real area-projected strip with two opposite operating
loads. The larger load has zero objective weight but still participates in the
original unweighted stress constraint and independent analysis. Material values,
large stress limit, and loose default endpoint-resolution allowance are software
fixtures, not calibrated engineering allowables. `--strict` requests a deliberately
unattainable endpoint tolerance within the two-round limit to exercise actual
local refinement and honest refusal. Neither mode promises acceptance or an optimum.

## What changes during an update

For each otherwise admissible coarse candidate, the accepted baseline and
candidate are independently solved on the **same** analysis grid and its one-level
enrichment. Both geometries stay fixed. A joint indicator map combines weighted
absolute cell contributions from both designs, so refinement does not favor only
one endpoint. The existing bulk marking and budgeted 2:1 closure select new cells.
The original design lattice and proposal sensitivity grid are not refined.

Every accepted comparison must satisfy all of these numerical screens: both
endpoint resolution tolerances; requested relative objective decrease on both
grids; non-cancelling per-case correction/gain changes within the specified
fraction of the coarse gain; and original target/pairwise material-area agreement
on both coarse and enriched grids. The final analysis does not replace the
optimizer's original objective, area, or stress history.

The first round must reproduce the snapshot, area, aggregate and **every**
original unweighted case compliance bit for bit, including zero-weight cases.
A failed numerical solve or unresolved comparison contracts the existing search.
Cell, round and level limits return the last complete pair, not a partial grid.
A refused bounded search is not a stationarity or optimality result.

## Bounds, cancellation, and continuation

With N cases and R maximum rounds, a complete candidate requires allowance for
N original equilibria plus 4*N*R adaptive solve slots. This worst-case family must
fit the original lifetime `max_solves` before proposal work. Only entered cases
are charged; each entered two-solve estimator conservatively spends two slots,
including a later numerical failure. Cancellation never refunds spent work.
An unfunded family starts no proposal or analysis callback.

The unchanged owner checkpoint retains accepted geometry, original independent
solutions through verified recovery, global ordinal, multiplier, stress and work
counter. The caller must retain `AdaptiveGainPolicy` alongside it and select the
same adaptive entry point on resume. The policy is not silently installed into
all owner methods. Cancellation polls bracket cases, rounds, refinement, and
publication as well as the original CG boundaries. A two-solve DWR case and its
integration are still indivisible; there is no hard latency guarantee.

This route supports weighted sums with one to sixteen cases, initial levels
1..=5, maximum analysis level 6, up to sixteen rounds and 4096 coarse leaves.
An estimator enriches once beyond those coarse leaves. Worst-weighted objectives
and stress restoration are rejected, not reinterpreted as weighted sums.

## Scope

Results remain **Estimated**. Shared bias can survive these screens; signed
DWR evidence is not a continuum-error bound. Original sampled stress is retained
but is neither re-evaluated nor certified on the adaptive grids. There is no
adaptive design-space optimization or new stress-adjoint method here.

This increment exposes the in-loop route through the Rust owner and runnable
example. Native `.fsim` declarations still use their existing fixed-grid gain
checks and optional adaptive **final** assessment; they do not automatically
select this new in-loop entry point. Native execution of the new regressions and
example must be checked separately from source review.
