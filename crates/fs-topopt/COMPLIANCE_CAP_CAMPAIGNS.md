# Repeated, resumable compliance-cap design

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --adapt-material 1000 1500 3 200 250000 0.02
```

Arguments are authored Y/Z compliance caps, maximum selective refits (0–4 in
this example), accepted updates per grid, one cumulative Krylov allowance, and
maximum coarse/fine change divided by each cap. The last three arguments default
to 200, 250000 and 0.02. Caps are illustrative dimensionless requirements, not
calibrated ratings. The original worst-load, minimum-material and one-refit
routes are unchanged.

`solve_reference_compliance_caps` solves the initial adaptive grid using the
same reference-load integration, feasible-start recovery, exact gradients,
constrained optimizer and fresh recheck used by finer-grid refits. Its result
retains both a restoration seed and a separate optimized incumbent on late
stops. The working study is consumed into this result; initial input/load
assembly errors occur before any usable design is admitted.

`limits::campaign::ComplianceCapCampaign3` takes a completed, rechecked result,
its background, the same reference laws and the original `SolveControl`. It
binds these to one geometry builder and immutable campaign settings. All calls
share the physical budget. The builder must preserve SDF, material, supports
and integration policy, honor the supplied checkpoint, and charge all builds
(including failures) against one cumulative geometry allowance. The example
does so and reserves a refinement level for the final diagnostic probe.

`run(n)` promotes at most n design grids. `run(0)` performs no work. Positive
segments probe their current grid before pausing or reaching the total refit
cap. Completed probes remain cached, so `run(1); run(1)` need not repeat the
intermediate solve performed by `run(2)`. Caps, load order and smoothing do not
change. The diagnostic grid inherits physical eroded stiffness; refits transfer
raw controls and protected regions and rebuild their design map.

A fully solved candidate remains pending until the promotion checkpoint.
Retrying cancellation at that boundary reuses the checked candidate. A failed
candidate remains separately inspectable and never replaces the installed
study. Retrying its build/refit consumes another candidate attempt; failed
attempts and physical work are not refunded. Cancellation after promotion
retains the complete new grid, whose probe is still required before a pass.
This is in-memory segmentation between grid refits, not a durable checkpoint
or continuation inside an interrupted candidate's optimizer.

Only `StableAndFeasible` passes the current coarse/fine checks for every cap,
both smooth envelopes and the separate cap-normalized change tolerance. A
refit count limit, attempt limit, missing marking signal or stopped candidate
is not an accuracy pass. Old operators and displacement families are released
on promotion; summaries and the latest full probe remain available, together
with complete fields for the installed and any pending/rejected candidate.

The example exports the installed densities and independent displacement fields
before propagating late errors; pending and rejected fields are not mislabeled
as installed. It returns nonzero unless the current two-grid check passes.
There is no continuum guarantee, volume discretization bound, stress/buckling
claim, manufactured length scale, global optimum or cross-grid descent claim.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::limits
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 campaign::tests
```

Native tests could not start locally because Cargo/rustc were absent. Static
source checks are not execution of these implementations or production CutFEM.
