# Feasibility-first thermal population fitting

`examples/scan_power_fit` now has a reusable `scan_power_fit::population`
driver and the `constrained_scan_fit` example. They connect the existing
`ValueConstrainedScanStudy` to `ConstrainedDifferentialEvolutionState`; neither
implements another thermal solver or evolutionary algorithm. The production
workspace dependency graph and the consumer's existing dependencies are unchanged.

## Execute the synthetic specimen

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit
cargo run --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit -- 45 0.5 2400
cargo run --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit -- 60 0.24 2400
cargo run --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit -- --workers 8
cargo run --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit -- --workers 4 45 0.5 2400
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --lib
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --example constrained_scan_fit
cargo test --manifest-path examples/scan_power_fit/Cargo.toml --test spatial_design value::batch
```

Arguments are a total absorbed-energy cap in joules, a dimensionless sampled
liquid-mass-fraction envelope cap, and lifetime search attempts. Defaults are
60 J, 0.5, and 2400 attempts. The run always reserves one additional physical
evaluation for replay. Parameters outside the admitted ranges refuse before
running the specimen.

`--workers 1..128` explicitly selects production TilePool execution. Omitting it
retains serial evaluation. A wave contains at most 24 candidates and stops at
the original generation boundary. The requested worker count and execution mode
are printed; these are configuration values, not measured parallel speedups.

The synthetic unit cube starts uniformly at 80 J/kg (308 K); its solidus is
100 J/kg (310 K), with latent heat 100 J/kg and reference density 1 kg/m^3.
The declared source footprint gives equal weights to the eight P1 vertices.
Targets are volume-mean liquid fractions of 0.1 and 0.3 after one and two
seconds. With this uniform unit-density specimen, the analytic target schedule
is 30 W then 20 W. That equation is an independent reference, not the program's
forward model: every candidate is evaluated by the original P1 enthalpy solver.

Five inequalities are checked on every successful trajectory: the two automatic
nonnegative-power rows, total authored energy, a 311 K first-endpoint temperature
cap, and a liquid-fraction envelope over every vertex at both endpoints.
The envelope has smoothing scale 0.001. Its actual sampled maximum and smoothing
gap are printed separately. The energy integral includes the full declared
source schedule; the final shutdown value remains exactly zero.

The 45 J case cannot meet the final 0.3 fraction target: only 25 J remain for
latent heat after the initial 20 J of sensible heating. Its analytic constrained
best fit is 30 W then 15 W, with fractions 0.1 and 0.25 and objective 0.125.
The executable reports the feasible mismatch and exits unsuccessfully rather
than calling it target recovery. Exhausting a heuristic search does not prove
that this is an optimum; the analytic reference is specific to this specimen.

## Driver semantics

`PopulationScanStudy::new` consumes a value-constrained physical study plus a
finite search box in its optimizer coordinates. It checks the exact constraint
count and both source-coordinate box corners without calling the model factory.
A one-watt reference scale with nonnegative-linear power coordinates makes those
coordinates watts; other coordinate maps retain their original meanings.

Each `advance` passes the original objective and every scaled inequality to the
existing feasibility-first ordering. All feasible candidates outrank infeasible
candidates. Among feasible candidates objective decides; among infeasible ones
maximum violation decides, then objective. No tunable penalty weights are used.
Simulation failures are returned as their original typed errors. Successful
infeasible simulations are still useful evaluated candidates.

The driver owns the requirements and search together; it exposes no mutable
handle that can replace the requirements during continuation. The caller's
physical factory must likewise keep non-design model data fixed. Runtime
source/readout identity checks are retained, but cannot prove arbitrary callback
semantics.

`extend_evaluation_budgets` checks both lifetime allowances before changing
either. Search attempts and physical reservations are distinct: a malformed
source input can fail before a physical reservation. Existing probes are not
refunded. `advance` always leaves one physical slot for replay. Extra replays can
reduce future search capacity, which is visible in the returned work report.

`replay_best` re-runs the best-ranked design, including every physical limit,
and verifies its objective, maximum violation and feasible flag against the
stored score before returning source declarations. It can return an infeasible
best-so-far; neither export nor replay implies feasibility. Its
`target_recovered()` requires both feasibility and the explicit objective target.

## Shared model, independent parallel trajectories

`PopulationScanStudy::advance_parallel(cx, pool, lease, additional, max_wave)`
uses the same owner and requirement set as `advance`. It additionally requires
a `Fn + Sync` factory and `Sync` observations. Mutable-only factories remain
valid for serial execution; they are not cloned, globally locked, or implicitly
shared. Every worker builds its own original physical run from its worker Cx.
There is no shared mutable trajectory or global physical-solver mutex.

The library also exposes `ScanPowerStudy::value_batch` and
`ValueConstrainedScanStudy::value_batch`. A batch exclusively borrows the original
study, sharing only immutable declarations. An atomic reservation check bounds
all worker starts by one allowance. The allowance is reserved in the lifetime
counter before dispatch; after worker drain, dropping the guard returns only
unstarted slots. Failed rebuilds and solves remain charged. Unwinding settles
the reservation too; forgetting a guard spends its full allowance rather than
refunding work. Serial and batch paths share the same forward producer and
constraint-row collector, including robust losses and explicit peak smoothing.

The driver limits each batch by both remaining search attempts and remaining
physical work, keeping one physical slot for serial replay. It settles actual
physical reservations before returning either a report or a typed error. The
underlying DE implementation keeps successful siblings cached, so serial or
parallel continuation retries only missing/failed results. Model inputs and
constraint scales cannot be replaced during a batch or between these modes.
Arbitrary interior-mutability in a caller's factory can still violate model
determinism; `Sync` alone does not prove numerical purity.

With deterministic physical callbacks and no early target stop, equal work
prefixes retain the original serial candidate order and population. Target
stopping occurs between waves; parallel execution may therefore evaluate more
candidates than serial execution before returning a recovered target. Every
started evaluation is charged, including this extra work. Successful replay
still has to reproduce the selected score bit-for-bit.

The example supplies a 64 MiB executor lease, separately from the physical
solver's spatial/work budgets. This does not cap total callback-owned heap or
process memory. Increasing workers increases simultaneous physical workspace.
No new snapshot format, KKT certificate, global optimum, continuous-peak coverage
or experimental calibration is added. The example is not a hardware controller.

## Verification boundary

Two driver budget regressions and four physical example regressions are authored.
Native Rust execution remains unverified. In the publication session, process
launches for DSR, RCH, Cargo and rustc were refused with PermissionError; the
recorded DSR fallback was absent. No compiler or native test pass is claimed.
The independent Python reference was rerun on three constrained cases and two
negative controls. It uses analytic uniform-P1 energy/phase balances and SciPy
constrained DE, not the Rust implementation or Rust concurrency.

Seven additional parallel/batch regressions are authored: serial/batch physical
value and constraint parity under robust loss; shared budget contention; failure,
cancellation, unwinding and forgotten-guard accounting; failed-member retry with
retained siblings and reserved replay; pre-dispatch refusal; one/four-worker and
pause-partition equivalence at fixed work; and CLI admission. They remain
unexecuted. In the parallel-integration session DSR, its documented fallback,
RCH, Cargo and rustc were absent. No native compilation, race-freedom test pass,
physical runtime pass, or wall-time improvement is claimed by this document.
