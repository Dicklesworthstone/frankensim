# Valves in reacting, molecular and wall-coupled shock networks

Behind `experimental-euler`, the existing finite-plenum valve model now shares
one `NetworkValveOperator` between the plain `ValvedNetwork1d`, molecular
`DissipatingNetwork1d`, and `ReactingNetwork1d` paths. It reuses the original
variable-Cp converging-nozzle, sonic-limit and finite adiabatic-blowdown solvers;
there is no second valve thermodynamic implementation.

## Prepare once, compose each accepted step

Prepare the molecular, wall and valve operators from the same raw
`MixtureNetwork1d`, then move that network into its reacting or molecular owner.
All numerical and receipt storage is prepared once. A source-only
`NetworkValveOperator::advance_for` changes no clock and stages every plenum
before publishing a complete interval.

`NetworkSources` is the same public type re-exported by both `network::reactive`
and `network::dissipating`. It borrows optional wall and valve operators with
separate controls, avoiding a different owner for each combination of physics.
For already prepared `flow`, `molecular`, `walls`, and `valves`:

```rust,ignore
use fs_gas::multicomponent::network::reactive::{NetworkSources, ReactiveStepControl};
use fs_gas::multicomponent::network::valve::ValveSolve;
use fs_gas::multicomponent::wall_exchange::WallControl;
use fs_gas::multicomponent::dissipation::DissipationControl;
use fs_qty::Time;

let report = flow.step_with_dissipation_and_sources(
    &mut molecular,
    Time::new(maximum_duration_s),
    ReactiveStepControl::default(),
    DissipationControl::default(),
    NetworkSources {
        wall: Some((&mut walls, WallControl::default())),
        valves: Some((&mut valves, ValveSolve::default())),
    },
    || false,
)?;
let valve_transfers = valves.accepted_exchanges().unwrap();
assert_eq!(valves.last_report().unwrap().elapsed, report.elapsed);
```

Without molecular transport, the reacting owner's `step_with_sources` takes
that same bundle. The dissipating owner's `step_with_sources` uses its already
owned molecular operator. All pre-existing `step`, `step_with_wall`, and
`step_with_dissipation_and_wall` APIs delegate to the common composition.

The full order is **W/2, V/2, D/2, C/2, A, C/2, D/2, V/2, W/2**; absent sources
are omitted. Valve links are swept in forward then reverse retained order.
Only A advances the clock. A shorter A duration, rejected source, or exhausted
source budget restarts every requested stage from the accepted gas state and
accepted finite-wall temperatures. Final valve or wall changes must still fit
the chemistry temperature interval. A final cancellation/admission point
precedes infallible publication of all gas fields and operator receipts.

## Energy, species and changing openings

Both halves' valve species/enthalpy transfers are retained per link. They are
internal and never enter the external boundary-energy account. Reaction changes
species-reference energy exactly once; no extra jet heat, reaction heat, or
stationary-drag energy is added to total gas energy. Finite-wall storage closes
the gas-plus-wall balance; prescribed-reservoir heat remains external.

Source-only valve residuals are summed separately across the two intervals.
Subtracting the first before-state from the second after-state would wrongly
attribute intervening chemistry, heat or advection to the valves. Complete
owner reports independently audit the actual final gas inventories.

`valves.set_orifice(index, admitted_orifice)` changes an opening exactly between
accepted calls, including zero-area closure and ideal check-valve direction.
Historical receipts continue to describe the preceding accepted interval. The
caller must limit a requested duration to the next scheduled opening-change
time; this API does not silently integrate across an actuation event.

## Checks and limitations

Four source-level tests cover actual shared-tank balances, intervening heat,
late cancellation/replay, zero duration and source/volume refusal. Six public-API
integration tests in `tests/valved_sources.rs` cover full per-plenum balances,
chemical formation energy, inert/reacting-path agreement, heat-driven complete
restarts, every operator's rollback/reuse, and final-valve cooling across a
narrow chemistry interval. Run:

```bash
cargo test -p fs-gas --features experimental-euler --test valved_sources
cargo test -p fs-gas --features experimental-euler multicomponent::network::valve
```

These Rust tests were authored but not compiled or run in the authoring session:
Rust/Cargo/DSR were unavailable. An independent constant-Cp two-tank calculation
checks the narrow-temperature regression's premise; it is not execution of the
production Rust or a native test verdict.

The model remains Estimated: compact rigid plenums, quasi-steady converging
orifices, caller-supplied discharge/opening data, ordered pair splitting and
first-order valve/chemistry/wall heat integration. Molecular transport retains
its constant-area continuum and explicit mixing-length assumptions. This does
not add resolved valve jets, lift/contact dynamics, connection inertia,
compliant-volume work, solid conduction between walls, real-gas phases, a
second-order proof, an entropy theorem, or experimental validation.
