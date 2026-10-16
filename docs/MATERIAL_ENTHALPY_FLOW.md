# Material enthalpy transport, electrical feedback, and flow reservoirs

`fs_thermal::enthalpy::advection::Advected<O>` adds mass-carried thermal energy to
an **existing** `DepositionSimulation` and its source/history observer. It does
not add a caloric law, thermal solver, current solver, or material-name branch.
The compatibility path is
`fs_conduction::transient::deposition::advection` (the same types, not a copy).

## What changes when material data change

The original immutable `ThermalLaw` still maps each cell's specific enthalpy to
temperature; an `EquilibriumEnthalpyPhaseCurve` retains sensible and latent heat
without replacing its plateau with an apparent heat capacity. An electrical or
conductive material resolver receives that same enthalpy and temperature at every
accepted explicit step. It can distinguish different liquid fractions at the
same temperature. The resulting resistivity, current, Joule source and link
conductance act on the same next state that the flow transport advances.

Use `advance_joule_enthalpy_with_material_transport_observed` or the existing
material-interface observed entry point with an `Advected` observer. No separate
advance or after-the-fact temperature correction is needed. Electrical material
payloads, selected-link laws, endpoint support checks and current budgets remain
owned by the existing electrical/material adapter. Other sources such as a
heater can be wrapped directly or combined through `CoupledObservers`.

The regression in `crates/fs-lbm/tests/joule_enthalpy_composed/advection.rs` uses
the actual electrical solve, original phase chart, material-resolved conduction,
and prescribed cold through-flow. Changing only latent heat or phase-dependent
resistivity data changes current while the core remains at its phase-transition
temperature. It is a synthetic implementation test, not measured metal data.

## Prescribed flow data

For a closed loop, create uniquely identified `MassFlow` connections. A positive
`kg_s` flows from `a` to `b`; negative values reverse it. The donor supplies its
current specific enthalpy, not `Cp * T`. At every active cell the incoming and
outgoing mass rates must balance. Different material IDs cannot exchange mass
through this lane: changing composition needs an explicit species/mixture model.

For a well-mixed through-flow volume, declare both an inlet and an outlet:

```rust,ignore
use fs_conduction::transient::deposition::{CellId, MaterialId};
use fs_conduction::transient::deposition::advection::{Advected, MassBoundary};

// `other_physics` is an existing deeply cloneable DepositionObserver (or ()).
let flow = Advected::new(other_physics, vec![], 0)?.with_boundaries(
    vec![
        MassBoundary::Inlet {
            id: 100, cell: CellId(1), kg_s: 0.01,
            material_id: MaterialId(1), specific_h_j_kg: inlet_enthalpy,
        },
        MassBoundary::Outlet { id: 101, cell: CellId(1), kg_s: 0.01 },
    ],
    2, // Declared ceiling across interior AND boundary connections.
)?;
// Pass &flow to the existing observed thermal/electrothermal advance.
```

For a flow-through chain, give the inlet to the first cell, connect the cells with
`MassFlow`, and put the outlet on the last cell. The same conservation check
includes all of these connections. Interior transfer count also respects the
thermal timeline's `max_links` budget; the contacts themselves need not coincide
with flow connections. A flowing face does not silently cover a thermal surface.

Inlet enthalpy is queried against the receiving cell's actual caloric law and
must use the same `MaterialId` and energy datum. There is no silently inferred
inlet phase, `Cp`, mixture, velocity, or pressure. An outlet carries the current
local enthalpy. Rates are constant for an advance; change them at an accepted
endpoint, reconstructing the wrapper (or replacing its boundary data). Select
those endpoints to resolve physical changes, not to conceal discontinuities.

## Conservation and integration

For fixed control-volume masses, the explicit step is

    m_i (h_i_new - h_i_old) / dt
      = conductive_i + ordinary_sources_i
        + sum(incoming mdot * h_donor) - sum(outgoing mdot * h_i)
        + sum(inlet mdot * h_inlet) - sum(outlet mdot * h_i).

Each internal connection contributes the identical power with opposite signs.
External reservoirs enter `EnergyAudit.external_heat_j` exactly once. Their net
contribution is also observed in `AdvectionAudit.external_enthalpy_j`; **do not
add that observation to the thermal audit again**. Joule/work observations still
cover only their own sources. Inlet/outlet mass counters describe throughput:
they neither change stored mass nor appear as deposited birth mass. Gross
absolute transfer is datum-dependent, not dissipation or heat generation.

The existing conduction stability restriction and the outgoing-mass restriction
share one step. If their maximum decay rates are `a` and `b`, the bound is at
most `0.8/(a+b)`, not the unsafe minimum of separately funded limits. Independent
observers tighten the proposal in sequence. The event/target distance clips the
result only after this restriction, avoiding asymptotic non-completion at an
event. Nonlinear sources can require smaller timesteps and independent refinement.

All physical state and observer history are tentative until the original thermal
and electrical admission succeeds. Invalid balances, unsupported inlet data,
nonfinite powers, cancellation, work exhaustion, and unsupported endpoint material
states publish neither candidate. Original boundaries, link templates, per-side
contact coverage, geometry, masses and material IDs are retained.

## Limits and verification

This is a first-order upwind, fixed-mass, prescribed-flow control-volume model.
It does not compute pressure or momentum, infer permeability through partially
molten matter, move interfaces, change density/volume, conserve species, account
for kinetic/potential fluxes, or solve mechanical/chemical work. Stream enthalpy
uses the caller's thermodynamic convention; do not add pressure-flow work already
included in it. A moving-conductor model needs conductor-frame electrical loss
and Lorentz work; this stationary Ohmic example does not supply those effects.
Cell-uniform thermal state and mixing need a physical resolution justification.
Numerical upwind smoothing is not a new material diffusivity or viscosity.

This increment contributes to MR19/MR43 and the transport prerequisite of MR29,
not completion of those broader beads. In particular MR29's computed free surface,
velocity, strain/stress, density-change continuity and experimental qualification
remain distinct work. No unsupported physical fidelity follows from passing an
energy balance alone.

Focused native commands (use the repository's pinned toolchain):

```sh
rustc --edition=2024 --test crates/fs-thermal/src/lib.rs \
  --crate-name fs_thermal -o /tmp/frankensim-thermal-tests
/tmp/frankensim-thermal-tests --test-threads=1
cargo test -p fs-lbm --features enthalpy-coupling \
  --test joule_enthalpy_composed -- --test-threads=1
```

The existing shared-thermal workflow executes these test seams. Native kernel
validation of the first, closed-flow increment (`52457eb93aad`) reported 45 passed,
0 failed, 1 ignored, including all eight new transport regressions. The boundary
increment adds four kernel and three real electrothermal regression sources;
its native result must be taken from its own CI revision, not inferred from the
first run. Independent reference equations cover 512 closed-flow cases/10,696
steps and 512 inlet/outlet-chain cases/12,288 steps; these are supplementary,
not substitutes for native or experimental validation.
