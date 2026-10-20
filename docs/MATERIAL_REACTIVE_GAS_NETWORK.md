# Reacting material gas networks and moving walls

`piston::exchange::graph::reactive::step_candidate` connects the original
shared-snapshot nozzle network, `GasReaction` and `GasPiston` in one physical
interval. `piston::material::graph::reactive_step_candidate` retains the original
source-resolved rigid bodies around that same implementation. Gas composition,
thermal state, pressure work and wall inertia are not unrelated coefficients.

## Flow, chemistry, motion: one interval

Every nozzle first samples the same accepted pressure/composition/caloric state.
The existing material transfer owner combines all inlet/outlet species masses
and upstream enthalpies before reconstructing one state per active volume.
Each selected volume then runs its existing declared reaction **once**, from
that actual mixed state, at fixed mass, volume and total internal energy.
Finally each existing piston advances the reacted gas and wall momentum.
Fixed-volume gas stores can react without becoming moving walls.

The stages act over the SAME duration, not three consecutive clock intervals.
Incoming reactants can participate now. Newly produced species leave in the
next flow interval. The reaction is not repeated for mechanical Newton probes.
Final pressure changes affect the next nozzle evaluation. This is an explicit
first-order operator split, not a simultaneous or stiff reactive-flow solver.
Refine physical time and split at changes to externally prescribed inputs.

Source formation energy is already included in gas internal energy. Do not add
reaction enthalpy as extra heat; do not deposit transported enthalpy or jet energy
a second time. The whole-network balance is measured from original to FINAL gas
and wall states:

`delta U_gases + delta K_walls + W_exterior - W_applied - Q_external = 0`.

The flow report retains its transport-only reference, each reaction report its
post-flow chemical reference, and each piston report its post-chemical reference.
The enclosing `MovingGasStep` gas storage and total energy cover all three.
Individual species are conserved during transport and changed only by the
selected reaction's original certified stoichiometry. Each reaction keeps its
existing source-consistency, consumed-fraction, species, energy and entropy gates.

## Data bindings

Reuse `ChamberReaction` from the existing finite thermal-network chemistry API;
it is re-exported by the new gas graph. Its `chamber` is the existing graph id,
including a finite fixed-volume gas when appropriate. There can be one declared
reaction per selected volume. Unlisted volumes remain nonreacting; duplicate or
foreign ids refuse. A catalog must include every incoming species and initially
absent product with consistent caloric, transport and pack sources.

```rust,ignore
use fs_couple::cavity_phs::piston::material::graph::{
    reactive_step_candidate, ChamberReaction,
};

// `walls` already owns MaterialGasPiston values resolved from solid cards.
// `reservoirs`, `links`, `endpoints`, `drives` and `limits` are the existing
// moving-network data. `reaction` and `kinetic_limits` are existing owners.
let chemistry = [ChamberReaction {
    chamber: 1,
    reaction: &reaction,
    controls: kinetic_limits,
}];
let proposal = reactive_step_candidate(
    &walls, &reservoirs, &links, dt_s, &endpoints, &drives,
    &chemistry, limits, || cancellation_requested(),
)?;
let (next_walls, next_reservoirs, interval, reactions) = proposal.into_parts();
// Publish every returned state together. The reported extents are already applied.
```

Changing supported reaction/material data changes the shared gas state and
subsequent flow and motion without modifying any solver. Changing a wall's
material density before construction changes the original body's canonical
inertia and the resulting gas trajectory. Solid source records survive the step
exactly. Gas chemistry does not silently heat, expand or melt the rigid wall.
No live mass/inertia exchange is inferred from a material substitution.

## Atomic publication and port observations

All flow, chemical and mechanical stages plus the independent final balance must
succeed before any state is returned. A late reaction refusal, a last-wall failure
or final cancellation leaves every input and port window untouched. Empty
reaction bindings use the identical existing nonreacting physical path.

Closed reacting walls have one accepted wall/heat port window including actual
chemical entropy change. Gross incoming OR outgoing mass still makes that
wall/heat-only audit unavailable: stream ports are needed to certify an open
chamber. Its real wall/heat observations and the enclosing graph balance remain.
An energy balance alone does not prove passivity or physical time accuracy.

## Focused checks and limits

```sh
cargo test -p fs-couple --test reactive_moving_gas
cargo test -p fs-couple --lib cavity_phs::piston::exchange::graph::
cargo test -p fs-couple --lib cavity_phs::piston::material::graph::
```

Ten new integration tests use synthetic balanced isomers and existing material
loaders. They cover original-owner parity, incoming-reactant activation, product
transport timing, symmetry/order, fixed reactors, source shutdown, density/rate
substitutions, refinement and late failures. They are authored but native
execution was unavailable: the authoring environment has no cargo/rustc.
Independent constant-cp Python equations over 480 steps show first-order velocity
refinement against DOP853 and whole-history energy residual below 4.30e-13 J.
These are reference equations, not Rust execution or measured-gas validation.

The physical scope remains uniform ideal gases, one declared reaction per
selected volume, quasi-steady restrictions and independent constant-area rigid
walls. This is not a complete combustion mechanism, ignition model, flame or
shock solver, implicit flow network, shared crankshaft, deforming solid or live
acoustic renderer. Updated gas sound speed, impedance and transport come from
the original material state and remain available to their existing consumers.
