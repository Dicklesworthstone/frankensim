# Species flow through the live material assembly

`MaterialExchangeState::advance_with_flow` connects the existing material
exchange state to a graph of finite gas transfers. No caller manually replaces
a chamber's species inventory, copies incoming enthalpy into heat, or advances
one connected state without the others.

```rust,ignore
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::{
    GasNetworkLink, GasNetworkSolve,
};

let report = state.advance_with_flow(
    &prepared_materials, &cx, target_s, &electrical_ports,
    &chamber_drives, &reactions, &gas_wall_patches,
    &nozzle_links, flow_controls, coupled_controls, &energy_guess_j,
)?;
// Circuit, thermal cells, all species inventories and all moving walls now
// belong to the same accepted instant. Reported transfers are already applied.
```

Each `GasNetworkLink` gives a stable ID, two existing chamber IDs, the original
`GasNozzle` geometry and `GasNozzleStepSolve`. Either endpoint can supply gas.
Shared nodes, parallel restrictions and cyclic graph topology are supported.
The lower-level `nozzle::network::step_candidate` is also usable directly with
`BTreeMap<u64, ClosedGasParcel>`, without pistons or electrical/thermal states.

## Physical ownership and time levels

The graph reuses the original source-resolved nozzle and two-parcel transfer
implementations. Choking, caloric inversion, species mixing, enthalpy transport,
source conflicts and pressure-crossing refusal remain there. Neither gamma nor
an air composition is introduced by the graph. Incoming species retain their
original caloric, transport and pack records. Both sides are finite reservoirs;
there is no implicit infinite ambient supply. The live assembly currently gives
each gas owner its existing translating-wall dynamics.

Edges are explicitly **Lie-split in increasing stable-ID order**, each over the
same physical duration on the result of preceding edges. Vector order is not
execution order. Changing stable IDs can change finite-step splitting error.
This is not a simultaneous pressure solution. Each original pair has its donor
fraction bound; the graph additionally limits gross withdrawal from each node
relative to its INITIAL mass, so incoming flow cannot hide excessive outflow.
Pairwise equilibrium-crossing admission is not an accuracy estimate for the
whole graph. Refine the interval, and split at source or valve discontinuities.

The assembly stages that graph once, then its original reaction stage once per
selected chamber, then the original coupled thermal/electrical/pressure-work
solve. Every stage approximates its operator over ONE common duration, derived
from the thermal clock. A nonlinear trial never transfers or reacts another
portion of accepted inventory. Arriving species must be supported by any active
reaction catalog; unsupported sources refuse instead of disappearing.

This remains a generally first-order composition of homogeneous finite volumes,
quasi-steady restrictions, declared reactions and thermal/motion models. It does
not add resolved pipe storage/waves, jet thrust, a valve actuator, check valves,
flames, vacuum endpoints or a stiff simultaneous kinetic-flow solve. Existing
radiation treats the gas as nonparticipating. New acoustic material states do
not by themselves rebuild an already-running renderer.

## Acceptance and observations

`flow_step_candidate` previews the complete advance. `advance_with_flow` moves
all resulting states into the live assembly only after final cancellation and
whole-interval conservation checks. The existing state assembly is reused with
staged chamber inputs, avoiding an extra copy of the accepted thermal timeline
and circuit merely to change gas inventories.

The graph checks total mass, every species and extensive internal energy from
its actual endpoints. With reactions enabled, individual molecular species may
then change under the reaction owner's conservation certificate. The assembly
checks total gas mass and an independent whole-interval energy balance including
electrical storage, thermal storage, gas energy, wall kinetic energy and external
work/heat. Flow enthalpy, internal wall heat, Joule heat and formation energy are
never counted twice. Extensive gas differences use the actual initial/final
masses, not the closed-gas formula m*delta_u.

`MaterialFlowAdvance.flow` retains each original throat and pair-transfer report;
`exchange` retains post-flow reference states in the existing reaction/thermal
and mechanical reports; `chambers` records full-interval extensive changes.
The top-level mass and energy residuals enclose the original live state through
the final endpoint. No report is an additional energy or species source.

A chamber with any nonzero stream transfer keeps one whole-interval port tick.
Its legacy heat/work ports cannot represent species streams, so their closed
boundary audit is explicitly unavailable, even at zero NET mass change. Actual
heat and mechanical transfers remain observable; use the enclosing flow balance.
No-flow chambers retain the original closed-step port records unchanged.

## Focused checks

```sh
cargo test -p fs-material --lib thermochemical_gas::energy::parcel::transfer::nozzle::network
cargo test -p fs-couple --test material_exchange_network
```

Eight authored Rust regressions cover original-owner reductions, branching and
transit, input permutation, reverse flow, aggregate withdrawal, original source
preservation, incoming-material changes through the complete live assembly,
reaction extent evaluated once after transfer, and failure/cancellation/retry.
Fixtures are synthetic. Native Cargo was unavailable during authoring; source
review and separate Python equation checks are not Rust compilation/execution.
