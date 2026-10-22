# Shared-state, material-backed reaction mechanisms

`GasReaction::combine` connects multiple declared homogeneous reaction channels
to the existing gas, piston, transport, electrical and thermal assembly paths.
It returns the existing `GasReaction` type: existing `ChamberReaction` bindings
continue to contain ONE chemistry owner per chamber, now with multiple channels.
No new species registry, gas state, equation of state or material-card schema is
introduced. The original single-channel constructor remains available.

## One inventory, one temperature reconstruction

Every channel sees the same accepted species amounts, temperature, pressure and
material conditions, after transport when the caller uses a transport adapter.
Its original mass-action owner supplies a signed rate `r_j` in mol/(m3 s), with
dimensionless activities `p_i/p0`. The combined explicit update is

```
dxi_j = duration * volume * r_j
n_i(new) = n_i(old) + sum_j nu_ij * dxi_j
U(new) = U(old)
```

All contributions are combined before the existing gas owner reconstructs the
new composition, temperature, pressure, transport and acoustic properties ONCE.
Formation energy is included in U; no reaction heat is deposited a second time.
Products created by one reaction cannot immediately feed another during this
explicit substep. The following physical interval samples the evolved state.

The consumption limit applies to the SUM of negative species contributions
across signed NET channel extents, divided by the ORIGINAL species inventory.
It is not a separate limit that each competing reaction can spend independently.
Production by another channel cannot conceal excessive withdrawal. This is not
a cap on the sum of microscopic forward and reverse equilibrium turnover.
Actual species amounts, mass, energy and total endpoint entropy retain the
original acceptance checks. No clipping or invented trace concentrations occur.

## Source data and construction

All channels must use one identical complete `ThermochemicalGasModel` catalog
(including its declared nominal composition and original transport/pack bindings),
select distinct reaction IDs from the SAME original `StoichiometricMatrix`, and
retain the SAME elemental/charge `ConservationCertificate`. This avoids combining
incompatible element definitions merely because species names and masses match.
No implicit catalog merge or inference of chemical accessibility is performed.
Each selected column already passes the original molar-mass and caloric-source
checks in `GasReaction::new`.

Bind each optional `MaterialRateScale` before combining:

```rust,ignore
let first = GasReaction::new(catalog.clone(), first_column_law)?
    .with_material_rate(first_rate_source)?;
let second = GasReaction::new(catalog.clone(), second_column_law)?
    .with_material_rate(second_rate_source)?;
let chemistry = GasReaction::combine(vec![second, first])?;
let proposal = chemistry.parcel_candidate(&gas, duration, controls, cancelled)?;
```

Each rate card resolves at the actual shared gas point, retaining its complete
sample. The original reverse-equilibrium closure uses that channel's SAME
forward scale. A fixed scale and a sourced scale may coexist. Missing source
data or an unsupported temperature refuses the complete proposal; no fallback
is selected. Changing a rate table or thermochemical source is data-only.
Binding one ambiguous rate to an already combined mechanism refuses.

Combination is bounded to 32 total channels, including flattened combinations.
Canonical reaction-ID ordering removes caller vector-order dependence; it does
not make evaluation order into sequential physical chemistry. Inputs are never
mutated, and cancellation is polled between channels and species operations.
A singleton combination returns its original object unchanged.

## Reports and existing consumers

Legacy `kinetics()`, `rate_material()`, and report fields `rate`, `rate_material`
and `extent_mol` refer to the first canonical channel. `primary_reaction` names
it. `additional_reactions` retains every other original rate, material sample
and extent. `reaction_extents()` visits ALL channels. Different reaction extents
must not be added as if they were one stoichiometric reaction.

The enclosing species, consumption, energy and entropy fields cover the ENTIRE
mechanism. `has_progress()` also notices an active secondary reaction when the
primary is idle. The original moving-gas graph uses it to retain the whole
chemical/pressure-work entropy and storage account, with one physical clock.
Open gas streams still require the enclosing network audit rather than a
closed two-port-only passivity claim.

Existing `GasPiston::reactive_candidate`, material-inertia wrappers, reacting
flow graphs, and the live electrical/thermal assembly consume the same parcel
entry point. Their nonlinear mechanical/thermal trials never repeat chemistry.
There is no additional kinetics integration path in those consumers.

## Checks and limits

```sh
cargo test -p fs-couple --test reactive_mechanism
cargo test -p fs-couple --lib cavity_phs::piston::reactive::
```

Seven new regression sources cover singleton/order equivalence, competing and
chained reactions, joint depletion despite another channel's production,
nonisomeric molecule-count changes, per-channel live material scales with the
original reverse closure, source/certificate refusal, cancellation and the
existing piston/network path with an idle primary channel.

Native compilation and Rust execution were unavailable in the authoring
environment (`cargo` and `rustc` absent). Independent manufactured Python
reference equations checked 512 multi-channel updates: maximum energy roundoff
3.41e-13 J and amount discrepancy 1.39e-17 mol. A 560-step trajectory study
against DOP853 gave error-reduction ratios 2.0019 and 2.0010 and maximum energy
drift 4.61e-12 J. These do NOT execute the Rust or validate measured chemistry.

This remains explicit first-order homogeneous ideal-gas evolution. It is not a
stiff-kinetics solver, ignition/flame model, complete validated combustion
mechanism or proof of arbitrary-step stability. Source discontinuities require
caller-selected interval boundaries. Small admitted steps and refinement remain
necessary. The existing formulation/reference convention is retained:
https://cantera.org/stable/reference/reactors/ideal-gas-reactor.html
