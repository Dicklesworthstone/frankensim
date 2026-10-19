# One material exchange network

`fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::exchange`
composes the existing electrical/contact network, bulk and contact conduction,
closed-enclosure radiation, and any supplied gas chambers on one thermal trial.
It introduces no material law, time integrator, circuit solver or heat store.

```rust,ignore
let assembly = MaterialExchangeNetwork {
    electrical: &electrical_contact_network,
    thermal_contacts: &thermal_contact_links,
    radiation: Some(&material_radiation),
};
let prepared = assembly.prepare(&cx, &thermal)?;
let proposal = prepared.step_candidate(
    &cx, &thermal, target_s, &circuit, &source_ports,
    &chambers, &gas_patches, controls, &energy_guess_j,
)?;
let (next_thermal, next_circuit, report, combined) = proposal.into_parts();
// Publish the circuit, thermal state and report.source.chambers together.
// All heat and pressure work are already accounted for. Do not apply them again.
```

Empty chamber and gas-patch slices give a finite electrical/thermal system.
Empty thermal-contact bindings omit those resistances; `radiation: None` omits
radiation. Neither choice infers missing physical data or supplies a default.
To compose other sources through the existing callback APIs, reuse
`prepared.resolve(&cx, trial_thermal)` and retain its complete source payload.

## Same state, separate physical mechanisms

The existing prepared electrical network resolves bulk electrical/thermal
properties and electrical interface resistance. Thermal contacts add their
resistance in series with freshly evaluated bulk-only link coefficients.
Radiation resolves each surface at the receiving cell's actual enthalpy/phase
and supplies signed power from the original radiosity solve. The original
circuit/enthalpy owner deposits bulk and contact Joule heat once. The original
multi-chamber owner solves each gas once with all its heat contacts and pressure
work. Every trial starts from the same accepted histories.

For a link shared by electrical and thermal contact transport, both bindings
must use the same immutable ordered interface card and complete physical
conditions. The two properties still have independent source selections and
units; neither resistance is inferred from the other. Radiating surfaces must
retain their bulk cell's exact material card and complete thermal-condition
point. Original caloric-law, quantity, phase, area and source-domain checks are
retained by their existing implementations.

Preparation pins the initially selected bulk, electrical-interface,
thermal-interface and emissivity claims. Numerical values and phase states are
not cached. All properties are evaluated again at each actual trial. A prepared
binding may survive accepted steps; explicitly preparing at a new starting
state permits new policy selections. Leaving a pinned domain refuses rather
than switching to a competing claim. Each radiosity solve builds its operator
from the current surface emissivities, not an old factorization.

## Energy accounting and source payload

`MaterialExchangeSample::bulk` contains the original bulk-only link samples.
`contacts` contains the actual series-corrected conductances and contact source
receipts. `radiation` contains the actual surface emissivities, cell powers,
phase states and radiosity residuals. Gas interfaces consume the original
complete bulk snapshots, not a new conductivity-only property bundle.

The top-level `NetworkStep` excludes the net roundoff from closed-enclosure
radiation from external supply and independently checks the resulting total
energy residual. Constituent circuit/thermal reports remain unchanged: their
source maps include internal transfers exactly as handed to those solvers.
Radiative heat is never clipped, redistributed or counted as Joule heat.
Convergence, local thermal balances, circuit storage/work, gas/momentum work,
radiosity equations and combined energy checks must all pass before returning
any endpoint. Cancellation of the supplied `Cx` returns no partial endpoint.

## Physical and numerical scope

This is fixed-geometry bulk/interface transport and an opaque gray-diffuse
closed enclosure on uniform thermal cells, coupled to fixed-topology circuits
and independent translating-wall, closed frozen-composition gas chambers.
Radiation treats those gases as nonparticipating. It does not model combustion,
participating-gas radiation, arcing, pressure-driven contact-area change, moving
melt, or mechanically linked crankshafts. Contact areas, thermal path lengths,
view factors and electrical heat partitions are explicit physical inputs.
Radiation and gas heat transfer may coexist on the same exposed surface;
each mechanism's own aggregate coverage checks still apply.

The underlying thermal endpoint/gas midpoint composition remains generally
first order. Energy balance is not a timestep-error estimate or empirical
qualification. The module adds no dependency or alternate numerical algorithm.

## Focused tests

```sh
cargo test -p fs-couple --test material_exchange_network
```

Four authored regressions cover an emissivity-only change through both bodies,
the circuit and two gas motions; reduction to the original solver when extra
mechanisms are omitted; exact coefficient-owner equivalence and latent-phase
feedback; and source/area/cancellation refusal without changing accepted states.
The fixtures reuse the original radiation and thermal-contact test sources.
Native compilation and execution were unavailable during authoring. The separate
Python/SciPy reference is an equation check, not execution of these Rust tests.
