# Material heat transfer and spatial electrical/beam feedback

`strip_enthalpy::bridge::electrical::step_heat_network_candidate` composes the
existing `PreparedHeatNetwork` with the original electrical/beam/enthalpy solve.
Thermal contact resistance and gray-diffuse radiation now change the temperature,
resistivity, relaxation, current and compatible beam deformation in that solve.
No additional thermal advance, constitutive law or nonlinear solver is introduced.

```rust,ignore
let prepared = heat_network.prepare(cx, &thermal)?;
let next = electrical::step_heat_network_candidate(
    cx, &thermal, &beam, &bindings, &circuit, &prepared, &ports,
    target_s, gravity_m_s2, &line_loads_n_m, controls, &energy_guesses_j,
)?;
let (thermal, beam, circuit, report) = next.into_parts();
```

Retain the prepared network across accepted intervals to preserve its material
claims. All its original thermal contact, bulk transport and radiation results
are retained in `report.material: HeatNetworkSample`. The report and candidate
now carry a defaulted source type; bulk-only and electrical-contact entry points
retain `MaterialNetworkSample` and their existing public signatures.

Thermal and electrical contact resistances remain independent material
properties. Each is evaluated by its original owner at the proposed endpoint.
The original heat network applies thermal resistance to fresh bulk conductances
once, and adds radiosity-derived signed cell powers once. Both electrical and
mechanical physical losses are deposited by the existing coupled residual.
Radiation-only beam cells must satisfy the mechanical material's source-context
checks even without a resistor or thermal-conductivity requirement.

Closed-enclosure radiation is internal exchange, not a heater. Its measured
roundoff-net energy is retained in `closed_radiation_net_j`; it is excluded from
`other_external_heat_j`. The combined electrical/elastic/thermal balance must
close without treating that net defect as external supply. The cell thermal
account still includes all deposited powers, including this small net, and no
redistribution is used to conceal it. Mechanical numerical loss is not heat.

Four focused regressions exercise original component replay, independent
emissivity-only and thermal-contact-only changes to actual sag/current, exact
contact-only behavior when extra transfer paths are omitted, and joint refusal
under cancellation/work exhaustion:

```sh
rch exec -- cargo test -p fs-couple --test bridge_heat_network \
  --test bridge_electrical_contacts --test bridge_electrothermal
```

Native Rust compilation and tests were not executed in the editing environment:
RCH and Cargo were absent (exit 127). Independent Python/SciPy equations passed
36 parameter cases and four load/unload intervals, with maximum combined-energy
discrepancy 3.34e-16 J. Those calculations are reference equations, not Rust
execution or experimental material validation.

Existing physical limits remain: quasistatic small-strain slender beams,
constant reference density, temperature-independent spring energy/free strain,
single solid phase, fixed electrical paths/topology/C/L and fixed thermal areas.
The radiation enclosure is closed, gray and diffuse with authored view factors;
no evolving geometry, contact detection, optical ray tracing, melting or field
force is inferred. Endpoint material/thermal rates and midpoint circuit current
are generally a first-order composition. Split intervals at scheduled events and
load/source discontinuities. Target events remain pending until joint acceptance.

## Scheduled deposition with the same complete heat network

`GrowingBridge::step_heat_network_candidate` invokes the same interval above,
then the existing scheduled birth transaction. Its defaulted report payload
preserves the original `HeatNetworkSample` through final acceptance. Bulk-only,
electrical-contact and full heat-network growth all share one completion function.

```rust,ignore
let (next, circuit_next, report) = specimen.step_heat_network_candidate(
    cx, &circuit, &prepared_heat_network, &ports, target_s, gravity_m_s2,
    &old_cell_loads_n_m, controls, &old_cell_energy_guesses_j,
)?;
```

No birth occurs inside a radiation, circuit or Newton trial. A final invalid
arrival rejects the solved interval, passive receivers' heat changes, finite
charge/flux proposal and all earlier appends. Incoming enthalpy is included once;
zero-duration arrivals perform no additional radiosity, mechanical or circuit
advance. The returned mechanical state still reports whether new material awaits
equilibrium under its enlarged load field.

The caller must supply complete bindings for the new topology before the next
interval. This method does not infer new radiative surfaces, change enclosure
view factors, migrate selected surface claims or add circuit branches. The
single-layer clamped-seed growth restrictions remain unchanged.

Three further regressions compare original interval/event/append replay, reject
a late invalid birth after a valid append, and retain passive incoming enthalpy
without debiting the circuit or adding a beam segment. Run:

```sh
rch exec -- cargo test -p fs-couple --test bridge_heat_growth \
  --test bridge_heat_network --test bridge_electrical_growth
```

These additional Rust event tests are authored but unexecuted here. Independent
component equations do not establish native event ordering or rollback behavior.
