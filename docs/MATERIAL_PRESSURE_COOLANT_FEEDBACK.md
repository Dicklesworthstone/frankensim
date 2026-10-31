# Pressure-driven circulation inside the original material/thermal transaction

`SourceConvectionHeatNetwork::step_pressure_driven_relaxation_with_links`
connects the existing series pressure law to the existing finite coolant,
convection, upwind enthalpy, circuit/field and material-memory transaction.
Every proposed thermal endpoint resolves the original bulk coolant materials,
solves the series flow, derives local velocities from original hydraulic areas,
evaluates declared bulk/film convection, and carries that SAME flow through the
original implicit caloric residual. No pressure-root search, new material law,
second thermal state or additional time integrator is introduced.

## Original inputs and one coherent trial

The new types live in
`fs_airflow::graph::thermal::medium::convection::deposition::source::heat::relaxation::pressure`.

`PressureCoolantSection` pairs an original finite coolant `CellId` with its
`MaterialFlowResistance`. Section order defines positive closed-loop flow.
`PressureConvectiveLink` names an existing wall/coolant thermal patch and a
hydraulic section, with explicit correlation, characteristic length and
property-temperature convention. It has no independent speed or flow-area knob.
The existing thermal patch remains authoritative for area and surface coverage.

```rust,ignore
let path = PressureCoolantLoop {
    sections: &hydraulic_sections,
    links: &wall_links,
    driving_pressure_pa: held_pressure,
    pressure_budget,
    maximum_links,
};
let report = convection.step_pressure_driven_relaxation_with_links(
    &cx, &mut owner, &mut thermal, &held_native_mechanical_forces,
    electrical_input, thermal_controls, &energy_guesses_j, &gate, path,
    |section, trial, requested_temperature| {
        original_coolant_laws[&section.cell].at_temperature(requested_temperature)
    },
)?;
let accepted = report.transport_source.as_ref().unwrap();
let mass_flow = accepted.hydraulics.mass_flow_kg_s();
```

Initialize the original thermal owner and prepared heat network with the
existing `*_with_link_conductances` APIs, as for prescribed finite convection.
Those initial link coefficients admit coverage; they are not frozen subsequent
conductances. The flow and heat source callback is repeatable and must retain
original material pins, caloric identity, pressure/composition and support.
A generic transport trait cannot prove that a caller selected the right source
card. The pressure and film results retain the actual typed source objects.

The shared lower-level
`FrameRelaxationThermalSystem::step_heat_network_with_resolved_advection`
accepts `RelaxationFlowCoefficients { heat, flow }` from one endpoint evaluation.
Existing prescribed-flow callers borrow their graph instead; no per-trial
clone of a prescribed graph is introduced. Every newly resolved graph is
checked for forbidden retained-solid inventories, caloric identity, enthalpy
storage, balanced fixed mass, native transport budgets and cancellation.
Changed graph rates are never accepted independently from temperatures.

The accepted payload, transported enthalpy and original source snapshots are
those of the actual final thermal trial. There is no callback after publication,
last-good-flow fallback, coefficient lag, or material-history reset. Late
source, graph, circuit, numerical-budget and cancellation failures publish none
of the mechanical state, material memory, fluid momenta, charge/flux, enthalpy,
source conditions or sample clock.

## Pressure work is not free cooling

All local `dp_i * Q_i` is deposited once into the corresponding finite coolant
cell. The prescribed pressure source supplies that hydraulic work; it is already
included in `external_thermal_heat_j`. The per-cell powers remain visible in
`transport_source.hydraulic_heat_w`. Do not add them again as wall heat or Joule
loss. Closed circulation carries enthalpy internally and has zero external
stream supply. Bulk/contact Joule heat, solid relaxation loss, field-reversible
heat and exported mechanical-fluid losses retain their original meanings.

The generic entropy account retains the pressure-source heat as an explicitly
unclassified at-cell contribution. This does not claim reservoir entropy,
stream entropy, a finite pump's dissipation partition, or photon entropy.

## Physical and computational boundaries

This is a quasi-steady, closed, fixed-volume/fixed-mass coolant loop. It requires
at least two distinct cells and one common hydraulic density matching each
cell's mass/volume within roundoff. Variable-density expansion, compression,
phase change, moving geometry, hydraulic inertial storage, pump curves and shaft
or motor efficiency are not inferred. Differential drive does not alter the
material's prescribed thermodynamic pressure. A scalar pressure solve is not a
full fluid momentum field. A source/correlation refusal, including unsupported
zero speed, is not regularized into another flow or heat law.

Bulk viscosity and density determine pressure loss; declared film properties
may differ for convection, but velocity is the actual bulk volume flow divided
by original hydraulic area. Both evaluations use the same retained source
resolver, not independent property presets. Mechanical exterior-fluid inertia
and drag are still separate, explicit bindings: this coolant update does not
silently retune their frozen coefficients.

The original first-order thermal/material split and upwind spatial error remain.
Mass-flow and heat coefficients are refreshed at every original Newton trial.
Hydraulic work/storage scale linearly with section count; graph construction,
material queries and native nonlinear/circuit preparation still cost work.
No measured hardware speedup, allocation-free stepping or real-time claim is made.

## Focused checks

```sh
rch exec -- cargo test -p fs-airflow --test relaxation_convection finite::advection::resolved
rch exec -- cargo test -p fs-airflow --test relaxation_convection finite::pressure
```

Five new native regression groups cover prescribed/dynamic parity, actual final
source payloads, changing rates within thermal trials, local/global pressure,
convection and RC accounts, hydraulic heat, data-only viscosity changes,
reverse drive, unpowered histories and complete late-failure/retry behavior.
These Rust tests, compilation and rustfmt remain UNEXECUTED: Rust/Cargo/rustfmt,
RCH and DSR are absent from the editing executor and network resolution failed.

Independent NumPy/SciPy manufactured equations ran six scenarios at 256, 512,
1024 and 2048 Hz against continuous DOP853. All assertions passed. Maximum
combined-energy discrepancy was 9.06e-11 J and maximum individual-cell discrepancy
was 1.10e-11 J. Refinement reduced temperature and trajectory error. Doubling
only viscosity changed the 256 Hz endpoint flow from 0.193663 to 0.100185 kg/s,
the heated coolant from 303.938221 to 305.974606 K, and capacitor charge from
7.49168e-5 to 7.48249e-5 C. These are checks of reference equations, NOT execution
of the Rust code, measured-material validation, or a full-workspace result.
