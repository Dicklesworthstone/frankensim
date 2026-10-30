# Electrical heat and relaxation on the same material trajectory

`FrameRelaxationThermalSystem::bind_material_circuit` connects the original
source-bound circuit to the same finite layer inventories already driving
multi-body relaxation and fluid-loaded motion. A charged capacitor, inductor,
or held electrical source can heat those layers. Their actual temperatures
change the selected relaxation times and subsequent reciprocal mechanical
response. A second temperature, heat capacity or mechanical integrator is not
created. Heating a fully relaxed body at rest does not generate motion in this
rate-only model: there is no thermal-expansion or field-force actuator here.

## Bind the original source network

Start with the existing single- or multi-frame relaxation thermal owner. Add
`CellResistor` bindings to a `MaterialNetwork` that borrows the SAME `cells` and
`links` slices, not copies. Use its explicit conductivity/resistivity choices to
prepare the initial circuit. Numeric resistor placeholders in the graph are
replaced by the original conductor law. Branch path length times area must
match its receiving cell's volume, under the original whole-cell model.

```rust,ignore
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::
    electrical::RelaxationElectricalInput;

let prepared = electrical_network.prepare_with_bases(&thermal, &bases, || false)?;
let initial = prepared.resolve(&thermal, || false)?;
let circuit = MaterialCircuit::new(
    graph, initial.resistors, initial_fluxes, initial_charges, &initial_ports,
)?;
owner.bind_material_circuit(
    &thermal, &electrical_network, circuit, &bases, &gate,
)?;
let report = owner.step_electrical(
    &mut thermal, &held_native_mechanical_forces,
    RelaxationElectricalInput { ports: &held_electrical_ports, solve: circuit_controls },
    thermal_controls, &cell_energy_guesses_j, &gate,
)?;
```

The input is an alias of the existing `FrameElectricalInput`; its electrical
work and power tolerances are independent of the thermal/Newton controls.
An empty port slice means finite stored charge/flux powers the circuit.
`electrical()` exposes read-only accepted charge, flux and source snapshots.
Ordinary `step` refuses an electrically bound owner, even with zero source
ports, because zero drive is not zero stored energy.

Cold binding does not change any original state or clock. Every resistor must
match the initial resolved material and geometry. Rebinding after stepping,
a different cell/link owner, incomplete bases or cancelled preparation refuses
without altering the owner. Electrical properties may be present only as
conductivity: no resistivity claim is fabricated.

## One proposed endpoint and separate energy mechanisms

The mechanical propagator first stages the interval at accepted relaxation
rates. Independently integrated Maxwell heat stays assigned to its original
layers. Every thermal residual evaluates the original material network at the
proposed temperatures, then runs the EXISTING native circuit from its SAME
accepted charge/flux under held source ports. Discarded Newton trials cannot
discharge the physical capacitor. Circuit solves count in the shared trial
budget, including finite-difference stencils and final fresh admission.

The shared caloric residual receives actual per-branch Joule heat and actual
per-layer Maxwell heat once. It checks combined mechanical/electrical storage,
thermal enthalpy, supplied work, external thermal exchange and exterior losses.
Only after finite-cell, source, rate-preparation and cancellation gates pass do
charge, flux, mechanical motion, fluid momenta, all material histories, layer
temperatures and the sample clock publish together.

`report.layers` and `deposited_material_heat_j` remain RELAXATION heat only.
`report.electrical` retains the native circuit account, including each branch's
actual midpoint current, voltage, source state and heat.
`deposited_joule_heat_j` is separate. The existing finite-cell entropy account's
`material_loss_generation_j_k` includes BOTH dissipative sources divided by
actual receiving endpoint Kelvin. Its temporal/algebraic defects are numerical,
not additional heat. Other mechanical/fluid losses remain external sinks.

R(T) and thermal conduction are implicit in the endpoint caloric residual.
Maxwell rates use the accepted start over the mechanical interval; updated
rates act NEXT interval. This retains the first-order thermal split. Refine time
separately from structural/fluid spectral resolution and thermal lumping.

## Explicit limits

Reference density, geometry, equilibrium and incremental moduli remain fixed.
Selected relaxation layers must retain wholly solid caloric charts and their
original source-domain validity. The bulk-only circuit binding holds authored
C/L fixed; material-dependent storage requires the field binding described below.
The bulk-only entry point refuses electrical contacts without their two-sided
source/heat binding; the interface path below provides that binding. There
is no moving wire geometry, motional EMF, magnetic force, nonlinear E/J response,
thermal expansion, melting, spatial temperature gradient or acoustic-heating
claim. Native circuit, source heat, thermal and preparation budgets all apply.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature::heat
```

Five electrical regression groups are authored beside the existing single/shared
heat tests. They cover analytical RC/RLC midpoint charge/flux, current-source
work, per-layer heat, conductivity-only sources, data-only resistivity changes,
complete energy, and late refusal/retry without charge or history loss.
Native Rust compilation, execution and formatting remain UNEXECUTED in the
editing environment: Rust/Cargo/rustfmt/RCH/DSR are absent and the toolchain
host cannot be resolved. Source inspection is not native qualification.

Independent NumPy/SciPy equations ran seven RC/RLC wet/dry cases at four time
resolutions against continuous DOP853. Maximum combined-energy discrepancy was
4.03e-18 J; time refinement reduced temperature and trajectory errors. Changing
only core resistivity changed the synthetic continuing displacement by 2.46e-8 m.
The fully relaxed stationary control heated electrically without moving. These
are manufactured reference equations, not execution of Rust, material measurements,
or proof of full-workspace or real-time performance.

## Material-dependent electric and magnetic storage

`bind_material_field_circuit` and `bind_material_contact_field_circuit` bind
the existing `FieldThermalBinding` objects to this same relaxation/heat owner.
Every capacitor and inductor must retain its original material law, physical
volume, finite caloric cell and charge/flux envelope, including constant laws.
Use `step_electrical` or `step_heat_network` normally after binding. The bulk-
only binding without fields still refuses unaccounted material field storage.

```rust,ignore
owner.bind_material_field_circuit(
    &thermal, &electrical_network, source_bound_circuit, &bases,
    field_thermal_bindings, &gate,
)?;
// Each FieldThermalBinding identifies branch, cell, original source law and
// maximum_abs_coordinate. There is no independent heat-capacity override.
```

The existing field owner evaluates C(T)/L(T) and its thermodynamic potential.
It retains charge or flux during a proposed coefficient replacement, then the
existing native circuit advances from those SAME accepted coordinates. The
field's internal energy is U = Psi + T*S, not just the isothermal circuit
Hamiltonian. Its signed reversible heat reaches the original field cell next
to, but separately from, Joule and material-relaxation loss. A negative source
slope can cool that cell while electrical resistance dissipates positive heat.
No energy residual is used as a physical heat source.

`report.field_storage` retains the original field source samples, actual
endpoint coordinates, internal energies, parameter work and reversible heat.
`report.field_entropy` records entropy storage, signed reversible Q/T and the
first-order splitting defect. Include its reversible transfer term along with
radiation and the existing caloric entropy terms; its signed temporal defect
is not irreversible entropy production. `external_thermal_heat_j` excludes
both field exchange and closed radiation. The combined energy gate uses
actual field U and original mechanical/caloric storage at both endpoints.

Field admission is prepared once per physical interval, not per nonlinear
stencil. Several independent fields may occupy one whole material cell: their
capacity contributions add while its caloric capacity is counted once. Missing
fields, foreign receivers, invalid source domains, exceeded coordinate limits
and nonpositive combined capacity refuse before publication. Late solver or
cancellation failure leaves all electrical, mechanical, material, fluid and
thermal histories available for retry.

Geometry, density and mechanical storage moduli remain fixed. Source field
laws obey the existing globally affine, solid-state model. This is not a
Maxwell-force, motional-EMF, hysteresis, thermal-expansion or phase-change model;
the native first-order thermal split and its refinement requirement remain.

Three native regression groups cover independently calculated RC/RL motion,
U/S and reversible heat, source-only coefficient changes, whole-cell/capacity
admission, and late rollback/retry. They are authored but UNEXECUTED: this
executor has no Rust/Cargo/rustfmt/RCH/DSR, and toolchain-host DNS is unavailable.
Independent NumPy/SciPy calculations passed nine wet/dry, capacitive/inductive
and complete-heat cases at four time resolutions, with maximum combined-energy
discrepancy 4.07e-18 J. Refinement reduced error against continuous DOP853.
These synthetic equations are not native execution, experimental validation or
a performance measurement.


## Original electrical interfaces and two-sided heat

`bind_material_contact_circuit` takes the existing
`circuit_enthalpy::contacts::material::ElectricalContactNetwork`. Its bulk source
must borrow the same original cell/link slices. Each interface retains its
ordered material pair, selected contact-resistance claim, fixed thermal-link
area, and explicit A/B heat partition. No bulk property substitutes for an
interface law, and no thermal partition is inferred from resistance.

```rust,ignore
let resolved = contact_network.prepare(&thermal, || false)?
    .resolve(&thermal, || false)?;
let circuit = MaterialCircuit::new_with_contacts(
    graph, resolved.bulk.resistors, resolved.contacts,
    initial_fluxes, initial_charges, &initial_ports,
)?;
owner.bind_material_contact_circuit(&thermal, &contact_network, circuit, &gate)?;
// Continue through the SAME step_electrical call shown above.
let routes = owner.electrical_contact_heat_destinations();
```

The original contact evaluator sees BOTH proposed side temperatures at every
caloric trial, preserving pinned sources and applicability. Native circuit bulk
and interface losses stay separate in `report.electrical`; their sum contributes
to `deposited_joule_heat_j`. The original contact-deposition helper delivers
interface loss to the two declared cells exactly once, beside their own bulk
Joule and Maxwell losses. A zero or unit partition is retained exactly.
Thermal conductance remains an independent source law, not inferred from the
electrical contact or mechanical connection.

The initial circuit must agree with original interface area, selected sources,
physical queries and resistance. A changed selection-policy receipt from source
pinning is not confused with a changed material. Source-domain, late nonlinear,
budget and cancellation refusals preserve circuit charge/flux/contact states,
all material/fluid histories, temperatures and clock. No accepted endpoint is
published before both the circuit and the finite thermal solve are admitted.

Interface area, context and heat fractions remain prescribed. This extension
adds no contact opening, pressure evolution, arcing, contact capacitance or
spreading-temperature model. The existing bulk/source and rate-only mechanical
limits still apply; electrical heating does not invent an expansion force.

Three additional Rust regression groups cover independent RC/contact work,
separate two-cell caloric balances at fractions 0, 0.35 and 1, interface-only
source substitution changing motion, initial area mismatch, and late interface-
domain/cancellation refusal with retry. These native tests remain UNEXECUTED.
Independent six-case wet/dry reference equations at four time resolutions gave
a maximum combined-energy discrepancy of 3.81e-18 J. Changing only interface
resistance changed continuing displacement by approximately -1.1723e-8 m in the
synthetic wet case. Refinement reduced trajectory and temperature errors against
continuous DOP853. These equation checks are not native execution or measured-
material validation.

## Complete heat transport on the same relaxation trajectory

`FrameRelaxationThermalSystem::step_heat_network` accepts the existing
`PreparedHeatNetwork`. Thermal interface resistance, source-integrated k(T)
layer conduction and gray-diffuse enclosure radiation now act on the same
finite cells as Maxwell loss and bulk/interface Joule heat. The same entry
supports unpowered material assemblies. It replaces bulk-only trial transport;
it does not add a parallel temperature state, contact law or heat integrator.

```rust,ignore
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::
    transport::RelaxationHeatInput;

let heat = material_heat_network.prepare(&cx, &thermal)?;
// Optionally use heat.with_integrated_conduction(...) with the original layer
// source preparation. A link cannot also have a separate thermal-contact owner.
let report = owner.step_heat_network(
    &mut thermal, &held_native_mechanical_forces,
    RelaxationHeatInput { network: &heat, cx: &cx },
    Some(electrical_input), // None only for an electrically unbound owner.
    thermal_controls, &cell_energy_guesses_j, &gate,
)?;
```

Keep the same preparation between intervals to retain its original source pins.
The heat assembly must borrow the original material-cell and thermal-link
slices. Before speculative motion, its resistor sources, accepted temperature
queries, conductor geometry, electrical interfaces, receiving cells and A/B
heat fractions are checked against the continuing circuit. An explicit new
heat preparation is a caller intervention, never an automatic source fallback.

Each proposed thermal endpoint resolves the complete network once. The native
circuit still starts every trial from its accepted charge/flux. Independently
integrated Maxwell loss, actual electrical loss, contact conduction and signed
radiative powers enter the original finite-cell residual exactly once. Source
rate updates and all mechanical, fluid, electrical and caloric histories still
publish at the existing common acceptance boundary. Cancellation through either
Cx or the mechanical gate refuses before publication; a Cx cancellation also
requests the native gate. Failed numerical/source trials cannot drain storage.

`report.source` retains the original bulk snapshots. `report.heat_transfer`
retains accepted thermal-contact, source-integrated layer and radiation results.
`closed_radiation_net_j` records the enclosure's signed roundoff net without
correcting individual cell powers. The physical total-energy gate excludes that
internal exchange from `external_thermal_heat_j`, under the original tolerance.
`radiation_exchange_entropy_j_k` is the signed material sum of Q/T at the actual
receivers; include it with the terms in `report.entropy` to reconstruct caloric
entropy storage. It is neither photon entropy nor a reversible-source claim.
Authored heat stays separately classified, even in a radiating assembly.

Geometry, view factors and contact areas remain fixed. The original solid,
rate-only storage model and first-order thermal split still apply. This does
not infer moving-enclosure radiation, contact opening, melting, expansion,
field forces or fluid heating from exported exterior losses.

Four additional native regression groups exercise continued powered/unpowered
motion, independent RC and two-cell transport equations, emissivity-only and
thermal-interface-only substitutions, integrated affine k(T), authored-source
classification, empty-path parity, electrical-interface/partition continuity,
and source-domain/budget/cancellation rollback with retry. They are authored but
UNEXECUTED here: Cargo, Rust, rustfmt, RCH and DSR are unavailable. Run the focused
heat test command above in the native workspace before treating this as qualified.

An independent NumPy/SciPy reference ran six wet/dry, powered/unpowered and
integrated-conduction cases at four time resolutions. Its maximum combined
energy discrepancy was 4.21e-18 J; refinement reduced temperature and trajectory
errors against continuous DOP853. Changing only emissivity from 0.2 to 0.8
changed the synthetic continuing displacement by approximately -8.24e-9 m.
These are checks of manufactured equations, NOT execution of Rust, experimental
material validation, a full-workspace pass or a performance measurement.
