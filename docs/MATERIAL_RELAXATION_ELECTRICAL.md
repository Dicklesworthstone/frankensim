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
original source-domain validity. C/L are explicitly authored fixed storage;
material-bound field storage refuses without its reciprocal thermal-work owner.
Electrical contacts refuse without their two-sided source/heat binding. There
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
