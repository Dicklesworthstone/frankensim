# Spatial piezoelectric plate, reactive shunts and acoustic sensing

Run the finite-plate example with the repository's pinned toolchain:

```bash
cargo run -p fs-couple --example piezoelectric_plate
cargo run -p fs-couple --example piezoelectric_plate -- --harmonic
cargo run -p fs-couple --example piezoelectric_plate -- --acoustic
```

The first command computes static actuation, the matching isolated-electrode
state, a voltage ramp, open-circuit ringing and resistively loaded motion.
The harmonic mode reports a small frequency sweep with open, resistive and
resistor–inductor loads. Voltage, current and deformation come from the coupled
field solution. The table also reports mean conductor heat and the residual
of the independently reconstructed power balance.
The acoustic mode advances the same plate and retained R/L/C histories under
an incident guide wave, including the fluid's reaction on the moving plate.

## Material and geometry

The example retains a 10 × 10 × 2 mm plate with P1 tetrahedra and complete
electrodes on its thickness faces. Six displacement constraints remove rigid
motion while allowing the free transverse and longitudinal strains. These
constraints are not a physically clamped face.

The original nominal PZT-5A coefficients come from the manufacturer's
[Materials Technical Data, typical values](https://info.piezo.com/hubfs/Data-Sheets/piezo-material-properties-data-sheet-20201112.pdf).
The example explicitly assumes a homogeneous, transversely isotropic poled
ceramic and fills its symmetry-related entries accordingly. It converts the
engineering shear convention to Mandel coordinates before the native
strain-charge ingress derives the stiffness, stress-charge coupling and
fully clamped permittivity. All selected source receipts remain attached.

This is an authored frozen small-signal research model. The source's nominal
values and listed measurement frequencies do not establish a measured
broadband tensor, uncertainty distribution or device-accuracy claim. The
example's frequency sweep explores this approximation; it is not an
experimental validation or a certified resonance prediction. Intrinsic
mechanical/dielectric losses are omitted. The external load and ideal winding
are explicitly declared models, not inferred properties of PZT-5A.

## Retained reactive loads in time

`piezoelectric::spatial::transient::load::reactive` adds a joint owner for the
original spatial body and external electrical storage:

| Load | Original source owner | Retained coordinate | Physical law |
| --- | --- | --- | --- |
| Series R–L | `ResolvedInductor`, optionally `ResolvedConductor` | Flux linkage λ [Wb] | λ̇ = V − RI, I = λ/L |
| Parallel capacitor | `ResolvedCapacitor` | Terminal charge q [C] | q̇ = I, V = q/C |

The winding's L is derived from the original permeability, magnetic path,
cross-section and turn count. The capacitor's C is derived from its original
permittivity, facing area and separation. A capacitor outside the ceramic is
additional physical storage; it does not duplicate the ceramic's dielectric
energy.

Construct `ReactivePiezoelectricTransient::new` from an accepted
`TetPiezoelectricTransient`, the source-bound `ReactiveElectrodeLoad` list,
matching `ReactiveElectrodeState` values, and any existing memoryless
resistors. A directly connected capacitor must already have charge consistent
with the body's terminal voltage. Construction does not project charge, flux,
displacement or velocity to make an inconsistent connection appear valid.

`step_candidate` eliminates the electrical midpoint equations into the
existing sparse solid/dielectric system. The same solve determines mechanical
motion, electrode charge, winding current and capacitor voltage. A successful
candidate owns all updated histories together. Rejection or cancellation
leaves the accepted histories unchanged.

Each accepted interval checks

```text
change(body kinetic + elastic + electric energy)
  + change(external capacitor + winding energy)
  + conductor Joule heat
  + signed acoustic outward work, when a guide is attached
  = mechanical work + external electrical source work.
```

The reports retain per-branch currents, transferred charge, storage at both
endpoints and heat. Heat belongs to the actual series or parallel conductor,
not automatically to the ceramic. These branches keep their source conditions
fixed and do not yet update the reactive elements' temperature.

Resonant shunting follows the established physical mechanism described by
Hagood and von Flotow, *Damping of structural vibrations with piezoelectric
materials and passive electrical networks* (1991),
[DOI 10.1016/0022-460X(91)90762-9](https://doi.org/10.1016/0022-460X(91)90762-9).
An inductor adds an electrical resonance; it does not improve vibration at
every frequency. Tuning must be assessed against the coupled response.

## Acoustic sensing with retained R/L/C storage

Run `cargo run -p fs-couple --example piezoelectric_plate -- --acoustic`.
This mode uses the same nominal PZT-5A tensor, mesh and six supports as the
other modes. The pressure load is integrated over the actual top-face P1
triangles; no fitted sensor gain or independently imposed mechanical force
replaces the reciprocal pressure/volume-flow map.

The declared circuit has the existing 1 kΩ conductor, ideal 1.256637 mH winding,
and a separate 1 nF capacitor derived from an authored permittivity and gap
geometry. The winding retains flux linkage and the capacitor retains charge.
Only the named external conductor produces heat; the ceramic and fluid remain
lossless in this demonstration. These are source-retaining authored circuit
models, not measured component specifications.

The fluid is explicitly authored with density 1000 kg/m³ and sound speed
1480 m/s. It occupies a lossless semi-infinite 10 × 10 mm rigid uniform guide.
The first transverse cutoff of that declared guide is 74 kHz; the incident
tone is 5 kHz with 1000 Pa peak pressure. This is not a measured water model,
free-field radiation, a finite cavity or a broadband transient prediction.

The initial displacement, velocity, electrode charge, winding flux and
capacitor charge come from the coupled periodic harmonic solution. Preparation
uses the implicit-midpoint frequency map `ωd = 2 tan(ω dt/2)/dt` and the
corresponding midpoint pressure scaling. This initializes the time-discrete
single-tone history without introducing a zero-start broadband excitation.
The declared interface bandwidth is 10 kHz; timestep size is not used as a
spectral certificate. There are 64 steps per cycle and two cycles; refine the
mesh and time sampling before drawing a design conclusion.

Eight CSV rows report time, electrode voltage, mean thickness motion, external
reactive energy, cumulative conductor heat, incoming/outgoing wave energies,
signed net outward work and the cumulative balance. Expected checks are:

- Voltage and motion oscillate; external reactive storage remains nonnegative.
- Heat and the separate incoming/outgoing wave energies accumulate positively.
- At 0.4 ms the prescribed incoming wave has supplied **1.351351351 × 10⁻⁸ J**
  to the guide section, from `duration × peak_pressure² × area / (2 ρ c)`.
  Most can be reflected; this is not the energy absorbed by the ceramic.
- Over a complete periodic cycle, net acoustic work is inward and pays the
  conductor heat. During a partial cycle either sign is possible as stored
  energy is exchanged with the guide.
- `balance_j` is `body energy + external reactive energy + cumulative heat
  + cumulative signed outward work − initial total energy`. It should remain
  small; it is printed without rescaling or zeroing. Each native step must
  satisfy the declared 10⁻¹⁸ J balance tolerance before publication.

The entry point is
`ReactivePiezoelectricTransient::step_with_acoustic_ports_candidate`.
It accepts optional moving-support velocities, retains original guide wave
reports, and publishes all body and electrical histories together. Empty port
lists preserve the existing non-acoustic behavior.

## Harmonic sensing, actuation and electrical loading

`TetPiezoelectricModel::solve_harmonic` uses the same prepared mechanical
stiffness K, consistent mass M, reciprocal coupling G and dielectric stiffness
Kₑ. Inputs and outputs are **peak phasors under exp(−iωt)**; ω is in rad/s,
not Hz. Convert cyclic frequency with ω = 2πf.

`HarmonicElectrodeCondition::Voltage` prescribes a potential phasor.
`Current` prescribes an external current into an electrode, before shunt
outflow. Zero current gives the AC open-circuit condition. It is a total
electrode condition, not a pointwise assertion that electric displacement
vanishes on each face. Non-electrode surfaces remain insulating.

`HarmonicElectrodeLoad` accepts the same original resistive and reactive
source definitions used by the transient path. For this time convention,
R has admittance 1/R, a series R–L has 1/(R − iωL), and a capacitor has −iωC.
Current, voltage and charge orientations all follow the ordered electrode pair.

The coupled complex symmetric pencil is represented as a sparse real
symmetric system for the native MINRES solver. No dense factorization or modal
truncation is introduced. Positive stiffness/inertia magnitudes supply the
coordinate scales, so a cancelling mechanical dynamic diagonal can still be
solved when the coupled electrical system is regular. Iteration and real
coordinate limits remain explicit. The returned numerical residual is not
a forward-error bound near an ill-conditioned resonance.

The returned response contains nodal motion and velocity, electric potentials,
cell strain/stress/field/displacement, support reactions, electrode free charge
and source current, and per-load current and heat. Mean powers use
Re(effort × conjugate(flow))/2 because the phasors are peak amplitudes.
Applied mechanical work includes any prescribed moving support. Electrical
work includes electrode and volume-charge sources. Potential differences
remove arbitrary common voltage from all work accounts.

Small-signal admission checks the maximum strain norm and field magnitude over
the **whole cycle**, rather than only each quadrature. Both quadratures also
pass the existing constitutive reconstruction, force/charge residual and
storage checks. The final mechanical plus electrical average supply must
equal the actual external load heat within the caller's power tolerance.
With intrinsic relaxation enabled, the same account also includes material heat.

## Intrinsic mechanical and dielectric relaxation

`TetPiezoelectricModel::solve_harmonic_with_relaxation` adds source-backed
ceramic loss and dispersion to that same spatial solve. It accepts one uniform
`RelaxingTensorPiezoelectric` or one per original tetrahedron. Every law must
retain the exact original base tensor identity and mesh frame.

The material owner is
`fs_material::piezoelectric::tensor::relaxation::RelaxingTensorPiezoelectric`.
Its `resolve` constructor takes the original base tensor and card, explicitly
named native `TensorBranchProperties`, directed `DirectedDebyeProperties`, and
the original material selection policy. All mechanical branch entries and all
three dielectric direction components must be provided; zero coefficients are
explicit source data. The direction is supplied in the base source frame and
co-rotates with the specimen.

**The base coefficients have specific limiting meanings:** mechanical C₀ is
the relaxed constant-field stiffness, and ε∞ is the instantaneous clamped
permittivity. Adding a Debye increment to an already relaxed permittivity
would count polarization twice. A single-frequency loss tangent or mechanical
quality factor does not identify the required spectrum. The nominal PZT-5A
example therefore keeps intrinsic loss disabled; its cited data do not supply
these relaxation parameters.

For peak phasors under exp(−iωt), each mechanical Maxwell branch adds
`C_j (-iωτ_j)/(1-iωτ_j)`. Each directed dielectric pole adds
`Δε_j n_j n_jᵀ/(1-iωτ_j)`. Native positive-semidefinite mechanical factors
and positive dielectric increments make their heat nonnegative. The original
real piezoelectric coupling remains reciprocal; this model does not infer
intrinsic coupling loss. Source frequency applicability is enforced without
freezing a frequency-dependent measured curve into a constant parameter.

The solve evaluates each supplied material's frequency coefficients once,
adds its increments through the original mesh gradients, and reuses the same
mass, coupling, boundary reduction and sparse solver. Constitutive response
evaluation allocates nothing. Independent stress and polarization recovery
supplies the actual forces, electrode charge and source currents.

Each response cell reports `average_mechanical_heat_w`,
`average_dielectric_heat_w` and `average_memory_energy_j`. These quantities are
already integrated over that cell. Memory energy is included in
`average_stored_energy_j`; it must not be added again. The response's
`average_material_heat_w` sums the two intrinsic heat channels, while
`average_load_heat_w` remains the external conductors' heat. The complete
average supply must pay both. Coefficients stay at their held material state;
these heat reports do not themselves evolve the ceramic's temperature.

The source-backed model and its analytic checks are documented by the primary
[piezoelectric-loss theory](https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_modeling.05.109.html)
and [Debye dispersion theory](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html).

## Retained intrinsic memories in time

The same source-bound material law also extends the existing
`TetPiezoelectricTransient` at initialization:

```rust,ignore
let body = TetPiezoelectricTransient::new(/* original mesh state and sources */)?
    .with_relaxation(
        vec![relaxing_material],
        maximum_internal_variables,
        controls,
        cx,
    )?;
```

Supply one uniform law or one per original tetrahedron. The explicit scalar
budget counts six native Gram coordinates per Maxwell branch and one directed
polarization per Debye pole, for every cell. Material identity, mesh frame and
base coefficients must match the original model. A frequency-limited source
fit cannot authorize an arbitrary transient. Initial memories are zero; they
are not inferred from the initial displacement or field, and they cannot be
reset after time has advanced.

`step_candidate`, `step_with_loads_candidate` and retained reactive shunts all
use these histories in the same sparse coupled midpoint solve. Each process
is eliminated locally, so no material-memory degrees of freedom are added to
the global matrix. The mechanical companion contributes
`tau/(tau+dt/2) * C_branch`; the polarization companion contributes
`(dt/2)/(tau+dt/2) * delta_epsilon * n n^T`. The original consistent mass,
reciprocal coupling, field maps and electrode topology remain intact.

The native step reports `mechanical_material_heat_j`,
`dielectric_material_heat_j`, their sum `material_heat_j`, and before/after
`material_memory_energy_j`. `material_relaxation_steps()` exposes the same
heat and storage accounts per original cell. The accepted state's
`material_memory_energy_j()` is a subaccount already included in
`stored_energy_j()` and `total_energy_j()`. Conductor Joule heat remains a
separate physical account. Complete cell stress and electric displacement
include the retained material response.

The acceptance identity is

```text
delta(kinetic + elastic + electric + material memory)
  + intrinsic material heat = mechanical work + net electrode work.
```

Enclosing resistor/reactive reports also pay conductor heat and reactive
storage against actual external source work. Both material heat channels are
computed from positive midpoint memory/friction norms, independently of this
identity. Rejected candidates leave all accepted mechanical, electrical and
material histories unchanged. Implicit midpoint is second order and passive,
but very fast unresolved relaxation modes can alternate; time refinement is
still needed to resolve their dynamics.

### Finite ceramic heating

`step_with_material_thermal_loads_candidate` deposits intrinsic ceramic losses
and external resistor heat through the existing finite enthalpy solver. Supply
one `ThermalCeramicCell { piezoelectric_cell, thermal_cell }` binding for every
original tetrahedron and the original absolute-temperature axis. Multiple
tetrahedra may share a thermal cell when their summed reference volume and mass
match its complete caloric inventory. Duplicate, missing and foreign routes
refuse; the resistor-only API continues to require explicit ceramic ownership
before it can accept a body with intrinsic memories.

Each receiver must retain the original ceramic card, density claim, complete
nonthermal source conditions and actual caloric law. Its chart remains fully
solid at the original reference density. The source owner checks every accepted
and trial temperature against all original base and relaxation claims. A changed
temperature refuses any coefficient with temperature dependence, including a
flat sampled segment whose wider constancy has not been established. This
models finite heating where the original claims support fixed coefficients.
Temperature-driven stiffness, coupling, permittivity, polarization rates, thermal
strain and depoling require their own source admission and material work.

Ceramic losses are internal heat transfer into those named cells, never an
extra external supply. The same enthalpy residual retains thermal conduction,
source-dependent resistor feedback, caloric storage and the combined physical
energy gate. The report preserves `ceramic_destinations` beside the original
resistor destinations; heat amounts remain in the native per-cell material
accounts. A rejected thermal candidate changes neither thermal nor material
history. With no resistor feedback, the thermal Newton trials share one
immutable spatial solve while still checking every trial's actual sources.
With endpoint-dependent thermal or resistor coefficients, the combined thermal
path remains generally first order; energy closure does not certify time accuracy.

### Thermal expansion, pyroelectricity and reversible heat

The separate full-tensor thermal path supports explicit thermal expansion and
pyroelectric response through
`ThermalTensorPiezoelectric::with_thermal_bias(reference_temperature_k, selection)`.
The same source card must declare all six
`piezoelectric_thermal_expansion_mandel_1..6` coefficients in 1/K and all three
`piezoelectric_pyroelectric_clamped_1..3` coefficients in C/(m² K), including
zeros. They are expressed in the original tensor's source frame and must remain
constant over the admitted temperature interval. The reference temperature
defines zero thermal strain and zero polarization increment.

The material potential uses total strain and clamped pyroelectric coefficients:

```text
sigma = C (S - alpha (T - Tref)) - transpose(e) E
D     = e S + epsilon E + p_clamped (T - Tref)
```

Free thermal expansion therefore contributes the secondary response `e alpha`;
a measured total pyroelectric coefficient must first be converted to the
clamped convention. The caloric chart supplies the stress-free, zero-field
reference heat capacity. Constraint-dependent field energy and entropy are
accounted by the coupled law, including reversible electrocaloric heating or
cooling. They are not assigned to a Joule-loss channel.

`TetPiezoelectricModel::assemble_thermal_materials` accepts resolved thermal
samples for static fields or fixed-temperature transients. For finite thermal
evolution, use `ThermalTensorPiezoelectricBody` with the original complete
caloric inventory. A bare initial spatial model is compatible only where the
actual thermal offsets vanish; a nonzero initial bias needs a matching thermal
model. Its `step_with_acoustic_ports_candidate` can retain an original guide
trace while temperature, polarization and thermal strain evolve. Support work,
electrical work, sound, resistor heat and reversible thermal exchange retain
separate accounts.

The nominal plate example does not invent the missing expansion/pyroelectric
source claims, so those effects remain disabled in its demonstration modes.
Synthetic analytic regressions exercise the new paths. Finite thermal feedback
remains generally first order and uses one uniform ceramic temperature.
Transient intrinsic memories with nonzero thermal bias still require an
accounted memory reference.

### Frequency response around a heated or DC-biased equilibrium

Prepare `model.harmonic_operating_point(...)` with the actual static supports,
electrode voltages/charges, mechanical forces and volume charge. The returned
immutable owner retains its solved fields and the exact original model. Call
its `solve_harmonic(...)` with angular frequency, acoustic ports, AC supports,
electrodes, loads, forces, charges and the original harmonic controls. The
support DOFs and electrode node sets/order stay fixed; DC and AC electrode
source kinds may differ. Use empty port/load slices when none are attached.

`point.state()` contains the DC fields. `response.harmonic` contains only the
oscillating perturbation. `average_total_stored_energy_j` adds the steady body
free energy once to the incremental mean storage. Each cell's total strain,
thermal elastic strain and electric field are checked over the full combined
cycle, including extrema between the real/imaginary quadrature samples.

The response also exposes the DC source currents needed to maintain that bias
against physical shunts, their DC Joule heat and external capacitor/inductor
storage. Those accounts remain separate from AC power and ceramic storage.
`solve_harmonic_with_inductor_bias(...)` requires `(load_index, current_a)` for
every ideal lossless winding, including explicit zero; nonzero DC voltage
across such a winding cannot be stationary and refuses. Finite-resistance RL
currents instead follow the actual bias voltage and winding resistance.

This is a held-temperature, small-strain affine material model. DC offsets do
not create geometric stress stiffening or nonlinear voltage tuning. The original
`model.solve_harmonic(...)` remains the zero-bias entry point and refuses a
thermally biased model unless its equilibrium has been prepared explicitly.

## Verification and boundaries

Run the focused native spatial suite:

```bash
cargo test -p fs-couple --test piezoelectric_spatial
```

The new tests use explicitly synthetic source data and independent RC/LC/RLC
and retained-coordinate equations. They exercise midpoint convergence,
electrical and mechanical resonance, Maxwell decay, Debye charging, joint
mechanical/polarization/shunt trajectories, separate ceramic/conductor heating,
complete thermal source and inventory admission, mean-power factors and signs, open/short
limits, orientation/gauge/phase invariance, full-cycle amplitude limits and
transactional refusal. They verify implementation behavior, not measured
material accuracy.

The spatial shunts currently connect electrodes within one admitted voltage
reference component. Arbitrary internal circuit nodes, transformers, mutual
inductance, nonlinear switching and reactive material thermal feedback remain
outside this path. The existing longitudinal `PiezoelectricCircuit` remains
the separate reduced rod adapter for the native general circuit graph.

These models do not supply ferroelectric hysteresis, saturation,
depoling, breakdown, fringing, free-field acoustic radiation or dissipative piezoelectric
cross coupling. Intrinsic relaxation is available in harmonic and transient paths.
The harmonic result is a periodic particular solution; it does not establish
startup decay or uniqueness of an unexcited lossless mode. The plate mesh and
time/frequency sampling must be refined for a concrete design question.

