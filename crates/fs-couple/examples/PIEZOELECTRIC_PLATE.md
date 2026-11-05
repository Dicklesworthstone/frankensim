# Spatial piezoelectric plate, reactive shunts and frequency response

Run the finite-plate example with the repository's pinned toolchain:

```bash
cargo run -p fs-couple --example piezoelectric_plate
cargo run -p fs-couple --example piezoelectric_plate -- --harmonic
```

The first command computes static actuation, the matching isolated-electrode
state, a voltage ramp, open-circuit ringing and resistively loaded motion.
The harmonic mode reports a small frequency sweep with open, resistive and
resistor–inductor loads. Voltage, current and deformation come from the coupled
field solution. The table also reports mean conductor heat and the residual
of the independently reconstructed power balance.

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

Intrinsic heat is reported but not automatically deposited into a temperature
field. The existing resistor-only thermal adapter refuses bodies with these
memories because it has no named ceramic heat destinations. Temperature-driven
coefficient changes require their own source admission and material work.

## Verification and boundaries

Run the focused native spatial suite:

```bash
cargo test -p fs-couple --test piezoelectric_spatial
```

The new tests use explicitly synthetic source data and independent RC/LC/RLC
and retained-coordinate equations. They exercise midpoint convergence,
electrical and mechanical resonance, Maxwell decay, Debye charging, joint
mechanical/polarization/shunt trajectories, mean-power factors and signs, open/short
limits, orientation/gauge/phase invariance, full-cycle amplitude limits and
transactional refusal. They verify implementation behavior, not measured
material accuracy.

The spatial shunts currently connect electrodes within one admitted voltage
reference component. Arbitrary internal circuit nodes, transformers, mutual
inductance, nonlinear switching and reactive material thermal feedback remain
outside this path. The existing longitudinal `PiezoelectricCircuit` remains
the separate reduced rod adapter for the native general circuit graph.

These models do not supply ferroelectric hysteresis, saturation,
depoling, breakdown, fringing, acoustic radiation or dissipative piezoelectric
cross coupling. Intrinsic relaxation is available in harmonic and transient paths.
The harmonic result is a periodic particular solution; it does not establish
startup decay or uniqueness of an unexcited lossless mode. The plate mesh and
time/frequency sampling must be refined for a concrete design question.

