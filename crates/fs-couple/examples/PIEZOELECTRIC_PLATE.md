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

## Verification and boundaries

Run the focused native spatial suite:

```bash
cargo test -p fs-couple --test piezoelectric_spatial
```

The new tests use explicitly synthetic source data and independent RC/LC/RLC
and retained-coordinate equations. They exercise midpoint convergence,
electrical and mechanical resonance, mean-power factors and signs, open/short
limits, orientation/gauge/phase invariance, full-cycle amplitude limits and
transactional refusal. They verify implementation behavior, not measured
material accuracy.

The spatial shunts currently connect electrodes within one admitted voltage
reference component. Arbitrary internal circuit nodes, transformers, mutual
inductance, nonlinear switching and reactive material thermal feedback remain
outside this path. The existing longitudinal `PiezoelectricCircuit` remains
the separate reduced rod adapter for the native general circuit graph.

Neither path added here supplies ferroelectric hysteresis, saturation,
depoling, breakdown, fringing, acoustic radiation or intrinsic ceramic loss.
The harmonic result is a periodic particular solution; it does not establish
startup decay or uniqueness of an unexcited lossless mode. The plate mesh and
time/frequency sampling must be refined for a concrete design question.
