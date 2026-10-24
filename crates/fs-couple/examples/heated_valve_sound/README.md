# Finite heating changes a material valve and its sound

`heated_valve_sound` joins the existing local enthalpy history, source-backed
three-dimensional solid, thermally prestressed valve reduction, acoustic
assembly and physical-pressure WAV encoder. Changing the material pack changes
that same calculation. It does not select an instrument preset.

## Run

The supplied demonstration data is manufactured, not measured cane or a
qualified material. From the repository root, with the admitted dependencies:

```bash
cargo run --release -p fs-couple --example heated_valve_sound -- \
  --demo 1 1 10 1 /tmp/heated-valve-cp1
cargo run --release -p fs-couple --example heated_valve_sound -- \
  --demo 2 1 10 1 /tmp/heated-valve-cp2
```

The five values after `--demo` are heat-capacity multiplier, total heating
power in watts, heating duration in seconds, PCM full scale in pascals, and
output prefix. Each run writes a stereo WAV (left: cold; right: heated) and a
CSV with both pressure histories. They contain the same 4096-sample, 48 kHz
blowing experiment at the same observation location. No peak normalization or
pitch correction is applied. Output paths must be new. Both complete physical
calculations and WAV encoding succeed before any file is created.

The second run changes only the material's heat capacity. Its reference mass,
reference elastic law, heating energy, geometry, acoustic gas and excitation
are unchanged. The same finite energy produces a different temperature,
expansion and stiffness, which changes the actual calculated pressure.

For an existing canonical normalized material-card pack:

```bash
cargo run --release -p fs-couple --example heated_valve_sound -- \
  --pack /path/to/material.matpack 1 10 1 /tmp/heated-valve-source
```

The selected source must support `T` from 300 to 400 K for density and specific
heat capacity. It must also supply one continuous, unambiguous source for
Young's modulus, Poisson ratio and the instantaneous linear thermal-expansion
coefficient along the actual heating path. Their canonical names are `density`,
`specific_heat_capacity`, `young_modulus`, `poisson_ratio` and
`linear_thermal_expansion_coefficient`. The user must select a card whose
300–400 K interval is fully solid: that caloric chart is an explicit model
assumption, and this example does not infer melting or transformation laws
from a material label. Missing data, a source domain exit or ambiguous claims
refuse; no substitute material is chosen.

## Physical calculation

The declared specimen is a 10 × 12 × 4 mm slab, with 27 reference nodes and
48 affine tetrahedra. Its root is clamped. A fixed plane lies 0.4 mm below its
lower surface. The tip half of that underside supplies the opening observation
and a geometric unit-pressure load; the original load assembler determines
the signed force-per-pascal pattern. All geometry is declared in the example.
The mesh is deliberately coarse and does not establish continuum accuracy.

Every tetrahedron has its own finite mass and caloric state. Uniform volumetric
heating on this homogeneous, insulated specimen gives the same temperature in
every cell and therefore zero internal conductive flux. This exact symmetry
is why the example needs no guessed thermal conductance. Surface heating,
nonuniform source power and transient cooling require an explicitly resolved
transport network or spatial conduction model; those are not silently inferred.

The original caloric owner integrates heat. The source observer checks material
support and solid phase at every accepted thermal substep and event. Expansion
integrates from the original stress-free reference temperature. The solid owner
solves compatible displacement and stress, and its signed geometric stiffness
enters the same modal pencil as current source elasticity. The actual displaced
surface supplies the new gap. Reference mass is preserved even if a current
density claim differs from its initial value.

The entire tangent must remain positive definite, and the retained first mode
must carry at least 90% of both opening and pressure compliance, with no more
than 10% relative error in static pressure-induced opening. Inadequate modal
reduction, instability, excessive strain/displacement or closed aperture
refuses. The example's structural loss factor 0.1 is declared research input;
it is not inferred from the material name or presented as measured damping.

The admitted valve enters the existing material acoustic compiler. It replaces
all template reed mechanics with its calculated mass, stiffness, gap, width,
pressure area and damping. The unchanged 0.5 m cylindrical duct, gas at 300 K
and 101325 Pa, 100 Pa blowing-pressure attack and 1 m listener then produce
physical pressure samples through the ordinary realizer.

## Scope and accounting

These are two **frozen-state acoustic observations**, each starting from rest.
They do not transfer a vibrating modal state while its temperature changes.
Heating is sequential and quasistatic: its audit covers caloric storage and
applied heat; equilibrium elastic energy is retained separately. Reversible
thermoelastic heat, acoustic work returned to the heating model, source-derived
loss spectra, follower pressure loads, moving contact, finite strain and phase
change need their own admitted closures.

The existing acoustic realizer's compact-source observation and finite-mode
approximations remain. This demonstration does not establish a calibrated
instrument, a full exterior wave solution or measured material fidelity.

Focused regressions check finite-energy temperature against `Q/(m cp)`, reference
mass conservation, source-only heat-capacity changes reaching actual pressure,
deterministic replay and source/cancellation refusal. The library's
`solid_enthalpy::valve` tests additionally compare the thermal tangent and both
static pressure compliances with independent six-tetrahedron analytic integrals.
