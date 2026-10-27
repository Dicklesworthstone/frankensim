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

## Heat the spatial valve from a finite source gas

Add `--finite-gas VOLUME_M3 GAS_DISTANCE_M STEPS` and the shared explicit gas
source options to use a finite thermal reservoir. The gas exchanges heat with
the actual exterior triangles of the same 48-cell solid used by the mechanical
reducer. The clamped root is insulated. Shared interior faces carry heat
between neighboring cells. Each cell retains its own temperature; a mean
temperature is never substituted for the structural field.

For the supplied hot-nitrogen record, compile the existing source packs and run:

```bash
run_dir="$(mktemp -d)"
demo=crates/fs-couple/examples/heated_valve_sound
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen/manifest.tsv \
  --out "$run_dir/nitrogen.species"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/manifest.tsv \
  --out "$run_dir/nitrogen.nasa9"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-transport-svehla1995/manifest.tsv \
  --out "$run_dir/nitrogen.transport"

cargo run --release -p fs-couple --example heated_valve_sound -- \
  --demo 1 0 10 1 "$run_dir/finite-valve" \
  --finite-gas 0.0001 0.001 32 \
  --gas-record "$demo/nitrogen.gas" \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen
```

This starts the solid at 300 K and 0.1 litre of nitrogen at 340 K and 100000 Pa.
The imposed solid power is zero: the gas cools as the solid warms. Change only
the volume to `0.0004` and choose a new output prefix to compare a reservoir
with four times the gas inventory. The initial acoustic gas is identical in
both runs; the subsequent temperature, source elasticity, aperture and sound
are calculated from each reservoir's actual energy history.

The original pressure CSV and stereo WAV are accompanied by
`PREFIX.thermal.csv`. Its initial row and each accepted thermal endpoint report
solid minimum/maximum temperature, gas temperature/pressure, cumulative external
heat, solid storage change, gas storage change, solid-to-gas heat and the actual
combined storage residual. With zero external power, the solid gains the
energy lost by the gas. Every complete physical calculation and sound render
still precedes file creation. File I/O itself is not an atomic multi-file
transaction.

The finite-gas mode explicitly admits a **fully solid 280–400 K interval**,
including a margin around the 300 K initial state for numerical Jacobian
probes. Packs must supply their own support across that interval; there is no
extrapolation. In addition to the mechanical/caloric properties, the bounded
native example requires a nonnegative scalar-valued `thermal_conductivity`
source at every actual local temperature and numerical probe. Constant claims,
temperature curves and supported fixed-context response surfaces use the same
material resolver. The source is selected once at the original state, then
refreshed at every trial endpoint. Internal links and gas-facing patches consume
the **same cell material snapshot**; conductivity is not frozen at 300 K or
reselected separately at a boundary. An explicit zero coefficient insulates;
missing or unsupported data refuses the entire candidate.

The manufactured demo defaults to 0.2 W/(m K). To compare a temperature-dependent
source with exactly the same initial conductivity, repeat the command above
with a new output prefix and add:

```text
--demo-conductivity 0.1 0.2 2.0
```

These three values specify conductivity at 280, 300 and 400 K, with the ordinary
source pack's linear-inside interpolation. The option requires `--demo` and
`--finite-gas`; it declares manufactured material data and does not override an
imported pack. All reference mechanical/caloric coefficients, geometry, gas,
heating and acoustic excitation remain the same. The calculation propagates
the changed conductivity through the spatial heat field, deformation, valve
frequency and actual pressure. `--pack` uses the supplied conductivity law
directly, with no demonstration override.

Interior conductances use each endpoint cell's actual source conductivity,
shared triangle area and the two centroid-to-face normal distances. The existing
prepared material-network owner lowers those two material layers in series
inside the **same implicit residual** as finite gas heat and local enthalpy
storage. Different cell temperatures therefore give different half-cell
coefficients even for a homogeneous specimen. This is an explicit coarse two-point
thermal network. The tetrahedral mesh is not orthogonal, so the network does
not claim a consistent continuum heat-equation discretization. It does conserve
each signed internal transfer. Boundary conductance places the source gas and
solid layers in series using the declared gas-layer distance and the cell's
centroid-to-face distance. A floating-point area summation allowance is kept
insulated; it creates no transfer patch.

Volume and layer distance must each be finite and in `(0,1]` SI units; steps
must be in `1..=1024`. Each interval uses the original backward-Euler thermal
owner and NASA9 gas internal-energy inversion in one bounded solve. Declared
numerical controls are: 4096 complete source trials, 24 Newton iterations,
`1e-7 J` finite-difference scale, `1e-9 J` absolute nonlinear tolerance,
`1e-8 J` per-cell and gas tolerances, `1e-7 J` thermal tolerance, and
`2e-7 J` combined tolerance. Gas inversion uses 280–400 K, `1e-7 J/kg`, at most
128 iterations, and an admitted temperature change of at most 20 K per step.
Domain exits or budget/energy failures refuse the run; more steps can resolve
a too-large physical interval. These controls do not certify time accuracy.

The cold sound uses the initial gas state; the heated sound uses the final gas
state. In each frozen experiment that same supplied medium reaches excitation,
bore, radiation and propagation. The finite volume determines **thermal storage**;
it is not a closed acoustic chamber geometry or a pressure load on the solid.
Gas mass/composition/volume and reference solid mass stay fixed. Moving surfaces,
gas flow and pressure-volume work need their own closures.

## Retained inertial motion during heating

The finite-gas mode can also advance the solid's actual displacement, velocity
and acceleration as the spatial temperature field changes. Add:

```text
--inertial-history MAX_STEP_S MASS_DAMPING_PER_S
```

`MAX_STEP_S` is a finite positive ceiling on every thermal/mechanical interval,
in seconds. `MASS_DAMPING_PER_S` is an explicitly authored nonnegative
mass-proportional damping coefficient, in inverse seconds. The option requires
`--finite-gas`. Choose `HEAT_SECONDS / STEPS <= MAX_STEP_S`; a larger interval
refuses and asks for more finite-gas steps or a shorter heating duration. The
ceiling is an admission control, not a temporal-error certificate: resolve the
fastest inertial behavior of interest and check timestep refinement.

For example, a short manufactured heating pulse can use `POWER_W=5000`,
`HEAT_SECONDS=0.002`, eight finite-gas steps, and:

```text
--inertial-history 0.0003 20
```

An optional explicit Rayleigh stiffness coefficient is now supported through
material softening and stiffening:

```text
--inertial-history 0.0003 20 --inertial-stiffness-damping 0.00001
```

The extra coefficient is in seconds, must be finite and nonnegative, and
requires `--inertial-history`. Either option order is accepted. Omitting it
retains zero stiffness damping and the original mechanical calculation.
These are authored mechanical model inputs, not source-calibrated loss spectra
or a conversion of the acoustic loss factor. Both coefficients appear in the
motion CSV header.

The original thermal solve first retains each accepted source state. The
existing sparse solid owner then advances one trajectory from stress-free
rest on those **same physical endpoint times**, using each cell's refreshed
stiffness and thermal strain. Displacement and velocity continue through the
whole trajectory; they are not reset to successive equilibria. The reference
mesh, solid mass and homogeneous clamped supports stay fixed. For changing
stiffness and nonzero stiffness damping, the native midpoint effort uses
`Cbar = eta*M + zeta*(K_old+K_next)/2`, so loss is
`dt*vbar^T*Cbar*vbar`. Physical endpoint acceleration is recovered from the
new endpoint operators. Do not assume its endpoint average equals the discrete
velocity increment in this variable-damping path. Fixed-stiffness and zero
stiffness-damping paths retain ordinary average-acceleration Newmark behavior.
No authored acoustic loss is copied into the mechanical history.

The additional `OUTPUT_PREFIX.motion.csv` contains a row for every node at
every accepted thermal time, including the initial state. Columns report
reference coordinates; displacement in metres; velocity in metres/second;
acceleration in metres/second²; kinetic and shifted elastic energy; cumulative
external, thermal-expansion and stiffness work; damping loss; and the cumulative
mechanical balance residual, all energies in joules. Global energy accounts are
repeated on each node row and must not be summed over nodes.

This is **prescribed one-way thermal loading**. Signed expansion/stiffness work
is accounted by the native mechanical owner, but is not debited from caloric
storage. Mechanical damping is not silently credited as heat. Reciprocal
thermoelastic heat needs an admitted closure. The existing pressure CSV and WAV
remain the two frozen-state acoustic experiments; this option does not claim
that a live vibrating modal state has been transferred into those recordings.
All thermal and mechanical admission and sound calculations finish before any
output file is created, and all output paths must be new.

## Physical calculation

The declared specimen is a 10 × 12 × 4 mm slab, with 27 reference nodes and
48 affine tetrahedra. Its root is clamped. A fixed plane lies 0.4 mm below its
lower surface. The tip half of that underside supplies the opening observation
and a geometric unit-pressure load; the original load assembler determines
the signed force-per-pascal pattern. All geometry is declared in the example.
The mesh is deliberately coarse and does not establish continuum accuracy.

Every tetrahedron has its own finite mass and caloric state. In the default
mode, uniform volumetric heating on the homogeneous insulated specimen gives the same temperature in
every cell and therefore zero internal conductive flux. This exact symmetry
is why the default mode needs no thermal conductance. The finite-gas mode uses
the explicit transport network above to produce a nonuniform thermal field
from the same finite source gas.

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

The finite-gas regressions change only reservoir volume and check nonuniform
solid heating, fixed gas mass through `p/T`, equal-and-opposite endpoint storage,
different valve frequency and actual pressure, deterministic replay and input
or cancellation refusal.

An additional material-swap regression changes only the conductivity's
temperature dependence while retaining its initial value. It checks unchanged
cold pressure, changed heated deformation/frequency/pressure, conserved combined
storage and deterministic replay. The library also covers heterogeneous original
material cards against independently written nonlinear cell/gas energy equations.

The inertial regressions exercise a finite heating pulse with changing source
stiffness, retained Newmark kinematics, identical thermal/mechanical clocks,
zero supports, explicit work and energy balance, actual inertial displacement,
deterministic replay, and timestep/control/cancellation refusal. Additional
regressions define explicit stiffness-damping input admission and check changed
motion/loss with unchanged prescribed heat and frozen acoustic observations.
