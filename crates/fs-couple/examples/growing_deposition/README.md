# Growing thermoviscoelastic overhang

```sh
cargo test -p fs-solid bridge::growth
cargo test -p fs-couple --example growing_deposition
cargo run -p fs-couple --example growing_deposition
```

A 10 mm clamped seed grows to 20 mm through four timed cell births. The
thermal owner activates each supplied mass and incoming enthalpy exactly at
its scheduled time. Its actual CFL-limited temperature substeps drive the
existing thermoviscoelastic beam. No future cell has weight, stiffness,
thermal capacity, conductance or mechanical history before birth.

The public `fs_solid::bridge::growth::CellBirth` and
`Bridge::append_cantilever_controlled` / `append_cantilever` APIs continue a
born single-layer clamped seed at its current tip tangent. Existing section
history and geometry are unchanged at birth. New material has zero curvature
and axial strain, its own material law and thermal reference, and the real
incoming mass. Its weight enters the next positive-time equilibrium. Simply
extending a pin/roller span would move a support; that operation is refused.

The example's `driver::GrowingRoad` composes this library API with the existing
`DepositionObserver`. It validates an already active seed, explicit future
material bindings, body and cell identities, reference coordinates, volumes,
masses and temperatures. Every birth of the selected body must be bound;
other thermal bodies are explicitly thermal-only. There is one allowance for
all Gauss-point initializations and updates in each call, not a renewed budget
at each birth or thermal substep. Heat and mechanics publish together only
after the complete thermal advance and all mechanical/material checks pass.
The driver is retained next to this executable example, not exported as a
new fs-couple library module; existing crate dependencies are sufficient.

Restart by cloning at an existing step/event boundary and using the same
thermal grid. Ending exactly at a birth returns `birth_equilibrium_pending`:
the new domain is an initial condition, not yet a solved equilibrium. A birth
also clears the last-equilibrium report so stale reactions cannot be mistaken
for the weight of the enlarged domain. The printed source temperatures are
enthalpy-backed, and the reported tip position is the actual beam geometry.

## Boundaries

Reference cell centrelines lie along local +x at the section centroid, in
metres. Thermal cell geometry remains material-reference geometry; the beam
exposes current coordinates separately. New material is assumed perfectly
joined to the current tip, not constrained to a prescribed nozzle trajectory.
The seed must already satisfy the existing slenderness screen. This does not
model a short initial extrusion nub, tension/slack, a free liquid jet, dynamic
settling, large rotations, automatic right-end anchoring, impact or adhesion.
Thermal coefficients and all fixture properties are synthetic. Uniform-cell
thermal applicability, reference-area boundaries and source properties remain
caller obligations. No viscous heating or other mechanical energy is inserted
into enthalpy. Mechanical work/storage/dissipation checks and the thermal
mass/energy audit are separate, not a whole-system passivity certificate.

The example uses a mechanical update on each thermal substep. The spatial mesh,
birth intervals, thermal accuracy and material/time accuracy remain distinct
resolution choices. More material cells do not automatically improve birth
localization, and smaller solver residuals do not establish prediction accuracy.

## Verification status

Six library regressions check stress-free birth, history/tip continuity,
thermal references, enlarged elastic-span limits and spatial refinement,
restart, future-support exclusion and batch refusal. Six example-driver
regressions cover actual birth/substep timing, temperature and mass transfer,
same-grid restart, late physical failure, shared budgets and cancellation.
Native Rust compilation/tests were not run in the authoring environment:
Cargo and rustc are absent. Independent Python/SciPy equations exercised
1280 growing mechanical endpoints and the timed thermal/mechanical example
with time refinement; those checks are not execution of this Rust source or
calibration against a printed part.


## Declared drawing tension and slack

The growing-road driver now carries an explicit `fs_solid::bridge::drawing::Drawing`
state through every positive-time mechanical update. `Slack` preserves the prior
zero-axial-resultant path. `TensionN(value)` applies the same positive tensile
force to each accepted section and the beam-column equilibrium, so tension changes
both Maxwell stress history and sag. New cells remain stress-free at their exact
thermal birth time; the drawing load first acts on their next positive-duration
mechanical solve. Compression is refused rather than repurposed as an unverified
buckling model.

The executable advances otherwise identical slack and 2 mN synthetic-tension
roads on the same thermal schedule and checks that their enthalpy/cell histories
remain identical while the tensioned geometry sags less. The force is an explicit
input, not inferred from nozzle speed, extrusion pressure, or filament feed.
