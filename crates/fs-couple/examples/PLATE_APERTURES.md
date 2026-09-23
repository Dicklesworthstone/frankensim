# Plate-derived moving valves and coupled tube motion

The `bernoulli_aperture::plate` path derives a moving slit from an actual flat
triangle plate chart, its supports, and its numeric or material-bound sections.
It then drives the **existing nonlinear midpoint aperture** and the existing
reciprocal characteristic tube/network. It does not prescribe an audio frequency,
fit an effective mass, or animate the surface after a quasistatic flow solve.

Run an authored numeric specimen, then change only its density:

```bash
cargo run --release -p fs-couple --example plate_aperture -- \
  4000000000 900 0.0003 > /tmp/plate-valve-light.csv
cargo run --release -p fs-couple --example plate_aperture -- \
  4000000000 1800 0.0003 > /tmp/plate-valve-heavy.csv
```

Arguments are Young's modulus [Pa], density [kg/m³] and thickness [m]. With no
arguments, the first line's values are used. The example explicitly constructs
an authored 25 × 10 mm plate, a clamped root and a 10 mm free-edge slit with a
0.2 mm rest gap. Poisson ratio is 0.3 and viscous damping ratio is 0.02. These are
**synthetic numeric inputs, not a measured cane/material card**. The estimated
one-coordinate lay has explicit stiffness, exponent and loss in the source;
none is secretly selected by the runtime. The same 5 Pa pressure for 1,024
steps, followed by release, drives every substitution over 4,096 steps at 100 kHz.
No pressure rescaling or output normalization conceals the change in mechanics.

The tube radius is 7 mm, requested length 250 mm and memoryless terminal
reflectance -0.8. Density and sound speed come from the existing moist-air model
at 293.15 K, 101325 Pa and zero relative humidity. The declared length-error
allowance is 2 mm; both requested and integer-transit represented lengths are
printed. Ambient gas state does **not** silently modify the solid specimen.
This example is a lossless tube interior with a passive terminal, not a claim
that a realistic instrument has no viscothermal loss. Existing network boundary
loss, compliant-wall and cavity owners can host the same bound moving valve.

CSV contains physical opening/velocity, reconstructed tip displacement/velocity,
maximum retained slope, swept/jet flows, inlet pressure and whole-system
energy/work/loss. Mechanical state is at step end; the explicitly labelled
pressure/flow quantities belong to the midpoint. The pressure is internal to
the tube, **not an exterior microphone or a radiated-audio prediction**.

## Geometry and material input

`PlateApertureReduction::from_chart` accepts an `fs_plate::PlateChart` with
arbitrary admitted flat triangles, explicit supports and per-triangle sections.
`from_material_chart` accepts the immutable `ResolvedPlateChart` produced by
`thin_plate::compile_plate_material_chart`. That path retains the original
regional material states, receipts, thickness/mass constraints and axis angles.
Both call the same existing plate assembly and modal eigensolver. Orthotropic,
rotated, tapered or regionally heterogeneous sections are not collapsed into a
uniform bounding rectangle.

Supply `PlateApertureOptions` explicitly: assembly/support choices, eigenvalue
search window and selected index, existing eigensolver budgets, node/triangle
caps, unique boundary edges forming the slit, rest gap, damping ratio, allowed
slit nonuniformity and maximum linear slope. Neither a support nor a slit is
inferred from an instrument name. Source state queries are the caller's choice;
a numeric chart alone acquires no material-identification authority.

The reduction keeps the complete source chart and material binding; it exposes
its mass, stiffness, pressure area, slit width, closing pressure and normalized
physical shape. The eigenvalue interval is evidence about the **discrete pencil**,
not a bound on omitted physical modes or a measured-reed fidelity certificate.
Damping is an explicit independent input. An elastic material card does not
provide a viscous loss law by implication.

## One work-conjugate coordinate

Let phi be the retained mode and s its length-weighted mean displacement on
the specified slit. For mean opening change x = s q:

```
M = (phiᵀ M_plate phi) / s²
K = M lambda
A = integral_over_plate(phi_displacement dA) / s
P_close = K H / A
F_pressure = -A deltaP
U_swept = -A x_dot
```

The pressure integral uses the plate's P1 displacement trace over its actual
triangles. The sign/normalization of the eigensolver's phi cancels. Geometry
sets width by the sum of the supplied edge lengths. The same area supplies
pressure work and swept flow; there is no legacy 25 mm effective-face fallback.
A non-closing or zero-slit projection refuses rather than being made positive
by an absolute value. Different densities at equal elasticity/geometry retain
static closure but change inertia and the **coupled** pressure/flow trajectory.

## Retained runtime and limits

Use `DynamicAperture::from_plate(reduction, rho, impedance, dt, max_steps,
initial_state, lay)` to keep the specimen attached. Then pass that valve to
`ApertureTube::new` or `ApertureNetwork::new`. No extra time integrator is
introduced. Scalar-only `DynamicAperture::new` remains a separate authored path;
`plate_reduction()` returns `None` there, so it cannot impersonate a bound specimen.
A `dynamic_spec()` export alone likewise does not retain a geometry guard.

The bound path checks both initial and candidate openings against the declared
nodal-rotation/P1-slope limit. A failed candidate changes neither valve state
nor composed propagation history. State/energy/clock and material bindings
survive cancellation and explicit budget extension. `nodal_motion()` reconstructs
physical displacement and rotation from the current opening and velocity; it
does not amplify or step the model. For this one linear mode the allowed opening
interval is convex, so accepted endpoints also bound the midpoint geometry.

This is a **single-mode linear plate, uniform-pressure loading, mean slit gap and
one generalized compliant lay contact**. It is not a distributed moving closure,
nonlinear reed geometry, wet-cane law, fluid-loaded eigensolve, viscoelastic memory,
a multimode truncation certificate or instrument reconstruction. The caller must
state the allowed slit-shape variation; a mean-gap approximation is not exact
when nodal gaps differ. The existing Bernoulli relation itself retains its model
limits. No physical limit is relaxed to make a material or mode pass.

Reduction is offline. Cancellation is polled around the existing assembly and
eigensolve, not inside those owners, and during projection. The admitted runtime
still uses the owners' allocations and nonlinear iterations; no allocation-free
or real-time throughput claim is made.

Focused native targets:

```bash
cargo test --release -p fs-couple --test plate_aperture \
  --test dynamic_aperture --test aperture_tube --example plate_aperture
```

The eleven new tests cover mechanical/pressure projection, material and thickness
changes, actual two-way tube trajectories, work accounting, physical nodal motion,
source-receipt retention, slope rejection with exact retry, cancellation and
budgeted continuation. They do not replace measured-instrument validation.
