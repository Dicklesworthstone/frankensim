# Deposited bridge sag and later-pass contact

`fs_solid::bridge` advances the spatial shape of a slender already-spanning road
under gravity and signed cellwise dead loads. Material history, thermal strain,
implicit tangents and discrete energy accounting remain owned by the existing
`deposition::LayeredStrip` / `fs_material::thermovisco` implementation. No second
Maxwell model, dense global matrix or finite-difference constitutive solver is
introduced. This is a bounded part of issue #1 R5, not completion of that issue.

## Physical scope

The local reference x axis is horizontal and z points up. Choose a left pin/right
roller or a left clamp/free right end. Axial resultant is zero: thermal contraction
is free, and no drawing tension, slack, membrane stiffening, inertia, shear or
molten-strand law is guessed. Supply an initially complete span; unsupported
moving-front deposition is not modeled. `deposit_bonded_layer` explicitly assumes
perfect bonding of a complete rectangular road onto the current section. Existing
stress/history is preserved, new points begin stress-free, and only born reference
mass contributes gravity. Disconnected or overlapping sections refuse.

Temperatures are keyed by `(axial_cell_id, layer_id)`. Each cell has constant
curvature and strain. Integrating curvature gives a C1 piecewise-quadratic shape;
rotation and deflection extrema include interior stationary points. The exact
cell-average equilibrium bending moment is the generalized load, not the midpoint
value. Consequently section mechanical work equals the integral of the declared
piecewise-constant dead load against this same displacement field. Global work
and storage are rechecked independently from the summed local material accounts.
The discrete balance separates material dissipation, backward-Euler numerical
dissipation, thermal-eigenstrain and modulus work. It is not a full entropy model.

Background: the Euler-Bernoulli weak form and its two integrations by parts are
summarized by TU Delft's *Computational Modelling*, section 4.1:
https://teachbooks.tudelft.nl/computational-modelling/structural_linear/euler_bernouilli.html
The moment sign here follows the existing section convention `epsilon=e0-z*kappa`
and `M=-integral(z*sigma dA)`, so `M=EI*w''` for a uniform symmetric section.
The implementation is an integral/force-based piecewise-curvature discretization,
not the textbook's cubic displacement element.

## Downstream contact geometry

`next_pass_contact` checks a centered straight later-pass footprint against the
actual top fiber. It inverts the deformed x map, including axial strain and fiber
rotation, and solves the resulting quadratic vertical gap using all polynomial
extrema and roots, not endpoint-only or midpoint-only sampling. An explicit
capture distance selects near-contact; it does not establish adhesion or force.
Excessive penetration refuses the whole query. Missing bridge coverage remains
reported as uncovered length, not proof that no other body exists there.

Each row keeps reference area and current horizontal projected area separately.
A reference-area thermal discretization must use the former, with independently
supplied conductance per area. Neither value is an exact curved surface area.
This geometry query supplies no force, nozzle impact prediction or strength
mapping. Separately declared mechanical stops can use the support solve below;
captured area is never silently converted to a support force.

```sh
cargo test -p fs-solid bridge
```

Thermal consumers may bind captured reference areas to their matching cell IDs
and an independently declared conductance per original area. Do not add a later
road to the perfectly bonded beam merely because it falls inside the geometric
capture distance. Contact-force transfer and strength need their own laws.

## Solved frictionless supports

`Bridge::advance_supported_controlled` and its executor-scoped companion
`advance_supported_to` solve upward reactions at declared `supports::Stop`
stations. An open stop carries zero force; loading can engage it and unloading
can lift the beam away. Reactions change actual section moments, material
history and integrated geometry, rather than clamping a displayed displacement.
The built-in pin/roller or clamp remains in place and is bilateral.

Each stop supplies a stable ID, reference x station and minimum allowed
reference-line displacement w in metres. For a contacted fiber at reference z,
subtract z from its support surface elevation. This is a material-reference
station with free tangential motion, not inversion of a fixed world-x obstacle.
Contact is enforced only at supplied stations: a beam can penetrate a continuous
bed between them. Dense sampling, geometry discovery, finite-area pressures,
friction, adhesion and impact need additional models.

`LayeredStrip::bending_response` condenses zero axial force from the same
history-aware implicit section matrix used by `advance_to`. It returns the
free curvature and positive moment compliance at the supplied temperature and
duration, without evolving history or subjecting a fictitious unit-load state
to the material limits. The beam's exact discrete curvature basis constructs
a Gram contact flexibility. The existing `fs_solver::orthant` owner solves
nonnegative reactions; no penalty stiffness or new optimization algorithm is
introduced. Duplicate IDs, invalid stations, exhausted work, and invalid final
strain/geometry refuse. Dependent contact rows do not acquire a uniqueness claim.

The support configuration bounds up to 64 stops, 4096 curvature cells,
stop-squared-times-cell setup products and coordinate sweeps. Its gradient
tolerance is in metres. Independent gap and reaction-times-gap gates are checked
on the final actual beam, not only on the condensed system. All physical point
updates publish together after the existing force/work/storage checks. Rejected
solves and cancellation leave the original beam untouched.

Reports keep built-in reactions, each stop force/gap, material energy accounts
and endpoint reaction work. A moved or newly engaged stop can do nonzero endpoint
quadrature work; this is not exact time-integrated impact work or a passivity
guarantee. Stop sets and bounds can change at explicit step boundaries. Clone
and replay with the same physical inputs preserve the accepted material history.
Thermal area, conductance and bond strength are not inferred from a point force.

## Coupled deposition example

```sh
cargo test -p fs-couple --example sagging_deposition
cargo run -p fs-couple --example sagging_deposition
```

The example advances the actual thermal solver and two separate beam histories:
freely sagging and supported by a thermally insulating midpoint roller. Both
see the same original cooling inputs. A later externally supported straight road
then contacts their different solved top surfaces. Captured reference area and
an explicitly supplied conductance determine subsequent reheating, compared
also with a flat-road counterfactual. The roller does not create a finite heat
contact, and the later road is not assumed perfectly bonded to the beam.
All coefficients are synthetic, not measured printer/material data.

## Verification and remaining limits

Seven beam regressions cover G1 elastic spatial refinement, creep time refinement
and a near-Maxwell viscous limit, G0/G3 load work and thermal bending, and G4/G5
birth/refusal/cancellation/replay. Four geometry regressions exercise capture
area, interior gaps, contraction and penetration. Spatial, material/time and
coupling resolution are separate.
Seven additional support regressions cover propped-beam reaction refinement,
support removal and lift-off, moved/multiple stops, creep engagement, replay,
budget/cancellation refusal, and agreement of the condensed section response
with the real material-history update. One coupled-example regression checks
the resulting contact-area and temperature changes, including the support case.

Native Rust compilation/tests have not run in the authoring environment:
DSR, its documented absolute fallback, RCH and Cargo were unavailable. Independent
Python/SciPy references test equations, not the Rust compiler or actual runtime.
A ten-to-one slenderness screen and declared strain/rotation/deflection gates
bound admission; they are not material calibration or discretization error bounds.
No finite rotations, fracture, dynamic settling, arbitrary body-to-body contact,
whole-print geometry discovery or qualified printed-part prediction is claimed.
