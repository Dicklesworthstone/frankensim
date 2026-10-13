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
There is no automatic force/contact iteration, nozzle impact prediction or
strength mapping; use an appropriate coupled/contact owner for those claims.

```sh
cargo test -p fs-solid bridge
```

Thermal consumers may bind captured reference areas to their matching cell IDs
and an independently declared conductance per original area. Do not add a later
road to the perfectly bonded beam merely because it falls inside the geometric
capture distance. Contact-force transfer and strength need their own laws.

## Verification and remaining limits

Seven beam regressions cover G1 elastic spatial refinement, creep time refinement
and a near-Maxwell viscous limit, G0/G3 load work and thermal bending, and G4/G5
birth/refusal/cancellation/replay. Four geometry regressions exercise capture
area, interior gaps, contraction and penetration. Spatial, material/time and
coupling resolution are separate.

Native Rust compilation/tests have not run in the authoring environment:
DSR, its documented absolute fallback, RCH and Cargo were unavailable. Independent
Python/SciPy references test equations, not the Rust compiler or actual runtime.
A ten-to-one slenderness screen and declared strain/rotation/deflection gates
bound admission; they are not material calibration or discretization error bounds.
No finite rotations, fracture, feedback contact forces, dynamic settling,
whole-print geometry discovery or qualified printed-part prediction is claimed.
