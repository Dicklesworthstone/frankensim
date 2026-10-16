# Bounded thermal-history-dependent deposited strips

Issue #1 now has a uniaxial thermoviscoelastic constitutive law and an axially
uniform, perfectly bonded section consumer. This is bounded R2 functionality,
not a general 3D print, adhesion, fracture or sag predictor.

`fs_material::thermovisco::ThermoMaxwell` accepts the reference modulus and
branch `(E, tau)` pairs used by `visco::GeneralizedMaxwell`, plus explicit
bounded temperature tables for common modulus/relaxation multipliers and free
thermal strain. Points start stress-free at their actual birth temperature
and total strain. Accepted internal variables are immutable inputs to a trial;
changing temperature never resets them. Existing fixed-temperature Prony APIs
are unchanged. This new update is backward Euler with an exact affine tangent,
not the older exact-exponential update. Do not equate a glass transition with
a latent heat plateau or infer crystallization kinetics from these tables.

`fs_solid::deposition::LayeredStrip` uses `epsilon(z)=epsilon0-z*kappa` and two
Gauss points per rectangular layer. Each point owns independent branch history.
It solves axial force/bending moment or prescribed strain/curvature, reports
reactions and independently checks section equilibrium and work conjugacy.
Layer birth uses the current strain and curvature, leaving old states untouched.
The deformed local position is derived from that same solved mechanical state.
Reference length, layer dimensions, loading, thermal histories, material ranges,
small-strain/rotation limits and work/state budgets are explicit SI inputs.
Thermal contact is not evidence of mechanical bonding: perfect bonding is a
separate caller assumption in this reduced consumer.

## Working example

```sh
cargo test -p fs-material thermovisco
cargo test -p fs-solid deposition
cargo test -p fs-solid --example deposited_thermovisco_strip
cargo run -p fs-solid --example deposited_thermovisco_strip -- 0.25
```

The example uses the actual std-only `fs-conduction` deposition producer by
source path, not a copied thermal solver. This development-only arrangement
avoids a new production dependency between physics crates. It activates three
passes at 0, 10 and 20 s; resolves reheating and constrained cooling; prescribes
bed release at 30 s and a two-second force/moment unloading ramp; then cools to
60 s. An interval publishes neither thermal nor mechanical state until both
candidates succeed. A layer born at the endpoint receives no preceding aging.
The example retains reference geometry, source IDs, forces, curvature, thermal
energy and before/after-release positions. The coefficients are synthetic,
not measured filament or interface data.

The endpoint work audit separates mechanical work, thermal-eigenstrain work,
modulus-change work, storage, physical dashpot dissipation and backward-Euler
numerical dissipation. The global section also measures work from the actual
resultants and strain/curvature increments. Numerical dissipation is not
physical viscous heating. A sudden release can retain algorithmic dissipation
even for a tiny dt, so resolve a loading ramp when release work matters.

## Validation scope

The `Deposition physics` workflow compiles the exact changed production modules
and coupled example in debug and optimized modes without substituting numerical
kernels. This avoids unrelated missing sibling crates, but does not substitute
for full workspace validation or the repository's pinned-toolchain gates.
Seven constitutive, six section and two coupled regressions cover independent
analytic limits, history retention, energy accounts, tangent differentiation,
refinement, birth, cancellation and regime refusal. Check each commit's native
job results; authored tests alone do not establish execution success.

Thermal subcycling, coupling cadence and section subdivision are separate
resolution choices. The example's cadence ladder refines thermal and mechanical
steps together; the section test checks isothermal subdivision invariance, not
3D spatial convergence. Energy residuals are not physical prediction accuracy.

## Remaining limitations

The consumer is axially uniform and uniaxial; layers span the complete reference
length. Slenderness, geometry and nonoverlapping volumes are caller obligations.
No shear/torsion, Poisson contraction, general material-frame anisotropy, inertia,
gravity, calibrated healing, cohesive damage, adhesion-based detachment time,
finite strain, geometry/contact feedback, or viscous-heating feedback is supplied.
Temperatures outside the supplied table, excessive strain/rate/rotation and
unrepresentable arithmetic refuse rather than silently escalating validity.
A prescribed release does not predict adhesion failure. Clone is an in-memory
checkpoint, not a persistent restart format.
