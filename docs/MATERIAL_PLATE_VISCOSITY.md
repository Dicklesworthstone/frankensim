# Material tensors in shared string-body dynamics

`fs_material::plane_stress_viscosity::PlaneStressViscosity` describes an explicitly
characterized **reduced plane-stress** Kelvin–Voigt tensor. This is not a scalar
loss factor or a three-dimensional viscosity condensed independently of elasticity.
The six original scalar properties are in Pa s, in material-frame order
11, 22, 66, 12, 16, 26. All six must exist, including explicitly zero couplings:

```text
plane_stress_viscosity_11
plane_stress_viscosity_22
plane_stress_viscosity_66
plane_stress_viscosity_12
plane_stress_viscosity_16
plane_stress_viscosity_26
```

The strain-rate vector is `[epsilon_dot_11, epsilon_dot_22, gamma_dot_12]` and
the stress vector is `[sigma_11, sigma_22, sigma_12]`. Their dot product is physical
dissipation density. V must be resolved positive definite or explicitly all zero;
nonzero singular tensors are currently refused, not repaired. Coefficients must
be constant at the selected state (temperature-only curves are allowed), with
explicit overlapping frequency support. No material name supplies missing data.

## Connect the existing regional body

The `SharedMaterialPlateChartBody` mesh, original material bundles, orientation,
thickness/mass constraints, support, bridge footprint, cavity and modal budgets
remain the inputs. Every region's material bundle must include the six selected
viscosity properties, and `body.reduction.damping_ratio` must be zero.

```rust
use fs_couple::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
use fs_couple::bowed_string::runtime::reciprocal::bank::ReciprocalStringBank;
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;

// body, string_sources, gas, cavity, friction, options, cx and gate are the
// existing physical assembly inputs. No notes or target frequencies are added.
let viscosity = PlateViscosity::from_chart(
    &body.specimen, PLANE_STRESS_VISCOSITY_PROPERTIES, &gate,
)?;
let limits = PlateViscosityLimits {
    max_modes: 16,
    max_projection_terms: 10_000_000,
    reference_resistance_n_s_m: 1.0,
};
let bank = ReciprocalStringBank::from_viscous_material_plate_chart(
    &cx, string_sources, gas, body, viscosity, limits, cavity,
    friction, options, &gate,
)?;
```

`bank.source().material_chart_viscosity()` retains the original tensors, material
states and full mass-normalized damping matrix. The same sources persist through
single- and multiple-bow scheduling. A cold material replacement uses the existing
`with_region_material` and then `viscosity.rebind(&new_chart, &gate)`: all coefficients
are re-resolved. Passing the old loss binding with a changed chart refuses. No
running mechanical history is reset or reinterpreted as a new material.

## One basis, full damping

The original DKT/modal reduction supplies force shapes, radiating areas, actual
modal masses and all normalized eigenvectors in **one** solve. Each element uses
`D_viscous = V_rotated h^3/12`; neither elastic stiffness nor geometric prestress
is automatically damped. Projection retains every entry of `C = Phi^T C_FE Phi`,
including signed off-diagonal terms. Arbitrary orientation and Orthotropic21's
ordered-axis exchange affect elasticity and viscosity consistently.

The existing linear-algebra owner factors C. Its columns become pure dashpots
in `CoupledModalSystem`, so the existing simultaneous bow/bridge/cavity solve sees
their reactions and accounts for their loss. The resistance value above is an
explicit factor-coordinate scale: changing it inversely rescales columns and
preserves C, but the original per-channel force-admission ceiling uses that scale.
Zero tensors add no connections and retain the legacy zero-loss path. Nonzero
rank-unresolved modal matrices refuse rather than losing channels.

Generic plate consumers can call `certified_viscous_radiators` or
`certified_viscous_radiators_with_force_ports`; returned radiators have **no scalar
material damping**. Connect `image.damping().connections(body_component)` exactly
once, to those same returned modes, rather than attaching it to another eigenbasis.

## Limits and checks

This is linear flat-plate bending in a truncated real basis, with full viscous
coupling inside that basis. It does not model carved shells, frequency-dependent
memory, thermal evolution, radiation back-loading or measured violin/cello/bass
fidelity. Spatial/modal and time refinement remain necessary, especially for
strong damping and unresolved modes. The original coupled step is not the exact
continuous damped-system exponential. Material/source bounds do not certify
nonlinear bow harmonics or complete instrument realism.

The existing native test seam selects the focused integration checks:

```bash
cargo test -p fs-material --test plane_stress_viscosity
cargo test -p fs-couple --lib thin_plate::material_loss::ports::tensor_tests -- --test-threads=1
```

The latter covers full DKT cross-mode reactions and heat, common force/eigenbasis,
zero-source parity, viscosity-only source feedback and pressure, scheduler callback
invariance, retained source identity, cold rebinding and numerical/physical refusal.
These tests were authored but not executed natively in the editing environment.
Independent matrix-equation checks are not a substitute for that native run.
