# Spatial material and thermal rod equilibrium

`fs_solid::rod::material_field::SpatialMaterialRod` connects source-resolved
concentric sections to the existing Cosserat rod statics. A shaft, member or
slender mechanism part may have different material layers in each segment.
Spatial order matters: moving a compliant section toward a clamp changes
bending even when the material inventory and arithmetic-average modulus stay
the same. Segments share nodal positions and frames; they are not independent
rod solves patched together afterwards.

## Material, thermal state and supports are separate inputs

Build each `ConcentricElasticSection` from its ordered material layers, then
call `with_shear_modulus(property_name, shear_correction)` to obtain the
`ConcentricShearSection`. The rod consumes the original density/E/G claims and
section integrals. It does not infer shear modulus from a material name or
silently impose a Poisson ratio. Static support must be inside any declared
frequency band. Current-temperature curves are evaluated at the same source
query coordinates used by their section.

For heating, integrate each layer's expansion with
`integrate_isotropic_thermal_expansion` and construct a
`ConcentricThermoelasticSection` from the SAME elastic section. Supply one
thermal section per spatial segment. This integrates the actual alpha(T) path,
not alpha at the endpoint times the total temperature difference.

```rust
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::supports::RodSupport;

// `mechanical_sections` and `thermal_sections` are source-resolved vectors
// in spatial order, with identical per-segment elastic bundles and geometry.
let gate = CancelGate::new();
let rod = SpatialMaterialRod::try_straight(length_m, mechanical_sections)?;
let heated = rod.with_thermal_sections(thermal_sections, 0.02, &gate)?;
let mut supported = heated.with_supports(&[RodSupport {
    node: segment_count,
    position_m: [Some(length_m), Some(0.0), Some(0.0)],
    frame: Some(So3::identity()),
}], &gate)?;
let load = TipLoad::default();
let residuals = supported.solve_static(&load, 1, 1e-7, &gate)?;
let reactions = supported.support_reactions(&load, &gate)?;
let layer_response = supported.thermal_axial_response(0)?;
```

Node zero remains clamped. Additional supports may fix any world-position
components and/or the complete body-to-world frame. Omitting a component leaves
it free. An empty support list restores the root-only boundary condition.
Constraints project the existing residual and tangent, not penalty springs.
`support_reactions` reports forces on the rod in world coordinates and moments
in the local body frame, after subtracting the supplied applied load. It also
reports the residual on free coordinates. It excludes the root reaction and is
not a whole-body balance certificate. Fully constrained members need no Newton
solve but can still have stresses, reactions and stored energy.

`positions_m`, `frames`, `segment_strains`, `segment_resultants` and `energy_j`
observe the same accepted configuration. Full thermal layer responses retain
self-equilibrated stresses: zero axial resultant does not imply zero layer
stress or zero stored energy. Failed source admission, an invalid imposed pose,
solver failure or cancellation does not alter the original specimen.

## Domain and verification

This is a bounded static rod model: 1..=64 equal-reference-length segments,
linear diagonal axial/shear/bending/torsional section laws, existing SO(3)
kinematics, dense finite-difference Newton, and the original body-frame tip
moment convention. Material/cross-section jumps idealize perfectly compatible
joints; 3D interface stresses, helical winding, slip, damage, dynamics and
stability qualification are not supplied. Thermal strain is axial and section
geometry is frozen. Circumferentially uniform concentric layers generate no
thermal bending moment in this model.

Thermal states are prescribed, not solved here. Binding a different state or
support returns an independent static specimen; it does not conserve the work,
heat or history of a dynamic transition. The strain limit constrains each
layer's thermal and axial elastic strain, not every bending-fiber strain or
yield envelope. FD probes must also lie in that domain. The stopping norm mixes
force and moment components, as in the original solver; it is not a physical
error bound or experimental certificate.

Targeted Rust tests are authored in `rod_material_field` and `rod_thermal_field`:

```sh
rch exec -- cargo test -p fs-solid --test rod_material_field --test rod_thermal_field
rch exec -- cargo test -p fs-solid --lib rod::tests
rch exec -- cargo test -p fs-material --test concentric_thermal
```

Native Rust execution was unavailable in the implementation environment.
Independent Python/SciPy equations checked series compliance, spatial bending,
heated extension, opposing layer stress, torsion and prescribed-support
reactions; these do not replace compilation or native regression execution.
