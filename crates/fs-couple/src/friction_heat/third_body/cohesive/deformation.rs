//! Heat through a surviving adhesive layer under bounded affine deformation.

use super::{
    CohesiveHeatError, CohesiveLayerHeatBinding, CohesiveLayerHeatControls,
    CohesiveLayerHeatSource, HeatedCohesiveLayerCandidate,
    propose_deformable_cohesive_layer_heat_inner,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsStep};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Explicit strain domain for a fixed-area, laterally uniform adhesive layer.
/// Normal jump is extension of the original layer, not a vacuum crack opening.
/// Tangential jump is bounded simple shear; edge overlap and lateral heat flux
/// are not resolved. The caller must choose bounds within the original joint's
/// small-deformation and material calibration domains.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveLayerDeformation {
    maximum_normal_strain: f64,
    maximum_shear_strain: f64,
}

impl CohesiveLayerDeformation {
    /// Maximum absolute normal and Euclidean shear jump divided by reference
    /// thickness. Both are finite and nonnegative; normal strain is below one
    /// so every admitted compression retains a positive physical thickness.
    /// A zero bound explicitly holds that channel fixed.
    pub fn new(
        maximum_normal_strain: f64,
        maximum_shear_strain: f64,
    ) -> Result<Self, DepositionError> {
        if !maximum_normal_strain.is_finite()
            || !(0.0..1.0).contains(&maximum_normal_strain)
            || !maximum_shear_strain.is_finite()
            || maximum_shear_strain < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "invalid adhesive layer strain domain",
            ));
        }
        Ok(Self {
            maximum_normal_strain,
            maximum_shear_strain,
        })
    }

    pub const fn maximum_normal_strain(&self) -> f64 {
        self.maximum_normal_strain
    }
    pub const fn maximum_shear_strain(&self) -> f64 {
        self.maximum_shear_strain
    }

    pub(super) fn normal_stretch(
        self,
        jump_m: [f64; 3],
        failed: bool,
        reference_thickness_m: f64,
    ) -> Result<f64, &'static str> {
        if failed {
            return Err("a ruptured adhesive needs an explicit crack/contact thermal transition");
        }
        if !reference_thickness_m.is_finite()
            || reference_thickness_m <= 0.0
            || jump_m.iter().any(|v| !v.is_finite())
        {
            return Err("invalid adhesive layer deformation");
        }
        let normal = jump_m[0] / reference_thickness_m;
        let shear = det::hypot(jump_m[1], jump_m[2]) / reference_thickness_m;
        if !normal.is_finite()
            || normal.abs() > self.maximum_normal_strain
            || !shear.is_finite()
            || shear > self.maximum_shear_strain
        {
            return Err("adhesive layer strain exceeds declared domain");
        }
        Ok(1.0 + normal)
    }
}

/// Advance native tensile/shear adhesion and its original finite thermal layer.
///
/// At accepted endpoint opening g, thickness is L=L0+g and the two thermal
/// paths are A/[d_body/k_body + L/(2 k_layer) + R_area]. Conductivities remain
/// the physical spatial source values; only layer distance changes. Equivalently
/// K0_nn=k/lambda follows J F^-1 k F^-T for fixed-area normal stretch lambda.
/// Reference mass and m*h storage remain unchanged, so no heat capacity or mass
/// is created by stretching. The report retains `heat.layer_normal_stretch`.
///
/// Both endpoint strains and surviving attachment are required. The original
/// straight constitutive path and irreversible damage then retain these bounds
/// throughout the step. Partial damage keeps its original nonthermal fracture
/// account; it is not inferred to be missing geometric area. The callback may
/// query explicitly calibrated damage-dependent conductivity/resistance. Complete
/// rupture refuses the whole proposal before any thermal source is evaluated.
///
/// The original source identity, native work, heat-location and energy gates
/// remain in force. A changing exposed side area is not modeled: the layer's
/// ambient heat-transfer coefficient must be zero. Explicit total source and
/// reservoir powers retain their own meaning. Mechanical coefficients are frozen
/// over the step and thermal transport uses endpoint backward Euler; resolve
/// load and temperature histories by timestep refinement.
#[allow(clippy::too_many_arguments)]
pub fn propose_deformable_cohesive_layer_heat_with_deformation<'s, 'a, 't, P, E>(
    motion: &CohesiveMotion<'s, 'a>,
    mean_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    dt_s: f64,
    binding: CohesiveLayerHeatBinding<'_>,
    deformation: CohesiveLayerDeformation,
    thermal: &DepositionSimulation<'t>,
    controls: CohesiveLayerHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        f64,
        &DynamicsStep<'s, 'a>,
    ) -> Result<CohesiveLayerHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedCohesiveLayerCandidate<'s, 'a, 't, P>, CohesiveHeatError<E>> {
    propose_deformable_cohesive_layer_heat_inner(
        motion,
        mean_force_n,
        endpoint_prescribed_jumps_m,
        dt_s,
        binding,
        Some(deformation),
        thermal,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )
}
