//! Physical fields and motion on the original Tet10/P1 mixed quadrature.
//!
//! The supplied pressure is interpolated in the original continuous P1 space;
//! the full first Piola stress comes from the SAME mixed material state as force
//! assembly. Exact incompressibility is weak: neither J nor Jdot is projected.
//! No pressure rate, constitutive history evolution or caloric heat is inferred.

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, MixedPointState};

use super::{
    MixedUpError, MixedUpProblem, PreparedMixedUpProblem, checkpoint, invalid, point_deformation,
};
use crate::hyper3::response::{motion_fields, stress_strain_fields};

/// One original mixed material point, in element-major Duffy quadrature order.
///
/// All tensors are row-major with physical tensor shear. Stress includes the
/// original pressure multiplier and reference/history adapter. Rates describe
/// total geometry, not an inferred elastic/plastic/thermal decomposition.
#[derive(Clone, Debug)]
pub struct MixedPointMotionResponse {
    pub element: usize,
    pub local_point: usize,
    /// Position of this quadrature point in the original reference mesh [m].
    pub reference_position_m: [f64; 3],
    /// Original quadrature volume, not the whole element volume [m³].
    pub reference_volume_m3: f64,
    pub current_volume_m3: f64,
    pub deformation_gradient: [f64; 9],
    /// Total Green--Lagrange strain (F^T F-I)/2.
    pub green_lagrange_strain: [f64; 9],
    /// Interpolated independent pressure [Pa], positive in tension for the
    /// native mixed potential. It need not equal mean Cauchy stress.
    pub pressure_pa: f64,
    /// Original full Piola stress, weak pressure residual and mixed potential.
    /// The saddle potential is not physical stored energy. A frozen inelastic
    /// remainder may be incremental; this observer advances no history.
    pub material_state: MixedPointState,
    /// Full physical P F^T/J [Pa], without stress symmetrization.
    pub cauchy_pa: [f64; 9],
    /// Actual compatible det(F), equal to material_state.volume_ratio.
    pub jacobian: f64,
    pub deformation_gradient_rate_per_s: [f64; 9],
    pub spatial_velocity_gradient_per_s: [f64; 9],
    /// D=sym(Fdot F^-1), transforming as Q D Q^T under moving spatial frames.
    pub spatial_strain_rate_per_s: [f64; 9],
    /// Edot=sym(F^T Fdot), invariant under moving spatial frames.
    pub green_lagrange_strain_rate_per_s: [f64; 9],
    /// Compatible Jdot=J tr(L), not a projected incompressibility constraint.
    pub jacobian_rate_per_s: f64,
    pub current_volume_rate_m3_s: f64,
    /// Quadrature weight times full physical P:Fdot [W].
    /// This is stress power, not the derivative of the saddle potential or heat.
    pub reference_stress_power_w: f64,
    /// Independently reconstructed current-volume sigma:L [W].
    pub spatial_stress_power_w: f64,
}

impl<L: MixedFiniteStrainPotential> MixedUpProblem<'_, L> {
    /// Observe native mixed fields and motion without solving or advancing history.
    ///
    /// Displacement and velocity are original-node xyz vectors in m and m/s.
    /// Pressure follows the sorted corner-node order of pressure_nodes() on the
    /// prepared problem. Supplied fields need not be an equilibrium solution.
    /// Preparation preserves the original geometry, material and budget gates.
    pub fn quadrature_motion_responses(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<MixedPointMotionResponse>, MixedUpError> {
        self.prepare(cx)?
            .quadrature_motion_responses(displacement_m, velocity_m_s, pressure_pa, cx)
    }
}

impl<L: MixedFiniteStrainPotential> PreparedMixedUpProblem<'_, L> {
    /// Reuse the original prepared Tet10 geometry and original point materials.
    ///
    /// Only the native mixed stress/state seam is queried, never its tangent or
    /// an extracted elastic backbone. Each point retains its own weak pressure
    /// residual and compatible volume rate, including for exact incompressibility.
    /// Cancellation, malformed fields, material-domain exits or unrepresentable
    /// rates/power return no partial response. Supports and external loads do not
    /// change this local observation, and force equilibrium is not asserted.
    pub fn quadrature_motion_responses(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<MixedPointMotionResponse>, MixedUpError> {
        checkpoint(cx)?;
        self.validate_state(displacement_m, pressure_pa)?;
        if velocity_m_s.len() != self.mesh.nodes_count {
            return Err(invalid("one xyz velocity per original mixed node required"));
        }
        for (node, velocity) in velocity_m_s.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            if velocity.iter().any(|v| !v.is_finite()) {
                return Err(invalid("mixed nodal velocity must be finite"));
            }
        }
        let mut responses =
            Vec::with_capacity(self.mesh.elements.len() * self.mesh.points_per_element);
        for (element_index, element) in self.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let refused = |what: &str| MixedUpError::Material {
                    element: element_index,
                    point: q,
                    error: MaterialError::State {
                        what: format!("mixed motion recovery: {what}"),
                    },
                };
                let f = point_deformation(element, point, displacement_m);
                let mut rate = [0.0; 9];
                for a in 1..10 {
                    for i in 0..3 {
                        let relative = velocity_m_s[element.connectivity[a]][i]
                            - velocity_m_s[element.connectivity[0]][i];
                        for j in 0..3 {
                            rate[3 * i + j] =
                                relative.mul_add(point.gradients[a][j], rate[3 * i + j]);
                        }
                    }
                }
                let p = (0..4)
                    .map(|a| pressure_pa[element.pressure_indices[a]] * point.pressure_shapes[a])
                    .sum();
                let material_state = self
                    .problem
                    .material
                    .at(element_index, q, self.mesh.points_per_element)
                    .mixed_state(&f, p, self.problem.volumetric)
                    .map_err(|error| MixedUpError::Material {
                        element: element_index,
                        point: q,
                        error,
                    })?;
                let jacobian = material_state.volume_ratio;
                let current_volume_m3 = point.weight_m3 * jacobian;
                let reference_position_m = core::array::from_fn(|i| {
                    (0..10).fold(0.0, |sum, a| {
                        point.shapes[a]
                            .mul_add(self.problem.nodes_m[element.connectivity[a]][i], sum)
                    })
                });
                let (cauchy_pa, green_lagrange_strain) =
                    stress_strain_fields(&f, &material_state.first_piola_pa, jacobian);
                if !(jacobian.is_finite()
                    && jacobian > 0.0
                    && current_volume_m3.is_finite()
                    && current_volume_m3 > 0.0)
                    || [
                        p,
                        material_state.potential_pa,
                        material_state.remainder_potential_pa,
                        material_state.bulk_stored_pa,
                        material_state.bulk_volume_ratio,
                        material_state.bulk_complementary_pa,
                        material_state.pressure_residual,
                    ]
                    .iter()
                    .any(|v| !v.is_finite())
                    || f.iter()
                        .chain(&rate)
                        .chain(&reference_position_m)
                        .chain(&material_state.first_piola_pa)
                        .chain(&cauchy_pa)
                        .chain(&green_lagrange_strain)
                        .any(|v| !v.is_finite())
                {
                    return Err(refused(
                        "stress, geometry or material state is unrepresentable",
                    ));
                }
                let fields = motion_fields(
                    f,
                    rate,
                    jacobian,
                    point.weight_m3,
                    current_volume_m3,
                    &material_state.first_piola_pa,
                    &cauchy_pa,
                )
                .map_err(refused)?;
                responses.push(MixedPointMotionResponse {
                    element: element_index,
                    local_point: q,
                    reference_position_m,
                    reference_volume_m3: point.weight_m3,
                    current_volume_m3,
                    deformation_gradient: f,
                    green_lagrange_strain,
                    pressure_pa: p,
                    material_state,
                    cauchy_pa,
                    jacobian,
                    deformation_gradient_rate_per_s: rate,
                    spatial_velocity_gradient_per_s: fields.spatial_velocity_gradient_per_s,
                    spatial_strain_rate_per_s: fields.spatial_strain_rate_per_s,
                    green_lagrange_strain_rate_per_s: fields.green_lagrange_strain_rate_per_s,
                    jacobian_rate_per_s: fields.jacobian_rate_per_s,
                    current_volume_rate_m3_s: fields.current_volume_rate_m3_s,
                    reference_stress_power_w: fields.reference_stress_power_w,
                    spatial_stress_power_w: fields.spatial_stress_power_w,
                });
            }
        }
        checkpoint(cx)?;
        Ok(responses)
    }
}
