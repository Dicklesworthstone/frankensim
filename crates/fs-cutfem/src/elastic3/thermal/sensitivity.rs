//! Temperature transpose of the original, fixed-geometry mechanical residual.
//!
//! Temperature changes the full constitutive tensor as well as integrated free
//! strain. The same change also reaches the symmetric Nitsche terms and both
//! sides of every ghost penalty. No density interpolation, finite difference,
//! constitutive fit or positive-definite surrogate for C' is used here.

use std::ops::ControlFlow;

use fs_material::elastic::thermal::ThermoelasticTemperatureDerivatives;
use fs_material::tensor::{contract, from_mandel, to_mandel};
use fs_solver::op::LinearOp;

use crate::elastic3::eigenstrain::{material_stress, physical_response, project_free};
use crate::elastic3::stress::mandel_strain;
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3, q1};

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

/// A signed derivative action, never an admitted elastic material.
#[derive(Clone, Copy)]
enum TemperatureTangent3<'a> {
    Tensor(&'a [[f64; 6]; 6]),
    Isotropic { bulk: f64, shear: f64 },
}

impl<'a> TemperatureTangent3<'a> {
    fn new(
        cell: &Cell3,
        derivative: &'a ThermoelasticTemperatureDerivatives,
    ) -> Result<Self, ElasticityError3> {
        if cell.kinematics.is_projected() {
            let (bulk, shear) =
                derivative
                    .isotropic_moduli_pa_k
                    .ok_or(ElasticityError3::Invalid(
                        "B-bar temperature derivative requires native K/G slopes",
                    ))?;
            Ok(Self::Isotropic { bulk, shear })
        } else {
            Ok(Self::Tensor(&derivative.stiffness_mandel_pa_k))
        }
    }

    /// Mandel strain to physical stress, following the native B-bar mode split.
    fn stress(self, strain: [f64; 6]) -> [f64; 6] {
        match self {
            Self::Tensor(matrix) => from_mandel(&std::array::from_fn(|i| {
                matrix[i].iter().zip(strain).map(|(c, e)| c * e).sum()
            })),
            Self::Isotropic { bulk, shear } => from_mandel(
                &fs_material::elastic::isotropic_stress_from_moduli(bulk, 2.0 * shear, &strain),
            ),
        }
    }
}

/// Differentiate max_i sum_j |C_ij| only on its differentiable branches.
/// An inactive row has a strict gap, so its zero entries do not affect the local
/// maximum. An active absolute-value kink or unequal active-row slopes does.
fn stiffness_bound_derivative(
    matrix: &[[f64; 6]; 6],
    derivative: &[[f64; 6]; 6],
) -> Result<f64, ElasticityError3> {
    let sums: [f64; 6] = std::array::from_fn(|i| matrix[i].iter().map(|v| v.abs()).sum());
    let maximum = sums.into_iter().fold(0.0_f64, f64::max);
    let mut active_slope = None;
    for i in 0..6 {
        if sums[i] != maximum {
            continue;
        }
        let mut terms = [0.0; 6];
        for j in 0..6 {
            if matrix[i][j] == 0.0 {
                if derivative[i][j] != 0.0 {
                    return Err(ElasticityError3::Invalid(
                        "temperature derivative of stiffness bound is nonsmooth at a zero entry",
                    ));
                }
            } else {
                terms[j] = matrix[i][j].signum() * derivative[i][j];
            }
        }
        // Permuting an isotropic row cannot change its exact derivative.
        // Canonical summation avoids manufacturing unequal tied slopes from
        // a different order of the same six floating-point contributions.
        terms.sort_by(f64::total_cmp);
        let slope: f64 = terms.into_iter().sum();
        if !slope.is_finite() {
            return Err(ElasticityError3::Invalid(
                "stiffness-bound temperature derivative overflow",
            ));
        }
        if active_slope.is_some_and(|previous| previous != slope) {
            return Err(ElasticityError3::Invalid(
                "tied stiffness bounds have different temperature derivatives",
            ));
        }
        active_slope = Some(slope);
    }
    active_slope.ok_or(ElasticityError3::Invalid(
        "invalid retained stiffness bound",
    ))
}

fn traction(stress: [f64; 6], normal: [f64; 3]) -> [f64; 3] {
    [
        stress[0] * normal[0] + stress[3] * normal[1] + stress[5] * normal[2],
        stress[3] * normal[0] + stress[1] * normal[1] + stress[4] * normal[2],
        stress[5] * normal[0] + stress[4] * normal[1] + stress[2] * normal[2],
    ]
}

/// Negative temperature derivative of the symmetric Nitsche motion load.
fn motion_temperature_work(
    motion: [f64; 3],
    motion_slope: [f64; 3],
    traction_slope: [f64; 3],
    source_traction: [f64; 3],
    dual_value: [f64; 3],
    penalty_slope: f64,
    penalty: f64,
) -> f64 {
    (0..3)
        .map(|c| {
            (traction_slope[c] - penalty_slope * dual_value[c]) * motion[c]
                + (source_traction[c] - penalty * dual_value[c]) * motion_slope[c]
        })
        .sum()
}

impl CutElasticity3 {
    pub(in crate::elastic3) fn admit_thermal_temperature_points(
        &self,
        max_points: usize,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(&mut checkpoint)?;
        let mut points = self.ghosts.len();
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            points = points
                .checked_add(cell.rules.bulk().len())
                .and_then(|n| {
                    n.checked_add(
                        self.embedded
                            .as_ref()
                            .map_or(0, |data| data.selected_points(id).len()),
                    )
                })
                .ok_or(ElasticityError3::Invalid(
                    "thermal temperature point count overflow",
                ))?;
            if points > max_points {
                return Err(ElasticityError3::Invalid(
                    "thermal temperature point allowance exhausted",
                ));
            }
            if cell.material.thermal_state().is_none() {
                return Err(ElasticityError3::Invalid(
                    "temperature pullback requires a native thermoelastic source in every cell",
                ));
            }
        }
        poll(&mut checkpoint)
    }

    /// Exact partial-temperature transpose of `R = K(T)u - b_free(T) - f`.
    ///
    /// Returns `dual^T dR/dT_cell` in original active-cell order. Temperatures
    /// are independent, cell-constant prescribed source states. Geometry,
    /// original reference density, stiffness scales and external force are
    /// held fixed. The caller supplies any explicit objective derivative and
    /// the negative sign required by its equilibrium-adjoint convention.
    /// This is not a heat-equation adjoint or an adiabatic tangent.
    ///
    /// All active cells must retain their native thermoelastic source. The
    /// integrated free strain and exact tensor/modulus slopes come from that
    /// same source; an unrelated prescribed eigenstrain cannot be substituted.
    /// Original actual and B-bar physical strains are admitted at retained bulk
    /// and selected Nitsche points, against the source and the separate total
    /// strain budget in `(0,1)`. The dual is algebraic and has no strain limit.
    ///
    /// Homogeneous strong and embedded supports are supported. Use
    /// [`Self::thermal_residual_temperature_pullback_with_motion`] for nonzero
    /// embedded motion and its optional temperature dependence. Strong rows have zero
    /// temperature derivative. The full tensor C', native B-bar K'/G', symmetric
    /// Nitsche consistency/penalty, and both ghost-side coefficients are kept.
    /// A relevant nonsmooth absolute-value or maximum-row branch in the native
    /// stabilization modulus is refused rather than assigned a gradient.
    ///
    /// `max_points` bounds bulk + selected boundary + ghost quadrature points
    /// before allocation. Scratch/output storage is O(active cells); there is
    /// no global displacement expansion, matrix assembly or source callback.
    /// Cancellation or failure publishes no partial result or operator change.
    pub fn thermal_residual_temperature_pullback(
        &self,
        u: &[f64],
        dual: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.thermal_residual_temperature_pullback_mapped(
            u,
            dual,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            None,
            None,
            |node, component, field| {
                if self.fixed[node] {
                    0.0
                } else {
                    field[3 * node + component]
                }
            },
            checkpoint,
        )
    }

    /// Partial-temperature transpose including the original embedded support lift.
    ///
    /// Returns `dual^T d(K(T)u-b_free(T)-b_g(T)-f)/dT_cell` for the
    /// displacement field supplied to `prescribed_displacement_load`.
    /// `prescribed` is the actual motion in metres at the retained support
    /// points. Its optional derivative is metres per kelvin with respect to
    /// that point's owning cell temperature; omission means fixed motion.
    /// Both callbacks must be pure and describe the same thermal experiment.
    ///
    /// This differentiates both constitutive consistency and penalty lifting,
    /// including their product with temperature-dependent boundary motion.
    /// Strong box supports remain homogeneous, so motion must be compatible
    /// with them. Geometry, reference inertia, scales and external forces stay
    /// fixed. Physical admission, nonsmooth stabilization refusal and bounded
    /// cancellation are inherited from
    /// [Self::thermal_residual_temperature_pullback].
    ///
    /// The original bulk/surface/ghost point budget is checked before either
    /// motion callback. No additional quadrature, global lifting vector or
    /// derivative matrix is allocated. An operator without an embedded support
    /// refuses this entry point.
    #[allow(clippy::too_many_arguments)]
    pub fn thermal_residual_temperature_pullback_with_motion(
        &self,
        u: &[f64],
        dual: &[f64],
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        prescribed_temperature_derivative: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        maximum_kinematic_strain: f64,
        max_points: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        self.thermal_residual_temperature_pullback_mapped(
            u,
            dual,
            self.n(),
            maximum_kinematic_strain,
            max_points,
            Some(prescribed),
            prescribed_temperature_derivative,
            |node, component, field| {
                if self.fixed[node] {
                    0.0
                } else {
                    field[3 * node + component]
                }
            },
            checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::elastic3) fn thermal_residual_temperature_pullback_mapped(
        &self,
        u: &[f64],
        dual: &[f64],
        n: usize,
        maximum_kinematic_strain: f64,
        max_points: usize,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        prescribed_temperature_derivative: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        coefficient: impl Fn(usize, usize, &[f64]) -> f64,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        if prescribed.is_some() && self.embedded.is_none() {
            return Err(ElasticityError3::Invalid(
                "thermal motion pullback requires an embedded Dirichlet support",
            ));
        }
        if u.len() != n || dual.len() != n {
            return Err(ElasticityError3::Invalid(
                "thermal temperature pullback field shape mismatch",
            ));
        }
        if !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(ElasticityError3::Invalid(
                "invalid total-strain temperature pullback budget",
            ));
        }
        for field in [u, dual] {
            for chunk in field.chunks(192) {
                poll(&mut checkpoint)?;
                if chunk.iter().any(|v| !v.is_finite()) {
                    return Err(ElasticityError3::Invalid(
                        "nonfinite thermal temperature pullback field",
                    ));
                }
            }
        }
        self.admit_thermal_temperature_points(max_points, &mut checkpoint)?;
        let mut ghost_needed = vec![false; self.cells.len()];
        for face in &self.ghosts {
            poll(&mut checkpoint)?;
            for id in face.cells {
                ghost_needed[id] = true;
            }
        }
        let mut ghost_slopes = vec![0.0; self.cells.len()];
        let mut result = vec![0.0; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let state = cell.material.thermal_state().expect("source preflight");
            let derivative = state.source().temperature_derivatives().map_err(|_| {
                ElasticityError3::Invalid("native thermoelastic temperature derivative unavailable")
            })?;
            let tangent = TemperatureTangent3::new(cell, &derivative)?;
            let local_u: [f64; 24] =
                std::array::from_fn(|i| coefficient(cell.nodes[i / 3], i % 3, u));
            let local_z: [f64; 24] =
                std::array::from_fn(|i| coefficient(cell.nodes[i / 3], i % 3, dual));
            if local_u.iter().chain(&local_z).any(|v| !v.is_finite()) {
                return Err(ElasticityError3::Invalid(
                    "thermal reconstructed field overflow",
                ));
            }
            let free_mandel = state.source().free_strain_mandel();
            let free = from_mandel(&free_mandel);
            let free_trace = free_mandel[..3].iter().sum::<f64>();
            let alpha = derivative.free_strain_mandel_per_k;
            let alpha = project_free(cell, alpha, alpha[..3].iter().sum());
            let expansion_stress = material_stress(cell, alpha);
            let selected = self
                .embedded
                .as_ref()
                .map_or(&[][..], |data| data.selected_points(id));
            let bound_slope =
                if !selected.is_empty() || (ghost_needed[id] && !cell.kinematics.is_projected()) {
                    stiffness_bound_derivative(
                        state.acoustic().stiffness_mandel_pa(),
                        &derivative.stiffness_mandel_pa_k,
                    )?
                } else {
                    0.0
                };
            if ghost_needed[id] {
                ghost_slopes[id] = match tangent {
                    TemperatureTangent3::Isotropic { shear, .. } => shear,
                    TemperatureTangent3::Tensor(_) => bound_slope,
                };
            }
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let physical = physical_response(
                    cell,
                    &local_u,
                    position,
                    free,
                    free_trace,
                    self.scales[id],
                    maximum_kinematic_strain,
                )?;
                let (_, gradients) = q1(cell.bounds, position);
                let strain_z =
                    cell.project_mandel_strain(mandel_strain(&gradients, &local_z), &local_z);
                let stiffness_stress =
                    tangent.stress(to_mandel(&physical.constitutive_elastic_strain));
                let residual_stress =
                    std::array::from_fn(|i| stiffness_stress[i] - expansion_stress[i]);
                result[id] +=
                    self.scales[id] * weight * contract(&from_mandel(&strain_z), &residual_stress);
            }
            if !selected.is_empty() {
                let data = self.embedded.as_ref().expect("selected support");
                let h = (0..3)
                    .map(|a| cell.bounds.hi()[a] - cell.bounds.lo()[a])
                    .fold(f64::INFINITY, f64::min);
                let gamma_slope = data.beta() * (bound_slope / h);
                let surface = cell.rules.surface().expect("retained support surface");
                for &index in selected {
                    poll(&mut checkpoint)?;
                    let point = &surface.points()[index];
                    let physical = physical_response(
                        cell,
                        &local_u,
                        point.position,
                        free,
                        free_trace,
                        self.scales[id],
                        maximum_kinematic_strain,
                    )?;
                    let (basis, gradients) = q1(cell.bounds, point.position);
                    let strain_z =
                        cell.project_mandel_strain(mandel_strain(&gradients, &local_z), &local_z);
                    let stiffness_stress =
                        tangent.stress(to_mandel(&physical.constitutive_elastic_strain));
                    let residual_traction = traction(
                        std::array::from_fn(|i| stiffness_stress[i] - expansion_stress[i]),
                        point.normal,
                    );
                    let dual_traction = traction(tangent.stress(strain_z), point.normal);
                    let value_u: [f64; 3] = std::array::from_fn(|c| {
                        (0..8).map(|a| basis[a] * local_u[3 * a + c]).sum()
                    });
                    let value_z: [f64; 3] = std::array::from_fn(|c| {
                        (0..8).map(|a| basis[a] * local_z[3 * a + c]).sum()
                    });
                    // b_free has NEGATIVE eigenstress traction: the first
                    // consistency term differentiates the elastic stress.
                    let mut work: f64 = (0..3)
                        .map(|c| {
                            -value_z[c] * residual_traction[c] - dual_traction[c] * value_u[c]
                                + gamma_slope * value_z[c] * value_u[c]
                        })
                        .sum();
                    if let Some(prescribed) = prescribed {
                        let motion = prescribed(point.position, point.normal);
                        poll(&mut checkpoint)?;
                        let motion_slope =
                            if let Some(derivative) = prescribed_temperature_derivative {
                                let slope = derivative(point.position, point.normal);
                                poll(&mut checkpoint)?;
                                slope
                            } else {
                                [0.0; 3]
                            };
                        if motion
                            .iter()
                            .chain(&motion_slope)
                            .any(|value| !value.is_finite())
                        {
                            return Err(ElasticityError3::Invalid(
                                "nonfinite prescribed motion or temperature derivative",
                            ));
                        }
                        // -d b_g/dT: the consistency load has the opposite sign
                        // from the positive penalty load. Include the product
                        // rule when the authored boundary itself moves with T.
                        let (source_traction, gamma) =
                            if prescribed_temperature_derivative.is_some() {
                                (
                                    traction(material_stress(cell, strain_z), point.normal),
                                    data.beta() * (cell.material.nitsche_modulus() / h),
                                )
                            } else {
                                ([0.0; 3], 0.0)
                            };
                        work += motion_temperature_work(
                            motion,
                            motion_slope,
                            dual_traction,
                            source_traction,
                            value_z,
                            gamma_slope,
                            gamma,
                        );
                    }
                    result[id] += self.scales[id] * point.weight * work;
                }
            }
            if !result[id].is_finite() {
                return Err(ElasticityError3::Invalid(
                    "thermal temperature pullback overflow",
                ));
            }
        }
        for face in &self.ghosts {
            poll(&mut checkpoint)?;
            let mut work = 0.0;
            for component in 0..3 {
                let mut jump_u = 0.0;
                let mut jump_z = 0.0;
                for (&node, &jump) in face.nodes.iter().zip(&face.jump) {
                    jump_u += jump * coefficient(node, component, u);
                    jump_z += jump * coefficient(node, component, dual);
                }
                work += jump_u * jump_z;
            }
            for id in face.cells {
                let slope = face.geometry.weight(ghost_slopes[id]);
                result[id] += (0.5 * slope * self.scales[id]) * work;
                if !result[id].is_finite() {
                    return Err(ElasticityError3::Invalid(
                        "thermal ghost temperature pullback overflow",
                    ));
                }
            }
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::stiffness_bound_derivative;

    #[test]
    fn moving_nitsche_load_temperature_work_matches_independent_differences() {
        let motion = [0.003, -0.002, 0.001];
        let traction = [7.0, -4.0, 2.0];
        let traction_slope = [-0.3, 0.4, 0.2];
        let dual_value = [1.5, -0.8, 0.7];
        let (penalty, penalty_slope) = (90.0, -1.2);
        for motion_slope in [[0.0; 3], [0.00004, -0.00002, 0.00003]] {
            let load_work = |delta: f64| -> f64 {
                (0..3)
                    .map(|c| {
                        (-(traction[c] + delta * traction_slope[c])
                            + (penalty + delta * penalty_slope) * dual_value[c])
                            * (motion[c] + delta * motion_slope[c])
                    })
                    .sum()
            };
            let h = 0.001;
            let numerical = -(load_work(h) - load_work(-h)) / (2.0 * h);
            let actual = super::motion_temperature_work(
                motion,
                motion_slope,
                traction_slope,
                traction,
                dual_value,
                penalty_slope,
                penalty,
            );
            assert!((actual - numerical).abs() < 2e-10 * actual.abs());
        }
    }

    #[test]
    fn bound_derivative_admits_only_smooth_active_rows() {
        let mut matrix = [[0.0; 6]; 6];
        let mut derivative = [[0.0; 6]; 6];
        matrix[0][0] = 3.0;
        matrix[1][1] = 2.0;
        derivative[0][0] = -0.4;
        derivative[1][2] = 0.7; // A strictly inactive kink is immaterial.
        assert_eq!(
            stiffness_bound_derivative(&matrix, &derivative).unwrap(),
            -0.4,
        );
        derivative[0][2] = 0.2;
        assert!(stiffness_bound_derivative(&matrix, &derivative).is_err());
        derivative[0][2] = 0.0;
        derivative[1][2] = 0.0;
        matrix[1][1] = 3.0;
        derivative[1][1] = 0.1;
        assert!(stiffness_bound_derivative(&matrix, &derivative).is_err());
        derivative[1][1] = -0.4;
        assert_eq!(
            stiffness_bound_derivative(&matrix, &derivative).unwrap(),
            -0.4,
        );
    }
}
