//! Frozen mixed temperature/mechanical residual block on the native cut rules.
//!
//! Preparation integrates direct local temperature-force columns, not dual
//! basis probes. Repeated forward and transpose actions share those columns,
//! the original ghost stencils and the original displacement constraint map.

use std::ops::ControlFlow;
use std::sync::Arc;

use fs_material::tensor::{from_mandel, to_mandel};

use super::sensitivity::{
    TemperatureTangent3, motion_temperature_work, stiffness_bound_derivative, traction,
};
use crate::elastic3::eigenstrain::{material_stress, physical_response, project_free};
use crate::elastic3::material::engineering_strain;
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3, q1};

/// Explicit preparation and reuse allowances; no unbounded default is supplied.
#[derive(Debug, Clone, Copy)]
pub struct ThermalJacobianLimits3 {
    /// Bulk, selected Nitsche and ghost points, counted before callbacks.
    pub max_points: usize,
    /// Peak cached/scratch scalar slots, each charged as eight bytes, excluding
    /// constant stack scratch and Vec/owner metadata. Boolean scratch is charged
    /// conservatively as a full slot. No displacement-sized scratch is needed.
    pub max_workspace_components: usize,
    /// Bound on scalar stencil work per action, including the original hanging
    /// rows. Linear input/output checks are additionally bounded by their sizes.
    pub max_action_terms: usize,
}

struct CellColumn3 {
    force: [f64; 24],
    scale_bits: u64,
}

struct GhostColumn3 {
    jump: [f64; 3],
    /// Already include each side's one-half, geometry and stiffness scale.
    weights: [f64; 2],
}

/// Exact partial `d(K(T)u-b_free(T)-b_g(T)-f)/dT` at a frozen native state.
///
/// Temperatures are independent per retained cell; geometry, reference mass,
/// scales, external forces and physical displacement are fixed. The optional
/// embedded support motion and its local temperature slope are sampled once.
/// Strong box supports remain homogeneous. This is a partial residual block,
/// not a solved-equilibrium derivative or a heat-equation Jacobian.
///
/// Local columns include both symmetric Nitsche terms and penalty derivatives;
/// ghost columns include both material sides. Forward/transpose actions use the
/// SAME cached coefficients. No positive material surrogate for signed C' or
/// finite-difference production kernel is introduced. A changed physical state
/// or boundary law needs new preparation even if the operator is unchanged.
/// Owner and exact scale checks reject stale geometry/material/space/scales.
pub struct PreparedThermalResidualJacobian3 {
    binding: Arc<()>,
    adaptive: bool,
    n: usize,
    cells: Vec<CellColumn3>,
    ghosts: Vec<GhostColumn3>,
    points: usize,
    workspace_components: usize,
    action_terms: usize,
}

impl core::fmt::Debug for PreparedThermalResidualJacobian3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedThermalResidualJacobian3")
            .field("displacement_dofs", &self.n)
            .field("temperatures", &self.cells.len())
            .field("points", &self.points)
            .field("workspace_components", &self.workspace_components)
            .field("action_terms", &self.action_terms)
            .finish_non_exhaustive()
    }
}

/// A borrowed view of the existing native map, never a second constraint owner.
#[derive(Clone, Copy)]
pub(in crate::elastic3) struct ThermalCoordinateMap3<'a> {
    pub(in crate::elastic3) rows: Option<&'a [Vec<(usize, f64)>]>,
    pub(in crate::elastic3) fixed: &'a [bool],
}

impl ThermalCoordinateMap3<'_> {
    pub(super) fn n(self) -> usize {
        self.fixed.len() * 3
    }

    pub(super) fn value(self, field: &[f64], node: usize, component: usize) -> f64 {
        if let Some(rows) = self.rows {
            rows[node]
                .iter()
                .map(|&(master, weight)| {
                    if self.fixed[master] {
                        0.0
                    } else {
                        weight * field[3 * master + component]
                    }
                })
                .sum()
        } else if self.fixed[node] {
            0.0
        } else {
            field[3 * node + component]
        }
    }

    pub(super) fn scatter(self, node: usize, component: usize, value: f64, output: &mut [f64]) {
        if let Some(rows) = self.rows {
            for &(master, weight) in &rows[node] {
                if !self.fixed[master] {
                    output[3 * master + component] += weight * value;
                }
            }
        } else if !self.fixed[node] {
            output[3 * node + component] += value;
        }
    }

    pub(super) fn terms(self, node: usize) -> usize {
        self.rows.map_or(1, |rows| rows[node].len())
    }
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(
    field: &[f64],
    checkpoint: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    for chunk in field.chunks(192) {
        poll(checkpoint)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite thermal Jacobian field or action"));
        }
    }
    Ok(())
}

// Direct derivative of the native u_a-u_0 strain reconstruction. Node zero is
// retained as the exact negative sum, including the original B-bar projection.
pub(super) fn strain_basis(cell: &Cell3, gradients: &[[f64; 3]; 8]) -> [[f64; 6]; 24] {
    let mut basis = [[0.0; 6]; 24];
    for a in 1..8 {
        for component in 0..3 {
            let mut strain = engineering_strain(gradients[a], component);
            for value in &mut strain[3..] {
                *value *= core::f64::consts::FRAC_1_SQRT_2;
            }
            if let Some((_, _, mean)) = cell.kinematics.projected_moduli() {
                let correction = (mean[a][component] - gradients[a][component]) / 3.0;
                for value in &mut strain[..3] {
                    *value += correction;
                }
            }
            basis[3 * a + component] = strain;
        }
    }
    for component in 0..3 {
        for k in 0..6 {
            basis[component][k] = -(1..8).map(|a| basis[3 * a + component][k]).sum::<f64>();
        }
    }
    basis
}

impl PreparedThermalResidualJacobian3 {
    /// Original independent displacement dimension, including zero clamp rows.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Number of independent native cell temperatures in original cell order.
    #[must_use]
    pub fn temperatures(&self) -> usize {
        self.cells.len()
    }

    /// Bulk, selected support and ghost points consumed during preparation.
    #[must_use]
    pub const fn points(&self) -> usize {
        self.points
    }

    /// Conservatively charged peak scalar storage during preparation.
    #[must_use]
    pub const fn workspace_components(&self) -> usize {
        self.workspace_components
    }

    /// Admitted bound on the original native stencil work per reusable action.
    #[must_use]
    pub const fn action_terms(&self) -> usize {
        self.action_terms
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::elastic3) fn prepare(
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermalJacobianLimits3,
        prescribed: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        motion_derivative: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        if u.len() != map.n()
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid("invalid thermal Jacobian displacement or strain budget"));
        }
        if prescribed.is_some() && source.embedded.is_none() {
            return Err(invalid("thermal motion Jacobian requires embedded support"));
        }
        finite(u, &mut checkpoint)?;
        let workspace_components = source.cells.len().checked_mul(27)
            .and_then(|n| source.ghosts.len().checked_mul(5).and_then(|g| n.checked_add(g)))
            .ok_or_else(|| invalid("thermal Jacobian workspace count overflow"))?;
        if workspace_components > limits.max_workspace_components {
            return Err(invalid("thermal Jacobian workspace allowance exhausted"));
        }
        let mut points = source.ghosts.len();
        let mut action_terms = source.ghosts.len().checked_mul(10)
            .ok_or_else(|| invalid("thermal Jacobian stencil count overflow"))?;
        for (id, cell) in source.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            if cell.material.thermal_state().is_none() {
                return Err(invalid("thermal Jacobian requires a native thermal source in every cell"));
            }
            points = points.checked_add(cell.rules.bulk().len())
                .and_then(|n| n.checked_add(source.embedded.as_ref().map_or(0, |data| data.selected_points(id).len())))
                .ok_or_else(|| invalid("thermal Jacobian point count overflow"))?;
            if points > limits.max_points {
                return Err(invalid("thermal Jacobian point allowance exhausted"));
            }
            for &node in &cell.nodes {
                action_terms = map.terms(node).checked_mul(6)
                    .and_then(|n| action_terms.checked_add(n))
                    .ok_or_else(|| invalid("thermal Jacobian stencil count overflow"))?;
                if action_terms > limits.max_action_terms {
                    return Err(invalid("thermal Jacobian action allowance exhausted"));
                }
            }
        }
        for face in &source.ghosts {
            poll(&mut checkpoint)?;
            for &node in &face.nodes {
                action_terms = map.terms(node).checked_mul(6)
                    .and_then(|n| action_terms.checked_add(n))
                    .ok_or_else(|| invalid("thermal Jacobian stencil count overflow"))?;
                if action_terms > limits.max_action_terms {
                    return Err(invalid("thermal Jacobian action allowance exhausted"));
                }
            }
        }
        if points > limits.max_points || action_terms > limits.max_action_terms {
            return Err(invalid("thermal Jacobian point or action allowance exhausted"));
        }

        // Every size/work check above precedes allocation and boundary callbacks.
        let mut ghost_needed = vec![false; source.cells.len()];
        for face in &source.ghosts {
            poll(&mut checkpoint)?;
            for id in face.cells {
                ghost_needed[id] = true;
            }
        }
        let mut ghost_slopes = vec![0.0; source.cells.len()];
        let mut cells = Vec::with_capacity(source.cells.len());
        for (id, cell) in source.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let state = cell.material.thermal_state().expect("thermal source preflight");
            let derivative = state.source().temperature_derivatives()
                .map_err(|_| invalid("native thermoelastic temperature derivative unavailable"))?;
            let tangent = TemperatureTangent3::new(cell, &derivative)?;
            let local_u: [f64; 24] = std::array::from_fn(|i| map.value(u, cell.nodes[i / 3], i % 3));
            finite(&local_u, &mut checkpoint)?;
            let free_mandel = state.source().free_strain_mandel();
            let free = from_mandel(&free_mandel);
            let free_trace = free_mandel[..3].iter().sum::<f64>();
            let alpha = derivative.free_strain_mandel_per_k;
            let expansion_stress = material_stress(cell, project_free(cell, alpha, alpha[..3].iter().sum()));
            let selected = source.embedded.as_ref().map_or(&[][..], |data| data.selected_points(id));
            let bound_slope = if !selected.is_empty() || (ghost_needed[id] && !cell.kinematics.is_projected()) {
                stiffness_bound_derivative(state.acoustic().stiffness_mandel_pa(), &derivative.stiffness_mandel_pa_k)?
            } else {
                0.0
            };
            ghost_slopes[id] = tangent.ghost_slope(cell, bound_slope);
            let scale = source.scales[id];
            let mut force = [0.0; 24];
            for &(position, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let physical = physical_response(cell, &local_u, position, free, free_trace, scale, maximum_kinematic_strain)?;
                let stiffness_stress = tangent.stress(to_mandel(&physical.constitutive_elastic_strain));
                let residual_stress = to_mandel(&std::array::from_fn(|i| stiffness_stress[i] - expansion_stress[i]));
                let (_, gradients) = q1(cell.bounds, position);
                let basis = strain_basis(cell, &gradients);
                for i in 3..24 {
                    force[i] += (scale * weight) * basis[i].iter().zip(residual_stress).map(|(b, stress)| b * stress).sum::<f64>();
                }
            }
            for component in 0..3 {
                force[component] = -(1..8).map(|a| force[3 * a + component]).sum::<f64>();
            }
            if !selected.is_empty() {
                let data = source.embedded.as_ref().expect("selected support");
                let h = (0..3).map(|a| cell.bounds.hi()[a] - cell.bounds.lo()[a]).fold(f64::INFINITY, f64::min);
                let gamma_slope = data.beta() * (bound_slope / h);
                let gamma = data.beta() * (cell.material.nitsche_modulus() / h);
                let surface = cell.rules.surface().expect("retained support surface");
                for &index in selected {
                    poll(&mut checkpoint)?;
                    let point = &surface.points()[index];
                    let physical = physical_response(cell, &local_u, point.position, free, free_trace, scale, maximum_kinematic_strain)?;
                    let stiffness_stress = tangent.stress(to_mandel(&physical.constitutive_elastic_strain));
                    let residual_traction = traction(std::array::from_fn(|i| stiffness_stress[i] - expansion_stress[i]), point.normal);
                    let (values, gradients) = q1(cell.bounds, point.position);
                    let basis = strain_basis(cell, &gradients);
                    let value_u: [f64; 3] = std::array::from_fn(|c| (0..8).map(|a| values[a] * local_u[3 * a + c]).sum());
                    let motion = if let Some(callback) = prescribed {
                        let value = callback(point.position, point.normal);
                        poll(&mut checkpoint)?;
                        value
                    } else {
                        [0.0; 3]
                    };
                    let motion_slope = if let Some(callback) = motion_derivative {
                        let value = callback(point.position, point.normal);
                        poll(&mut checkpoint)?;
                        value
                    } else {
                        [0.0; 3]
                    };
                    finite(&motion, &mut checkpoint)?;
                    finite(&motion_slope, &mut checkpoint)?;
                    for i in 0..24 {
                        let c = i % 3;
                        let value = values[i / 3];
                        let traction_slope = traction(tangent.stress(basis[i]), point.normal);
                        let mut work = -value * residual_traction[c]
                            - traction_slope.iter().zip(value_u).map(|(t, v)| t * v).sum::<f64>()
                            + gamma_slope * value * value_u[c];
                        if prescribed.is_some() {
                            let source_traction = if motion_derivative.is_some() {
                                traction(material_stress(cell, basis[i]), point.normal)
                            } else {
                                [0.0; 3]
                            };
                            work += motion_temperature_work(motion, motion_slope, traction_slope, source_traction,
                                std::array::from_fn(|component| if component == c { value } else { 0.0 }), gamma_slope, gamma);
                        }
                        force[i] += (scale * point.weight) * work;
                    }
                }
            }
            finite(&force, &mut checkpoint)?;
            cells.push(CellColumn3 { force, scale_bits: scale.to_bits() });
        }
        let mut ghosts = Vec::with_capacity(source.ghosts.len());
        for face in &source.ghosts {
            poll(&mut checkpoint)?;
            let jump: [f64; 3] = std::array::from_fn(|component| face.nodes.iter().zip(&face.jump)
                .map(|(&node, &coefficient)| coefficient * map.value(u, node, component)).sum());
            let weights = face.cells.map(|id| 0.5 * face.geometry.weight(ghost_slopes[id]) * source.scales[id]);
            finite(&jump, &mut checkpoint)?;
            finite(&weights, &mut checkpoint)?;
            ghosts.push(GhostColumn3 { jump, weights });
        }
        poll(&mut checkpoint)?;
        Ok(Self {
            binding: Arc::clone(&source.material_body_binding),
            adaptive: map.rows.is_some(),
            n: map.n(), cells, ghosts, points, workspace_components, action_terms,
        })
    }

    fn validate_for(
        &self,
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if !Arc::ptr_eq(&self.binding, &source.material_body_binding)
            || self.adaptive != map.rows.is_some() || self.n != map.n()
        {
            return Err(invalid("thermal Jacobian belongs to another operator or space"));
        }
        for (id, (cell, scale)) in self.cells.iter().zip(&source.scales).enumerate() {
            if id % 64 == 0 { poll(checkpoint)?; }
            if cell.scale_bits != scale.to_bits() {
                return Err(invalid("thermal Jacobian stiffness scales changed"));
            }
        }
        poll(checkpoint)
    }

    pub(in crate::elastic3) fn apply_mapped(
        &self,
        source: &CutElasticity3,
        map: ThermalCoordinateMap3<'_>,
        input: &[f64],
        output: &mut [f64],
        transpose: bool,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.validate_for(source, map, &mut checkpoint)?;
        let (rows, columns) = if transpose { (self.cells.len(), self.n) } else { (self.n, self.cells.len()) };
        if input.len() != columns || output.len() != rows {
            return Err(invalid("thermal Jacobian action shape mismatch"));
        }
        finite(input, &mut checkpoint)?;
        poll(&mut checkpoint)?;
        for chunk in output.chunks_mut(192) {
            poll(&mut checkpoint)?;
            chunk.fill(0.0);
        }
        for (id, (cell, column)) in source.cells.iter().zip(&self.cells).enumerate() {
            poll(&mut checkpoint)?;
            for i in 0..24 {
                let node = cell.nodes[i / 3];
                if transpose {
                    output[id] += column.force[i] * map.value(input, node, i % 3);
                } else {
                    map.scatter(node, i % 3, column.force[i] * input[id], output);
                }
            }
        }
        for (face, column) in source.ghosts.iter().zip(&self.ghosts) {
            poll(&mut checkpoint)?;
            if transpose {
                let mut work = 0.0;
                for component in 0..3 {
                    let jump: f64 = face.nodes.iter().zip(&face.jump)
                        .map(|(&node, &coefficient)| coefficient * map.value(input, node, component)).sum();
                    work += column.jump[component] * jump;
                }
                for (id, weight) in face.cells.into_iter().zip(column.weights) {
                    output[id] += weight * work;
                }
            } else {
                let weight: f64 = face.cells.into_iter().zip(column.weights).map(|(id, weight)| weight * input[id]).sum();
                let weighted_jump = column.jump.map(|value| weight * value);
                for (&node, &jump) in face.nodes.iter().zip(&face.jump) {
                    for component in 0..3 {
                        map.scatter(node, component, jump * weighted_jump[component], output);
                    }
                }
            }
        }
        finite(output, &mut checkpoint)?;
        poll(&mut checkpoint)
    }
}

impl CutElasticity3 {
    /// Freeze the native temperature residual block at this physical u.
    /// Domain, derivative-branch, point, storage and stencil-work admission occur
    /// once. Repeated actions integrate no quadrature and query no sources.
    pub fn prepare_thermal_residual_jacobian(
        &self,
        u: &[f64],
        maximum_kinematic_strain: f64,
        limits: ThermalJacobianLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermalResidualJacobian3, ElasticityError3> {
        PreparedThermalResidualJacobian3::prepare(self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, None, None, checkpoint)
    }

    /// Include nonzero embedded support motion and its optional LOCAL owning-cell
    /// temperature derivative, in metres per kelvin. Both pure callbacks must
    /// describe the same motion as `prescribed_displacement_load`. Point/work/
    /// storage preflight precedes callbacks; each callback is followed by a poll.
    /// Geometry, external forcing and strong homogeneous box clamps stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_thermal_residual_jacobian_with_motion(
        &self,
        u: &[f64],
        prescribed: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        prescribed_temperature_derivative: Option<&dyn Fn([f64; 3], [f64; 3]) -> [f64; 3]>,
        maximum_kinematic_strain: f64,
        limits: ThermalJacobianLimits3,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedThermalResidualJacobian3, ElasticityError3> {
        PreparedThermalResidualJacobian3::prepare(self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            u, maximum_kinematic_strain, limits, Some(prescribed), prescribed_temperature_derivative, checkpoint)
    }

    /// Allocation-free `J_T * delta_temperature` in native displacement rows.
    /// No physical strain/temperature-domain bound applies to algebraic directions.
    /// Ownership, exact scales, dimensions and finite input are checked before
    /// writing. Mid-action cancellation/overflow may leave partial output; discard
    /// it on error. All strong-support output rows are zero.
    pub fn thermal_residual_temperature_apply_into(
        &self,
        prepared: &PreparedThermalResidualJacobian3,
        temperature_direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            temperature_direction, output, false, checkpoint)
    }

    /// Allocation-free `J_T^T * dual`, reusing the exact forward coefficients.
    /// The finite dual is algebraic and has no physical strain ceiling. Refusal
    /// and partial-output semantics match `thermal_residual_temperature_apply_into`.
    pub fn thermal_residual_temperature_transpose_into(
        &self,
        prepared: &PreparedThermalResidualJacobian3,
        dual: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        prepared.apply_mapped(self, ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            dual, output, true, checkpoint)
    }
}
