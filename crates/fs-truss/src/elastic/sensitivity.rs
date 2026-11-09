//! Implicit design derivatives of compatible truss displacement functionals.
//!
//! For J = sum(g dot u), the adjoint solves K_ff lambda_f = g_f with homogeneous
//! support motion through the same native matrix-free owner as the primal.
//! No derivative is taken through CG iterations. Member derivatives include
//! both stiffness and the imposed-state RHS, using the actual elastic extension.
//!
//! Geometry, support flags, objective weights and source material states are
//! fixed. Positive retained areas may vary smoothly; removed members have no
//! activation derivative. This does not differentiate catalog snapping, material
//! state resolution, buckling, geometry, or a self-weight fixed-point iteration.

use super::{
    ElasticControls, ElasticError, ElasticImposedState, ElasticResponse, checkpoint, invalid,
};
use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;

/// Partial derivatives for one retained original member at fixed geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(clippy::struct_field_names)] // Each field is a named physical-input derivative.
pub struct MemberSensitivity {
    /// dJ/dA with E, stress-free strain and physical nodal loads fixed [J-units/m²].
    pub area_derivative: f64,
    /// dJ/dE with A and the original material state otherwise frozen [J-units/Pa].
    pub young_modulus_derivative: f64,
    /// dJ/d(epsilon_0), with original length fixed [J-units].
    pub stress_free_strain_derivative: f64,
}

/// One admitted primal/adjoint pair and its original-order design derivatives.
///
/// "J-units" means the units of the caller's weighted displacement functional.
/// Choosing physical loads for the weights gives J in joules (load compliance).
#[derive(Debug, Clone)]
pub struct ElasticSensitivity {
    /// Actual selected-area response, including imposed strain and support motion.
    pub response: ElasticResponse,
    /// sum(g dot u), including weighted supported displacements.
    pub value: f64,
    /// Original member order; None means removed (zero area), not zero derivative
    /// for adding that member back. Retained-member derivatives preserve topology.
    pub members: Vec<Option<MemberSensitivity>>,
    /// dJ/df at fixed objective weights, in original node/component order [J-units/N].
    ///
    /// If compliance weights track the loads, its total load derivative is
    /// lambda + u, not this fixed-weight lambda alone. It reduces to 2u only
    /// for zero eigenstrain and homogeneous supports.
    /// Supported components are zero. On a mechanism-containing structure, only
    /// load variations that retain solvable equilibrium are admissible; this
    /// covector is not a claim that every independent load can be applied.
    pub load_gradient: Vec<[f64; 2]>,
    /// dJ/du_p = g_p - K_pf lambda_f on supported coordinates; zero on free ones.
    /// Units are the functional's units per metre.
    pub prescribed_displacement_gradient: Vec<[f64; 2]>,
    /// Recomputed physical equilibrium residual of the adjoint solve.
    pub adjoint_relative_equilibrium_residual: f64,
    /// Adjoint native CG iterations. Together with response.iterations this
    /// never exceeds the supplied maximum_iterations.
    pub adjoint_iterations: usize,
}

/// Differentiate J = sum(displacement_weights dot u) by one native adjoint.
///
/// The input weights have the complete original node ordering. For compliance,
/// pass case.loads(); those weights and the actual external loads remain fixed
/// during each reported partial derivative. For a retained member, writing
/// e_el = B^T u - L epsilon_0 and z = B_f^T lambda_f gives
/// dJ/dA = -(E/L) e_el z, dJ/dE = -(A/L) e_el z, and dJ/d(epsilon_0) = EA z.
/// The elastic extension includes imposed strain and support motion: replacing
/// it by total extension would give incorrect compliance gradients.
///
/// The iteration budget is shared across the primal and adjoint solves. Each
/// uses the existing checked-CG acceptance and bounded defect correction; no
/// partial derivative result is published if either physical equilibrium fails.
///
/// These are frozen-geometry, frozen-source partial derivatives. If a parameter
/// also changes external loads, add load_gradient dot df/dtheta. If it changes
/// the objective weights, also add u dot dg/dtheta. Thus load-dependent
/// compliance requires both terms, including when loads come from self-weight.
/// Source-temperature, density, geometry and catalog derivatives are not inferred.
///
/// # Errors
/// Malformed/nonfinite weights, primal or adjoint refusal, unrepresentable
/// derivatives, exhausted shared work and cancellation return no accepted result.
/// An objective that loads an unsupported free mechanism refuses even if the
/// original physical loads alone admitted a compatible response.
#[allow(clippy::too_many_arguments)] // One physical solve plus its objective covector.
pub fn displacement_functional(
    ground: &GroundStructure,
    case: &LayoutCase,
    areas_m2: &[f64],
    young_pa: &[f64],
    imposed: ElasticImposedState<'_>,
    displacement_weights: &[[f64; 2]],
    controls: ElasticControls,
    cx: &Cx<'_>,
) -> Result<ElasticSensitivity, ElasticError> {
    checkpoint(cx)?;
    if displacement_weights.len() != ground.nodes().len() {
        return Err(TrussConstructionError::VectorLength {
            field: "elastic objective weights",
            expected: ground.nodes().len(),
            actual: displacement_weights.len(),
        }
        .into());
    }
    for (node, weights) in displacement_weights.iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        if weights.iter().any(|value| !value.is_finite()) {
            return Err(invalid("elastic objective weights", "must be finite").into());
        }
    }
    let response =
        super::solve_with_imposed(ground, case, areas_m2, young_pa, imposed, controls, cx)?;
    let adjoint_case = LayoutCase::try_new(
        case.supported().to_vec(),
        displacement_weights.to_vec(),
        ground.nodes().len(),
    )?;
    let adjoint = super::solve(
        ground,
        &adjoint_case,
        areas_m2,
        young_pa,
        ElasticControls {
            maximum_iterations: controls
                .maximum_iterations
                .saturating_sub(response.iterations),
            ..controls
        },
        cx,
    )?;
    let mut members = Vec::with_capacity(ground.members().len());
    for (index, &area) in areas_m2.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if area == 0.0 {
            members.push(None);
            continue;
        }
        let length = ground.lengths()[index];
        let adjoint_extension = length * adjoint.member_strain[index];
        let derivative = MemberSensitivity {
            area_derivative: -response.member_stress_pa[index] * adjoint_extension,
            young_modulus_derivative: -(area * response.member_elastic_strain[index])
                * adjoint_extension,
            stress_free_strain_derivative: adjoint.member_force_n[index] * length,
        };
        if [
            derivative.area_derivative,
            derivative.young_modulus_derivative,
            derivative.stress_free_strain_derivative,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(invalid("elastic design derivatives", "must remain finite").into());
        }
        members.push(Some(derivative));
    }
    let mut value = 0.0;
    let mut support_gradient = adjoint.support_reactions_n;
    for (node, gradient) in support_gradient.iter_mut().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (axis, entry) in gradient.iter_mut().enumerate() {
            value += displacement_weights[node][axis] * response.displacements_m[node][axis];
            *entry = if case.supported()[node][axis] {
                -*entry
            } else {
                0.0
            };
        }
        if !value.is_finite() {
            return Err(invalid("elastic objective value", "must remain finite").into());
        }
    }
    checkpoint(cx)?;
    Ok(ElasticSensitivity {
        response,
        value,
        members,
        load_gradient: adjoint.displacements_m,
        prescribed_displacement_gradient: support_gradient,
        adjoint_relative_equilibrium_residual: adjoint.relative_equilibrium_residual,
        adjoint_iterations: adjoint.iterations,
    })
}
