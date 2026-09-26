//! Two-field goal comparison using the same residual and tangent as the
//! production steady solve. Nonlinear remainder is observed, not discarded.

use super::{
    ConductionError, ConductionProblem, Cx, DofMap, LinearConfig, RobinResponse, ThermalInterfaces,
    add, admit_linear, assemble_operator_scaled_with_interfaces, checked, failed, invalid,
    nonlinear, poll, reduce, temperature_dependent, true_residual, vector,
};

mod algebraic;
pub use algebraic::{
    LinearGoalAnalysis, LinearGoalAnalysisConfig, LinearGoalAnalyzer, LinearGoalSolve,
    LinearGoalSolveConfig, LinearGoalStop, LinearMaximumAnalysis, LinearMaximumSolve,
    analyze_linear_goal, analyze_linear_maximum,
};

mod feedback;
pub use feedback::{RobinGoalFeedback, RobinGoalLinearization};

/// An observed linear-goal difference on one discrete mesh. An enriched-mesh
/// caller may supply a prolonged coarse field as the approximation. This is
/// neither a continuum bound nor a material/shape derivative.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscreteGoalComparison {
    /// Sum of lambda_i [R(T_reference) - R(T_approximate)]_i, in goal units.
    pub signed_residual_change: f64,
    /// Observed goal change minus the dual-weighted residual change. Includes
    /// the finite nonlinear linearization remainder and remaining dual-solve
    /// or floating-point error; it must not be omitted from an error estimate.
    pub linearization_remainder: f64,
    /// Directly measured weights . (T_reference - T_approximate).
    pub signed_goal_change: f64,
    /// Signed residual contributions in FULL nodal order, zero on fixed nodes.
    /// These prioritize refinement; they are not local continuum error bounds.
    pub nodal_contributions: Vec<f64>,
    /// Independently checked residual of the supplied reference primal.
    pub primal_relative_residual: f64,
    /// Independently checked residual of the actual transposed tangent solve.
    pub dual_relative_residual: f64,
    /// Inner Krylov iterations consumed by this comparison, excluding primals.
    pub dual_iterations: usize,
    /// True when the actual nonsymmetric K'(T) material tangent was used.
    pub uses_nonlinear_jacobian: bool,
}

/// Compare a linear temperature goal against a residual-checked reference.
///
/// With R(T) = A(T) T - b and delta = T_reference - T_approximate, solve
/// J(T_reference)^T lambda = weights. The finite difference R(T_reference) -
/// R(T_approximate) retains both material states and any matching-P1 contact
/// blocks. The explicitly observed remainder closes the difference to
/// weights . delta; freezing k(T) and ignoring that remainder is not allowed.
/// The reference's small algebraic residual is included in the subtraction,
/// rather than presumed zero. No additional primal solve is performed.
///
/// Geometry, loads, Dirichlet values and fixed contact resistance are shared
/// by both fields. Heterogeneous material assignments take precedence over
/// the fallback material, exactly as in the primal. Approximate fields may
/// cross table segments, but an active kink at the reference refuses a
/// unique tangent. A successful comparison certifies neither continuum
/// accuracy nor the correctness of a material law.
///
/// # Errors
/// Refuses invalid budgets, field lengths, nonfinite arithmetic, changed
/// prescribed temperatures, material extrapolation, undeclared/mismatched
/// contact, a reference failing the primal residual gate, a failed dual,
/// a nonsmooth reference tangent, or cancellation at the existing assembly
/// and Krylov checkpoints.
#[allow(clippy::too_many_arguments)]
pub fn compare_discrete_goal(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig,
    reference_temperature: &[f64],
    approximate_temperature: &[f64],
    full_nodal_weights: &[f64],
) -> Result<DiscreteGoalComparison, ConductionError> {
    poll(cx, 0)?;
    admit_linear(linear)?;
    let n = problem.mesh.vertex_count();
    for values in [
        reference_temperature,
        approximate_temperature,
        full_nodal_weights,
    ] {
        vector(cx, values, n)?;
    }
    let dofs = DofMap::new(problem.boundary, n)?;
    if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
        return Err(ConductionError::SingularPureNeumann);
    }
    for (index, &vertex) in dofs.fixed().iter().enumerate() {
        if index % 512 == 0 {
            poll(cx, index)?;
        }
        if reference_temperature[vertex] != dofs.prescribed()[vertex]
            || approximate_temperature[vertex] != dofs.prescribed()[vertex]
        {
            return Err(invalid(
                "discrete goal comparison requires the same prescribed temperatures in both fields",
            ));
        }
    }
    let dependent = temperature_dependent(cx, problem)?;
    if dependent && linear.restart == 0 {
        return Err(invalid(
            "nonlinear discrete goal comparison requires a positive FGMRES restart",
        ));
    }
    let assemble = |temperature| {
        assemble_operator_scaled_with_interfaces(
            cx,
            problem.mesh,
            problem.boundary,
            problem.material,
            problem.source,
            temperature,
            None,
            interfaces,
            problem.element_materials,
        )
    };
    let reference_system = assemble(reference_temperature)?;
    let (matrix, rhs) = reduce(&reference_system, &dofs);
    let primal_relative_residual =
        true_residual(&matrix, &dofs.gather(reference_temperature), &rhs)?;
    if primal_relative_residual >= linear.tolerance {
        return Err(failed(0, primal_relative_residual, linear));
    }
    let reference_residual =
        crate::assemble::residual(&reference_system, &dofs, reference_temperature);
    let approximate_system = assemble(approximate_temperature)?;
    let approximate_residual =
        crate::assemble::residual(&approximate_system, &dofs, approximate_temperature);
    let (matrix, nonlinear) = if dependent {
        let (jacobian, tangent) =
            nonlinear::prepare(cx, problem, interfaces, reference_temperature, &dofs)?;
        (jacobian, Some(tangent))
    } else {
        (matrix, None)
    };
    // Reuse the exact normalized PCG/FGMRES paths and true-residual gates of
    // the physical response. This private temporary publishes no flux data.
    let response = RobinResponse {
        temperature: reference_temperature.to_vec(),
        robin_fluxes: Vec::new(),
        matrix,
        dofs,
        ports: Vec::new(),
        linear,
        nonlinear,
    };
    let (lambda, dual_relative_residual, dual_iterations) =
        response.solve_rhs(cx, full_nodal_weights, true)?;
    let mut nodal_contributions = vec![0.0; n];
    let mut signed_residual_change = 0.0;
    for (index, &vertex) in response.dofs.free().iter().enumerate() {
        if index % 512 == 0 {
            poll(cx, index)?;
        }
        let difference = checked(reference_residual[index] - approximate_residual[index])?;
        let contribution = checked(lambda[vertex] * difference)?;
        nodal_contributions[vertex] = contribution;
        add(&mut signed_residual_change, contribution)?;
    }
    let mut signed_goal_change = 0.0;
    for vertex in 0..n {
        if vertex % 512 == 0 {
            poll(cx, vertex)?;
        }
        let difference = checked(reference_temperature[vertex] - approximate_temperature[vertex])?;
        add(
            &mut signed_goal_change,
            checked(full_nodal_weights[vertex] * difference)?,
        )?;
    }
    let linearization_remainder = checked(signed_goal_change - signed_residual_change)?;
    poll(cx, dual_iterations)?;
    Ok(DiscreteGoalComparison {
        signed_residual_change,
        linearization_remainder,
        signed_goal_change,
        nodal_contributions,
        primal_relative_residual,
        dual_relative_residual,
        dual_iterations,
        uses_nonlinear_jacobian: dependent,
    })
}
