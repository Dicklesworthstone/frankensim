//! Stored-system goal analysis on actual conduction assembly. The nonlinear,
//! coupled-airflow, radiation and maximum-relocation gaps are not hidden by
//! freezing them and changing the meaning of a certificate.

use fs_solver::goal::{
    GoalResidualError, GoalResidualLimits, GoalResidualReport, enclose_goal_error,
};

use super::{
    ConductionError, ConductionProblem, Cx, DofMap, LinearConfig, RobinResponse,
    ThermalInterfaces, admit_linear, assemble_operator_scaled_with_interfaces,
    invalid, poll, reduce, temperature_dependent, vector,
};

mod solve;
pub use solve::{LinearGoalSolve, LinearGoalSolveConfig, LinearGoalStop};

mod maximum;
pub use maximum::{LinearMaximumAnalysis, analyze_linear_maximum};

/// Additional work admitted for a discrete thermal goal analysis. The dual
/// separately uses the caller's existing [`LinearConfig`] iteration budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LinearGoalAnalysisConfig {
    /// Limits for every outward residual pass, not a replacement for the
    /// assembly/runtime memory and work budgets.
    pub residual_limits: GoalResidualLimits,
    /// Extra CG iterations allowed to propose a positive stability scaling.
    /// Zero disables that proposal, not the inverse verification itself.
    pub max_stability_iterations: usize,
}

/// Algebraic error for one declared linear temperature goal on one retained
/// discrete operator. A bound on a selected node is NOT a bound on a maximum
/// whose active node may move. Boundary coefficients and references are fixed.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearGoalAnalysis {
    /// Outward arithmetic, dual-error and checked inverse evidence.
    pub enclosure: GoalResidualReport,
    /// Recomputed relative residual of the physical transposed solve. It may
    /// exceed the requested Krylov target when the explicit budget ends; the
    /// actual dual defect is retained in the full goal enclosure.
    pub dual_relative_residual: f64,
    /// Inner Krylov work used for the dual.
    pub dual_iterations: usize,
    /// Inner Krylov work actually spent proposing a stability scaling.
    pub stability_iterations: usize,
    /// Recomputed proposal residual, which need not meet its loose target.
    pub stability_relative_residual: Option<f64>,
    /// The positive proposal checked against the exact retained matrix, in
    /// free-dof order. Presence alone does not establish an inverse bound.
    pub stability_scaling: Option<Vec<f64>>,
}

impl LinearGoalAnalysis {
    /// Goal-oriented stopping test, in the declared linear goal's units.
    /// A missing full algebraic bound or an invalid tolerance cannot pass.
    /// This says nothing about discretization, assembly or material error.
    #[must_use]
    pub fn meets_absolute_tolerance(&self, tolerance: f64) -> bool {
        tolerance.is_finite()
            && tolerance > 0.0
            && self.enclosure.goal_error().is_some_and(|bound| {
                bound.magnitude_upper() <= tolerance
            })
    }
}

// A finite inexact dual is useful: its omitted term is bounded separately.
// Requiring an over-solved dual would defeat goal-oriented primal stopping.
fn bounded_solve(
    cx: &Cx<'_>, matrix: &fs_sparse::Csr, rhs: &[f64],
    linear: LinearConfig, transpose: bool,
) -> Result<(Vec<f64>, f64, usize), ConductionError> {
    use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, checked_cg};
    poll(cx, 0)?;
    // Retain the ACTUAL stored transpose, including assembly roundoff
    // asymmetry; an approximately symmetric matrix is not its own transpose.
    let matrix = if transpose { fs_sparse::ops::transpose(matrix) } else { matrix.clone() };
    let op = CsrOp::symmetric(matrix);
    let pre = crate::solve::spd_preconditioner(op.matrix());
    poll(cx, 0)?;
    let result = checked_cg(&op, &pre, rhs, CheckedCgConfig {
        tolerance: linear.tolerance, max_iterations: linear.max_iterations,
        max_corrections: 0,
    }, |iteration| poll(cx, iteration)).map_err(|error| match error {
        CheckedCgError::InvalidInput(what) => invalid(what),
        CheckedCgError::Interrupted(error) => error,
    })?;
    Ok((result.x, result.report.rel_residual, result.report.iters))
}

fn map_enclosure(error: GoalResidualError) -> ConductionError {
    match error {
        GoalResidualError::Cancelled => ConductionError::Cancelled {
            stage: "linear-goal-enclosure",
            at: 0,
        },
        error => ConductionError::Config {
            parameter: "linear-goal-enclosure",
            what: error.to_string(),
        },
    }
}

/// A fixed linear thermal operator, goal and computed dual that can assess
/// many iterates without repeating a dual or stability solve. Immutable
/// borrows bind material/mesh/boundary identity for the analyzer's lifetime.
/// The outward inverse check still runs against the stored matrix each time.
pub struct LinearGoalAnalyzer<'m> {
    problem: ConductionProblem<'m>,
    response: RobinResponse,
    rhs: Vec<f64>,
    weights: Vec<f64>,
    free_dual: Vec<f64>,
    config: LinearGoalAnalysisConfig,
    dual_relative_residual: f64,
    dual_iterations: usize,
    stability_iterations: usize,
    stability_relative_residual: Option<f64>,
    stability_scaling: Option<Vec<f64>>,
}

fn validate_field(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    dofs: &DofMap,
    temperature: &[f64],
) -> Result<(), ConductionError> {
    vector(cx, temperature, problem.mesh.vertex_count())?;
    for (index, &vertex) in dofs.fixed().iter().enumerate() {
        if index % 512 == 0 { poll(cx, index)?; }
        if temperature[vertex] != dofs.prescribed()[vertex] {
            return Err(invalid("linear goal analysis requires the declared prescribed temperatures"));
        }
    }
    for element in 0..problem.mesh.element_count() {
        if element % 512 == 0 { poll(cx, element)?; }
        let model = match problem.element_materials {
            Some(materials) => materials.model_for(element)?,
            None => problem.material,
        };
        model.temperature_span().check(crate::assemble::element_temperature(
            problem.mesh, element, temperature,
        ))?;
    }
    Ok(())
}

impl<'m> LinearGoalAnalyzer<'m> {
    /// Prepare actual production assembly and a bounded transposed dual solve.
    ///
    /// Supports heterogeneous constant tensors, the nonzero Dirichlet lift,
    /// fixed Robin boundaries and matching-P1 contact. An unconverged
    /// reference primal is intentional; no primal convergence is presumed.
    ///
    /// If unscaled dominance cannot prove an inverse bound, the optional
    /// extra solve `A w = 1` proposes a positive scaling. Its loose 1% residual
    /// target only controls proposal cost. Outward verification of EVERY row,
    /// not that target, establishes the inverse bound. An insufficient
    /// proposal leaves the full bound absent. No stability constant is supplied
    /// by the caller or invented from a residual tolerance.
    ///
    /// # Errors
    /// Refuses malformed inputs, changed prescribed temperatures, nonlinear
    /// conductivity, material extrapolation, assembly failure, exhausted
    /// structural limits, nonfinite arithmetic and cancellation. A finite
    /// inexact dual is retained with its error term, not upgraded to an exact
    /// dual. An insufficient stability proposal is only a no-bound state.
    pub fn new(
        cx: &Cx<'_>,
        problem: ConductionProblem<'m>,
        interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig,
        reference_temperature: &[f64],
        full_nodal_weights: &[f64],
        config: LinearGoalAnalysisConfig,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        admit_linear(linear)?;
        let n = problem.mesh.vertex_count();
        vector(cx, full_nodal_weights, n)?;
        if temperature_dependent(cx, problem)? {
            return Err(invalid(
                "an algebraic goal certificate requires a linear conductivity operator; \
                 a frozen nonlinear tangent is not the original problem",
            ));
        }
        let dofs = DofMap::new(problem.boundary, n)?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        if dofs.n() > config.residual_limits.max_rows {
            return Err(invalid("linear goal analysis exceeds its declared row budget"));
        }
        validate_field(cx, problem, &dofs, reference_temperature)?;
        let system = assemble_operator_scaled_with_interfaces(
            cx, problem.mesh, problem.boundary, problem.material, problem.source,
            reference_temperature, None, interfaces, problem.element_materials,
        )?;
        let (matrix, rhs) = reduce(&system, &dofs);
        if matrix.nnz() > config.residual_limits.max_nonzeros {
            return Err(invalid("linear goal analysis exceeds its declared nonzero budget"));
        }
        let free_temperature = dofs.gather(reference_temperature);
        let weights = dofs.gather(full_nodal_weights);
        let response = RobinResponse {
            temperature: reference_temperature.to_vec(),
            robin_fluxes: Vec::new(),
            matrix,
            dofs,
            ports: Vec::new(),
            linear,
            nonlinear: None,
        };
        let (free_dual, dual_relative_residual, dual_iterations) =
            bounded_solve(cx, &response.matrix, &weights, linear, true)?;
        let evaluate = |response: &RobinResponse, scaling: Option<&[f64]>| {
            enclose_goal_error(
                &response.matrix, &rhs, &free_temperature, &weights, &free_dual,
                scaling, config.residual_limits, || cx.checkpoint().is_ok(),
            ).map_err(map_enclosure)
        };
        let enclosure = evaluate(&response, None)?;
        let mut stability_iterations = 0;
        let mut stability_relative_residual = None;
        let mut stability_scaling = None;
        if enclosure.goal_error().is_none() && config.max_stability_iterations > 0 {
            let one = vec![1.0; response.dofs.n()];
            let proposal_config = LinearConfig {
                tolerance: 0.01, max_iterations: config.max_stability_iterations,
                restart: linear.restart,
            };
            let (scaling, residual, iterations) =
                bounded_solve(cx, &response.matrix, &one, proposal_config, false)?;
            stability_iterations = iterations;
            stability_relative_residual = Some(residual);
            if scaling.iter().all(|w| w.is_finite() && *w > 0.0) {
                // Neither positivity nor a small proposal residual is enough.
                // Check the stored rows now and again on every assessed field.
                let _ = evaluate(&response, Some(&scaling))?;
                stability_scaling = Some(scaling);
            }
        }
        poll(cx, dual_iterations.saturating_add(stability_iterations))?;
        Ok(Self {
            problem, response, rhs, weights, free_dual, config,
            dual_relative_residual, dual_iterations, stability_iterations,
            stability_relative_residual, stability_scaling,
        })
    }

    /// Degree-of-freedom identity of the retained reduced operator.
    #[must_use]
    pub const fn dofs(&self) -> &DofMap {
        &self.response.dofs
    }

    /// Actual approximate dual, in free-dof order, for independent replay.
    #[must_use]
    pub fn free_dual(&self) -> &[f64] {
        &self.free_dual
    }

    /// Positive scaling proposal, if one was admitted for independent checking.
    /// Presence alone is not an inverse bound; inspect the analysis result.
    #[must_use]
    pub fn stability_scaling(&self) -> Option<&[f64]> {
        self.stability_scaling.as_deref()
    }

    /// Assess another admissible field on the SAME operator and fixed goal.
    /// No primal, dual, or witness solve is repeated. Prescribed temperatures
    /// and every material's temperature support are rechecked before use.
    ///
    /// # Errors
    /// Field/temperature/Dirichlet, outward-arithmetic and cancellation refusals.
    pub fn analyze(
        &self,
        cx: &Cx<'_>,
        approximate_temperature: &[f64],
    ) -> Result<LinearGoalAnalysis, ConductionError> {
        poll(cx, 0)?;
        validate_field(cx, self.problem, &self.response.dofs, approximate_temperature)?;
        let free_temperature = self.response.dofs.gather(approximate_temperature);
        let enclosure = enclose_goal_error(
            &self.response.matrix, &self.rhs, &free_temperature, &self.weights,
            &self.free_dual, self.stability_scaling.as_deref(),
            self.config.residual_limits, || cx.checkpoint().is_ok(),
        ).map_err(map_enclosure)?;
        poll(cx, 0)?;
        Ok(LinearGoalAnalysis {
            enclosure,
            dual_relative_residual: self.dual_relative_residual,
            dual_iterations: self.dual_iterations,
            stability_iterations: self.stability_iterations,
            stability_relative_residual: self.stability_relative_residual,
            stability_scaling: self.stability_scaling.clone(),
        })
    }
}

/// One-shot version of [`LinearGoalAnalyzer`], without another primal solve.
///
/// Weights use full nodal order. Prescribed-node weights add the same constant
/// to both fields and are eliminated. A selected-node goal is NOT a relocating
/// maximum. This does not differentiate airflow feedback, admit nonlinear k(T)
/// or radiation, certify assembly error, or carry a continuum-error claim.
///
/// # Errors
/// Every preparation and field-assessment refusal of [`LinearGoalAnalyzer`].
pub fn analyze_linear_goal(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig,
    approximate_temperature: &[f64],
    full_nodal_weights: &[f64],
    config: LinearGoalAnalysisConfig,
) -> Result<LinearGoalAnalysis, ConductionError> {
    LinearGoalAnalyzer::new(
        cx, problem, interfaces, linear, approximate_temperature,
        full_nodal_weights, config,
    )?.analyze(cx, approximate_temperature)
}
