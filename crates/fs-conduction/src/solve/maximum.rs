//! Bounded goal correction of an already accepted physical solve. A smaller
//! maximum-error bound alone never authorizes replacing its physical report.

use super::{
    ConductionError, ConductionProblem, ConductionSolution, Cx, DofMap, LinearConfig,
    PhysicalReport, StopReason, ThermalInterfaces, assemble_operator_scaled_with_interfaces,
    physical_report,
};
use crate::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalAnalyzer, LinearGoalSolveConfig, LinearGoalStop,
    LinearMaximumAnalysis,
};

/// A goal-improving candidate failed the original physical residual criterion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaximumPhysicalGateRefusal {
    /// Independently recomputed residual of the proposed temperature field.
    pub candidate_residual_w: f64,
    /// Original accepted solve's declared residual threshold.
    pub residual_threshold_w: f64,
}

/// Accepted physical field, its matching maximum analysis, and all correction
/// work. The underlying goal stop is distinct from physical acceptance.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearMaximumPolish {
    /// The adopted candidate or the retained original field, with freshly
    /// computed residual, energy, contact and Robin reporting.
    pub solution: ConductionSolution,
    /// Maximum analysis of exactly `solution.temperature`.
    pub analysis: LinearMaximumAnalysis,
    /// Original field's bound before any correction; absent is not zero.
    pub initial_bound_k: Option<f64>,
    /// Caller-declared absolute maximum-error allowance, K.
    pub requested_tolerance_k: f64,
    /// Whether the ACCEPTED field meets the requested maximum allowance.
    pub goal_met: bool,
    /// A changed, physically accepted correction replaced the original field.
    pub candidate_accepted: bool,
    /// Why goal correction stopped, even if physical acceptance later refused.
    pub stop: LinearGoalStop,
    /// All correction PCG iterations, including rejected proposals.
    pub primal_iterations: usize,
    /// Additional defect solves actually started.
    pub defect_corrections: usize,
    /// Completed goal checks inside the correction driver.
    pub goal_checks: usize,
    /// Present only when a changed candidate failed the physical residual gate.
    pub physical_gate_refusal: Option<MaximumPhysicalGateRefusal>,
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config {
        parameter: "maximum goal polish",
        what: what.to_string(),
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "maximum-goal-polish",
        at: 0,
    })
}

fn physical(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    dofs: &DofMap,
    temperature: &[f64],
) -> Result<PhysicalReport, ConductionError> {
    poll(cx)?;
    let system = assemble_operator_scaled_with_interfaces(
        cx,
        problem.mesh,
        problem.boundary,
        problem.material,
        problem.source,
        temperature,
        None,
        interfaces,
        problem.element_materials,
    )?;
    let report = physical_report(problem, interfaces, &system, dofs, temperature)?;
    let energy = report.energy;
    if ![
        report.final_residual,
        energy.source_w,
        energy.neumann_out_w,
        energy.robin_out_w,
        energy.dirichlet_in_w,
        energy.closure_w,
        energy.scale_w,
        energy.relative_closure(),
    ]
    .iter()
    .all(|value| value.is_finite())
        || energy.scale_w <= 0.0
    {
        return Err(invalid(
            "physical residual or energy scale is nonfinite or invalid",
        ));
    }
    poll(cx)?;
    Ok(report)
}

fn refresh(solution: &mut ConductionSolution, report: PhysicalReport) {
    solution.report.final_residual = report.final_residual;
    solution.report.energy = report.energy;
    solution.report.interface_fluxes = report.interface_fluxes;
    solution.report.robin_fluxes = report.robin_fluxes;
}

/// Correct a fixed linear physical solve to an explicit regional maximum
/// accuracy, preserving its original residual acceptance criterion.
///
/// Reassembles the original field and refuses it if its true residual exceeds
/// its finite declared threshold. The cached maximum analyzer then performs
/// bounded correction. A changed best candidate, including an improvement
/// retained at work exhaustion, is adopted only after independent reassembly
/// passes that same physical threshold. Otherwise the original field and its
/// analysis survive, with an explicit physical-gate refusal. All physical
/// report values are recomputed from the returned field using the same helper
/// as the ordinary solver; original iteration history and linear evidence
/// remain historical, while correction work is counted separately.
///
/// The supplied original report's threshold is the caller's declared policy,
/// not authenticated source evidence. Constant heterogeneous materials and
/// fixed matching contact/Robin rows are supported. Frozen radiation or air
/// feedback is not a bound on the coupled nonlinear equations. Missing inverse
/// evidence returns an unresolved goal with the original physical solution.
///
/// # Errors
/// Invalid policy/threshold, inconsistent mesh/report dimensions, failed
/// baseline residual, nonlinear conductivity, material/assembly/arithmetic
/// refusal and cancellation. No partially refreshed solution is published.
#[allow(clippy::too_many_arguments)]
pub fn polish_linear_maximum(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig,
    original: &ConductionSolution,
    vertices: &[usize],
    analysis_config: LinearGoalAnalysisConfig,
    control_config: LinearGoalSolveConfig,
) -> Result<LinearMaximumPolish, ConductionError> {
    poll(cx)?;
    let threshold = original.report.residual_threshold;
    if !(threshold.is_finite() && threshold >= 0.0)
        || !(original.report.final_residual.is_finite() && original.report.final_residual >= 0.0)
        || !(original.report.energy.scale_w.is_finite() && original.report.energy.scale_w > 0.0)
        || !(control_config.absolute_tolerance.is_finite()
            && control_config.absolute_tolerance > 0.0)
        || !(1..=32).contains(&control_config.check_every)
    {
        return Err(invalid(
            "finite nonnegative residual threshold, positive goal tolerance and check cadence 1..=32 required",
        ));
    }
    let analyzer = LinearGoalAnalyzer::new_for_maximum(
        cx,
        problem,
        interfaces,
        linear,
        &original.temperature,
        analysis_config,
    )?;
    if original.report.elements != problem.mesh.element_count()
        || original.report.free_dofs != analyzer.dofs().n()
    {
        return Err(invalid(
            "original solution report does not match the supplied mesh",
        ));
    }
    let material = problem.element_materials;
    if original.report.material_provenance
        != material.map_or_else(
            || problem.material.provenance(),
            crate::material::ElementMaterials::provenance,
        )
        || original.report.material_receipts
            != material.map_or_else(
                || problem.material.receipts().len(),
                |assigned| assigned.receipts().len(),
            )
        || original.report.element_material_identity
            != material.map(crate::material::ElementMaterials::identity)
    {
        return Err(invalid(
            "original material report does not match the supplied assignment",
        ));
    }
    let baseline_physical = physical(
        cx,
        problem,
        interfaces,
        analyzer.dofs(),
        &original.temperature,
    )?;
    if baseline_physical.final_residual > threshold {
        return Err(ConductionError::NotConverged {
            iterations: original.report.iterations,
            residual: baseline_physical.final_residual,
            threshold,
        });
    }
    let baseline = analyzer.analyze_maximum(cx, &original.temperature, vertices)?;
    let initial_bound_k = baseline.algebraic_half_width_k();
    let correction =
        analyzer.solve_maximum_to_goal(cx, &original.temperature, vertices, control_config)?;
    let mut solution = original.clone();
    refresh(&mut solution, baseline_physical);
    let mut analysis = baseline;
    let mut candidate_accepted = false;
    let mut physical_gate_refusal = None;
    if correction
        .temperature
        .iter()
        .zip(&original.temperature)
        .any(|(a, b)| a.to_bits() != b.to_bits())
    {
        let checked = physical(
            cx,
            problem,
            interfaces,
            analyzer.dofs(),
            &correction.temperature,
        )?;
        if checked.final_residual <= threshold {
            solution.temperature = correction.temperature;
            refresh(&mut solution, checked);
            solution.report.stop_reason = StopReason::ResidualTolerance;
            analysis = correction.analysis;
            candidate_accepted = true;
        } else {
            physical_gate_refusal = Some(MaximumPhysicalGateRefusal {
                candidate_residual_w: checked.final_residual,
                residual_threshold_w: threshold,
            });
        }
    }
    let goal_met = analysis.meets_absolute_tolerance(control_config.absolute_tolerance);
    poll(cx)?;
    Ok(LinearMaximumPolish {
        solution,
        analysis,
        initial_bound_k,
        requested_tolerance_k: control_config.absolute_tolerance,
        goal_met,
        candidate_accepted,
        stop: correction.stop,
        primal_iterations: correction.primal_iterations,
        defect_corrections: correction.defect_corrections,
        goal_checks: correction.goal_checks,
        physical_gate_refusal,
    })
}
