//! Missing solid stability is prepared once before complete solid/air correction.

use fs_conduction::adjoint::{
    LinearGoalSolveConfig, LinearRobinFeedbackAnalyzer, LinearRobinMaximumAnalysis,
    SpectralInverseLimits, SpectralPreparation, SpectralStop,
};

use super::{
    AirPath, ConductionProblem, Cx, LinearAirMaximumSolve, LinearConfig,
    LinearGoalAnalysisConfig, Result, RobinFeedbackAnalysisConfig, ThermalInterfaces,
    bad, poll, prepare_linear_maximum, solve,
};

/// Explicit additional work to prove a solid inverse before coupled correction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralMaximumControl {
    /// Positive finite shift proposal in the stored matrix's coefficient units.
    /// It is NOT an assumed eigenvalue bound; the Gram residual must prove it.
    pub initial_shift: f64,
    /// One preparation allowance, shared by all shifts, factor fill and proof
    /// work. System limits are intersected with the existing solid envelope.
    /// The ordinary feedback allowance remains per residual/goal assessment.
    pub limits: SpectralInverseLimits,
}

/// Read-only-by-convention diagnostics, not independently reusable inverse proof.
/// No solver accepts these scalars in place of the matrix-owning certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralPreparationSummary {
    /// Actual preparation outcome, including an unsuccessful bounded attempt.
    pub stop: SpectralStop,
    /// All preparation visits, including retries; never charged per correction.
    pub work_entries: usize,
    /// Peak admitted logical records, not allocator bytes or measured RSS.
    pub peak_storage_entries: usize,
    /// Actual attempted shifts under the single preparation allowance.
    pub shift_attempts: usize,
    /// Checked shift, absent when no certificate was produced.
    pub shift: Option<f64>,
    /// Full stored-matrix Gram defect upper bound, if proved.
    pub defect_upper: Option<f64>,
    /// Positive stored-system coercivity lower bound, if proved.
    pub coercivity_lower: Option<f64>,
}
impl SpectralPreparationSummary {
    fn from_preparation(preparation: &SpectralPreparation) -> Self {
        let certificate = preparation.certificate.as_ref();
        Self {
            stop: preparation.stop,
            work_entries: preparation.work_entries,
            peak_storage_entries: preparation.peak_storage_entries,
            shift_attempts: preparation.shift_attempts,
            shift: certificate.map(|proof| proof.shift()),
            defect_upper: certificate.map(|proof| proof.defect_upper()),
            coercivity_lower: certificate.map(|proof| proof.coercivity_lower()),
        }
    }
}

/// A complete candidate bundle plus the preliminary check and preparation cost.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralAirMaximumSolve {
    /// Best checked temperature and air states rebuilt from that exact field.
    /// Inspect the retained stop; neither a budget stop nor missing proof is success.
    pub solution: LinearAirMaximumSolve,
    /// The ONE preliminary full-system check that decided whether preparation
    /// was needed. This check is additional to `solution.solid.goal_checks`.
    pub initial_analysis: LinearRobinMaximumAnalysis,
    /// Paid preparation, even if it failed. None means existing inverse evidence
    /// sufficed (or the requested maximum contained no free vertices).
    pub preparation: Option<SpectralPreparationSummary>,
}

/// Correct a linear cooling model whose solid inverse needs a sparse proof.
///
/// Use the ordinary physical AirPath/Robin admission and first check the
/// requested maximum. Only a missing solid inverse triggers ONE bounded
/// preparation. The immutable analyzer then reuses its checked certificate
/// throughout the existing complete-operator FGMRES correction loop; neither
/// factorization nor response solves are repeated for successive fields.
///
/// Preparation has its explicit lifetime work/storage allowance. Each outward
/// field/response assessment retains the existing feedback limits; all primal
/// corrections/restarts share `control.max_primal_iterations`. The additional
/// initial check is retained separately, not hidden in the solve's check count.
/// Failed preparation yields the ordinary honest BoundUnavailable result,
/// retaining its original field and freshly rebuilt production air marches.
///
/// This is opt-in: `solve_linear_maximum` retains its previous preparation and
/// numerical path. Neither function creates a ConductionSolution, validates a
/// physical model, or substitutes a stored-affine bound for continuum error,
/// coefficient-lowering error, nonlinear conductivity, radiation or hydraulics.
///
/// # Errors
/// Invalid controls/proposals, physical or field admission, proof input faults,
/// arithmetic/allocation failures and cancellation. No partial field/air bundle
/// or preparation success is returned after cancellation.
#[allow(clippy::too_many_arguments)]
pub fn solve_linear_maximum_with_spectral(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath],
    linear: LinearConfig,
    initial_temperature: &[f64],
    region_vertices: &[usize],
    solid_config: LinearGoalAnalysisConfig,
    feedback_config: RobinFeedbackAnalysisConfig,
    control: LinearGoalSolveConfig,
    spectral: SpectralMaximumControl,
) -> Result<SpectralAirMaximumSolve> {
    let PreparedMaximum { analyzer, initial_analysis, preparation } = prepare(
        cx, problem, interfaces, paths, linear, initial_temperature, region_vertices,
        solid_config, feedback_config, control, spectral,
    )?;
    let solution = solve::finish(cx, &analyzer, paths, initial_temperature, region_vertices, control)?;
    poll(cx)?;
    Ok(SpectralAirMaximumSolve { solution, initial_analysis, preparation })
}

// Shared by numerical-only correction and physical in-loop admission. No
// second preparation or response solve is needed after a physical rejection.
pub(super) struct PreparedMaximum<'m> {
    pub(super) analyzer: LinearRobinFeedbackAnalyzer<'m>,
    pub(super) initial_analysis: LinearRobinMaximumAnalysis,
    pub(super) preparation: Option<SpectralPreparationSummary>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare<'m>(
    cx: &Cx<'_>, problem: ConductionProblem<'m>, interfaces: Option<&ThermalInterfaces>,
    paths: &[AirPath], linear: LinearConfig, initial_temperature: &[f64],
    region_vertices: &[usize], solid_config: LinearGoalAnalysisConfig,
    feedback_config: RobinFeedbackAnalysisConfig, control: LinearGoalSolveConfig,
    spectral: SpectralMaximumControl,
) -> Result<PreparedMaximum<'m>> {
    solve::admit_control(cx, control)?;
    if !spectral.initial_shift.is_finite() || spectral.initial_shift <= 0.0
        || spectral.limits.max_shift_attempts == 0
    {
        return Err(bad("spectral maximum control needs a positive finite shift and nonzero shift attempts"));
    }
    let mut analyzer = prepare_linear_maximum(
        cx, problem, interfaces, paths, linear, initial_temperature,
        solid_config, feedback_config,
    )?;
    let initial_analysis = analyzer.analyze_maximum(cx, initial_temperature, region_vertices)?;
    if initial_analysis.free_vertices() > 0
        && initial_analysis.coupled().solid_inverse_infinity_upper().is_none()
    {
        analyzer = analyzer.with_spectral_inverse(cx, spectral.initial_shift, spectral.limits)?;
    }
    let preparation = analyzer.spectral_preparation().map(SpectralPreparationSummary::from_preparation);
    poll(cx)?;
    Ok(PreparedMaximum { analyzer, initial_analysis, preparation })
}
