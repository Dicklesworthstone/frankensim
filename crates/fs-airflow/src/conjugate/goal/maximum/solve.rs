//! One-shot coupled maximum correction with air states rebuilt from its field.

use fs_conduction::adjoint::{LinearGoalSolve, LinearGoalSolveConfig, LinearRobinMaximumAnalysis};
use crate::conjugate::AirMarch;
use super::{AirPath, ConductionProblem, Cx, LinearConfig, LinearGoalAnalysisConfig,
    Result, RobinFeedbackAnalysisConfig, ThermalInterfaces, bad, poll, prepare_linear_maximum};

/// Coupled temperature candidate and air marches belonging to that exact field.
/// This is not a whole-model `ConductionSolution`, an energy audit, or a
/// continuum/assembly error certificate. A consumer must inspect `solid.stop`
/// and revalidate its physical gates before replacing a published simulation.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearAirMaximumSolve {
    /// Best field and independently checked STORED coupled-system goal bound.
    /// A budget or missing-bound outcome is retained, never changed to success.
    pub solid: LinearGoalSolve<LinearRobinMaximumAnalysis>,
    /// Wall means in branch-major, stream-wise port order.
    pub wall_temperatures_k: Vec<f64>,
    /// Recomputed production air marches in input path order. No inlet,
    /// outlet, heat rate or reference is copied from the initial field.
    pub air: Vec<AirMarch>,
}

/// Correct the full linear solid/air system to an absolute regional-max goal.
///
/// Prepare the existing production AirPath/Robin binding and use its new
/// nonsymmetric FGMRES correction path. Upstream wall dependencies and separate
/// inlets stay active throughout. All correction work shares `control`'s
/// iteration cap; response and inverse preparation keep their own existing
/// explicitly supplied caps and are not repeated between corrections.
///
/// Recompute every wall mean and original exponential-law air march from the
/// returned best temperature, even on a non-success outcome. The existing
/// bound concerns the rounded stored affine model, not the difference between
/// coefficient lowering and the production march, nonlinear material laws,
/// radiation, flow redistribution, or physical/continuum accuracy. This routine
/// does not manufacture a new physical convergence or energy-balance report.
///
/// # Errors
/// Preparation/control/field/selection, air-law, finite-range/allocation and
/// cancellation refusals. No partial temperature/air bundle is published.
#[allow(clippy::too_many_arguments)]
pub fn solve_linear_maximum(
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
) -> Result<LinearAirMaximumSolve> {
    poll(cx)?;
    // Reject invalid correction controls before paying for preparation.
    if !(control.absolute_tolerance.is_finite() && control.absolute_tolerance > 0.0)
        || !(1..=32).contains(&control.check_every)
    {
        return Err(bad("air maximum solve needs a positive finite tolerance and check_every in 1..=32"));
    }
    let analyzer = prepare_linear_maximum(
        cx, problem, interfaces, paths, linear, initial_temperature,
        solid_config, feedback_config,
    )?;
    let solid = analyzer.solve_maximum_to_goal(cx, initial_temperature, region_vertices, control)?;
    let wall_temperatures_k = analyzer.wall_mean_temperatures(cx, &solid.temperature)?;
    let mut air = Vec::new();
    air.try_reserve_exact(paths.len()).map_err(|_| bad("air maximum result allocation refused"))?;
    let mut start = 0;
    for path in paths {
        poll(cx)?;
        // The preparation established the exact total and coordinate order.
        let end = start + path.segments().len();
        air.push(path.march(&wall_temperatures_k[start..end])?);
        start = end;
    }
    poll(cx)?;
    Ok(LinearAirMaximumSolve { solid, wall_temperatures_k, air })
}
