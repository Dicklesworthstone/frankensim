//! Selected component means through the original whole-assembly primal and dual.
//! A region indicator is cellwise constant and generally discontinuous. Assemble
//! its exact P1 load per selected element; turning it into a nodal source would
//! smear the functional across shared vertices and change the requested mean.
use super::*;
use crate::{DofMap, InterfaceSurface, LinearConfig, ThermalInterfaces};
use crate::adjoint::{LinearGoalAnalysis, LinearGoalAnalysisConfig, LinearGoalAnalyzer};
pub use fs_solver::goal::GoalResidualLimits;
pub use fs_verify::tet::{RegionMeanBound, RegionMeanSelection};

/// Independent physical-dual, matrix-analysis and flux-reconstruction limits.
/// No stability-proposal solves are needed for this continuum majorant.
#[derive(Debug, Clone, Copy)]
pub struct RegionMeanConfig {
    pub dual: LinearConfig,
    pub residual_limits: GoalResidualLimits,
    pub flux: FluxBudget,
}

/// Actual approximate adjoint and its production analysis, with a continuum
/// regional bound. No primal solve or report is fabricated for a supplied field.
#[derive(Debug)]
pub struct RegionMeanFieldBound {
    /// Full nodal order, with the homogeneous Dirichlet lift. Not normalized
    /// by region volume: the adjoint's functional is the regional INTEGRAL.
    pub dual_temperature: Vec<f64>,
    /// Original stored-system report, in integral units. This is NOT substituted
    /// for the complete continuum bound below; a finite inexact dual is allowed.
    pub dual_analysis: LinearGoalAnalysis,
    pub bound: RegionMeanBound,
}

/// The unchanged physical solver's result and its component-mean analysis.
#[derive(Debug)]
pub struct RegionMeanSolution {
    pub primal: ConductionSolution,
    pub region: RegionMeanFieldBound,
}

fn selection(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, cells: &[usize], config: RegionMeanConfig,
) -> Result<RegionMeanSelection> {
    poll(cx)?;
    if !config.dual.tolerance.is_finite() || config.dual.tolerance <= 0.0
        || config.dual.tolerance >= 1.0 || config.dual.max_iterations == 0 {
        return Err(TetError::Invalid("regional dual needs a positive iteration budget and tolerance in (0,1)").into());
    }
    let selected = RegionMeanSelection::new(problem.mesh.element_count(), cells,
        config.flux, || cx.checkpoint().is_ok())?;
    if DofMap::new(problem.boundary, problem.mesh.vertex_count())?.n() > config.residual_limits.max_rows {
        return Err(TetError::Budget.into());
    }
    Ok(selected)
}

fn bound_field(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, admitted: &Admitted,
    interfaces: &ThermalInterfaces, temperature: &[f64],
    selected: &RegionMeanSelection, config: RegionMeanConfig,
) -> Result<RegionMeanFieldBound> {
    poll(cx)?;
    let n = problem.mesh.vertex_count();
    if temperature.len() != n { return Err(TetError::Invalid("candidate temperature length").into()); }
    for chunk in temperature.chunks(256) {
        poll(cx)?;
        if chunk.iter().any(|v| !v.is_finite()) {
            return Err(TetError::Invalid("finite candidate temperature required").into());
        }
    }
    for &(vertex, value) in problem.boundary.dirichlet() {
        poll(cx)?;
        if temperature[vertex] != value {
            return Err(TetError::Invalid("candidate does not match prescribed Dirichlet trace").into());
        }
    }
    let mut weights = vec![0.0; n];
    for &cell in selected.cells() {
        poll(cx)?;
        // Integral_cell lambda_i = volume/4. Use the geometry owner already
        // consumed by native assembly, not equal vertex/cell averages. These
        // f64 coefficients merely propose a dual; the verifier independently
        // bounds its error for the EXACT cell indicator and outward volume.
        let mass = problem.mesh.element_volume(cell) / 4.0;
        if !mass.is_finite() || mass <= 0.0 { return Err(TetError::Unbounded.into()); }
        for &vertex in &problem.mesh.complex().tets[cell] {
            let value = &mut weights[vertex as usize];
            *value += mass;
            if !value.is_finite() { return Err(TetError::Unbounded.into()); }
        }
    }
    let analyzer = LinearGoalAnalyzer::new(cx, problem, Some(interfaces), config.dual,
        temperature, &weights, LinearGoalAnalysisConfig {
            residual_limits: config.residual_limits, max_stability_iterations: 0,
        })?;
    let mut dual_temperature = vec![0.0; n];
    for (i, (&vertex, &value)) in analyzer.dofs().free().iter().zip(analyzer.free_dual()).enumerate() {
        if i % 256 == 0 { poll(cx)?; }
        dual_temperature[vertex] = value;
    }
    let dual_analysis = analyzer.analyze(cx, temperature)?;
    drop(analyzer);
    let bound = admitted.problem(problem.mesh.positions()).region_mean_bound(
        temperature, &dual_temperature, selected, config.flux, || cx.checkpoint().is_ok())?;
    poll(cx)?;
    Ok(RegionMeanFieldBound { dual_temperature, dual_analysis, bound })
}

/// Bound the selected component's mean without repeating the supplied field's
/// primal solve. Pass every original matching contact surface (or `&[]` for a
/// no-contact model); unselected solids still carry their original heat paths.
/// Cells use ORIGINAL mesh-element indices, not surface slots or vertex IDs.
///
/// The original contact/material/source/boundary admission is retained. A
/// finite inexact dual is reported as such and its continuum error remains in
/// the bound; a small discrete residual is not itself a continuum certificate.
///
/// # Errors
/// Invalid cells, original physical/geometry/solver refusals, matrix or flux
/// budget exhaustion, nonfinite arithmetic, or cancellation return no partial
/// bound. Nominal linear tensor/P1-source/matching-contact scope only; no point
/// maximum, nonlinear, CAD or physical-uncertainty certificate is inferred.
pub fn bound_temperature_region_mean(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, surfaces: &[InterfaceSurface],
    temperature: &[f64], cells: &[usize], config: RegionMeanConfig,
) -> Result<RegionMeanFieldBound> {
    let selected = selection(cx, problem, cells, config)?;
    let (admitted, interfaces) = contact::admit_contacts(cx, problem, surfaces, config.flux)?;
    bound_field(cx, problem, &admitted, &interfaces, temperature, &selected, config)
}

/// Run the original thermal primal once, then its regional-integral adjoint
/// and full-domain majorants. No new boundary is introduced around the region.
/// The returned primal is the actual `solve_with_interfaces` result, unchanged.
/// All source/material/card and mapped contact values are retained together.
///
/// # Errors
/// The same complete-model refusals as `bound_temperature_region_mean`, plus
/// physical primal failure. Selection and known structural limits are admitted
/// before primal work; the actual matrix nonzero limit is checked at assembly.
pub fn solve_with_region_mean_bound(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, surfaces: &[InterfaceSurface],
    cells: &[usize], primal_config: SolveConfig, config: RegionMeanConfig,
) -> Result<RegionMeanSolution> {
    let selected = selection(cx, problem, cells, config)?;
    let (admitted, interfaces) = contact::admit_contacts(cx, problem, surfaces, config.flux)?;
    let primal = crate::solve_with_interfaces(cx, problem, &interfaces, primal_config)?;
    let region = bound_field(cx, problem, &admitted, &interfaces,
        &primal.temperature, &selected, config)?;
    Ok(RegionMeanSolution { primal, region })
}
