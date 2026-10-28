//! Closed-enclosure exchange with original external boundary heat preserved.
//!
//! Positive Robin patches use h(T-(T_ref-q_rad/h)); prescribed-flux patches
//! add q_rad to their authored outward flux. Both representations share the
//! same exact patch binding, reflected exchange and physical heat accounts.

use super::{
    CoupledRadiationConfig, CoupledRadiationReport, CoupledRadiationSolution, GrayDiffuseEnclosure,
    RadiosityReport, enclosure_exchange as exchange, solve_enclosure,
};
use crate::{
    ConductionError, ConductionProblem, ConductionSolution, RobinFlux, SolveConfig,
    ThermalBoundary, ThermalInterfaces,
};
use fs_exec::Cx;

/// Steady enclosure equilibrium with original external heat kept separate.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvectiveEnclosureSolution {
    /// Accepted field and report for the combined boundary actually solved.
    /// Robin and Neumann totals may both include applied internal radiation.
    pub conduction: ConductionSolution,
    /// Applied radiosity, fresh final-state radiation and the full residual gate.
    pub radiation: CoupledRadiationReport,
    /// Exact frozen-load boundary of the accepted inner solve.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in the original region-declaration order.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convective heat leaving the solid [W]; only this heats the air.
    pub convective_out_w: f64,
    /// Authored outward Neumann heat [W], excluding internal enclosure transfer.
    pub prescribed_neumann_out_w: f64,
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "radiation-convection-binding",
        at,
    })
}

/// Solve steady conduction with reflected radiation on exact named boundaries.
///
/// Every enclosure surface covers one complete same-named uniform positive
/// Robin or uniform Neumann region. Zero outward flux expresses a vacuum
/// patch without introducing convection. Other exterior laws, heterogeneous
/// k(T), matching contact and prescribed-temperature baths remain unchanged.
///
/// The existing initialization requires a solvable conduction problem before
/// radiation is applied. An isolated steady enclosure has no external thermal
/// anchor and retains the ordinary singularity/balance refusal. Radiation-only
/// anchoring of disconnected solid components is not inferred.
///
/// Surface-update convergence alone is insufficient: fresh final radiosity is
/// assembled into the complete nonlinear residual and must pass the original
/// solid threshold. Derivative consumers use the complete enclosure tangent;
/// differentiating only the returned frozen boundary would omit reflection.
///
/// # Errors
/// Refuses mismatched traces, unsupported boundary laws, invalid physical
/// coefficients, source-card violations, cancellation or exhausted solver work.
pub fn solve_with_convective_gray_diffuse_enclosure(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    enclosure: &GrayDiffuseEnclosure,
    conduction_config: SolveConfig,
    coupling_config: CoupledRadiationConfig,
) -> Result<ConvectiveEnclosureSolution, ConductionError> {
    poll(cx, 0)?;
    let patches = exchange::bind_traces(cx, problem, enclosure)?;
    let prescribed_neumann_out_w = exchange::prescribed_neumann_heat(cx, problem)?;
    let apply_boundary = |radiosity: &RadiosityReport| {
        exchange::shifted_boundary(problem, &patches, radiosity)
    };
    let (
        CoupledRadiationSolution { conduction, radiation },
        combined_boundary,
    ) = solve_enclosure(
        cx, problem, interfaces, enclosure, conduction_config, coupling_config,
        &apply_boundary, true,
    )?;
    let (convective_robin_fluxes, convective_out_w) =
        exchange::convection(cx, problem, &conduction.report.robin_fluxes)?;
    poll(cx, convective_robin_fluxes.len())?;
    Ok(ConvectiveEnclosureSolution {
        conduction, radiation, combined_boundary, convective_robin_fluxes,
        convective_out_w, prescribed_neumann_out_w,
    })
}
