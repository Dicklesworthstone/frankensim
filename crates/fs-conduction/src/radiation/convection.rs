//! Coincident convection and closed-enclosure radiation on the original trace.
//!
//! For one frozen radiosity iterate, h(T-T_air)+q_rad is exactly represented
//! by h(T-(T_air-q_rad/h)). The shifted reference is algebra, not an air
//! temperature. Only the separately reconstructed ORIGINAL convection heats
//! the air; radiation remains internal exchange between the enclosure patches.

use super::{
    CoupledRadiationConfig, CoupledRadiationReport, CoupledRadiationSolution, GrayDiffuseEnclosure,
    RadiosityReport, radiation_error, require_temperature, solve_enclosure,
};
use crate::{
    ConductionError, ConductionProblem, ConductionSolution, RobinFlux, ScalarField, SolveConfig,
    ThermalBc, ThermalBoundary, ThermalInterfaces,
};
use fs_exec::Cx;

/// Steady enclosure equilibrium with convection kept separate from radiation.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvectiveEnclosureSolution {
    /// Accepted field and report for the combined boundary actually solved.
    /// Its Robin heat includes applied radiation; do not pass that total to air.
    pub conduction: ConductionSolution,
    /// Applied radiosity, fresh final-state radiation and the full residual gate.
    pub radiation: CoupledRadiationReport,
    /// Exact frozen-load boundary of the accepted inner solve. Its shifted
    /// Robin references are not physical reservoir temperatures.
    pub combined_boundary: ThermalBoundary,
    /// Original convection only, in the original region-declaration order.
    /// Unrelated nonuniform Robin regions retain their original quadrature.
    pub convective_robin_fluxes: Vec<RobinFlux>,
    /// Original convective heat leaving the solid [W]; only this heats the air.
    pub convective_out_w: f64,
}

struct Patch {
    region: usize,
    htc: f64,
    reference: f64,
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "radiation-convection-binding",
        at,
    })
}

/// Solve steady conduction with convection and reflected radiation on the
/// SAME named surface, using the existing bounded radiosity/solid fixed point.
///
/// Every enclosure surface must exactly cover one original uniform positive
/// Robin region with the same name. This explicit entry point declares the two
/// physical channels; the radiation-only solver still refuses occupied faces.
/// Other boundary regions, heterogeneous k(T), and contact remain unchanged.
/// The base problem must already admit a conduction solution without radiation.
///
/// All radiative temperatures remain subject to the original material-card
/// limits. Surface-update convergence alone is insufficient: the complete
/// nonlinear residual is reassembled using fresh final-field radiosity and
/// must satisfy the original solid threshold. No radiation adjoint is implied;
/// differentiating the returned shifted Robin rows freezes radiation feedback
/// and is not a derivative of this coupled model.
///
/// # Errors
/// Mismatched/incomplete traces, nonuniform or non-Robin radiating rows,
/// nonphysical reservoir temperatures, nonfinite composition, cancellation,
/// original material/solver refusals, or an exhausted outer budget.
pub fn solve_with_convective_gray_diffuse_enclosure(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    enclosure: &GrayDiffuseEnclosure,
    conduction_config: SolveConfig,
    coupling_config: CoupledRadiationConfig,
) -> Result<ConvectiveEnclosureSolution, ConductionError> {
    poll(cx, 0)?;
    enclosure.validate_for(problem.mesh)?;
    let mut patches = Vec::with_capacity(enclosure.surfaces().len());
    for surface in enclosure.surfaces() {
        poll(cx, patches.len())?;
        let invalid = |what: &str| {
            radiation_error(
                surface.name(),
                what,
                "bind each enclosure patch to its exact same-named uniform convective region",
            )
        };
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == surface.name())
            .ok_or_else(|| invalid("enclosure patch has no same-named boundary region"))?;
        let ThermalBc::Robin {
            htc: ScalarField::Uniform(htc),
            t_ref: ScalarField::Uniform(reference),
        } = problem.boundary.conditions()[region]
        else {
            return Err(invalid(
                "enclosure convection requires a uniform Robin region",
            ));
        };
        if !htc.is_finite() || htc <= 0.0 {
            return Err(invalid("enclosure convection requires finite positive htc"));
        }
        require_temperature(surface.name(), "convective reservoir", reference)?;
        let mut faces = 0usize;
        for slot in 0..problem.mesh.boundary().len() {
            if slot % crate::assemble::ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            faces += usize::from(problem.boundary.region_for(slot) == Some(region));
        }
        if faces != surface.face_slots().len()
            || surface
                .face_slots()
                .iter()
                .any(|&slot| problem.boundary.region_for(slot) != Some(region))
        {
            return Err(invalid(
                "enclosure patch must exactly cover its convective region",
            ));
        }
        patches.push(Patch {
            region,
            htc,
            reference,
        });
    }
    let apply_boundary = |radiosity: &RadiosityReport| {
        let rows = patches
            .iter()
            .zip(&radiosity.net_outward_flux_w_m2)
            .map(|(patch, &flux)| (patch.region, patch.htc, patch.reference - flux / patch.htc))
            .collect::<Vec<_>>();
        problem.boundary.with_uniform_robin_replacements(&rows)
    };
    let (
        CoupledRadiationSolution {
            conduction,
            radiation,
        },
        combined_boundary,
    ) = solve_enclosure(
        cx,
        problem,
        interfaces,
        enclosure,
        conduction_config,
        coupling_config,
        &apply_boundary,
        true,
    )?;
    let mut convective_robin_fluxes = conduction.report.robin_fluxes.clone();
    let mut convective_out_w = 0.0;
    for (ordinal, row) in convective_robin_fluxes.iter_mut().enumerate() {
        poll(cx, ordinal)?;
        if let Some(patch) = patches
            .iter()
            .find(|patch| problem.boundary.region_names()[patch.region] == row.region)
        {
            row.mean_reference_temperature_k = patch.reference;
            row.mean_htc_w_per_m2_k = patch.htc;
            row.heat_rate_w =
                patch.htc * row.area_m2 * (row.mean_wall_temperature_k - patch.reference);
        }
        convective_out_w += row.heat_rate_w;
        if !row.heat_rate_w.is_finite() || !convective_out_w.is_finite() {
            return Err(radiation_error(
                &row.region,
                "nonfinite original convective heat",
                "retain finite supported physical coefficients and temperatures",
            ));
        }
    }
    poll(cx, convective_robin_fluxes.len())?;
    Ok(ConvectiveEnclosureSolution {
        conduction,
        radiation,
        combined_boundary,
        convective_robin_fluxes,
        convective_out_w,
    })
}
