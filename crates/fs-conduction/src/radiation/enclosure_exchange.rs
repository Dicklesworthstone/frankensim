//! Shared trace binding and physical heat accounts for closed enclosure steps.
//! Storage owners supply their own accepted fields; reflection, shifted loads
//! and original external heat retain one implementation for every trajectory.

use fs_exec::Cx;
use crate::{
    AmbientRadiationConfig, ConductionError, ConductionProblem, RobinFlux,
    ScalarField, ThermalBc, ThermalBoundary,
    assemble::ASSEMBLY_TILE,
};
use super::{GrayDiffuseEnclosure, RadiosityReport};

/// Accepted exchange and actual work over all enclosure trials.
#[derive(Debug, Clone)]
pub struct EnclosureExchangeReport {
    /// Fresh radiosity evaluated at the accepted endpoint temperatures.
    pub radiosity: RadiosityReport,
    /// Signed patch watts actually supplied to the accepted inner solve.
    pub applied_heat_w: Vec<f64>,
    /// Maximum applied-versus-fresh patch heat difference, W.
    pub max_heat_mismatch_w: f64,
    /// Number of inner solid solves; none advances history independently.
    pub iterations: usize,
    /// Sum of actual Newton attempts across all inner solves.
    pub solid_iterations: usize,
    /// Sum of actual Newton linear iterations across all inner solves.
    pub krylov_iterations: usize,
}

pub(crate) struct Patch {
    pub(crate) region: usize,
    pub(crate) law: PatchLaw,
}

#[derive(Clone, Copy)]
pub(crate) enum PatchLaw {
    Convection { htc: f64, reference: f64 },
    PrescribedFlux { outward_flux: f64 },
}

pub(crate) fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
) -> Result<Vec<Patch>, ConductionError> {
    let patches = bind_traces(cx, problem, enclosure)?;
    for (index, surface) in enclosure.surfaces().iter().enumerate() {
        poll(cx, index)?;
        surface.mean_temperature(problem.mesh, temperature)?;
    }
    Ok(patches)
}

/// Bind the exact boundary partition without inventing a temperature seed.
/// Steady coupling owns its physical initialization; storage owners separately
/// validate their supplied candidate through bind and means.
pub(crate) fn bind_traces(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    enclosure: &GrayDiffuseEnclosure,
) -> Result<Vec<Patch>, ConductionError> {
    enclosure.validate_for(problem.mesh)?;
    if enclosure.surfaces().is_empty() || enclosure.surfaces().len() > 64 {
        return Err(invalid("closed enclosure requires 1..=64 patches"));
    }
    let mut patches = Vec::with_capacity(enclosure.surfaces().len());
    for surface in enclosure.surfaces() {
        poll(cx, patches.len())?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| name == surface.name())
            .ok_or_else(|| invalid("unknown closed enclosure region"))?;
        let law = match problem.boundary.conditions()[region] {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(htc),
                t_ref: ScalarField::Uniform(reference),
            } if htc.is_finite() && htc > 0.0 => {
                super::require_temperature(surface.name(), "original convective reference", reference)?;
                PatchLaw::Convection { htc, reference }
            },
            ThermalBc::Neumann { outward_flux: ScalarField::Uniform(outward_flux) } => {
                PatchLaw::PrescribedFlux { outward_flux: finite(outward_flux)? }
            }
            _ => {
                return Err(invalid("enclosure requires uniform Robin or Neumann patches"));
            }
        };
        let mut count = 0_usize;
        for slot in 0..problem.mesh.boundary().len() {
            if slot % ASSEMBLY_TILE == 0 {
                poll(cx, slot)?;
            }
            if problem.boundary.region_for(slot) == Some(region) {
                count += 1;
            }
        }
        if count != surface.face_slots().len()
            || surface
                .face_slots()
                .iter()
                .any(|&slot| problem.boundary.region_for(slot) != Some(region))
        {
            return Err(invalid("enclosure patch must exactly cover its named boundary region"));
        }
        patches.push(Patch { region, law });
    }
    Ok(patches)
}

pub(crate) fn means(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
) -> Result<Vec<f64>, ConductionError> {
    enclosure
        .surfaces()
        .iter()
        .enumerate()
        .map(|(i, surface)| {
            poll(cx, i)?;
            Ok(surface.mean_temperature(problem.mesh, temperature)?)
        })
        .collect()
}

pub(crate) fn heat_tolerance(config: AmbientRadiationConfig, a: f64, b: f64) -> Result<f64, ConductionError> {
    finite(
        config
            .balance_relative_tolerance
            .mul_add(a.abs().max(b.abs()), config.balance_tolerance_w),
    )
}

pub(crate) fn radiosity(
    cx: &Cx<'_>,
    enclosure: &GrayDiffuseEnclosure,
    temperature: &[f64],
    config: AmbientRadiationConfig,
) -> Result<RadiosityReport, ConductionError> {
    poll(cx, 0)?;
    let report = enclosure.solve(cx, temperature)?;
    for &value in report.radiosity_w_m2.iter().chain(&report.irradiation_w_m2) {
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("nonphysical enclosure radiosity"));
        }
    }
    for &value in report
        .net_outward_flux_w_m2
        .iter()
        .chain(&report.net_outward_heat_w)
    {
        finite(value)?;
    }
    let area = enclosure
        .surfaces()
        .iter()
        .map(|surface| surface.area_m2())
        .fold(0.0_f64, f64::max);
    let tolerance = heat_tolerance(config, report.enclosure_energy_scale_w, 0.0)?;
    if finite(report.linear_residual_max_w_m2 * area)? > tolerance
        || finite(report.enclosure_energy_closure_w)?.abs() > tolerance
    {
        return Err(invalid("enclosure radiosity equations or closed energy sum exceed their allowance"));
    }
    Ok(report)
}

pub(crate) fn shifted_boundary(
    problem: ConductionProblem<'_>,
    patches: &[Patch],
    radiosity: &RadiosityReport,
) -> Result<ThermalBoundary, ConductionError> {
    let mut robin = Vec::new();
    let mut neumann = Vec::new();
    for (patch, &flux) in patches.iter().zip(&radiosity.net_outward_flux_w_m2) {
        match patch.law {
            PatchLaw::Convection { htc, reference } => {
                robin.push((patch.region, htc, finite(reference - flux / htc)?));
            }
            PatchLaw::PrescribedFlux { .. } => neumann.push((patch.region, flux)),
        }
    }
    Ok(problem.boundary.with_uniform_robin_replacements(&robin)?
        .with_uniform_neumann_additions(&neumann)?)
}

// Integrate the authored boundary independently. A mixed Robin/Neumann
// enclosure generally transfers nonzero net radiative heat through either
// subset, so the inner solve's Neumann watts cannot be treated as external.
pub(crate) fn prescribed_neumann_heat(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
) -> Result<f64, ConductionError> {
    let mut heat = 0.0;
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % ASSEMBLY_TILE == 0 { poll(cx, slot)?; }
        if let Some(ThermalBc::Neumann { outward_flux }) = problem.boundary.condition_for(slot) {
            let mean = face.vertices.iter().map(|&v| outward_flux.at(v as usize)).sum::<f64>() / 3.0;
            heat = finite(heat + face.area * mean)?;
        }
    }
    poll(cx, problem.mesh.boundary().len())?;
    Ok(heat)
}

pub(crate) fn convection(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    fluxes: &[RobinFlux],
) -> Result<(Vec<RobinFlux>, f64), ConductionError> {
    let mut total = 0.0;
    let mut rows = Vec::with_capacity(fluxes.len());
    for (ordinal, row) in fluxes.iter().enumerate() {
        poll(cx, ordinal)?;
        let region = problem
            .boundary
            .region_names()
            .iter()
            .position(|name| *name == row.region)
            .ok_or_else(|| invalid("unknown original convection region"))?;
        let mut original = row.clone();
        match &problem.boundary.conditions()[region] {
            ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(reference),
            } => {
                original.mean_htc_w_per_m2_k = *h;
                original.mean_reference_temperature_k = *reference;
                original.heat_rate_w =
                    finite(h * row.area_m2 * (row.mean_wall_temperature_k - reference))?;
            }
            // Only radiating rows are changed. Unrelated nonuniform Robin
            // rows retain their exact inner-solve integrals.
            ThermalBc::Robin { .. } => {}
            _ => {
                return Err(invalid("enclosure report contains a non-Robin row"));
            }
        }
        total = finite(total + original.heat_rate_w)?;
        rows.push(original);
    }
    Ok((rows, total))
}

fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("nonfinite enclosure arithmetic")) }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "enclosure-exchange", at,
    })
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Radiation {
        surface: "enclosure-exchange".into(), what: what.into(),
        fix: "retain exact named traces, physical sourced emissivities and explicit heat/work allowances".into(),
    }
}
