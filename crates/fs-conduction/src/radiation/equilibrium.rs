//! Reassemble the coupled discrete equations at the proposed final field.
//! A temperature update is not a substitute for the actual radiative load.
mod state;

use super::{
    ConductionError, ConductionProblem, ConductionSolution, Cx, DofMap,
    RadiationSurface, RadiosityReport, ThermalInterfaces,
    radiation_error,
};
use crate::assemble::{assemble_operator_scaled_with_interfaces, residual};

/// Final-state evidence for the actual nonlinear enclosure equations.
///
/// The field has not been advanced or corrected here. Radiation is evaluated
/// at its actual surface means, and the original conduction/contact assembly
/// is rerun at that SAME field with those fluxes. Only free equilibrium rows
/// enter the norm; prescribed rows are reactions, not unresolved equations.
/// This is a discrete residual check, not a temperature-error certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct EnclosureEquilibriumReport {
    /// Radiation recomputed at the returned conduction field, not its driver.
    pub nonlinear_radiosity: RadiosityReport,
    /// Euclidean norm of the reassembled free solid residual, W.
    pub residual_w: f64,
    /// The original final inner solve's threshold, retained without rescaling.
    pub residual_threshold_w: f64,
    /// Largest absolute difference between applied and recomputed patch heat.
    /// This diagnostic includes constrained trace vertices; it is not an
    /// independent per-patch tolerance or a continuum heat-error bound.
    pub max_patch_heat_mismatch_w: f64,
}

fn poll(cx: &Cx<'_>) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "radiation-equilibrium-residual", at: 0,
    })
}

/// Called only with fresh radiosity at the actual candidate's temperatures.
/// The driver reuses its factorization, never a previous temperature RHS.
#[allow(clippy::too_many_arguments)]
pub(super) fn check<F>(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    interfaces: Option<&ThermalInterfaces>,
    candidate: &ConductionSolution,
    nonlinear: RadiosityReport,
    applied: &RadiosityReport,
    apply_boundary: &F,
    original_convection: Option<&[RadiationSurface]>,
) -> Result<EnclosureEquilibriumReport, ConductionError>
where
    F: Fn(&RadiosityReport) -> Result<crate::ThermalBoundary, ConductionError>,
{
    poll(cx)?;
    let mut max_patch_heat_mismatch_w = 0.0f64;
    for (&new, &old) in nonlinear.net_outward_heat_w.iter().zip(&applied.net_outward_heat_w) {
        let mismatch = (new - old).abs();
        if !mismatch.is_finite() {
            return Err(radiation_error("enclosure-equilibrium", "nonfinite patch heat mismatch",
                "retain finite physical radiation and conduction inputs"));
        }
        max_patch_heat_mismatch_w = max_patch_heat_mismatch_w.max(mismatch);
    }
    if nonlinear.radiosity_w_m2.iter().chain(&nonlinear.irradiation_w_m2)
        .chain(&nonlinear.net_outward_flux_w_m2).chain(&nonlinear.net_outward_heat_w)
        .any(|v| !v.is_finite())
        || !nonlinear.enclosure_energy_closure_w.is_finite()
        || !nonlinear.enclosure_energy_scale_w.is_finite()
        || !nonlinear.linear_residual_max_w_m2.is_finite()
    {
        return Err(radiation_error("enclosure-equilibrium", "nonfinite nonlinear radiation report",
            "retain finite physical radiation and conduction inputs"));
    }
    let boundary = match original_convection {
        Some(_) => problem.boundary.clone(),
        None => apply_boundary(&nonlinear)?,
    };
    let dofs = DofMap::new(&boundary, problem.mesh.vertex_count())?;
    // Crucially retain the original heterogeneous k(T), contact and source.
    // The matrix is assembled at T, not at the preceding radiation driver.
    let system = assemble_operator_scaled_with_interfaces(cx, problem.mesh, &boundary,
        problem.material, problem.source, &candidate.temperature, None, interfaces,
        problem.element_materials)?;
    poll(cx)?;
    let mut values = residual(&system, &dofs, &candidate.temperature);
    if let Some(surfaces) = original_convection {
        // Reassemble the ORIGINAL convection and add actual radiative loads
        // directly. Reconstructing a shifted Robin reference here could hide
        // radiation when q/h rounds away against the reservoir temperature.
        // Uniform outward flux contributes +area*q/3 to every face vertex.
        for (surface, &flux) in surfaces.iter().zip(&nonlinear.net_outward_flux_w_m2) {
            for (ordinal, &slot) in surface.face_slots.iter().enumerate() {
                if ordinal % crate::assemble::ASSEMBLY_TILE == 0 { poll(cx)?; }
                let face = &problem.mesh.boundary()[slot];
                let nodal_heat = face.area * (flux / 3.0);
                if !nodal_heat.is_finite() {
                    return Err(radiation_error(surface.name(), "nonfinite radiative nodal load",
                        "use supported finite heat flux and surface geometry"));
                }
                for &vertex in &face.vertices {
                    if let Some(index) = dofs.slot_of(vertex as usize) {
                        values[index] += nodal_heat;
                    }
                }
            }
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(radiation_error("enclosure-equilibrium", "nonfinite reassembled residual",
            "inspect the original physical model and numerical range"));
    }
    let residual_w = fs_solver::norm2(&values);
    let threshold = candidate.report.residual_threshold;
    if !residual_w.is_finite() || !threshold.is_finite() || threshold < 0.0 {
        return Err(radiation_error("enclosure-equilibrium", "invalid nonlinear residual or threshold",
            "retain the original finite nonnegative solid residual policy"));
    }
    poll(cx)?;
    Ok(EnclosureEquilibriumReport {
        nonlinear_radiosity: nonlinear, residual_w, residual_threshold_w: threshold,
        max_patch_heat_mismatch_w,
    })
}
