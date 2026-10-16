//! Natural-convection derivatives at the accepted capacity-storage state.
//!
//! The existing card owns the signed wall and film-density physics. This
//! adapter binds its actual coefficient to the original Robin partition and
//! contracts independent ambient controls after the complete temperature dual.
use std::collections::BTreeSet;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};
use fs_conduction::adjoint::RobinGradient;
use crate::solve::natural::{NaturalLaw, differentiated_coefficient};
use super::{Cx, Instant, SolveRefusal, bad, finite, lower, poll};
use super::super::super::{json_string, number};

struct NaturalPort {
    target: String,
    port: usize,
    htc: f64,
    ambient_slope: f64,
}

/// Original convection with every natural coefficient evaluated at this state.
/// Slopes cover all caller-selected ports; prescribed rows have zero slope.
pub(super) struct Bound {
    pub(super) boundary: ThermalBoundary,
    pub(super) h_mean_slopes: Vec<f64>,
    nodes: usize,
    ports: usize,
    controls: Vec<NaturalPort>,
}

/// One active complete-feedback transpose plus natural-law/endpoint scratch.
/// Boundary and ordinary sparse physical assembly are admitted by the caller.
pub(super) fn storage_bytes(
    nodes: usize, ports: usize, natural_ports: usize,
) -> Result<usize, SolveRefusal> {
    if nodes == 0 || ports == 0 || ports > 64 || natural_ports == 0 || natural_ports > ports {
        return Err(bad("capacity natural adjoints require one to 64 distinct physical Robin ports"));
    }
    nodes.checked_mul(ports).and_then(|v| v.checked_mul(96))
        .and_then(|v| ports.checked_mul(ports).and_then(|p| p.checked_mul(48))
            .and_then(|p| v.checked_add(p)))
        .and_then(|v| nodes.checked_mul(128).and_then(|n| v.checked_add(n)))
        .and_then(|v| ports.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .and_then(|v| natural_ports.checked_mul(2048).and_then(|p| v.checked_add(p)))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| bad("capacity natural adjoint workspace overflows"))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn bind(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, accepted: &[f64],
    laws: &[NaturalLaw], pressure_pa: f64, regions: &[&str],
    deadline: Option<(Instant, f64)>,
) -> Result<Bound, SolveRefusal> {
    poll(cx, deadline)?;
    let nodes = problem.mesh.vertex_count();
    storage_bytes(nodes, regions.len(), laws.len())?;
    if accepted.len() != nodes || !pressure_pa.is_finite() || pressure_pa <= 0.0 {
        return Err(bad("capacity natural binding lost its accepted field or absolute pressure"));
    }
    for (vertex, &temperature) in accepted.iter().enumerate() {
        if vertex % 256 == 0 { poll(cx, deadline)?; }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(bad("capacity natural adjoints require positive accepted absolute temperatures"));
        }
    }
    let mut seen = BTreeSet::new();
    for &region in regions {
        if !seen.insert(region) {
            return Err(bad("capacity natural adjoint selected Robin controls repeat a region"));
        }
    }
    seen.clear();
    let mut selected_ports = Vec::with_capacity(laws.len());
    let mut boundary_regions = Vec::with_capacity(laws.len());
    for law in laws {
        poll(cx, deadline)?;
        if !seen.insert(law.target.as_str()) {
            return Err(bad("capacity natural adjoint laws repeat a boundary target"));
        }
        let port = regions.iter().position(|&name| name == law.target)
            .ok_or_else(|| bad("a natural law has no selected physical Robin control"))?;
        let region = problem.boundary.region_names().iter().position(|name| name == &law.target)
            .ok_or_else(|| bad("a natural law has no original physical boundary"))?;
        let Some(ThermalBc::Robin {
            htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
        }) = problem.boundary.conditions().get(region)
        else { return Err(bad("capacity natural adjoints require original uniform Robin rows")); };
        if !h.is_finite() || *h <= 0.0 || !reference.is_finite()
            || reference.to_bits() != law.ambient_k.to_bits()
        {
            return Err(bad("capacity natural adjoint lost its original coefficient or ambient identity"));
        }
        selected_ports.push(port);
        boundary_regions.push(region);
    }
    let mut areas = vec![0.0; laws.len()];
    let mut means = vec![0.0; laws.len()];
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 { poll(cx, deadline)?; }
        let Some(index) = problem.boundary.region_for(slot)
            .and_then(|region| boundary_regions.iter().position(|&found| found == region))
        else { continue; };
        // EXACT Robin report arithmetic, which the forward natural iteration
        // uses to query the card: ordered nodal sum, /3, fused area accumulate.
        let mut wall_sum = 0.0;
        for &vertex in &face.vertices {
            wall_sum = finite(wall_sum + accepted[vertex as usize])?;
        }
        areas[index] = finite(areas[index] + face.area)?;
        means[index] = finite(face.area.mul_add(wall_sum / 3.0, means[index]))?;
    }
    let mut replacements = Vec::with_capacity(laws.len());
    let mut h_mean_slopes = vec![0.0; regions.len()];
    let mut controls = Vec::with_capacity(laws.len());
    for (index, law) in laws.iter().enumerate() {
        poll(cx, deadline)?;
        if areas[index] <= 0.0 {
            return Err(bad("capacity natural adjoint has an empty exterior trace"));
        }
        let mean = finite(means[index] * (1.0 / areas[index]))?;
        // Both wall and ambient neighborhoods pass through the actual card;
        // equilibrium, unsupported cards and domain endpoints retain refusal.
        let (coefficient, partials) = differentiated_coefficient(law, mean, pressure_pa)?;
        let port = selected_ports[index];
        h_mean_slopes[port] = partials[0];
        replacements.push((boundary_regions[index], coefficient.htc_w_m2_k, law.ambient_k));
        controls.push(NaturalPort {
            target: law.target.clone(), port, htc: coefficient.htc_w_m2_k,
            ambient_slope: partials[1],
        });
    }
    let boundary = problem.boundary.with_uniform_robin_replacements(&replacements).map_err(lower)?;
    poll(cx, deadline)?;
    Ok(Bound { boundary, h_mean_slopes, nodes, ports: regions.len(), controls })
}

impl Bound {
    pub(super) fn owns_region(&self, name: &str) -> bool {
        self.controls.iter().any(|control| control.target == name)
    }

    /// The caller supplies ORIGINAL convection bars after any radiation remap.
    /// Wall feedback is already present in that complete temperature dual.
    /// An ambient control changes its explicit reference and film-dependent h.
    pub(super) fn ambient_partials(
        &self, cx: &Cx<'_>, convection: &RobinGradient, deadline: Option<(Instant, f64)>,
    ) -> Result<Vec<f64>, SolveRefusal> {
        poll(cx, deadline)?;
        if convection.references.len() != self.ports || convection.log_htc.len() != self.ports
            || convection.nodal_load.len() != self.nodes
            || !convection.relative_residual.is_finite() || convection.relative_residual < 0.0
        {
            return Err(bad("capacity natural controls need the complete original convection gradient"));
        }
        let mut values = Vec::with_capacity(self.controls.len());
        for control in &self.controls {
            poll(cx, deadline)?;
            let reference_bar = finite(convection.references[control.port])?;
            let coefficient_bar = finite(convection.log_htc[control.port])?;
            values.push(finite(reference_bar
                + coefficient_bar * finite(control.ambient_slope / control.htc)?)?);
        }
        poll(cx, deadline)?;
        Ok(values)
    }
}

/// Natural targets remain independent controls even when their cards match.
/// Query pressure and transport properties remain fixed throughout the history.
pub(super) fn fragment(
    laws: &[NaturalLaw], derivatives: &[f64], step_count: usize,
) -> Result<String, SolveRefusal> {
    if laws.is_empty() || laws.len() > 64 || derivatives.len() != laws.len() {
        return Err(bad("capacity natural report lost its physical reservoir controls"));
    }
    let mut seen = BTreeSet::new();
    let mut rows = Vec::with_capacity(laws.len());
    for (law, &derivative) in laws.iter().zip(derivatives) {
        if !seen.insert(law.target.as_str()) {
            return Err(bad("capacity natural report repeats a reservoir target"));
        }
        rows.push(format!(
            "{{\"target\":{},\"ambient_temperature_k\":{},\"dtemperature_dambient_temperature_k_per_k\":{}}}",
            json_string(&law.target), number(law.ambient_k)?, number(derivative)?,
        ));
    }
    Ok(format!(
        ",\"natural_control_steps\":{step_count},\"natural_reservoirs\":[{}]",
        rows.join(","),
    ))
}
