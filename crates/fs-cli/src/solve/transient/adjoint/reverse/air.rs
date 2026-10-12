//! Actual streamwise air feedback at an accepted capacity-storage endpoint.
//!
//! The physical path march supplies the Robin references and air heat gain.
//! Its analytic affine wall map supplies every upstream reference derivative;
//! the differently rounded affine replay is never used for endpoint admission.

use std::collections::BTreeSet;

use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::affine_reference_law};
use fs_conduction::transient::backward_euler::EndpointEvaluation;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};

use super::{Cx, Instant, SolveRefusal, bad, finite, lower, poll};
use crate::solve::conjugate::ConjugatePath;

/// Original convection evaluated at this field, before optional radiation.
/// The square feedback is in the caller's complete Robin-control order.
pub(super) struct Bound {
    pub(super) boundary: ThermalBoundary,
    pub(super) reference_feedback: Vec<f64>,
    air_regions: BTreeSet<String>,
    air_heat_w: f64,
}

/// One active complete-feedback transpose and physical air reconstruction.
/// The caller separately admits ordinary boundary and sparse assembly storage.
pub(super) fn storage_bytes(
    nodes: usize,
    ports: usize,
    air_ports: usize,
) -> Result<usize, SolveRefusal> {
    if nodes == 0 || ports == 0 || ports > 64 || air_ports == 0 || air_ports > ports {
        return Err(bad(
            "capacity air adjoints require one to 64 distinct physical Robin ports",
        ));
    }
    nodes
        .checked_mul(ports)
        .and_then(|v| v.checked_mul(96))
        .and_then(|v| {
            ports
                .checked_mul(ports)
                .and_then(|p| p.checked_mul(48))
                .and_then(|p| v.checked_add(p))
        })
        .and_then(|v| nodes.checked_mul(128).and_then(|n| v.checked_add(n)))
        .and_then(|v| ports.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .and_then(|v| air_ports.checked_mul(8192).and_then(|p| v.checked_add(p)))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| bad("capacity air adjoint workspace overflows"))
}

fn air_error(error: CoupledGoalError) -> SolveRefusal {
    match error {
        CoupledGoalError::Interrupted => crate::solve::conduction_error(
            "cli-solve-cancelled",
            "capacity air-reference differentiation was cancelled",
            "resume the accepted pipeline prefix",
        ),
        other => bad(format!(
            "capacity air-reference derivative refused: {other}"
        )),
    }
}

/// Keep any already evaluated natural rows, replacing only forced-air rows.
/// Natural and air ownership is disjoint in the admitted forward problem.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn bind(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    accepted: &[f64],
    lowering: &ConjugatePath,
    regions: &[&str],
    deadline: Option<(Instant, f64)>,
) -> Result<Bound, SolveRefusal> {
    poll(cx, deadline)?;
    let nodes = problem.mesh.vertex_count();
    let count = lowering.segments.len();
    storage_bytes(nodes, regions.len(), count)?;
    if accepted.len() != nodes {
        return Err(bad("capacity air binding lost its accepted field layout"));
    }
    for (vertex, &temperature) in accepted.iter().enumerate() {
        if vertex % 256 == 0 {
            poll(cx, deadline)?;
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(bad(
                "capacity air adjoints require positive accepted absolute temperatures",
            ));
        }
    }
    let mut selected = BTreeSet::new();
    for &region in regions {
        if !selected.insert(region) {
            return Err(bad(
                "capacity air adjoint selected Robin controls repeat a region",
            ));
        }
    }
    let paths = lowering.air_paths();
    let mut names = Vec::with_capacity(count);
    let mut ports = Vec::with_capacity(count);
    let mut boundary_regions = Vec::with_capacity(count);
    let mut air_regions = BTreeSet::new();
    for path in &paths {
        poll(cx, deadline)?;
        for segment in path.segments() {
            let ordinal = names.len();
            let derivation = lowering
                .segments
                .get(ordinal)
                .ok_or_else(|| bad("capacity air path has extra physical segments"))?;
            let name = segment.region();
            if name != derivation.target.as_str()
                || segment.htc_w_per_m2_k().to_bits() != derivation.htc_w_m2_k.to_bits()
                || segment.area_m2().to_bits() != derivation.wetted_area_m2.to_bits()
                || !air_regions.insert(name.to_string())
            {
                return Err(bad(
                    "capacity air path and original segment identities differ",
                ));
            }
            let port = regions
                .iter()
                .position(|&found| found == name)
                .ok_or_else(|| bad("capacity air segment has no selected Robin control"))?;
            let region = problem
                .boundary
                .region_names()
                .iter()
                .position(|found| found == name)
                .ok_or_else(|| bad("capacity air segment has no original physical boundary"))?;
            let Some(ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(reference),
            }) = problem.boundary.conditions().get(region)
            else {
                return Err(bad(
                    "capacity air adjoints require original uniform convection",
                ));
            };
            if !h.is_finite()
                || *h <= 0.0
                || !reference.is_finite()
                || *reference <= 0.0
                || h.to_bits() != segment.htc_w_per_m2_k().to_bits()
            {
                return Err(bad(
                    "capacity air and solid must share the same original positive coefficient",
                ));
            }
            names.push(name.to_string());
            ports.push(port);
            boundary_regions.push(region);
        }
    }
    if names.len() != count || paths.is_empty() {
        return Err(bad("capacity air branch ordering omits physical segments"));
    }
    let mut areas = vec![0.0; count];
    let mut means = vec![0.0; count];
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % 256 == 0 {
            poll(cx, deadline)?;
        }
        let Some(index) = problem
            .boundary
            .region_for(slot)
            .and_then(|region| boundary_regions.iter().position(|&found| found == region))
        else {
            continue;
        };
        // EXACT Robin report arithmetic used by the forward air exchange:
        // ordered nodal sum, division by three, fused area accumulation and
        // final multiplication by reciprocal total area.
        let mut wall_sum = 0.0;
        for &vertex in &face.vertices {
            wall_sum = finite(wall_sum + accepted[vertex as usize])?;
        }
        areas[index] = finite(areas[index] + face.area)?;
        means[index] = finite(face.area.mul_add(wall_sum / 3.0, means[index]))?;
    }
    for (index, mean) in means.iter_mut().enumerate() {
        if areas[index] <= 0.0 {
            return Err(bad(
                "capacity air derivative has an empty physical wall trace",
            ));
        }
        *mean = finite(*mean * (1.0 / areas[index]))?;
    }
    let law = affine_reference_law(cx, &paths, 64, count * count).map_err(air_error)?;
    if law.regions() != names.as_slice() || law.wall_matrix().len() != count * count {
        return Err(bad(
            "capacity air derivative lost its branch-major wall coordinates",
        ));
    }
    let mut reference_feedback = vec![0.0; regions.len() * regions.len()];
    for (row, &port) in ports.iter().enumerate() {
        poll(cx, deadline)?;
        for (column, &source) in ports.iter().enumerate() {
            reference_feedback[port * regions.len() + source] =
                finite(law.wall_matrix()[row * count + column])?;
        }
    }
    let mut replacements = Vec::with_capacity(count);
    let mut start = 0;
    let mut air_heat_w = 0.0;
    for path in &paths {
        poll(cx, deadline)?;
        let end = start + path.segments().len();
        // The physical march deliberately remains separate from affine
        // coefficient replay, so revalidation uses the accepted arithmetic.
        let march = path
            .march(&means[start..end])
            .map_err(|error| bad(format!("accepted capacity air march refused: {error}")))?;
        if march.segments.len() != path.segments().len() {
            return Err(bad("capacity air march lost its physical segment layout"));
        }
        for (local, (segment, state)) in path.segments().iter().zip(&march.segments).enumerate() {
            let index = start + local;
            let area = segment.area_m2();
            if state.region != names[index]
                || !area.is_finite()
                || area <= 0.0
                || (areas[index] - area).abs() > 128.0 * f64::EPSILON * areas[index].max(area)
                || !state.reference_temperature_k.is_finite()
                || state.reference_temperature_k <= 0.0
            {
                return Err(bad(
                    "capacity air and solid must share names, wetted area and physical references",
                ));
            }
            replacements.push((
                boundary_regions[index],
                segment.htc_w_per_m2_k(),
                state.reference_temperature_k,
            ));
        }
        air_heat_w = finite(air_heat_w + march.total_heat_rate_w)?;
        start = end;
    }
    let boundary = problem
        .boundary
        .with_uniform_robin_replacements(&replacements)
        .map_err(lower)?;
    poll(cx, deadline)?;
    Ok(Bound {
        boundary,
        reference_feedback,
        air_regions,
        air_heat_w,
    })
}

impl Bound {
    pub(super) fn owns_region(&self, name: &str) -> bool {
        self.air_regions.contains(name)
    }

    /// Replace only on-path convection by actual air enthalpy gain. Independent
    /// natural reservoirs, gray radiation and prescribed-node capacity reaction
    /// retain their own terms. The caller separately gates the full free nodal
    /// residual against the ORIGINAL accepted forward tolerance.
    pub(super) fn check_air_energy(
        &self,
        physical: &EndpointEvaluation,
        dt_s: f64,
        tolerance_j: f64,
    ) -> Result<(), SolveRefusal> {
        if !dt_s.is_finite() || dt_s <= 0.0 || !tolerance_j.is_finite() || tolerance_j <= 0.0 {
            return Err(bad(
                "capacity air energy audit needs positive time and joule tolerance",
            ));
        }
        let mut off_path_w = 0.0;
        let mut matched = BTreeSet::new();
        for flux in &physical.convective_robin_fluxes {
            if self.air_regions.contains(&flux.region) {
                if !matched.insert(flux.region.as_str()) {
                    return Err(bad("capacity adjoint physical air audit repeats a port"));
                }
            } else {
                off_path_w = finite(off_path_w + flux.heat_rate_w)?;
            }
        }
        if matched.len() != self.air_regions.len() {
            return Err(bad("capacity adjoint physical air audit lost a port"));
        }
        let net_w = finite(
            physical.source_w + physical.dirichlet_in_w
                - physical.neumann_out_w
                - off_path_w
                - self.air_heat_w
                - physical.radiation_out_w,
        )?;
        let residual_j = finite(physical.stored_energy_change_j - dt_s * net_w)?;
        if residual_j.abs() > tolerance_j {
            return Err(bad(format!(
                "accepted capacity adjoint actual-air energy residual {residual_j:e} J exceeds {tolerance_j:e} J",
            )));
        }
        Ok(())
    }
}
