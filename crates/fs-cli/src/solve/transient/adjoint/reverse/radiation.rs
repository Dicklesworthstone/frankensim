//! Bind native capacity-storage radiation to its accepted physical endpoint.
//!
//! Radiation augments the original convection row. The complete mean-state
//! feedback belongs to the temperature transpose; parameter bars are then
//! mapped back to the two independent reservoirs and the sourced emissivity.
use std::collections::BTreeSet;
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};
use fs_conduction::adjoint::RobinGradient;
use fs_conduction::radiation::AmbientRadiationGradient;
use fs_conduction::transient::backward_euler::{EndpointEvaluation, StepConfig};
use crate::solve::radiation::LoweredRadiation;
use super::{Cx, Instant, SolveRefusal, bad, finite, lower, poll};
use super::super::super::{json_string, number};

#[derive(Clone, Copy)]
struct PhysicalPort {
    convection_h: f64,
    convection_reference: f64,
    combined_h: f64,
    combined_reference: f64,
}

struct RadiativePort {
    port: usize,
    region: String,
    secant: f64,
    reservoir_partial: f64,
    reservoir: f64,
    emissivity: f64,
}

/// One independently audited endpoint; no forward or transpose solve is run here.
pub(super) struct Bound {
    pub(super) boundary: ThermalBoundary,
    pub(super) h_mean_slopes: Vec<f64>,
    pub(super) local_reference_slopes: Vec<f64>,
    pub(super) reference_feedback: Vec<f64>,
    pub(super) max_feedback_entries: usize,
    nodes: usize,
    dual_tolerance: f64,
    regions: Vec<String>,
    physical: Vec<PhysicalPort>,
    radiative: Vec<RadiativePort>,
}

/// Additional mean-feedback factors, endpoint scratch and parameter metadata.
/// The native caller separately admits the ordinary physical assembly and
/// boundary storage. No factor or matrix allocation precedes this checked cap.
pub(in crate::solve::transient) fn storage_bytes(
    nodes: usize, ports: usize, patches: usize,
) -> Result<u64, SolveRefusal> {
    if nodes == 0 || ports == 0 || ports > 64 || patches == 0 || patches > ports {
        return Err(bad("capacity radiation adjoint requires one to 64 bound Robin ports and patches"));
    }
    let nodes = u64::try_from(nodes).map_err(|_| bad("capacity radiation node count overflows"))?;
    let ports = ports as u64;
    let patches = patches as u64;
    nodes.checked_mul(ports).and_then(|v| v.checked_mul(96))
        .and_then(|v| ports.checked_mul(ports).and_then(|m| m.checked_mul(48))
            .and_then(|m| v.checked_add(m)))
        .and_then(|v| nodes.checked_mul(128).and_then(|n| v.checked_add(n)))
        .and_then(|v| ports.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .and_then(|v| patches.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .ok_or_else(|| bad("capacity radiation feedback workspace overflows"))
}

/// Enclosure exchange has its own patch basis, including Neumann-only walls.
/// Reserve the Robin workspace plus the independent node-by-patch factors,
/// their assembly scratch and dense patch-space radiosity/transpose work.
/// A closed vacuum enclosure therefore remains admissible with zero ports.
pub(in crate::solve::transient) fn enclosure_storage_bytes(
    nodes: usize, ports: usize, patches: usize,
) -> Result<u64, SolveRefusal> {
    if nodes == 0 || ports > 64 || patches == 0 || patches > 64 {
        return Err(bad("capacity enclosure adjoint requires one to 64 patches and at most 64 Robin ports"));
    }
    let nodes = u64::try_from(nodes).map_err(|_| bad("capacity enclosure node count overflows"))?;
    let ports = ports as u64;
    let patches = patches as u64;
    nodes.checked_mul(ports).and_then(|v| v.checked_mul(96))
        .and_then(|v| ports.checked_mul(ports).and_then(|m| m.checked_mul(48))
            .and_then(|m| v.checked_add(m)))
        .and_then(|v| nodes.checked_mul(patches).and_then(|p| p.checked_mul(64))
            .and_then(|p| v.checked_add(p)))
        .and_then(|v| patches.checked_mul(patches).and_then(|p| p.checked_mul(64))
            .and_then(|p| v.checked_add(p)))
        .and_then(|v| nodes.checked_mul(128).and_then(|n| v.checked_add(n)))
        .and_then(|v| ports.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .and_then(|v| patches.checked_mul(4096).and_then(|p| v.checked_add(p)))
        .ok_or_else(|| bad("capacity enclosure feedback workspace overflows"))
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn bind(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, accepted: &[f64],
    config: StepConfig, lowering: &LoweredRadiation, regions: &[&str],
    convection_mean_slopes: Option<&[f64]>, convection_reference_feedback: Option<&[f64]>,
    audit: EndpointEvaluation, deadline: Option<(Instant, f64)>,
) -> Result<Bound, SolveRefusal> {
    poll(cx, deadline)?;
    let n = problem.mesh.vertex_count();
    storage_bytes(n, regions.len(), lowering.patches.len())?;
    if accepted.len() != n
        || lowering.surface_labels().len() != lowering.patches.len()
    {
        return Err(bad("capacity radiation binding lost its field, source or physical target"));
    }
    if let Some(slopes) = convection_mean_slopes {
        if slopes.len() != regions.len() {
            return Err(bad("capacity radiation lost its natural coefficient slope layout"));
        }
        for &slope in slopes { finite(slope)?; }
    }
    if let Some(feedback) = convection_reference_feedback {
        if feedback.len() != regions.len() * regions.len() {
            return Err(bad("capacity radiation lost its air-reference feedback layout"));
        }
        for &value in feedback { finite(value)?; }
    }
    let mut reference_feedback = convection_reference_feedback.map_or_else(
        || vec![0.0; regions.len() * regions.len()], |feedback| feedback.to_vec(),
    );
    let mut seen = BTreeSet::new();
    for &region in regions {
        if !seen.insert(region) {
            return Err(bad("capacity radiation controls repeat a Robin region"));
        }
    }
    seen.clear();
    for ((name, target), patch) in lowering.surface_labels().zip(&lowering.patches) {
        poll(cx, deadline)?;
        if name.trim().is_empty() || target != patch.region()
            || !regions.contains(&target) || !seen.insert(target)
        {
            return Err(bad("capacity radiation surface identities differ from bound Robin ports"));
        }
    }
    for (vertex, &temperature) in accepted.iter().enumerate() {
        if vertex % 256 == 0 { poll(cx, deadline)?; }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(bad("capacity radiation needs positive accepted absolute temperatures"));
        }
    }
    // The shared reverse owner audited the actual gray heat, full free nodal
    // balance and any transported air energy at this unchanged endpoint.
    let boundary = audit.combined_boundary;
    if boundary.region_names() != problem.boundary.region_names() {
        return Err(bad("capacity radiation changed physical boundary identities"));
    }
    let mut physical = Vec::with_capacity(regions.len());
    let mut h_mean_slopes = convection_mean_slopes.map_or_else(
        || vec![0.0; regions.len()], |slopes| slopes.to_vec(),
    );
    let mut local_reference_slopes = vec![0.0; regions.len()];
    for &name in regions {
        poll(cx, deadline)?;
        let region = problem.boundary.region_names().iter().position(|value| value == name)
            .ok_or_else(|| bad("capacity radiation has an unknown control region"))?;
        let (convection_h, convection_reference) = uniform_row(problem.boundary, region)?;
        let (combined_h, combined_reference) = uniform_row(&boundary, region)?;
        if convection_h < 0.0 || combined_h <= 0.0 {
            return Err(bad("capacity radiation requires nonnegative convection and positive combined coefficients"));
        }
        physical.push(PhysicalPort {
            convection_h, convection_reference, combined_h, combined_reference,
        });
    }
    let mut radiative = Vec::with_capacity(lowering.patches.len());
    for patch in &lowering.patches {
        poll(cx, deadline)?;
        let port = regions.iter().position(|&name| name == patch.region())
            .ok_or_else(|| bad("capacity radiation patch lost its selected control"))?;
        let region = problem.boundary.region_names().iter().position(|name| name == patch.region())
            .ok_or_else(|| bad("capacity radiation patch lost its physical boundary"))?;
        let mut area = 0.0;
        let mut integral = 0.0;
        for (slot, face) in problem.mesh.boundary().iter().enumerate() {
            if slot % 256 == 0 { poll(cx, deadline)?; }
            if boundary.region_for(slot) != problem.boundary.region_for(slot) {
                return Err(bad("capacity radiation changed exterior face ownership"));
            }
            if problem.boundary.region_for(slot) != Some(region) { continue; }
            area = finite(area + face.area)?;
            // EXACT capacity-forward arithmetic: divide each vertex by three,
            // sum the face, multiply by area, and finally divide by total area.
            let mean = face.vertices.iter()
                .map(|&vertex| accepted[vertex as usize] / 3.0).sum::<f64>();
            integral = finite(integral + face.area * mean)?;
        }
        if area <= 0.0 {
            return Err(bad("capacity radiation patch has no positive exterior area"));
        }
        let mean = finite(integral / area)?;
        let secant = patch.secant_coefficient_w_m2_k(mean).map_err(lower)?;
        let [wall_partial, reservoir_partial] = patch.secant_partials_w_m2_k2(mean).map_err(lower)?;
        let reservoir = patch.ambient_temperature_k();
        if reservoir.next_down() <= 0.0 || !reservoir.next_up().is_finite() {
            return Err(bad("capacity radiation reservoir lacks a positive two-sided temperature neighborhood"));
        }
        let point = physical[port];
        let combined_h = finite(point.convection_h + secant)?;
        let combined_reference = finite(
            (point.convection_h / combined_h) * point.convection_reference
                + (secant / combined_h) * reservoir,
        )?;
        if combined_h.to_bits() != point.combined_h.to_bits()
            || combined_reference.to_bits() != point.combined_reference.to_bits()
        {
            return Err(bad("capacity radiation derivative differs from the actual accepted boundary"));
        }
        // The combined reference is hc/H times the original air reference.
        // Scale the whole row, including upstream walls, exactly once.
        for source in 0..regions.len() {
            let entry = port * regions.len() + source;
            reference_feedback[entry] = finite(
                reference_feedback[entry] * (point.convection_h / combined_h),
            )?;
        }
        let convection_wall_partial = h_mean_slopes[port];
        h_mean_slopes[port] = finite(convection_wall_partial + wall_partial)?;
        local_reference_slopes[port] = finite(
            (wall_partial / combined_h) * (reservoir - combined_reference),
        )?;
        if convection_wall_partial != 0.0 {
            local_reference_slopes[port] = finite(local_reference_slopes[port]
                + (convection_wall_partial / combined_h)
                    * (point.convection_reference - combined_reference))?;
        }
        radiative.push(RadiativePort {
            port, region: patch.region().to_owned(), secant, reservoir_partial,
            reservoir, emissivity: patch.emissivity().value(),
        });
    }
    let max_feedback_entries = n.checked_mul(regions.len()).and_then(|v| v.checked_mul(4))
        .ok_or_else(|| bad("capacity radiation feedback entry count overflows"))?;
    poll(cx, deadline)?;
    Ok(Bound {
        boundary, h_mean_slopes, local_reference_slopes,
        reference_feedback,
        max_feedback_entries, nodes: n, dual_tolerance: config.linear.tolerance,
        regions: regions.iter().map(|name| (*name).to_owned()).collect(), physical, radiative,
    })
}

fn uniform_row(boundary: &ThermalBoundary, region: usize) -> Result<(f64, f64), SolveRefusal> {
    let Some(ThermalBc::Robin {
        htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
    }) = boundary.conditions().get(region)
    else { return Err(bad("capacity radiation controls require original uniform Robin rows")); };
    Ok((finite(*h)?, finite(*reference)?))
}

impl Bound {
    /// Map one COMPLETE residual-checked Robin gradient to physical controls.
    /// The watt-coordinate nodal load multiplier passes through unchanged;
    /// neither this chain nor the radiation coefficient receives another dt.
    pub(super) fn remap_gradient(
        &self, cx: &Cx<'_>, mut gradient: RobinGradient, deadline: Option<(Instant, f64)>,
    ) -> Result<AmbientRadiationGradient, SolveRefusal> {
        poll(cx, deadline)?;
        if gradient.references.len() != self.regions.len()
            || gradient.log_htc.len() != self.regions.len()
            || gradient.nodal_load.len() != self.nodes
            || !gradient.relative_residual.is_finite()
            || gradient.relative_residual >= self.dual_tolerance
        {
            return Err(bad("capacity radiation parameter chain needs the complete checked Robin gradient"));
        }
        for (index, &value) in gradient.references.iter().chain(&gradient.log_htc)
            .chain(&gradient.nodal_load).enumerate()
        {
            if index % 256 == 0 { poll(cx, deadline)?; }
            finite(value)?;
        }
        let mut emissivities = Vec::with_capacity(self.radiative.len());
        let mut reservoir_temperatures = Vec::with_capacity(self.radiative.len());
        for patch in &self.radiative {
            poll(cx, deadline)?;
            let port = patch.port;
            let values = parameter_partials(
                self.physical[port], patch, gradient.references[port], gradient.log_htc[port],
            )?;
            gradient.references[port] = values[0];
            gradient.log_htc[port] = values[1];
            emissivities.push(values[2]);
            reservoir_temperatures.push(values[3]);
        }
        poll(cx, deadline)?;
        Ok(AmbientRadiationGradient {
            regions: self.regions.clone(), convection: gradient,
            radiation_regions: self.radiative.iter().map(|patch| patch.region.clone()).collect(),
            reservoir_temperatures, emissivities,
        })
    }
}

fn parameter_partials(
    point: PhysicalPort, patch: &RadiativePort, reference_bar: f64, log_coefficient_bar: f64,
) -> Result<[f64; 4], SolveRefusal> {
    let coefficient_bar = finite(
        log_coefficient_bar + reference_bar * (patch.reservoir - point.combined_reference),
    )?;
    Ok([
        finite(reference_bar * (point.convection_h / point.combined_h))?,
        finite((point.convection_h / point.combined_h)
            * (log_coefficient_bar + reference_bar
                * (point.convection_reference - point.combined_reference)))?,
        finite(coefficient_bar * (patch.secant / point.combined_h) / patch.emissivity)?,
        finite(reference_bar * (patch.secant / point.combined_h)
            + coefficient_bar * (patch.reservoir_partial / point.combined_h))?,
    ])
}

/// Immutable source identities stay distinct even when surfaces share a card.
pub(super) fn fragment(
    lowering: &LoweredRadiation, emissivities: &[f64], reservoirs: &[f64], step_count: usize,
) -> Result<String, SolveRefusal> {
    let count = lowering.patches.len();
    if count == 0 || count > 64 || emissivities.len() != count || reservoirs.len() != count
        || lowering.surface_labels().len() != count
    {
        return Err(bad("capacity radiation report lost its complete surface-control layout"));
    }
    let mut rows = Vec::with_capacity(count);
    for (index, ((name, target), patch)) in lowering.surface_labels()
        .zip(&lowering.patches).enumerate()
    {
        if target != patch.region() {
            return Err(bad("capacity radiation report changed a surface identity"));
        }
        let emissivity = patch.emissivity();
        let directions = if emissivity.value() == 1.0 {
            "nonpositive-at-emissivity-one"
        } else {
            "two-sided-within-zero-exclusive-one-inclusive"
        };
        rows.push(format!(
            "{{\"name\":{},\"target\":{},\"card\":{},\"material_state\":{},\"emissivity_receipt\":{},\"emissivity\":{},\"query_temperature_k\":{},\"ambient_temperature_k\":{},\"dtemperature_demissivity_k\":{},\"dtemperature_dreservoir_temperature_k_per_k\":{},\"emissivity_feasible_directions\":{}}}",
            json_string(name), json_string(target), json_string(&emissivity.card_identity().to_hex()),
            json_string(emissivity.material_state()),
            json_string(&emissivity.receipt().content_hash().to_hex()),
            number(emissivity.value())?, number(emissivity.temperature_k())?,
            number(patch.ambient_temperature_k())?, number(emissivities[index])?,
            number(reservoirs[index])?, json_string(directions),
        ));
    }
    Ok(format!(
        ",\"radiation\":{{\"method\":\"complete-area-mean-gray-capacity-tangent\",\"patch_count\":{count},\"control_steps\":{step_count},\"parameters\":\"absolute-per-surface-emissivity-and-reservoir-temperature-throughout-trajectory\",\"residual_target\":\"original-accepted-forward-physical-tolerance\",\"surfaces\":[{}]}}",
        rows.join(","),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enclosure_workspace_admits_independent_flux_patches_without_relaxing_ambient() {
        let vacuum = enclosure_storage_bytes(32, 0, 6).unwrap();
        let mixed = enclosure_storage_bytes(32, 1, 6).unwrap();
        let all_robin = enclosure_storage_bytes(32, 6, 6).unwrap();
        assert!(vacuum > 0 && mixed > vacuum && all_robin > mixed);
        // The ambient secant is carried by an existing Robin row; it cannot
        // use the enclosure's separate mean-flux factors.
        assert!(storage_bytes(32, 0, 6).is_err());
        assert!(storage_bytes(32, 1, 6).is_err());
        assert!(storage_bytes(32, 6, 6).is_ok());
        assert!(all_robin > storage_bytes(32, 6, 6).unwrap());
        for (ports, patches) in [(0, 1), (0, 64), (1, 6), (64, 1), (64, 64)] {
            let bytes = enclosure_storage_bytes(32, ports, patches).unwrap();
            let retained_factor_bytes = 8 * (4 * 32 * ports + 2 * 32 * patches) as u64;
            assert!(bytes >= retained_factor_bytes);
        }
    }

    #[test]
    fn enclosure_workspace_refuses_unbounded_layouts_and_checked_overflow() {
        for (nodes, ports, patches) in [(0, 0, 1), (32, 65, 1), (32, 0, 0), (32, 0, 65)] {
            assert!(enclosure_storage_bytes(nodes, ports, patches).is_err());
        }
        if usize::BITS == 64 {
            assert!(enclosure_storage_bytes(usize::MAX, 0, 1).is_err());
            assert!(enclosure_storage_bytes(usize::MAX, 1, 64).is_err());
        }
    }

    fn point(wall: f64, reservoir: f64, emissivity: f64, convection_h: f64, reference: f64)
        -> (PhysicalPort, RadiativePort)
    {
        let sigma = 5.670_374_419e-8;
        let secant = emissivity * sigma * (wall + reservoir)
            * wall.mul_add(wall, reservoir * reservoir);
        let reservoir_partial = emissivity * sigma
            * (wall * wall + 2.0 * wall * reservoir + 3.0 * reservoir * reservoir);
        let combined_h = convection_h + secant;
        (PhysicalPort {
            convection_h, convection_reference: reference, combined_h,
            combined_reference: (convection_h / combined_h) * reference
                + (secant / combined_h) * reservoir,
        }, RadiativePort {
            port: 0, region: "wall".into(), secant, reservoir_partial, reservoir, emissivity,
        })
    }

    #[test]
    fn physical_control_chain_matches_fresh_combined_rows_for_both_heat_signs() {
        let reference_bar = -0.4;
        let log_bar = 2.5;
        for reservoir in [270.0, 300.0, 350.0] {
            let wall = 300.0;
            let emissivity = 0.72;
            let h = 12.3;
            let reference = 293.15;
            let (physical, patch) = point(wall, reservoir, emissivity, h, reference);
            let actual = parameter_partials(physical, &patch, reference_bar, log_bar).unwrap();
            let objective = |ambient, epsilon, coefficient, reference| {
                let (physical, _) = point(wall, ambient, epsilon, coefficient, reference);
                reference_bar * physical.combined_reference + log_bar * physical.combined_h.ln()
            };
            let delta = 1.0e-4;
            let expected = [
                (objective(reservoir, emissivity, h, reference + delta)
                    - objective(reservoir, emissivity, h, reference - delta)) / (2.0 * delta),
                (objective(reservoir, emissivity, h * delta.exp(), reference)
                    - objective(reservoir, emissivity, h * (-delta).exp(), reference)) / (2.0 * delta),
                (objective(reservoir, emissivity + delta, h, reference)
                    - objective(reservoir, emissivity - delta, h, reference)) / (2.0 * delta),
                (objective(reservoir + delta, emissivity, h, reference)
                    - objective(reservoir - delta, emissivity, h, reference)) / (2.0 * delta),
            ];
            for (actual, expected) in actual.into_iter().zip(expected) {
                assert!((actual - expected).abs() <= 2.0e-7 * expected.abs().max(1.0),
                    "physical control {actual} differs from fresh-row derivative {expected}");
            }
        }
    }
}

