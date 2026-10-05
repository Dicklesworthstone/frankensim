//! Independently controlled native fan banks through the COMPLETE thermal dual.
//! Fixed quadratic vent/leakage losses keep passive branch-flow fractions fixed.
//! A single bank retains q(s)/s scaling; heterogeneous series/parallel banks use
//! the nominal hydraulic derivative of each independently declared speed.

use std::collections::BTreeMap;
use fs_airflow::conjugate::goal::{CoupledGoalError, maximum::pullback_transport_controls};
use fs_conduction::{ScalarField, ThermalBc};
use fs_convection::{CorrelationInputs, NusseltEvaluation, evaluate};
#[cfg(test)]
use fs_convection::CorrelationId;
use fs_project::fansystem::FanSystemTopology;
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, conduction_error,
    contractions, finite, poll, row, unsupported, zeros};
use super::super::conjugate::{self, AIR_DYNAMIC_VISCOSITY_PA_S, AIR_PRANDTL,
    AIR_SPECIFIC_HEAT_J_KG_K, AIR_THERMAL_CONDUCTIVITY_W_M_K};

mod hydraulic;

/// Append supported controls, or explicit per-bank reasons. Unsupported
/// parameter maps do not invalidate the already complete thermal state dual.
pub(super) fn append(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, lambda: &[f64],
    rows: &mut Vec<String>, missing: &mut Vec<String>,
) -> Result<(), SolveRefusal> {
    poll(cx)?;
    let Some(system) = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref()) else { return Ok(()); };
    system.validate().map_err(|e| bad(format!("fan-speed declaration refused: {}", e.detail)))?;
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("fan derivative has no retained thermal operator"))?;
    let unavailable = |missing: &mut Vec<String>, reason: &str| {
        for bank in &system.banks { missing.push(unsupported("fan-speed-ratio", &bank.bank_id, reason)); }
    };
    if system.banks.len() > 64 { return Err(bad("fan-speed derivative exceeds the 64-bank budget")); }
    if system.banks.iter().any(|bank| bank.speed_ratio <= bank.speed_ratio_domain.0
        || bank.speed_ratio >= bank.speed_ratio_domain.1) {
        unavailable(missing, "a two-sided fan-speed derivative requires interior points of the declared speed domains");
        return Ok(());
    }
    if data.air_paths.is_empty() {
        unavailable(missing, "no native airflow-convection path connects this fan to the thermal field");
        return Ok(());
    }
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("fan derivative has no declared thermal setup"))?;
    let laws = conjugate::airflow_laws(setup)?;
    let m = data.air_paths.iter().try_fold(0_usize, |n, p| n.checked_add(p.segments().len()))
        .ok_or_else(|| bad("fan derivative air-trace count overflow"))?;
    if m == 0 || m > 64 || laws.len() != m || lambda.len() != solved.mesh.vertex_count() {
        return Err(bad("fan derivative requires every retained air trace within the 64-port budget"));
    }
    let ports: Vec<_> = data.air_paths.iter().flat_map(|p| p.segments()).collect();
    let mut regions = BTreeMap::new();
    let mut names = Vec::with_capacity(m);
    let mut slopes = Vec::with_capacity(m);
    let mut at = 0;
    for path in &data.air_paths {
        poll(cx)?;
        let branch = &laws[at].branch;
        // Transport properties are the SAME frozen native values. Recover Re
        // from retained mass capacity, not a second hydraulic operating solve.
        let mass_flow = finite(path.capacity_rate_w_per_k()/AIR_SPECIFIC_HEAT_J_KG_K)?;
        for segment in path.segments() {
            poll(cx)?;
            let law = &laws[at];
            if &law.branch != branch || law.target != segment.region()
                || law.inlet_temperature_k.to_bits() != path.inlet_temperature_k().to_bits() {
                return Err(bad("fan derivative air-branch ordering differs from the declared producer"));
            }
            let region = data.boundary.region_names().iter().position(|n| n == segment.region())
                .ok_or_else(|| bad("fan derivative has no original convective trace"))?;
            if regions.insert(region, at).is_some() { return Err(bad("fan derivative repeats an air trace")); }
            let ThermalBc::Robin { htc: ScalarField::Uniform(h), .. } = &data.boundary.conditions()[region]
                else { return Err(bad("fan derivative requires original uniform convective rows")); };
            if *h != segment.htc_w_per_m2_k() { return Err(bad("fan derivative cannot substitute the radiative or another coefficient for convection")); }
            let re = finite((mass_flow/law.flow_area_m2)*law.hydraulic_diameter_m/AIR_DYNAMIC_VISCOSITY_PA_S)?;
            let inputs = |r| CorrelationInputs::forced(r, AIR_PRANDTL)
                .with_length_ratio(law.channel_length_m/law.hydraulic_diameter_m);
            let Ok(evaluation) = evaluate(law.correlation, inputs(re)) else {
                unavailable(missing, "the retained flow point has no admitted convection-card differential; no extrapolation is supplied");
                return Ok(());
            };
            let expected_h = finite(evaluation.evidence().value*AIR_THERMAL_CONDUCTIVITY_W_M_K/law.hydraulic_diameter_m)?;
            if (expected_h-h).abs() > 512.0*f64::EPSILON*expected_h.abs().max(h.abs()) {
                return Err(bad("retained convection coefficient differs from its fan-speed correlation law"));
            }
            // The local derivative cannot extend the card domain. Checking
            // neighbouring representable Re is admission, not finite differences.
            let interior = [re.next_down(), re.next_up()].into_iter().all(|r|
                evaluate(law.correlation, inputs(r)).is_ok_and(|v| v.evidence().model.in_domain));
            let Some(slope) = reynolds_elasticity(&evaluation)?.filter(|_| interior) else {
                unavailable(missing, "a convection card has no admitted smooth Reynolds derivative at this operating point; use explicit coordinate-secant calibration");
                return Ok(());
            };
            names.push(segment.region()); slopes.push(slope); at += 1;
        }
    }
    let mut areas = zeros(m)?;
    let mut means = zeros(m)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let Some(&i) = data.boundary.region_for(slot).and_then(|r| regions.get(&r)) else { continue; };
        areas[i] = finite(areas[i]+face.area)?;
        for &v in &face.vertices { means[i] = finite(means[i]+(face.area/3.0)*solved.solution.temperature[v as usize])?; }
    }
    for i in 0..m {
        let area = ports[i].area_m2();
        if areas[i] <= 0.0 || (areas[i]-area).abs() > 128.0*f64::EPSILON*areas[i].max(area) {
            return Err(bad("fan derivative requires matching physical wetted areas"));
        }
        means[i] = finite(means[i]/areas[i])?;
    }
    let mut references = Vec::with_capacity(m);
    let mut start = 0;
    for path in &data.air_paths {
        poll(cx)?;
        let end = start+path.segments().len();
        references.extend(path.march(&means[start..end])
            .map_err(|e| bad(format!("fan derivative air march refused: {e}")))?.reference_temperatures_k());
        start = end;
    }
    let mut reference_bars = zeros(m)?;
    let mut direct_h_bars = zeros(m)?;
    for (slot, face) in solved.mesh.boundary().iter().enumerate() {
        if slot % 512 == 0 { poll(cx)?; }
        let Some(&i) = data.boundary.region_for(slot).and_then(|r| regions.get(&r)) else { continue; };
        let vertices = face.vertices.map(|v| v as usize);
        let h = ports[i].htc_w_per_m2_k();
        let bars = contractions::boundary(face.area, vertices.map(|v| lambda[v]),
            vertices.map(|v| solved.solution.temperature[v]), h, references[i]);
        reference_bars[i] = finite(reference_bars[i]+bars[1])?;
        direct_h_bars[i] = finite(direct_h_bars[i]+finite(h*bars[0])?)?;
    }
    let controls = pullback_transport_controls(cx, &data.air_paths, &names, &means,
        &reference_bars, &direct_h_bars, 64).map_err(|e| match e {
            CoupledGoalError::Interrupted => conduction_error("cli-solve-cancelled",
                "fan-speed adjoint interrupted", "resume the accepted pipeline prefix"),
            e => bad(format!("complete fan thermal controls refused: {e}")),
        })?;
    let mut total = 0.0;
    for value in &controls.log_capacity_rates { total = finite(total+value)?; }
    for (&bar, slope) in controls.log_htc.iter().zip(slopes) { total = finite(total+finite(bar*slope)?)?; }
    let weights = if matches!(&system.topology, FanSystemTopology::Single) {
        vec![1.0] // Preserve the existing single-bank affinity calculation.
    } else {
        let Some(weights) = hydraulic::flow_weights(cx, spec, &data.air_paths, &laws)? else {
            unavailable(missing, "a member curve is at a nonsmooth knot, nonunique parallel inverse or validity endpoint; no two-sided hydraulic derivative is supplied");
            return Ok(());
        };
        weights
    };
    if weights.len() != system.banks.len() { return Err(bad("hydraulic derivatives lost bank identity")); }
    // Each project field is an independently controlled ABSOLUTE ratio, not
    // ln(speed), and not a common speed change applied to all fan banks.
    for (ordinal, (bank, weight)) in system.banks.iter().zip(weights).enumerate() {
        poll(cx)?;
        rows.push(row("fan-speed-ratio", &bank.bank_id, ordinal, "1",
            finite(finite(total*weight)?/bank.speed_ratio)?)?);
    }
    poll(cx)
}

// Keep the native binding/domain checks above, but use the correlation owner's
// derivative rather than maintaining another copy of its physical formulas.
fn reynolds_elasticity(evaluation: &NusseltEvaluation) -> Result<Option<f64>, SolveRefusal> {
    evaluation.reynolds_elasticity().map_err(|error|
        bad(format!("fan-convection derivative refused: {error}")))
}

#[cfg(test)]
mod tests;
