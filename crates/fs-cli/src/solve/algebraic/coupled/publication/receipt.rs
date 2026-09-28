//! Rebuild live air results while retaining the initial exchange as history.
//! Geometry/correlation/flow metadata are copied only after segment identities
//! agree; no old temperature or heat rate survives in the live result fields.

use fs_airflow::conjugate::goal::maximum::physical::AcceptedLinearCooling;
use fs_conduction::{ScalarField, ThermalBc};
use crate::json_read::JsonValue as Json;
use super::{AirPath, Cx, SolveRefusal, cancelled, error, json_string, number};

pub(super) fn set(value: &mut Json, key: &str, next: Json) -> Result<(), SolveRefusal> {
    let Json::Object(rows) = value else { return Err(error("air receipt member is not an object")); };
    if let Some((_, old)) = rows.iter_mut().find(|(name, _)| name == key) { *old = next; }
    else { rows.push((key.into(), next)); }
    Ok(())
}
pub(super) fn numeric(value: f64) -> Result<Json, SolveRefusal> {
    Ok(Json::Number { value, raw: number(value)? })
}
fn text(value: &str) -> Json { Json::Str(value.into()) }
fn checked_sum(values: impl Iterator<Item = f64>) -> Result<f64, SolveRefusal> {
    let mut total = 0.0;
    for value in values { total += value; number(total)?; }
    Ok(total)
}

fn render(cx: &Cx<'_>, value: &Json, out: &mut String) -> Result<(), SolveRefusal> {
    encode(value, out, &mut || cx.checkpoint().map_err(|_| cancelled()))
}

pub(super) fn encode(
    value: &Json, out: &mut String, checkpoint: &mut impl FnMut() -> Result<(), SolveRefusal>,
) -> Result<(), SolveRefusal> {
    checkpoint()?;
    match value {
        Json::Null => out.push_str("null"),
        Json::Bool(value) => out.push_str(if *value { "true" } else { "false" }),
        Json::Number { raw, .. } => out.push_str(raw),
        Json::Str(value) => out.push_str(&json_string(value)),
        Json::Array(values) => {
            out.push('[');
            for (i, value) in values.iter().enumerate() {
                if i != 0 { out.push(','); }
                encode(value, out, checkpoint)?;
            }
            out.push(']');
        }
        Json::Object(rows) => {
            out.push('{');
            for (i, (key, value)) in rows.iter().enumerate() {
                if i != 0 { out.push(','); }
                out.push_str(&json_string(key)); out.push(':');
                encode(value, out, checkpoint)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

pub(super) fn rebuild(
    cx: &Cx<'_>, original: &str, paths: &[AirPath], accepted: &AcceptedLinearCooling,
    iterations: usize,
) -> Result<String, SolveRefusal> {
    cx.checkpoint().map_err(|_| cancelled())?;
    // A bounded internal receipt, not an unbounded external JSON rewrite.
    if original.len() > 1024 * 1024 || paths.is_empty()
        || paths.len() != accepted.air.len() || paths.len() != accepted.branches.len()
    { return Err(error("air publication has invalid branch arity or an oversized initial receipt")); }
    let history = Json::parse(original).map_err(|e| error(e.to_string()))?;
    let old_branches: Vec<&Json> = if paths.len() == 1 { vec![&history] } else {
        let branches = history.get("branches").and_then(Json::as_array)
            .ok_or_else(|| error("multi-branch history has no branch array"))?;
        if branches.len() != paths.len() { return Err(error("initial air branch count changed")); }
        branches.iter().collect()
    };
    let mut branches = Vec::with_capacity(paths.len());
    let mut solid_totals = Vec::with_capacity(paths.len());
    let decomposition = checked_sum(accepted.solid.report.robin_fluxes.iter().map(|f| f.heat_rate_w))?
        - accepted.solid.report.energy.robin_out_w;
    for (i, ((path, air), old)) in paths.iter().zip(&accepted.air).zip(old_branches).enumerate() {
        cx.checkpoint().map_err(|_| cancelled())?;
        if old.f64_field("inlet_k") != Some(path.inlet_temperature_k()) {
            return Err(error("initial receipt and retained air inlet differ"));
        }
        let old_segments = old.get("segments").and_then(Json::as_array)
            .ok_or_else(|| error("initial branch has no segment array"))?;
        if old_segments.len() != air.segments.len() || old_segments.len() != path.segments().len() {
            return Err(error("initial and accepted air segment counts differ"));
        }
        let mut segments = Vec::with_capacity(old_segments.len());
        let mut solid_total = 0.0;
        let mut max_imbalance = 0.0_f64;
        for ((old, state), segment) in old_segments.iter().zip(&air.segments).zip(path.segments()) {
            cx.checkpoint().map_err(|_| cancelled())?;
            if old.str_field("target") != Some(segment.region()) || state.region != segment.region() {
                return Err(error("initial and accepted segment identities differ"));
            }
            let flux = accepted.solid.report.robin_fluxes.iter().find(|flux| flux.region == state.region)
                .ok_or_else(|| error("accepted air segment has no physical Robin flux"))?;
            let region = accepted.boundary.region_names().iter().position(|name| *name == state.region)
                .ok_or_else(|| error("accepted air segment has no boundary region"))?;
            let reference = match &accepted.boundary.conditions()[region] {
                ThermalBc::Robin { t_ref: ScalarField::Uniform(value), .. } => *value,
                _ => return Err(error("accepted air boundary is not uniform Robin")),
            };
            let mut row = old.clone();
            for (key, value) in [
                ("air_in_k", state.inlet_temperature_k), ("air_out_k", state.outlet_temperature_k),
                ("reference_k", reference), ("marched_reference_k", state.reference_temperature_k),
                ("wall_temperature_k", flux.mean_wall_temperature_k), ("ntu", state.ntu),
                ("effectiveness", state.effectiveness), ("solid_heat_rate_w", flux.heat_rate_w),
                ("air_heat_rate_w", state.heat_rate_w), ("imbalance_w", flux.heat_rate_w - state.heat_rate_w),
            ] { set(&mut row, key, numeric(value)?)?; }
            solid_total += flux.heat_rate_w; number(solid_total)?;
            max_imbalance = max_imbalance.max((flux.heat_rate_w - state.heat_rate_w).abs());
            segments.push(row);
        }
        solid_totals.push(solid_total);
        let mut branch = old.clone();
        set(&mut branch, "segments", Json::Array(segments))?;
        for (key, value) in [
            ("outlet_k", air.outlet_temperature_k), ("solid_total_w", solid_total),
            ("air_total_w", air.total_heat_rate_w), ("interface_imbalance_w", solid_total - air.total_heat_rate_w),
            ("max_region_imbalance_w", max_imbalance), ("balance_tolerance_w", accepted.branches[i].watt_limit),
            ("reference_delta_k", accepted.branches[i].reference_delta_k),
            ("enthalpy_imbalance_w", accepted.branches[i].enthalpy_imbalance_w),
        ] { set(&mut branch, key, numeric(value)?)?; }
        set(&mut branch, "decomposition_residual_w", if paths.len() == 1 { numeric(decomposition)? } else { Json::Null })?;
        stamp(&mut branch, iterations)?;
        branches.push(branch);
    }
    let mut output = if branches.len() == 1 { branches.remove(0) } else {
        let mut output = history.clone();
        let solid_total = checked_sum(solid_totals.into_iter())?;
        let air_total = checked_sum(accepted.air.iter().map(|air| air.total_heat_rate_w))?;
        set(&mut output, "branches", Json::Array(branches))?;
        for (key, value) in [("solid_total_w", solid_total), ("air_total_w", air_total),
            ("interface_imbalance_w", solid_total - air_total), ("decomposition_residual_w", decomposition),
            ("balance_tolerance_w", checked_sum(accepted.branches.iter().map(|b| b.watt_limit))?)]
        { set(&mut output, key, numeric(value)?)?; }
        stamp(&mut output, iterations)?;
        output
    };
    set(&mut output, "publication_schema", text("frankensim.cli.accepted-cooling.v1"))?;
    set(&mut output, "initial_exchange", history)?;
    let mut result = String::new();
    render(cx, &output, &mut result)?;
    cx.checkpoint().map_err(|_| cancelled())?;
    Ok(result)
}

fn stamp(value: &mut Json, iterations: usize) -> Result<(), SolveRefusal> {
    set(value, "iterations", Json::Number { value: iterations as f64, raw: iterations.to_string() })?;
    set(value, "iteration_scope", text("shared coupled correction; initial exchange is historical"))?;
    set(value, "worst_recorded_imbalance_w", Json::Null)?;
    set(value, "acceleration", Json::Object(vec![("method".into(), text("fgmres-physical-goal-polish"))]))?;
    set(value, "authority", text("physically reassembled accepted cooling field with refreshed air; initial partitioned exchange retained separately"))?;
    Ok(())
}
