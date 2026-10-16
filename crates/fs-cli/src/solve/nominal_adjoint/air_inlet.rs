//! Bind each full-system inlet contraction to exactly one declared air branch.

use std::collections::{BTreeMap, BTreeSet};
use fs_airflow::conjugate::{AirPath, goal::{CoupledGoalError, maximum::pullback_inlet_temperatures}};
use fs_project::ConductionSetup;
use super::{Cx, SolveRefusal, bad, conduction_error, poll, row};

pub(super) fn rows(
    cx: &Cx<'_>, setup: &ConductionSetup, paths: &[AirPath],
    port_names: &[&str], reference_bars: &[f64],
) -> Result<Vec<String>, SolveRefusal> {
    poll(cx)?;
    let laws = super::super::conjugate::airflow_laws(setup)?;
    if laws.len() != port_names.len() { return Err(bad("declared air segments differ from the complete adjoint ports")); }
    let targets: BTreeMap<_, _> = laws.iter().map(|law| (law.target.as_str(), law)).collect();
    let mut branches: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for law in &laws {
        poll(cx)?;
        branches.entry(law.branch.as_str()).or_default().push(law.target.as_str());
    }
    if branches.len() != paths.len() { return Err(bad("one retained air path is required per declared inlet branch")); }
    let derivatives = pullback_inlet_temperatures(cx, paths, port_names, reference_bars, 64)
        .map_err(|error| match error {
            CoupledGoalError::Interrupted => conduction_error("cli-solve-cancelled",
                "air-inlet adjoint interrupted", "resume the accepted pipeline prefix"),
            error => bad(format!("complete air-inlet contraction refused: {error}")),
        })?;
    let mut seen = BTreeSet::new();
    let mut result = Vec::with_capacity(paths.len());
    for (ordinal, (path, derivative)) in paths.iter().zip(derivatives).enumerate() {
        poll(cx)?;
        let first = path.segments().first().ok_or_else(|| bad("empty retained inlet path"))?;
        let law = targets.get(first.region()).ok_or_else(|| bad("retained inlet has no declared segment"))?;
        let branch = law.branch.as_str();
        if !seen.insert(branch) { return Err(bad("a declared air inlet was split into multiple derivative paths")); }
        let declared = &branches[branch];
        if declared.len() != path.segments().len() { return Err(bad("retained inlet path omits a declared segment")); }
        for (segment, &target) in path.segments().iter().zip(declared) {
            poll(cx)?;
            if segment.region() != target
                || targets[target].inlet_temperature_k.to_bits() != path.inlet_temperature_k().to_bits()
            { return Err(bad("retained air inlet or stream-wise segment order differs from its declaration")); }
        }
        // One physical coordinate controls ALL segments of this branch. A row
        // per segment would be ambiguous to uncertainty/optimization consumers.
        result.push(row("air-inlet-temperature", branch, ordinal, "K", derivative)?);
    }
    poll(cx)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_airflow::conjugate::AirSegment;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_project::{ThermalBoundary, ThermalBoundaryCondition as B};
    use fs_project::spec::dims;
    use fs_qty::QtyAny;

    #[test]
    fn inlet_rows_aggregate_segments_and_reject_changed_physical_binding() {
        let law = |target: &str, order| ThermalBoundary { target: target.into(), condition: B::AirflowConvection {
            branch: "supply".into(), order, inlet_temperature: QtyAny::new(300.0, dims::TEMPERATURE),
            hydraulic_diameter: QtyAny::new(0.02, dims::LENGTH), flow_area: QtyAny::new(0.004, dims::AREA),
            channel_length: QtyAny::new(0.3, dims::LENGTH), correlation: "convection.gnielinski".into(),
        } };
        let setup = ConductionSetup { regions: vec![], boundaries: vec![law("out", 1), law("in", 0)],
            adiabatic_remainder: true, radiation: None, transient: None };
        let segments = vec![AirSegment::new("in", 0.1, 10.0).unwrap(), AirSegment::new("out", 0.1, 20.0).unwrap()];
        let paths = [AirPath::new(300.0, 0.01, 1000.0, segments.clone()).unwrap()];
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey { seed: 2, kernel_id: 832, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            let result = rows(&cx, &setup, &paths, &["in", "out"], &[0.5, 0.8]).unwrap();
            assert_eq!(result.len(), 1);
            let json = crate::json_read::JsonValue::parse(&result[0]).unwrap();
            assert_eq!(json.str_field("entity"), Some("supply"));
            assert_eq!(json.str_field("parameter_unit"), Some("K"));
            let changed = [AirPath::new(301.0, 0.01, 1000.0, segments.clone()).unwrap()];
            assert!(rows(&cx, &setup, &changed, &["in", "out"], &[0.5, 0.8]).is_err());
            let reversed = [AirPath::new(300.0, 0.01, 1000.0, segments.into_iter().rev().collect()).unwrap()];
            assert!(rows(&cx, &setup, &reversed, &["out", "in"], &[0.8, 0.5]).is_err());
        });
    }
}
