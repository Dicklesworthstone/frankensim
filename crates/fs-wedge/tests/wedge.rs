//! Battery for go-to-market wedge selection (addendum Proposal 7). Verifies
//! historical-score supersession, evidence-complete measured inputs, workspace
//! evidence drift, candidate rankings, and the cycle-time kill criterion.

use fs_wedge::{
    CHT_BASELINE, EvidenceKind, InputAxis, Measurement, Readiness, STRONG_THRESHOLD, ScoreUse,
    WEDGE_DOCTRINE, WedgeCriterion, audit, chosen_wedge, four_criteria, measured_inputs_for,
    measured_wedge_inputs, to_json, verticals,
};
use std::path::Path;

#[test]
fn the_historical_beachhead_is_conjugate_heat_transfer() {
    let w = chosen_wedge();
    assert_eq!(w.name, "conjugate-heat-transfer");
    assert_eq!(w.rank, 1);
    // it exercises incremental re-solve (2), adjoints (1), the ladder (3),
    // and the evidence package (12).
    assert!(w.exercises.contains(&"2") && w.exercises.contains(&"3"));
}

#[test]
fn historical_scores_are_preserved_but_superseded_for_decisions() {
    let w = chosen_wedge();
    assert_eq!(w.score_use, ScoreUse::SupersededForDecisionUse);
    for c in four_criteria() {
        // Replay retains the plan's values; the decision API refuses them.
        assert!(
            w.score(c) >= STRONG_THRESHOLD,
            "historical {} score changed on {}",
            w.name,
            c.label()
        );
        assert_eq!(w.decision_score(c), None);
    }
    assert!(w.weakest_criterion_score() >= STRONG_THRESHOLD);
    assert!(
        verticals()
            .iter()
            .all(|vertical| !vertical.score_use.permits_decision())
    );
}

#[test]
fn every_candidate_has_complete_measured_inputs_on_all_four_axes() {
    let inputs = measured_wedge_inputs();
    assert_eq!(inputs.len(), verticals().len());
    for vertical in verticals() {
        let measured = measured_inputs_for(vertical.name)
            .unwrap_or_else(|| panic!("missing measured inputs for {}", vertical.name));
        assert!(measured.is_complete(), "incomplete: {}", measured.vertical);
        assert!(!measured.kernels.is_empty());
        assert!(!measured.validation_data.is_empty());
        assert!(!measured.cad_burden.is_empty());
        assert!(!measured.compute_cost.is_empty());
        for measurement in measured.measurements() {
            assert!(measurement.is_complete(), "{measurement:?}");
            assert!(!measurement.evidence.is_empty());
            assert!(
                measurement
                    .evidence
                    .iter()
                    .all(|pointer| pointer.is_complete())
            );
        }
    }
}

#[test]
fn absent_inputs_cannot_carry_strong_scores() {
    for inputs in measured_wedge_inputs() {
        for measurement in inputs.measurements() {
            assert!(
                measurement.score <= measurement.readiness.score_ceiling(),
                "{} has score {} above {:?} ceiling {}",
                inputs.vertical,
                measurement.score,
                measurement.readiness,
                measurement.readiness.score_ceiling()
            );
            if measurement.readiness == Readiness::Absent {
                assert!(
                    measurement.score < STRONG_THRESHOLD,
                    "absent {} input scored {}",
                    inputs.vertical,
                    measurement.score
                );
            }
        }
    }
}

fn check_workspace_measurement(
    root: &Path,
    vertical: &str,
    axis: InputAxis,
    label: &str,
    measurement: Measurement,
) -> Vec<String> {
    let mut failures = Vec::new();
    for pointer in measurement
        .evidence
        .iter()
        .filter(|pointer| pointer.kind == EvidenceKind::WorkspacePath)
    {
        let path = root.join(pointer.reference);
        let result = std::fs::read_to_string(&path);
        let (passed, detail) = match result {
            Ok(contents) if contents.contains(pointer.locator) => {
                (true, "marker-found".to_string())
            }
            Ok(_) => (false, format!("missing marker {:?}", pointer.locator)),
            Err(error) => (false, format!("read failed: {error}")),
        };
        eprintln!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            if passed { "PASS" } else { "FAIL" },
            vertical,
            axis.label(),
            label,
            pointer.reference,
            detail
        );
        if !passed {
            failures.push(format!(
                "{} {} {}: {} ({detail})",
                vertical,
                axis.label(),
                label,
                pointer.reference
            ));
        }
    }
    failures
}

#[test]
fn workspace_evidence_paths_and_markers_have_not_drifted() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let root = manifest
        .parent()
        .and_then(Path::parent)
        .expect("fs-wedge lives at <workspace>/crates/fs-wedge");
    let mut failures = Vec::new();

    eprintln!("RESULT\tVERTICAL\tAXIS\tENTRY\tPATH\tDETAIL");
    for inputs in measured_wedge_inputs() {
        for entry in inputs.kernels {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::KernelReadiness,
                entry.capability,
                entry.measurement,
            ));
        }
        for entry in inputs.validation_data {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::ValidationDataAccess,
                entry.dataset,
                entry.measurement,
            ));
        }
        for entry in inputs.cad_burden {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::CadBurden,
                entry.required_geometry,
                entry.measurement,
            ));
        }
        for entry in inputs.compute_cost {
            failures.extend(check_workspace_measurement(
                root,
                inputs.vertical,
                InputAxis::ComputeCost,
                entry.rung,
                entry.measurement,
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "evidence drift:\n{}",
        failures.join("\n")
    );
}

#[test]
fn three_verticals_are_ranked_with_proposal_mappings() {
    let vs = verticals();
    assert_eq!(vs.len(), 3);
    let mut ranks: Vec<u8> = vs.iter().map(|v| v.rank).collect();
    ranks.sort_unstable();
    assert_eq!(ranks, vec![1, 2, 3]);
    // second vertical exercises Proposal 1; third exercises 11 and 4.
    let aero = vs
        .iter()
        .find(|v| v.name == "aeroelastic-screening")
        .unwrap();
    assert_eq!(aero.rank, 2);
    assert!(aero.exercises.contains(&"1"));
    let am = vs
        .iter()
        .find(|v| v.name == "additive-manufacturing-distortion")
        .unwrap();
    assert_eq!(am.rank, 3);
    assert!(am.exercises.contains(&"11") && am.exercises.contains(&"4"));
    // every vertical names at least one exercised proposal.
    assert!(vs.iter().all(|v| !v.exercises.is_empty()));
}

#[test]
fn the_cycle_time_kill_criterion_is_measurable() {
    assert!((CHT_BASELINE.baseline_days - 5.0).abs() < 1e-12);
    assert!((CHT_BASELINE.target_reduction - 3.0).abs() < 1e-12);
    assert_eq!(CHT_BASELINE.kill_within_quarters, 2);
    // a 1.5-day cycle is a 3.33x reduction -> meets the criterion.
    assert!(CHT_BASELINE.meets_kill_criterion(1.5));
    // a 2-day cycle is only 2.5x -> does not.
    assert!(!CHT_BASELINE.meets_kill_criterion(2.0));
    // guard against divide-by-zero.
    assert!(!CHT_BASELINE.meets_kill_criterion(0.0));
}

#[test]
fn the_audit_is_complete() {
    let a = audit();
    assert!(a.ok(), "gaps: {:?}", a.gaps);
    assert!(a.passed("historic-scores-superseded"));
    assert!(a.passed("measured-inputs-complete"));
    assert!(a.passed("no-absent-strong-scores"));
    assert!(a.passed("ranks-complete"));
    assert!(a.passed("all-exercise-proposals"));
    assert!(a.passed("kill-criterion-measurable"));
    assert_eq!(a.checks.len(), 6);
}

#[test]
fn the_negative_doctrine_is_stated() {
    // the load-bearing anti-pattern: don't sell against peak single-physics.
    assert!(
        WEDGE_DOCTRINE
            .to_lowercase()
            .contains("peak single-physics")
    );
    // criterion labels are unique.
    let labels: Vec<&str> = WedgeCriterion::ALL.iter().map(|c| c.label()).collect();
    let mut sorted = labels.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), labels.len());
}

#[test]
fn json_is_well_formed_and_deterministic() {
    let j = to_json();
    assert_eq!(j, to_json());
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains("conjugate-heat-transfer"));
    assert!(j.contains("\"score_use\":\"superseded-for-decision-use\""));
    assert!(j.contains("\"measured_inputs\":"));
    assert!(j.contains("\"validation_data\":"));
    assert!(j.contains("NIST Additive Manufacturing Benchmark Test Series"));
    assert!(j.contains("\"target_reduction\":3"));
    assert_eq!(j.matches("\"rank\":").count(), 3);
}
