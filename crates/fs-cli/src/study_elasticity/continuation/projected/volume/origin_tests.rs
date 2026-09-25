use super::*;
use crate::study::elasticity::{canonical as canonical_study, parse as parse_study};
use fs_topols::evaluated::DesignEvaluationStage;
use fs_topols::refinement::projected::VolumeRefinementStage;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn policy(spec: &ElasticitySpec) -> &Controls {
    let Some(ProjectedControls::Volume(policy)) = &spec.projected else { panic!("volume policy") };
    policy
}
fn read(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, OptimizeReport, VolumeEvidence) {
    let receipt = json(out);
    let design = document(&linked(ledger, &receipt, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &receipt, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &receipt, &design, &rows).unwrap();
    let history = VolumeEvidence::read(receipt.path(&["continuation", "constraints"]).unwrap(), &report, policy(spec)).unwrap();
    (phi, report, history)
}
fn refined(source: &ElasticitySpec, pointer: &str) -> ElasticitySpec {
    let mut target = source.clone();
    target.base.physics.as_mut().unwrap().mesh_level += 1;
    let Some(ProjectedControls::Volume(policy)) = &mut target.projected else { panic!("volume") };
    policy.refine_from = Some(ContentHash::from_hex(pointer.strip_prefix("study-").unwrap()).unwrap());
    parse_study(&canonical_study(&target)).unwrap()
}
fn coarse(ledger: &Ledger) -> (ElasticitySpec, Outcome) {
    let spec = parse_study(BASE).unwrap();
    let out = drive(&spec, ledger, Some(1), &gate(), None).unwrap();
    assert_eq!(integer(&json(&out), "iterations_completed").unwrap(), 1);
    (spec, out)
}

#[test]
fn refinement_declaration_is_canonical_and_bad_or_unknown_sources_refuse() {
    let source = parse_study(BASE).unwrap();
    assert!(policy(&source).refine_from.is_none());
    let pointer = format!("study-{}", "ab".repeat(32));
    let fine = refined(&source, &pointer);
    assert_eq!(parse_study(&fine.canonical).unwrap().id, fine.id);
    assert_ne!(source.id, fine.id);
    for bad in ["\"not-a-study\"".to_string(), "3".into(), "\"study-123\"".into()] {
        assert!(parse_study(&fine.canonical.replace(&format!("\"{pointer}\""), &bad)).is_err());
    }
    let line = format!("    :refine-from \"{pointer}\"\n");
    assert!(parse_study(&fine.canonical.replace(&line, &format!("{line}{line}"))).is_err());
    let ledger = Ledger::open(":memory:").unwrap();
    assert!(drive(&fine, &ledger, None, &gate(), None).is_err());
    assert_eq!(ledger.table_count("ops").unwrap(), 0);
}

#[test]
fn refined_native_baseline_uses_the_accepted_design_and_preserves_real_lineage() {
    let ledger = Ledger::open(":memory:").unwrap();
    let (spec, old) = coarse(&ledger);
    let old_bytes = old.receipt.clone();
    let (phi, _, old_history) = read(&ledger, &old, &spec);
    let fine = refined(&spec, &old.pointer);
    let out = drive(&fine, &ledger, Some(0), &gate(), None).unwrap();
    let (actual, rows, history) = read(&ledger, &out, &fine);
    assert!(rows.rows.is_empty());
    assert_eq!(actual.n(), 2 * phi.n());
    assert_ne!(actual.nodes(), initial_phi(&fine).nodes());
    let oracle = fs_topols::evaluate_compliance_design(&actual, fixture(&fine), settings(&fine, fine.steps)).unwrap();
    assert!(Measured::from(oracle).same(history.baseline));
    let constraints = json(&out).path(&["continuation", "constraints"]).unwrap().clone();
    assert_eq!(constraints.f64_field("relative_reduction"), Some(0.0));
    assert_eq!(constraints.path(&["refinement_origin", "source_updates"]).and_then(JsonValue::as_f64), Some(1.0));
    assert_eq!(constraints.path(&["refinement_origin", "coarse_endpoint"]), Some(&document(old_history.current().json().as_bytes()).unwrap()));
    let source_hash = json(&out).str_field("source").and_then(ContentHash::from_hex).unwrap();
    let source_op = ledger.artifact_output_seal(&source_hash).unwrap().unwrap();
    let old_hash = ContentHash::from_hex(old.pointer.strip_prefix("study-").unwrap()).unwrap();
    assert!(ledger.edge_exists(source_op, &old_hash, EdgeRole::In).unwrap());
    assert!(ledger.edge_exists(source_op, &source_hash, EdgeRole::Out).unwrap());
    assert_eq!(load(&ledger, &old.pointer).unwrap().bytes, old_bytes);
    // A second handoff reconstructs the inherited pin support through ancestry.
    let finer = refined(&fine, &out.pointer);
    let next = drive(&finer, &ledger, Some(0), &gate(), None).unwrap();
    let (next_phi, _, _) = read(&ledger, &next, &finer);
    let ControlFlow::Continue(pins) = origin::fixed(&finer, &ledger, |_| ControlFlow::<()>::Continue(())).unwrap()
        else { panic!("fixed-node reconstruction stopped") };
    for (node, value) in pins { assert_eq!(next_phi.nodes()[node].to_bits(), value.to_bits()); }
}

#[test]
fn fine_grid_resume_keeps_the_new_baseline_and_does_not_replay_coarse_mechanics() {
    let ledger = Ledger::open(":memory:").unwrap();
    let (spec, old) = coarse(&ledger);
    let fine = refined(&spec, &old.pointer);
    let seed = drive(&fine, &ledger, Some(0), &gate(), None).unwrap();
    let before = json(&seed).path(&["continuation", "constraints", "baseline"]).unwrap().clone();
    let loaded = load(&ledger, &seed.pointer).unwrap();
    let resumed = drive_observed(&fine, &ledger, Some(1), &gate(), Some(&loaded), |stage| {
        assert!(!matches!(stage, VolumeStage::Origin(origin::Stage::SourceSetup(_) | origin::Stage::Refine(_))),
            "normal resume must not redo the source or refinement solve");
    }).unwrap();
    let fresh = drive(&fine, &ledger, Some(1), &gate(), None).unwrap();
    assert_eq!(fresh.status, resumed.status);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(), linked(&ledger, &json(&fresh), key, kind).unwrap());
    }
    assert_eq!(json(&resumed).path(&["continuation", "constraints", "baseline"]), Some(&before));
    assert_eq!(load(&ledger, &old.pointer).unwrap().bytes, old.receipt);
}

#[test]
fn refinement_cancellation_discards_fine_work_without_writing_a_child_study() {
    let ledger = Ledger::open(":memory:").unwrap();
    let (spec, old) = coarse(&ledger);
    let fine = refined(&spec, &old.pointer);
    let count = ledger.table_count("ops").unwrap();
    for final_stage in [false, true] {
        let cancel = gate();
        let mut reached = false;
        let error = drive_observed(&fine, &ledger, None, &cancel, None, |stage| {
            let should_stop = if final_stage {
                matches!(stage, VolumeStage::Origin(origin::Stage::Refine(VolumeRefinementStage::Publish)))
            } else {
                matches!(stage, VolumeStage::Origin(origin::Stage::Refine(VolumeRefinementStage::Baseline(
                    ProjectedSetupStage::Evaluation(DesignEvaluationStage::Solve(n))))) if n > 0)
            };
            if should_stop { reached = true; cancel.request(); }
        }).unwrap_err();
        assert!(reached);
        assert_eq!(error.exit, exit::CANCELLED);
        assert_eq!(ledger.table_count("ops").unwrap(), count);
        assert_eq!(load(&ledger, &old.pointer).unwrap().bytes, old.receipt);
    }
    assert!(drive(&fine, &ledger, Some(0), &gate(), None).is_ok());
}

#[test]
fn refinement_cannot_change_physics_or_use_a_same_level_restart() {
    let ledger = Ledger::open(":memory:").unwrap();
    let (spec, old) = coarse(&ledger);
    let fine = refined(&spec, &old.pointer);
    let before = ledger.table_count("ops").unwrap();
    for which in 0..4 {
        let mut changed = fine.clone();
        match which {
            0 => changed.youngs_pa *= 2.0,
            1 => changed.load_traction_pa *= 2.0,
            2 => {
                changed.base.constraints.as_mut().unwrap().volume_fraction = 0.7;
                let Some(ProjectedControls::Volume(policy)) = &mut changed.projected else { panic!() };
                policy.area.target = 0.7;
            }
            _ => changed.base.physics.as_mut().unwrap().mesh_level -= 1,
        }
        let changed = parse_study(&canonical_study(&changed)).unwrap();
        assert!(drive_observed(&changed, &ledger, None, &gate(), None, |stage| {
            assert!(!matches!(stage, VolumeStage::Origin(origin::Stage::SourceSetup(_) | origin::Stage::Refine(_))));
        }).is_err());
    }
    assert_eq!(ledger.table_count("ops").unwrap(), before);
}

#[test]
fn public_refined_study_and_report_use_the_same_retained_origin() {
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-refined-study-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let db = dir.join("study.db");
    let coarse_path = dir.join("coarse.fsim");
    let fine_path = dir.join("fine.fsim");
    fs::write(&coarse_path, BASE).unwrap();
    let invoke = |args: &[&str]| crate::run(args.iter().map(|s| s.to_string()));
    let old = invoke(&["--json", "study", coarse_path.to_str().unwrap(), db.to_str().unwrap(), "--budget", "1"]);
    assert_eq!(old.exit_code, exit::BUDGET, "{}", old.stderr);
    let old = JsonValue::parse(&old.stdout).unwrap();
    let fine = refined(&parse_study(BASE).unwrap(), old.str_field("run_id").unwrap());
    fs::write(&fine_path, &fine.canonical).unwrap();
    let out = invoke(&["--json", "study", fine_path.to_str().unwrap(), db.to_str().unwrap(), "--budget", "1"]);
    assert!(matches!(out.exit_code, exit::BUDGET | exit::REFUSED), "{}", out.stderr);
    let out = JsonValue::parse(&out.stdout).unwrap();
    assert!(matches!(out.str_field("status"), Some("budget-exhausted" | "no-feasible-descent")));
    let exported = invoke(&["--json", "report", out.str_field("run_id").unwrap(), db.to_str().unwrap()]);
    assert_eq!(exported.exit_code, exit::SUCCESS, "{}", exported.stderr);
    let paths = JsonValue::parse(&exported.stdout).unwrap();
    let report = document(&fs::read(paths.str_field("report_json").unwrap()).unwrap()).unwrap();
    assert_eq!(report.path(&["constraints", "refinement_origin"]), out.path(&["receipt", "continuation", "constraints", "refinement_origin"]));
    assert!(fs::read_to_string(paths.str_field("report_html").unwrap()).unwrap().contains("new work budget"));
}
