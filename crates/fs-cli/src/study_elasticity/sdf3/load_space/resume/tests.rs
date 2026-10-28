use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-load-space.fsim"));

#[test]
fn g5_native_load_space_restore_preserves_the_accepted_trajectory() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let reference = drive(&spec, &ledger, Some(4), &gate, None).unwrap();
    let reference = load(&ledger, &reference.pointer).unwrap();
    let prefix = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let prefix = load(&ledger, &prefix.pointer).unwrap();
    let restored = drive(&spec, &ledger, Some(3), &gate, Some(&prefix)).unwrap();
    let restored = load(&ledger, &restored.pointer).unwrap();
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &restored.value, key, kind).unwrap(),
            linked(&ledger, &reference.value, key, kind).unwrap());
    }
    assert_eq!(restored.value.str_field("trace_hash"), reference.value.str_field("trace_hash"));
    let restored_work = field(&restored.value, "optimizer_work").unwrap();
    let reference_work = field(&reference.value, "optimizer_work").unwrap();
    let extra = integer(restored_work, "restoration_evaluations").unwrap();
    assert!(extra >= 1);
    assert_eq!(integer(restored_work, "evaluations").unwrap(),
        integer(reference_work, "evaluations").unwrap() + extra);
    let prefix_report = document(&linked(&ledger, &prefix.value, "report_json", "study-report-json").unwrap()).unwrap();
    let restored_report = document(&linked(&ledger, &restored.value, "report_json", "study-report-json").unwrap()).unwrap();
    assert_eq!(field(&prefix_report, "gradient_check").unwrap(), field(&restored_report, "gradient_check").unwrap());
    assert_eq!(load(&ledger, &format!("study-{}", prefix.hash.to_hex())).unwrap().bytes, prefix.bytes);
}

#[test]
fn g3_load_space_checkpoint_pins_original_physics_and_executable() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let out = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let old = load(&ledger, &out.pointer).unwrap();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    recover(&spec, &ledger, &old, producer).unwrap();
    for source in [
        FIXTURE.replace(":compliance-limit-j 10.0", ":compliance-limit-j 11.0"),
        FIXTURE.replace("(0.0 -1.0 -0.5)", "(0.0 -1.0 -0.6)"),
        FIXTURE.replace(":youngs-modulus-pa 1.0", ":youngs-modulus-pa 1.1"),
        FIXTURE.replace(":threshold-offset 0.1", ":threshold-offset 0.2"),
        FIXTURE.replace(":max-optimizer-evaluations 2000", ":max-optimizer-evaluations 3000"),
        FIXTURE.replace(":linear-iterations 250000", ":linear-iterations 300000"),
    ] {
        let changed = spec::parse(&source).unwrap();
        assert_ne!(changed.id, spec.id);
        assert!(recover(&changed, &ledger, &old, producer).is_err());
    }
    assert!(recover(&spec, &ledger, &old, hash_bytes(b"different executable")).is_err());
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, old.bytes);
}

#[test]
fn g4_resume_refuses_spent_allowances_and_precancellation_without_losing_results() {
    let spec = spec::parse(&FIXTURE.replace(":max-optimizer-evaluations 2000",
        ":max-optimizer-evaluations 1")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let out = drive(&spec, &ledger, Some(0), &gate, None).unwrap();
    let old = load(&ledger, &out.pointer).unwrap();
    let error = drive(&spec, &ledger, Some(1), &gate, Some(&old)).unwrap_err();
    assert_eq!(error.exit, exit::BUDGET);
    assert_eq!(error.code, "cli-study-sdf3-resume-budget");
    gate.request();
    let error = drive(&spec, &ledger, Some(1), &gate, Some(&old)).unwrap_err();
    assert_eq!(error.exit, exit::CANCELLED);
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, old.bytes);
}
