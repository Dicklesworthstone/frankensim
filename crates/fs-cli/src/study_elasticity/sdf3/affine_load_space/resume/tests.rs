use super::*;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-load-space.fsim"
));

#[test]
fn g5_affine_restore_continues_the_same_accepted_trajectory_and_charges_reconstruction() {
    let spec = spec::parse(&FIXTURE.replace(":max-updates 80", ":max-updates 4")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let reference = drive(&spec, &ledger, Some(4), &gate, None).unwrap();
    let reference = load(&ledger, &reference.pointer).unwrap();
    let prefix = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let prefix = load(&ledger, &prefix.pointer).unwrap();
    let restored = super::super::resume(&ledger, &prefix, Some(100), &gate).unwrap();
    let restored = load(&ledger, &restored.pointer).unwrap();
    for (key, kind) in [
        ("design", "study-design"),
        ("iterations", "study-iterations"),
    ] {
        assert_eq!(
            linked(&ledger, &restored.value, key, kind).unwrap(),
            linked(&ledger, &reference.value, key, kind).unwrap(),
        );
    }
    assert_eq!(
        restored.value.str_field("trace_hash"),
        reference.value.str_field("trace_hash")
    );
    assert_eq!(integer(&restored.value, "iterations_completed").unwrap(), 4);
    let restored_work = field(&restored.value, "optimizer_work").unwrap();
    let reference_work = field(&reference.value, "optimizer_work").unwrap();
    let extra = integer(restored_work, "restoration_evaluations").unwrap();
    assert!((1..=2).contains(&extra));
    assert_eq!(
        integer(restored_work, "evaluations").unwrap(),
        integer(reference_work, "evaluations").unwrap() + extra,
    );
    let prefix_report =
        document(&linked(&ledger, &prefix.value, "report_json", "study-report-json").unwrap())
            .unwrap();
    let restored_report =
        document(&linked(&ledger, &restored.value, "report_json", "study-report-json").unwrap())
            .unwrap();
    assert_eq!(
        field(&prefix_report, "gradient_check").unwrap(),
        field(&restored_report, "gradient_check").unwrap()
    );
    let before = Spent::read(&prefix.value).unwrap();
    let full = Spent::read(&reference.value).unwrap();
    let after = Spent::read(&restored.value).unwrap();
    assert_eq!(
        after.geometry.boxes,
        before.geometry.boxes + full.geometry.boxes
    );
    assert_eq!(
        after.geometry.points,
        before.geometry.points + full.geometry.points
    );
    assert!(after.linear.linear_iterations > full.linear.linear_iterations);
    assert!(after.linear.linear_solves > full.linear.linear_solves);
    assert!(after.wall_s >= before.wall_s);
    after.validate(&spec).unwrap();
    let finished = super::super::resume(&ledger, &restored, Some(1), &gate).unwrap();
    assert_eq!(
        finished.pointer,
        format!("study-{}", restored.hash.to_hex())
    );
    assert_eq!(finished.receipt, restored.bytes);
    assert_eq!(
        load(&ledger, &format!("study-{}", prefix.hash.to_hex()))
            .unwrap()
            .bytes,
        prefix.bytes
    );
}

#[test]
fn g3_affine_checkpoint_pins_self_weight_uncertainty_source_and_executable() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let out = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let old = load(&ledger, &out.pointer).unwrap();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    let recovery = recover(&spec, &ledger, &old, producer).unwrap();
    let checkpoint = document(&linked(&ledger, &old.value, "checkpoint", KIND).unwrap()).unwrap();
    assert_eq!(
        checkpoint.str_field("binding"),
        Some(recovery.checkpoint.binding.to_hex().as_str())
    );
    for source in [
        FIXTURE.replace(
            ":compliance-limit-j 0.000009765625",
            ":compliance-limit-j 0.00001953125",
        ),
        FIXTURE.replace("(0.0 -1.0 -0.5)", "(0.0 -1.0 -0.6)"),
        FIXTURE.replace(
            ":youngs-modulus-pa 1024000.0",
            ":youngs-modulus-pa 1126400.0",
        ),
        FIXTURE.replace(":density-kg-m3 0.2", ":density-kg-m3 0.3"),
        FIXTURE.replace(
            ":acceleration-m-s2 (0.0 0.0 -1.0)",
            ":acceleration-m-s2 (0.0 0.0 -1.1)",
        ),
        FIXTURE.replace(":threshold-offset 0.1", ":threshold-offset 0.2"),
        FIXTURE.replace(
            ":max-optimizer-evaluations 2000",
            ":max-optimizer-evaluations 3000",
        ),
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
fn g4_affine_restore_refuses_spent_evaluations_and_precancellation_without_losing_fields() {
    let spec = spec::parse(&FIXTURE.replace(
        ":max-optimizer-evaluations 2000",
        ":max-optimizer-evaluations 1",
    ))
    .unwrap();
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
    assert!(
        !linked(&ledger, &old.value, "design", "study-design")
            .unwrap()
            .is_empty()
    );
}

#[test]
fn g3_historical_affine_receipts_without_optimizer_state_cannot_replay_a_fake_resume() {
    let spec = spec::parse(FIXTURE).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let out = drive(&spec, &ledger, Some(0), &gate, None).unwrap();
    let mut old = load(&ledger, &out.pointer).unwrap();
    let JsonValue::Object(fields) = &mut old.value else {
        panic!("receipt object")
    };
    fields.retain(|(key, _)| key != "resume_mode" && key != "checkpoint");
    let (_, supported) = fields
        .iter_mut()
        .find(|(key, _)| key == "resume_supported")
        .unwrap();
    *supported = JsonValue::Bool(false);
    let producer = checkpoint::producer_identity(&gate).unwrap();
    let error = recover(&spec, &ledger, &old, producer).err().unwrap();
    assert_eq!(error.code, "cli-study-sdf3-resume-unsupported");
    assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, old.bytes);
}
