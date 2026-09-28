use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));

fn source() -> String {
    FIXTURE.replace(":updates-per-stage 3", ":updates-per-stage 1")
        .replace(":schedule ((1.0 1.0) (2.0 2.0) (3.0 8.0))", ":schedule ((1.0 1.0) (2.0 2.0))")
}
fn spec() -> Spec { spec::parse(&source()).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { JsonValue::parse(&out.receipt).unwrap() }
fn keep(_: &Outcome, _: usize) -> Result<()> { Ok(()) }
fn bytes(ledger: &Ledger, value: &JsonValue, key: &str) -> Vec<u8> {
    let kind = match key {
        "design" => "study-design", "iterations" => "study-iterations",
        "checkpoint" => KIND, _ => panic!("unsupported test artifact"),
    };
    linked(ledger, value, key, kind).unwrap()
}
fn set(value: &mut JsonValue, key: &str, replacement: JsonValue) {
    let JsonValue::Object(members) = value else { panic!("object expected") };
    members.iter_mut().find(|(name, _)| name == key).unwrap().1 = replacement;
}

#[test]
fn g1_g5_split_stages_replay_exact_fields_and_charge_all_repeated_work() {
    let spec = spec();
    let reference_ledger = Ledger::open(":memory:").unwrap();
    let full = drive(&spec, &reference_ledger, None, &gate(), None, keep).unwrap();
    assert_eq!(full.status, "completed");
    let reference = json(&full);
    let ledger = Ledger::open(":memory:").unwrap();
    let first = drive(&spec, &ledger, Some(1), &gate(), None, keep).unwrap();
    assert_eq!(first.status, "budget-exhausted");
    let first_json = json(&first);
    assert_eq!(integer(&first_json, "stages_completed").unwrap(), 1);
    assert_eq!(first_json.get("resume_supported"), Some(&JsonValue::Bool(true)));
    let original = bytes(&ledger, &first_json, "checkpoint");
    let old = load(&ledger, &first.pointer).unwrap();
    let mut seen = Vec::new();
    let resumed = drive(&spec, &ledger, Some(1), &gate(), Some(&old), |out, stage| {
        seen.push(stage);
        assert_eq!(load(&ledger, &out.pointer).unwrap().bytes, out.receipt,
            "observer sees a fully committed and sealed receipt");
        Ok(())
    }).unwrap();
    assert_eq!(seen, vec![2]);
    assert_eq!(resumed.status, "completed");
    let resumed_json = json(&resumed);
    assert_eq!(integer(&resumed_json, "stages_replayed").unwrap(), 1);
    assert_eq!(resumed_json.str_field("predecessor"), Some(old.hash.to_hex().as_str()));
    for key in ["design", "iterations", "checkpoint"] {
        assert_eq!(bytes(&reference_ledger, &reference, key), bytes(&ledger, &resumed_json, key));
    }
    assert_eq!(original, bytes(&ledger, &first_json, "checkpoint"));
    let initial = Spent::read(&first_json).unwrap();
    let uninterrupted = Spent::read(&reference).unwrap();
    let repeated = Spent::read(&resumed_json).unwrap();
    assert_eq!(repeated.linear, initial.add(uninterrupted.linear, uninterrupted.geometry, 0.0).unwrap().linear);
    assert_eq!(repeated.geometry, initial.add(uninterrupted.linear, uninterrupted.geometry, 0.0).unwrap().geometry);
    assert!(repeated.wall_s >= initial.wall_s);
    repeated.validate(&spec).unwrap();
    let complete = load(&ledger, &resumed.pointer).unwrap();
    let again = drive(&spec, &ledger, None, &gate(), Some(&complete), |_, _| {
        panic!("completed resume must not publish or solve another stage")
    }).unwrap();
    assert_eq!(again.pointer, resumed.pointer);
    assert_eq!(again.receipt, resumed.receipt);
}

#[test]
fn g4_interruption_after_commit_leaves_a_real_resumable_checkpoint() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let mut pointer = None;
    let failure = drive(&spec, &ledger, None, &gate(), None, |out, stage| {
        assert_eq!(stage, 1);
        pointer = Some(out.pointer.clone());
        Err(fail("test-interruption", "stop after first durable stage"))
    }).unwrap_err();
    let pointer = pointer.unwrap();
    assert_eq!(failure.code, "test-interruption");
    assert!(failure.message.contains(&pointer));
    let old = load(&ledger, &pointer).unwrap();
    assert_eq!(old.value.str_field("status"), Some("checkpointed"));
    let resumed = drive(&spec, &ledger, None, &gate(), Some(&old), keep).unwrap();
    assert_eq!(resumed.status, "completed");
    assert_eq!(integer(&json(&resumed), "stages_completed").unwrap(), 2);
}

#[test]
fn g4_replay_cannot_renew_the_original_krylov_allowance() {
    let reference_ledger = Ledger::open(":memory:").unwrap();
    let first = drive(&spec(), &reference_ledger, Some(1), &gate(), None, keep).unwrap();
    let cost = Spent::read(&json(&first)).unwrap().linear.linear_iterations;
    assert!(cost > 1);
    // Enough for the first accepted stage, NOT enough to pay for its replay.
    let restricted_source = source().replace(":linear-iterations 250000",
        &format!(":linear-iterations {}", 2 * cost - 1));
    let restricted = spec::parse(&restricted_source).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = drive(&restricted, &ledger, Some(1), &gate(), None, keep).unwrap();
    assert_eq!(first.status, "budget-exhausted");
    let old = load(&ledger, &first.pointer).unwrap();
    let before = bytes(&ledger, &old.value, "checkpoint");
    let error = drive(&restricted, &ledger, None, &gate(), Some(&old), |_, _| {
        panic!("no publication before full replay verification")
    }).unwrap_err();
    assert_eq!(error.code, "cli-study-sdf3-replay-incomplete");
    assert_eq!(error.exit, exit::BUDGET);
    assert!(error.message.contains(&first.pointer));
    assert_eq!(before, bytes(&ledger, &old.value, "checkpoint"));
    assert_eq!(load(&ledger, &first.pointer).unwrap().bytes, old.bytes);
}

/// Manufacture a new structurally sealed receipt through the public ledger API.
/// Neither that seal nor a well-shaped state is authority for the native field:
/// the actual numerical replay must independently reproduce its complete hash.
fn forged_state(ledger: &Ledger, spec: &Spec, old: &Loaded) -> Loaded {
    let original_hash = old.value.str_field("checkpoint").unwrap();
    let original = String::from_utf8(bytes(ledger, &old.value, "checkpoint")).unwrap();
    let marker = "\"native_displacement_bits\":[[\"";
    let digit = original.find(marker).unwrap() + marker.len() + 15;
    let mut altered = original.as_bytes().to_vec();
    altered[digit] = if altered[digit] == b'0' { b'1' } else { b'0' };
    ledger.begin().unwrap();
    let state = ledger.put_artifact(KIND, &altered, None).unwrap();
    let receipt = old.bytes.replace(original_hash, &state.hash.to_hex());
    let ordinal = integer(&old.value, "iterations_completed").unwrap();
    let op = ledger.begin_op(Some(spec.id.as_bytes()), &ir_for(SDF3_DRIVER, spec.id, ordinal),
        &FiveExplicits { seed: &spec.seed.to_le_bytes(), versions: "{}", budget: "{}", capability: "{}" }, 0).unwrap();
    for (key, role) in [
        ("source", EdgeRole::In), ("design", EdgeRole::Out), ("iterations", EdgeRole::Out),
        ("report_html", EdgeRole::Out), ("report_json", EdgeRole::Out), ("package", EdgeRole::Out),
    ] {
        let hash = ContentHash::from_hex(old.value.str_field(key).unwrap()).unwrap();
        ledger.link(op, &hash, role).unwrap();
    }
    ledger.link(op, &state.hash, EdgeRole::Out).unwrap();
    let output = ledger.put_artifact(RECEIPT_KIND, receipt.as_bytes(), None).unwrap();
    ledger.link(op, &output.hash, EdgeRole::Out).unwrap();
    ledger.seal_artifact_output(&output.hash, op).unwrap();
    ledger.finish_op(op, OpOutcome::Ok, None, 1).unwrap();
    ledger.commit().unwrap();
    load(ledger, &format!("study-{}", output.hash.to_hex())).unwrap()
}

#[test]
fn g3_one_changed_native_field_bit_cannot_authorize_a_resumed_solve() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = drive(&spec, &ledger, Some(1), &gate(), None, keep).unwrap();
    let old = load(&ledger, &first.pointer).unwrap();
    let forged = forged_state(&ledger, &spec, &old);
    expected(&spec, &ledger, &forged, producer_identity(&gate()).unwrap()).unwrap();
    let error = drive(&spec, &ledger, None, &gate(), Some(&forged), |_, _| {
        panic!("a forged native field must not publish a successor")
    }).unwrap_err();
    assert_eq!(error.code, "cli-study-sdf3-replay-mismatch");
    assert_eq!(load(&ledger, &first.pointer).unwrap().bytes, old.bytes);
}

#[test]
fn g3_g4_admission_rejects_changed_producer_source_and_precancellation() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let first = drive(&spec, &ledger, Some(1), &gate(), None, keep).unwrap();
    for key in ["producer", "study_id"] {
        let mut old = load(&ledger, &first.pointer).unwrap();
        let replacement = JsonValue::parse(&quoted(&ContentHash([7; 32]).to_hex())).unwrap();
        set(&mut old.value, key, replacement);
        let error = drive(&spec, &ledger, None, &gate(), Some(&old), |_, _| {
            panic!("admission failure must not publish")
        }).unwrap_err();
        assert_eq!(error.code, "cli-study-sdf3-checkpoint");
    }
    let mut old = load(&ledger, &first.pointer).unwrap();
    set(&mut old.value, "resume_supported", JsonValue::Bool(false));
    assert_eq!(drive(&spec, &ledger, None, &gate(), Some(&old), keep).unwrap_err().code,
        "cli-study-sdf3-resume-unsupported");
    let cancelled = gate();
    cancelled.request();
    assert_eq!(drive(&spec, &ledger, None, &cancelled, None, |_, _| {
        panic!("pre-cancellation must not publish")
    }).unwrap_err().exit, exit::CANCELLED);
    assert_eq!(load(&ledger, &first.pointer).unwrap().bytes, first.receipt);
}

#[test]
fn g4_every_work_counter_is_retained_and_overflow_cannot_reset_a_budget() {
    let spec = spec();
    let spent = Spent {
        wall_s: 0.5,
        linear: SolveWork { linear_iterations: 7, linear_solves: 3,
            preconditioner_operator_applications: 11, preconditioner_galerkin_products: 13 },
        geometry: QuadratureWork3 { boxes: 17, points: 19, field_evaluations: 23 },
    };
    let json = JsonValue::parse(&format!("{{\"consumed_wall_s\":0.5,\"work\":{}}}", spent.json())).unwrap();
    assert_eq!(Spent::read(&json).unwrap(), spent);
    let twice = spent.add(spent.linear, spent.geometry, spent.wall_s).unwrap();
    assert_eq!(twice.linear.preconditioner_galerkin_products, 26);
    assert_eq!(twice.geometry.field_evaluations, 46);
    twice.require_remaining(&spec).unwrap();
    for exhausted in [
        Spent { wall_s: spec.wall_s, ..spent },
        Spent { linear: SolveWork { linear_iterations: spec.linear, ..spent.linear }, ..spent },
        Spent { geometry: QuadratureWork3 { boxes: spec.boxes, ..spent.geometry }, ..spent },
        Spent { geometry: QuadratureWork3 { points: spec.points, ..spent.geometry }, ..spent },
    ] { assert_eq!(exhausted.require_remaining(&spec).unwrap_err().exit, exit::BUDGET); }
    for invalid in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(spent.add(SolveWork::default(), QuadratureWork3::default(), invalid).is_err());
    }
    let overflow = SolveWork { preconditioner_galerkin_products: usize::MAX, ..Default::default() };
    assert!(spent.add(overflow, QuadratureWork3::default(), 0.0).is_err());
    let overflow = QuadratureWork3 { field_evaluations: usize::MAX, ..Default::default() };
    assert!(spent.add(SolveWork::default(), overflow, 0.0).is_err());
}
