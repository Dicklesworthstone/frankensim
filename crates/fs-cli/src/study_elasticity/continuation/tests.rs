use super::*;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-2d.fsim"));

fn spec() -> ElasticitySpec {
    parse(&FIXTURE.replace(":mesh-level 4", ":mesh-level 2")
        .replace(":max-iterations 8", ":max-iterations 3")
        .replace(":steps 8", ":steps 3")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":nucleation-period 4", ":nucleation-period 2"))
        .expect("admitted small native study")
}

fn json(output: &Outcome) -> JsonValue { JsonValue::parse(&output.receipt).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn artifact_bytes(ledger: &Ledger, receipt: &JsonValue, key: &str) -> Vec<u8> {
    let kind = if key == "design" { "study-design" } else { "study-iterations" };
    linked(ledger, receipt, key, kind).unwrap()
}
fn set(value: &mut JsonValue, key: &str, replacement: JsonValue) {
    let JsonValue::Object(members) = value else { panic!("expected object") };
    let slot = members.iter_mut().find(|(name, _)| name == key).unwrap();
    slot.1 = replacement;
}

#[test]
fn native_chunks_restore_exact_design_and_trace_without_prefix_replay() {
    let spec = spec();
    let full_ledger = Ledger::open(":memory:").unwrap();
    let full = drive(&spec, &full_ledger, None, &gate(), None).unwrap();
    assert_eq!(full.status, "completed");
    let full_json = json(&full);
    assert_eq!(full_json.path(&["continuation", "updates_this_invocation"]).and_then(JsonValue::as_f64), Some(3.0));

    let ledger = Ledger::open(":memory:").unwrap();
    let prefix = drive(&spec, &ledger, Some(1), &gate(), None).unwrap();
    assert_eq!(prefix.status, "budget-exhausted");
    let prefix_json = json(&prefix);
    let untouched = artifact_bytes(&ledger, &prefix_json, "design");
    let loaded = load(&ledger, &prefix.pointer).unwrap();
    let resumed = drive(&spec, &ledger, None, &gate(), Some(&loaded)).unwrap();
    assert_eq!(resumed.status, "completed");
    let resumed_json = json(&resumed);
    assert_eq!(resumed_json.path(&["continuation", "updates_this_invocation"]).and_then(JsonValue::as_f64), Some(2.0));
    assert_eq!(resumed_json.path(&["continuation", "legacy_prefix_updates_replayed"]).and_then(JsonValue::as_f64), Some(0.0));
    for key in ["design", "iterations"] {
        assert_eq!(artifact_bytes(&full_ledger, &full_json, key), artifact_bytes(&ledger, &resumed_json, key));
    }
    assert_eq!(full_json.str_field("trace_hash"), resumed_json.str_field("trace_hash"));
    assert_eq!(artifact_bytes(&ledger, &prefix_json, "design"), untouched);
    let oracle = run_prefix(&spec, spec.steps).unwrap();
    let expected = document(&artifact_bytes(&ledger, &resumed_json, "iterations")).unwrap();
    let design = document(&artifact_bytes(&ledger, &resumed_json, "design")).unwrap();
    let decoded = decode(&spec, &resumed_json, &design, &expected).unwrap();
    assert_eq!(decoded.1.rows, oracle.1.rows);
    assert_eq!(decoded.0.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        oracle.0.nodes().iter().map(|v| v.to_bits()).collect::<Vec<_>>());
}

#[test]
fn cancellation_before_first_solve_retains_a_resumable_seed() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let cancelled = gate();
    cancelled.request();
    let output = drive(&spec, &ledger, None, &cancelled, None).unwrap();
    assert_eq!(output.status, "cancelled");
    assert_eq!(integer(&json(&output), "iterations_completed").unwrap(), 0);
    let loaded = load(&ledger, &output.pointer).unwrap();
    let continued = drive(&spec, &ledger, Some(1), &gate(), Some(&loaded)).unwrap();
    assert_eq!(integer(&json(&continued), "iterations_completed").unwrap(), 1);
    assert_eq!(continued.status, "budget-exhausted");
}

#[test]
fn final_accepted_checkpoint_can_finalize_without_another_geometry_update() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let output = drive(&spec, &ledger, None, &gate(), None).unwrap();
    let receipt = json(&output);
    // Simulate losing the final publication: the predecessor is the already
    // committed, numerically complete state, still labelled running.
    let pointer = format!("study-{}", receipt.str_field("predecessor").unwrap());
    let accepted = load(&ledger, &pointer).unwrap();
    assert_eq!(accepted.value.str_field("status"), Some("running"));
    assert_eq!(integer(&accepted.value, "iterations_completed").unwrap(), spec.steps);
    let finalized = drive(&spec, &ledger, None, &gate(), Some(&accepted)).unwrap();
    assert_eq!(finalized.status, "completed");
    assert_eq!(json(&finalized).path(&["continuation", "updates_this_invocation"]).and_then(JsonValue::as_f64), Some(0.0));
    assert_eq!(artifact_bytes(&ledger, &receipt, "design"), artifact_bytes(&ledger, &json(&finalized), "design"));
    let complete = load(&ledger, &finalized.pointer).unwrap();
    let again = drive(&spec, &ledger, None, &gate(), Some(&complete)).unwrap();
    assert_eq!(again.pointer, finalized.pointer);
    assert_eq!(again.receipt, finalized.receipt);
}

#[test]
fn retained_state_decoder_rejects_changed_geometry_multiplier_and_lattice() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let output = drive(&spec, &ledger, Some(1), &gate(), None).unwrap();
    let receipt = json(&output);
    let mut design = document(&artifact_bytes(&ledger, &receipt, "design")).unwrap();
    let iterations = document(&artifact_bytes(&ledger, &receipt, "iterations")).unwrap();
    let JsonValue::Object(ref mut members) = design else { panic!("design object") };
    let (_, JsonValue::Array(bits)) = members.iter_mut().find(|(name, _)| name == "phi_bits").unwrap()
        else { panic!("bits") };
    bits[0] = JsonValue::Str("7ff8000000000000".into());
    assert!(decode(&spec, &receipt, &design, &iterations).is_err());

    let mut design = document(&artifact_bytes(&ledger, &receipt, "design")).unwrap();
    set(&mut design, "n", JsonValue::Number { value: 1e9, raw: "1000000000".into() });
    assert!(decode(&spec, &receipt, &design, &iterations).is_err());

    let design = document(&artifact_bytes(&ledger, &receipt, "design")).unwrap();
    let mut iterations = iterations;
    let JsonValue::Object(ref mut members) = iterations else { panic!("iterations object") };
    let (_, JsonValue::Array(rows)) = members.iter_mut().find(|(name, _)| name == "iterations").unwrap()
        else { panic!("rows") };
    set(&mut rows[0], "ell", JsonValue::Number { value: 9876.0, raw: "9876".into() });
    assert!(decode(&spec, &receipt, &design, &iterations).is_err());
}

#[test]
fn changed_executable_binding_refuses_without_extending_the_study() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let output = drive(&spec, &ledger, Some(1), &gate(), None).unwrap();
    let mut loaded = load(&ledger, &output.pointer).unwrap();
    // Exercise the internal continuation admission, not the public seal loader.
    let JsonValue::Object(ref mut members) = loaded.value else { panic!("receipt object") };
    let (_, binding) = members.iter_mut().find(|(name, _)| name == "continuation").unwrap();
    set(binding, "producer", JsonValue::Str("00".repeat(32)));
    let error = drive(&spec, &ledger, None, &gate(), Some(&loaded)).unwrap_err();
    assert!(error.message.contains("identical executable"));
    assert_eq!(integer(&load(&ledger, &output.pointer).unwrap().value, "iterations_completed").unwrap(), 1);
}
