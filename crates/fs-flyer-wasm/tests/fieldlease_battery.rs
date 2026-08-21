//! E7.1-ii battery (bead wf-root-guzez.8.1.2): sample_field over the
//! E5.0 snapshot LEASE. Staleness labeling at the threshold AND one
//! past it (labeled HISTORICAL under its OWN unchanged id, never
//! mislabeled current); the raw-tick-as-state refusal EXECUTED
//! (Round-4 S-02); forbidden-claims propagation through the lease
//! layer; payload-derived bound system (airborne publishes it,
//! on-rail does not); JSON export bounds at cap AND cap+1; empty-ring
//! refusal propagation; determinism.
//! Repro: cd crates/fs-flyer-wasm && cargo test --test fieldlease_battery

use fs_flyer::fieldsvc::{
    C_BOUND_CIRCULATION, C_GROUND_IMAGES, C_MEAN_ATMO, C_TURB_ATMO, GridSpec, claim_total_flow,
};
use fs_flyer_wasm::fieldlease::{
    MAX_JSON_POINTS, STALENESS_THRESHOLD_TICKS, field_selftest_json, leased_sample_json,
    sample_field_by_tick, sample_field_leased,
};
use fs_flyer_wasm::ring::SnapshotRing;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-wasm-fieldlease\",\"case\":\"{case}\",{payload}}}");
}

const AIRBORNE: [f64; 12] = [
    50.0, 3.0, 12.0, 0.8, 0.01, 0.05, 0.02, 0.0, 35.0, 0.0, 0.0, 1.0,
];
const ON_RAIL: [f64; 12] = [
    10.0, 0.1, 5.0, 0.0, 0.0, 0.02, 0.01, 0.0, 30.0, 0.0, 0.0, 0.0,
];

fn ring_with(tick: u64, payload: &[f64; 12]) -> SnapshotRing {
    let mut r = SnapshotRing::new(3, 12, 0x1903, 0xd17).unwrap();
    r.publish(tick, payload).unwrap();
    r
}

fn grid() -> GridSpec {
    GridSpec {
        origin_m: [45.0, -4.0, 0.5],
        dx_m: 2.0,
        nx: 3,
        ny: 3,
        nz: 2,
    }
}

#[test]
fn staleness_labeled_at_threshold_and_beyond() {
    let mut ring = ring_with(100, &AIRBORNE);
    let mask = C_MEAN_ATMO | C_TURB_ATMO;
    // Fresh.
    let fresh = sample_field_leased(&mut ring, 110, 1903, 8.0, 1.294, &grid(), mask).unwrap();
    assert!(!fresh.historical);
    assert_eq!(fresh.staleness_ticks, 10);
    let id = fresh.sample.meta.field_source_snapshot_id.clone();
    // AT the threshold: still current (cap).
    let at = sample_field_leased(
        &mut ring,
        100 + STALENESS_THRESHOLD_TICKS,
        1903,
        8.0,
        1.294,
        &grid(),
        mask,
    )
    .unwrap();
    assert!(!at.historical, "AT the threshold is still current");
    // One past: HISTORICAL, under the SAME id — never relabeled.
    let past = sample_field_leased(
        &mut ring,
        100 + STALENESS_THRESHOLD_TICKS + 1,
        1903,
        8.0,
        1.294,
        &grid(),
        mask,
    )
    .unwrap();
    assert!(past.historical, "one past the threshold is historical");
    assert_eq!(
        past.sample.meta.field_source_snapshot_id, id,
        "a stale snapshot keeps its OWN id"
    );
    assert_eq!(past.sample.meta.source_tick, 100);
    jlog(
        "staleness",
        &format!("\"id\":\"{id}\",\"at_threshold_historical\":false,\"past_historical\":true"),
    );
}

#[test]
fn raw_tick_as_state_is_forbidden() {
    let err = sample_field_by_tick(1234).unwrap_err();
    assert_eq!(err.code, "tick-as-state-forbidden");
    jlog("tick-forbidden", &format!("\"code\":\"{}\"", err.code));
}

#[test]
fn payload_derives_the_bound_system_and_claims_propagate() {
    // Airborne payload: the bound system exists, the full supported
    // sum satisfies the total-flow claim…
    let mut ring = ring_with(100, &AIRBORNE);
    let full = sample_field_leased(
        &mut ring,
        110,
        1903,
        8.0,
        1.294,
        &grid(),
        C_MEAN_ATMO | C_TURB_ATMO | C_BOUND_CIRCULATION | C_GROUND_IMAGES,
    )
    .unwrap();
    claim_total_flow(&full.sample).unwrap();
    // …and an atmosphere-only sum REFUSES the claim through the lease
    // layer (the falsifier, executed).
    let partial = sample_field_leased(
        &mut ring,
        110,
        1903,
        8.0,
        1.294,
        &grid(),
        C_MEAN_ATMO | C_TURB_ATMO,
    )
    .unwrap();
    let err = claim_total_flow(&partial.sample).unwrap_err();
    assert_eq!(err.code, "forbidden-claim-total-flow");
    // On-rail payload: no bound system — requesting it refuses TYPED.
    let mut rail_ring = ring_with(100, &ON_RAIL);
    let err = sample_field_leased(
        &mut rail_ring,
        110,
        1903,
        8.0,
        1.294,
        &grid(),
        C_MEAN_ATMO | C_BOUND_CIRCULATION,
    )
    .unwrap_err();
    assert_eq!(err.code, "component-unsupported");
    jlog("claims", "\"propagated\":true");
}

#[test]
fn json_export_bounds_and_selftest() {
    let mut ring = ring_with(100, &AIRBORNE);
    // AT the JSON cap admits (mean-only keeps it cheap): 4096 points.
    let flat = GridSpec {
        origin_m: [0.0, 0.0, 1.0],
        dx_m: 0.01,
        nx: MAX_JSON_POINTS,
        ny: 1,
        nz: 1,
    };
    let at_cap = sample_field_leased(&mut ring, 110, 1903, 8.0, 1.294, &flat, C_MEAN_ATMO).unwrap();
    let j = leased_sample_json(&at_cap).unwrap();
    assert!(j.contains("\"n_points\":4096"));
    // One more point refuses.
    let over = GridSpec {
        nx: MAX_JSON_POINTS + 1,
        ..flat
    };
    let too_big =
        sample_field_leased(&mut ring, 110, 1903, 8.0, 1.294, &over, C_MEAN_ATMO).unwrap();
    assert_eq!(
        leased_sample_json(&too_big).unwrap_err().code,
        "json-points-exceeded"
    );
    // Whole-path selftest (the wasm-lane surface): well-formed, not a
    // refusal envelope, and carries the honesty fields.
    let st = field_selftest_json();
    assert!(st.contains("\"schema\":\"wf-field-sample-v1\""));
    assert!(!st.contains("refusal"), "selftest refused: {st}");
    assert!(st.contains("\"omitted_components\":["));
    assert!(st.contains("\"historical\":false"));
    // Determinism: the selftest is a pure function.
    assert_eq!(st, field_selftest_json(), "bit-identical twice");
    jlog("json", &format!("\"selftest_len\":{}", st.len()));
}

#[test]
fn empty_ring_refusal_propagates() {
    let mut ring = SnapshotRing::new(3, 12, 0x1903, 0xd17).unwrap();
    let err =
        sample_field_leased(&mut ring, 10, 1903, 8.0, 1.294, &grid(), C_MEAN_ATMO).unwrap_err();
    // Whatever the ring's typed code is, it must PROPAGATE — the
    // lease layer never fabricates a field without a snapshot.
    assert!(!err.code.is_empty());
    jlog("empty-ring", &format!("\"code\":\"{}\"", err.code));
}
