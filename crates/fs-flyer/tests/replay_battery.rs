//! E3.2-ii replay battery (bead wf-root-guzez.4.2.2): ring-window
//! semantics with caps, record→replay bit-identity, LOCALIZING divergence
//! falsifier (a tampered event is pinpointed to its applied tick), and the
//! frozen InputTraceId law (extent-in-domain, ordinal sensitivity,
//! hand-computed preimage).
//! Repro: cargo test -p fs-flyer --test replay_battery

use fs_flyer::replay::{
    AppliedEvent, InputTrace, MAX_RING_CAPACITY, StateRing, record, run_recorded, verify_replay,
};
use fs_flyer::spine::{Loads, RigidBody, SixDofState};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-replay\",\"case\":\"{case}\",{payload}}}");
}

const DT: f64 = 1.0 / 120.0;

fn body() -> RigidBody {
    RigidBody {
        mass_kg: 340.17,
        inertia_kgm2: [1787.0, 367.4, 1820.9],
    }
}

fn rest() -> SixDofState {
    SixDofState {
        pos_m: [0.0; 3],
        vel_mps: [0.0; 3],
        quat: [1.0, 0.0, 0.0, 0.0],
        omega_body: [0.0; 3],
    }
}

// Control-driven loads: channel 0 = pitch-moment command, channel 1 =
// thrust command (a deterministic stand-in for the real aero model).
fn loads(t: f64, s: &SixDofState, c: &[f64]) -> Loads {
    Loads {
        force_n: [
            200.0 * c[1] + 5.0 * (2.0 * t).sin(),
            0.0,
            3336.0 - 30.0 * s.vel_mps[2],
        ],
        moment_nm: [0.0, 400.0 * c[0] - 50.0 * s.omega_body[1], 0.0],
    }
}

fn demo_trace() -> InputTrace {
    InputTrace {
        end_tick_exclusive: 360,
        events: vec![
            AppliedEvent {
                channel: 1,
                applied_tick: 0,
                ordinal_within_tick: 0,
                quantized_value: 0.75,
            },
            AppliedEvent {
                channel: 0,
                applied_tick: 120,
                ordinal_within_tick: 0,
                quantized_value: 0.25,
            },
            AppliedEvent {
                channel: 1,
                applied_tick: 120,
                ordinal_within_tick: 1,
                quantized_value: 1.0,
            },
            AppliedEvent {
                channel: 0,
                applied_tick: 240,
                ordinal_within_tick: 0,
                quantized_value: -0.25,
            },
        ],
    }
}

#[test]
fn ring_window_semantics_and_caps() {
    let mut ring = StateRing::new(8).unwrap();
    for tick in 0..20u32 {
        let mut s = rest();
        s.pos_m[0] = f64::from(tick);
        ring.push(tick, s);
    }
    // Live window is the trailing 8 ticks: 12..=19.
    assert!(ring.get(11).is_none(), "evicted tick must be gone");
    assert!(ring.get(12).is_some() && ring.get(19).is_some());
    assert_eq!(ring.get(15).unwrap().pos_m[0], 15.0);
    assert_eq!(ring.pushes(), 20);
    // Capacity caps at cap AND cap+1 (plus zero).
    assert!(StateRing::new(MAX_RING_CAPACITY).is_ok());
    assert_eq!(
        StateRing::new(MAX_RING_CAPACITY + 1).unwrap_err().code,
        "ring-capacity-invalid"
    );
    assert_eq!(StateRing::new(0).unwrap_err().code, "ring-capacity-invalid");
    jlog(
        "ring",
        "\"window\":\"trailing-8 verified, caps at cap/cap+1\"",
    );
}

#[test]
fn record_replay_bit_identity() {
    let rec = record(&body(), &rest(), DT, demo_trace(), loads).unwrap();
    assert_eq!(rec.tick_digests.len(), 360);
    verify_replay(&rec, loads).expect("same-artifact replay must be bit-identical");
    jlog(
        "bit-identity",
        &format!("\"ticks\":360,\"trace_id\":\"{}\"", rec.input_trace_id),
    );
}

#[test]
fn divergence_localizes_to_the_tampered_tick() {
    let rec = record(&body(), &rest(), DT, demo_trace(), loads).unwrap();
    // FALSIFIER: tamper one mid-run event (tick 240) AND recompute the trace
    // id so the tamper is not caught by the id gate — the digest walk must
    // localize the divergence to exactly tick 240.
    let mut tampered = rec.clone();
    tampered.trace.events[3].quantized_value = -0.26;
    tampered.input_trace_id = tampered.trace.trace_id();
    let refusal = verify_replay(&tampered, loads).unwrap_err();
    assert_eq!(refusal.code, "replay-digest-mismatch");
    assert!(
        refusal.message.contains("first divergence at tick 240"),
        "must localize to tick 240: {}",
        refusal.message
    );
    // And the id gate alone catches silent trace edits.
    let mut silent = rec.clone();
    silent.trace.events[3].quantized_value = -0.26;
    assert_eq!(
        verify_replay(&silent, loads).unwrap_err().code,
        "replay-trace-id-mismatch"
    );
    jlog("localization", "\"tampered_tick\":240,\"localized\":true");
}

#[test]
fn input_trace_id_frozen_law() {
    let trace = demo_trace();
    let id = trace.trace_id();
    // Deterministic and hand-verifiable: recompute the preimage here,
    // independently, per the frozen formula.
    let mut payload = Vec::new();
    payload.extend_from_slice(b"org.frankensim.wright-flyer.input-trace.v1");
    payload.extend_from_slice(&360u32.to_le_bytes());
    for e in &trace.events {
        payload.extend_from_slice(&e.channel.to_le_bytes());
        payload.extend_from_slice(&e.applied_tick.to_le_bytes());
        payload.extend_from_slice(&e.ordinal_within_tick.to_le_bytes());
        payload.extend_from_slice(&e.quantized_value.to_bits().to_le_bytes());
    }
    let hand = fs_blake3::hash_domain("fs-flyer/applied-input-trace/v1", &payload).to_hex();
    assert_eq!(
        id, hand,
        "trace id must match the hand-computed frozen preimage"
    );
    // EXTENT law: identical events, different end_tick_exclusive → different
    // id (two event-free runs stopped at different ticks differ).
    let mut longer = trace.clone();
    longer.end_tick_exclusive = 361;
    assert_ne!(
        id,
        longer.trace_id(),
        "the trace EXTENT is part of the domain"
    );
    // Ordinal sensitivity: swapping same-tick ordinals changes the id.
    let mut swapped = trace.clone();
    swapped.events.swap(1, 2);
    assert_ne!(id, swapped.trace_id());
    jlog("trace-id", &format!("\"id\":\"{id}\""));
}

#[test]
fn trace_admission_refusals() {
    // Event outside the extent.
    let mut out = demo_trace();
    out.events.push(AppliedEvent {
        channel: 0,
        applied_tick: 360,
        ordinal_within_tick: 0,
        quantized_value: 0.0,
    });
    assert_eq!(
        run_recorded(&body(), &rest(), DT, &out, 8, loads)
            .unwrap_err()
            .code,
        "trace-extent-invalid"
    );
    // Broken canonical order (regressing tick).
    let mut disorder = demo_trace();
    disorder.events.swap(1, 3);
    assert_eq!(
        run_recorded(&body(), &rest(), DT, &disorder, 8, loads)
            .unwrap_err()
            .code,
        "trace-order-invalid"
    );
    // Channel above the cap.
    let mut wide = demo_trace();
    wide.events[0].channel = 16;
    assert_eq!(
        run_recorded(&body(), &rest(), DT, &wide, 8, loads)
            .unwrap_err()
            .code,
        "channel-outside-domain"
    );
    jlog("admission", "\"gates\":\"extent, order, channel\"");
}
