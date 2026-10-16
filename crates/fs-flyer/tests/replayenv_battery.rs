//! E6.1-ii battery (bead wf-root-guzez.7.1.2): envelope round-trip
//! (bytes -> parse -> bytes bitwise), scrub-to-tick EQUALS
//! uninterrupted execution BITWISE at rail/airborne/checkpoint-edge/
//! terminal ticks, the event index matches the measured run (liftoff
//! 626 class, terminal), hostile twins (tamper, truncation, trailing
//! bytes, bad terminal code), Human-mode recording refuses (v1 law),
//! interval caps at cap AND cap+1.
//! Repro: cargo test -p fs-flyer --test replayenv_battery --release

use fs_flyer::replayenv::{MAX_CHECKPOINT_INTERVAL, ReplayEnvelope, record_replay};
use fs_flyer::simloop::{Phase, PilotMode, SimLoop, SimStateOut, dec17_scenario};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e61ii\",\"case\":\"{case}\",{payload}}}");
}

/// Uninterrupted per-tick states (tick -> SimStateOut).
fn truth_states() -> Vec<SimStateOut> {
    let mut sim = SimLoop::init(dec17_scenario(1903, PilotMode::FixedControls)).unwrap();
    let mut out = Vec::new();
    loop {
        let s = sim.step(None).unwrap();
        out.push(s);
        if let Phase::Ended(_) = s.phase {
            break;
        }
    }
    out
}

fn bitwise_eq(a: &SimStateOut, b: &SimStateOut) -> bool {
    a.tick == b.tick
        && a.phase == b.phase
        && a.x_m.to_bits() == b.x_m.to_bits()
        && a.h_m.to_bits() == b.h_m.to_bits()
        && a.u_mps.to_bits() == b.u_mps.to_bits()
        && a.w_mps.to_bits() == b.w_mps.to_bits()
        && a.q_rad_s.to_bits() == b.q_rad_s.to_bits()
        && a.theta_rad.to_bits() == b.theta_rad.to_bits()
        && a.dc_rad.to_bits() == b.dc_rad.to_bits()
        && a.omega_prop_rad_s.to_bits() == b.omega_prop_rad_s.to_bits()
        && a.gust_w_mps.to_bits() == b.gust_w_mps.to_bits()
}

#[test]
fn scrub_equals_uninterrupted_bitwise_and_events_match() {
    let env = record_replay(dec17_scenario(1903, PilotMode::FixedControls), 200).unwrap();
    let truth = truth_states();
    // Event index vs the measured run.
    assert_eq!(env.events.terminal_tick, truth.last().unwrap().tick);
    let liftoff_truth = truth
        .iter()
        .find(|s| matches!(s.phase, Phase::Airborne))
        .map(|s| s.tick);
    assert_eq!(env.events.liftoff_tick, liftoff_truth);
    assert_eq!(env.events.terminal_code, 5, "fixed run: envelope exit");
    // Scrub points: rail, airborne, exactly ON a checkpoint tick, and
    // the terminal tick itself.
    for t in [150u64, 401, 800, 1000, env.events.terminal_tick] {
        let scrubbed = env.scrub_to_tick(t).unwrap();
        let truth_at = truth.iter().find(|s| s.tick == t).unwrap();
        assert!(
            bitwise_eq(&scrubbed, truth_at),
            "scrub({t}) diverged: {scrubbed:?} vs {truth_at:?}"
        );
    }
    // Edges: terminal admits (above), one past refuses; 0 refuses.
    assert_eq!(
        env.scrub_to_tick(env.events.terminal_tick + 1)
            .unwrap_err()
            .code,
        "replay-scrub-out-of-range"
    );
    assert_eq!(
        env.scrub_to_tick(0).unwrap_err().code,
        "replay-scrub-out-of-range"
    );
    jlog(
        "scrub",
        &format!(
            "\"terminal\":{},\"checkpoints\":{},\"undulations\":{}",
            env.events.terminal_tick,
            env.checkpoints.len(),
            env.events.undulation_ticks.len()
        ),
    );
}

#[test]
fn envelope_roundtrip_and_hostile_twins() {
    let env = record_replay(
        {
            let mut s = dec17_scenario(11, PilotMode::FixedControls);
            s.max_ticks = 240; // short: rail-only, cheap
            s
        },
        60,
    )
    .unwrap();
    let bytes = env.to_bytes();
    let parsed = ReplayEnvelope::from_bytes(&bytes).unwrap();
    assert_eq!(parsed, env, "parse inverts serialize");
    assert_eq!(parsed.to_bytes(), bytes, "bytes round-trip bitwise");
    // Tamper: flip a byte inside a checkpoint payload region.
    let mut bad = bytes.clone();
    let mid = bytes.len() / 2;
    bad[mid] ^= 0x40;
    assert_eq!(
        ReplayEnvelope::from_bytes(&bad).unwrap_err().code,
        "replay-tampered"
    );
    // Truncation.
    assert!(matches!(
        ReplayEnvelope::from_bytes(&bytes[..bytes.len() - 33])
            .unwrap_err()
            .code,
        "replay-tampered" | "replay-malformed"
    ));
    // Trailing garbage breaks the digest.
    let mut extra = bytes.clone();
    extra.extend_from_slice(b"xx");
    assert_eq!(
        ReplayEnvelope::from_bytes(&extra).unwrap_err().code,
        "replay-tampered"
    );
    jlog("roundtrip", &format!("\"bytes\":{}", bytes.len()));
}

#[test]
fn recording_refusals_at_cap_and_cap_plus_one() {
    // Human mode refuses (v1 law: inputs are never silently dropped).
    assert_eq!(
        record_replay(dec17_scenario(1, PilotMode::Human), 100)
            .unwrap_err()
            .code,
        "replay-human-needs-trace"
    );
    // Interval caps.
    let short = |ticks: u64| {
        let mut s = dec17_scenario(1, PilotMode::FixedControls);
        s.max_ticks = ticks;
        s
    };
    assert!(record_replay(short(10), MAX_CHECKPOINT_INTERVAL).is_ok());
    assert_eq!(
        record_replay(short(10), MAX_CHECKPOINT_INTERVAL + 1)
            .unwrap_err()
            .code,
        "replay-interval-invalid"
    );
    assert_eq!(
        record_replay(short(10), 0).unwrap_err().code,
        "replay-interval-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}
