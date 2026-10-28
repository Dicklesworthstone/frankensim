//! E6.1-ii battery (bead wf-root-guzez.7.1.2): envelope round-trip
//! (bytes -> parse -> bytes bitwise), scrub-to-tick EQUALS
//! uninterrupted execution BITWISE at rail/airborne/checkpoint-edge/
//! terminal ticks, the event index matches the measured explicit catapult
//! fixture (liftoff and terminal), hostile twins (tamper, truncation, trailing
//! bytes, bad terminal code), Human-mode recording refuses (v1 law),
//! interval caps at cap AND cap+1.
//! Repro: cargo test -p fs-flyer --test replayenv_battery --release

use fs_atmo::FlatSiteLogLaw;
use fs_flyer::replayenv::{
    MAX_CHECKPOINT_INTERVAL, ReplayEnvelope, record_replay,
    record_replay_with_distributed_spatial_wind, record_replay_with_spatial_wind,
};
use fs_flyer::simloop::{
    CATAPULT_1904_V1, Phase, PilotMode, ScenarioSpec, SimLoop, SimStateOut, SpatialWindSpec,
    dec17_scenario,
};
use fs_math::det;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e61ii\",\"case\":\"{case}\",{payload}}}");
}

/// Uninterrupted per-tick states (tick -> SimStateOut).
fn truth_states() -> Vec<SimStateOut> {
    let mut sim = SimLoop::init(assisted_fixture()).unwrap();
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

// Explicit test launch: the corrected unassisted Dec-17 reference reaches
// its rail end before liftoff. Keep actual airborne scrub coverage without
// changing the production scenario or tuning its gains/initial attitude.
fn assisted_fixture() -> ScenarioSpec {
    let mut spec = dec17_scenario(1903, PilotMode::FixedControls);
    spec.catapult = Some(CATAPULT_1904_V1);
    spec
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
        && a.p_rad_s.to_bits() == b.p_rad_s.to_bits()
        && a.phi_rad.to_bits() == b.phi_rad.to_bits()
        && a.r_rad_s.to_bits() == b.r_rad_s.to_bits()
        && a.psi_rad.to_bits() == b.psi_rad.to_bits()
        && a.dc_rad.to_bits() == b.dc_rad.to_bits()
        && a.warp_rad.to_bits() == b.warp_rad.to_bits()
        && a.omega_prop_rad_s.to_bits() == b.omega_prop_rad_s.to_bits()
        && a.gust_w_mps.to_bits() == b.gust_w_mps.to_bits()
        && a.assist_active == b.assist_active
        && a.assist_dc_rad.to_bits() == b.assist_dc_rad.to_bits()
}

#[test]
fn scrub_equals_uninterrupted_bitwise_and_events_match() {
    let env = record_replay(assisted_fixture(), 200).unwrap();
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
    for t in [50u64, 74, 75, 200, 401, 600, env.events.terminal_tick] {
        let scrubbed = env.scrub_to_tick(t).unwrap();
        let truth_at = truth.iter().find(|s| s.tick == t).unwrap();
        assert!(
            bitwise_eq(&scrubbed, truth_at),
            "scrub({t}) diverged: {scrubbed:?} vs {truth_at:?}"
        );
    }
    // Edges: terminal and the true initialization admit; one past refuses.
    assert_eq!(
        env.scrub_to_tick(env.events.terminal_tick + 1)
            .unwrap_err()
            .code,
        "replay-scrub-out-of-range"
    );
    let initial = env.scrub_to_tick(0).unwrap();
    assert_eq!(initial.tick, 0);
    assert_eq!(initial.phase, Phase::OnRail);
    assert_eq!(initial.x_m, 0.0);
    assert!(
        env.checkpoints.iter().all(|(tick, _)| *tick > 0),
        "legacy recording keeps its original checkpoint sequence and bytes"
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
    // The lifecycle's explicit damage-model boundary is a recorded terminal,
    // not a malformed envelope (phase_code already emits this value).
    let mut damage_boundary = env.clone();
    damage_boundary.events.terminal_code = 6;
    assert_eq!(
        ReplayEnvelope::from_bytes(&damage_boundary.to_bytes()).unwrap(),
        damage_boundary
    );
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

fn spatial_wind() -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: 2.4,
            reference_speed_mps: 11.0,
        },
        sigma_mps: 0.3,
        length_m: 20.0,
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 8],
    }
}

#[test]
fn spatial_recordings_replay_zero_early_checkpoint_and_terminal_states() {
    for distributed in [false, true] {
        let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
        scenario.max_ticks = 12;
        let mut sim = if distributed {
            SimLoop::init_with_distributed_spatial_wind(scenario.clone(), spatial_wind())
        } else {
            SimLoop::init_with_spatial_wind(scenario.clone(), spatial_wind())
        }
        .unwrap();
        let original_initial = sim.save_checkpoint().unwrap();
        let wind_at_zero = sim.sample_spatial_wind([0.0, 0.0, 2.4]).unwrap();
        let theta = sim.tick0().trim.alpha_rad;
        let canard = sim.tick0().mech.delta_rad;
        let mut truth = Vec::new();
        loop {
            let state = sim.step(None).unwrap();
            truth.push(state);
            if matches!(state.phase, Phase::Ended(_)) {
                break;
            }
        }
        assert!(
            truth.len() >= 9,
            "short wind fixture must reach its later checkpoint"
        );
        let recorded = if distributed {
            record_replay_with_distributed_spatial_wind(scenario, spatial_wind(), 4)
        } else {
            record_replay_with_spatial_wind(scenario, spatial_wind(), 4)
        }
        .unwrap();
        assert_eq!(recorded.checkpoints[0], (0, original_initial));
        let bytes = recorded.to_bytes();
        let replay = ReplayEnvelope::from_bytes(&bytes).unwrap();
        assert_eq!(replay.to_bytes(), bytes);
        assert_eq!(replay.final_digest, sim.digest_hex());
        let zero = replay.scrub_to_tick(0).unwrap();
        assert_eq!(zero.tick, 0);
        assert_eq!(zero.phase, Phase::OnRail);
        assert_eq!((zero.x_m, zero.h_m, zero.q_rad_s), (0.0, 2.4, 0.0));
        assert_eq!(zero.theta_rad.to_bits(), theta.to_bits());
        assert_eq!(zero.dc_rad.to_bits(), canard.to_bits());
        assert_eq!(zero.gust_w_mps.to_bits(), wind_at_zero.u[2].to_bits());
        let u = det::cos(theta) * -wind_at_zero.u[0] + det::sin(theta) * -wind_at_zero.u[2];
        let w = det::sin(theta) * -wind_at_zero.u[0] - det::cos(theta) * -wind_at_zero.u[2];
        assert_eq!(zero.u_mps.to_bits(), u.to_bits());
        assert_eq!(zero.w_mps.to_bits(), w.to_bits());
        assert!(!zero.assist_active);
        for tick in [1, 3, 4, 5, 9, replay.events.terminal_tick] {
            let actual = replay.scrub_to_tick(tick).unwrap();
            assert!(
                bitwise_eq(&actual, &truth[(tick - 1) as usize]),
                "distributed={distributed} scrub={tick}"
            );
        }
    }
}

#[test]
fn missing_or_downgraded_wind_initialization_refuses_even_after_a_checkpoint() {
    let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
    scenario.max_ticks = 8;
    let original =
        record_replay_with_distributed_spatial_wind(scenario.clone(), spatial_wind(), 4).unwrap();
    assert!(original.events.terminal_tick > 4);
    let mut missing = original.clone();
    missing.checkpoints.retain(|(tick, _)| *tick != 0);
    for tick in [0, 1, 5] {
        assert_eq!(
            missing.scrub_to_tick(tick).unwrap_err().code,
            "replay-initial-state-missing"
        );
    }
    let legacy = SimLoop::init(scenario.clone())
        .unwrap()
        .save_checkpoint()
        .unwrap();
    let cg = SimLoop::init_with_spatial_wind(scenario, spatial_wind())
        .unwrap()
        .save_checkpoint()
        .unwrap();
    for wrong_mode in [legacy, cg] {
        let mut downgraded = original.clone();
        downgraded.checkpoints[0].1 = wrong_mode;
        // Re-serialization gives these bytes a valid outer checksum; the
        // original mode/intent mismatch still cannot become a legacy replay.
        let downgraded = ReplayEnvelope::from_bytes(&downgraded.to_bytes()).unwrap();
        for tick in [0, 1, 5] {
            assert_eq!(
                downgraded.scrub_to_tick(tick).unwrap_err().code,
                "replay-identity-mismatch"
            );
        }
    }
    let mut duplicate_zero = original;
    duplicate_zero
        .checkpoints
        .insert(0, duplicate_zero.checkpoints[0].clone());
    assert_eq!(
        ReplayEnvelope::from_bytes(&duplicate_zero.to_bytes())
            .unwrap_err()
            .code,
        "replay-malformed"
    );
}
