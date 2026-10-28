//! E5.1 wasm-surface battery (bead wf-root-guzez.6.2): the envelope
//! API over the REAL engine, exercised natively (the wasm exports are
//! 1:1 wrappers over `EngineSlot`). Every documented refusal code at
//! this surface is EXECUTED; the short lifecycle is run twice
//! bit-identically; envelopes are well-formed for ok AND refusal.
//! Repro: cd crates/fs-flyer-wasm && cargo test --test engine_battery

use fs_flyer::simloop::{PilotMode, ScenarioSpec, SimLoop, SpatialWindSpec};
use fs_flyer_wasm::engine::{EngineSlot, MODE_FIXED, MODE_HISTORICAL, MODE_HUMAN};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-wasm-e51\",\"case\":\"{case}\",{payload}}}");
}

const DEC17: (f64, f64, f64, u64) = (1.294, 11.0, 18.3, 40);

fn init_short(slot: &mut EngineSlot, seed: u64, mode: u32, member: u32) -> String {
    let (rho, wind, rail, ticks) = DEC17;
    slot.init(seed, rho, wind, mode, member, rail, ticks, false, false)
}

fn explicit_wind() -> SpatialWindSpec {
    SpatialWindSpec {
        mean: fs_atmo::FlatSiteLogLaw {
            scenario_effective_z0_m: 0.004,
            displacement_height_m: 0.02,
            reference_height_m: 2.6,
            reference_speed_mps: 11.0,
        },
        sigma_mps: 0.08,
        length_m: 18.0,
        advection_mps: 9.5,
        correlation_times_s: vec![0.8, 1.1, 1.7, 2.4, 3.0, 3.5, 4.0, 4.8],
    }
}

fn init_spatial(slot: &mut EngineSlot, distributed: bool, wind: &SpatialWindSpec) -> String {
    let (rho, headwind, rail, ticks) = DEC17;
    slot.init_with_spatial_wind(
        1903,
        rho,
        headwind,
        MODE_FIXED,
        0,
        rail,
        ticks,
        false,
        false,
        distributed,
        wind.mean.scenario_effective_z0_m,
        wind.mean.displacement_height_m,
        wind.mean.reference_height_m,
        wind.sigma_mps,
        wind.length_m,
        wind.advection_mps,
        &wind.correlation_times_s,
    )
}

#[test]
fn spatial_initializer_matches_native_modes_and_live_checkpoints() {
    let mut identities = Vec::new();
    for distributed in [false, true] {
        let wind = explicit_wind();
        let mut slot = EngineSlot::default();
        let envelope = init_spatial(&mut slot, distributed, &wind);
        assert!(envelope.starts_with("{\"ok\":{"), "{envelope}");
        let (rho, headwind, rail, ticks) = DEC17;
        let scenario = ScenarioSpec {
            seed: 1903,
            rho_kg_m3: rho,
            headwind_mps: headwind,
            pilot_mode: PilotMode::FixedControls,
            rail_length_m: rail,
            max_ticks: ticks,
            assist: None,
            catapult: None,
        };
        let mut native = if distributed {
            SimLoop::init_with_distributed_spatial_wind(scenario, wind)
        } else {
            SimLoop::init_with_spatial_wind(scenario, wind)
        }
        .unwrap();
        assert!(envelope.contains(&native.run_intent_id));
        assert!(envelope.contains(&native.tick0().digest));
        identities.push(native.run_intent_id.clone());
        for tick in 0..=4 {
            if tick > 0 {
                let step = slot.step(false, 0.0, 0.0);
                assert!(step.starts_with("{\"ok\":{"), "{step}");
                native.step(None).unwrap();
            }
            let bytes = native.save_checkpoint().unwrap();
            assert_eq!(
                u32::from_le_bytes(bytes[..4].try_into().unwrap()),
                if distributed { 4 } else { 3 }
            );
            let hex: String = bytes.iter().map(|b| format!("{b:02x}")).collect();
            assert_eq!(
                slot.checkpoint(),
                format!(
                    "{{\"ok\":{{\"run_intent_id\":\"{}\",\"checkpoint_hex\":\"{hex}\"}}}}",
                    native.run_intent_id
                )
            );
            assert_eq!(
                slot.digest(),
                format!("{{\"ok\":{{\"digest\":\"{}\"}}}}", native.digest_hex())
            );
        }
    }
    assert_ne!(identities[0], identities[1]);
}

#[test]
fn spatial_admission_refuses_without_legacy_fallback_or_stale_run() {
    let mut slot = EngineSlot::default();
    assert!(init_short(&mut slot, 1903, MODE_FIXED, 0).starts_with("{\"ok\":{"));
    for case in 0..6 {
        let mut wind = explicit_wind();
        let code = match case {
            0 => {
                wind.mean.scenario_effective_z0_m = 0.0;
                "z0-outside-domain"
            }
            1 => {
                wind.correlation_times_s.clear();
                "mode-count-invalid"
            }
            2 => {
                wind.correlation_times_s[0] = 0.0;
                "ou-params-invalid"
            }
            3 => {
                wind.sigma_mps = f64::NAN;
                "non-finite-input"
            }
            4 => {
                wind.correlation_times_s = vec![1.0; fs_atmo::MAX_MODES + 1];
                "mode-count-invalid"
            }
            _ => {
                wind.length_m = -1.0;
                "turbulence-params-invalid"
            }
        };
        let refusal = init_spatial(&mut slot, case % 2 == 0, &wind);
        assert!(
            refusal.contains(&format!("\"code\":\"{code}\"")),
            "{refusal}"
        );
        assert!(
            slot.checkpoint()
                .contains("\"code\":\"engine-not-initialized\"")
        );
        assert!(
            slot.step(false, 0.0, 0.0)
                .contains("\"code\":\"engine-not-initialized\"")
        );
    }
}

#[test]
fn init_envelope_carries_run_identity_and_trim() {
    let mut slot = EngineSlot::default();
    let env = init_short(&mut slot, 1903, MODE_FIXED, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    for key in [
        "run_intent_id",
        "tick0_digest",
        "trim_v_mps",
        "trim_alpha_rad",
        "trim_dc_rad",
        "trim_omega_rad_s",
    ] {
        assert!(env.contains(key), "missing {key}: {env}");
    }
    jlog("init", &format!("\"len\":{}", env.len()));
}

#[test]
fn step_envelope_walks_the_rail_and_max_ticks_terminal() {
    let mut slot = EngineSlot::default();
    let env = init_short(&mut slot, 1903, MODE_FIXED, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let mut last = String::new();
    for _ in 0..40 {
        last = slot.step(false, 0.0, 0.0);
        assert!(last.starts_with("{\"ok\":{"), "{last}");
    }
    // 40 ticks on an 18.3 m rail: still constrained, then the tick
    // budget (40) ends the run in the SAME step envelope.
    assert!(
        last.contains("\"phase\":\"ended:max-ticks\""),
        "terminal in-band: {last}"
    );
    // Stepping past the terminal is the typed refusal.
    let past = slot.step(false, 0.0, 0.0);
    assert!(past.contains("\"code\":\"run-ended\""), "{past}");
    jlog("lifecycle", &format!("\"last\":{}", last.len()));
}

#[test]
fn every_documented_surface_refusal_code_executes() {
    // engine-not-initialized (step and digest).
    let mut fresh = EngineSlot::default();
    let e1 = fresh.step(false, 0.0, 0.0);
    assert!(e1.contains("\"code\":\"engine-not-initialized\""), "{e1}");
    let e2 = fresh.digest();
    assert!(e2.contains("\"code\":\"engine-not-initialized\""), "{e2}");
    let e_checkpoint = fresh.checkpoint();
    assert!(
        e_checkpoint.contains("\"code\":\"engine-not-initialized\""),
        "{e_checkpoint}"
    );
    // mode-invalid (3 is one past the last mode word — cap AND cap+1).
    let e3 = fresh.init(1, 1.294, 11.0, 3, 0, 18.3, 40, false, false);
    assert!(e3.contains("\"code\":\"mode-invalid\""), "{e3}");
    // scenario-invalid (headwind past the cap; refuses before any
    // heavy equilibration work).
    let e4 = fresh.init(
        1,
        1.294,
        20.0_f64.next_up(),
        MODE_FIXED,
        0,
        18.3,
        40,
        false,
        false,
    );
    assert!(e4.contains("\"code\":\"scenario-invalid\""), "{e4}");
    // control-input-missing: Human mode without input;
    // control-input-invalid: non-finite lever force.
    let mut human = EngineSlot::default();
    let env = init_short(&mut human, 7, MODE_HUMAN, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let e5 = human.step(false, 0.0, 0.0);
    assert!(e5.contains("\"code\":\"control-input-missing\""), "{e5}");
    let e6 = human.step(true, f64::NAN, 0.0);
    assert!(e6.contains("\"code\":\"control-input-invalid\""), "{e6}");
    // run-ended is executed in the lifecycle test above; execute it
    // here too so THIS test alone covers the documented list.
    let mut ended = EngineSlot::default();
    let (rho, wind, rail, _) = DEC17;
    let env = ended.init(1, rho, wind, MODE_FIXED, 0, rail, 2, false, false);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    ended.step(false, 0.0, 0.0);
    ended.step(false, 0.0, 0.0);
    let e7 = ended.step(false, 0.0, 0.0);
    assert!(e7.contains("\"code\":\"run-ended\""), "{e7}");
    jlog("refusals", "\"codes\":6");
}

#[test]
fn checkpoint_envelope_carries_live_simloop_bytes_and_refuses_after_terminal() {
    let mut slot = EngineSlot::default();
    let env = init_short(&mut slot, 1903, MODE_FIXED, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    slot.step(false, 0.0, 0.0);
    let checkpoint = slot.checkpoint();
    assert!(checkpoint.starts_with("{\"ok\":{"), "{checkpoint}");
    assert!(checkpoint.contains("\"run_intent_id\":\""), "{checkpoint}");
    let hex_prefix = "\"checkpoint_hex\":\"";
    let hex_start = checkpoint.find(hex_prefix).expect("checkpoint hex field") + hex_prefix.len();
    let hex = &checkpoint[hex_start..checkpoint.len() - 3];
    assert!(!hex.is_empty() && hex.len() % 2 == 0, "{checkpoint}");
    assert!(
        hex.bytes().all(|byte| byte.is_ascii_hexdigit()),
        "{checkpoint}"
    );

    let mut terminal = EngineSlot::default();
    let (rho, wind, rail, _) = DEC17;
    let init = terminal.init(1903, rho, wind, MODE_FIXED, 0, rail, 1, false, false);
    assert!(init.starts_with("{\"ok\":{"), "{init}");
    terminal.step(false, 0.0, 0.0);
    let after_terminal = terminal.checkpoint();
    assert!(
        after_terminal.contains("\"code\":\"checkpoint-after-terminal\""),
        "{after_terminal}"
    );
    jlog("checkpoint", &format!("\"bytes_hex\":{}", hex.len()));
}

#[test]
fn failed_reinit_invalidates_the_previous_lifecycle_before_refusing() {
    let mut slot = EngineSlot::default();
    let valid = init_short(&mut slot, 1903, MODE_FIXED, 0);
    assert!(valid.starts_with("{\"ok\":{"), "{valid}");
    assert!(slot.checkpoint().starts_with("{\"ok\":{"));

    // Mode 3 is deliberately invalid. The former implementation kept the
    // valid A lifecycle live when this B request refused.
    let refused = slot.init(1904, 1.294, 11.0, 3, 0, 18.3, 40, false, false);
    assert!(refused.contains("\"code\":\"mode-invalid\""), "{refused}");
    let stale = slot.checkpoint();
    assert!(
        stale.contains("\"code\":\"engine-not-initialized\""),
        "failed B must not expose A checkpoint: {stale}"
    );
}

#[test]
fn human_mode_steps_with_input_and_historical_member_binds() {
    let mut human = EngineSlot::default();
    let env = init_short(&mut human, 7, MODE_HUMAN, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let s = human.step(true, 25.0, 0.01);
    assert!(s.starts_with("{\"ok\":{"), "{s}");
    assert!(s.contains("\"phase\":\"on-rail\""), "{s}");
    for key in ["p_rad_s", "phi_rad", "r_rad_s", "psi_rad"] {
        assert!(s.contains(key), "missing lateral field {key}: {s}");
    }
    // Historical member is intent-bearing: run_intent_id differs.
    let mut h0 = EngineSlot::default();
    let mut h1 = EngineSlot::default();
    let e0 = init_short(&mut h0, 7, MODE_HISTORICAL, 0);
    let e1 = init_short(&mut h1, 7, MODE_HISTORICAL, 1);
    let intent = |e: &str| {
        let k = "\"run_intent_id\":\"";
        let i = e.find(k).unwrap() + k.len();
        e[i..i + 32].to_string()
    };
    assert_ne!(intent(&e0), intent(&e1), "member must bind into intent");
    jlog("modes", "\"human_and_member_binding\":true");
}

#[test]
fn assist_visibility_and_counterfactual_receipt() {
    // E5.3c (bead guzez.6.6): assist ON — the envelope reports the
    // active flag AND the bounded contribution; assist OFF over the
    // SAME scenario is the model counterfactual and its chained digest
    // DIFFERS (on this short rail-only run the divergence enters via
    // the flagged snapshot: the assist flag is run identity, honestly
    // labeled; the airborne physics counterfactual is the app-level
    // replay path, E6.1's ABComparisonReceipt scope).
    let (rho, wind, rail, _) = DEC17;
    let mut on = EngineSlot::default();
    let env = on.init(1903, rho, wind, MODE_FIXED, 0, rail, 30, true, false);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let mut saw_active = false;
    for _ in 0..30 {
        let s = on.step(false, 0.0, 0.0);
        assert!(s.starts_with("{\"ok\":{"), "{s}");
        if s.contains("\"assist_active\":true") {
            saw_active = true;
            // Authority bound: |dc_assist| <= 0.3 * stop (0.5236) rad.
            let key = "\"assist_dc_rad\":";
            let i = s.find(key).unwrap() + key.len();
            let end = s[i..].find([',', '}']).unwrap() + i;
            let v: f64 = s[i..end].parse().unwrap();
            assert!(v.abs() <= 0.3 * 0.5236 + 1e-12, "authority bound: {v}");
        }
    }
    assert!(saw_active, "assist visibility: the flag must surface");
    let d_on = on.digest();
    let mut off = EngineSlot::default();
    let env = off.init(1903, rho, wind, MODE_FIXED, 0, rail, 30, false, false);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    for _ in 0..30 {
        let s = off.step(false, 0.0, 0.0);
        assert!(s.contains("\"assist_active\":false"), "{s}");
        assert!(
            s.contains("\"assist_dc_rad\":0"),
            "no-assist reports zero: {s}"
        );
    }
    let d_off = off.digest();
    assert_ne!(d_on, d_off, "the counterfactual MUST diverge");
    jlog(
        "assist-counterfactual",
        "\"label\":\"model counterfactual\",\"diverges\":true",
    );
}

#[test]
fn selftest_matches_its_lane_golden_and_the_falsifier_fires() {
    // The shipped self-test (per-lane golden) passes on this build...
    let ok = fs_flyer_wasm::engine::selftest();
    assert!(ok.contains("\"matched\":true"), "{ok}");
    // ...and the DELIBERATELY PERTURBED golden fires the determinism
    // badge refusal (the DONE-WHEN falsifier, executed).
    let bad = fs_flyer_wasm::engine::selftest_against(&"0".repeat(64));
    assert!(
        bad.contains("\"code\":\"determinism-selftest-failed\""),
        "{bad}"
    );
    jlog("selftest", "\"lane_golden\":true,\"falsifier\":true");
}

#[test]
fn short_lifecycle_is_bit_identical_twice() {
    let run = || {
        let mut slot = EngineSlot::default();
        let env = init_short(&mut slot, 1903, MODE_FIXED, 0);
        assert!(env.starts_with("{\"ok\":{"), "{env}");
        let mut transcript = String::new();
        for _ in 0..40 {
            transcript.push_str(&slot.step(false, 0.0, 0.0));
        }
        transcript.push_str(&slot.digest());
        transcript
    };
    let a = run();
    let b = run();
    assert_eq!(a, b, "byte-identical envelope transcripts");
    jlog("determinism", &format!("\"transcript_bytes\":{}", a.len()));
}
