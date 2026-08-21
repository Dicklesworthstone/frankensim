//! E5.1 wasm-surface battery (bead wf-root-guzez.6.2): the envelope
//! API over the REAL engine, exercised natively (the wasm exports are
//! 1:1 wrappers over `EngineSlot`). Every documented refusal code at
//! this surface is EXECUTED; the short lifecycle is run twice
//! bit-identically; envelopes are well-formed for ok AND refusal.
//! Repro: cd crates/fs-flyer-wasm && cargo test --test engine_battery

use fs_flyer_wasm::engine::{EngineSlot, MODE_FIXED, MODE_HISTORICAL, MODE_HUMAN};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-wasm-e51\",\"case\":\"{case}\",{payload}}}");
}

const DEC17: (f64, f64, f64, u64) = (1.294, 11.0, 18.3, 40);

fn init_short(slot: &mut EngineSlot, seed: u64, mode: u32, member: u32) -> String {
    let (rho, wind, rail, ticks) = DEC17;
    slot.init(seed, rho, wind, mode, member, rail, ticks, false, false)
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
    // control-input-missing: Human mode without input, then with a
    // non-finite lever force.
    let mut human = EngineSlot::default();
    let env = init_short(&mut human, 7, MODE_HUMAN, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let e5 = human.step(false, 0.0, 0.0);
    assert!(e5.contains("\"code\":\"control-input-missing\""), "{e5}");
    let e6 = human.step(true, f64::NAN, 0.0);
    assert!(e6.contains("\"code\":\"control-input-missing\""), "{e6}");
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
fn human_mode_steps_with_input_and_historical_member_binds() {
    let mut human = EngineSlot::default();
    let env = init_short(&mut human, 7, MODE_HUMAN, 0);
    assert!(env.starts_with("{\"ok\":{"), "{env}");
    let s = human.step(true, 25.0, 0.01);
    assert!(s.starts_with("{\"ok\":{"), "{s}");
    assert!(s.contains("\"phase\":\"on-rail\""), "{s}");
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
