//! E5.1 native-core battery (bead wf-root-guzez.6.2): the full
//! lifecycle (init → rail → airborne → terminal) EXECUTED twice with
//! bit-identical digests; RunIntentId minted after (and sensitive to
//! intent, insensitive backward — tick-0 digest never changes with
//! intent); phase transitions with receipts; human-mode input law;
//! scenario caps at cap AND cap+1; explicit catapult-assisted flight
//! fixtures; the original unaided launch shortfall; deterministic golden.
//! Repro: cargo test -p fs-flyer --test simloop_battery

use fs_flyer::simloop::{
    ControlInput, MAX_HEADWIND_MPS, MAX_RAIL_M, MAX_TICKS, MIN_HEADWIND_MPS, Phase, PilotMode,
    SNAPSHOT_LEN, SNAPSHOT_PHASE_INDEX, ScenarioSpec, SimLoop, TerminalEvent, dec17_scenario,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e51\",\"case\":\"{case}\",{payload}}}");
}

// Test-only assisted launch, explicitly distinct from historical Dec-17.
// Correct rail kinematics no longer lift the unchanged 18.3 m reference.
// Use the existing registered tow law to exercise airborne consumers without
// retuning pitch, controls, aerodynamic parameters or pilot gains.
fn assisted_launch_scenario(mode: PilotMode) -> ScenarioSpec {
    let mut spec = dec17_scenario(1903, mode);
    spec.catapult = Some(fs_flyer::simloop::CATAPULT_1904_V1);
    spec
}

fn run_assisted_lifecycle(mode: PilotMode) -> (SimLoop, u64, Phase, f64, f64) {
    let mut sim = SimLoop::init(assisted_launch_scenario(mode)).unwrap();
    let mut last_phase = Phase::OnRail;
    let mut x_end = 0.0;
    let mut max_forward_m = 0.0f64;
    let mut last_out: Option<fs_flyer::simloop::SimStateOut> = None;
    loop {
        let input = match mode {
            PilotMode::Human => Some(ControlInput {
                lever_force_n: 0.0,
                warp_cmd_rad: 0.0,
            }),
            _ => None,
        };
        match sim.step(input) {
            Ok(out) => {
                // Per-second trajectory receipts (JSONL log doctrine).
                if out.tick % 120 == 0 || matches!(out.phase, Phase::Ended(_)) {
                    jlog(
                        "traj",
                        &format!(
                            "\"mode\":\"{mode:?}\",\"tick\":{},\"phase\":\"{:?}\",\"x\":{:.2},\"h\":{:.2},\"u\":{:.2},\"w\":{:.2},\"q\":{:.3},\"theta\":{:.3},\"dc\":{:.3},\"omega\":{:.1}",
                            out.tick,
                            out.phase,
                            out.x_m,
                            out.h_m,
                            out.u_mps,
                            out.w_mps,
                            out.q_rad_s,
                            out.theta_rad,
                            out.dc_rad,
                            out.omega_prop_rad_s
                        ),
                    );
                }
                last_phase = out.phase;
                x_end = out.x_m;
                max_forward_m = max_forward_m.max(out.x_m);
                last_out = Some(out);
                if let Phase::Ended(_) = out.phase {
                    break;
                }
            }
            Err(e) => panic!("lifecycle refusal after {last_out:?}: {e:?}"),
        }
    }
    let tick = last_out.expect("lifecycle has a terminal output").tick;
    assert_eq!(sim.step(None).unwrap_err().code, "run-ended");
    (sim, tick, last_phase, x_end, max_forward_m)
}

#[test]
fn fixed_control_assisted_lifecycle_preserves_its_numerical_refusal() {
    let (sim, tick, phase, x_end, max_forward_m) = run_assisted_lifecycle(PilotMode::FixedControls);
    // Same-tick release retains its free vertical/pitch impulse. The aircraft
    // climbs and turns back before the native outer coupling reaches
    // its correction budget. This is a numerical refusal, not proof of a physical
    // flight-envelope boundary. Keep its exact code/tick and real travel gate.
    assert_eq!(phase, Phase::Ended(TerminalEvent::EnvelopeExceeded));
    assert_eq!(tick, 655);
    let r = sim
        .envelope_refusal()
        .expect("numerical terminal keeps its receipt");
    assert_eq!(r.code, "PropAirframeCouplingDidNotConverge");
    assert!(r.message.starts_with("cap 12 exhausted:"));
    jlog("envelope-receipt", &format!("\"code\":\"{}\"", r.code));
    assert!(
        max_forward_m > 5.0,
        "must travel before turning back: {max_forward_m} m"
    );
    jlog(
        "fixed-lifecycle",
        &format!(
            "\"terminal\":\"{phase:?}\",\"end_tick\":{tick},\"distance_m\":{x_end},\"maximum_forward_m\":{max_forward_m},\"digest\":\"{}\"",
            &sim.digest_hex()[..16]
        ),
    );
}

#[test]
fn lifecycle_is_bit_identical_twice() {
    // The DONE-WHEN determinism clause, native form.
    let (a, _, _, ax, _) = run_assisted_lifecycle(PilotMode::FixedControls);
    let (b, _, _, bx, _) = run_assisted_lifecycle(PilotMode::FixedControls);
    assert_eq!(a.digest_hex(), b.digest_hex(), "bit-identical lifecycles");
    assert_eq!(ax.to_bits(), bx.to_bits());
    // And the historical pilot too (remnant streams are tick-addressed).
    let (c, _, _, _, _) = run_assisted_lifecycle(PilotMode::Historical(2));
    let (d, _, _, _, _) = run_assisted_lifecycle(PilotMode::Historical(2));
    assert_eq!(c.digest_hex(), d.digest_hex());
    assert_ne!(a.digest_hex(), c.digest_hex(), "modes must differ");
    jlog(
        "determinism",
        &format!("\"digest\":\"{}\"", &a.digest_hex()[..16]),
    );
}

#[test]
fn historical_pilot_assisted_fixture_reports_its_terminal() {
    let (sim, _, phase, x_end, _) = run_assisted_lifecycle(PilotMode::Historical(2));
    assert!(matches!(
        phase,
        Phase::Ended(TerminalEvent::GroundContact) | Phase::Ended(TerminalEvent::EnvelopeExceeded)
    ));
    assert!(x_end > 5.0);
    jlog(
        "historical",
        &format!(
            "\"terminal\":\"{phase:?}\",\"distance_m\":{x_end},\"digest\":\"{}\"",
            &sim.digest_hex()[..16]
        ),
    );
}

#[test]
fn run_intent_id_is_minted_after_and_downstream_of_tick0() {
    let a = SimLoop::init(dec17_scenario(1903, PilotMode::FixedControls)).unwrap();
    let b = SimLoop::init(dec17_scenario(1903, PilotMode::Historical(0))).unwrap();
    // Same seed/scenario physics: the tick-0 digest is IDENTICAL (intent
    // cannot reach backward)...
    assert_eq!(a.tick0().digest, b.tick0().digest);
    // ...while RunIntentId differs (intent binds forward only).
    assert_ne!(a.run_intent_id, b.run_intent_id);
    // Headwind is intent-bearing scenario state: different headwind =
    // different RunIntentId.
    let mut s = dec17_scenario(1903, PilotMode::FixedControls);
    s.headwind_mps = 9.0;
    let c = SimLoop::init(s).unwrap();
    assert_ne!(a.run_intent_id, c.run_intent_id);
    jlog(
        "run-intent",
        &format!("\"tick0\":\"{}\"", &a.tick0().digest[..16]),
    );
}

#[test]
fn human_input_refusals_do_not_consume_physical_time() {
    let spec = dec17_scenario(7, PilotMode::Human);
    let mut retried = SimLoop::init(spec.clone()).unwrap();
    let checkpoint_before = retried.save_checkpoint().unwrap();
    let digest_before = retried.digest_hex();

    // Missing and invalid input are distinct typed refusals. Neither may
    // advance the tick-addressed atmosphere or any other checkpoint state.
    assert_eq!(
        retried.step(None).unwrap_err().code,
        "control-input-missing"
    );
    assert_eq!(retried.save_checkpoint().unwrap(), checkpoint_before);
    assert_eq!(retried.digest_hex(), digest_before);

    for invalid in [
        ControlInput {
            lever_force_n: f64::NAN,
            warp_cmd_rad: 0.0,
        },
        ControlInput {
            lever_force_n: 0.0,
            warp_cmd_rad: f64::INFINITY,
        },
    ] {
        assert_eq!(
            retried.step(Some(invalid)).unwrap_err().code,
            "control-input-invalid"
        );
        assert_eq!(retried.save_checkpoint().unwrap(), checkpoint_before);
        assert_eq!(retried.digest_hex(), digest_before);
    }

    // A corrected retry is the same physical tick as a first-attempt valid
    // control, not the next atmosphere sample after a hidden skipped tick.
    let valid = ControlInput {
        lever_force_n: 10.0,
        warp_cmd_rad: 0.02,
    };
    let retried_out = retried.step(Some(valid)).unwrap();
    let mut fresh = SimLoop::init(spec).unwrap();
    let fresh_out = fresh.step(Some(valid)).unwrap();
    assert_eq!(retried_out, fresh_out);
    assert_eq!(retried.digest_hex(), fresh.digest_hex());
    assert_eq!(
        retried.save_checkpoint().unwrap(),
        fresh.save_checkpoint().unwrap()
    );
    jlog(
        "human-input",
        "\"required_every_tick\":true,\"refusals_are_atomic\":true",
    );
}

#[test]
fn human_warp_drives_real_lateral_state_into_payload_v2() {
    let mut spec = assisted_launch_scenario(PilotMode::Human);
    spec.max_ticks = 1_200;
    let mut sim = SimLoop::init(spec).unwrap();
    let mut last = None;
    for _ in 0..1_200 {
        let out = sim
            .step(Some(ControlInput {
                lever_force_n: 0.0,
                warp_cmd_rad: 0.08,
            }))
            .unwrap();
        last = Some(out);
        if out.phi_rad.abs() > 1.0e-5 && out.psi_rad.abs() > 1.0e-8 {
            break;
        }
        if matches!(out.phase, Phase::Ended(_)) {
            break;
        }
    }
    let out = last.expect("at least one step");
    assert!(out.phi_rad > 1.0e-5, "positive warp banks right: {out:?}");
    assert!(out.p_rad_s > 0.0, "bank has a real roll rate: {out:?}");
    assert!(
        out.psi_rad.abs() > 1.0e-8,
        "linked yaw changes heading: {out:?}"
    );
    assert!(
        out.r_rad_s.abs() > 0.0,
        "heading has a real yaw rate: {out:?}"
    );
    let payload = sim.snapshot_payload(&out);
    assert_eq!(payload.len(), SNAPSHOT_LEN);
    assert_eq!(
        payload[SNAPSHOT_PHASE_INDEX], 1.0,
        "airborne phase stays in the v1 slot"
    );
    assert_eq!(payload[12].to_bits(), out.phi_rad.to_bits());
    assert_eq!(payload[13].to_bits(), out.psi_rad.to_bits());
    jlog(
        "payload-v2-lateral",
        &format!(
            "\"tick\":{},\"phi_rad\":{},\"psi_rad\":{}",
            out.tick, out.phi_rad, out.psi_rad
        ),
    );
}

#[test]
fn scenario_caps_at_cap_and_cap_plus_one() {
    let mk = |f: &dyn Fn(&mut ScenarioSpec)| -> Result<SimLoop, fs_flyer::Refusal> {
        let mut s = dec17_scenario(1, PilotMode::FixedControls);
        f(&mut s);
        SimLoop::init(s)
    };
    assert!(mk(&|s| s.headwind_mps = MAX_HEADWIND_MPS).is_ok());
    assert!(matches!(
        mk(&|s| s.headwind_mps = MAX_HEADWIND_MPS.next_up()),
        Err(e) if e.code == "scenario-invalid"
    ));
    // The FLOOR too (E5.4): at the floor admits, one below refuses —
    // a sub-floor start would refuse inside the coupled solve at tick
    // 1 instead of at admission (measured on the calm Huffman twin).
    assert!(mk(&|s| s.headwind_mps = MIN_HEADWIND_MPS).is_ok());
    assert!(matches!(
        mk(&|s| s.headwind_mps = MIN_HEADWIND_MPS.next_down()),
        Err(e) if e.code == "scenario-invalid"
    ));
    assert!(mk(&|s| s.rail_length_m = MAX_RAIL_M).is_ok());
    assert!(matches!(
        mk(&|s| s.rail_length_m = MAX_RAIL_M.next_up()),
        Err(e) if e.code == "scenario-invalid"
    ));
    assert!(mk(&|s| s.max_ticks = MAX_TICKS).is_ok());
    assert!(matches!(
        mk(&|s| s.max_ticks = MAX_TICKS + 1),
        Err(e) if e.code == "scenario-invalid"
    ));
    assert!(matches!(
        mk(&|s| s.max_ticks = 0),
        Err(e) if e.code == "scenario-invalid"
    ));
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn rail_phase_produces_receipted_transition() {
    let mut sim = SimLoop::init(assisted_launch_scenario(PilotMode::FixedControls)).unwrap();
    let mut saw_rail = false;
    let mut transition_tick = None;
    for _ in 0..2400 {
        let out = sim.step(None).unwrap();
        match out.phase {
            Phase::OnRail => saw_rail = true,
            Phase::Airborne => {
                transition_tick = Some(out.tick);
                break;
            }
            Phase::Ended(e) => panic!("ended before liftoff: {e:?}"),
        }
    }
    assert!(saw_rail, "must start on the rail");
    let t = transition_tick.expect("the explicit catapult fixture must lift off");
    assert!(t > 1, "liftoff is not instantaneous");
    jlog("rail", &format!("\"liftoff_tick\":{t}"));
}

#[test]
fn member3_assisted_fixture_has_bounded_undulating_flight() {
    // Registered gains are unchanged. This explicit catapult fixture yields
    // 11.41 s airborne, five undulations and 44.90 m before the unchanged
    // outer coupling solver refuses midpoint residual growth at tick 1443.
    // It is a bounded plant regression, not a historical calibration claim.
    // The original unaided reference is checked separately below.
    let run = || {
        let mut sim = SimLoop::init(assisted_launch_scenario(PilotMode::Historical(3))).unwrap();
        let mut liftoff = None;
        let mut undulation_flips = 0u32;
        let mut last_sign = 0i8;
        let mut end = None;
        loop {
            let out = sim.step(None).unwrap();
            if matches!(out.phase, Phase::Airborne) {
                if liftoff.is_none() {
                    liftoff = Some(out.tick);
                }
                let s = if out.q_rad_s > 1e-3 {
                    1i8
                } else if out.q_rad_s < -1e-3 {
                    -1i8
                } else {
                    0
                };
                if s != 0 && last_sign != 0 && s != last_sign {
                    undulation_flips += 1;
                }
                if s != 0 {
                    last_sign = s;
                }
            }
            if let Phase::Ended(e) = out.phase {
                end = Some((e, out.tick, out.x_m));
                break;
            }
        }
        let digest = sim.digest_hex();
        let refusal = sim.envelope_refusal().map(|r| r.code);
        (
            liftoff.unwrap(),
            undulation_flips / 2,
            end.unwrap(),
            digest,
            refusal,
        )
    };
    let (liftoff, undulations, (terminal, end_tick, x_end), digest, refusal) = run();
    assert_eq!(
        terminal,
        TerminalEvent::EnvelopeExceeded,
        "the explicit assisted fixture has a measured numerical solver refusal"
    );
    assert_eq!(refusal, Some("nonlinear-did-not-converge"));
    assert_eq!(
        end_tick, 1_444,
        "registered controls and the corrected plant must reproduce the terminal tick"
    );
    let airborne_s = (end_tick - liftoff) as f64 / 120.0;
    assert!(
        airborne_s >= 8.0,
        "controlled flight duration: {airborne_s} s"
    );
    assert!(undulations >= 3, "undulating flight: {undulations}");
    assert!(
        (25.0..90.0).contains(&x_end),
        "ground distance in the assisted-flight band: {x_end} m"
    );
    // Deterministic registration evidence.
    let (_, _, _, digest2, _) = run();
    assert_eq!(digest, digest2, "bit-identical twice");
    jlog(
        "member3-flight",
        &format!(
            "\"airborne_s\":{airborne_s},\"undulations\":{undulations},\"x_end_m\":{x_end},\"terminal\":\"{terminal:?}\",\"end_tick\":{end_tick},\"digest\":\"{}\"",
            &digest[..16]
        ),
    );
}

#[test]
fn original_dec17_reference_reports_the_corrected_rail_shortfall() {
    for mode in [
        PilotMode::FixedControls,
        PilotMode::Historical(2),
        PilotMode::Historical(3),
    ] {
        let spec = dec17_scenario(1903, mode);
        let mut sim = SimLoop::init(spec.clone()).unwrap();
        let end = loop {
            let out = sim.step(None).unwrap();
            // Independently reconstruct ground vertical velocity: the pitched
            // aircraft remains constrained to the rail, including its endpoint.
            let (s, c) = out.theta_rad.sin_cos();
            assert!((s * out.u_mps - c * out.w_mps).abs() < 2e-12);
            assert_ne!(
                out.phase,
                Phase::Airborne,
                "do not fabricate a launch: {out:?}"
            );
            if let Phase::Ended(_) = out.phase {
                break out;
            }
        };
        assert_eq!(end.phase, Phase::Ended(TerminalEvent::RailEndWithoutLift));
        assert_eq!(end.tick, 975);
        assert!((spec.rail_length_m..spec.rail_length_m + 0.04).contains(&end.x_m));
        assert!(sim.envelope_refusal().is_none());
        if mode == PilotMode::FixedControls {
            // The corrected actuator-disk wash changes the settled trim and
            // digest while preserving this same tick-975 rail shortfall.
            assert_eq!(
                sim.digest_hex(),
                "949fca76eb2879406d5616bec83955a83b7cf0e63f485595125095d55d56d36d"
            );
        }
        jlog(
            "corrected-rail-shortfall",
            &format!(
                "\"mode\":\"{mode:?}\",\"end_tick\":{},\"x_end_m\":{}",
                end.tick, end.x_m
            ),
        );
    }
}

#[test]
fn catapult_mechanics_and_the_huffman_low_speed_block() {
    // E5.4 (bead guzez.6.7). Part 1 — catapult mechanics INSIDE the
    // admitted numerical domain: explicit 60 m test rails with the same KDH wind
    // let both runs lift; the registered 1904 catapult launches earlier
    // (executed comparative discriminator), the tow binds into
    // RunIntentId but never the tick-0 digest, and both catapult axes
    // cap at cap AND cap+1.
    use fs_flyer::simloop::{
        CATAPULT_1904_V1, CatapultSpec, MAX_CATAPULT_FORCE_N, MAX_CATAPULT_PULL_M, huffman_scenario,
    };
    let liftoff_x = |catapult: Option<CatapultSpec>| -> f64 {
        let mut spec = dec17_scenario(1903, PilotMode::FixedControls);
        spec.rail_length_m = 60.0; // test-only extension; original shortfall is tested above
        spec.catapult = catapult;
        let mut sim = SimLoop::init(spec).unwrap();
        loop {
            let out = sim.step(None).unwrap();
            if matches!(out.phase, Phase::Airborne) {
                return out.x_m;
            }
            assert!(!matches!(out.phase, Phase::Ended(_)), "must lift: {out:?}");
        }
    };
    let x_unaided = liftoff_x(None);
    let x_cat = liftoff_x(Some(CATAPULT_1904_V1));
    assert!(
        x_cat < 0.6 * x_unaided,
        "the catapult must shorten the run substantially: {x_cat} vs {x_unaided} m"
    );
    jlog(
        "catapult",
        &format!("\"liftoff_x_m\":{x_cat},\"unaided_liftoff_x_m\":{x_unaided}"),
    );
    // Caps at cap AND cap+1 on both catapult axes.
    let mk = |c: CatapultSpec| {
        let mut s = dec17_scenario(1, PilotMode::FixedControls);
        s.catapult = Some(c);
        SimLoop::init(s)
    };
    assert!(
        mk(CatapultSpec {
            pull_force_n: MAX_CATAPULT_FORCE_N,
            pull_length_m: MAX_CATAPULT_PULL_M
        })
        .is_ok()
    );
    for bad in [
        CatapultSpec {
            pull_force_n: MAX_CATAPULT_FORCE_N.next_up(),
            pull_length_m: 10.0,
        },
        CatapultSpec {
            pull_force_n: 1000.0,
            pull_length_m: MAX_CATAPULT_PULL_M.next_up(),
        },
        CatapultSpec {
            pull_force_n: 0.0,
            pull_length_m: 10.0,
        },
    ] {
        assert!(
            matches!(mk(bad), Err(e) if e.code == "scenario-invalid"),
            "{bad:?}"
        );
    }
    // Intent binding: catapult present/absent changes RunIntentId,
    // never the tick-0 digest.
    let mut with_spec = dec17_scenario(1, PilotMode::FixedControls);
    with_spec.catapult = Some(CATAPULT_1904_V1);
    let with = SimLoop::init(with_spec).unwrap();
    let without = SimLoop::init(dec17_scenario(1, PilotMode::FixedControls)).unwrap();
    assert_ne!(with.run_intent_id, without.run_intent_id);
    assert_eq!(with.tick0().digest, without.tick0().digest);
    assert_eq!(CATAPULT_1904_V1.pull_force_n, 2073.0);
    assert_eq!(CATAPULT_1904_V1.pull_length_m, 14.6);

    // Part 2 — an executed low-speed numerical limitation (bead guzez.5.7.1).
    // The registered near-calm scenario still refuses with the unchanged
    // catapult and solver budgets. This receipt is not evidence of a physical
    // launch boundary. A future solver improvement may admit this same case;
    // it must retain the physical inputs and independently demonstrate liftoff.
    let mut sim = SimLoop::init(huffman_scenario(1904, PilotMode::FixedControls)).unwrap();
    let end = loop {
        let out = sim.step(None).unwrap();
        assert!(
            !matches!(out.phase, Phase::Airborne),
            "low-speed solve improved: preserve the inputs and register the liftoff receipt"
        );
        if let Phase::Ended(e) = out.phase {
            break e;
        }
    };
    assert_eq!(end, TerminalEvent::EnvelopeExceeded);
    let receipt = sim.envelope_refusal().expect("receipt kept");
    jlog(
        "huffman-blocked",
        &format!("\"code\":\"{}\",\"tracked\":\"guzez.5.7.1\"", receipt.code),
    );

    // Part 3 — unchanged low-speed inputs on both sides of the current
    // numerical limit. Correct disk wash and the native strip-cycle safeguard
    // retain the 4 m/s liftoff gate. At 3 m/s the wing now solves, but the
    // propeller station budget is exhausted later; keep that actual receipt.
    let outcome_at = |headwind: f64| -> (bool, u64, Option<String>) {
        let mut spec = huffman_scenario(1904, PilotMode::FixedControls);
        spec.headwind_mps = headwind;
        let mut sim = SimLoop::init(spec.clone()).unwrap();
        loop {
            let out = sim.step(None).unwrap();
            if matches!(out.phase, Phase::Airborne) {
                assert!(out.x_m < spec.rail_length_m, "must lift within the rail");
                return (true, out.tick, None);
            }
            if let Phase::Ended(e) = out.phase {
                assert_eq!(e, TerminalEvent::EnvelopeExceeded);
                return (
                    false,
                    out.tick,
                    Some(
                        sim.envelope_refusal()
                            .expect("numerical receipt")
                            .code
                            .to_owned(),
                    ),
                );
            }
        }
    };
    let (lift3, tick3, refusal3) = outcome_at(3.0);
    assert!(!lift3, "the unresolved 3 m/s solve must keep its refusal");
    assert_eq!(tick3, 8);
    assert_eq!(refusal3.as_deref(), Some("station-did-not-converge"));
    let (lift4, _, _) = outcome_at(4.0);
    assert!(lift4, "4 m/s must lift with the declared catapult");
    jlog(
        "huffman-floor-bracket",
        "\"refuses_at\":3.0,\"lifts_at\":4.0",
    );
}

#[test]
fn golden_digest() {
    let (sim, _, _, _, _) = run_assisted_lifecycle(PilotMode::FixedControls);
    let digest = sim.digest_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    // Golden bump 2026-08-21 (bead guzez.7.2.1): ground-contact
    // friction moved off platform tanh onto det::tanh (and the
    // Prandtl factors off platform acos), unifying the native
    // and wasm lanes; the lifecycle trajectory shifted by ulps.
    // Golden bump 2026-08-24 (bead frankensim-4pa2k): snapshot v2
    // appends the simulated roll and heading states to every chained
    // payload; the v1 prefix and longitudinal trajectory are unchanged.
    // Golden bump 2026-10-28: converged passive hinge friction changes
    // the settled initial state and the ensuing deterministic trajectory.
    // Golden bump 2026-10-28: airborne translation now advances world
    // momentum with actual vertical gravity and midpoint held-load motion;
    // explicit body-frame transport no longer creates kinetic energy.
    // Golden bump 2026-10-28: coupled translation/pitch/rotor midpoint loads
    // and a physically stationary rail state replace tick-start flight loads
    // and the spurious launch velocity. This flight golden explicitly uses
    // the existing catapult fixture; the unmodified reference's rail refusal
    // has its own exact digest in the preceding regression.
    // Golden bump 2026-10-28: already-unloaded rail states retain the full
    // same-tick free-flight impulse; positive midpoint stages decide contact
    // from their accepted endpoint. At that revision, station nonconvergence
    // at 673 was retained as a numerical receipt, not a physical limit.
    // Golden bump 2026-10-28: stable propeller station momentum arithmetic
    // and residual-certified local finishing solve the former tick-673 case.
    // The same evaluation budget still refuses an unresolved station at 676;
    // the equation, controls and substantive flight gates are unchanged.
    // Golden bump 2026-10-29: actuator-disk vi is no longer halved again
    // before entering the spatial wash law. The unchanged assisted fixture
    // now reaches the outer correction cap at tick 655; no physical gate,
    // pilot gain, tolerance or iteration cap was relaxed for this digest.
    assert_eq!(
        digest, "73d1fbb87f3df4b830cf92fafc1f4603d6839742dac485d4b7166a43e43ac5f7",
        "lifecycle golden moved — determinism regression or an \
         intentional physics change requiring the golden-bump protocol"
    );
}
