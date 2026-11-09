use super::*;
use crate::perception::PERCEPTION_HZ;
use crate::simloop::{CATAPULT_1904_V1, ControlInput, PilotMode, dec17_scenario};

fn simulation(mode: PilotMode) -> SimLoop {
    let mut spec = dec17_scenario(1903, mode);
    spec.catapult = Some(CATAPULT_1904_V1);
    SimLoop::init(spec).unwrap()
}

#[test]
fn specific_force_limits_distinguish_gravity_support_and_load_factor() {
    let (g, mass) = (motion::GRAVITY_MPS2, 317.0);
    for alpha in [0.0_f64, 0.3, -0.4] {
        // Assembly includes level-flight gravity at incidence alpha. Pure
        // weight must yield zero accelerometer specific force at every angle.
        let gravity_force = [-mass * g * alpha.sin(), mass * g * alpha.cos()];
        let free = motion::nongravity_acceleration(gravity_force, alpha, mass);
        assert!(free[0].abs() < 2e-14 && free[1].abs() < 2e-14);
    }
    for load_factor in [0.0, 1.0, 2.0] {
        let total_force = [0.0, mass * g * (1.0 - load_factor)];
        let specific = motion::nongravity_acceleration(total_force, 0.0, mass)[1];
        assert!((specific + load_factor * g).abs() < 2e-14);
    }
    assert_eq!(rail_specific_force(0.0, 0.0), -g);
    let theta = 0.3;
    assert!(
        (rail_specific_force(theta, 2.0) - (2.0 * theta.sin() - g * theta.cos())).abs() < 2e-14
    );
}

#[test]
fn real_rail_impulse_reaches_perception_at_six_ticks_without_a_new_solve() {
    let mut sim = simulation(PilotMode::FixedControls);
    for cue in &mut sim.perception.cues {
        cue.remnant_sigma = 0.0;
    }
    assert!(sim.heave_cue_sample().is_none());
    assert!(sim.perceived_cues().unwrap().is_none());
    let ground_speed = |sim: &SimLoop| {
        motion::ground_velocity([sim.u, sim.w], sim.theta, [-sim.spec.headwind_mps, 0.0])[0]
    };
    let mut first = 0.0;
    for tick in 0..=6 {
        let start_v = ground_speed(&sim);
        let theta = sim.theta;
        let out = sim.step(None).unwrap();
        let sample = sim.heave_cue_sample().unwrap();
        assert_eq!(sample.tick, tick);
        if tick == 0 {
            assert_eq!(out.phase, Phase::OnRail);
            let end_v = ground_speed(&sim);
            let expected = (theta.sin() * (end_v - start_v) * PERCEPTION_HZ
                - motion::GRAVITY_MPS2 * theta.cos())
                - sample.reference_specific_force_mps2;
            assert!((sample.relative_specific_force_mps2 - expected).abs() < 2e-11);
            first = sample.relative_specific_force_mps2;
            assert!(
                first.abs() > 1e-4,
                "tow acceleration must produce a real cue at pitched trim"
            );
        }
        let perceived = sim.perceived_cues().unwrap().unwrap();
        assert_eq!(perceived.tick, tick);
        if tick < 6 {
            assert_eq!(perceived.values[2], 0.0);
        } else {
            let expected = first * (1.0 - (-1.0 / (PERCEPTION_HZ * 0.08)).exp());
            assert!((perceived.values[2] - expected).abs() < 2e-13);
        }
    }
}

#[test]
fn airborne_sample_uses_the_native_start_load_with_its_own_incidence() {
    let mut sim = simulation(PilotMode::Human);
    sim.phase = Phase::Airborne;
    sim.h_m = 20.0;
    sim.theta = 0.2;
    sim.q = 0.08;
    sim.u = 13.0;
    sim.w = 0.9;
    let mut ou = sim.ou.clone();
    ou.advance_to(1).unwrap();
    let gust: f64 = ou
        .amplitudes()
        .iter()
        .enumerate()
        .map(|(i, a)| a * 0.12 / (1.0 + i as f64))
        .sum();
    let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
    let alpha = (det::atan2(sim.w, sim.u) + gust / speed).clamp(-0.55, 0.55);
    let dc = sim
        .mech_spec
        .step(sim.mech, 0.0, 0.0, 1.0 / PERCEPTION_HZ)
        .unwrap()
        .0
        .delta_rad;
    let initial = sim
        .sampled_loads_current_wind(speed, alpha, dc)
        .unwrap()
        .forces;
    // Independent force balance removes the force assembly's gravity exactly
    // once; current attitude must not be substituted for its incidence.
    let expected = initial.force_n[2] / sim.design.gross_mass_kg
        - motion::GRAVITY_MPS2 * alpha.cos()
        - sim.heave_reference();
    let out = sim
        .step(Some(ControlInput {
            lever_force_n: 0.0,
            warp_cmd_rad: 0.0,
        }))
        .unwrap();
    assert_eq!(out.phase, Phase::Airborne);
    let sample = sim.heave_cue_sample().unwrap();
    assert!((sample.relative_specific_force_mps2 - expected).abs() < 2e-10);
}

#[test]
fn checkpoint_preserves_force_acquisition_and_perceived_noise_without_new_state() {
    let mut sim = simulation(PilotMode::FixedControls);
    for _ in 0..9 {
        sim.step(None).unwrap();
    }
    let checkpoint = sim.save_checkpoint().unwrap();
    let raw = sim.heave_cue_sample();
    let cues = sim.perceived_cues().unwrap();
    assert_eq!(
        checkpoint,
        sim.save_checkpoint().unwrap(),
        "observation is read-only"
    );
    let mut restored = SimLoop::restore_checkpoint(sim.spec.clone(), &checkpoint).unwrap();
    assert_eq!(raw, restored.heave_cue_sample());
    assert_eq!(cues, restored.perceived_cues().unwrap());
    for _ in 0..12 {
        assert_eq!(sim.step(None).unwrap(), restored.step(None).unwrap());
        assert_eq!(sim.heave_cue_sample(), restored.heave_cue_sample());
        assert_eq!(
            sim.perceived_cues().unwrap(),
            restored.perceived_cues().unwrap()
        );
        assert_eq!(
            sim.save_checkpoint().unwrap(),
            restored.save_checkpoint().unwrap()
        );
    }
}

#[test]
fn unsupported_acquisition_and_terminal_load_failure_publish_no_false_sample() {
    let mut sim = simulation(PilotMode::FixedControls);
    sim.perception.cues[2].delay_ticks = 0;
    let before = sim.save_checkpoint().unwrap();
    assert_eq!(
        sim.step(None).unwrap_err().code,
        "perception-deferred-input-invalid"
    );
    assert_eq!(before, sim.save_checkpoint().unwrap());
    sim.perception.cues[2].delay_ticks = 6;
    sim.phase = Phase::Airborne;
    sim.u = 0.0;
    sim.w = 0.0;
    let out = sim.step(None).unwrap();
    assert_eq!(out.phase, Phase::Ended(TerminalEvent::EnvelopeExceeded));
    assert!(sim.heave_cue_sample().is_none());
}
