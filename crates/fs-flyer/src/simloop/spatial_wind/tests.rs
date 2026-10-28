use super::*;
use crate::simloop::{ControlInput, Phase, PilotMode, TerminalEvent, dec17_scenario};

fn wind_spec(sigma_mps: f64) -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: 2.4,
            reference_speed_mps: 11.0,
        },
        sigma_mps,
        length_m: 20.0,
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 8],
    }
}

fn simulation(seed: u64, sigma_mps: f64) -> SimLoop {
    SimLoop::init_with_spatial_wind(
        dec17_scenario(seed, PilotMode::FixedControls),
        wind_spec(sigma_mps),
    )
    .unwrap()
}

#[test]
fn spatial_gust_changes_actual_force_and_propeller_torque() {
    let mut calm = simulation(1903, 0.0);
    let mut turbulent = simulation(1903, 0.8);
    assert_eq!(turbulent.spatial_ground_velocity_mps(), Some([0.0; 2]));
    let wind = turbulent.sample_spatial_wind([0.0, 0.0, 2.4]).unwrap();
    let (s, c) = (det::sin(turbulent.theta), det::cos(turbulent.theta));
    assert_eq!(turbulent.u, -c * wind.u[0] - s * wind.u[2]);
    assert_eq!(turbulent.w, -s * wind.u[0] + c * wind.u[2]);
    let a = calm.step(None).unwrap();
    let b = turbulent.step(None).unwrap();
    assert!(!matches!(a.phase, Phase::Ended(_)), "{a:?}");
    assert!(!matches!(b.phase, Phase::Ended(_)), "{b:?}");
    // Changing only reported air-relative velocity would not change these:
    // the gust must reach both the real force solve and rotor torque balance.
    assert_ne!(
        calm.spatial_ground_velocity_mps(),
        turbulent.spatial_ground_velocity_mps()
    );
    assert_ne!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    let at_end = turbulent.sample_spatial_wind([b.x_m, 0.0, b.h_m]).unwrap();
    assert_eq!(b.gust_w_mps.to_bits(), at_end.u[2].to_bits());
    assert_eq!(
        turbulent.spatial_wind.as_ref().unwrap().field.tick(),
        b.tick as i64
    );
}

#[test]
fn zero_strength_has_no_hidden_legacy_gust_or_seed_dependent_load() {
    let mut a = simulation(1903, 0.0);
    let mut b = simulation(7, 0.0);
    assert_ne!(a.run_intent_id, b.run_intent_id);
    for _ in 0..24 {
        let x = a.step(None).unwrap();
        let y = b.step(None).unwrap();
        assert_eq!(
            x, y,
            "zero spatial amplitude must not add the legacy OU gust"
        );
        assert_eq!(x.gust_w_mps, 0.0);
        assert_eq!(
            a.spatial_ground_velocity_mps(),
            b.spatial_ground_velocity_mps()
        );
        if matches!(x.phase, Phase::Ended(_)) {
            panic!("fixture should remain live through 24 ticks: {x:?}");
        }
    }
}

#[test]
fn observers_share_physics_field_without_changing_it() {
    let mut sim = simulation(1903, 0.8);
    let before = sim.save_checkpoint().unwrap();
    let near = sim.sample_spatial_wind([0.0, 0.0, 2.4]).unwrap();
    let far = sim.sample_spatial_wind([17.0, -3.0, 2.4]).unwrap();
    assert_ne!(near.u, far.u);
    for i in (0..200).rev() {
        let p = [i as f64 * 0.2, -2.0, 0.0];
        assert_eq!(sim.sample_spatial_wind(p).unwrap().u[2], 0.0);
    }
    assert_eq!(sim.sample_spatial_wind([0.0, 0.0, 2.4]).unwrap(), near);
    assert_eq!(sim.save_checkpoint().unwrap(), before);
    let mut untouched = SimLoop::restore_checkpoint(sim.spec.clone(), &before).unwrap();
    assert_eq!(sim.step(None).unwrap(), untouched.step(None).unwrap());
    assert_eq!(sim.digest_hex(), untouched.digest_hex());
}

#[test]
fn checkpoint_resumes_rail_and_airborne_wind_and_forces_bit_for_bit() {
    let mut sim = simulation(1903, 0.3);
    let initial = sim.save_checkpoint().unwrap();
    assert_eq!(&initial[..4], &3u32.to_le_bytes());
    let mut resumed = SimLoop::restore_checkpoint(sim.spec.clone(), &initial).unwrap();
    assert_eq!(resumed.save_checkpoint().unwrap(), initial);
    for _ in 0..24 {
        let a = sim.step(None).unwrap();
        let b = resumed.step(None).unwrap();
        assert_eq!(a, b);
        assert_eq!(sim.digest_hex(), resumed.digest_hex());
        assert_eq!(a.phase, Phase::OnRail, "rail fixture: {a:?}");
    }

    // Replay tests need both integration branches, without assuming that
    // this synthetic gust seed and fixed controls yield a successful launch.
    // Construct a live airborne state at the admitted trim airspeed; ground
    // velocity includes the current wind and u,w are synchronized exactly as
    // in production before checkpointing. Subsequent forces/OU steps are real.
    sim.phase = Phase::Airborne;
    sim.h_m = 20.0;
    sim.theta = sim.tick0.trim.alpha_rad;
    let wind = sim.sample_spatial_wind([sim.x_m, 0.0, sim.h_m]).unwrap();
    sim.spatial_wind.as_mut().unwrap().ground_mps = [sim.tick0.trim.v_mps + wind.u[0], wind.u[2]];
    sim.sync_spatial_air_velocity().unwrap();
    let checkpoint = sim.save_checkpoint().unwrap();
    resumed = SimLoop::restore_checkpoint(sim.spec.clone(), &checkpoint).unwrap();
    assert_eq!(resumed.save_checkpoint().unwrap(), checkpoint);
    for _ in 0..16 {
        let a = sim.step(None).unwrap();
        let b = resumed.step(None).unwrap();
        assert_eq!(a, b);
        assert_eq!(sim.digest_hex(), resumed.digest_hex());
        assert_eq!(a.phase, Phase::Airborne, "airborne fixture: {a:?}");
    }
    assert_eq!(
        sim.save_checkpoint().unwrap(),
        resumed.save_checkpoint().unwrap()
    );
}

#[test]
fn gust_does_not_rotate_gravity_or_impulse_ground_velocity() {
    let mut sim = simulation(1903, 0.0);
    sim.phase = Phase::Airborne;
    sim.h_m = 20.0;
    sim.theta = 0.3;
    let weight = sim.design.gross_mass_kg * 9.80665;
    let mut loads = sim.tick0.trim.buildup.clone();
    loads.moment_y_nm = 0.0;
    loads.torque_imbalance_nm = 0.0;
    for alpha in [-0.2, 0.15] {
        sim.spatial_wind.as_mut().unwrap().ground_mps = [0.0; 2];
        // A force build-up containing only its documented gravity term.
        loads.force_n = [-weight * det::sin(alpha), 0.0, weight * det::cos(alpha)];
        sim.advance_spatial_motion(&loads, alpha, 1.0 / 120.0);
        let ground = sim.spatial_ground_velocity_mps().unwrap();
        assert_eq!(ground[0], 0.0, "gravity must never accelerate horizontally");
        assert_eq!(ground[1], -9.80665 / 120.0);
        let before = ground;
        sim.tick += 1;
        sim.finish_spatial_tick().unwrap();
        assert_eq!(sim.spatial_ground_velocity_mps().unwrap(), before);
    }
    // An envelope exit and ground crossing in the same integration step
    // still produce a terminal output, not a below-surface sampling error.
    sim.phase = Phase::Ended(TerminalEvent::EnvelopeExceeded);
    sim.h_m = -0.01;
    sim.tick += 1;
    assert_eq!(sim.finish_spatial_tick().unwrap(), Some(0.0));
    assert_eq!(sim.h_m, 0.0);
    assert_eq!(sim.phase, Phase::Ended(TerminalEvent::EnvelopeExceeded));
}

#[test]
fn refused_input_leaves_wind_and_checkpoint_unchanged() {
    let mut spec = dec17_scenario(1903, PilotMode::Human);
    let mut sim = SimLoop::init_with_spatial_wind(spec.clone(), wind_spec(0.3)).unwrap();
    let checkpoint = sim.save_checkpoint().unwrap();
    assert_eq!(sim.step(None).unwrap_err().code, "control-input-missing");
    let bad = ControlInput {
        lever_force_n: f64::NAN,
        warp_cmd_rad: 0.0,
    };
    assert_eq!(
        sim.step(Some(bad)).unwrap_err().code,
        "control-input-invalid"
    );
    assert_eq!(sim.save_checkpoint().unwrap(), checkpoint);
    spec.headwind_mps = 10.0;
    assert!(SimLoop::restore_checkpoint(spec, &checkpoint).is_err());
    let mut corrupt = checkpoint;
    let last = corrupt.len() - 1;
    corrupt[last] ^= 1;
    assert_eq!(
        SimLoop::restore_checkpoint(sim.spec.clone(), &corrupt)
            .err()
            .unwrap()
            .code,
        "checkpoint-tampered"
    );
    assert_ne!(sim.phase(), Phase::Ended(TerminalEvent::EnvelopeExceeded));
}

#[test]
fn distributed_geometry_loads_share_one_field_and_v4_replays_exactly() {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let mut sim =
        SimLoop::init_with_distributed_spatial_wind(scenario.clone(), wind_spec(0.3)).unwrap();
    let cg = SimLoop::init_with_spatial_wind(scenario.clone(), wind_spec(0.3)).unwrap();
    assert_ne!(sim.run_intent_id, cg.run_intent_id);
    let checkpoint = sim.save_checkpoint().unwrap();
    assert_eq!(&checkpoint[..4], &4u32.to_le_bytes());
    let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
    let alpha = det::atan2(sim.w, sim.u).clamp(-0.55, 0.55);
    let loads = sim
        .force_buildup_current_wind(speed, alpha, sim.mech.delta_rad)
        .unwrap();
    let (panels, strips) = sim.design.layout(sim.mech.delta_rad).unwrap();
    let local = loads
        .coupled
        .strip_velocity_mps
        .as_ref()
        .expect("actual distributed force path");
    assert_eq!(local.len(), strips.len());
    for (s, strip) in strips.iter().enumerate() {
        let panel = panels[strip.panel_indices[0]];
        let point = [
            0.5 * (panel.a[0] + panel.b[0]),
            0.5 * (panel.a[1] + panel.b[1]),
            0.5 * (panel.a[2] + panel.b[2]),
        ];
        let incident = sim.incident_at_body_point(point).unwrap();
        // At q=0 prop wash changes only u; every strip's vertical wind
        // comes from the same current physical field, not a per-query RNG.
        assert_eq!(local[s][2].to_bits(), incident[2].to_bits());
    }
    assert_eq!(
        sim.save_checkpoint().unwrap(),
        checkpoint,
        "sampling loads cannot consume OU time"
    );
    let mut resumed = SimLoop::restore_checkpoint(scenario.clone(), &checkpoint).unwrap();
    assert_eq!(resumed.save_checkpoint().unwrap(), checkpoint);
    for _ in 0..8 {
        let a = sim.step(None).unwrap();
        let b = resumed.step(None).unwrap();
        assert_eq!(a, b);
        assert!(
            !matches!(a.phase, Phase::Ended(_)),
            "short distributed fixture must remain live"
        );
        assert_eq!(sim.digest_hex(), resumed.digest_hex());
    }
    // Rehashing a V4 payload as V3 cannot silently select CG-only loading.
    let mut downgraded = checkpoint;
    downgraded[..4].copy_from_slice(&3u32.to_le_bytes());
    let end = downgraded.len() - 32;
    let digest = hash_domain("org.frankensim.wf.checkpoint-state.v3", &downgraded[..end]);
    downgraded[end..].copy_from_slice(digest.as_bytes());
    assert!(SimLoop::restore_checkpoint(scenario, &downgraded).is_err());
}

#[test]
fn distributed_surface_below_ground_refuses_instead_of_using_cg_flow() {
    let mut sim = SimLoop::init_with_distributed_spatial_wind(
        dec17_scenario(1903, PilotMode::FixedControls),
        wind_spec(0.0),
    )
    .unwrap();
    sim.h_m = 0.1; // CG remains above the plane, some real surfaces do not.
    let out = sim.step(None).unwrap();
    assert_eq!(out.phase, Phase::Ended(TerminalEvent::EnvelopeExceeded));
    assert_eq!(sim.envelope_refusal().unwrap().code, "below-surface-query");
}
