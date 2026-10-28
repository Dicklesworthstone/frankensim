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
fn ordinary_rail_initialization_is_stationary_in_the_world_frame() {
    let sim = SimLoop::init(dec17_scenario(1903, PilotMode::FixedControls)).unwrap();
    let (s, c) = sim.theta.sin_cos();
    assert!(
        sim.theta.abs() > 0.001,
        "pitched fixture must expose the frame error"
    );
    assert!((c * sim.u + s * sim.w - sim.spec.headwind_mps).abs() < 1e-12);
    assert!((s * sim.u - c * sim.w).abs() < 1e-12);
    assert!((sim.w - s * sim.spec.headwind_mps).abs() < 1e-12);
    assert_eq!(sim.phase, Phase::OnRail);
}

#[test]
fn releasing_the_rail_reaction_adds_no_horizontal_or_vertical_impulse() {
    for mode in 0..3 {
        let spec = dec17_scenario(1903, PilotMode::FixedControls);
        let mut sim = match mode {
            0 => SimLoop::init(spec).unwrap(),
            1 => SimLoop::init_with_spatial_wind(spec, wind_spec(0.0)).unwrap(),
            _ => SimLoop::init_with_distributed_spatial_wind(spec, wind_spec(0.0)).unwrap(),
        };
        let speed = 3.7;
        if let Some(wind) = &mut sim.spatial_wind {
            wind.ground_mps = [speed, 0.0];
            sim.sync_spatial_air_velocity().unwrap();
        } else {
            [sim.u, sim.w] =
                motion::body_air_velocity([speed, 0.0], sim.theta, [-sim.spec.headwind_mps, 0.0]);
        }
        let start = sim.current_flight_motion();
        let alpha = 0.12;
        let weight = sim.design.gross_mass_kg * motion::GRAVITY_MPS2;
        let upward = 1.1 * weight;
        let (s, c) = (det::sin(sim.theta), det::cos(sim.theta));
        let mut loads = sim.tick0.trim.buildup.clone();
        // A purely vertical prescribed non-gravity force releases the rail;
        // it supplies no horizontal impulse during this constrained tick.
        loads.force_n = [
            s * upward - weight * det::sin(alpha),
            0.0,
            -c * upward + weight * det::cos(alpha),
        ];
        loads.moment_y_nm = 0.0;
        loads.torque_imbalance_nm = 0.0;
        let dt = 1.0 / 120.0;
        sim.advance_rail_motion(&loads, alpha, dt);
        let next = sim.current_flight_motion();
        assert_eq!(sim.phase, Phase::Airborne);
        assert!((next.longitudinal.ground_mps[0] - speed).abs() < 1e-12);
        assert!(next.longitudinal.ground_mps[1].abs() < 1e-12);
        assert!((next.longitudinal.position_m[0] - dt * speed).abs() < 1e-12);
        assert_eq!(
            next.longitudinal.position_m[1],
            start.longitudinal.position_m[1]
        );
        assert_eq!(next.longitudinal.theta_rad, start.longitudinal.theta_rad);
        assert_eq!(next.longitudinal.q_rad_s, start.longitudinal.q_rad_s);
        assert_eq!(next.rotor_rad_s, start.rotor_rad_s);
    }
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
        let start = sim.current_flight_motion();
        let acceleration = sim.flight_acceleration(start, &loads, alpha);
        let stage = motion::midpoint_stage(start, acceleration, 1.0 / 120.0);
        sim.publish_flight_motion(motion::finish_midpoint(
            start,
            stage,
            acceleration,
            1.0 / 120.0,
        ));
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

fn airborne_fixture(mode: usize) -> SimLoop {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let mut sim = match mode {
        0 => SimLoop::init(scenario).unwrap(),
        1 => SimLoop::init_with_spatial_wind(scenario, wind_spec(0.3)).unwrap(),
        _ => SimLoop::init_with_distributed_spatial_wind(scenario, wind_spec(0.3)).unwrap(),
    };
    sim.phase = Phase::Airborne;
    sim.h_m = 20.0;
    sim.q = 0.15;
    sim.omega += 2.0;
    let wind = sim
        .spatial_wind
        .as_ref()
        .map_or([-sim.spec.headwind_mps, 0.0], |w| {
            let sample = w.sample([sim.x_m, 0.0, sim.h_m]).unwrap();
            [sample.u[0], sample.u[2]]
        });
    let ground = [sim.tick0.trim.v_mps + wind[0], wind[1]];
    if let Some(w) = &mut sim.spatial_wind {
        w.ground_mps = ground;
    }
    [sim.u, sim.w] = motion::body_air_velocity(ground, sim.theta, wind);
    sim
}

fn motion_values(state: motion::FlightMotion) -> [f64; 7] {
    [
        state.longitudinal.position_m[0],
        state.longitudinal.position_m[1],
        state.longitudinal.ground_mps[0],
        state.longitudinal.ground_mps[1],
        state.longitudinal.theta_rad,
        state.longitudinal.q_rad_s,
        state.rotor_rad_s,
    ]
}

#[test]
fn actual_airborne_impulses_use_midpoint_aerodynamics_and_rotor_speed() {
    for mode in 0..3 {
        let mut sim = airborne_fixture(mode);
        let gust = if mode == 0 { 0.17 } else { 0.0 };
        let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
        let alpha = (det::atan2(sim.w, sim.u) + gust / speed).clamp(-0.55, 0.55);
        let dc = sim.mech.delta_rad;
        let first = sim.force_buildup_current_wind(speed, alpha, dc).unwrap();
        let start = sim.current_flight_motion();
        let dt = 1.0 / 60.0;
        let stage =
            motion::midpoint_stage(start, sim.flight_acceleration(start, &first, alpha), dt);
        let before_wind = sim.spatial_wind.as_ref().map(|w| w.field.checkpoint());
        let middle = sim
            .advance_airborne_motion(&first, alpha, dc, gust, dt)
            .unwrap()
            .unwrap();
        let next = sim.current_flight_motion();
        assert_ne!(middle.forces.torque_imbalance_nm, first.torque_imbalance_nm);
        assert_ne!(middle.forces.force_n, first.force_n);
        assert_ne!(stage.rotor_rad_s, start.rotor_rad_s);
        let mass = sim.design.gross_mass_kg;
        let weight = mass * motion::GRAVITY_MPS2;
        let body = [
            middle.forces.force_n[0] + weight * middle.alpha_rad.sin(),
            middle.forces.force_n[2] - weight * middle.alpha_rad.cos(),
        ];
        let (s, c) = stage.longitudinal.theta_rad.sin_cos();
        let expected_force = [
            c * body[0] + s * body[1],
            s * body[0] - c * body[1] - weight,
        ];
        for axis in 0..2 {
            let actual =
                mass * (next.longitudinal.ground_mps[axis] - start.longitudinal.ground_mps[axis]);
            assert!((actual - dt * expected_force[axis]).abs() < 2e-10);
        }
        assert!(
            (super::super::IYY_KG_M2 * (next.longitudinal.q_rad_s - start.longitudinal.q_rad_s)
                - dt * middle.forces.moment_y_nm)
                .abs()
                < 2e-12
        );
        assert!(
            (1.6 * (next.rotor_rad_s - start.rotor_rad_s) - dt * middle.forces.torque_imbalance_nm)
                .abs()
                < 2e-12
        );
        assert_eq!(sim.warm_slip, Some(middle.forces.coupled.w_slip));
        assert_eq!(
            sim.spatial_wind.as_ref().map(|w| w.field.checkpoint()),
            before_wind,
            "midpoint load sampling must not consume an OU innovation"
        );
    }
}

#[test]
fn refused_midpoint_keeps_motion_and_warm_start_unchanged() {
    let mut sim = airborne_fixture(0);
    sim.q = 0.0;
    let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
    let alpha = det::atan2(sim.w, sim.u);
    let dc = sim.mech.delta_rad;
    let mut first = sim.force_buildup_current_wind(speed, alpha, dc).unwrap();
    let start = sim.current_flight_motion();
    let warm = sim.warm_slip;
    let dt = 1.0 / 120.0;
    // A finite prescribed start load predicts almost zero airspeed. The
    // ORIGINAL aero owner then refuses that midpoint (its floor is 1 m/s).
    let wind = [-sim.spec.headwind_mps, 0.0];
    let acceleration = [
        (wind[0] - start.longitudinal.ground_mps[0]) / (0.5 * dt),
        (wind[1] - start.longitudinal.ground_mps[1]) / (0.5 * dt) + motion::GRAVITY_MPS2,
    ];
    let body = motion::body_air_velocity(acceleration, sim.theta, [0.0; 2]);
    let mass = sim.design.gross_mass_kg;
    let weight = mass * motion::GRAVITY_MPS2;
    first.force_n = [
        mass * body[0] - weight * det::sin(alpha),
        0.0,
        mass * body[1] + weight * det::cos(alpha),
    ];
    first.moment_y_nm = 0.0;
    first.torque_imbalance_nm = 0.0;
    assert!(
        sim.advance_airborne_motion(&first, alpha, dc, 0.0, dt)
            .is_err()
    );
    assert_eq!(
        motion_values(sim.current_flight_motion()),
        motion_values(start)
    );
    assert_eq!(sim.warm_slip, warm);
    assert_eq!(sim.tick, 0);
    assert_eq!(sim.phase, Phase::Airborne);
}

#[test]
fn actual_ground_crossing_uses_terminal_branch_and_one_endpoint_wind_advance() {
    for spatial in [false, true] {
        let mut sim = if spatial {
            simulation(1903, 0.0)
        } else {
            SimLoop::init(dec17_scenario(1903, PilotMode::FixedControls)).unwrap()
        };
        sim.phase = Phase::Airborne;
        sim.h_m = 0.001;
        sim.q = 0.0;
        let wind = if spatial {
            let sample = sim.sample_spatial_wind([sim.x_m, 0.0, sim.h_m]).unwrap();
            [sample.u[0], sample.u[2]]
        } else {
            [-sim.spec.headwind_mps, 0.0]
        };
        let ground = [sim.tick0.trim.v_mps + wind[0], -0.4 + wind[1]];
        if let Some(w) = &mut sim.spatial_wind {
            w.ground_mps = ground;
        }
        [sim.u, sim.w] = motion::body_air_velocity(ground, sim.theta, wind);
        let mut expected = sim
            .spatial_wind
            .as_ref()
            .map(|w| EvolvingTurbulenceField::resume(w.field.checkpoint()).unwrap());
        if let Some(field) = &mut expected {
            field.advance_to(1).unwrap();
        }
        let out = sim.step(None).unwrap();
        assert_eq!(
            out.phase,
            Phase::Ended(TerminalEvent::GroundContact),
            "{out:?}"
        );
        assert_eq!(out.h_m, 0.0);
        assert_eq!(out.tick, 1);
        assert_eq!(
            sim.spatial_wind.as_ref().map(|w| w.field.checkpoint()),
            expected.map(|w| w.checkpoint())
        );
    }
}
