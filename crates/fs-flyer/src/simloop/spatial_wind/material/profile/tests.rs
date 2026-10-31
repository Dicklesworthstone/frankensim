use super::*;
use crate::replayenv::{
    ReplayEnvelope, record_replay_with_profile_distributed_spatial_wind,
    record_replay_with_profile_spatial_wind,
};
use crate::simloop::{ControlInput, Phase, PilotMode, dec17_scenario, motion};
use fs_atmo::FlatSiteLogLaw;
use fs_math::det;

fn wind(sigma_mps: f64) -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: RAIL_HEIGHT_M,
            reference_speed_mps: 11.0,
        },
        sigma_mps,
        length_m: 20.0,
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 8],
    }
}

fn initialize(distributed: bool, profile: MaterialAirProfile, pilot: PilotMode) -> SimLoop {
    let scenario = dec17_scenario(1903, pilot);
    if distributed {
        SimLoop::init_with_profile_distributed_spatial_wind(scenario, wind(0.3), profile).unwrap()
    } else {
        SimLoop::init_with_profile_spatial_wind(scenario, wind(0.3), profile).unwrap()
    }
}

fn airborne(sim: &mut SimLoop, height: f64, climb: f64) {
    sim.phase = Phase::Airborne;
    sim.h_m = height;
    let flow = sim.sample_spatial_wind([sim.x_m, 0.0, height]).unwrap();
    sim.theta = sim.tick0.trim.alpha_rad + det::atan2(climb - flow.u[2], 13.0);
    sim.spatial_wind.as_mut().unwrap().ground_mps = [13.0 + flow.u[0], climb];
    sim.sync_spatial_air_velocity().unwrap();
}

fn motion_values(state: motion::FlightMotion) -> [f64; 7] {
    let p = state.longitudinal;
    [
        p.position_m[0],
        p.position_m[1],
        p.ground_mps[0],
        p.ground_mps[1],
        p.theta_rad,
        p.q_rad_s,
        state.rotor_rad_s,
    ]
}

#[test]
fn profile_density_drives_both_native_stages_warp_and_lateral_motion() {
    for distributed in [false, true] {
        let profile = MaterialAirProfile::standard1976(0.0).unwrap();
        let mut sim = initialize(distributed, profile, PilotMode::Human);
        airborne(&mut sim, 2_000.0, 2.0);
        sim.q = 0.08;
        sim.lateral_state.p_rad_s = 0.03;
        sim.lateral_state.r_rad_s = -0.02;
        let before = sim.save_checkpoint().unwrap();
        let dt = 1.0 / 120.0;
        let speed = (sim.u * sim.u + sim.w * sim.w).sqrt().max(0.1);
        let alpha = det::atan2(sim.w, sim.u).clamp(-0.55, 0.55);
        let dc = sim
            .mech_spec
            .step(sim.mech, 0.0, 0.0, dt)
            .unwrap()
            .0
            .delta_rad;
        let first = sim.sampled_loads_current_wind(speed, alpha, dc).unwrap();
        assert_eq!(
            first.rho_kg_m3,
            profile.gas_at_height(sim.h_m).unwrap().density
        );
        assert!(first.rho_kg_m3 < 0.85 * sim.scenario().rho_kg_m3);
        let start = sim.current_flight_motion();
        let stage = motion::midpoint_stage(
            start,
            sim.flight_acceleration(start, &first.forces, alpha),
            dt,
        );
        let [x, h] = stage.longitudinal.position_m;
        let air = sim.sample_material_air([x, 0.0, h]).unwrap();
        assert_eq!(air.gas, profile.gas_at_height(h).unwrap());
        assert_ne!(air.gas.density, first.rho_kg_m3);
        let [u, w] = motion::body_air_velocity(
            stage.longitudinal.ground_mps,
            stage.longitudinal.theta_rad,
            [air.wind.u[0], air.wind.u[2]],
        );
        let speed = det::sqrt(u * u + w * w).max(0.1);
        let alpha = det::atan2(w, u).clamp(-0.55, 0.55);
        // Call the original aircraft owner with independently sourced density.
        // This oracle detects freezing the rail/start density at the midpoint.
        let loads = if distributed {
            sim.design.force_buildup_with_incident(
                speed,
                alpha,
                dc.clamp(-0.55, 0.55),
                stage.rotor_rad_s.clamp(10.5, 119.5),
                stage.longitudinal.q_rad_s.clamp(-2.9, 2.9),
                air.gas.density,
                Some(first.forces.coupled.w_slip),
                &|point| sim.incident_at_motion(stage, point),
            )
        } else {
            sim.design.force_buildup_warm(
                speed,
                alpha,
                dc.clamp(-0.55, 0.55),
                stage.rotor_rad_s.clamp(10.5, 119.5),
                stage.longitudinal.q_rad_s.clamp(-2.9, 2.9),
                air.gas.density,
                Some(first.forces.coupled.w_slip),
            )
        }
        .unwrap();
        let expected = motion::finish_midpoint(
            start,
            stage,
            sim.flight_acceleration(stage, &loads, alpha),
            dt,
        );
        let loaded = sim
            .warp_model
            .evaluate(0.08, 0.5 * air.gas.density * speed * speed, alpha)
            .unwrap();
        let (projection, norm) = loaded
            .strips
            .iter()
            .zip(&sim.warp_model.basis)
            .fold((0.0, 0.0), |(p, n), (strip, basis)| {
                (p + strip.loaded_rad * basis, n + basis * basis)
            });
        let mut lateral = sim.lateral_state;
        if distributed {
            let moments = sim
                .design
                .distributed_lateral_loads(&loads, dc.clamp(-0.55, 0.55), air.gas.density)
                .unwrap()
                .net();
            sim.lateral_model
                .step_with_moments(
                    &mut lateral,
                    projection / norm,
                    speed,
                    air.gas.density,
                    dt,
                    moments,
                )
                .unwrap();
        } else {
            sim.lateral_model
                .step(&mut lateral, projection / norm, speed, air.gas.density, dt)
                .unwrap();
        }
        assert_eq!(
            sim.save_checkpoint().unwrap(),
            before,
            "stage queries are read-only"
        );
        let out = sim
            .step(Some(ControlInput {
                lever_force_n: 0.0,
                warp_cmd_rad: 0.08,
            }))
            .unwrap();
        assert_eq!(out.phase, Phase::Airborne);
        assert_eq!(
            motion_values(sim.current_flight_motion()),
            motion_values(expected)
        );
        assert_eq!(sim.lateral_state, lateral);
        assert_eq!(sim.warm_slip, Some(loads.coupled.w_slip));
        assert_eq!(
            sim.spatial_wind.as_ref().unwrap().field.tick(),
            out.tick as i64
        );
        let final_air = sim.sample_material_air([out.x_m, 0.0, out.h_m]).unwrap();
        assert_eq!(final_air.gas, profile.gas_at_height(out.h_m).unwrap());
        assert!(final_air.gas.pressure < profile.gas_at_height(2_000.0).unwrap().pressure);
    }
}

#[test]
fn changing_altitude_changes_actual_flight_while_uniform_profile_preserves_native_motion() {
    for distributed in [false, true] {
        let profile = MaterialAirProfile::standard1976(0.0).unwrap();
        let gas = profile.gas_at_height(RAIL_HEIGHT_M).unwrap();
        let uniform = MaterialAirProfile::uniform(gas).unwrap();
        let mut varying = initialize(distributed, profile, PilotMode::FixedControls);
        let mut fixed = initialize(distributed, uniform, PilotMode::FixedControls);
        let mut native = if distributed {
            SimLoop::init_with_material_distributed_spatial_wind(
                fixed.scenario().clone(),
                wind(0.3),
                gas,
            )
        } else {
            SimLoop::init_with_material_spatial_wind(fixed.scenario().clone(), wind(0.3), gas)
        }
        .unwrap();
        assert_eq!(varying.tick0().trim, fixed.tick0().trim);
        for sim in [&mut varying, &mut fixed, &mut native] {
            airborne(sim, 2_000.0, 2.0);
        }
        for _ in 0..4 {
            let a = varying.step(None).unwrap();
            let b = fixed.step(None).unwrap();
            assert_eq!(b, native.step(None).unwrap());
            assert_eq!(a.phase, Phase::Airborne);
            assert_eq!(b.phase, Phase::Airborne);
            assert_ne!(a.h_m, b.h_m);
            assert_ne!(a.omega_prop_rad_s, b.omega_prop_rad_s);
        }
    }
}

#[test]
fn evolving_profile_checkpoint_and_replay_preserve_the_model_not_one_density() {
    for distributed in [false, true] {
        let profile = MaterialAirProfile::standard1976(300.0).unwrap();
        let mut sim = initialize(distributed, profile, PilotMode::FixedControls);
        airborne(&mut sim, 600.0, 1.0);
        for _ in 0..3 {
            sim.step(None).unwrap();
        }
        let checkpoint = sim.save_checkpoint().unwrap();
        assert_eq!(&checkpoint[..4], &6_u32.to_le_bytes());
        let mut resumed = SimLoop::restore_checkpoint(sim.scenario().clone(), &checkpoint).unwrap();
        assert_eq!(resumed.material_air_profile(), Some(&profile));
        assert_eq!(resumed.save_checkpoint().unwrap(), checkpoint);
        for _ in 0..4 {
            assert_eq!(sim.step(None).unwrap(), resumed.step(None).unwrap());
            assert_eq!(sim.digest_hex(), resumed.digest_hex());
        }
        assert_ne!(
            resumed.sample_material_air([0.0, 0.0, 600.0]).unwrap().gas,
            resumed.sample_material_air([0.0, 0.0, 700.0]).unwrap().gas,
        );
        let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
        scenario.max_ticks = 6;
        let envelope = if distributed {
            record_replay_with_profile_distributed_spatial_wind(scenario, wind(0.3), profile, 3)
        } else {
            record_replay_with_profile_spatial_wind(scenario, wind(0.3), profile, 3)
        }
        .unwrap();
        let decoded = ReplayEnvelope::from_bytes(&envelope.to_bytes()).unwrap();
        let mut direct =
            SimLoop::restore_checkpoint(decoded.spec.clone(), &decoded.checkpoints[0].1).unwrap();
        assert_eq!(
            decoded.scrub_to_tick(0).unwrap(),
            direct.replay_initial_state().unwrap()
        );
        for tick in 1..=decoded.events.terminal_tick {
            assert_eq!(
                decoded.scrub_to_tick(tick).unwrap(),
                direct.step(None).unwrap()
            );
        }
        assert_eq!(decoded.final_digest, direct.digest_hex());
    }
}

#[test]
fn profile_domain_refusal_does_not_publish_predicted_motion_or_warm_start() {
    let profile = MaterialAirProfile::standard1976(0.0).unwrap();
    let mut sim = SimLoop::init_with_profile_spatial_wind(
        dec17_scenario(1903, PilotMode::FixedControls),
        wind(0.0),
        profile,
    )
    .unwrap();
    for height in [79_999.99, 79_999.98] {
        // The first case crosses the profile top at the force point; the second
        // has an admitted midpoint but an out-of-domain accepted endpoint.
        airborne(&mut sim, height, 3.0);
        let before = sim.save_checkpoint().unwrap();
        let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
        let alpha = det::atan2(sim.w, sim.u);
        let dc = sim.mech.delta_rad;
        let first = sim.force_buildup_current_wind(speed, alpha, dc).unwrap();
        let error = sim
            .advance_airborne_motion(&first, alpha, dc, 0.0, 1.0 / 120.0)
            .err()
            .unwrap();
        assert_eq!(error.code, "material-atmosphere-invalid");
        assert_eq!(sim.save_checkpoint().unwrap(), before);
    }
}

#[test]
fn profile_checkpoint_rejects_changed_definition_downgrade_and_invalid_restored_height() {
    let profile = MaterialAirProfile::standard1976(300.0).unwrap();
    let sim = initialize(false, profile, PilotMode::FixedControls);
    let checkpoint = sim.save_checkpoint().unwrap();
    let body = &checkpoint[..checkpoint.len() - 32];
    let mut wind_bytes = Vec::new();
    sim.spatial_wind
        .as_ref()
        .unwrap()
        .encode(&sim.run_intent_id, &mut wind_bytes);
    let profile_start = body.len() - wind_bytes.len() - 9;
    let seal = |mut body: Vec<u8>, domain| {
        let digest = hash_domain(domain, &body);
        body.extend_from_slice(digest.as_bytes());
        body
    };
    let restore = |bytes: &[u8]| SimLoop::restore_checkpoint(sim.scenario().clone(), bytes);
    let domain = "org.frankensim.wf.checkpoint-state.v6";
    let mut changed = body.to_vec();
    changed[profile_start + 1..profile_start + 9].copy_from_slice(&301.0_f64.to_le_bytes());
    assert!(restore(&seal(changed, domain)).is_err());
    let mut unknown = body.to_vec();
    unknown[profile_start] = 2;
    assert!(restore(&seal(unknown, domain)).is_err());
    let mut height = body.to_vec();
    height[21..29].copy_from_slice(&80_000.0_f64.to_le_bytes());
    assert!(restore(&seal(height, domain)).is_err());
    let mut downgraded = body.to_vec();
    downgraded.drain(profile_start - 1..profile_start + 9);
    downgraded[..4].copy_from_slice(&3_u32.to_le_bytes());
    assert!(restore(&seal(downgraded, "org.frankensim.wf.checkpoint-state.v3")).is_err());
}

#[test]
fn profile_codec_retains_complete_mixture_and_the_geometric_datum() {
    let moist = fs_material::gas::GasState::try_new_moist_air(293.15, 101_325.0, 0.6).unwrap();
    for profile in [
        MaterialAirProfile::uniform(moist).unwrap(),
        MaterialAirProfile::standard1976(1_500.0).unwrap(),
    ] {
        let mut bytes = Vec::new();
        encode_profile(&profile, &mut bytes);
        let (decoded, used) = decode_profile(&bytes).unwrap();
        assert_eq!(decoded, profile);
        assert_eq!(used, bytes.len());
        for height in [0.0, 600.0, 5_000.0] {
            assert_eq!(decoded.gas_at_height(height), profile.gas_at_height(height));
        }
        assert!(decode_profile(&[]).is_err());
        assert!(decode_profile(&bytes[..bytes.len() - 1]).is_err());
    }
}
