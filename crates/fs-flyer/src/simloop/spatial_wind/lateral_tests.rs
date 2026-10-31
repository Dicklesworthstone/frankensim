use super::*;
use crate::simloop::{ControlInput, Phase, PilotMode, dec17_scenario};

fn simulation(distributed: bool) -> SimLoop {
    let scenario = dec17_scenario(1903, PilotMode::Human);
    let wind = SpatialWindSpec {
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
    };
    if distributed {
        SimLoop::init_with_distributed_spatial_wind(scenario, wind).unwrap()
    } else {
        SimLoop::init_with_spatial_wind(scenario, wind).unwrap()
    }
}

fn airborne(sim: &mut SimLoop) {
    sim.phase = Phase::Airborne;
    sim.h_m = 20.0;
    sim.sync_spatial_air_velocity().unwrap();
}

fn controls() -> Option<ControlInput> {
    Some(ControlInput { lever_force_n: 0.0, warp_cmd_rad: 0.0 })
}

#[test]
fn airborne_step_consumes_the_accepted_gust_torque_once_and_replays() {
    let mut sim = simulation(true);
    airborne(&mut sim);
    let dt = 1.0 / super::super::PERCEPTION_HZ;
    let speed = (sim.u * sim.u + sim.w * sim.w).sqrt().max(0.1);
    let alpha = det::atan2(sim.w, sim.u).clamp(-0.55, 0.55);
    // Independently prepare the physical controls and two native load stages
    // of the next human-input step. Neither query advances the wind field.
    let dc = sim.mech_spec.step(sim.mech, 0.0, 0.0, dt).unwrap().0.delta_rad;
    let first = sim.force_buildup_current_wind(speed, alpha, dc).unwrap();
    let start = sim.current_flight_motion();
    let stage = motion::midpoint_stage(start, sim.flight_acceleration(start, &first, alpha), dt);
    let [x, h] = stage.longitudinal.position_m;
    let wind = sim.sample_spatial_wind([x, 0.0, h]).unwrap();
    let [u, w] = motion::body_air_velocity(stage.longitudinal.ground_mps,
        stage.longitudinal.theta_rad, [wind.u[0], wind.u[2]]);
    let speed = det::sqrt(u * u + w * w).max(0.1);
    let alpha = det::atan2(w, u).clamp(-0.55, 0.55);
    let loads = sim.force_buildup_at_motion(stage, speed, alpha, dc,
        Some(first.coupled.w_slip)).unwrap();
    let moments = sim.design.distributed_lateral_loads(
        &loads, dc.clamp(-0.55, 0.55), sim.spec.rho_kg_m3,
    ).unwrap().net();
    assert!(moments.roll_nm.abs() + moments.yaw_nm.abs() > 1e-8,
        "the fixture must actually resolve asymmetric gust loading");
    let loaded = sim.warp_model.evaluate(0.0,
        0.5 * sim.spec.rho_kg_m3 * speed * speed, alpha).unwrap();
    let (projection, norm) = loaded.strips.iter().zip(&sim.warp_model.basis)
        .fold((0.0, 0.0), |(p, n), (strip, basis)| {
            (p + strip.loaded_rad * basis, n + basis * basis)
        });
    let mut expected = sim.lateral_state;
    sim.lateral_model.step_with_moments(&mut expected, projection / norm,
        speed, sim.spec.rho_kg_m3, dt, moments).unwrap();
    let before = sim.save_checkpoint().unwrap();
    let mut replay = SimLoop::restore_checkpoint(sim.spec.clone(), &before).unwrap();
    let out = sim.step(controls()).unwrap();
    assert_eq!(out.phase, Phase::Airborne);
    assert_eq!(sim.lateral_state, expected);
    assert_eq!(out.p_rad_s, expected.p_rad_s);
    assert_eq!(out.r_rad_s, expected.r_rad_s);
    assert_eq!(out.phi_rad, expected.phi_rad);
    assert_eq!(out.psi_rad, expected.psi_rad);
    assert_eq!(out, replay.step(controls()).unwrap());
    assert_eq!(sim.digest_hex(), replay.digest_hex());
    // Restore after the torques have produced nonzero physical angular state.
    let checkpoint = sim.save_checkpoint().unwrap();
    replay = SimLoop::restore_checkpoint(sim.spec.clone(), &checkpoint).unwrap();
    assert_eq!(replay.save_checkpoint().unwrap(), checkpoint);
    for _ in 0..8 {
        let a = sim.step(controls()).unwrap();
        let b = replay.step(controls()).unwrap();
        assert_eq!(a, b);
        assert_eq!(sim.digest_hex(), replay.digest_hex());
        assert_eq!(a.phase, Phase::Airborne);
    }
}

#[test]
fn rail_constraints_and_cg_only_lateral_behavior_stay_unchanged() {
    let mut rail = simulation(true);
    assert_eq!(rail.phase, Phase::OnRail);
    let out = rail.step(controls()).unwrap();
    assert!(!matches!(out.phase, Phase::Ended(_)));
    assert_eq!([out.phi_rad, out.p_rad_s, out.psi_rad, out.r_rad_s], [0.0; 4]);
    let mut cg = simulation(false);
    airborne(&mut cg);
    let out = cg.step(controls()).unwrap();
    assert_eq!(out.phase, Phase::Airborne);
    assert!(out.p_rad_s.abs() < 1e-12 && out.r_rad_s.abs() < 1e-12,
        "CG-only gusts must not gain distributed moment forcing");
}

#[test]
fn missing_distributed_force_velocities_cannot_silently_disable_gust_torque() {
    let mut sim = simulation(true);
    airborne(&mut sim);
    let speed = (sim.u * sim.u + sim.w * sim.w).sqrt();
    let alpha = det::atan2(sim.w, sim.u).clamp(-0.55, 0.55);
    let dc = sim.mech.delta_rad;
    let mut loads = sim.force_buildup_current_wind(speed, alpha, dc).unwrap();
    loads.coupled.strip_velocity_mps = None;
    let before = sim.save_checkpoint().unwrap();
    let error = sim.step_lateral_current_wind(0.0, speed, sim.spec.rho_kg_m3, dc, &loads, 1.0 / 120.0)
        .unwrap_err();
    assert_eq!(error.code, "distributed-loads-invalid");
    assert_eq!(sim.save_checkpoint().unwrap(), before);
}

#[test]
fn pre_lateral_distributed_checkpoint_identity_is_not_silently_reused() {
    let sim = simulation(true);
    let cg = simulation(false);
    let old = hash_domain("org.frankensim.wf.distributed-wind-intent.v1",
        cg.run_intent_id.as_bytes()).to_hex();
    assert_ne!(sim.run_intent_id, old);
    let mut checkpoint = sim.save_checkpoint().unwrap();
    assert_eq!(&checkpoint[..4], &4u32.to_le_bytes());
    // The unchanged V4 layout ends with the 64-byte run identity and digest.
    // Use a valid old integrity envelope, not merely corrupt bytes.
    let end = checkpoint.len() - 32;
    checkpoint[end - 64..end].copy_from_slice(old.as_bytes());
    let digest = hash_domain("org.frankensim.wf.checkpoint-state.v4", &checkpoint[..end]);
    checkpoint[end..].copy_from_slice(digest.as_bytes());
    let error = SimLoop::restore_checkpoint(sim.spec.clone(), &checkpoint).err().unwrap();
    assert_eq!(error.code, "checkpoint-malformed");
    assert!(error.message.contains("different scenario"));
    // A fresh matching V4 checkpoint still restores byte-for-byte.
    let current = sim.save_checkpoint().unwrap();
    let replay = SimLoop::restore_checkpoint(sim.spec.clone(), &current).unwrap();
    assert_eq!(replay.save_checkpoint().unwrap(), current);
}
