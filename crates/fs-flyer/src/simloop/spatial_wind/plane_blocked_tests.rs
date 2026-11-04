//! Plane-blocked Mann wind in the original flight, load and replay owners.

use super::*;
use crate::simloop::{ControlInput, Phase, PilotMode, RAIL_HEIGHT_M, dec17_scenario};
use fs_atmo::mann::MANN_TARGET_V1;

fn wind() -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: RAIL_HEIGHT_M,
            reference_speed_mps: 11.0,
        },
        spectrum: SpatialWindSpectrum::MannPlaneBlocked {
            parameters: MannParams {
                alpha_eps23: 0.02,
                ..MANN_TARGET_V1
            },
            min_kl: 0.03,
            max_kl: 30.0,
        },
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 32],
    }
}

fn native(spec: &SpatialWindSpec) -> HomogeneousMannField {
    let SpatialWindSpectrum::MannPlaneBlocked {
        parameters,
        min_kl,
        max_kl,
    } = spec.spectrum
    else {
        unreachable!()
    };
    let model = MannSpectralModel::build(
        1903,
        parameters,
        MannQuadrature {
            n_modes: spec.correlation_times_s.len(),
            min_kl,
            max_kl,
        },
        [spec.advection_mps, 0.0, 0.0],
    )
    .unwrap();
    let mut field =
        HomogeneousMannField::at_anchor(model, 1903, spec.correlation_times_s.clone()).unwrap();
    field.advance_to(0).unwrap();
    field
}

#[test]
fn ground_blocking_retains_native_velocity_gradient_and_complete_covariance() {
    let spec = wind();
    let field = native(&spec);
    let sim = SimLoop::init_with_distributed_spatial_wind(
        dec17_scenario(1903, PilotMode::FixedControls),
        spec.clone(),
    )
    .unwrap();
    let before = sim.save_checkpoint().unwrap();
    let signs = [-1.0, 1.0, 1.0];
    for [x, y, h] in [[0.7, -0.3, 0.0], [1.2, -2.3, 0.2], [-0.5, 4.0, 3.0]] {
        let mut expected = field
            .sample_wall_blocked([-x, y, h], ground_plane())
            .unwrap();
        expected.u[0] += spec.mean.speed(h);
        expected.grad[0][2] += spec.mean.dspeed_dh(h);
        for i in 0..3 {
            expected.u[i] *= signs[i];
            for j in 0..3 {
                expected.grad[i][j] *= signs[i] * signs[j];
            }
        }
        let actual = sim.sample_spatial_wind([x, y, h]).unwrap();
        assert_eq!(actual, expected);
        assert!((actual.grad[0][0] + actual.grad[1][1] + actual.grad[2][2]).abs() < 1e-12);
        if h == 0.0 {
            assert_eq!(actual.u[2], 0.0);
            assert_eq!(actual.grad[2][0], 0.0);
            assert_eq!(actual.grad[2][1], 0.0);
        }
    }
    let a = [0.3, -0.2, 0.0];
    let b = [0.8, 0.7, 0.6];
    let mut expected = field
        .covariance_wall_blocked(
            [-a[0], a[1], a[2]],
            [-b[0], b[1], b[2]],
            0.2,
            ground_plane(),
        )
        .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            expected[i][j] *= signs[i] * signs[j];
        }
    }
    let actual = sim.mann_wind_covariance(a, b, 0.2).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(actual[2], [0.0; 3]);
    assert!(actual[0][2].abs() > 1e-10);
    let elevated = sim
        .mann_wind_covariance([a[0], a[1], 2.0], [b[0], b[1], 2.6], 0.2)
        .unwrap();
    assert_ne!(elevated, actual); // Equal separation, different wall distances.
    assert!(sim.homogeneous_wind_covariance(a, b, 0.2).is_err());
    for height in [-1e-12, f64::INFINITY, f64::NAN] {
        assert_eq!(
            sim.sample_spatial_wind([0.0, 0.0, height])
                .unwrap_err()
                .code,
            APPLICABILITY_REFUSAL
        );
    }
    assert_eq!(sim.save_checkpoint().unwrap(), before);
}

#[test]
fn plane_correction_changes_actual_strip_rotor_and_distributed_motion() {
    let spec = wind();
    let mut homogeneous = spec.clone();
    let SpatialWindSpectrum::MannPlaneBlocked {
        parameters,
        min_kl,
        max_kl,
    } = spec.spectrum
    else {
        unreachable!()
    };
    homogeneous.spectrum = SpatialWindSpectrum::MannHomogeneousAbove {
        parameters,
        min_kl,
        max_kl,
        minimum_height_m: 0.01,
    };
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let mut blocked = SimLoop::init_with_distributed_spatial_wind(scenario.clone(), spec).unwrap();
    let mut unblocked = SimLoop::init_with_distributed_spatial_wind(scenario, homogeneous).unwrap();
    // Only the spatial observation changes; both use the same random history
    // and the same physical aircraft position, attitude and ground velocity.
    assert_eq!(blocked.tick0.ou_amplitudes, unblocked.tick0.ou_amplitudes);
    for sim in [&mut blocked, &mut unblocked] {
        sim.phase = Phase::Airborne;
        sim.h_m = 3.0;
        sim.spatial_wind.as_mut().unwrap().ground_mps = [2.0, 0.0];
        sim.sync_spatial_air_velocity().unwrap();
    }
    let loads = |sim: &SimLoop| {
        sim.force_buildup_current_wind(
            det::hypot(sim.u, sim.w),
            det::atan2(sim.w, sim.u),
            sim.mech.delta_rad,
        )
        .unwrap()
    };
    let a = loads(&blocked);
    let b = loads(&unblocked);
    assert!((a.lift_n - b.lift_n).abs() > 1e-7);
    assert!((a.moment_y_nm - b.moment_y_nm).abs() > 1e-7);
    assert!((a.torque_imbalance_nm - b.torque_imbalance_nm).abs() > 1e-7);
    let a = blocked.step(None).unwrap();
    let b = unblocked.step(None).unwrap();
    assert_eq!(a.phase, Phase::Airborne);
    assert_eq!(b.phase, Phase::Airborne);
    assert_ne!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    assert_ne!(
        blocked.spatial_ground_velocity_mps(),
        unblocked.spatial_ground_velocity_mps()
    );
    assert!((a.p_rad_s - b.p_rad_s).abs() + (a.r_rad_s - b.r_rad_s).abs() > 1e-10);
}

#[test]
fn v9_restores_four_paths_for_every_sampling_and_material_choice() {
    for distributed in [false, true] {
        for material in 0..3 {
            let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
            scenario.max_ticks = 24;
            let gas = GasState::try_new_moist_air(288.15, 98000.0, 0.4).unwrap();
            let profile = MaterialAirProfile::standard1976(240.0).unwrap();
            let mut sim = match (distributed, material) {
                (false, 0) => SimLoop::init_with_spatial_wind(scenario, wind()),
                (true, 0) => SimLoop::init_with_distributed_spatial_wind(scenario, wind()),
                (false, 1) => SimLoop::init_with_material_spatial_wind(scenario, wind(), gas),
                (true, 1) => {
                    SimLoop::init_with_material_distributed_spatial_wind(scenario, wind(), gas)
                }
                (false, _) => SimLoop::init_with_profile_spatial_wind(scenario, wind(), profile),
                (true, _) => {
                    SimLoop::init_with_profile_distributed_spatial_wind(scenario, wind(), profile)
                }
            }
            .unwrap();
            for _ in 0..3 {
                sim.step(None).unwrap();
            }
            let checkpoint = sim.save_checkpoint().unwrap();
            assert_eq!(&checkpoint[..4], &9u32.to_le_bytes());
            let mut resumed = SimLoop::restore_checkpoint(sim.spec.clone(), &checkpoint).unwrap();
            assert_eq!(resumed.save_checkpoint().unwrap(), checkpoint);
            for _ in 0..4 {
                assert_eq!(sim.step(None).unwrap(), resumed.step(None).unwrap());
                assert_eq!(sim.digest_hex(), resumed.digest_hex());
            }
            assert_eq!(
                sim.save_checkpoint().unwrap(),
                resumed.save_checkpoint().unwrap()
            );
        }
    }
}

#[test]
fn below_ground_strip_refusal_restores_flight_controls_and_all_random_paths() {
    for mode in [PilotMode::Human, PilotMode::Historical(1)] {
        let mut sim =
            SimLoop::init_with_distributed_spatial_wind(dec17_scenario(1903, mode), wind())
                .unwrap();
        sim.phase = Phase::Airborne;
        sim.h_m = 0.01;
        sim.theta = 0.0;
        let flow = sim.sample_spatial_wind([0.0, 0.0, sim.h_m]).unwrap();
        sim.spatial_wind.as_mut().unwrap().ground_mps = [13.0 + flow.u[0], flow.u[2]];
        sim.sync_spatial_air_velocity().unwrap();
        let before = sim.save_checkpoint().unwrap();
        let input = Some(ControlInput {
            lever_force_n: 30.0,
            warp_cmd_rad: 0.04,
        });
        for _ in 0..2 {
            assert_eq!(sim.step(input).unwrap_err().code, APPLICABILITY_REFUSAL);
            assert_eq!(sim.save_checkpoint().unwrap(), before);
        }
    }
}

#[test]
fn plane_blocked_recording_scrubs_and_binds_its_spectral_inputs() {
    use crate::replayenv::{ReplayEnvelope, record_replay_with_distributed_spatial_wind};
    let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
    scenario.max_ticks = 12;
    let mut direct = SimLoop::init_with_distributed_spatial_wind(scenario.clone(), wind()).unwrap();
    let checkpoint = direct.save_checkpoint().unwrap();
    let begin = checkpoint
        .windows(PLANE_BLOCKED_MARKER.len())
        .position(|s| s == PLANE_BLOCKED_MARKER)
        .unwrap()
        + PLANE_BLOCKED_MARKER.len();
    for field in 0..5 {
        let mut changed = checkpoint.clone();
        let offset = begin + 8 * field;
        let value = f64::from_le_bytes(changed[offset..offset + 8].try_into().unwrap());
        changed[offset..offset + 8].copy_from_slice(&(value * 1.01).to_le_bytes());
        let body = changed.len() - 32;
        let digest = hash_domain("org.frankensim.wf.checkpoint-state.v9", &changed[..body]);
        changed[body..].copy_from_slice(digest.as_bytes());
        let error = SimLoop::restore_checkpoint(scenario.clone(), &changed)
            .err()
            .unwrap();
        assert_eq!(error.code, "checkpoint-malformed");
        assert!(error.message.contains("different scenario"));
    }
    let envelope = record_replay_with_distributed_spatial_wind(scenario, wind(), 4).unwrap();
    let restored = ReplayEnvelope::from_bytes(&envelope.to_bytes()).unwrap();
    for tick in 1..=12 {
        let expected = direct.step(None).unwrap();
        if [1, 5, 9, 12].contains(&tick) {
            assert_eq!(restored.scrub_to_tick(tick).unwrap(), expected);
        }
    }
    assert_eq!(restored.final_digest, direct.digest_hex());
}
