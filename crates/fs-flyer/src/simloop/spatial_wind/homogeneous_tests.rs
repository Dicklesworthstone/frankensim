//! Actual flight loads, domain admission and replay for homogeneous Mann wind.

use super::*;
use crate::simloop::{ControlInput, Phase, PilotMode, RAIL_HEIGHT_M, dec17_scenario};
use fs_atmo::mann::MANN_TARGET_V1;

fn wind(minimum_height_m: f64) -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: RAIL_HEIGHT_M,
            reference_speed_mps: 11.0,
        },
        spectrum: SpatialWindSpectrum::MannHomogeneousAbove {
            parameters: MannParams {
                alpha_eps23: 0.02,
                ..MANN_TARGET_V1
            },
            min_kl: 0.03,
            max_kl: 30.0,
            minimum_height_m,
        },
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 32],
    }
}

fn airborne(sim: &mut SimLoop, height: f64, climb: f64) {
    sim.phase = Phase::Airborne;
    sim.h_m = height;
    let flow = sim.sample_spatial_wind([sim.x_m, 0.0, height]).unwrap();
    sim.spatial_wind.as_mut().unwrap().ground_mps = [13.0 + flow.u[0], climb];
    sim.sync_spatial_air_velocity().unwrap();
}

#[test]
fn full_tensor_velocity_gradient_and_covariance_reach_the_flight_frame() {
    let spec = wind(0.5);
    let sim = SimLoop::init_with_distributed_spatial_wind(
        dec17_scenario(1903, PilotMode::FixedControls),
        spec.clone(),
    )
    .unwrap();
    let SpatialWindSpectrum::MannHomogeneousAbove {
        parameters,
        min_kl,
        max_kl,
        ..
    } = spec.spectrum
    else {
        unreachable!()
    };
    let model = MannSpectralModel::build(
        1903,
        parameters,
        MannQuadrature {
            n_modes: 32,
            min_kl,
            max_kl,
        },
        [11.0, 0.0, 0.0],
    )
    .unwrap();
    let mut field =
        HomogeneousMannField::at_anchor(model, 1903, spec.correlation_times_s.clone()).unwrap();
    field.advance_to(0).unwrap();
    let before = sim.save_checkpoint().unwrap();
    assert!(sim.spatial_wind_stress_fit().is_none());
    for [x, y, h] in [[0.0, 0.0, 0.5], [1.2, -2.3, 3.0], [-0.5, 4.0, 5.0]] {
        let mut expected = field.sample([-x, y, h]).unwrap();
        expected.u[0] += spec.mean.speed(h);
        expected.grad[0][2] += spec.mean.dspeed_dh(h);
        let signs = [-1.0, 1.0, 1.0];
        for i in 0..3 {
            expected.u[i] *= signs[i];
            for j in 0..3 {
                expected.grad[i][j] *= signs[i] * signs[j];
            }
        }
        assert_eq!(sim.sample_spatial_wind([x, y, h]).unwrap(), expected);
        assert!((expected.grad[0][0] + expected.grad[1][1] + expected.grad[2][2]).abs() < 1e-12);
    }
    let a = [0.3, -0.2, 2.0];
    let b = [0.8, 0.7, 2.6];
    let mut expected = field
        .covariance([a[0] - b[0], b[1] - a[1], b[2] - a[2]], 0.2)
        .unwrap();
    let signs = [-1.0, 1.0, 1.0];
    for i in 0..3 {
        for j in 0..3 {
            expected[i][j] *= signs[i] * signs[j];
        }
    }
    assert_eq!(
        sim.homogeneous_wind_covariance(a, b, 0.2).unwrap(),
        expected
    );
    assert_eq!(
        sim.sample_spatial_wind([0.0, 0.0, 0.49]).unwrap_err().code,
        APPLICABILITY_REFUSAL
    );
    assert_eq!(
        sim.homogeneous_wind_covariance(a, [0.0, 0.0, 0.1], 0.0)
            .unwrap_err()
            .code,
        APPLICABILITY_REFUSAL
    );
    assert_eq!(sim.save_checkpoint().unwrap(), before);
}

#[test]
fn homogeneous_wind_changes_actual_rotor_and_distributed_motion() {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let mut quiet = wind(0.5);
    quiet.spectrum = SpatialWindSpectrum::Analytic {
        sigma_mps: 0.0,
        length_m: 1.77,
    };
    let mut calm = SimLoop::init_with_distributed_spatial_wind(scenario.clone(), quiet).unwrap();
    let mut turbulent = SimLoop::init_with_distributed_spatial_wind(scenario, wind(0.5)).unwrap();
    for sim in [&mut calm, &mut turbulent] {
        airborne(sim, 20.0, 0.0);
    }
    let a = calm.step(None).unwrap();
    let b = turbulent.step(None).unwrap();
    assert_eq!(a.phase, Phase::Airborne);
    assert_eq!(b.phase, Phase::Airborne);
    assert_ne!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    assert_ne!(
        calm.spatial_ground_velocity_mps(),
        turbulent.spatial_ground_velocity_mps()
    );
    assert!((a.p_rad_s - b.p_rad_s).abs() + (a.r_rad_s - b.r_rad_s).abs() > 1e-10);
}

#[test]
fn v8_restores_all_sampling_material_combinations_and_four_quadratures() {
    for distributed in [false, true] {
        for material in 0..3 {
            let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
            scenario.max_ticks = 24;
            let gas = GasState::try_new_moist_air(288.15, 98000.0, 0.4).unwrap();
            let profile = MaterialAirProfile::standard1976(240.0).unwrap();
            let mut sim = match (distributed, material) {
                (false, 0) => SimLoop::init_with_spatial_wind(scenario, wind(0.5)),
                (true, 0) => SimLoop::init_with_distributed_spatial_wind(scenario, wind(0.5)),
                (false, 1) => SimLoop::init_with_material_spatial_wind(scenario, wind(0.5), gas),
                (true, 1) => {
                    SimLoop::init_with_material_distributed_spatial_wind(scenario, wind(0.5), gas)
                }
                (false, _) => SimLoop::init_with_profile_spatial_wind(scenario, wind(0.5), profile),
                (true, _) => SimLoop::init_with_profile_distributed_spatial_wind(
                    scenario,
                    wind(0.5),
                    profile,
                ),
            }
            .unwrap();
            assert_eq!(sim.tick0.ou_amplitudes.len(), 4 * 32);
            for _ in 0..3 {
                sim.step(None).unwrap();
            }
            let checkpoint = sim.save_checkpoint().unwrap();
            assert_eq!(&checkpoint[..4], &8u32.to_le_bytes());
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
fn resealed_v8_cannot_substitute_target_cutoffs_or_applicability_height() {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let sim = SimLoop::init_with_spatial_wind(scenario.clone(), wind(0.5)).unwrap();
    let checkpoint = sim.save_checkpoint().unwrap();
    let begin = checkpoint
        .windows(HOMOGENEOUS_MARKER.len())
        .position(|s| s == HOMOGENEOUS_MARKER)
        .unwrap()
        + HOMOGENEOUS_MARKER.len();
    for field in 0..6 {
        let mut changed = checkpoint.clone();
        let offset = begin + 8 * field;
        let value = f64::from_le_bytes(changed[offset..offset + 8].try_into().unwrap());
        changed[offset..offset + 8].copy_from_slice(&(value * 1.01).to_le_bytes());
        let body = changed.len() - 32;
        let digest = hash_domain("org.frankensim.wf.checkpoint-state.v8", &changed[..body]);
        changed[body..].copy_from_slice(digest.as_bytes());
        let error = SimLoop::restore_checkpoint(scenario.clone(), &changed)
            .err()
            .expect("model mismatch");
        assert_eq!(error.code, "checkpoint-malformed");
        assert!(error.message.contains("different scenario"));
    }
}

#[test]
fn strip_and_midpoint_domain_refusals_preserve_all_flight_control_and_random_state() {
    // CG is admitted, but the lower wing is below the authored height. Both
    // pilot modes mutate perception/control state before the real strip query.
    for mode in [PilotMode::Human, PilotMode::Historical(1)] {
        let mut sim = SimLoop::init_with_distributed_spatial_wind(
            dec17_scenario(1903, mode),
            wind(RAIL_HEIGHT_M),
        )
        .unwrap();
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
    // A CG-only path can leave the domain at a predicted midpoint or at the
    // accepted endpoint. Neither case may consume a tick or an OU innovation.
    for offset in [1e-6, 0.006] {
        let mut sim = SimLoop::init_with_spatial_wind(
            dec17_scenario(1903, PilotMode::FixedControls),
            wind(0.5),
        )
        .unwrap();
        airborne(&mut sim, 0.5 + offset, -1.0);
        let before = sim.save_checkpoint().unwrap();
        assert_eq!(sim.step(None).unwrap_err().code, APPLICABILITY_REFUSAL);
        assert_eq!(sim.save_checkpoint().unwrap(), before);
        airborne(&mut sim, 20.0, 0.0);
        let mut twin =
            SimLoop::restore_checkpoint(sim.spec.clone(), &sim.save_checkpoint().unwrap()).unwrap();
        assert_eq!(sim.step(None).unwrap(), twin.step(None).unwrap());
        assert_eq!(
            sim.save_checkpoint().unwrap(),
            twin.save_checkpoint().unwrap()
        );
    }
}

#[test]
fn homogeneous_recording_scrubs_without_external_spectral_state() {
    use crate::replayenv::{ReplayEnvelope, record_replay_with_distributed_spatial_wind};
    let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
    scenario.max_ticks = 12;
    let mut direct =
        SimLoop::init_with_distributed_spatial_wind(scenario.clone(), wind(0.5)).unwrap();
    let envelope = record_replay_with_distributed_spatial_wind(scenario, wind(0.5), 4).unwrap();
    let restored = ReplayEnvelope::from_bytes(&envelope.to_bytes()).unwrap();
    for tick in 1..=12 {
        let expected = direct.step(None).unwrap();
        if [1, 5, 9, 12].contains(&tick) {
            assert_eq!(restored.scrub_to_tick(tick).unwrap(), expected);
        }
    }
    assert_eq!(restored.final_digest, direct.digest_hex());
    for floor in [0.0, -1.0, f64::NAN, RAIL_HEIGHT_M + 0.1] {
        assert!(
            SimLoop::init_with_spatial_wind(
                dec17_scenario(1903, PilotMode::FixedControls),
                wind(floor)
            )
            .is_err()
        );
    }
}
