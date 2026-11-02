//! G0/G3/G5: the calibrated field reaches physical loads and exact replay.
//! These are finite-mode consistency tests, not historical wind validation.

use super::*;
use crate::simloop::{Phase, PilotMode, dec17_scenario};
use fs_atmo::mann::MANN_TARGET_V1;

fn wind() -> SpatialWindSpec {
    SpatialWindSpec {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.0,
            reference_height_m: 2.4,
            reference_speed_mps: 11.0,
        },
        spectrum: SpatialWindSpectrum::MannStress {
            parameters: MannParams {
                alpha_eps23: 0.02,
                ..MANN_TARGET_V1
            },
            reference_height_m: 3.0,
        },
        advection_mps: 11.0,
        correlation_times_s: vec![1.8; 32],
    }
}

#[test]
fn fitted_wind_reaches_flight_sampling_with_the_same_tensor_and_clock() {
    let spec = wind();
    let seed = 1903;
    let sim = SimLoop::init_with_distributed_spatial_wind(
        dec17_scenario(seed, PilotMode::FixedControls),
        spec.clone(),
    )
    .unwrap();
    let SpatialWindSpectrum::MannStress {
        parameters,
        reference_height_m,
    } = spec.spectrum
    else {
        unreachable!()
    };
    let fit = fit_amplitudes(seed, 32, &parameters, reference_height_m).unwrap();
    assert_eq!(sim.spatial_wind_stress_fit(), Some(&fit));
    assert!(fit.uw_fraction > 0.0);
    assert!(fit.target_uw_m2_s2 < 0.0);
    let basis = build_fitted_field(seed, &fit, &parameters, spec.advection_mps).unwrap();
    let mut field =
        EvolvingTurbulenceField::at_anchor(basis, seed, spec.correlation_times_s).unwrap();
    field.advance_to(0).unwrap();
    let before = sim.digest_hex();
    for [x, y, h] in [[0.0, 0.0, 0.0], [1.2, -2.3, 3.0], [-0.5, 4.0, 5.0]] {
        let mut expected = field.sample(-x, y, h).unwrap();
        expected.u[0] += spec.mean.speed(h);
        expected.grad[0][2] += spec.mean.dspeed_dh(h);
        let signs = [-1.0, 1.0, 1.0];
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
        }
    }
    assert_eq!(
        sim.digest_hex(),
        before,
        "observations cannot advance the atmosphere"
    );
}

#[test]
fn fitted_turbulence_changes_ground_rotor_and_distributed_lateral_motion() {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let mut quiet = wind();
    quiet.spectrum = SpatialWindSpectrum::Analytic {
        sigma_mps: 0.0,
        length_m: 1.77,
    };
    let mut calm = SimLoop::init_with_distributed_spatial_wind(scenario.clone(), quiet).unwrap();
    let mut turbulent = SimLoop::init_with_distributed_spatial_wind(scenario, wind()).unwrap();
    for sim in [&mut calm, &mut turbulent] {
        sim.phase = Phase::Airborne;
        sim.h_m = 20.0;
        sim.sync_spatial_air_velocity().unwrap();
    }
    let a = calm.step(None).unwrap();
    let b = turbulent.step(None).unwrap();
    assert_eq!(a.phase, Phase::Airborne);
    assert_eq!(b.phase, Phase::Airborne);
    assert_ne!(
        calm.spatial_ground_velocity_mps(),
        turbulent.spatial_ground_velocity_mps()
    );
    assert_ne!(a.omega_prop_rad_s.to_bits(), b.omega_prop_rad_s.to_bits());
    assert!(
        (a.p_rad_s - b.p_rad_s).abs() + (a.r_rad_s - b.r_rad_s).abs() > 1e-10,
        "the same calibrated field must drive the distributed load moments"
    );
}

#[test]
fn fitted_checkpoint_restores_cg_distributed_and_material_profile_runs() {
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
            for _ in 0..4 {
                sim.step(None).unwrap();
            }
            let bytes = sim.save_checkpoint().unwrap();
            assert_eq!(&bytes[..4], &7u32.to_le_bytes());
            let mut resumed = SimLoop::restore_checkpoint(sim.spec.clone(), &bytes).unwrap();
            assert_eq!(resumed.save_checkpoint().unwrap(), bytes);
            assert_eq!(
                resumed.spatial_wind_stress_fit(),
                sim.spatial_wind_stress_fit()
            );
            for _ in 0..5 {
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
fn a_resealed_checkpoint_cannot_substitute_different_mann_parameters() {
    let scenario = dec17_scenario(1903, PilotMode::FixedControls);
    let sim = SimLoop::init_with_spatial_wind(scenario.clone(), wind()).unwrap();
    let original = sim.save_checkpoint().unwrap();
    let marker = b"mann-stress-v1\0";
    let begin = original
        .windows(marker.len())
        .position(|v| v == marker)
        .unwrap()
        + marker.len();
    for field in 0..4 {
        let mut changed = original.clone();
        let offset = begin + 8 * field;
        let old = f64::from_le_bytes(changed[offset..offset + 8].try_into().unwrap());
        changed[offset..offset + 8].copy_from_slice(&(1.01 * old).to_le_bytes());
        let body_len = changed.len() - 32;
        let digest = hash_domain(
            "org.frankensim.wf.checkpoint-state.v7",
            &changed[..body_len],
        );
        changed[body_len..].copy_from_slice(digest.as_bytes());
        let error = SimLoop::restore_checkpoint(scenario.clone(), &changed)
            .err()
            .expect("model mismatch");
        assert_eq!(error.code, "checkpoint-malformed");
        assert!(error.message.contains("different scenario"));
    }
}

#[test]
fn recorded_fitted_wind_replays_without_an_external_fit_or_wind_specification() {
    use crate::replayenv::{ReplayEnvelope, record_replay_with_distributed_spatial_wind};
    let mut scenario = dec17_scenario(1903, PilotMode::FixedControls);
    scenario.max_ticks = 12;
    let mut direct = SimLoop::init_with_distributed_spatial_wind(scenario.clone(), wind()).unwrap();
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
