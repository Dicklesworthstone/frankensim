//! WebGPU field-compute experiment battery (bead `frankensim-wf-root-guzez.8.8`, E7.6).

use fs_atmo::{Atmosphere, DEC17_AIR, FlatSiteLogLaw, TurbulenceField};
use fs_flyer::fieldsvc::{BoundSystem, FieldSourceStateV1, GridSpec};
use fs_flyer::webgpu_experiment::{
    run_webgpu_field_experiment, FieldBackendId, WebGpuExperimentConfig,
};

fn test_state() -> FieldSourceStateV1 {
    FieldSourceStateV1 {
        tick: 120,
        source_state_digest: "webgpu-test-state".into(),
        atmosphere: Atmosphere {
            mean: FlatSiteLogLaw {
                scenario_effective_z0_m: 5.0e-3,
                displacement_height_m: 0.02,
                reference_height_m: 10.0,
                reference_speed_mps: 8.0,
            },
            turbulence: TurbulenceField::build(1903, 8, 0.9, 20.0, 8.0).unwrap(),
            air: DEC17_AIR,
        },
        bound: Some(BoundSystem {
            gamma_m2ps: 8.0,
            tip_left_m: [0.0, -6.0, 2.5],
            tip_right_m: [0.0, 6.0, 2.5],
            trail_m: 40.0,
            core_m: 0.05,
        }),
        images_active: true,
    }
}

#[test]
fn webgpu_experiment_runs_and_evaluates_promotion_criteria() {
    let state = test_state();
    let config = WebGpuExperimentConfig {
        grid: GridSpec {
            origin_m: [-4.0, -4.0, 1.0],
            dx_m: 1.0,
            nx: 8,
            ny: 8,
            nz: 4,
        },
        speedup_threshold: 2.0,
        allow_render_readback: false,
        simulate_device_loss: false,
    };

    let receipt = run_webgpu_field_experiment(&config, &state).expect("experiment runs");

    assert_eq!(receipt.backend_evaluated, FieldBackendId::WebGpuExperiment);
    assert_eq!(receipt.points_sampled, 8 * 8 * 4);
    assert!(receipt.parity_verified);
    assert!(receipt.measured_speedup >= 2.0);
    assert!(receipt.promoted);
    assert!(receipt.no_physics_consumer_asserted);
    assert!(!receipt.receipt_digest.is_empty());
}

#[test]
fn webgpu_experiment_device_loss_fallback_drill() {
    let state = test_state();
    let config = WebGpuExperimentConfig {
        simulate_device_loss: true,
        ..Default::default()
    };

    let receipt = run_webgpu_field_experiment(&config, &state).expect("experiment with fallback runs");
    assert!(receipt.device_loss_fallback_verified);
}

#[test]
fn webgpu_experiment_refuses_promotion_if_readback_enabled() {
    let state = test_state();
    let config = WebGpuExperimentConfig {
        allow_render_readback: true, // Violates no-readback rule
        ..Default::default()
    };

    let receipt = run_webgpu_field_experiment(&config, &state).expect("experiment runs");
    assert!(!receipt.promoted, "must not promote when readback is enabled on render path");
}
