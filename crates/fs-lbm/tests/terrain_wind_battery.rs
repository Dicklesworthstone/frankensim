//! V-04c wind-over-terrain battery (bead `frankensim-wf-root-guzez.11.7`, E10.3).

use fs_exec::{CancelGate, TilePool};
use fs_lbm::terrain_wind::{
    TerrainWindConfig, V04C_RECEIPT_SCHEMA_V2, run_v04c_terrain_wind,
    run_v04c_terrain_wind_with_heightfield,
};

fn small_config() -> TerrainWindConfig {
    TerrainWindConfig {
        domain_size_m: [80.0, 40.0, 40.0],
        resolution: [16, 8, 8],
        num_steps: 40,
        tau: 0.8,
        declared_band_mps: 0.5,
        ..Default::default()
    }
}

#[test]
fn v04c_terrain_wind_runs_and_evaluates_discrepancy_bands() {
    let config = small_config();
    let receipt = run_v04c_terrain_wind(&config).expect("simulation succeeds");

    assert_eq!(receipt.schema_version, V04C_RECEIPT_SCHEMA_V2);
    assert!(receipt.mean_wind_rmse_mps.is_finite());
    assert!(receipt.maximum_speed_mps > 0.0);
    assert!(receipt.mean_velocity_mps[0] > 0.0);
    assert_eq!(receipt.cell_size_m, 5.0);
    assert_eq!(receipt.elapsed_time_s, 1.0);
    assert!((receipt.kinematic_viscosity_m2_s - 100.0).abs() < 1e-12);
    assert_eq!(
        receipt.within_declared_band,
        receipt.mean_wind_rmse_mps <= config.declared_band_mps
            && receipt.max_discrepancy_mps <= config.declared_band_mps
    );
    assert!(!receipt.run_manifest_digest.is_empty());
    assert!(!receipt.receipt_digest.is_empty());
}

#[test]
fn actual_time_viscosity_and_voxel_geometry_change_the_velocity() {
    let mut config = small_config();
    let flat = vec![0.0; 16 * 8];
    let pool = TilePool::for_host(1, 0);
    let gate = CancelGate::new();
    config.num_steps = 0;
    let rest = run_v04c_terrain_wind_with_heightfield(&config, &flat, &pool, &gate).unwrap();
    assert_eq!(rest.mean_velocity_mps, [0.0; 3]);
    assert_eq!(rest.maximum_speed_mps, 0.0);
    assert!(!rest.recirculation_detected);
    assert!(!rest.within_declared_band);

    config.num_steps = 40;
    let moving = run_v04c_terrain_wind_with_heightfield(&config, &flat, &pool, &gate).unwrap();
    assert!(moving.mean_velocity_mps[0] > 0.1);
    assert_ne!(rest.run_manifest_digest, moving.run_manifest_digest);

    config.tau = 1.1;
    let viscous = run_v04c_terrain_wind_with_heightfield(&config, &flat, &pool, &gate).unwrap();
    assert!((viscous.mean_velocity_mps[0] - moving.mean_velocity_mps[0]).abs() > 1e-5);
    assert_ne!(viscous.run_manifest_digest, moving.run_manifest_digest);

    let mut ridge = flat;
    for y in 0..8 {
        for x in 6..10 {
            ridge[y * 16 + x] = 10.0;
        }
    }
    let obstructed = run_v04c_terrain_wind_with_heightfield(&config, &ridge, &pool, &gate).unwrap();
    assert_eq!(viscous.sampled_cells - obstructed.sampled_cells, 4 * 8 * 2);
    assert!((viscous.mean_velocity_mps[0] - obstructed.mean_velocity_mps[0]).abs() > 1e-5);
    assert_ne!(viscous.run_manifest_digest, obstructed.run_manifest_digest);
}

#[test]
fn physical_velocity_scaling_leaves_the_same_native_trajectory() {
    let mut config = small_config();
    let flat = vec![0.0; 16 * 8];
    let pool = TilePool::for_host(1, 0);
    let gate = CancelGate::new();
    let original = run_v04c_terrain_wind_with_heightfield(&config, &flat, &pool, &gate).unwrap();
    config.inflow_speed_mps *= 2.0;
    config.timestep_s *= 0.5;
    let doubled = run_v04c_terrain_wind_with_heightfield(&config, &flat, &pool, &gate).unwrap();
    assert_eq!(doubled.elapsed_time_s, 0.5 * original.elapsed_time_s);
    assert_eq!(
        doubled.kinematic_viscosity_m2_s,
        2.0 * original.kinematic_viscosity_m2_s
    );
    for axis in 0..3 {
        assert!(
            (doubled.mean_velocity_mps[axis] - 2.0 * original.mean_velocity_mps[axis]).abs()
                < 1e-12
        );
    }
    assert!((doubled.maximum_speed_mps - 2.0 * original.maximum_speed_mps).abs() < 1e-12);
    assert_eq!(
        doubled.recirculation_detected,
        original.recirculation_detected
    );
}

#[test]
fn refuses_unsupported_grid_scaling_terrain_and_cancellation() {
    let base = small_config();
    let mut config = base.clone();
    config.domain_size_m[2] = 30.0;
    assert!(
        run_v04c_terrain_wind(&config)
            .unwrap_err()
            .contains("cubic cells")
    );
    config = base.clone();
    config.max_cells = 100;
    assert!(
        run_v04c_terrain_wind(&config)
            .unwrap_err()
            .contains("max_cells")
    );
    config = base.clone();
    config.timestep_s = 1.0;
    assert!(
        run_v04c_terrain_wind(&config)
            .unwrap_err()
            .contains("low-Mach")
    );
    config = base.clone();
    config.roughness_z0_m = f64::NAN;
    assert!(run_v04c_terrain_wind(&config).is_err());

    let pool = TilePool::for_host(1, 0);
    let gate = CancelGate::new();
    let mut heights = vec![0.0; 16 * 8];
    heights[1] = 10.0;
    assert!(
        run_v04c_terrain_wind_with_heightfield(&base, &heights, &pool, &gate)
            .unwrap_err()
            .contains("first interior")
    );
    heights[1] = f64::NAN;
    assert!(run_v04c_terrain_wind_with_heightfield(&base, &heights, &pool, &gate).is_err());
    gate.request();
    assert!(
        run_v04c_terrain_wind_with_heightfield(&base, &[], &pool, &gate)
            .unwrap_err()
            .contains("cancelled")
    );
}

#[test]
fn v04c_simulation_refuses_unstable_tau() {
    let config = TerrainWindConfig {
        tau: 0.50, // Unstable in LBM (tau must be > 0.5)
        ..Default::default()
    };

    let err = run_v04c_terrain_wind(&config).expect_err("must refuse unstable tau");
    assert!(err.contains("tau must be strictly greater than 0.5"));
}

#[test]
fn v04c_simulation_refuses_invalid_inflow_speed() {
    let config = TerrainWindConfig {
        inflow_speed_mps: 0.0,
        ..Default::default()
    };

    let err = run_v04c_terrain_wind(&config).expect_err("must refuse non-positive inflow");
    assert!(err.contains("inflow_speed_mps must be strictly positive"));
}
