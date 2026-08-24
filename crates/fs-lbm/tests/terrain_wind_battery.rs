//! V-04c wind-over-terrain battery (bead `frankensim-wf-root-guzez.11.7`, E10.3).

use fs_lbm::terrain_wind::{run_v04c_terrain_wind, TerrainWindConfig, V04C_RECEIPT_SCHEMA_V1};

#[test]
fn v04c_terrain_wind_runs_and_evaluates_discrepancy_bands() {
    let config = TerrainWindConfig {
        inflow_speed_mps: 8.0,
        roughness_z0_m: 0.01,
        resolution: [16, 16, 12],
        num_steps: 50,
        tau: 0.60,
        declared_band_mps: 0.5,
        ..Default::default()
    };

    let receipt = run_v04c_terrain_wind(&config).expect("simulation succeeds");

    assert_eq!(receipt.schema_version, V04C_RECEIPT_SCHEMA_V1);
    assert!(receipt.mean_wind_rmse_mps < config.declared_band_mps);
    assert!(receipt.within_declared_band);
    assert!(!receipt.run_manifest_digest.is_empty());
    assert!(!receipt.receipt_digest.is_empty());
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
