//! SmoothedTangentPlane test suite (bead `frankensim-wf-root-guzez.5.11`, E4.4b, V-06c).

use fs_wing::smoothed_plane::{
    SmoothedTangentPlane, SmoothedTangentPlaneConfig, DEFAULT_MAX_SLOPE_RAD,
};

#[test]
fn smoothed_plane_slope_gate_admits_at_cap_and_refuses_at_cap_plus_epsilon() {
    let mut plane = SmoothedTangentPlane::new(SmoothedTangentPlaneConfig::default(), 0.0);

    // Flat horizontal terrain admits
    let normal_flat = [0.0, 0.0, -1.0];
    let receipt = plane.update(1.0, normal_flat, 500.0, 0.05).expect("flat admits");
    assert_eq!(receipt.claim_class, "EstimateOnly");
    assert!(receipt.slope_rad < 1e-6);

    // Exactly at slope cap admits (nx = sin(cap), nz = -cos(cap))
    let cap = DEFAULT_MAX_SLOPE_RAD;
    let normal_at_cap = [cap.sin(), 0.0, -cap.cos()];
    let res_at_cap = plane.update(1.0, normal_at_cap, 500.0, 0.05);
    assert!(res_at_cap.is_ok(), "slope at cap should admit");

    // Slope beyond cap refuses with typed code
    let beyond_cap = DEFAULT_MAX_SLOPE_RAD + 0.02;
    let normal_beyond = [beyond_cap.sin(), 0.0, -beyond_cap.cos()];
    let err = plane.update(1.0, normal_beyond, 500.0, 0.05).expect_err("slope beyond cap must refuse");
    assert_eq!(err.code, "smoothed-plane-slope-exceeded");
}

#[test]
fn smoothed_plane_hysteresis_and_filter_smoothing() {
    let config = SmoothedTangentPlaneConfig {
        hysteresis_band_m: 0.1,
        filter_cutoff_hz: 1.0,
        ..Default::default()
    };
    let mut plane = SmoothedTangentPlane::new(config, 0.0);

    // Micro-perturbation within hysteresis deadband does not move the plane
    let normal = [0.0, 0.0, -1.0];
    let r1 = plane.update(0.05, normal, 100.0, 0.01).expect("update");
    assert_eq!(r1.plane_z_m, 0.0, "within deadband height should not change");
    assert_eq!(r1.artificial_boundary_power_w, 0.0);

    // Significant step changes plane smoothly
    let mut last_z = 0.0;
    for _ in 0..50 {
        let r = plane.update(2.0, normal, 100.0, 0.01).expect("step update");
        assert!(r.plane_z_m >= last_z, "plane height increases monotonically toward target");
        assert!(r.artificial_boundary_power_w > 0.0 || r.plane_z_m > 1.9);
        last_z = r.plane_z_m;
    }
    assert!(last_z > 1.0 && last_z <= 2.0, "approaches target elevation");
}

#[test]
fn smoothed_plane_state_checkpoint_and_replay() {
    let mut plane = SmoothedTangentPlane::new(SmoothedTangentPlaneConfig::default(), 0.0);
    let normal = [0.0, 0.0, -1.0];

    for _ in 0..10 {
        plane.update(1.5, normal, 200.0, 0.01).expect("update");
    }

    let checkpoint = *plane.state();
    let r_cont = plane.update(2.5, normal, 200.0, 0.01).expect("cont update");

    // Rewind and replay from checkpoint
    let mut replayed_plane = SmoothedTangentPlane::new(SmoothedTangentPlaneConfig::default(), 0.0);
    replayed_plane.restore_state(checkpoint);
    let r_replayed = replayed_plane.update(2.5, normal, 200.0, 0.01).expect("replayed update");

    assert_eq!(
        r_cont.plane_z_m.to_bits(),
        r_replayed.plane_z_m.to_bits(),
        "replayed plane height must match bitwise"
    );
    assert_eq!(
        r_cont.cumulative_artificial_work_j.to_bits(),
        r_replayed.cumulative_artificial_work_j.to_bits(),
        "replayed artificial work must match bitwise"
    );
}
