//! Nonlinear fiber forces determine acceleration even below drift resolution.

use fs_frame::{GroundMotion, HistoryLimits, StoryFrame, StoryParams};

#[test]
fn g1_tiny_step_retains_fiber_equilibrium_and_velocity() {
    let params = StoryParams::default();
    let mut frame = StoryFrame::new(params);
    frame.x = 1e-3;
    let restoring = frame.restoring(frame.x).0;
    let damping = 2.0 * params.zeta * (frame.initial_stiffness() * params.mass).sqrt();
    frame.a = -restoring / params.mass;
    let before = frame.clone();
    let h = 1e-12;
    let ground = 0.75;
    let expected_acceleration =
        (-ground - restoring / params.mass - damping / params.mass * (0.5 * h * before.a))
            / (1.0 + damping / params.mass * 0.5 * h);
    let expected_velocity = 0.5 * h * (before.a + expected_acceleration);
    let response = frame
        .run_checked(
            GroundMotion::new(&[ground], h),
            HistoryLimits::new(1, 8, 1e-20, 1e-8),
        )
        .expect("finite physical acceleration must survive a sub-ULP drift increment");
    assert_eq!(frame.x, before.x);
    assert!((frame.a - expected_acceleration).abs() < 2e-15);
    assert!((frame.v / expected_velocity - 1.0).abs() < 2e-14);
    assert_eq!(response.restoring_shear_n[0], restoring);
    assert!(response.max_abs_equilibrium_residual_n < 1e-8);
}

/// The viscous law is based on virgin stiffness, not the unloading tangent
/// reached at an arbitrary record or checkpoint boundary.
#[test]
fn g3_hysteretic_history_keeps_initial_damping_across_record_chunks() {
    let params = StoryParams::default();
    let dt = 0.02;
    // Exercise retained material history within the concrete strain domain.
    let motion: Vec<_> = (0..160)
        .map(|i| 0.5 * (1.2 * f64::from(i) * dt * std::f64::consts::TAU).sin())
        .collect();
    let limits = HistoryLimits::new(motion.len(), 40, 1e-12, 1e-2);
    let mut complete = StoryFrame::new(params);
    let initial = complete.initial_stiffness();
    let response = complete
        .run_checked(GroundMotion::new(&motion, dt), limits)
        .unwrap();
    assert!(response.peak_abs_displacement_m / params.h > 0.002);
    assert_eq!(complete.initial_stiffness().to_bits(), initial.to_bits());

    for chunk_size in [1, 7, 40] {
        let mut resumed = StoryFrame::new(params);
        let mut displacement = Vec::new();
        let mut shear = Vec::new();
        for chunk in motion.chunks(chunk_size) {
            let part = resumed
                .run_checked(GroundMotion::new(chunk, dt), limits)
                .unwrap();
            displacement.extend(part.displacement_m);
            shear.extend(part.restoring_shear_n);
            assert_eq!(resumed.initial_stiffness().to_bits(), initial.to_bits());
        }
        assert_eq!(displacement, response.displacement_m);
        assert_eq!(shear, response.restoring_shear_n);
        assert_eq!(
            [
                resumed.x.to_bits(),
                resumed.v.to_bits(),
                resumed.a.to_bits()
            ],
            [
                complete.x.to_bits(),
                complete.v.to_bits(),
                complete.a.to_bits()
            ]
        );
    }
}

#[test]
fn g3_legacy_history_and_cloned_restart_keep_the_same_viscous_law() {
    let params = StoryParams::default();
    let dt = 0.02;
    let motion: Vec<_> = (0..160)
        .map(|i| 0.5 * (1.2 * f64::from(i) * dt * std::f64::consts::TAU).sin())
        .collect();
    let mut uninterrupted = StoryFrame::new(params);
    let expected = uninterrupted.run(&motion, dt);
    let mut prefix = StoryFrame::new(params);
    let mut actual = prefix.run(&motion[..73], dt);
    let mut resumed = prefix.clone();
    actual.extend(resumed.run(&motion[73..], dt));
    assert_eq!(actual, expected);
    assert_eq!(resumed.v.to_bits(), uninterrupted.v.to_bits());
    assert_eq!(resumed.a.to_bits(), uninterrupted.a.to_bits());
}
