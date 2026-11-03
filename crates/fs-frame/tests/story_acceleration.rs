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
