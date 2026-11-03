use super::*;

fn controls() -> TurningControls {
    TurningControls {
        minimum_step_s: 1e-10,
        maximum_angle_rad: core::f64::consts::FRAC_PI_4,
        maximum_refinements: 16,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

/// Analytical original midpoint solve for scalar mass 1 and Coulomb capacity 1.
/// vbar = v0 + h(F + f)/2 with f in -partial|vbar| is a disk resolvent.
fn midpoint(v0: [f64; 2], force: [f64; 2], h: f64) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let w = [v0[0] + 0.5 * h * force[0], v0[1] + 0.5 * h * force[1]];
    let length = w[0].hypot(w[1]);
    let scale = if length == 0.0 {
        0.0
    } else {
        (1.0 - 0.5 * h / length).max(0.0)
    };
    let mean = [scale * w[0], scale * w[1]];
    let after = [2.0 * mean[0] - v0[0], 2.0 * mean[1] - v0[1]];
    let friction = [
        (after[0] - v0[0]) / h - force[0],
        (after[1] - v0[1]) / h - force[1],
    ];
    (after, [h * mean[0], h * mean[1]], friction)
}

#[test]
fn original_midpoint_turn_preserves_both_momentum_and_coulomb_work() {
    let before = [1.2, 0.0];
    let force = [-0.6, 2.4];
    let (after, slip, friction) = midpoint(before, force, 1.0);
    close(after[0], 0.0, 1e-15);
    close(after[1], 1.6, 1e-15);
    close(slip[0], 0.6, 1e-15);
    close(slip[1], 0.8, 1e-15);
    close(friction[0], -0.6, 1e-15);
    close(friction[1], -0.8, 1e-15);
    let change = measure(before, after, 1e-10).unwrap().unwrap();
    assert!(change.transverse_resolved);
    close(change.angle_rad, core::f64::consts::FRAC_PI_2, 1e-15);
    let work = force[0] * slip[0] + force[1] * slip[1];
    let loss = slip[0].hypot(slip[1]);
    let kinetic = 0.5 * (after[0].powi(2) + after[1].powi(2) - before[0].powi(2));
    close(work, 1.56, 1e-15);
    close(loss, 1.0, 1e-15);
    close(kinetic, 0.56, 1e-15);
    close(work, kinetic + loss, 1e-15);
}

#[test]
fn bounded_refinement_keeps_the_actual_nonzero_turning_velocity() {
    let before = [1.2, 0.0];
    let force = [-0.6, 2.4];
    let c = controls();
    let mut h = 2.0;
    let mut refinements = 0;
    let (accepted, angle) = loop {
        let (after, _, _) = midpoint(before, force, h);
        let change = measure(before, after, 1e-10).unwrap().unwrap();
        if change.angle_rad <= c.maximum_angle_rad {
            break (after, change.angle_rad);
        }
        assert!(change.transverse_resolved);
        h = shorter_duration(h, c, refinements).unwrap();
        refinements += 1;
    };
    assert_eq!(refinements, 3);
    assert_eq!(h, 0.25);
    assert!(angle > 0.0 && angle < core::f64::consts::FRAC_PI_4);
    assert!(accepted[0] > 0.0 && accepted[1] > 0.0);
    assert!(accepted[0].hypot(accepted[1]) > 0.9);
    assert_eq!(accepted, midpoint(before, force, h).0);
}

#[test]
fn intermediate_original_response_exposes_turn_hidden_by_antiparallel_endpoints() {
    let before = [1.0, 0.0];
    let force = [-1.0, 1.0];
    let (upper, slip, friction) = midpoint(before, force, 2.0);
    assert_eq!(upper, [-1.0, 0.0]);
    assert_eq!(slip, [0.0, 0.0]);
    assert_eq!(friction, [0.0, -1.0]);
    let endpoint = measure(before, upper, 1e-10).unwrap().unwrap();
    assert!(!endpoint.transverse_resolved);
    assert_eq!(endpoint.angle_rad, core::f64::consts::PI);
    // A native stopping probe must inspect its actual full velocity, not only
    // its scalar projection. It can now hand this trial to angular refinement.
    let (intermediate, _, _) = midpoint(before, force, 1.0);
    close(intermediate[0], -core::f64::consts::FRAC_1_SQRT_2, 1e-15);
    close(
        intermediate[1],
        1.0 - core::f64::consts::FRAC_1_SQRT_2,
        1e-15,
    );
    let change = measure(before, intermediate, 1e-10).unwrap().unwrap();
    assert!(change.transverse_resolved);
    assert!(change.angle_rad > controls().maximum_angle_rad);
    assert!(intermediate[0].hypot(intermediate[1]) > 0.7);
}

#[test]
fn relative_frame_rotation_and_velocity_units_leave_the_angle_unchanged() {
    let before = [3.0, 0.0];
    let after = [1.2, 1.6];
    let expected = measure(before, after, 1e-10).unwrap().unwrap();
    for scale in [1e-200, 1.0, 1e200] {
        // A 3-4-5 orthogonal frame change; products of raw large speeds overflow.
        let rotated = |v: [f64; 2]| {
            [
                scale * (0.6 * v[0] - 0.8 * v[1]),
                scale * (0.8 * v[0] + 0.6 * v[1]),
            ]
        };
        let actual = measure(rotated(before), rotated(after), scale * 1e-10)
            .unwrap()
            .unwrap();
        close(actual.angle_rad, expected.angle_rad, 1e-15);
        assert_eq!(actual.transverse_resolved, expected.transverse_resolved);
    }
}

#[test]
fn both_physical_velocity_bands_protect_stops_and_collinear_reversals() {
    assert_eq!(measure([1.0, 0.0], [0.0, 0.0], 1e-8).unwrap(), None);
    assert_eq!(measure([1e-8, 0.0], [1.0, 0.0], 1e-8).unwrap(), None);
    assert_eq!(measure([1.0, 0.0], [0.0, -1e-8], 1e-8).unwrap(), None);
    let collinear = measure([1.0, 0.0], [-2.0, -0.0], 1e-8).unwrap().unwrap();
    assert_eq!(collinear.angle_rad, core::f64::consts::PI);
    assert!(!collinear.transverse_resolved);
    // The actual collinear reversal must remain collinear even when raw speed
    // products overflow and the physical tolerance is far below one speed ulp.
    let large_collinear = measure([3e200, 4e200], [-6e200, -8e200], 1e-8)
        .unwrap()
        .unwrap();
    assert_eq!(large_collinear.angle_rad, core::f64::consts::PI);
    assert!(!large_collinear.transverse_resolved);
    let unresolved = measure([1.0, 0.0], [-1.0, 1.5e-8], 1e-8).unwrap().unwrap();
    assert!(!unresolved.transverse_resolved);
    let resolved = measure([1.0, 0.0], [-1.0, 3e-8], 1e-8).unwrap().unwrap();
    assert!(resolved.transverse_resolved);
}

#[test]
fn native_trial_duration_and_refinement_limits_are_explicit() {
    let mut c = controls();
    c.minimum_step_s = 0.1;
    c.maximum_refinements = 2;
    assert_eq!(shorter_duration(0.3, c, 0).unwrap(), 0.15);
    assert_eq!(shorter_duration(0.15, c, 1).unwrap(), 0.1);
    assert!(matches!(
        shorter_duration(0.3, c, 2),
        Err(TurningError::Refinement {
            refinements: 2,
            maximum: 2
        })
    ));
    assert!(matches!(
        shorter_duration(0.1, c, 0),
        Err(TurningError::MinimumStep { .. })
    ));
    c.minimum_step_s = f64::from_bits(1);
    assert_eq!(
        shorter_duration(f64::from_bits(2), c, 0).unwrap(),
        c.minimum_step_s
    );
    assert!(matches!(
        shorter_duration(c.minimum_step_s, c, 1),
        Err(TurningError::MinimumStep { .. })
    ));
}

#[test]
fn malformed_angles_and_unrepresentable_speeds_refuse() {
    for angle in [0.0, -0.1, core::f64::consts::FRAC_PI_2, f64::NAN] {
        let mut c = controls();
        c.maximum_angle_rad = angle;
        assert!(matches!(validate(c), Err(TurningError::Invalid(_))));
    }
    assert!(matches!(
        measure([f64::INFINITY, 0.0], [1.0, 0.0], 1e-8),
        Err(TurningError::NonFinite(_))
    ));
    assert!(matches!(
        measure([f64::MAX, f64::MAX], [1.0, 0.0], 1e-8),
        Err(TurningError::NonFinite(_))
    ));
    assert!(matches!(
        measure([1.0, 0.0], [1.0, 0.0], 0.0),
        Err(TurningError::Invalid(_))
    ));
}
