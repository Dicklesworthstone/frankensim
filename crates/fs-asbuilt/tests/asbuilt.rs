//! Battery for as-built ingestion (addendum Proposal 11). Covers rigid
//! registration recovery (exact + noisy), fiducial well-posedness (too-few,
//! collinear), the R8 signal-vs-noise gate, and the validated-color as-built
//! δ anchored to the calibration certificate.

use fs_asbuilt::{
    Color, Fiducial, Point2, RegError, Registration, as_built_diff, register, well_posed,
};

/// Apply a ground-truth rigid transform to a design point (for building scans).
fn xform(p: Point2, theta: f64, tx: f64, ty: f64) -> Point2 {
    let (s, c) = theta.sin_cos();
    Point2::new(c * p.x - s * p.y + tx, s * p.x + c * p.y + ty)
}

fn triangle() -> [Point2; 3] {
    [
        Point2::new(0.0, 0.0),
        Point2::new(2.0, 0.0),
        Point2::new(0.0, 2.0),
    ]
}

#[test]
fn registration_recovers_a_known_rigid_transform() {
    let (theta, tx, ty) = (std::f64::consts::FRAC_PI_6, 5.0, 2.0); // 30 degrees
    let fids: Vec<Fiducial> = triangle()
        .iter()
        .map(|&d| Fiducial::new(d, xform(d, theta, tx, ty)))
        .collect();
    let reg = register(&fids).unwrap();
    assert!(
        (reg.rotation_rad - theta).abs() < 1e-9,
        "theta {}",
        reg.rotation_rad
    );
    assert!((reg.tx - tx).abs() < 1e-9 && (reg.ty - ty).abs() < 1e-9);
    assert!(reg.residual_rms < 1e-9, "residual {}", reg.residual_rms);
    // and it maps a design point onto its scanned location.
    let p = Point2::new(1.0, 1.0);
    let mapped = reg.apply(p);
    let truth = xform(p, theta, tx, ty);
    assert!((mapped.x - truth.x).abs() < 1e-9 && (mapped.y - truth.y).abs() < 1e-9);
}

#[test]
fn noisy_measurements_carry_a_positive_residual() {
    let (theta, tx, ty) = (0.2, 1.0, -1.0);
    let noise = [(0.01, -0.02), (-0.015, 0.01), (0.02, 0.005)];
    let fids: Vec<Fiducial> = triangle()
        .iter()
        .zip(noise)
        .map(|(&d, (nx, ny))| {
            let m = xform(d, theta, tx, ty);
            Fiducial::new(d, Point2::new(m.x + nx, m.y + ny))
        })
        .collect();
    let reg = register(&fids).unwrap();
    // the registration error is carried forward, not discarded.
    assert!(reg.residual_rms > 0.0 && reg.residual_rms < 0.1);
}

#[test]
fn registration_is_ill_posed_without_enough_non_collinear_fiducials() {
    // too few.
    let two = [
        Fiducial::new(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0)),
        Fiducial::new(Point2::new(1.0, 0.0), Point2::new(2.0, 1.0)),
    ];
    assert!(matches!(
        register(&two),
        Err(RegError::TooFewFiducials { have: 2, need: 3 })
    ));
    // collinear design points (all on the x-axis) are rank-deficient.
    let collinear: Vec<Fiducial> = [0.0, 1.0, 2.0]
        .iter()
        .map(|&x| Fiducial::new(Point2::new(x, 0.0), Point2::new(x + 0.3, 5.0)))
        .collect();
    assert_eq!(register(&collinear), Err(RegError::CollinearFiducials));
}

#[test]
fn the_r8_gate_rejects_registration_below_the_noise_floor() {
    // signal (certified deviation) 0.5 above the registration residual 0.01 -> ok.
    let sharp = Registration {
        rotation_rad: 0.0,
        tx: 0.0,
        ty: 0.0,
        residual_rms: 0.01,
    };
    assert!(well_posed(&sharp, 0.5));
    // registration residual 0.6 exceeds the 0.5 deviation being certified -> R8 kill.
    let blurry = Registration {
        rotation_rad: 0.0,
        tx: 0.0,
        ty: 0.0,
        residual_rms: 0.6,
    };
    assert!(!well_posed(&blurry, 0.5));
    // a non-positive certified deviation is never well-posed.
    assert!(!well_posed(&sharp, 0.0));
}

#[test]
fn the_as_built_diff_is_validated_and_anchored_to_the_cert() {
    let reg = Registration {
        rotation_rad: 0.0,
        tx: 0.0,
        ty: 0.0,
        residual_rms: 0.0,
    };
    let design = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 1.0)];
    let scanned = vec![Point2::new(0.0, 0.1), Point2::new(1.0, 1.0)];
    let diff = as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "metrology-cal-2026").unwrap();
    assert!((diff.max_deviation - 0.1).abs() < 1e-12);
    assert!(diff.within_tolerance); // 0.1 <= 0.2
    assert!(diff.above_noise_floor); // 0.1 > 0.05 (distinguishable from noise)
    // validated color anchored to the calibration certificate.
    match &diff.color {
        Color::Validated { dataset, .. } => assert_eq!(dataset, "metrology-cal-2026"),
        other => panic!("expected validated, got {other:?}"),
    }
}

#[test]
fn a_deviation_below_the_noise_floor_is_flagged() {
    let reg = Registration {
        rotation_rad: 0.0,
        tx: 0.0,
        ty: 0.0,
        residual_rms: 0.0,
    };
    let design = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 1.0)];
    let scanned = vec![Point2::new(0.0, 0.01), Point2::new(1.0, 1.0)];
    // deviation 0.01 is below the 0.05 measurement noise floor.
    let diff = as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "cal").unwrap();
    assert!(!diff.above_noise_floor);
}

#[test]
fn as_built_diff_rejects_malformed_input() {
    let reg = Registration {
        rotation_rad: 0.0,
        tx: 0.0,
        ty: 0.0,
        residual_rms: 0.0,
    };
    assert_eq!(
        as_built_diff(&reg, &[], &[], 0.1, 0.01, "c"),
        Err(RegError::Empty)
    );
    assert!(matches!(
        as_built_diff(&reg, &[Point2::new(0.0, 0.0)], &[], 0.1, 0.01, "c"),
        Err(RegError::LengthMismatch { .. })
    ));
}

#[test]
fn registration_is_deterministic() {
    let fids: Vec<Fiducial> = triangle()
        .iter()
        .map(|&d| Fiducial::new(d, xform(d, 0.3, 1.0, 2.0)))
        .collect();
    assert_eq!(register(&fids), register(&fids));
}
