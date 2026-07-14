//! Battery for as-built ingestion (addendum Proposal 11). Covers rigid
//! registration recovery (exact + noisy), fiducial well-posedness (too-few,
//! collinear), the R8 signal-vs-noise gate, and the as-built δ retained as an
//! estimated candidate until calibration authority exists.

use fs_asbuilt::{
    AS_BUILT_POLL_POLICY_VERSION, AS_BUILT_POLL_STRIDE_POINTS, Color, Fiducial, Point2, RegError,
    Registration, as_built_diff, register, well_posed,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(cancelled: bool, mode: ExecMode, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xA5B0_117,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            budget,
            mode,
        );
        f(&cx)
    })
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(false, ExecMode::Deterministic, Budget::INFINITE, f)
}

/// Apply a ground-truth rigid transform to a design point (for building scans).
fn xform(p: Point2, theta: f64, tx: f64, ty: f64) -> Point2 {
    let (s, c) = theta.sin_cos();
    point(c * p.x() - s * p.y() + tx, s * p.x() + c * p.y() + ty)
}

fn point(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).expect("finite test point")
}

fn registration(rotation: f64, tx: f64, ty: f64, residual: f64) -> Registration {
    Registration::new(rotation, tx, ty, residual).expect("valid test registration")
}

fn triangle() -> [Point2; 3] {
    [point(0.0, 0.0), point(2.0, 0.0), point(0.0, 2.0)]
}

#[test]
fn registration_recovers_a_known_rigid_transform() {
    let (theta, tx, ty) = (std::f64::consts::FRAC_PI_6, 5.0, 2.0); // 30 degrees
    let fids: Vec<Fiducial> = triangle()
        .iter()
        .map(|&d| Fiducial::new(d, xform(d, theta, tx, ty)))
        .collect();
    let reg = with_default_cx(|cx| register(&fids, cx)).unwrap();
    assert!(
        (reg.rotation_rad() - theta).abs() < 1e-9,
        "theta {}",
        reg.rotation_rad()
    );
    assert!((reg.tx() - tx).abs() < 1e-9 && (reg.ty() - ty).abs() < 1e-9);
    assert!(reg.residual_rms() < 1e-9, "residual {}", reg.residual_rms());
    // and it maps a design point onto its scanned location.
    let p = point(1.0, 1.0);
    let mapped = reg.apply(p).expect("finite mapped point");
    let truth = xform(p, theta, tx, ty);
    assert!((mapped.x() - truth.x()).abs() < 1e-9 && (mapped.y() - truth.y()).abs() < 1e-9);
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
            Fiducial::new(d, point(m.x() + nx, m.y() + ny))
        })
        .collect();
    let reg = with_default_cx(|cx| register(&fids, cx)).unwrap();
    // the registration error is carried forward, not discarded.
    assert!(reg.residual_rms() > 0.0 && reg.residual_rms() < 0.1);
}

#[test]
fn registration_is_ill_posed_without_enough_non_collinear_fiducials() {
    // too few.
    let two = [
        Fiducial::new(point(0.0, 0.0), point(1.0, 1.0)),
        Fiducial::new(point(1.0, 0.0), point(2.0, 1.0)),
    ];
    assert!(matches!(
        with_default_cx(|cx| register(&two, cx)),
        Err(RegError::TooFewFiducials { have: 2, need: 3 })
    ));
    // collinear design points (all on the x-axis) are rank-deficient.
    let collinear: Vec<Fiducial> = [0.0, 1.0, 2.0]
        .iter()
        .map(|&x| Fiducial::new(point(x, 0.0), point(x + 0.3, 5.0)))
        .collect();
    assert_eq!(
        with_default_cx(|cx| register(&collinear, cx)),
        Err(RegError::CollinearFiducials)
    );
}

#[test]
fn registration_rank_gate_is_scale_invariant_for_small_geometry() {
    let scale = 1.0e-9;
    let design = [point(0.0, 0.0), point(scale, 0.0), point(0.0, scale)];
    let fiducials: Vec<_> = design
        .iter()
        .copied()
        .map(|datum| Fiducial::new(datum, datum))
        .collect();

    let registration = with_default_cx(|cx| register(&fiducials, cx))
        .expect("small non-collinear triangle is rank two");
    assert_eq!(registration.rotation_rad().to_bits(), 0.0f64.to_bits());
    assert_eq!(registration.tx().to_bits(), 0.0f64.to_bits());
    assert_eq!(registration.ty().to_bits(), 0.0f64.to_bits());
    assert!(registration.residual_rms() <= f64::MIN_POSITIVE);
}

#[test]
fn the_r8_gate_rejects_registration_below_the_noise_floor() {
    // signal (certified deviation) 0.5 above the registration residual 0.01 -> ok.
    let sharp = registration(0.0, 0.0, 0.0, 0.01);
    assert!(well_posed(&sharp, 0.5));
    // registration residual 0.6 exceeds the 0.5 deviation being certified -> R8 kill.
    let blurry = registration(0.0, 0.0, 0.0, 0.6);
    assert!(!well_posed(&blurry, 0.5));
    // a non-positive certified deviation is never well-posed.
    assert!(!well_posed(&sharp, 0.0));
}

#[test]
fn the_as_built_diff_is_an_estimated_candidate_with_a_proposed_regime() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = vec![point(0.0, 0.0), point(1.0, 1.0)];
    let scanned = vec![point(0.0, 0.1), point(1.0, 1.0)];
    let diff = with_default_cx(|cx| {
        as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "metrology-cal-2026", cx)
    })
    .unwrap();
    assert!((diff.max_deviation() - 0.1).abs() < 1e-12);
    assert!(diff.within_tolerance()); // 0.1 + 0.05 <= 0.2 (advisory one-dispersion screen)
    assert!(diff.above_noise_floor()); // 0.1 > combined dispersion 0.05
    assert_eq!(
        diff.proposed_regime().bound("measurement_noise"),
        Some((0.0, 0.05))
    );
    assert_eq!(
        diff.proposed_regime().bound("design_tolerance"),
        Some((0.0, 0.2))
    );
    match diff.color() {
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            assert!(estimator.starts_with("asbuilt-diff-v3:"));
            assert_eq!(dispersion.to_bits(), 0.05f64.to_bits());
        }
        other => panic!("expected estimated candidate, got {other:?}"),
    }
}

#[test]
fn a_deviation_below_the_noise_floor_is_flagged() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = vec![point(0.0, 0.0), point(1.0, 1.0)];
    let scanned = vec![point(0.0, 0.01), point(1.0, 1.0)];
    // deviation 0.01 is below the 0.05 measurement noise floor.
    let diff =
        with_default_cx(|cx| as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "cal", cx)).unwrap();
    assert!(!diff.above_noise_floor());
}

#[test]
fn registration_residual_is_included_in_advisory_deviation_screens() {
    let reg = registration(0.0, 0.0, 0.0, 0.08);
    let design = [point(0.0, 0.0)];
    let scanned = [point(0.0, 0.1)];
    let diff = with_default_cx(|cx| as_built_diff(&reg, &design, &scanned, 0.2, 0.06, "cal", cx))
        .expect("finite as-built inputs");

    assert!((diff.max_deviation() - 0.1).abs() < 1e-12);
    assert!(!diff.within_tolerance());
    assert!(!diff.above_noise_floor());
    match diff.color() {
        Color::Estimated { dispersion, .. } => {
            assert!((*dispersion - 0.14).abs() < 1e-12);
        }
        other => panic!("expected estimated candidate, got {other:?}"),
    }
}

#[test]
fn as_built_diff_rejects_malformed_input() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    assert_eq!(
        with_default_cx(|cx| as_built_diff(&reg, &[], &[], 0.1, 0.01, "c", cx)),
        Err(RegError::Empty)
    );
    assert!(matches!(
        with_default_cx(|cx| { as_built_diff(&reg, &[point(0.0, 0.0)], &[], 0.1, 0.01, "c", cx) }),
        Err(RegError::LengthMismatch { .. })
    ));
}

#[test]
fn registration_is_deterministic() {
    let fids: Vec<Fiducial> = triangle()
        .iter()
        .map(|&d| Fiducial::new(d, xform(d, 0.3, 1.0, 2.0)))
        .collect();
    let first = with_default_cx(|cx| register(&fids, cx));
    let replay = with_default_cx(|cx| register(&fids, cx));
    assert_eq!(first, replay);
}

fn estimator_identity(diff: &fs_asbuilt::AsBuiltDiff) -> &str {
    match diff.color() {
        Color::Estimated { estimator, .. } => estimator,
        other => panic!("default diff must remain estimated, got {other:?}"),
    }
}

fn fixture_diff(mode: ExecMode, budget: Budget) -> fs_asbuilt::AsBuiltDiff {
    let reg = registration(0.1, 2.0, -3.0, 0.01);
    let design = [point(0.0, 0.0), point(1.0, 1.0), point(-2.0, 3.0)];
    let scanned = [
        reg.apply(design[0]).unwrap(),
        point(reg.apply(design[1]).unwrap().x() + 0.02, 0.0),
        reg.apply(design[2]).unwrap(),
    ];
    with_cx(false, mode, budget, |cx| {
        as_built_diff(
            &reg,
            &design,
            &scanned,
            10.0,
            0.05,
            "identity-context-fixture",
            cx,
        )
    })
    .expect("valid context-bound diff")
}

fn numeric_signature(diff: &fs_asbuilt::AsBuiltDiff) -> (Vec<u64>, u64, bool, bool, u64) {
    let dispersion = match diff.color() {
        Color::Estimated { dispersion, .. } => dispersion.to_bits(),
        other => panic!("default diff must remain estimated, got {other:?}"),
    };
    (
        diff.deviations()
            .iter()
            .map(|value| value.to_bits())
            .collect(),
        diff.max_deviation().to_bits(),
        diff.within_tolerance(),
        diff.above_noise_floor(),
        dispersion,
    )
}

#[test]
fn g4_pre_cancelled_entry_points_report_exact_zero_progress() {
    let fiducials: Vec<_> = triangle()
        .into_iter()
        .map(|datum| Fiducial::new(datum, datum))
        .collect();
    assert_eq!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            register(&fiducials, cx)
        },),
        Err(RegError::Cancelled {
            phase: "register.initial",
            completed_work: 0,
            planned_work: 12,
        })
    );

    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = [point(0.0, 0.0), point(1.0, 1.0)];
    assert_eq!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            as_built_diff(&reg, &design, &design, 0.1, 0.01, "cal", cx)
        },),
        Err(RegError::Cancelled {
            phase: "as-built-diff.initial",
            completed_work: 0,
            planned_work: 6,
        })
    );
}

#[test]
fn g5_identity_binds_mode_and_every_budget_field_without_changing_numerics() {
    let baseline = fixture_diff(ExecMode::Deterministic, Budget::new());
    let baseline_identity = estimator_identity(&baseline).to_owned();
    let baseline_numeric = numeric_signature(&baseline);
    let variants = [
        fixture_diff(ExecMode::Fast, Budget::new()),
        fixture_diff(ExecMode::Deterministic, Budget::with_deadline_at_ns(17)),
        fixture_diff(ExecMode::Deterministic, Budget::new().with_poll_quota(31)),
        fixture_diff(ExecMode::Deterministic, Budget::new().with_cost_quota(47)),
        fixture_diff(ExecMode::Deterministic, Budget::new().with_priority(199)),
    ];

    for variant in variants {
        assert_ne!(estimator_identity(&variant), baseline_identity);
        assert_eq!(numeric_signature(&variant), baseline_numeric);
    }

    for (left, right) in [
        (
            fixture_diff(ExecMode::Deterministic, Budget::with_deadline_at_ns(17)),
            fixture_diff(ExecMode::Deterministic, Budget::with_deadline_at_ns(18)),
        ),
        (
            fixture_diff(ExecMode::Deterministic, Budget::new().with_poll_quota(31)),
            fixture_diff(ExecMode::Deterministic, Budget::new().with_poll_quota(32)),
        ),
        (
            fixture_diff(ExecMode::Deterministic, Budget::new().with_cost_quota(47)),
            fixture_diff(ExecMode::Deterministic, Budget::new().with_cost_quota(48)),
        ),
        (
            fixture_diff(ExecMode::Deterministic, Budget::new().with_priority(199)),
            fixture_diff(ExecMode::Deterministic, Budget::new().with_priority(200)),
        ),
    ] {
        assert_ne!(estimator_identity(&left), estimator_identity(&right));
        assert_eq!(numeric_signature(&left), baseline_numeric);
        assert_eq!(numeric_signature(&right), baseline_numeric);
    }
}

#[test]
fn g5_retained_identity_declares_the_v1_256_point_poll_policy() {
    assert_eq!(AS_BUILT_POLL_POLICY_VERSION, 1);
    assert_eq!(AS_BUILT_POLL_STRIDE_POINTS, 256);
    let diff = fixture_diff(ExecMode::Deterministic, Budget::INFINITE);
    assert!(estimator_identity(&diff).starts_with("asbuilt-diff-v3:"));
}

#[test]
fn public_geometry_and_registration_construction_refuse_non_finite_values() {
    for (x, y, field) in [
        (f64::NAN, 0.0, "point.x"),
        (f64::INFINITY, 0.0, "point.x"),
        (0.0, f64::NEG_INFINITY, "point.y"),
    ] {
        assert_eq!(Point2::new(x, y), Err(RegError::NonFinite { field }));
    }

    for (rotation, tx, ty, residual, field) in [
        (f64::NAN, 0.0, 0.0, 0.0, "registration.rotation_rad"),
        (0.0, f64::INFINITY, 0.0, 0.0, "registration.tx"),
        (0.0, 0.0, f64::NEG_INFINITY, 0.0, "registration.ty"),
        (0.0, 0.0, 0.0, f64::NAN, "registration.residual_rms"),
    ] {
        assert_eq!(
            Registration::new(rotation, tx, ty, residual),
            Err(RegError::NonFinite { field })
        );
    }
    assert_eq!(
        Registration::new(0.0, 0.0, 0.0, -0.01),
        Err(RegError::Negative {
            field: "registration.residual_rms"
        })
    );

    let overflowing = Registration::new(0.0, f64::MAX, 0.0, 0.0).unwrap();
    assert!(matches!(
        overflowing.apply(point(f64::MAX, 0.0)),
        Err(RegError::NonFinite { field: "point.x" })
    ));
}

#[test]
fn public_numeric_values_canonicalize_signed_zero() {
    let point = point(-0.0, -0.0);
    assert_eq!(point.x().to_bits(), 0.0f64.to_bits());
    assert_eq!(point.y().to_bits(), 0.0f64.to_bits());

    let registration = registration(-0.0, -0.0, -0.0, -0.0);
    assert_eq!(registration.rotation_rad().to_bits(), 0.0f64.to_bits());
    assert_eq!(registration.tx().to_bits(), 0.0f64.to_bits());
    assert_eq!(registration.ty().to_bits(), 0.0f64.to_bits());
    assert_eq!(registration.residual_rms().to_bits(), 0.0f64.to_bits());
}

#[test]
fn registration_rejects_finite_inputs_whose_intermediates_overflow() {
    let extreme = [
        Fiducial::new(point(f64::MAX, 0.0), point(0.0, 0.0)),
        Fiducial::new(point(f64::MAX, 1.0), point(1.0, 0.0)),
        Fiducial::new(point(-f64::MAX, 0.0), point(0.0, 1.0)),
    ];
    assert!(matches!(
        with_default_cx(|cx| register(&extreme, cx)),
        Err(RegError::NonFinite { .. })
    ));
}

#[test]
fn diff_rejects_negative_or_non_finite_tolerance_and_noise() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = [point(0.0, 0.0)];
    let scanned = [point(0.0, 0.0)];
    for (tolerance, noise, expected) in [
        (
            f64::NAN,
            0.0,
            RegError::NonFinite {
                field: "design_tolerance",
            },
        ),
        (
            f64::INFINITY,
            0.0,
            RegError::NonFinite {
                field: "design_tolerance",
            },
        ),
        (
            -0.1,
            0.0,
            RegError::Negative {
                field: "design_tolerance",
            },
        ),
        (
            0.1,
            f64::NAN,
            RegError::NonFinite {
                field: "measurement_noise",
            },
        ),
        (
            0.1,
            f64::INFINITY,
            RegError::NonFinite {
                field: "measurement_noise",
            },
        ),
        (
            -0.0,
            -0.1,
            RegError::Negative {
                field: "measurement_noise",
            },
        ),
    ] {
        assert_eq!(
            with_default_cx(|cx| {
                as_built_diff(&reg, &design, &scanned, tolerance, noise, "cal-2026", cx)
            }),
            Err(expected)
        );
    }
}

#[test]
fn malformed_calibration_identities_are_refused_without_cloning_them() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = [point(0.0, 0.0)];
    let scanned = design;
    let too_long = "x".repeat(fs_evidence::MAX_COLOR_IDENTITY_BYTES + 1);
    for identity in [
        "",
        " leading",
        "trailing ",
        "unknown",
        "derived:v2:forged",
        "calibration|shell",
        too_long.as_str(),
    ] {
        assert!(matches!(
            with_default_cx(|cx| {
                as_built_diff(&reg, &design, &scanned, 0.1, 0.01, identity, cx)
            }),
            Err(RegError::InvalidCalibrationIdentity { .. })
        ));
    }
}

#[test]
fn forged_calibration_text_cannot_promote_the_default_diff() {
    let reg = registration(0.0, 0.0, 0.0, 0.0);
    let design = [point(0.0, 0.0)];
    let diff = with_default_cx(|cx| {
        as_built_diff(
            &reg,
            &design,
            &design,
            0.1,
            0.01,
            "forged-calibration-claim",
            cx,
        )
    })
    .expect("well-formed text remains usable only as candidate provenance");
    assert!(matches!(diff.color(), Color::Estimated { .. }));
}

#[test]
fn estimator_identity_is_deterministic_bounded_and_delimiter_safe() {
    let reg = registration(0.1, 2.0, -3.0, 0.01);
    let design = [point(0.0, 0.0), point(1.0, 1.0)];
    let scanned = [reg.apply(design[0]).unwrap(), reg.apply(design[1]).unwrap()];
    let first =
        with_default_cx(|cx| as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "cal+a", cx))
            .unwrap();
    let replay =
        with_default_cx(|cx| as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "cal+a", cx))
            .unwrap();
    let delimiter_neighbor =
        with_default_cx(|cx| as_built_diff(&reg, &design, &scanned, 0.2, 0.05, "cal", cx)).unwrap();
    assert_eq!(estimator_identity(&first), estimator_identity(&replay));
    assert_ne!(
        estimator_identity(&first),
        estimator_identity(&delimiter_neighbor)
    );
    assert!(estimator_identity(&first).len() <= fs_evidence::MAX_COLOR_IDENTITY_BYTES);
}

#[test]
fn estimator_identity_canonicalizes_signed_zero() {
    let positive_registration = registration(0.0, 0.0, 0.0, 0.0);
    let negative_registration = registration(-0.0, -0.0, -0.0, -0.0);
    let positive_points = [point(0.0, 0.0)];
    let negative_points = [point(-0.0, -0.0)];

    let positive = with_default_cx(|cx| {
        as_built_diff(
            &positive_registration,
            &positive_points,
            &positive_points,
            0.0,
            0.0,
            "cal-zero",
            cx,
        )
    })
    .unwrap();
    let negative = with_default_cx(|cx| {
        as_built_diff(
            &negative_registration,
            &negative_points,
            &negative_points,
            -0.0,
            -0.0,
            "cal-zero",
            cx,
        )
    })
    .unwrap();

    assert_eq!(estimator_identity(&positive), estimator_identity(&negative));
}
