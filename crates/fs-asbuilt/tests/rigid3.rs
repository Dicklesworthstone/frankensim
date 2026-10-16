//! Battery for 3-D rigid/similarity registration and the calibrated 6-dof
//! pose covariance. Naming: `g0_` analytic oracles, `g3_` metamorphic
//! invariances, `g4_` cancellation, `g5_` identity/determinism.

#![allow(clippy::needless_range_loop)] // Fixed 3x3/6x6 indices mirror the source modules.
#![allow(clippy::float_cmp)] // G5 bitwise-replay assertions compare exact bits on purpose.

use fs_asbuilt::rigid3::{
    Covariance3, CrossFiducialModel3, DegeneracyDiagnosis, Fiducial3, MetrologyModel3, Point3,
    Rigid3Error, estimate_calibrated_rigid3, register3, register3_similarity,
};
use fs_asbuilt::uncertainty::{HuberPolicy, OutlierDisposition};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<R>(cancelled: bool, mode: ExecMode, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = fs_exec::VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0A5B_0117,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            budget,
            mode,
        )
        .with_time_source(&clock);
        f(&cx)
    });
    let stats = pool.stats();
    assert!(
        stats.quiescent(),
        "Cx arena must be quiescent after scope: {}",
        stats.to_json()
    );
    result
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(false, ExecMode::Deterministic, Budget::INFINITE, f)
}

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).expect("finite test point")
}

/// Rodrigues rotation matrix about a unit axis.
fn rotation_matrix(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let a = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
    let (s, c) = angle.sin_cos();
    let v = 1.0 - c;
    [
        [
            c + v * a[0] * a[0],
            v * a[0] * a[1] - s * a[2],
            v * a[0] * a[2] + s * a[1],
        ],
        [
            v * a[1] * a[0] + s * a[2],
            c + v * a[1] * a[1],
            v * a[1] * a[2] - s * a[0],
        ],
        [
            v * a[2] * a[0] - s * a[1],
            v * a[2] * a[1] + s * a[0],
            c + v * a[2] * a[2],
        ],
    ]
}

fn apply_pose(rotation: &[[f64; 3]; 3], translation: [f64; 3], scale: f64, p: Point3) -> Point3 {
    let x = p.x();
    let y = p.y();
    let z = p.z();
    p3(
        scale * (rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z) + translation[0],
        scale * (rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z) + translation[1],
        scale * (rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z) + translation[2],
    )
}

fn general_design() -> Vec<Point3> {
    vec![
        p3(0.0, 0.0, 0.0),
        p3(1.0, 0.0, 0.0),
        p3(0.0, 2.0, 0.0),
        p3(0.0, 0.0, 3.0),
        p3(2.0, 1.0, 0.5),
        p3(0.5, 1.5, 2.5),
    ]
}

fn transformed_fiducials(
    design: &[Point3],
    rotation: &[[f64; 3]; 3],
    translation: [f64; 3],
    scale: f64,
) -> Vec<Fiducial3> {
    design
        .iter()
        .map(|p| Fiducial3::new(*p, apply_pose(rotation, translation, scale, *p)))
        .collect()
}

/// Deterministic test-local RNG: xorshift64* with Box-Muller normals. Test
/// infrastructure only; production randomness lives in fs-rand.
struct TestRng(u64);

impl TestRng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn uniform(&mut self) -> f64 {
        // Strictly inside (0, 1) so the Box-Muller logarithm is finite.
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 + f64::EPSILON
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (core::f64::consts::TAU * u2).cos()
    }
}

const TRUTH_AXIS: [f64; 3] = [1.0, 2.0, 3.0];
const TRUTH_ANGLE: f64 = 0.7;
const TRUTH_TRANSLATION: [f64; 3] = [0.4, -1.1, 2.2];

#[test]
fn g0_exact_recovery_of_a_general_rotation() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0);
    let registration = with_default_cx(|cx| register3(&fiducials, cx)).expect("clean recovery");
    for row in 0..3 {
        for column in 0..3 {
            assert!(
                (registration.rotation()[row][column] - rotation[row][column]).abs() < 1e-10,
                "rotation[{row}][{column}] drifted"
            );
        }
    }
    for axis in 0..3 {
        assert!((registration.translation()[axis] - TRUTH_TRANSLATION[axis]).abs() < 1e-9);
    }
    assert!(registration.residual_rms() < 1e-9);
    assert!((registration.rotation_angle_rad() - TRUTH_ANGLE).abs() < 1e-10);
    let condition = registration.condition();
    assert!(!condition.coplanar_design());
    assert!(!condition.reflection_preferred());
    let mapped = registration
        .apply(p3(0.25, 0.5, 0.75))
        .expect("finite apply");
    let expected = apply_pose(&rotation, TRUTH_TRANSLATION, 1.0, p3(0.25, 0.5, 0.75));
    assert!((mapped.x() - expected.x()).abs() < 1e-9);
    assert!((mapped.y() - expected.y()).abs() < 1e-9);
    assert!((mapped.z() - expected.z()).abs() < 1e-9);
}

#[test]
fn g0_coplanar_design_recovers_exactly_and_is_flagged() {
    let design = vec![
        p3(0.0, 0.0, 0.0),
        p3(4.0, 0.0, 0.0),
        p3(4.0, 3.0, 0.0),
        p3(0.0, 3.0, 0.0),
        p3(1.0, 2.0, 0.0),
    ];
    let rotation = rotation_matrix([0.0, 0.0, 1.0], 0.5);
    let fiducials = transformed_fiducials(&design, &rotation, [1.0, -2.0, 0.25], 1.0);
    let registration = with_default_cx(|cx| register3(&fiducials, cx)).expect("coplanar admitted");
    assert!(registration.condition().coplanar_design());
    assert!(registration.condition().coplanar_cross());
    assert!(!registration.condition().reflection_preferred());
    for row in 0..3 {
        for column in 0..3 {
            assert!((registration.rotation()[row][column] - rotation[row][column]).abs() < 1e-10);
        }
    }
    assert!(registration.residual_rms() < 1e-9);
}

#[test]
fn degenerate_configurations_refuse_with_geometric_diagnosis() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);

    let two = &transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0)[..2];
    assert!(matches!(
        with_default_cx(|cx| register3(two, cx)),
        Err(Rigid3Error::TooFewFiducials { have: 2, need: 3 })
    ));

    let coincident: Vec<Fiducial3> = (0..4)
        .map(|i| Fiducial3::new(p3(1.0, 1.0, 1.0), p3(f64::from(i), 0.0, 0.0)))
        .collect();
    assert!(matches!(
        with_default_cx(|cx| register3(&coincident, cx)),
        Err(Rigid3Error::DegenerateDesign {
            diagnosis: DegeneracyDiagnosis::Coincident
        })
    ));

    let collinear_design: Vec<Fiducial3> = (0..4)
        .map(|i| {
            let p = p3(f64::from(i), 0.0, 0.0);
            Fiducial3::new(p, apply_pose(&rotation, TRUTH_TRANSLATION, 1.0, p))
        })
        .collect();
    assert!(matches!(
        with_default_cx(|cx| register3(&collinear_design, cx)),
        Err(Rigid3Error::DegenerateDesign {
            diagnosis: DegeneracyDiagnosis::Collinear
        })
    ));

    // Full-rank design, measured collapsed onto a line.
    let collinear_measured: Vec<Fiducial3> = general_design()
        .iter()
        .map(|p| {
            let along = p.x() + p.y() + p.z();
            Fiducial3::new(*p, p3(along, 0.0, 0.0))
        })
        .collect();
    assert!(matches!(
        with_default_cx(|cx| register3(&collinear_measured, cx)),
        Err(Rigid3Error::DegenerateMeasured {
            diagnosis: DegeneracyDiagnosis::Collinear
        })
    ));
}

#[test]
fn mirrored_generic_data_reports_reflection_preference() {
    let fiducials: Vec<Fiducial3> = general_design()
        .iter()
        .map(|p| Fiducial3::new(*p, p3(-p.x(), p.y(), p.z())))
        .collect();
    let registration =
        with_default_cx(|cx| register3(&fiducials, cx)).expect("mirrored generic data still fits");
    assert!(registration.condition().reflection_preferred());
    assert!(registration.residual_rms() > 0.1);
}

#[test]
fn mirrored_symmetric_data_refuses_as_ambiguous() {
    // Centered unit-cube corners have an isotropic scatter, so the mirrored
    // cross-covariance has coincident trailing singular values.
    let mut fiducials = Vec::new();
    for sx in [-0.5f64, 0.5] {
        for sy in [-0.5f64, 0.5] {
            for sz in [-0.5f64, 0.5] {
                fiducials.push(Fiducial3::new(p3(sx, sy, sz), p3(-sx, sy, sz)));
            }
        }
    }
    assert!(matches!(
        with_default_cx(|cx| register3(&fiducials, cx)),
        Err(Rigid3Error::AmbiguousRotation)
    ));
}

#[test]
fn g3_rigid_pretransform_of_both_clouds_preserves_fit_quality() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let mut rng = TestRng(0x5EED_0001);
    let noisy: Vec<Fiducial3> = general_design()
        .iter()
        .map(|p| {
            let q = apply_pose(&rotation, TRUTH_TRANSLATION, 1.0, *p);
            Fiducial3::new(
                *p,
                p3(
                    q.x() + 1e-3 * rng.normal(),
                    q.y() + 1e-3 * rng.normal(),
                    q.z() + 1e-3 * rng.normal(),
                ),
            )
        })
        .collect();
    let baseline = with_default_cx(|cx| register3(&noisy, cx)).expect("baseline");

    let g = rotation_matrix([-2.0, 1.0, 0.4], 1.1);
    let g_translation = [3.0, 4.0, -5.0];
    let moved: Vec<Fiducial3> = noisy
        .iter()
        .map(|f| {
            Fiducial3::new(
                apply_pose(&g, g_translation, 1.0, f.design()),
                apply_pose(&g, g_translation, 1.0, f.measured()),
            )
        })
        .collect();
    let conjugated = with_default_cx(|cx| register3(&moved, cx)).expect("conjugated");

    assert!((baseline.residual_rms() - conjugated.residual_rms()).abs() < 1e-9);
    // Conjugation by a rigid motion preserves the rotation angle.
    assert!((baseline.rotation_angle_rad() - conjugated.rotation_angle_rad()).abs() < 1e-9);
}

#[test]
fn similarity_reports_seeded_unit_error_with_nearest_conversion() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 25.4);
    let similarity = with_default_cx(|cx| register3_similarity(&fiducials, 0.01, cx))
        .expect("similarity recovery");
    let scale = similarity.scale();
    assert!((scale.estimate() - 25.4).abs() < 1e-8 * 25.4);
    assert!(scale.standard_error() < 1e-8);
    let suspicion = scale.suspicion().expect("scale 25.4 must raise suspicion");
    assert_eq!(suspicion.conversion_name(), "25.4 (inch to millimetre)");
    assert!(suspicion.log_gap() < 1e-8);
    for row in 0..3 {
        for column in 0..3 {
            assert!((similarity.rotation()[row][column] - rotation[row][column]).abs() < 1e-9);
        }
    }
}

#[test]
fn similarity_near_unity_scale_raises_no_suspicion() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0004);
    let similarity = with_default_cx(|cx| register3_similarity(&fiducials, 0.01, cx))
        .expect("similarity recovery");
    assert!((similarity.scale().estimate() - 1.0004).abs() < 1e-8);
    assert!(similarity.scale().suspicion().is_none());
}

#[test]
fn similarity_rejects_invalid_tolerance() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0);
    for bad in [f64::NAN, f64::INFINITY, -0.5] {
        assert!(matches!(
            with_default_cx(|cx| register3_similarity(&fiducials, bad, cx)),
            Err(Rigid3Error::InvalidScalar {
                field: "scale_tolerance",
                ..
            })
        ));
    }
}

fn isotropic_model(count: usize, variance: f64, huber: HuberPolicy) -> MetrologyModel3 {
    let covariance =
        Covariance3::new(variance, 0.0, 0.0, variance, 0.0, variance).expect("iso covariance");
    MetrologyModel3::new(
        vec![covariance; count],
        CrossFiducialModel3::Independent,
        huber,
        "sha256:1111111111111111111111111111111111111111111111111111111111111111",
    )
    .expect("valid model")
}

#[test]
fn g0_calibrated_covariance_matches_the_analytic_axis_configuration() {
    // Centered +/- axis points make every block diagonal and hand-checkable:
    // cov(t) = (sigma^2 / n) I, cov(r) = sigma^2 diag(1/(2(b^2+c^2)), ...).
    let (a, b, c) = (1.0f64, 2.0f64, 3.0f64);
    let design = [
        p3(a, 0.0, 0.0),
        p3(-a, 0.0, 0.0),
        p3(0.0, b, 0.0),
        p3(0.0, -b, 0.0),
        p3(0.0, 0.0, c),
        p3(0.0, 0.0, -c),
    ];
    let fiducials: Vec<Fiducial3> = design.iter().map(|p| Fiducial3::new(*p, *p)).collect();
    let variance = 0.04;
    let model = isotropic_model(fiducials.len(), variance, HuberPolicy::Disabled);
    let calibrated = with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &model, cx))
        .expect("calibrated fit");

    assert_eq!(calibrated.degrees_of_freedom(), 12);
    assert!(!calibrated.robust_conditional());
    let covariance = calibrated.covariance();
    let expected_translation = variance / 6.0;
    let expected_rotation = [
        variance / (2.0 * (b * b + c * c)),
        variance / (2.0 * (a * a + c * c)),
        variance / (2.0 * (a * a + b * b)),
    ];
    for axis in 0..3 {
        assert!(
            (covariance[axis][axis] - expected_translation).abs() < 1e-12,
            "translation variance axis {axis}"
        );
        assert!(
            (covariance[3 + axis][3 + axis] - expected_rotation[axis]).abs()
                < 1e-12 * expected_rotation[axis].max(1.0),
            "rotation variance axis {axis}"
        );
    }
    for row in 0..6 {
        for column in 0..6 {
            if row != column {
                assert!(
                    covariance[row][column].abs() < 1e-12,
                    "off-diagonal [{row}][{column}] should vanish for this configuration"
                );
            }
        }
    }
    let leverage_sum: f64 = calibrated.leverage().iter().sum();
    assert!((leverage_sum - 6.0).abs() < 1e-9);
    for diagnostic in calibrated.outlier_diagnostics() {
        assert_eq!(diagnostic.disposition(), OutlierDisposition::NotEvaluated);
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_calibrated_covariance_agrees_with_monte_carlo() {
    let rotation = rotation_matrix([0.3, -1.0, 0.2], 0.4);
    let translation = [1.0, 2.0, -0.5];
    let design = [
        p3(0.0, 0.0, 0.0),
        p3(2.0, 0.0, 0.0),
        p3(2.0, 3.0, 0.0),
        p3(0.0, 3.0, 0.0),
        p3(0.0, 0.0, 1.0),
        p3(2.0, 0.0, 1.0),
        p3(2.0, 3.0, 1.0),
        p3(0.5, 1.0, 0.6),
    ];
    let sigma = 0.01f64;
    let model = isotropic_model(design.len(), sigma * sigma, HuberPolicy::Disabled);

    // The reported covariance is for (centroid image, rotation vector); the
    // design centroid is the pivot for equal weights.
    let mut centroid = [0.0f64; 3];
    for p in &design {
        centroid[0] += p.x() / design.len() as f64;
        centroid[1] += p.y() / design.len() as f64;
        centroid[2] += p.z() / design.len() as f64;
    }

    let replicates = 1500usize;
    let mut rng = TestRng(0xDEAD_BEEF_0042);
    let mut samples: Vec<[f64; 6]> = Vec::with_capacity(replicates);
    let mut reported = [[0.0f64; 6]; 6];
    for _ in 0..replicates {
        let fiducials: Vec<Fiducial3> = design
            .iter()
            .map(|p| {
                let q = apply_pose(&rotation, translation, 1.0, *p);
                Fiducial3::new(
                    *p,
                    p3(
                        q.x() + sigma * rng.normal(),
                        q.y() + sigma * rng.normal(),
                        q.z() + sigma * rng.normal(),
                    ),
                )
            })
            .collect();
        let calibrated = with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &model, cx))
            .expect("replicate fit");
        let fit = calibrated.registration();
        // Small-angle rotation-vector delta of R_hat * R*^T.
        let mut relative = [[0.0f64; 3]; 3];
        for row in 0..3 {
            for column in 0..3 {
                let mut value = 0.0;
                for inner in 0..3 {
                    value += fit.rotation()[row][inner] * rotation[column][inner];
                }
                relative[row][column] = value;
            }
        }
        let omega = [
            (relative[2][1] - relative[1][2]) / 2.0,
            (relative[0][2] - relative[2][0]) / 2.0,
            (relative[1][0] - relative[0][1]) / 2.0,
        ];
        // Predicted centroid image under the fitted and true poses.
        let fitted_centroid = {
            let r = fit.rotation();
            [
                r[0][0] * centroid[0]
                    + r[0][1] * centroid[1]
                    + r[0][2] * centroid[2]
                    + fit.translation()[0],
                r[1][0] * centroid[0]
                    + r[1][1] * centroid[1]
                    + r[1][2] * centroid[2]
                    + fit.translation()[1],
                r[2][0] * centroid[0]
                    + r[2][1] * centroid[1]
                    + r[2][2] * centroid[2]
                    + fit.translation()[2],
            ]
        };
        let true_centroid = {
            let r = &rotation;
            [
                r[0][0] * centroid[0]
                    + r[0][1] * centroid[1]
                    + r[0][2] * centroid[2]
                    + translation[0],
                r[1][0] * centroid[0]
                    + r[1][1] * centroid[1]
                    + r[1][2] * centroid[2]
                    + translation[1],
                r[2][0] * centroid[0]
                    + r[2][1] * centroid[1]
                    + r[2][2] * centroid[2]
                    + translation[2],
            ]
        };
        samples.push([
            fitted_centroid[0] - true_centroid[0],
            fitted_centroid[1] - true_centroid[1],
            fitted_centroid[2] - true_centroid[2],
            omega[0],
            omega[1],
            omega[2],
        ]);
        for row in 0..6 {
            for column in 0..6 {
                reported[row][column] += calibrated.covariance()[row][column] / replicates as f64;
            }
        }
    }

    let mut mean = [0.0f64; 6];
    for sample in &samples {
        for axis in 0..6 {
            mean[axis] += sample[axis] / replicates as f64;
        }
    }
    for axis in 0..6 {
        let mut empirical = 0.0f64;
        for sample in &samples {
            let delta = sample[axis] - mean[axis];
            empirical += delta * delta / (replicates as f64 - 1.0);
        }
        let predicted = reported[axis][axis];
        assert!(
            (empirical - predicted).abs() < 0.2 * predicted,
            "axis {axis}: empirical {empirical:.3e} vs predicted {predicted:.3e}"
        );
    }
}

#[test]
fn huber_downweights_a_seeded_outlier_and_improves_the_fit() {
    let rotation = rotation_matrix([0.0, 0.0, 1.0], 0.3);
    let translation = [0.5, -0.25, 1.5];
    let design = [
        p3(0.0, 0.0, 0.0),
        p3(2.0, 0.0, 0.0),
        p3(2.0, 3.0, 0.0),
        p3(0.0, 3.0, 0.0),
        p3(0.0, 0.0, 1.0),
        p3(2.0, 3.0, 1.0),
        p3(1.0, 1.5, 0.5),
        p3(0.5, 2.0, 0.8),
    ];
    let sigma = 0.01f64;
    let mut rng = TestRng(0x0BAD_F00D_7777);
    let mut fiducials: Vec<Fiducial3> = design
        .iter()
        .map(|p| {
            let q = apply_pose(&rotation, translation, 1.0, *p);
            Fiducial3::new(
                *p,
                p3(
                    q.x() + sigma * rng.normal(),
                    q.y() + sigma * rng.normal(),
                    q.z() + sigma * rng.normal(),
                ),
            )
        })
        .collect();
    // Corrupt one measurement far beyond the noise floor.
    let corrupted = fiducials[3].measured();
    fiducials[3] = Fiducial3::new(
        fiducials[3].design(),
        p3(corrupted.x() + 0.5, corrupted.y(), corrupted.z()),
    );

    let plain_model = isotropic_model(fiducials.len(), sigma * sigma, HuberPolicy::Disabled);
    let robust_model = isotropic_model(
        fiducials.len(),
        sigma * sigma,
        HuberPolicy::new(2.0, 8).expect("huber policy"),
    );
    let (plain, robust) = with_default_cx(|cx| {
        (
            estimate_calibrated_rigid3(&fiducials, &plain_model, cx).expect("plain fit"),
            estimate_calibrated_rigid3(&fiducials, &robust_model, cx).expect("robust fit"),
        )
    });

    assert!(robust.robust_conditional());
    let outlier = &robust.outlier_diagnostics()[3];
    assert_eq!(outlier.disposition(), OutlierDisposition::Downweighted);
    assert!(outlier.robust_weight() < 0.2);
    // Good points may brush the 2-sigma Huber threshold (a chi-3 norm
    // exceeds 2 about a quarter of the time), so the locked property is
    // SEPARATION: every clean point keeps at least five times the outlier's
    // weight, and the outlier carries the strict minimum.
    for (index, diagnostic) in robust.outlier_diagnostics().iter().enumerate() {
        if index != 3 {
            assert!(
                diagnostic.robust_weight() > 5.0 * outlier.robust_weight(),
                "clean point {index} weight {:.3} not separated from outlier {:.3}",
                diagnostic.robust_weight(),
                outlier.robust_weight()
            );
        }
    }

    let angle_error = |fit: &fs_asbuilt::rigid3::CalibratedRigid3Registration| {
        let mut relative = [[0.0f64; 3]; 3];
        let r = fit.registration().rotation();
        for row in 0..3 {
            for column in 0..3 {
                let mut value = 0.0;
                for inner in 0..3 {
                    value += r[row][inner] * rotation[column][inner];
                }
                relative[row][column] = value;
            }
        }
        let trace = relative[0][0] + relative[1][1] + relative[2][2];
        (((trace - 1.0) / 2.0).clamp(-1.0, 1.0)).acos()
    };
    assert!(
        angle_error(&robust) < angle_error(&plain),
        "robust fit should sit closer to the seeded truth"
    );
}

#[test]
fn calibrated_model_refusals_are_typed() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0);

    let covariance = Covariance3::new(1e-4, 0.0, 0.0, 1e-4, 0.0, 1e-4).expect("iso");
    let unknown = MetrologyModel3::new(
        vec![covariance; fiducials.len()],
        CrossFiducialModel3::Unknown,
        HuberPolicy::Disabled,
        "sha256:1111111111111111111111111111111111111111111111111111111111111111",
    )
    .expect("constructible");
    assert!(matches!(
        with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &unknown, cx)),
        Err(Rigid3Error::UnknownDependence)
    ));

    let short = MetrologyModel3::new(
        vec![covariance; fiducials.len() - 1],
        CrossFiducialModel3::Independent,
        HuberPolicy::Disabled,
        "sha256:1111111111111111111111111111111111111111111111111111111111111111",
    )
    .expect("constructible");
    assert!(matches!(
        with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &short, cx)),
        Err(Rigid3Error::LengthMismatch { .. })
    ));

    assert!(matches!(
        Covariance3::new(1.0, 2.0, 0.0, 1.0, 0.0, 1.0),
        Err(Rigid3Error::NonPositiveDefiniteCovariance { .. })
    ));
    assert!(matches!(
        Covariance3::new(f64::NAN, 0.0, 0.0, 1.0, 0.0, 1.0),
        Err(Rigid3Error::InvalidScalar { .. })
    ));
    assert!(matches!(
        Covariance3::new(-1.0, 0.0, 0.0, 1.0, 0.0, 1.0),
        Err(Rigid3Error::NonPositiveDefiniteCovariance { .. })
    ));

    assert!(matches!(
        MetrologyModel3::new(
            vec![covariance; 2],
            CrossFiducialModel3::Independent,
            HuberPolicy::Disabled,
            "sha256:1111111111111111111111111111111111111111111111111111111111111111",
        ),
        Err(Rigid3Error::TooFewFiducials { .. })
    ));
    assert!(matches!(
        MetrologyModel3::new(
            vec![covariance; 4],
            CrossFiducialModel3::Independent,
            HuberPolicy::Enabled {
                threshold: -1.0,
                iterations: 4
            },
            "sha256:1111111111111111111111111111111111111111111111111111111111111111",
        ),
        Err(Rigid3Error::InvalidScalar {
            field: "huber.threshold",
            ..
        })
    ));
    assert!(matches!(
        MetrologyModel3::new(
            vec![covariance; 4],
            CrossFiducialModel3::Independent,
            HuberPolicy::Disabled,
            "not a valid identity !!",
        ),
        Err(Rigid3Error::InvalidCalibrationIdentity { .. })
    ));
}

#[test]
fn oversized_input_refuses_before_scanning() {
    let fiducial = Fiducial3::new(p3(0.0, 0.0, 0.0), p3(1.0, 0.0, 0.0));
    let oversized = vec![fiducial; 1_000_001];
    assert!(matches!(
        with_default_cx(|cx| register3(&oversized, cx)),
        Err(Rigid3Error::TooManyPoints {
            have: 1_000_001,
            max: 1_000_000
        })
    ));
}

#[test]
fn g4_cancellation_is_structured_and_publishes_nothing() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let fiducials = transformed_fiducials(&general_design(), &rotation, TRUTH_TRANSLATION, 1.0);
    assert!(matches!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            register3(&fiducials, cx)
        }),
        Err(Rigid3Error::Cancelled {
            phase: "rigid3.design-extent"
        })
    ));
    let model = isotropic_model(fiducials.len(), 1e-4, HuberPolicy::Disabled);
    assert!(matches!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            estimate_calibrated_rigid3(&fiducials, &model, cx)
        }),
        Err(Rigid3Error::Cancelled {
            phase: "rigid3.whitening"
        })
    ));
}

#[test]
fn g5_replay_is_bitwise_and_identity_moves_with_inputs() {
    let rotation = rotation_matrix(TRUTH_AXIS, TRUTH_ANGLE);
    let mut rng = TestRng(0x1DEA_0002);
    let fiducials: Vec<Fiducial3> = general_design()
        .iter()
        .map(|p| {
            let q = apply_pose(&rotation, TRUTH_TRANSLATION, 1.0, *p);
            Fiducial3::new(
                *p,
                p3(
                    q.x() + 1e-3 * rng.normal(),
                    q.y() + 1e-3 * rng.normal(),
                    q.z() + 1e-3 * rng.normal(),
                ),
            )
        })
        .collect();
    let model = isotropic_model(
        fiducials.len(),
        1e-6,
        HuberPolicy::new(2.5, 4).expect("huber"),
    );

    let (first, second) = with_default_cx(|cx| {
        (
            estimate_calibrated_rigid3(&fiducials, &model, cx).expect("first"),
            estimate_calibrated_rigid3(&fiducials, &model, cx).expect("second"),
        )
    });
    assert_eq!(first.model_identity(), second.model_identity());
    assert_eq!(
        first.registration().rotation(),
        second.registration().rotation()
    );
    assert_eq!(first.covariance(), second.covariance());
    assert_eq!(first.robust_weights(), second.robust_weights());

    let mut perturbed = fiducials.clone();
    let moved = perturbed[0].measured();
    perturbed[0] = Fiducial3::new(
        perturbed[0].design(),
        p3(moved.x() + 1e-9, moved.y(), moved.z()),
    );
    let third = with_default_cx(|cx| estimate_calibrated_rigid3(&perturbed, &model, cx))
        .expect("perturbed fit");
    assert_ne!(first.model_identity(), third.model_identity());

    let renamed = MetrologyModel3::new(
        model.fiducial_covariances().to_vec(),
        CrossFiducialModel3::Independent,
        model.huber(),
        "sha256:2222222222222222222222222222222222222222222222222222222222222222",
    )
    .expect("renamed model");
    let fourth = with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &renamed, cx))
        .expect("renamed fit");
    assert_ne!(first.model_identity(), fourth.model_identity());
    // The calibration identity is bound into the identity but must not move
    // the numbers.
    assert_eq!(first.covariance(), fourth.covariance());
}
