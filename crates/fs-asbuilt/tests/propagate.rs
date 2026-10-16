//! Battery for pose-covariance → geometry-budget propagation. Naming: `g0_`
//! analytic oracles, `g4_` cancellation, `g5_` identity/determinism. The e2e
//! case drives register → calibrate → propagate → full 8-term budget and
//! logs the per-stage table with the correlation structure.

#![allow(clippy::needless_range_loop)] // Fixed 6-dof indices mirror the source module.
#![allow(clippy::float_cmp)] // G5 determinism assertions compare exact bits on purpose.

use fs_asbuilt::propagate::{
    CoveragePolicy, GEOMETRY_PROPAGATION_SCHEMA_VERSION, PropagateError, PropagationMethod,
    QoiEvaluator, QoiSensitivity, propagate_pose_covariance,
};
use fs_asbuilt::rigid3::{
    Covariance3, CrossFiducialModel3, Fiducial3, MetrologyModel3, Point3,
    estimate_calibrated_rigid3,
};
use fs_asbuilt::uncertainty::HuberPolicy;
use fs_evidence::uncertainty::{
    EngineeringUncertaintyBudget, EngineeringUncertaintyKind, EngineeringUncertaintyTerm,
    TermValue, UncertaintyArtifactRef,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use std::fmt::Write as _;

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
                seed: 0x0A5B_0119,
                kernel_id: 5,
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

/// The hand-checkable axis configuration from the rigid3 battery: centered
/// +/- axis points at identity pose with isotropic noise, so the pose
/// covariance is diagonal and known in closed form.
fn axis_fixture(variance: f64) -> fs_asbuilt::rigid3::CalibratedRigid3Registration {
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
    let covariance =
        Covariance3::new(variance, 0.0, 0.0, variance, 0.0, variance).expect("iso covariance");
    let model = MetrologyModel3::new(
        vec![covariance; fiducials.len()],
        CrossFiducialModel3::Independent,
        HuberPolicy::Disabled,
        "sha256:1111111111111111111111111111111111111111111111111111111111111111",
    )
    .expect("model");
    with_default_cx(|cx| estimate_calibrated_rigid3(&fiducials, &model, cx)).expect("fit")
}

fn coverage() -> CoveragePolicy {
    CoveragePolicy::new(0.95, 2.0).expect("coverage policy")
}

struct LinearEvaluator {
    gradients: Vec<[f64; 6]>,
}

impl QoiEvaluator for LinearEvaluator {
    fn evaluate(&self, pose_delta: &[f64; 6]) -> Result<Vec<f64>, String> {
        Ok(self
            .gradients
            .iter()
            .map(|gradient| {
                let mut value = 0.0;
                for axis in 0..6 {
                    value += gradient[axis] * pose_delta[axis];
                }
                value
            })
            .collect())
    }
}

struct QuadraticEvaluator {
    gradients: Vec<[f64; 6]>,
    curvature: f64,
}

impl QoiEvaluator for QuadraticEvaluator {
    fn evaluate(&self, pose_delta: &[f64; 6]) -> Result<Vec<f64>, String> {
        let mut norm_squared = 0.0;
        for axis in 0..6 {
            norm_squared += pose_delta[axis] * pose_delta[axis];
        }
        Ok(self
            .gradients
            .iter()
            .map(|gradient| {
                let mut value = self.curvature * norm_squared;
                for axis in 0..6 {
                    value += gradient[axis] * pose_delta[axis];
                }
                value
            })
            .collect())
    }
}

#[test]
fn g0_translation_gradient_recovers_the_analytic_marginal() {
    let variance = 0.04f64;
    let calibrated = axis_fixture(variance);
    let qoi = QoiSensitivity::new("gap-z", "millimetre", [0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        .expect("sensitivity");
    let propagation = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &[qoi],
            coverage(),
            None,
            0.0,
            "no QoI evaluator wired in this unit fixture",
            cx,
        )
    })
    .expect("propagation");
    // cov(t) = (variance / 6) I for this fixture.
    let expected = (variance / 6.0).sqrt();
    let deviation = propagation.standard_deviations()[0];
    assert!(
        (deviation - expected).abs() < 1e-12,
        "expected {expected}, saw {deviation}"
    );
    assert!(matches!(
        propagation.method(),
        PropagationMethod::LinearizedUnchecked { .. }
    ));
    let term = propagation.geometry_term(0).expect("term");
    assert_eq!(term.kind(), EngineeringUncertaintyKind::Geometry);
    match term.value() {
        TermValue::Distribution(distribution) => {
            assert_eq!(distribution.mean, 0.0);
            assert!((distribution.standard_deviation - expected).abs() < 1e-12);
            assert!((distribution.conservative_half_width - 2.0 * expected).abs() < 1e-12);
            assert_eq!(distribution.level, 0.95);
            assert_eq!(distribution.replay.digest(), propagation.record_identity());
        }
        other => panic!("expected a distribution term, got {other:?}"),
    }
}

#[test]
fn g0_shared_pose_error_produces_exact_correlation() {
    let calibrated = axis_fixture(0.01);
    let base = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let doubled = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let negated = [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let qois = vec![
        QoiSensitivity::new("qoi-a", "millimetre", base).expect("a"),
        QoiSensitivity::new("qoi-b", "millimetre", doubled).expect("b"),
        QoiSensitivity::new("qoi-c", "millimetre", negated).expect("c"),
    ];
    let propagation = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            None,
            0.0,
            "correlation oracle fixture",
            cx,
        )
    })
    .expect("propagation");
    // One pose error moves every QoI: parallel gradients correlate exactly.
    assert!((propagation.correlation(0, 1).expect("rho ab") - 1.0).abs() < 1e-12);
    assert!((propagation.correlation(0, 2).expect("rho ac") + 1.0).abs() < 1e-12);
    // Covariance scales bilinearly with the gradients.
    let var_a = propagation.cross_covariance()[0];
    let cov_ab = propagation.cross_covariance()[1];
    assert!((cov_ab - 2.0 * var_a).abs() < 1e-12 * var_a.max(1.0));
}

#[test]
fn spot_check_accepts_a_linear_map_and_records_the_method() {
    let calibrated = axis_fixture(0.01);
    let gradients = vec![
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.5, -0.25, 0.0, 0.0, 0.0, 1.5],
    ];
    let qois = vec![
        QoiSensitivity::new("gap-z", "millimetre", gradients[0]).expect("a"),
        QoiSensitivity::new("tilt-mix", "millimetre", gradients[1]).expect("b"),
    ];
    let evaluator = LinearEvaluator {
        gradients: gradients.clone(),
    };
    let propagation = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            Some(&evaluator),
            1e-9,
            "",
            cx,
        )
    })
    .expect("propagation");
    match propagation.method() {
        PropagationMethod::LinearizedSpotChecked {
            samples,
            max_relative_gap,
            tolerance,
        } => {
            assert_eq!(*samples, 12);
            assert!(*max_relative_gap < 1e-10);
            assert_eq!(*tolerance, 1e-9);
        }
        other => panic!("expected a spot-checked method, got {other:?}"),
    }
    assert!(matches!(
        propagation.geometry_term(1).expect("term").value(),
        TermValue::Distribution(_)
    ));
}

#[test]
fn spot_check_rejection_downgrades_terms_to_unknown() {
    let calibrated = axis_fixture(0.01);
    let gradients = vec![[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]];
    let qois = vec![QoiSensitivity::new("gap-z", "millimetre", gradients[0]).expect("a")];
    // Curvature strong enough that the one-sigma probes leave the linear
    // regime by far more than the declared tolerance.
    let evaluator = QuadraticEvaluator {
        gradients,
        curvature: 500.0,
    };
    let propagation = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            Some(&evaluator),
            0.05,
            "",
            cx,
        )
    })
    .expect("propagation still publishes a record");
    match propagation.method() {
        PropagationMethod::LinearizationRejected {
            max_relative_gap,
            tolerance,
            worst_qoi,
            ..
        } => {
            assert!(*max_relative_gap > *tolerance);
            assert_eq!(*worst_qoi, 0);
        }
        other => panic!("expected rejection, got {other:?}"),
    }
    let term = propagation.geometry_term(0).expect("term");
    match term.value() {
        TermValue::Unknown { reason } => {
            assert!(reason.contains("linearization"), "reason: {reason}");
            assert!(
                reason.contains("gap-z"),
                "reason names the worst QoI: {reason}"
            );
        }
        other => panic!("expected an unknown term after rejection, got {other:?}"),
    }
}

/// Deterministic test-local RNG (xorshift64* + Box-Muller).
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
        (self.next_u64() >> 11) as f64 / 9_007_199_254_740_992.0 + f64::EPSILON
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (core::f64::consts::TAU * u2).cos()
    }
}

#[test]
fn g0_synthetic_truth_variance_is_recovered_by_sampling() {
    let calibrated = axis_fixture(0.0025);
    let gradient = [0.7, -0.2, 1.1, 0.4, -0.9, 0.3];
    let qois = vec![QoiSensitivity::new("mixed", "millimetre", gradient).expect("q")];
    let propagation = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            None,
            0.0,
            "synthetic-truth sampling oracle",
            cx,
        )
    })
    .expect("propagation");
    let predicted_variance = propagation.cross_covariance()[0];

    // Sample pose errors from the reported covariance via its Cholesky
    // factor and push them through the linear QoI: the empirical variance
    // must match the propagated one.
    let covariance = propagation.pose_covariance();
    let mut lower = [[0.0f64; 6]; 6];
    for row in 0..6 {
        for column in 0..=row {
            let mut value = covariance[row][column];
            for prior in 0..column {
                value -= lower[row][prior] * lower[column][prior];
            }
            if row == column {
                lower[row][column] = value.sqrt();
            } else {
                lower[row][column] = value / lower[column][column];
            }
        }
    }
    let mut rng = TestRng(0x5EED_0003);
    let replicates = 4000usize;
    let mut sum = 0.0f64;
    let mut sum_squares = 0.0f64;
    for _ in 0..replicates {
        let mut standard = [0.0f64; 6];
        for axis in 0..6 {
            standard[axis] = rng.normal();
        }
        let mut delta = [0.0f64; 6];
        for row in 0..6 {
            for column in 0..=row {
                delta[row] += lower[row][column] * standard[column];
            }
        }
        let mut qoi_delta = 0.0f64;
        for axis in 0..6 {
            qoi_delta += gradient[axis] * delta[axis];
        }
        sum += qoi_delta;
        sum_squares += qoi_delta * qoi_delta;
    }
    let mean = sum / replicates as f64;
    let empirical_variance = sum_squares / replicates as f64 - mean * mean;
    assert!(
        (empirical_variance - predicted_variance).abs() < 0.1 * predicted_variance,
        "empirical {empirical_variance:.3e} vs predicted {predicted_variance:.3e}"
    );
}

#[test]
fn declaration_refusals_are_typed() {
    let calibrated = axis_fixture(0.01);
    let qoi = QoiSensitivity::new("gap-z", "millimetre", [0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
        .expect("sensitivity");

    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            &[],
            coverage(),
            None,
            0.0,
            "reason",
            cx
        )),
        Err(PropagateError::EmptyQois)
    ));

    let duplicated = vec![qoi.clone(), qoi.clone()];
    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            &duplicated,
            coverage(),
            None,
            0.0,
            "reason",
            cx
        )),
        Err(PropagateError::DuplicateQoi { index: 1 })
    ));

    let oversized = vec![qoi.clone(); 65];
    // Names are identical, but the cap check fires first by design order?
    // No: the cap check runs before the duplicate scan, so this exercises
    // the cap refusal deterministically.
    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            &oversized,
            coverage(),
            None,
            0.0,
            "reason",
            cx
        )),
        Err(PropagateError::TooManyQois { have: 65, max: 64 })
    ));

    let evaluator = LinearEvaluator {
        gradients: vec![[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]],
    };
    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            core::slice::from_ref(&qoi),
            coverage(),
            Some(&evaluator),
            f64::NAN,
            "",
            cx
        )),
        Err(PropagateError::InvalidScalar {
            field: "spot_check_tolerance",
            ..
        })
    ));
    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            core::slice::from_ref(&qoi),
            coverage(),
            None,
            0.0,
            "   ",
            cx
        )),
        Err(PropagateError::InvalidScalar {
            field: "unchecked_reason",
            ..
        })
    ));

    let short_evaluator = LinearEvaluator { gradients: vec![] };
    assert!(matches!(
        with_default_cx(|cx| propagate_pose_covariance(
            &calibrated,
            core::slice::from_ref(&qoi),
            coverage(),
            Some(&short_evaluator),
            0.1,
            "",
            cx
        )),
        Err(PropagateError::EvaluatorLengthMismatch {
            expected: 1,
            found: 0
        })
    ));

    assert!(matches!(
        QoiSensitivity::new("bad name with spaces", "millimetre", [0.0; 6]),
        Err(PropagateError::InvalidName { .. })
    ));
    assert!(matches!(
        QoiSensitivity::new("ok", "millimetre", [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0]),
        Err(PropagateError::InvalidScalar { .. })
    ));
    assert!(CoveragePolicy::new(1.5, 2.0).is_err());
    assert!(CoveragePolicy::new(0.95, 0.0).is_err());
}

#[test]
fn g5_record_identity_is_deterministic_and_input_sensitive() {
    let calibrated = axis_fixture(0.01);
    let qois = vec![
        QoiSensitivity::new("gap-z", "millimetre", [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]).expect("a"),
    ];
    let (first, second) = with_default_cx(|cx| {
        (
            propagate_pose_covariance(&calibrated, &qois, coverage(), None, 0.0, "reason", cx)
                .expect("first"),
            propagate_pose_covariance(&calibrated, &qois, coverage(), None, 0.0, "reason", cx)
                .expect("second"),
        )
    });
    assert_eq!(first.record_identity(), second.record_identity());
    assert_eq!(first.cross_covariance(), second.cross_covariance());

    let renamed = vec![
        QoiSensitivity::new("gap-z2", "millimetre", [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]).expect("a"),
    ];
    let third = with_default_cx(|cx| {
        propagate_pose_covariance(&calibrated, &renamed, coverage(), None, 0.0, "reason", cx)
    })
    .expect("third");
    assert_ne!(first.record_identity(), third.record_identity());

    let reworded = with_default_cx(|cx| {
        propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            None,
            0.0,
            "other reason",
            cx,
        )
    })
    .expect("fourth");
    assert_ne!(first.record_identity(), reworded.record_identity());
    assert_eq!(GEOMETRY_PROPAGATION_SCHEMA_VERSION, 1);
}

#[test]
fn g4_cancellation_is_structured() {
    let calibrated = axis_fixture(0.01);
    let qois = vec![
        QoiSensitivity::new("gap-z", "millimetre", [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]).expect("a"),
    ];
    assert!(matches!(
        with_cx(true, ExecMode::Deterministic, Budget::INFINITE, |cx| {
            propagate_pose_covariance(&calibrated, &qois, coverage(), None, 0.0, "reason", cx)
        }),
        Err(PropagateError::Cancelled {
            phase: "propagate.declarations"
        })
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn e2e_registration_to_full_budget_with_correlated_geometry() {
    // Register a noisy scan, calibrate, propagate three QoIs, and assemble
    // the complete eight-term budget with the geometry slot populated.
    let sigma = 1e-3f64;
    let mut rng = TestRng(0xE2E_0002);
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
    let fiducials: Vec<Fiducial3> = design
        .iter()
        .map(|p| {
            Fiducial3::new(
                *p,
                p3(
                    p.x() + sigma * rng.normal(),
                    p.y() + sigma * rng.normal(),
                    p.z() + sigma * rng.normal(),
                ),
            )
        })
        .collect();
    let covariance = Covariance3::new(sigma * sigma, 0.0, 0.0, sigma * sigma, 0.0, sigma * sigma)
        .expect("covariance");
    let model = MetrologyModel3::new(
        vec![covariance; fiducials.len()],
        CrossFiducialModel3::Independent,
        HuberPolicy::Disabled,
        "cmm-calibration-2026-07/block-fixture@rev3",
    )
    .expect("model");

    let gradients = vec![
        [0.0, 0.0, 1.0, 0.2, -0.1, 0.0],
        [0.0, 0.0, 1.0, -0.2, 0.1, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.8],
    ];
    let qois = vec![
        QoiSensitivity::new("t-junction-gap", "millimetre", gradients[0]).expect("a"),
        QoiSensitivity::new("contact-plane-gap", "millimetre", gradients[1]).expect("b"),
        QoiSensitivity::new("connector-offset-x", "millimetre", gradients[2]).expect("c"),
    ];
    let evaluator = LinearEvaluator {
        gradients: gradients.clone(),
    };

    let (calibrated, propagation) = with_default_cx(|cx| {
        let calibrated = estimate_calibrated_rigid3(&fiducials, &model, cx).expect("calibrated");
        let propagation = propagate_pose_covariance(
            &calibrated,
            &qois,
            coverage(),
            Some(&evaluator),
            1e-6,
            "",
            cx,
        )
        .expect("propagation");
        (calibrated, propagation)
    });

    // Correlation structure is retained: the two gap QoIs share the pose
    // translation almost fully.
    let rho = propagation.correlation(0, 1).expect("correlation");
    assert!(
        rho > 0.5,
        "shared translation must correlate the gaps: {rho}"
    );

    // Build the full 8-term budget: geometry from the propagation, the
    // numerical terms as certified intervals, everything else honestly
    // declared.
    let artifact = |role: &str| {
        UncertaintyArtifactRef::new(role, propagation.record_identity()).expect("artifact")
    };
    let mut terms: Vec<EngineeringUncertaintyTerm> = Vec::new();
    for kind in EngineeringUncertaintyKind::ALL {
        let term = match kind {
            EngineeringUncertaintyKind::Geometry => {
                propagation.geometry_term(0).expect("geometry term")
            }
            EngineeringUncertaintyKind::Roundoff
            | EngineeringUncertaintyKind::SolverAlgebraic
            | EngineeringUncertaintyKind::Discretization => EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::interval(0.0, 1e-9).expect("interval"),
                artifact("numerical-certificate-placeholder"),
            )
            .expect("numerical term"),
            EngineeringUncertaintyKind::Parameters
            | EngineeringUncertaintyKind::BoundaryConditions => {
                EngineeringUncertaintyTerm::try_new(
                    kind,
                    TermValue::unknown("not propagated in this as-built fixture").expect("unknown"),
                    artifact("declared-gap"),
                )
                .expect("unknown term")
            }
            EngineeringUncertaintyKind::ModelForm | EngineeringUncertaintyKind::Measurement => {
                EngineeringUncertaintyTerm::try_new(
                    kind,
                    TermValue::negligible("synthetic fixture with declared noise only")
                        .expect("negligible"),
                    artifact("fixture-declaration"),
                )
                .expect("negligible term")
            }
            // The kind enum is non-exhaustive upstream; any future source
            // this fixture does not understand is an honest evidence gap.
            _ => EngineeringUncertaintyTerm::try_new(
                kind,
                TermValue::unknown("source kind unknown to this fixture").expect("unknown"),
                artifact("declared-gap"),
            )
            .expect("wildcard term"),
        };
        terms.push(term);
    }
    let budget = EngineeringUncertaintyBudget::try_new("t-junction-gap", "millimetre", terms)
        .expect("budget");
    let geometry = budget.term(EngineeringUncertaintyKind::Geometry);
    assert!(matches!(geometry.value(), TermValue::Distribution(_)));
    assert_eq!(
        geometry.provenance().digest(),
        propagation.record_identity()
    );

    // Forensic log: per-stage table with the correlation structure.
    let mut log = String::new();
    let _ = writeln!(
        log,
        "registration: dof={} identity={}",
        calibrated.degrees_of_freedom(),
        calibrated.model_identity()
    );
    let _ = writeln!(
        log,
        "propagation: record={} method={:?} sd={:?}",
        propagation.record_identity(),
        propagation.method(),
        propagation.standard_deviations()
    );
    let _ = writeln!(
        log,
        "correlation: rho01={:.6} rho02={:.6} rho12={:.6}",
        propagation.correlation(0, 1).expect("rho01"),
        propagation.correlation(0, 2).expect("rho02"),
        propagation.correlation(1, 2).expect("rho12")
    );
    let _ = writeln!(log, "budget: {}", budget.render_report());
    for marker in [
        "registration: dof=18",
        "propagation: record=",
        "correlation: rho01=",
    ] {
        assert!(log.contains(marker), "forensic log lost {marker:?}\n{log}");
    }
    println!("{log}");
}
