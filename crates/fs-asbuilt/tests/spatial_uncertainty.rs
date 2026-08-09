//! G0/G3/G5 battery for calibrated as-built spatial uncertainty.

use fs_asbuilt::uncertainty::{
    BiasBound, Covariance2, CrossFiducialModel, DecisionReason, DecisionState,
    EvidenceAuthenticationError, EvidenceReceipt, EvidenceVerification, EvidenceVerifier,
    HuberPolicy, InspectionRelation, MetrologyModel, NoEvidenceVerifier, OutlierDisposition,
    REGISTRATION_UNCERTAINTY_SCHEMA_VERSION, SPATIAL_EVIDENCE_SCHEMA_VERSION,
    SpatialUncertaintyError, assess_calibrated_as_built, estimate_calibrated_registration,
};
use fs_asbuilt::{Fiducial, Point2};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn point(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).expect("finite fixture point")
}

fn covariance(xx: f64, xy: f64, yy: f64) -> Covariance2 {
    Covariance2::new(xx, xy, yy).expect("positive-definite fixture covariance")
}

fn cardinal_fiducials() -> Vec<Fiducial> {
    [
        point(1.0, 0.0),
        point(-1.0, 0.0),
        point(0.0, 1.0),
        point(0.0, -1.0),
    ]
    .into_iter()
    .map(|value| Fiducial::new(value, value))
    .collect()
}

fn rigid_map(value: Point2, angle: f64, tx: f64, ty: f64) -> Point2 {
    let (sine, cosine) = angle.sin_cos();
    point(
        cosine * value.x() - sine * value.y() + tx,
        sine * value.x() + cosine * value.y() + ty,
    )
}

fn model(
    count: usize,
    variance: f64,
    cross: CrossFiducialModel,
    huber: HuberPolicy,
    bias: BiasBound,
    identity: &str,
) -> MetrologyModel {
    MetrologyModel::new(
        vec![covariance(variance, 0.0, variance); count],
        cross,
        huber,
        bias,
        identity,
    )
    .expect("valid fixture model")
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5A71_A1,
                kernel_id: 94,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    });
    assert!(pool.stats().quiescent(), "fixture context leaked its arena");
    result
}

fn assert_near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} within {tolerance:.3e}"
    );
}

#[test]
fn g0_cardinal_geometry_matches_the_full_equicorrelation_oracle() {
    let fiducials = cardinal_fiducials();
    let independent = model(
        fiducials.len(),
        1.0,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "cardinal-independent-v1",
    );
    let independent_fit = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &independent, cx)
    })
    .expect("rank-two cardinal fit");
    assert_eq!(independent_fit.degrees_of_freedom(), 5);
    for diagonal in 0..3 {
        assert_near(
            independent_fit.covariance()[diagonal][diagonal],
            0.25,
            1e-12,
        );
    }
    assert_near(independent_fit.leverage().iter().sum(), 3.0, 1e-12);

    let rho = 0.2;
    let correlated = model(
        fiducials.len(),
        1.0,
        CrossFiducialModel::EquicorrelatedStandardized { rho },
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "cardinal-rho-v1",
    );
    let correlated_fit = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &correlated, cx)
    })
    .expect("equicorrelated cardinal fit");
    assert_near(
        correlated_fit.covariance()[0][0],
        (1.0 + 3.0 * rho) / 4.0,
        1e-12,
    );
    assert_near(
        correlated_fit.covariance()[1][1],
        (1.0 + 3.0 * rho) / 4.0,
        1e-12,
    );
    assert_near(correlated_fit.covariance()[2][2], (1.0 - rho) / 4.0, 1e-12);
    for (row, column) in [(0, 1), (0, 2), (1, 2)] {
        assert_near(correlated_fit.covariance()[row][column], 0.0, 1e-12);
    }
    assert_near(correlated_fit.leverage().iter().sum(), 3.0, 1e-12);
    assert_ne!(
        independent_fit.model_identity(),
        correlated_fit.model_identity(),
        "correlation is identity-semantic"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One paired unit/order metamorphic fixture shares exact inputs.
fn g3_global_heteroscedastic_fit_is_unit_equivariant_and_order_stable() {
    let angle = 0.37;
    let tx = 0.8;
    let ty = -0.45;
    let design = [
        point(-2.0, -1.0),
        point(0.5, -1.5),
        point(2.0, 0.25),
        point(-0.75, 2.5),
        point(3.0, 1.5),
    ];
    let fiducials: Vec<_> = design
        .into_iter()
        .map(|value| Fiducial::new(value, rigid_map(value, angle, tx, ty)))
        .collect();
    let covariances = vec![
        covariance(0.04, 0.012, 0.09),
        covariance(0.01, -0.004, 0.03),
        covariance(0.08, 0.02, 0.05),
        covariance(0.03, -0.009, 0.07),
        covariance(0.02, 0.006, 0.04),
    ];
    let calibrated = MetrologyModel::new(
        covariances.clone(),
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "heteroscedastic-offdiagonal-v1",
    )
    .unwrap();
    let fit = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &calibrated, cx)
    })
    .expect("global constrained GLS fit");
    assert_near(fit.registration().rotation_rad(), angle, 1e-12);
    assert_near(fit.registration().tx(), tx, 1e-12);
    assert_near(fit.registration().ty(), ty, 1e-12);
    for row in 0..3 {
        for column in 0..3 {
            assert_eq!(
                fit.covariance()[row][column].to_bits(),
                fit.covariance()[column][row].to_bits(),
                "published covariance must be bit-symmetric"
            );
        }
    }
    assert!(
        fit.outlier_diagnostics()
            .iter()
            .all(|diagnostic| diagnostic.disposition() == OutlierDisposition::NotEvaluated)
    );

    let distant_design = [point(1.0e6, -2.0e6)];
    let distant_scan = [rigid_map(distant_design[0], angle, tx, ty)];
    let distant_inspection = [covariance(0.01, 0.003, 0.02)];
    let distant = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &distant_design,
            &distant_scan,
            &distant_inspection,
            InspectionRelation::DisjointFromRegistration,
            f64::MAX,
            0.75,
            cx,
        )
    })
    .unwrap();
    let distant_bound = distant.point_bounds()[0];
    let naive_upper = distant_bound.observed_deviation()
        + (distant_bound.total_covariance().trace() / (1.0 - 0.75)).sqrt();
    let distant_boundary = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &distant_design,
            &distant_scan,
            &distant_inspection,
            InspectionRelation::DisjointFromRegistration,
            naive_upper,
            0.75,
            cx,
        )
    })
    .unwrap();
    assert_ne!(
        distant_boundary.decision().state(),
        DecisionState::WithinTolerance,
        "nonzero-angle affine/trigonometric rounding must be outward"
    );

    let scale = 1_000.0;
    let scaled_fiducials: Vec<_> = fiducials
        .iter()
        .map(|fiducial| {
            Fiducial::new(
                point(scale * fiducial.design().x(), scale * fiducial.design().y()),
                point(
                    scale * fiducial.measured().x(),
                    scale * fiducial.measured().y(),
                ),
            )
        })
        .collect();
    let scaled_covariances: Vec<_> = covariances
        .iter()
        .map(|entry| {
            covariance(
                scale * scale * entry.xx(),
                scale * scale * entry.xy(),
                scale * scale * entry.yy(),
            )
        })
        .collect();
    let scaled_model = MetrologyModel::new(
        scaled_covariances,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "heteroscedastic-offdiagonal-scaled-v1",
    )
    .unwrap();
    let scaled = with_cx(false, |cx| {
        estimate_calibrated_registration(&scaled_fiducials, &scaled_model, cx)
    })
    .expect("length-unit rescaling remains admissible");
    assert_near(scaled.registration().rotation_rad(), angle, 1e-12);
    assert_near(scaled.registration().tx(), scale * tx, 1e-9);
    assert_near(scaled.registration().ty(), scale * ty, 1e-9);
    assert_near(
        scaled.covariance()[0][0],
        scale * scale * fit.covariance()[0][0],
        1e-8,
    );
    assert_near(scaled.covariance()[2][2], fit.covariance()[2][2], 1e-12);

    let mut reversed_fiducials = fiducials;
    let mut reversed_covariances = covariances;
    reversed_fiducials.reverse();
    reversed_covariances.reverse();
    let reversed_model = MetrologyModel::new(
        reversed_covariances,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "heteroscedastic-offdiagonal-v1",
    )
    .unwrap();
    let reversed = with_cx(false, |cx| {
        estimate_calibrated_registration(&reversed_fiducials, &reversed_model, cx)
    })
    .unwrap();
    assert_near(
        reversed.registration().rotation_rad(),
        fit.registration().rotation_rad(),
        1e-12,
    );
    assert_near(reversed.registration().tx(), fit.registration().tx(), 1e-12);
    assert_near(reversed.registration().ty(), fit.registration().ty(), 1e-12);
}

#[test]
#[allow(clippy::too_many_lines)] // One ordered fail-closed matrix covers interacting model domains.
fn g0_covariance_correlation_and_geometry_domains_fail_closed() {
    assert!(matches!(
        Covariance2::new(1.0, 1.0, 1.0),
        Err(SpatialUncertaintyError::NonPositiveDefiniteCovariance { .. })
    ));
    let count = 4;
    for rho in [-1.0 / 3.0, 1.0, f64::NAN] {
        assert!(matches!(
            MetrologyModel::new(
                vec![covariance(1.0, 0.0, 1.0); count],
                CrossFiducialModel::EquicorrelatedStandardized { rho },
                HuberPolicy::Disabled,
                BiasBound::Bounded(0.0),
                "bad-rho-fixture",
            ),
            Err(SpatialUncertaintyError::InvalidScalar {
                field: "cross_fiducial.rho",
                ..
            })
        ));
    }
    let huber = HuberPolicy::new(1.5, 4).expect("valid Huber policy");
    assert_eq!(
        MetrologyModel::new(
            vec![covariance(1.0, 0.0, 1.0); count],
            CrossFiducialModel::EquicorrelatedStandardized { rho: 0.1 },
            huber,
            BiasBound::Bounded(0.0),
            "unsupported-robust-correlation",
        ),
        Err(SpatialUncertaintyError::RobustCorrelationUnsupported)
    );
    for invalid_huber in [
        HuberPolicy::Enabled {
            threshold: f64::NAN,
            iterations: 1,
        },
        HuberPolicy::Enabled {
            threshold: -1.0,
            iterations: 1,
        },
        HuberPolicy::Enabled {
            threshold: 1.0,
            iterations: 0,
        },
        HuberPolicy::Enabled {
            threshold: 1.0,
            iterations: 33,
        },
    ] {
        assert!(matches!(
            MetrologyModel::new(
                vec![covariance(1.0, 0.0, 1.0); count],
                CrossFiducialModel::Independent,
                invalid_huber,
                BiasBound::Bounded(0.0),
                "invalid-direct-huber",
            ),
            Err(SpatialUncertaintyError::InvalidScalar { .. })
        ));
    }

    let unknown_model = model(
        count,
        1.0,
        CrossFiducialModel::Unknown,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "unknown-dependence-v1",
    );
    assert_eq!(
        with_cx(false, |cx| estimate_calibrated_registration(
            &cardinal_fiducials(),
            &unknown_model,
            cx,
        )),
        Err(SpatialUncertaintyError::UnknownDependence)
    );

    let collinear: Vec<_> = [-1.0, 0.0, 1.0]
        .into_iter()
        .map(|x| {
            let value = point(x, 0.0);
            Fiducial::new(value, value)
        })
        .collect();
    let collinear_model = model(
        collinear.len(),
        1.0,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "collinear-model-v1",
    );
    assert_eq!(
        with_cx(false, |cx| estimate_calibrated_registration(
            &collinear,
            &collinear_model,
            cx,
        )),
        Err(SpatialUncertaintyError::SingularInformation)
    );

    let rotation_ambiguous: Vec<_> = cardinal_fiducials()
        .into_iter()
        .map(|fiducial| Fiducial::new(fiducial.design(), point(0.0, 0.0)))
        .collect();
    let ambiguous_model = model(
        rotation_ambiguous.len(),
        1.0,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "ambiguous-rotation-v1",
    );
    assert_eq!(
        with_cx(false, |cx| estimate_calibrated_registration(
            &rotation_ambiguous,
            &ambiguous_model,
            cx,
        )),
        Err(SpatialUncertaintyError::AmbiguousGlobalMinimum)
    );

    let near_hard = vec![
        Fiducial::new(point(1.0, 0.0), point(0.1, 0.2)),
        Fiducial::new(point(-1.0, 0.0), point(-0.1, -0.2)),
        Fiducial::new(point(0.0, 1.0), point(0.0, 0.0)),
        Fiducial::new(point(0.0, -1.0), point(0.0, 0.0)),
    ];
    let near_hard_model = MetrologyModel::new(
        vec![
            covariance(0.68, -0.24, 0.32),
            covariance(0.68, -0.24, 0.32),
            covariance(0.32, 0.24, 0.68),
            covariance(0.32, 0.24, 0.68),
        ],
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "near-hard-rotation-v1",
    )
    .unwrap();
    assert!(matches!(
        with_cx(false, |cx| estimate_calibrated_registration(
            &near_hard,
            &near_hard_model,
            cx,
        )),
        Err(SpatialUncertaintyError::AmbiguousGlobalMinimum
            | SpatialUncertaintyError::SingularInformation)
    ));
}

#[test]
fn robust_refit_downweights_outlier_but_keeps_coverage_conditional() {
    let mut fiducials = cardinal_fiducials();
    fiducials.push(Fiducial::new(point(2.0, 2.0), point(20.0, -15.0)));
    let unweighted_model = model(
        fiducials.len(),
        0.01,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "unweighted-outlier-v1",
    );
    let unweighted = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &unweighted_model, cx)
    })
    .unwrap();
    let one_pass_model = model(
        fiducials.len(),
        0.01,
        CrossFiducialModel::Independent,
        HuberPolicy::new(1.5, 1).expect("one bounded weight refresh"),
        BiasBound::Bounded(0.0),
        "one-pass-outlier-v1",
    );
    let one_pass = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &one_pass_model, cx)
    })
    .expect("last refreshed weights receive a final global solve");
    let pose_error = |fit: &fs_asbuilt::uncertainty::CalibratedRegistration| {
        fit.registration().tx().hypot(fit.registration().ty())
            + fit.registration().rotation_rad().abs()
    };
    assert!(one_pass.weights()[4] < 1.0);
    assert!(
        pose_error(&one_pass) < pose_error(&unweighted),
        "one-pass output must not publish the pre-refresh unweighted transform"
    );

    let robust_model = model(
        fiducials.len(),
        0.01,
        CrossFiducialModel::Independent,
        HuberPolicy::new(1.5, 12).expect("bounded robust policy"),
        BiasBound::Bounded(0.0),
        "robust-outlier-v1",
    );
    let fit = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &robust_model, cx)
    })
    .expect("robust fit completes");
    assert!(fit.weights()[4] < 0.1, "gross outlier must be downweighted");
    assert_eq!(
        fit.outlier_diagnostics()[4].disposition(),
        OutlierDisposition::Downweighted
    );
    assert!(fit.outlier_diagnostics()[4].standardized_residual_norm() > 1.5);
    assert!(fit.robust_conditional());

    let design = [point(0.0, 0.0)];
    let evidence = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &design,
            &design,
            &[covariance(0.01, 0.0, 0.01)],
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.9,
            cx,
        )
    })
    .expect("conditional assessment is represented honestly");
    assert_eq!(evidence.decision().state(), DecisionState::Indeterminate);
    assert_eq!(
        evidence.decision().reason(),
        DecisionReason::AdaptiveWeightsConditional
    );
    assert!(evidence.point_bounds().is_empty());
}

fn narrow_fit(bias: BiasBound) -> fs_asbuilt::uncertainty::CalibratedRegistration {
    let fiducials = cardinal_fiducials();
    let noise = model(
        fiducials.len(),
        1e-4,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        bias,
        "narrow-cardinal-v1",
    );
    with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &noise, cx)
    })
    .expect("narrow cardinal fit")
}

#[test]
#[allow(clippy::too_many_lines)] // One propagation oracle shares the cardinal covariance baseline.
fn propagation_uses_pose_plus_inspection_once_and_family_size_widens() {
    let fit = narrow_fit(BiasBound::Bounded(0.0));
    assert_eq!(
        fit.registration().residual_rms().to_bits(),
        0.0f64.to_bits(),
        "perfect-fit residual is not inflated into covariance"
    );
    let one_design = [point(0.0, 0.0)];
    let one_covariance = [covariance(1e-4, 0.0, 1e-4)];
    let one = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &one_design,
            &one_design,
            &one_covariance,
            InspectionRelation::DisjointFromRegistration,
            0.1,
            0.75,
            cx,
        )
    })
    .expect("one-point simultaneous assessment");
    assert_eq!(one.decision().state(), DecisionState::WithinTolerance);
    let bound = one.point_bounds()[0];
    assert_near(bound.total_covariance().xx(), 1.25e-4, 1e-15);
    assert_near(bound.total_covariance().yy(), 1.25e-4, 1e-15);
    let plain_radius = (bound.total_covariance().trace() / (1.0 - 0.75)).sqrt();
    assert!(bound.simultaneous_radius() >= plain_radius);
    assert!(bound.lower() <= 0.0);
    assert!(bound.upper() >= plain_radius);

    let boundary = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &one_design,
            &one_design,
            &one_covariance,
            InspectionRelation::DisjointFromRegistration,
            plain_radius,
            0.75,
            cx,
        )
    })
    .expect("outward-rounded boundary assessment");
    assert_ne!(
        boundary.decision().state(),
        DecisionState::WithinTolerance,
        "round-to-nearest equality must not create a false acceptance"
    );

    let four_design = [point(0.0, 0.0); 4];
    let four_covariance = [covariance(1e-4, 0.0, 1e-4); 4];
    let four = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &four_design,
            &four_design,
            &four_covariance,
            InspectionRelation::DisjointFromRegistration,
            0.2,
            0.75,
            cx,
        )
    })
    .expect("four-point simultaneous assessment");
    assert_near(
        four.point_bounds()[0].simultaneous_radius(),
        2.0 * bound.simultaneous_radius(),
        1e-14,
    );

    let far_design = [point(10.0, 0.0)];
    let far = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &far_design,
            &far_design,
            &one_covariance,
            InspectionRelation::DisjointFromRegistration,
            10.0,
            0.75,
            cx,
        )
    })
    .unwrap();
    assert_near(
        far.point_bounds()[0].total_covariance().xx(),
        1.25e-4,
        1e-15,
    );
    assert_near(
        far.point_bounds()[0].total_covariance().yy(),
        2.625e-3,
        1e-15,
    );

    let biased_fit = narrow_fit(BiasBound::Bounded(0.2));
    let biased = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &biased_fit,
            &one_design,
            &one_design,
            &one_covariance,
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.75,
            cx,
        )
    })
    .unwrap();
    assert!(biased.point_bounds()[0].simultaneous_radius() >= bound.simultaneous_radius() + 0.2);
}

#[test]
fn calibrated_bounds_produce_all_three_decision_states() {
    let fit = narrow_fit(BiasBound::Bounded(0.0));
    let design = [point(0.0, 0.0)];
    let covariance = [covariance(1e-4, 0.0, 1e-4)];
    let assess = |scan_y: f64, tolerance: f64| {
        with_cx(false, |cx| {
            assess_calibrated_as_built(
                &fit,
                &design,
                &[point(0.0, scan_y)],
                &covariance,
                InspectionRelation::DisjointFromRegistration,
                tolerance,
                0.75,
                cx,
            )
        })
        .expect("finite decision fixture")
    };
    assert_eq!(
        assess(0.0, 0.1).decision().state(),
        DecisionState::WithinTolerance
    );
    assert_eq!(
        assess(1.0, 0.5).decision().state(),
        DecisionState::ExceedsTolerance
    );
    assert_eq!(
        assess(0.05, 0.05).decision().state(),
        DecisionState::Indeterminate
    );
}

#[test]
fn unknown_overlap_and_unbounded_bias_are_explicit_no_claims() {
    let design = [point(0.0, 0.0)];
    let covariance = [covariance(1e-4, 0.0, 1e-4)];
    let overlap_fit = narrow_fit(BiasBound::Bounded(0.0));
    let overlap = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &overlap_fit,
            &design,
            &design,
            &covariance,
            InspectionRelation::UnknownOrOverlapping,
            1.0,
            0.9,
            cx,
        )
    })
    .expect("overlap becomes an explicit candidate");
    assert_eq!(
        overlap.decision().reason(),
        DecisionReason::RegistrationInspectionDependence
    );
    assert_eq!(overlap.decision().lower(), None);
    assert_eq!(overlap.decision().upper(), None);

    let unbounded_fit = narrow_fit(BiasBound::Unbounded);
    let unbounded = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &unbounded_fit,
            &design,
            &design,
            &covariance,
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.9,
            cx,
        )
    })
    .expect("unbounded bias becomes an explicit candidate");
    assert_eq!(unbounded.decision().reason(), DecisionReason::UnboundedBias);
    assert!(unbounded.point_bounds().is_empty());
}

#[test]
fn g3_permutation_preserves_numerics_and_g5_semantics_move_identities() {
    let fiducials = cardinal_fiducials();
    let baseline_model = model(
        fiducials.len(),
        0.01,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "identity-metamorphic-v1",
    );
    let baseline = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &baseline_model, cx)
    })
    .unwrap();
    let mut reversed = fiducials.clone();
    reversed.reverse();
    let permuted = with_cx(false, |cx| {
        estimate_calibrated_registration(&reversed, &baseline_model, cx)
    })
    .unwrap();
    assert_eq!(baseline.covariance(), permuted.covariance());
    assert_eq!(baseline.registration(), permuted.registration());
    assert_ne!(baseline.model_identity(), permuted.model_identity());

    let changed_model = model(
        fiducials.len(),
        0.02,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "identity-metamorphic-v1",
    );
    let changed = with_cx(false, |cx| {
        estimate_calibrated_registration(&fiducials, &changed_model, cx)
    })
    .unwrap();
    assert_ne!(baseline.model_identity(), changed.model_identity());

    let design = [point(0.0, 0.0)];
    let inspection = [covariance(0.01, 0.0, 0.01)];
    let first = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &baseline,
            &design,
            &design,
            &inspection,
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.75,
            cx,
        )
    })
    .unwrap();
    let replay = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &baseline,
            &design,
            &design,
            &inspection,
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.75,
            cx,
        )
    })
    .unwrap();
    let confidence_changed = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &baseline,
            &design,
            &design,
            &inspection,
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.8,
            cx,
        )
    })
    .unwrap();
    assert_eq!(first, replay);
    assert_ne!(
        first.evidence_identity(),
        confidence_changed.evidence_identity()
    );
    assert_eq!(REGISTRATION_UNCERTAINTY_SCHEMA_VERSION, 1);
}

struct AcceptExactPolicy(fs_blake3::ContentHash);

impl EvidenceVerifier for AcceptExactPolicy {
    fn verify(
        &self,
        _evidence: &fs_asbuilt::uncertainty::CalibratedAsBuiltEvidence,
        _receipt: &EvidenceReceipt,
    ) -> EvidenceVerification {
        EvidenceVerification::accept(self.0)
    }
}

#[test]
fn g5_receipt_authentication_fails_closed_on_every_independent_binding() {
    let fit = narrow_fit(BiasBound::Bounded(0.0));
    let design = [point(0.0, 0.0)];
    let evidence = with_cx(false, |cx| {
        assess_calibrated_as_built(
            &fit,
            &design,
            &design,
            &[covariance(1e-4, 0.0, 1e-4)],
            InspectionRelation::DisjointFromRegistration,
            1.0,
            0.75,
            cx,
        )
    })
    .unwrap();
    let policy = fs_blake3::hash_domain("fs-asbuilt.test.policy", b"exact");
    let receipt = EvidenceReceipt::from_parts(
        evidence.evidence_identity(),
        SPATIAL_EVIDENCE_SCHEMA_VERSION,
        policy,
    );
    assert!(matches!(
        evidence.clone().authenticate(receipt, &NoEvidenceVerifier),
        Err(EvidenceAuthenticationError::Refused { .. })
    ));
    let authenticated = evidence
        .clone()
        .authenticate(receipt, &AcceptExactPolicy(policy))
        .expect("matching injected authority admits exact lineage");
    assert_eq!(
        authenticated.evidence().evidence_identity(),
        evidence.evidence_identity()
    );

    let other_root = fs_blake3::hash_domain("fs-asbuilt.test.evidence", b"other");
    assert_eq!(
        evidence.clone().authenticate(
            EvidenceReceipt::from_parts(other_root, SPATIAL_EVIDENCE_SCHEMA_VERSION, policy),
            &AcceptExactPolicy(policy),
        ),
        Err(EvidenceAuthenticationError::EvidenceMismatch)
    );
    assert!(matches!(
        evidence.clone().authenticate(
            EvidenceReceipt::from_parts(evidence.evidence_identity(), 99, policy),
            &AcceptExactPolicy(policy),
        ),
        Err(EvidenceAuthenticationError::SchemaMismatch { .. })
    ));
    let other_policy = fs_blake3::hash_domain("fs-asbuilt.test.policy", b"other");
    assert!(matches!(
        evidence.authenticate(receipt, &AcceptExactPolicy(other_policy)),
        Err(EvidenceAuthenticationError::PolicyMismatch { .. })
    ));
}

#[test]
fn cancellation_refuses_before_publication() {
    let fiducials = cardinal_fiducials();
    let model = model(
        fiducials.len(),
        1.0,
        CrossFiducialModel::Independent,
        HuberPolicy::Disabled,
        BiasBound::Bounded(0.0),
        "cancel-spatial-v1",
    );
    assert_eq!(
        with_cx(true, |cx| estimate_calibrated_registration(
            &fiducials, &model, cx,
        )),
        Err(SpatialUncertaintyError::Cancelled {
            phase: "registration rank scale",
        })
    );
}

// sj31i.7.2: the uncertainty lane's covariance dimensions bind to the
// shared inference core — planar covariance entries are length-squared,
// information entries inverse-length-squared, derived mechanically.
#[test]
fn uncertainty_covariance_dims_match_the_shared_core() {
    use fs_asbuilt::dimensioned::planar_covariance_schema;
    let covariance = planar_covariance_schema().expect("planar schema");
    let entry = covariance.entry_dims(0, 1).expect("entry dims");
    assert_eq!(entry, fs_qty::Dims([2, 0, 0, 0, 0, 0]));
    let information = covariance
        .information_entry_dims(1, 1)
        .expect("information dims");
    assert_eq!(information, fs_qty::Dims([-2, 0, 0, 0, 0, 0]));
    // A calibrated Covariance2 lives in the squared reading dimensions of a
    // length observation, so its entries are exactly the core's planar
    // covariance dimensions.
    let reading = fs_qty::inference::ObservationSchema::new(
        fs_qty::semantic::QuantitySpec::dimensional(fs_qty::Dims([1, 0, 0, 0, 0, 0])),
    );
    assert_eq!(
        reading.noise_variance_dims().expect("variance dims"),
        entry,
        "Covariance2 entries carry squared reading dimensions"
    );
}
