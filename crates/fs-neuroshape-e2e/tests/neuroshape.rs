//! End-to-end battery: a learned neural implicit whose Lipschitz bound, safe
//! rendering, and existence of an enclosed negative component are PROVEN.
//! The global component count deliberately remains inexact.

use fs_exec::BudgetRefusal;
use fs_neuroshape_e2e::{
    CampaignError, CampaignParameter, CancellationKind, ComponentCountEvidence,
    LocalizationDiagnostic, LocalizationStage, NEUROSHAPE_COMPONENT_EVIDENCE_SCHEMA_VERSION,
    NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION, StageDetail, SurfaceLocalization,
    SurfaceLocalizationStatus, blob_sdf_net, iso_contour_resource_code, run_campaign,
    try_run_campaign,
};
use fs_rep_neural::{Layer, MlpSdf, SAFE_STEP_POLICY_VERSION, SafeStepStatus};
use fs_viz::{Grid2Error, IsoContourError, IsoContourResource};

#[test]
fn component_evidence_schema_versions_the_lower_bound_semantics() {
    assert_eq!(NEUROSHAPE_COMPONENT_EVIDENCE_SCHEMA_VERSION, 1);
}

#[test]
fn an_enclosed_component_lower_bound_is_certified() {
    let net = blob_sdf_net();
    let report = run_campaign(&net, 2.5, 0.3);
    // a finite certified Lipschitz bound underwrites everything.
    assert!(
        report.lipschitz.is_finite() && report.lipschitz > 0.0,
        "L {}",
        report.lipschitz
    );
    // sound sphere tracing: the origin is negative and the safe step is a
    // positive, finite, non-tunneling distance.
    assert!(report.origin_value < 0.0, "origin {}", report.origin_value);
    assert_eq!(report.safe_step.status(), SafeStepStatus::SignSeparated);
    assert!(report.safe_step.magnitude_lower_bound() > 0.0);
    assert!(report.safe_step.radius() > 0.0 && report.safe_step.radius().is_finite());
    assert_eq!(
        report.safe_step.enclosure(),
        net.eval_interval(&[0.0, 0.0], &[0.0, 0.0])
    );
    assert_eq!(report.safe_step_policy_version, SAFE_STEP_POLICY_VERSION);
    assert_eq!(report.safe_step.policy_version(), SAFE_STEP_POLICY_VERSION);
    assert_eq!(report.safe_step.policy(), report.safe_step_policy);
    assert_eq!(report.field_identity, net.identity());
    // Sampled localization is a useful independent corroboration, not the
    // theorem. The no-tunnel authority comes from the interval sign margin and
    // Lipschitz bound above.
    assert!(
        report.safe_step.radius() < report.nearest_surface_radius,
        "safe {} !< nearest surface {}",
        report.safe_step.radius(),
        report.nearest_surface_radius
    );
    assert!(report.nearest_surface_radius <= report.max_crossing_radius);
    // TOPOLOGY, PROVEN: a certified-inside interior enclosed by a CLOSED,
    // fully-certified boundary frame implies at least one enclosed component.
    // It does not exclude additional components inside or outside the frame.
    assert!(
        report.certified_inside,
        "inside interval {:?}",
        report.inside_interval
    );
    assert_eq!(report.boundary_segments, 4);
    assert_eq!(report.boundary_certified, report.boundary_segments);
    assert!(report.boundary_frame_certified);
    assert_eq!(report.component_count_evidence.lower_bound(), 1);
    assert_eq!(report.component_count_evidence.exact_count(), None);
    let ComponentCountEvidence::LowerBound(enclosed) = &report.component_count_evidence else {
        panic!("closed frame must produce a typed lower-bound witness");
    };
    assert_eq!(
        enclosed.central_box_half_width().to_bits(),
        0.3_f64.to_bits()
    );
    let enclosed_interval = enclosed.central_box_interval();
    assert_eq!(
        enclosed_interval.0.to_bits(),
        report.inside_interval.0.to_bits()
    );
    assert_eq!(
        enclosed_interval.1.to_bits(),
        report.inside_interval.1.to_bits()
    );
    assert_eq!(
        enclosed.boundary_frame_outer_half_width().to_bits(),
        2.5_f64.to_bits()
    );
    assert_eq!(
        enclosed.boundary_frame_inner_half_width().to_bits(),
        2.1_f64.to_bits()
    );
    assert!(
        enclosed
            .boundary_strip_intervals()
            .iter()
            .all(|(lo, hi)| lo.is_finite() && hi.is_finite() && *lo > 0.0 && lo <= hi)
    );
    // The origin Hessian is positive definite, but no zero-gradient witness is
    // present. This curvature check must never promote the lower bound to an
    // exact count or even claim a critical point.
    assert!(report.origin_hessian_positive_definite);
    // the visualization localizes a closed surface, all inside the ring.
    assert!(report.surface_crossings > 0);
    assert!(
        report.max_crossing_radius < 2.5,
        "surface escaped the ring: {}",
        report.max_crossing_radius
    );
    println!(
        "{{\"campaign\":\"neuroshapecert\",\"L\":{:.3},\"origin\":{:.3},\"safe_radius\":{:.3},\
         \"inside\":[{:.3},{:.3}],\"boundary\":{}/{},\"component_count_lower_bound\":{},\
         \"exact_component_count\":null,\"origin_hessian_positive_definite\":{},\
         \"crossings\":{},\"max_r\":{:.3}}}",
        report.lipschitz,
        report.origin_value,
        report.safe_step.radius(),
        report.inside_interval.0,
        report.inside_interval.1,
        report.boundary_certified,
        report.boundary_segments,
        report.component_count_evidence.lower_bound(),
        report.origin_hessian_positive_definite,
        report.surface_crossings,
        report.max_crossing_radius,
    );
}

#[test]
fn an_open_ring_yields_no_topology_certificate() {
    // too small a box: its boundary frame overlaps the surface → not certified.
    let net = blob_sdf_net();
    let report = run_campaign(&net, 0.3, 0.3);
    assert!(!report.boundary_frame_certified || !report.certified_inside);
    assert!(matches!(
        &report.component_count_evidence,
        ComponentCountEvidence::Unknown
    ));
    assert_eq!(report.component_count_evidence.lower_bound(), 0);
    assert_eq!(report.component_count_evidence.exact_count(), None);
}

#[test]
fn the_campaign_is_deterministic() {
    let net = blob_sdf_net();
    let a = run_campaign(&net, 2.5, 0.3);
    let b = run_campaign(&net, 2.5, 0.3);
    assert_eq!(a.lipschitz.to_bits(), b.lipschitz.to_bits());
    assert_eq!(a.surface_crossings, b.surface_crossings);
    assert_eq!(a.field_identity, b.field_identity);
    assert_eq!(a.safe_step.status(), b.safe_step.status());
    assert_eq!(
        a.safe_step.radius().to_bits(),
        b.safe_step.radius().to_bits()
    );
    assert_eq!(
        a.safe_step.magnitude_lower_bound().to_bits(),
        b.safe_step.magnitude_lower_bound().to_bits()
    );
    assert_eq!(a.component_count_evidence, b.component_count_evidence);
}

#[test]
fn campaign_admission_refuses_wrong_dimension_and_invalid_geometry() {
    let one_dimensional = MlpSdf::new(vec![Layer::new(vec![vec![1.0]], vec![0.0])], 1.0);
    // `NeuroShapeReport` is deliberately not `PartialEq` (a certificate report
    // is not a comparable value), so admission refusals are compared through
    // the error side of the `Result`.
    assert_eq!(
        try_run_campaign(&one_dimensional, 2.5, 0.3).err(),
        Some(CampaignError::InputDimension {
            expected: 2,
            actual: 1,
        })
    );

    let net = blob_sdf_net();
    for (ring_r, inner, expected) in [
        (
            f64::NAN,
            0.3,
            CampaignError::NonFiniteParameter(CampaignParameter::RingRadius),
        ),
        (
            f64::INFINITY,
            0.3,
            CampaignError::NonFiniteParameter(CampaignParameter::RingRadius),
        ),
        (
            2.5,
            f64::NEG_INFINITY,
            CampaignError::NonFiniteParameter(CampaignParameter::InnerHalfWidth),
        ),
        (
            0.0,
            0.3,
            CampaignError::OutOfRangeParameter(CampaignParameter::RingRadius),
        ),
        (
            2.5,
            -f64::MIN_POSITIVE,
            CampaignError::OutOfRangeParameter(CampaignParameter::InnerHalfWidth),
        ),
    ] {
        assert_eq!(try_run_campaign(&net, ring_r, inner).err(), Some(expected));
    }
}

#[test]
fn localization_schema_version_pins_the_wire_codes() {
    assert_eq!(NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION, 1);
    assert_eq!(SurfaceLocalizationStatus::Localized.code(), 1);
    assert_eq!(SurfaceLocalizationStatus::ValidEmpty.code(), 2);
    assert_eq!(SurfaceLocalizationStatus::InvalidInput.code(), 3);
    assert_eq!(SurfaceLocalizationStatus::Unrepresentable.code(), 4);
    assert_eq!(SurfaceLocalizationStatus::ResourceRefused.code(), 5);
    assert_eq!(SurfaceLocalizationStatus::Cancelled.code(), 6);
    assert_eq!(SurfaceLocalizationStatus::AllocationRefused.code(), 7);
    assert_eq!(SurfaceLocalizationStatus::InternalFault.code(), 8);
    assert_eq!(LocalizationStage::GridConstruction.code(), 1);
    assert_eq!(LocalizationStage::IsoContourExtraction.code(), 2);
}

/// Every `Grid2Error` variant must cross the boundary as the documented
/// typed outcome with its bounded detail — never as an erased `None`.
#[test]
fn every_grid_error_maps_to_a_typed_localization_outcome() {
    let cases: [(Grid2Error, SurfaceLocalizationStatus, LocalizationDiagnostic); 7] = [
        (
            Grid2Error::InvalidDimensions { dimensions: [1, 9] },
            SurfaceLocalizationStatus::InvalidInput,
            LocalizationDiagnostic::GridInvalidDimensions,
        ),
        (
            Grid2Error::NodeCountOverflow {
                dimensions: [usize::MAX, 3],
            },
            SurfaceLocalizationStatus::InvalidInput,
            LocalizationDiagnostic::GridNodeCountOverflow,
        ),
        (
            Grid2Error::NodeBudgetExceeded {
                required: 100,
                limit: 64,
            },
            SurfaceLocalizationStatus::ResourceRefused,
            LocalizationDiagnostic::GridNodeBudgetExceeded,
        ),
        (
            Grid2Error::InvalidBounds {
                axis: 1,
                lower: 0.0,
                upper: f64::NAN,
            },
            SurfaceLocalizationStatus::InvalidInput,
            LocalizationDiagnostic::GridInvalidBounds,
        ),
        (
            Grid2Error::UnrepresentableCoordinates {
                axis: 0,
                first_index: 3,
                first: 0.5,
                second_index: 4,
                second: 0.5,
            },
            SurfaceLocalizationStatus::Unrepresentable,
            LocalizationDiagnostic::GridUnrepresentableCoordinates,
        ),
        (
            Grid2Error::NonFiniteValue {
                index: 17,
                value: f64::NAN,
            },
            SurfaceLocalizationStatus::Unrepresentable,
            LocalizationDiagnostic::GridNonFiniteValue,
        ),
        (
            Grid2Error::AllocationFailed { nodes: 4096 },
            SurfaceLocalizationStatus::AllocationRefused,
            LocalizationDiagnostic::GridAllocationFailed,
        ),
    ];
    for (error, expected_status, expected_diagnostic) in cases {
        let outcome = SurfaceLocalization::from(error);
        assert_eq!(outcome.status(), expected_status, "{error}");
        let stage = outcome.stage().expect("refusals name their stage");
        assert_eq!(stage, LocalizationStage::GridConstruction, "{error}");
        match &outcome {
            SurfaceLocalization::InvalidInput(detail)
            | SurfaceLocalization::Unrepresentable(detail)
            | SurfaceLocalization::ResourceRefused(detail)
            | SurfaceLocalization::AllocationRefused(detail)
            | SurfaceLocalization::InternalFault(detail) => {
                assert_eq!(detail.diagnostic, expected_diagnostic, "{error}");
            }
            other => panic!("expected refusal detail for {error}, got {other:?}"),
        }
    }
}

/// Every non-cancellation `IsoContourError` variant maps to its documented
/// typed outcome; plan overflow keeps a stable per-resource auxiliary code.
#[test]
fn every_contour_error_maps_to_a_typed_localization_outcome() {
    // Resource codes are frozen by NEUROSHAPE_LOCALIZATION_SCHEMA_VERSION = 1.
    assert_eq!(
        iso_contour_resource_code(IsoContourResource::Cells),
        1
    );
    assert_eq!(
        iso_contour_resource_code(IsoContourResource::WorkUnits),
        13
    );

    let overflow = IsoContourError::PlanOverflow {
        resource: IsoContourResource::EdgeVisits,
    };
    match SurfaceLocalization::from(overflow) {
        SurfaceLocalization::InternalFault(detail) => {
            assert_eq!(detail.diagnostic, LocalizationDiagnostic::IsoPlanOverflow);
            assert_eq!(
                detail.aux,
                iso_contour_resource_code(IsoContourResource::EdgeVisits)
            );
        }
        other => panic!("plan overflow must be an internal fault, got {other:?}"),
    }

    let coincident = IsoContourError::CoincidentLevelEdge {
        first: [2, 3],
        second: [2, 4],
    };
    match SurfaceLocalization::from(coincident) {
        SurfaceLocalization::Unrepresentable(detail) => {
            assert_eq!(detail.diagnostic, LocalizationDiagnostic::IsoCoincidentLevelEdge);
            assert_eq!(detail.first_index, (2u64 << 32) | 3);
            assert_eq!(detail.second_index, (2u64 << 32) | 4);
        }
        other => panic!("coincident edge must be unrepresentable, got {other:?}"),
    }

    let budget = IsoContourError::OperationBudgetExceeded {
        resource: IsoContourResource::LiveBytes,
        required: u128::from(u64::MAX) + 7,
        limit: 1024,
    };
    match SurfaceLocalization::from(budget) {
        SurfaceLocalization::ResourceRefused(detail) => {
            assert_eq!(
                detail.diagnostic,
                LocalizationDiagnostic::IsoOperationBudgetExceeded
            );
            assert_eq!(detail.required, u64::MAX, "u128 requirements saturate");
            assert_eq!(detail.limit, 1024);
        }
        other => panic!("operation budgets are resource refusals, got {other:?}"),
    }

    let unrepresentable = IsoContourError::UnrepresentableIntersection {
        first: [1, 1],
        second: [1, 2],
        first_point_bits: [0, 0],
        second_point_bits: [0, 0],
        first_value_bits: 1,
        second_value_bits: 2,
        iso_bits: 3,
        first_distance_bits: 4,
        second_distance_bits: 5,
        interpolation_bits: 6,
        point_bits: [7, 8],
        collapsed_axis: 1,
    };
    match SurfaceLocalization::from(unrepresentable) {
        SurfaceLocalization::Unrepresentable(detail) => {
            assert_eq!(
                detail.diagnostic,
                LocalizationDiagnostic::IsoUnrepresentableIntersection
            );
            assert_eq!(detail.first_index, (1u64 << 32) | 1);
            assert_eq!(detail.second_index, (1u64 << 32) | 2);
            assert_eq!(detail.scalar_bits, 3, "level bits retained");
            assert_eq!(detail.second_bits, 6, "interpolation bits retained");
            assert_eq!(detail.axis, 1, "collapse axis retained");
        }
        other => panic!("unrepresentable intersections keep detail, got {other:?}"),
    }

    let allocation = IsoContourError::AllocationFailed { required: 33 };
    match SurfaceLocalization::from(allocation) {
        SurfaceLocalization::AllocationRefused(detail) => {
            assert_eq!(detail.required, 33);
        }
        other => panic!("allocation failures stay allocation refusals, got {other:?}"),
    }
}

/// Every fs-exec budget/cancellation refusal kind maps to the Cancelled
/// outcome with its exact stable phase and numeric context.
#[test]
fn every_cancellation_refusal_kind_is_retained() {
    let refusals = [
        BudgetRefusal::Cancelled { phase: "extract" },
        BudgetRefusal::DeadlineExpiredAtAdmission {
            deadline_ns: 10,
            observed_ns: 11,
        },
        BudgetRefusal::DeadlineWithoutClock { deadline_ns: 12 },
        BudgetRefusal::CostPlanExceedsQuota {
            planned: 13,
            quota: 14,
        },
        BudgetRefusal::DeadlineExpired {
            phase: "sample",
            deadline_ns: 15,
            observed_ns: 16,
        },
        BudgetRefusal::PollsExhausted {
            phase: "seal",
            quota: 17,
        },
        BudgetRefusal::CostExhausted {
            phase: "identity",
            requested: 18,
            remaining: 19,
            quota: 20,
        },
    ];
    let kinds = [
        CancellationKind::CancelledAtCheckpoint,
        CancellationKind::DeadlineExpiredAtAdmission,
        CancellationKind::DeadlineWithoutClock,
        CancellationKind::CostPlanExceedsQuota,
        CancellationKind::DeadlineExpiredMidRun,
        CancellationKind::PollsExhausted,
        CancellationKind::CostQuotaExhausted,
    ];
    for (refusal, kind) in refusals.into_iter().zip(kinds) {
        let outcome = SurfaceLocalization::from(IsoContourError::ExecutionBudgetRefused {
            refusal: refusal.clone(),
        });
        assert_eq!(outcome.status(), SurfaceLocalizationStatus::Cancelled, "{refusal}");
        match &outcome {
            SurfaceLocalization::Cancelled(detail) => {
                assert_eq!(detail.kind, kind, "{refusal}");
                assert_eq!(detail.stage, LocalizationStage::IsoContourExtraction);
            }
            other => panic!("budget refusals are cancellations here, got {other:?}"),
        }
    }
}

/// The default campaign localizes: the typed record agrees bit-for-bit with
/// the derived legacy views.
#[test]
fn localized_campaign_record_agrees_with_derived_views() {
    let net = blob_sdf_net();
    let report = run_campaign(&net, 2.5, 0.3);
    match report.surface_localization {
        SurfaceLocalization::Localized {
            crossings,
            max_radius,
            nearest_radius,
        } => {
            assert!(crossings > 0);
            assert_eq!(crossings, report.surface_crossings);
            assert_eq!(max_radius.to_bits(), report.max_crossing_radius.to_bits());
            assert_eq!(
                nearest_radius.to_bits(),
                report.nearest_surface_radius.to_bits()
            );
            assert!(report.surface_localization.stage().is_none());
        }
        other => panic!("default campaign must localize, got {other:?}"),
    }
}

/// A field that is identically on the level makes every grid edge a
/// coincident segment: the campaign must report the typed Unrepresentable
/// outcome naming the exact contour diagnostic, with NaN legacy radii.
#[test]
fn flat_zero_field_reports_coincident_edge_not_silent_none() {
    let flat = MlpSdf::new(
        vec![
            Layer::new(vec![vec![0.0, 0.0]], vec![0.0]),
            Layer::new(vec![vec![0.0]], vec![0.0]),
        ],
        1.0,
    );
    let report = try_run_campaign(&flat, 2.5, 0.3).expect("flat net admits");
    match report.surface_localization {
        SurfaceLocalization::Unrepresentable(detail) => {
            assert_eq!(detail.stage, LocalizationStage::IsoContourExtraction);
            assert_eq!(
                detail.diagnostic,
                LocalizationDiagnostic::IsoCoincidentLevelEdge
            );
        }
        other => panic!("flat zero field must refuse as coincident, got {other:?}"),
    }
    assert_eq!(report.surface_crossings, 0);
    assert!(report.max_crossing_radius.is_nan());
    assert!(report.nearest_surface_radius.is_nan());
}

/// A strictly positive field is a VALID EMPTY grid: no crossings, finite
/// sentinels (`0` / `+inf`), and no refusal anywhere.
#[test]
fn all_positive_field_reports_valid_empty() {
    let positive = MlpSdf::new(
        vec![
            Layer::new(vec![vec![0.0, 0.0]], vec![0.0]),
            Layer::new(vec![vec![0.0]], vec![20.0]),
        ],
        1.0,
    );
    let report = try_run_campaign(&positive, 2.5, 0.3).expect("positive net admits");
    assert_eq!(report.surface_localization, SurfaceLocalization::ValidEmpty);
    assert_eq!(report.surface_crossings, 0);
    assert_eq!(report.max_crossing_radius.to_bits(), 0.0f64.to_bits());
    assert!(report.nearest_surface_radius.is_infinite());
    assert!(report.nearest_surface_radius > 0.0);
    assert!(report.surface_localization.stage().is_none());
}

/// A network whose evaluation produces NaN must surface the exact first
/// offending node instead of collapsing to the same `None` as valid-empty.
#[test]
fn non_finite_samples_report_the_first_offender() {
    let poisoned = MlpSdf::new(
        vec![
            Layer::new(vec![vec![f64::NAN, 0.0]], vec![0.0]),
            Layer::new(vec![vec![1.0]], vec![0.0]),
        ],
        1.0,
    );
    let report = try_run_campaign(&poisoned, 2.5, 0.3).expect("poisoned net admits");
    match report.surface_localization {
        SurfaceLocalization::Unrepresentable(detail) => {
            assert_eq!(detail.stage, LocalizationStage::GridConstruction);
            assert_eq!(detail.diagnostic, LocalizationDiagnostic::GridNonFiniteValue);
            assert_eq!(detail.first_index, 0, "first sampled node offends");
            assert!(f64::from_bits(detail.scalar_bits).is_nan());
        }
        other => panic!("NaN samples must be unrepresentable, got {other:?}"),
    }
    assert!(report.max_crossing_radius.is_nan());
}
