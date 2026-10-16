//! G0/G3/G4 coverage for the tolerance-aware import census and promotion policy
//! (`frankensim-extreal-program-f85xj.11.3`).

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{
    ImportCensusPolicy, ImportPromotionError, ImportPromotionPolicy, ImportRefusalThresholds,
    IntersectionInspection, census_with_policy, promote_with_policy,
};
use fs_rep_mesh::Soup;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xc3_11_03,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn exhaustive_policy(tolerance: f64, pair_budget: usize) -> ImportCensusPolicy {
    ImportCensusPolicy::try_new(
        tolerance,
        IntersectionInspection::ExhaustiveF64 {
            max_pair_tests: pair_budget,
        },
        1,
    )
    .expect("fixture policy is valid")
}

fn tetra(short_edge: f64) -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(short_edge, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ],
        triangles: vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
    }
}

fn census(soup: &Soup, policy: ImportCensusPolicy) -> fs_io::ImportCensusReport {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        census_with_policy(soup, policy, cx).expect("fixture census succeeds")
    })
}

#[test]
fn g0_clean_mesh_has_no_extended_findings() {
    let report = census(&tetra(1.0), exhaustive_policy(1.0e-6, 16));
    for class in [
        "small-edge",
        "sliver-face",
        "near-boundary-loop-gap",
        "shell-overlap-or-self-intersection",
    ] {
        assert_eq!(report.count(class), 0, "unexpected {class}: {report:?}");
    }
    assert!(report.intersection.complete);
    assert_eq!(report.geometry_budget.tolerance_sensitive_residuals, 0);
    assert_eq!(
        report.geometry_budget.authority,
        "diagnostic-input-not-a-spatial-error-bound"
    );
}

#[test]
fn g0_small_edge_and_sliver_are_relative_to_declared_tolerance() {
    let soup = Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.001, 0.0, 0.0),
            Point3::new(1.0, 0.000_1, 0.0),
        ],
        triangles: vec![[0, 1, 2]],
    };
    let report = census(&soup, exhaustive_policy(0.002, 1));
    assert_eq!(report.count("small-edge"), 1);
    assert_eq!(report.count("sliver-face"), 1);
    assert_eq!(report.geometry_budget.tolerance_sensitive_residuals, 2);
    assert!(report.smallest_edge.is_some_and(|edge| edge <= 0.001));
    assert!(
        report
            .smallest_altitude
            .is_some_and(|altitude| altitude <= 0.002)
    );
}

#[test]
fn g0_near_boundary_loops_report_one_gap() {
    let soup = Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.001, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2], [3, 4, 5]],
    };
    let report = census(&soup, exhaustive_policy(0.01, 1));
    assert_eq!(report.boundary_loops, 2);
    assert_eq!(report.count("near-boundary-loop-gap"), 1);
    assert!(
        report
            .largest_detected_boundary_gap
            .is_some_and(|gap| (gap - 0.001).abs() < 1.0e-12)
    );
}

#[test]
fn g0_crossing_and_coplanar_overlapping_faces_are_detected() {
    let crossing = Soup {
        positions: vec![
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -0.5, -1.0),
            Point3::new(0.0, 0.5, 1.0),
            Point3::new(0.0, 1.0, -1.0),
        ],
        triangles: vec![[0, 1, 2], [3, 4, 5]],
    };
    let crossing_report = census(&crossing, exhaustive_policy(1.0e-6, 1));
    assert_eq!(
        crossing_report.count("shell-overlap-or-self-intersection"),
        1
    );

    let overlap = Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(0.25, 0.25, 0.0),
            Point3::new(1.25, 0.25, 0.0),
            Point3::new(0.25, 1.25, 0.0),
        ],
        triangles: vec![[0, 1, 2], [3, 4, 5]],
    };
    let overlap_report = census(&overlap, exhaustive_policy(1.0e-6, 1));
    assert_eq!(
        overlap_report.count("shell-overlap-or-self-intersection"),
        1
    );
    assert_eq!(
        overlap_report.intersection.authority,
        "f64-geometric-filter-no-exact-predicate-certificate"
    );
}

#[test]
fn g3_sample_level_and_coverage_are_explicit_and_deterministic() {
    let policy = ImportCensusPolicy::try_new(
        1.0e-6,
        IntersectionInspection::DeterministicSampleF64 { sample_count: 1 },
        1,
    )
    .expect("sample policy");
    let first = census(&tetra(1.0), policy);
    let second = census(&tetra(1.0), policy);
    assert_eq!(first, second);
    assert!(!first.intersection.complete);
    assert_eq!(
        first.intersection.requested_level,
        "deterministic-even-sample-f64-filter"
    );
    assert!(first.to_json().contains("\"complete\":false"));
}

#[test]
fn g0_exhaustive_pair_budget_caps_raw_visits_including_adjacency_skips() {
    let mut soup = tetra(1.0);
    soup.triangles = vec![[0, 1, 2]; 128];
    let report = census(&soup, exhaustive_policy(1.0e-6, 3));
    assert_eq!(report.intersection.visited_raw_pairs, 3);
    assert_eq!(report.intersection.inspected_pairs, 0);
    assert_eq!(report.intersection.shared_vertex_pairs_skipped, 3);
    assert!(!report.intersection.complete);
}

#[test]
fn g4_census_observes_pre_requested_cancellation() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    let refusal = with_cx(&gate, |cx| {
        census_with_policy(&tetra(1.0), exhaustive_policy(1.0e-6, 16), cx)
            .expect_err("pre-requested cancellation must refuse")
    });
    assert!(refusal.cancelled);
    assert_eq!(refusal.stage, "census-start");
}

#[test]
fn g0_invalid_structure_is_receipted_and_refused_before_repair() {
    let soup = Soup {
        positions: vec![
            Point3::new(f64::NAN, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 7]],
    };
    let policy = ImportPromotionPolicy {
        profile: "hostile-\"profile\nv1",
        max_hole_edges: 0,
        census: exhaustive_policy(1.0e-6, 1),
        thresholds: ImportRefusalThresholds::validation_grade(),
    };
    let gate = CancelGate::new_clock_free();
    let error = with_cx(&gate, |cx| {
        promote_with_policy(
            fs_io::quarantine::quarantine(soup, "obj\"", b"invalid-fixture"),
            policy,
            cx,
        )
        .expect_err("unsafe structure must refuse before repair indexes it")
    });
    let ImportPromotionError::Refused(refusal) = error else {
        panic!("expected structural refusal")
    };
    assert!(refusal.receipt_json.contains("\"invalid-face-index\""));
    assert!(refusal.receipt_json.contains("\"non-finite-vertex\""));
    assert!(refusal.receipt_json.contains("\"format\":\"obj\\\"\""));
    assert!(
        refusal
            .receipt_json
            .contains("\"profile\":\"hostile-\\\"profile\\nv1\"")
    );
    assert!(refusal.receipt_json.contains("\"trust\":\"refused\""));
}

#[test]
fn g0_residual_slivers_can_promote_only_under_receipted_thresholds() {
    let soup = tetra(0.000_1);
    let mut scoping_thresholds = ImportRefusalThresholds::validation_grade();
    scoping_thresholds.small_edges = 8;
    scoping_thresholds.sliver_faces = 8;
    let scoping = ImportPromotionPolicy {
        profile: "supplier-scoping-v1",
        max_hole_edges: 0,
        census: exhaustive_policy(0.001, 16),
        thresholds: scoping_thresholds,
    };
    let gate = CancelGate::new_clock_free();
    let (evidence, receipt) = with_cx(&gate, |cx| {
        promote_with_policy(
            fs_io::quarantine::quarantine(soup.clone(), "obj", b"sliver-fixture"),
            scoping,
            cx,
        )
        .expect("scoping policy accepts receipted residual slivers")
    });
    assert_eq!(evidence.value.triangles.len(), 4);
    assert!(receipt.after.count("sliver-face") > 0);
    let json = receipt.to_json();
    assert!(json.contains("\"profile\":\"supplier-scoping-v1\""));
    assert!(json.contains("\"sliver-face\""));
    assert!(json.contains("\"tolerance_sensitive_residuals\":"));
    assert!(json.contains("\"trust\":\"promoted\""));

    let validation = ImportPromotionPolicy {
        profile: "supplier-validation-v1",
        thresholds: ImportRefusalThresholds::validation_grade(),
        ..scoping
    };
    let gate = CancelGate::new_clock_free();
    let refusal = with_cx(&gate, |cx| {
        promote_with_policy(
            fs_io::quarantine::quarantine(soup, "obj", b"sliver-fixture"),
            validation,
            cx,
        )
        .expect_err("validation policy refuses the same residuals")
    });
    let ImportPromotionError::Refused(refusal) = refusal else {
        panic!("expected threshold refusal")
    };
    assert!(
        refusal
            .blocking
            .iter()
            .any(|blocker| blocker.contains("sliver-face"))
    );
    assert!(
        refusal
            .receipt_json
            .contains("\"profile\":\"supplier-validation-v1\"")
    );
    assert!(refusal.receipt_json.contains("\"trust\":\"refused\""));
}

#[test]
fn g0_receipt_carries_repair_operations_and_class_deltas() {
    let mut dirty = tetra(1.0);
    dirty.triangles.push(dirty.triangles[0]);
    let policy = ImportPromotionPolicy {
        profile: "repair-history-v1",
        max_hole_edges: 0,
        census: exhaustive_policy(1.0e-6, 32),
        thresholds: ImportRefusalThresholds::validation_grade(),
    };
    let gate = CancelGate::new_clock_free();
    let (_, receipt) = with_cx(&gate, |cx| {
        promote_with_policy(
            fs_io::quarantine::quarantine(dirty, "obj", b"duplicate-fixture"),
            policy,
            cx,
        )
        .expect("duplicate repair promotes")
    });
    assert_eq!(receipt.before.count("duplicate-face"), 1);
    assert_eq!(receipt.after.count("duplicate-face"), 0);
    let json = receipt.to_json();
    assert!(json.contains("\"repair_history\""));
    assert!(json.contains("\"defect\":\"duplicate-face\""));
    assert!(json.contains("\"class\":\"duplicate-face\",\"before\":1,\"after\":0"));
    assert!(json.contains("\"tolerance_based_repair_performed\":false"));
}

#[test]
fn g0_complete_intersection_requirement_refuses_a_sampled_census() {
    let policy = ImportPromotionPolicy {
        profile: "sampled-but-validation-v1",
        max_hole_edges: 0,
        census: ImportCensusPolicy::try_new(
            1.0e-6,
            IntersectionInspection::DeterministicSampleF64 { sample_count: 1 },
            1,
        )
        .expect("sample policy"),
        thresholds: ImportRefusalThresholds::validation_grade(),
    };
    let gate = CancelGate::new_clock_free();
    let error = with_cx(&gate, |cx| {
        promote_with_policy(
            fs_io::quarantine::quarantine(tetra(1.0), "obj", b"sample-fixture"),
            policy,
            cx,
        )
        .expect_err("validation requires complete intersection coverage")
    });
    let ImportPromotionError::Refused(refusal) = error else {
        panic!("expected policy refusal")
    };
    assert!(
        refusal
            .blocking
            .iter()
            .any(|blocker| blocker.contains("intersection-census-incomplete"))
    );
}
