//! G0/G3/G4 battery for selector-anchored revision association.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{
    AssignmentLimits, AssignmentReport, AssignmentRequest, AssociationPolicy, AssociationVerdict,
    MeshSelector, MigrationAction, RegionChange, RigidTransform3, associate_mesh_assignments,
    resolve_mesh_assignments,
};
use fs_rep_mesh::Soup;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xa5_50_c1_a7,
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

fn square() -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2], [0, 2, 3]],
    }
}

fn refined_square() -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
        ],
        triangles: vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    }
}

fn two_squares() -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 1.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
        ],
        triangles: vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
    }
}

fn box_selector() -> MeshSelector {
    MeshSelector::Box {
        min: [-0.1, -0.1, -0.1],
        max: [1.1, 1.1, 0.1],
        tolerance: 0.0,
    }
}

fn explicit(faces: &[u32]) -> MeshSelector {
    MeshSelector::ExplicitFaceSet {
        faces: faces.to_vec(),
        fragility_acknowledged: true,
    }
}

fn request(subject: &str, selector: MeshSelector, allow_overlap: bool) -> AssignmentRequest {
    AssignmentRequest {
        subject: subject.to_string(),
        selector,
        allow_overlap,
    }
}

fn resolve(
    soup: &Soup,
    identity: &str,
    unit: &str,
    requests: &[AssignmentRequest],
    cx: &Cx<'_>,
) -> AssignmentReport {
    resolve_mesh_assignments(
        soup,
        identity,
        unit,
        &[],
        requests,
        AssignmentLimits::DEFAULT,
        cx,
    )
    .expect("assignment fixture resolves")
}

fn policy() -> AssociationPolicy {
    AssociationPolicy {
        stable_relative_area: 1.0e-12,
        stable_distance: 1.0e-12,
        stable_orientation: 1.0e-12,
        stable_relative_extent: 1.0e-12,
        moved_relative_area: 0.5,
        moved_distance: 1.0,
        moved_orientation: 0.5,
        moved_relative_extent: 0.5,
        ambiguity_score_gap: 1.0e-12,
        ..AssociationPolicy::engineering(1.0e-12, 1.0)
    }
}

#[test]
fn g3_retessellation_is_stable_without_face_count_identity() {
    let source_soup = square();
    let target_soup = refined_square();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let source = resolve(
            &source_soup,
            "blake3:panel-r1",
            "m",
            &[request("surface:panel", box_selector(), false)],
            cx,
        );
        let target = resolve(
            &target_soup,
            "blake3:panel-r2",
            "m",
            &[request("surface:panel", box_selector(), false)],
            cx,
        );
        let report = associate_mesh_assignments(
            &source_soup,
            &source,
            &target_soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("retessellation associates");

        let decision = &report.decisions[0];
        assert_eq!(decision.verdict, AssociationVerdict::Stable);
        assert_eq!(decision.change, RegionChange::Unchanged);
        assert_eq!(decision.migration, MigrationAction::AutoApply);
        assert_eq!(decision.source_fingerprint.face_count, 2);
        assert_eq!(
            decision
                .target_fingerprint
                .as_ref()
                .expect("target fingerprint")
                .face_count,
            4
        );
        assert_eq!(
            decision.source_fingerprint.topology,
            decision
                .target_fingerprint
                .as_ref()
                .expect("target fingerprint")
                .topology
        );
        assert!(report.added.is_empty());
        let json = report.to_json();
        assert!(json.contains("\"association\":\"stable\""));
        assert!(json.contains("\"authority\":\"finite-tessellation-association-diagnostic\""));
        assert!(json.contains("\"no_claim\":"));
        assert!(!json.contains('\n'));
        let markdown = report.render_markdown();
        assert!(
            markdown
                .contains("| surface:panel | surface:panel | stable | unchanged | auto-apply |")
        );
    });
}

#[test]
fn g3_declared_rigid_frame_transform_preserves_stability_and_is_receipted() {
    let source_soup = square();
    let mut target_soup = square();
    for point in &mut target_soup.positions {
        let source = [point.x, point.y, point.z];
        point.x = -source[1] + 3.0;
        point.y = source[0] - 2.0;
        point.z = source[2] + 5.0;
    }
    let transform = RigidTransform3::new(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [3.0, -2.0, 5.0],
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let requests = [request("surface:panel", explicit(&[0, 1]), false)];
        let source = resolve(&source_soup, "blake3:panel-r1", "m", &requests, cx);
        let target = resolve(&target_soup, "blake3:panel-r2", "m", &requests, cx);
        let report = associate_mesh_assignments(
            &source_soup,
            &source,
            &target_soup,
            &target,
            transform,
            policy(),
            cx,
        )
        .expect("declared rigid transform aligns");

        assert_eq!(report.decisions[0].verdict, AssociationVerdict::Stable);
        assert_eq!(report.receipt.source_to_target(), transform);
        assert!(report.to_json().contains("\"translation\":[3,-2,5]"));
    });
}

#[test]
fn g0_exact_subject_distinguishes_motion_deformation_and_topology_change() {
    let source_soup = square();
    let mut moved_soup = square();
    for point in &mut moved_soup.positions {
        point.z += 0.2;
    }
    let mut deformed_soup = square();
    for point in &mut deformed_soup.positions {
        point.x *= 1.2;
    }
    let topology_soup = two_squares();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let source_requests = [request("surface:panel", explicit(&[0, 1]), false)];
        let source = resolve(&source_soup, "blake3:panel-r1", "m", &source_requests, cx);

        let moved = resolve(&moved_soup, "blake3:panel-moved", "m", &source_requests, cx);
        let moved_report = associate_mesh_assignments(
            &source_soup,
            &source,
            &moved_soup,
            &moved,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("motion classified");
        assert_eq!(
            (
                moved_report.decisions[0].verdict,
                moved_report.decisions[0].change,
                moved_report.decisions[0].migration,
            ),
            (
                AssociationVerdict::Moved,
                RegionChange::Moved,
                MigrationAction::Propose,
            )
        );

        let deformed = resolve(
            &deformed_soup,
            "blake3:panel-deformed",
            "m",
            &source_requests,
            cx,
        );
        let deformed_report = associate_mesh_assignments(
            &source_soup,
            &source,
            &deformed_soup,
            &deformed,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("deformation classified");
        assert_eq!(deformed_report.decisions[0].change, RegionChange::Deformed);
        assert_eq!(
            deformed_report.decisions[0].migration,
            MigrationAction::Propose
        );

        let topology_requests = [request("surface:panel", explicit(&[0, 1, 2, 3]), false)];
        let topology = resolve(
            &topology_soup,
            "blake3:panel-topology",
            "m",
            &topology_requests,
            cx,
        );
        let topology_report = associate_mesh_assignments(
            &source_soup,
            &source,
            &topology_soup,
            &topology,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("topology change is a report, not a silent rebind");
        assert_eq!(
            (
                topology_report.decisions[0].verdict,
                topology_report.decisions[0].change,
                topology_report.decisions[0].migration,
            ),
            (
                AssociationVerdict::Lost,
                RegionChange::TopologyChanged,
                MigrationAction::Refuse,
            )
        );
    });
}

#[test]
fn g0_missing_subject_with_two_equal_candidates_is_ambiguous_not_guessed() {
    let source_soup = square();
    let target_soup = square();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let source = resolve(
            &source_soup,
            "blake3:panel-r1",
            "m",
            &[request("surface:old-panel", box_selector(), false)],
            cx,
        );
        let target = resolve(
            &target_soup,
            "blake3:panel-r2",
            "m",
            &[
                request("surface:candidate-b", box_selector(), true),
                request("surface:candidate-a", box_selector(), true),
            ],
            cx,
        );
        let report = associate_mesh_assignments(
            &source_soup,
            &source,
            &target_soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("ambiguity is typed output");

        let decision = &report.decisions[0];
        assert_eq!(decision.verdict, AssociationVerdict::Ambiguous);
        assert_eq!(decision.change, RegionChange::Ambiguous);
        assert_eq!(decision.migration, MigrationAction::Refuse);
        assert_eq!(
            decision.candidates,
            ["surface:candidate-a", "surface:candidate-b"]
        );
        assert_eq!(report.added.len(), 2);
    });
}

#[test]
fn g0_unique_selector_anchored_fallback_is_only_a_proposal() {
    let source_soup = square();
    let target_soup = square();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let source = resolve(
            &source_soup,
            "blake3:panel-r1",
            "m",
            &[request("surface:old-panel", box_selector(), false)],
            cx,
        );
        let target = resolve(
            &target_soup,
            "blake3:panel-r2",
            "m",
            &[request("surface:renamed-panel", box_selector(), false)],
            cx,
        );
        let report = associate_mesh_assignments(
            &source_soup,
            &source,
            &target_soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("unique fallback resolves");

        assert_eq!(
            report.decisions[0].target_subject.as_deref(),
            Some("surface:renamed-panel")
        );
        assert_eq!(report.decisions[0].verdict, AssociationVerdict::Moved);
        assert_eq!(
            report.decisions[0].migration,
            MigrationAction::Propose,
            "renaming can never auto-apply merely because geometry matches"
        );
        assert!(report.added.is_empty());
    });
}

#[test]
fn g0_receipt_binding_and_one_to_one_fallback_prevent_authority_laundering() {
    let soup = square();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let source = resolve(
            &soup,
            "blake3:panel-r1",
            "m",
            &[request("surface:panel", box_selector(), false)],
            cx,
        );
        let target = resolve(
            &soup,
            "blake3:panel-r2",
            "m",
            &[request("surface:panel", box_selector(), false)],
            cx,
        );

        let mut wrong_soup = soup.clone();
        wrong_soup.positions[0].x = 0.125;
        let soup_error = associate_mesh_assignments(
            &wrong_soup,
            &source,
            &soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect_err("index-compatible different soup cannot borrow a receipt");
        assert_eq!(soup_error.code, "mesh-association-soup-receipt-mismatch");

        let mut mutated_report = source.clone();
        mutated_report.assignments[0].subject = "surface:mutated".to_string();
        let report_error = associate_mesh_assignments(
            &soup,
            &mutated_report,
            &soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect_err("mutated assignment rows cannot borrow a receipt");
        assert_eq!(
            report_error.code,
            "mesh-association-assignment-receipt-mismatch"
        );

        let competing_source = resolve(
            &soup,
            "blake3:panel-r1",
            "m",
            &[
                request("surface:old-a", box_selector(), true),
                request("surface:old-b", box_selector(), true),
            ],
            cx,
        );
        let renamed_target = resolve(
            &soup,
            "blake3:panel-r2",
            "m",
            &[request("surface:renamed", box_selector(), false)],
            cx,
        );
        let conflict = associate_mesh_assignments(
            &soup,
            &competing_source,
            &soup,
            &renamed_target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect("many-to-one fallback becomes typed ambiguity");
        assert!(
            conflict
                .decisions
                .iter()
                .all(|decision| decision.verdict == AssociationVerdict::Ambiguous)
        );
        assert!(
            conflict
                .decisions
                .iter()
                .all(|decision| decision.migration == MigrationAction::Refuse)
        );
        assert_eq!(conflict.added[0].subject, "surface:renamed");
    });
}

#[test]
fn g4_invalid_frame_units_work_cap_and_cancellation_refuse_atomically() {
    let soup = square();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let requests = [request("surface:panel", explicit(&[0, 1]), false)];
        let source = resolve(&soup, "blake3:panel-r1", "m", &requests, cx);
        let millimetres = resolve(&soup, "blake3:panel-r2", "mm", &requests, cx);
        let unit_error = associate_mesh_assignments(
            &soup,
            &source,
            &soup,
            &millimetres,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect_err("unit mismatch refuses");
        assert_eq!(unit_error.code, "mesh-association-unit-mismatch");

        let target = resolve(&soup, "blake3:panel-r2", "m", &requests, cx);
        let reflected = RigidTransform3::new(
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [0.0; 3],
        );
        let frame_error =
            associate_mesh_assignments(&soup, &source, &soup, &target, reflected, policy(), cx)
                .expect_err("reflection is not a proper rigid frame");
        assert_eq!(frame_error.code, "mesh-association-invalid-frame");

        let work_error = associate_mesh_assignments(
            &soup,
            &source,
            &soup,
            &target,
            RigidTransform3::IDENTITY,
            AssociationPolicy {
                max_face_references: 3,
                ..policy()
            },
            cx,
        )
        .expect_err("face-reference cap refuses before fingerprinting");
        assert_eq!(work_error.code, "mesh-association-resource-bound");
    });

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let live_gate = CancelGate::new_clock_free();
        let (source, target) = with_cx(&live_gate, |live_cx| {
            let requests = [request("surface:panel", explicit(&[0, 1]), false)];
            (
                resolve(&soup, "blake3:panel-r1", "m", &requests, live_cx),
                resolve(&soup, "blake3:panel-r2", "m", &requests, live_cx),
            )
        });
        let error = associate_mesh_assignments(
            &soup,
            &source,
            &soup,
            &target,
            RigidTransform3::IDENTITY,
            policy(),
            cx,
        )
        .expect_err("pre-requested cancellation refuses");
        assert_eq!(error.code, "mesh-association-cancelled");
    });
}
