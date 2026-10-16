//! G0/G3/G4 coverage for mesh-index-free assignment selectors.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{
    AssignmentLimits, AssignmentRequest, HalfSpaceSide, MeshSelector, NamedFaceGroup,
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
                seed: 0xa5_51_6e,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn cube() -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        triangles: vec![
            [4, 5, 6],
            [4, 6, 7],
            [0, 2, 1],
            [0, 3, 2],
            [0, 1, 5],
            [0, 5, 4],
            [3, 7, 6],
            [3, 6, 2],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
        ],
    }
}

fn retessellated_cube() -> Soup {
    let mut soup = cube();
    soup.positions.push(Point3::new(0.5, 0.5, 1.0));
    soup.triangles
        .splice(0..2, [[4, 5, 8], [5, 6, 8], [6, 7, 8], [7, 4, 8]]);
    soup
}

fn request(subject: &str, selector: MeshSelector) -> AssignmentRequest {
    AssignmentRequest {
        subject: subject.to_string(),
        selector,
        allow_overlap: false,
    }
}

fn top_half_space() -> MeshSelector {
    MeshSelector::HalfSpace {
        normal: [0.0, 0.0, 1.0],
        offset: 1.0,
        side: HalfSpaceSide::AtLeast,
        tolerance: 0.0,
    }
}

fn resolve_one(soup: &Soup, selector: MeshSelector) -> fs_io::AssignmentReport {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        resolve_mesh_assignments(
            soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request("surface:top", selector)],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect("selector resolves")
    })
}

#[test]
fn g0_named_geometric_and_acknowledged_selectors_resolve_with_stats() {
    let soup = cube();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let named = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[NamedFaceGroup {
                name: "TOP".to_string(),
                faces: vec![0, 1],
            }],
            &[request(
                "surface:top",
                MeshSelector::NamedGroup {
                    name: "TOP".to_string(),
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect("named group resolves");
        assert_eq!(named.assignments[0].faces, vec![0, 1]);

        for selector in [
            top_half_space(),
            MeshSelector::Box {
                min: [0.0, 0.0, 1.0],
                max: [1.0, 1.0, 1.0],
                tolerance: 0.0,
            },
            MeshSelector::Cylinder {
                origin: [0.5, 0.5, 0.0],
                axis: [0.0, 0.0, 2.0],
                radius: 0.8,
                axial_min: 1.0,
                axial_max: 1.0,
                tolerance: 0.0,
            },
            MeshSelector::NearestDatum {
                point: [0.5, 0.5, 1.0],
                max_distance: 0.0,
                tolerance: 0.0,
            },
        ] {
            let report = resolve_mesh_assignments(
                &soup,
                "blake3:fixture-cube",
                "m",
                &[],
                &[request("surface:top", selector)],
                AssignmentLimits::DEFAULT,
                cx,
            )
            .expect("geometric selector resolves");
            let assignment = &report.assignments[0];
            assert_eq!(assignment.faces, vec![0, 1]);
            assert_eq!(assignment.stats.surface_area, 1.0);
            assert_eq!(assignment.stats.enclosed_volume, None);
            assert_eq!(assignment.stats.bounds_min, [0.0, 0.0, 1.0]);
            assert_eq!(assignment.stats.bounds_max, [1.0, 1.0, 1.0]);
        }

        let all_faces: Vec<u32> = (0..12).collect();
        let volume = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "region:enclosure",
                MeshSelector::ExplicitFaceSet {
                    faces: all_faces,
                    fragility_acknowledged: true,
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect("acknowledged full boundary resolves");
        assert_eq!(volume.assignments[0].stats.surface_area, 6.0);
        let enclosed_volume = volume.assignments[0]
            .stats
            .enclosed_volume
            .expect("closed oriented cube reports volume");
        assert!((enclosed_volume - 1.0).abs() <= 8.0 * f64::EPSILON);
        let receipt = volume.to_json();
        assert!(receipt.contains("\"authority\":\"finite-tessellation-selection\""));
        assert!(receipt.contains("\"enclosed_volume\":"));
        assert!(!receipt.contains("\"enclosed_volume\":null"));
        assert!(!receipt.contains('\n'));
    });
}

#[test]
fn g3_rigid_translation_and_retessellation_preserve_geometric_assignment() {
    let base = resolve_one(&cube(), top_half_space());
    let refined = resolve_one(&retessellated_cube(), top_half_space());
    assert_eq!(base.assignments[0].subject, refined.assignments[0].subject);
    assert_eq!(base.assignments[0].stats.surface_area, 1.0);
    assert_eq!(refined.assignments[0].stats.surface_area, 1.0);
    assert_eq!(base.assignments[0].stats.face_count, 2);
    assert_eq!(refined.assignments[0].stats.face_count, 4);
    assert_eq!(
        base.receipt.requests_fingerprint(),
        refined.receipt.requests_fingerprint(),
        "same identity and selector remain the same request"
    );
    assert_ne!(
        base.receipt.source_mesh_fingerprint(),
        refined.receipt.source_mesh_fingerprint(),
        "retessellation remains visible in provenance"
    );

    let mut translated = cube();
    for point in &mut translated.positions {
        point.x += 2.0;
        point.y -= 3.0;
        point.z += 5.0;
    }
    let moved = resolve_one(
        &translated,
        MeshSelector::HalfSpace {
            normal: [0.0, 0.0, 1.0],
            offset: 6.0,
            side: HalfSpaceSide::AtLeast,
            tolerance: 0.0,
        },
    );
    assert_eq!(moved.assignments[0].stats.surface_area, 1.0);
    assert_eq!(moved.assignments[0].stats.bounds_min, [2.0, -3.0, 6.0]);
    assert_eq!(moved.assignments[0].stats.bounds_max, [3.0, -2.0, 6.0]);
}

#[test]
fn g0_overlap_empty_group_and_fragile_face_sets_refuse_with_fixes() {
    let soup = cube();
    let groups = [NamedFaceGroup {
        name: "TOP".to_string(),
        faces: vec![0, 1],
    }];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let overlap = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &groups,
            &[
                request(
                    "surface:named-top",
                    MeshSelector::NamedGroup {
                        name: "TOP".to_string(),
                    },
                ),
                request("surface:geometric-top", top_half_space()),
            ],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("undeclared overlap refuses");
        assert_eq!(overlap.code, "mesh-assignment-overlap");
        assert!(!overlap.fix.is_empty());

        let mut allowed = [
            request(
                "surface:named-top",
                MeshSelector::NamedGroup {
                    name: "TOP".to_string(),
                },
            ),
            request("surface:geometric-top", top_half_space()),
        ];
        for row in &mut allowed {
            row.allow_overlap = true;
        }
        assert_eq!(
            resolve_mesh_assignments(
                &soup,
                "blake3:fixture-cube",
                "m",
                &groups,
                &allowed,
                AssignmentLimits::DEFAULT,
                cx,
            )
            .expect("explicit overlap resolves")
            .assignments
            .len(),
            2
        );

        let fragile = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "surface:fragile",
                MeshSelector::ExplicitFaceSet {
                    faces: vec![0],
                    fragility_acknowledged: false,
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("unacknowledged mesh ordinals refuse");
        assert_eq!(fragile.code, "mesh-assignment-fragility-unacknowledged");

        let empty = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "surface:absent",
                MeshSelector::Box {
                    min: [10.0, 10.0, 10.0],
                    max: [11.0, 11.0, 11.0],
                    tolerance: 0.0,
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("empty selector refuses");
        assert_eq!(empty.code, "mesh-assignment-empty-selection");

        let overflow = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "surface:overflow",
                MeshSelector::HalfSpace {
                    normal: [1.0, 0.0, 0.0],
                    offset: f64::MAX,
                    side: HalfSpaceSide::AtMost,
                    tolerance: f64::MAX,
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("overflowing admitted threshold refuses");
        assert_eq!(overflow.code, "mesh-assignment-invalid-selector");

        let unnormalizable = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "surface:unnormalizable",
                MeshSelector::Cylinder {
                    origin: [0.0; 3],
                    axis: [f64::MAX, f64::MAX, 0.0],
                    radius: 1.0,
                    axial_min: 0.0,
                    axial_max: 1.0,
                    tolerance: 0.0,
                },
            )],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("non-normalizable finite axis refuses");
        assert_eq!(unnormalizable.code, "mesh-assignment-invalid-selector");

        let mut tight_limits = AssignmentLimits::DEFAULT;
        tight_limits.max_selected_faces = 1;
        let oversized_explicit = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request(
                "surface:oversized-explicit",
                MeshSelector::ExplicitFaceSet {
                    faces: vec![0, 1],
                    fragility_acknowledged: true,
                },
            )],
            tight_limits,
            cx,
        )
        .expect_err("explicit selectors respect the publication cap");
        assert_eq!(oversized_explicit.code, "mesh-assignment-resource-limit");
    });
}

#[test]
fn g4_cancellation_and_work_caps_refuse_before_publication() {
    let soup = cube();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let error = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request("surface:top", top_half_space())],
            AssignmentLimits::DEFAULT,
            cx,
        )
        .expect_err("pre-requested cancellation refuses");
        assert_eq!(error.code, "mesh-assignment-cancelled");
    });

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let error = resolve_mesh_assignments(
            &soup,
            "blake3:fixture-cube",
            "m",
            &[],
            &[request("surface:top", top_half_space())],
            AssignmentLimits {
                max_predicate_tests: 11,
                ..AssignmentLimits::DEFAULT
            },
            cx,
        )
        .expect_err("work cap refuses");
        assert_eq!(error.code, "mesh-assignment-work-limit");
        assert!(!error.fix.is_empty());
    });
}
