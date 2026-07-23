//! G0/G3/G4 coverage for the L6 project-to-fs-io assignment adapter.

use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{AssignmentLimits, HalfSpaceSide, MeshSelector, NamedFaceGroup};
use fs_project::assignment::{
    InterfaceAuditLimits, InterfaceDeclarationAudit, audit_interface_declarations,
};
use fs_project::{
    EntityDecl, GEOMETRY_ASSIGNMENT_REPORT_DOMAIN, GeometryArtifact, GeometryAssignment,
    GeometryResolution, ImportedMeshLibrary, ProjectSpec, geometry_source_identity,
    resolve_geometry_assignments,
};
use fs_rep_mesh::Soup;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x6a_03,
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

fn artifact() -> GeometryArtifact {
    GeometryArtifact {
        role: "enclosure".to_string(),
        format: "stl".to_string(),
        source_hash: 0x0123_4567_89ab_cdef,
        parser_version: "fs-io/stl/v1".to_string(),
    }
}

fn assembly() -> Vec<EntityDecl> {
    vec![
        EntityDecl::Assembly {
            name: "rig".to_string(),
            display: "Rig".to_string(),
            expect_id: None,
        },
        EntityDecl::Part {
            parent: "rig".to_string(),
            name: "stack".to_string(),
            display: "Stack".to_string(),
            expect_id: None,
        },
        EntityDecl::Region {
            parent: "stack".to_string(),
            name: "top".to_string(),
            display: "Top region".to_string(),
            expect_id: None,
        },
        EntityDecl::Region {
            parent: "stack".to_string(),
            name: "bottom".to_string(),
            display: "Bottom region".to_string(),
            expect_id: None,
        },
        EntityDecl::Interface {
            parent: "rig".to_string(),
            name: "front-interface".to_string(),
            display: "Front interface".to_string(),
            from: "top".to_string(),
            to: "bottom".to_string(),
            expect_id: None,
        },
    ]
}

fn named_spec() -> ProjectSpec {
    ProjectSpec {
        geometry: Some(vec![artifact()]),
        assignments: Some(vec![
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "top".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "TOP".to_string(),
                },
                allow_overlap: false,
            },
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "bottom".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "BOTTOM".to_string(),
                },
                allow_overlap: false,
            },
            GeometryAssignment {
                artifact: "enclosure".to_string(),
                target: "front-interface".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "FRONT".to_string(),
                },
                allow_overlap: false,
            },
        ]),
        assembly: Some(assembly()),
        ..ProjectSpec::default()
    }
}

fn named_groups() -> Vec<NamedFaceGroup> {
    vec![
        NamedFaceGroup {
            name: "TOP".to_string(),
            faces: vec![0, 1],
        },
        NamedFaceGroup {
            name: "BOTTOM".to_string(),
            faces: vec![2, 3],
        },
        NamedFaceGroup {
            name: "FRONT".to_string(),
            faces: vec![4, 5],
        },
    ]
}

fn spec_without_interface() -> ProjectSpec {
    let mut spec = named_spec();
    spec.assembly
        .as_mut()
        .expect("assembly")
        .retain(|entity| !matches!(entity, EntityDecl::Interface { .. }));
    spec.assignments
        .as_mut()
        .expect("assignments")
        .retain(|assignment| assignment.target != "front-interface");
    spec
}

fn one_top_spec() -> ProjectSpec {
    ProjectSpec {
        geometry: Some(vec![artifact()]),
        assignments: Some(vec![GeometryAssignment {
            artifact: "enclosure".to_string(),
            target: "top".to_string(),
            length_unit: "m".to_string(),
            selector: MeshSelector::HalfSpace {
                normal: [0.0, 0.0, 1.0],
                offset: 1.0,
                side: HalfSpaceSide::AtLeast,
                tolerance: 0.0,
            },
            allow_overlap: false,
        }]),
        assembly: Some(vec![
            EntityDecl::Assembly {
                name: "rig".to_string(),
                display: "Rig".to_string(),
                expect_id: None,
            },
            EntityDecl::Part {
                parent: "rig".to_string(),
                name: "part".to_string(),
                display: "Part".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "part".to_string(),
                name: "top".to_string(),
                display: "Top".to_string(),
                expect_id: None,
            },
        ]),
        ..ProjectSpec::default()
    }
}

#[test]
fn g0_project_names_bind_to_entity_ids_and_exact_lower_reports() {
    let spec = named_spec();
    let artifact = spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut library = ImportedMeshLibrary::new();
    let source_identity = library.insert(&artifact, cube(), "m", named_groups());
    assert_eq!(source_identity, geometry_source_identity(&artifact));

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let resolution =
            resolve_geometry_assignments(&spec, &library, AssignmentLimits::DEFAULT, cx);
        assert!(
            resolution.admissible(),
            "assignment must resolve: {:?}",
            resolution.violations
        );
        assert_eq!(resolution.artifacts.len(), 1);
        let retained = &resolution.artifacts[0];
        assert_eq!(retained.source_identity, source_identity);
        assert_eq!(retained.entities.len(), 3);
        assert_eq!(retained.report.assignments.len(), 3);
        assert_eq!(retained.report_bytes, retained.report.to_json().as_bytes());
        assert_eq!(
            retained.report_hash,
            hash_domain(GEOMETRY_ASSIGNMENT_REPORT_DOMAIN, &retained.report_bytes).to_hex()
        );

        let mut entity_findings = Vec::new();
        let ids = spec.resolve_entities(&mut entity_findings);
        assert!(entity_findings.is_empty());
        for (bound, lower) in retained.entities.iter().zip(&retained.report.assignments) {
            assert_eq!(bound.entity_id, ids[&bound.declared_target]);
            assert_eq!(lower.subject, bound.entity_id.token());
            assert!(lower.stats.face_count > 0);
        }
        let table = resolution.render_table();
        assert!(table.contains("top | entity region:"));
        assert!(table.contains("front-interface | entity interface:"));
        assert!(table.contains(&retained.report_hash));
        assert!(GeometryResolution::no_claim().contains("does not authenticate"));
    });
}

#[test]
fn g3_retessellation_preserves_the_persistent_subject_and_surface_stats() {
    let spec = one_top_spec();
    let artifact = spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut original_library = ImportedMeshLibrary::new();
    original_library.insert(&artifact, cube(), "m", Vec::new());
    let mut retessellated_spec = spec.clone();
    retessellated_spec.geometry.as_mut().expect("geometry")[0].source_hash ^= 1;
    let retessellated_artifact = retessellated_spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut retessellated_library = ImportedMeshLibrary::new();
    retessellated_library.insert(
        &retessellated_artifact,
        retessellated_cube(),
        "m",
        Vec::new(),
    );

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let original =
            resolve_geometry_assignments(&spec, &original_library, AssignmentLimits::DEFAULT, cx);
        let retessellated = resolve_geometry_assignments(
            &retessellated_spec,
            &retessellated_library,
            AssignmentLimits::DEFAULT,
            cx,
        );
        assert!(original.admissible(), "{:?}", original.violations);
        assert!(retessellated.admissible(), "{:?}", retessellated.violations);
        assert_eq!(
            original.artifacts[0].entities[0].entity_id,
            retessellated.artifacts[0].entities[0].entity_id
        );
        assert_eq!(
            original.artifacts[0].report.assignments[0]
                .stats
                .surface_area,
            1.0
        );
        assert_eq!(
            retessellated.artifacts[0].report.assignments[0]
                .stats
                .surface_area,
            1.0
        );
        assert_ne!(
            original.artifacts[0].report.assignments[0].stats.face_count,
            retessellated.artifacts[0].report.assignments[0]
                .stats
                .face_count
        );
    });
}

#[test]
fn g0_dangling_kind_unit_and_lower_selector_refusals_publish_nothing() {
    let base = one_top_spec();
    let artifact = base.geometry.as_ref().expect("geometry")[0].clone();
    let mut library = ImportedMeshLibrary::new();
    library.insert(&artifact, cube(), "m", Vec::new());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut dangling = base.clone();
        dangling.assignments.as_mut().expect("assignments")[0].target = "ghost".to_string();
        let result =
            resolve_geometry_assignments(&dangling, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(
            result.violations[0].code,
            "project-assignment-target-unknown"
        );
        assert!(result.artifacts.is_empty());

        let mut wrong_kind = base.clone();
        wrong_kind.assignments.as_mut().expect("assignments")[0].target = "part".to_string();
        let result =
            resolve_geometry_assignments(&wrong_kind, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(result.violations[0].code, "project-assignment-target-kind");
        assert!(result.artifacts.is_empty());

        let mut wrong_unit = base.clone();
        wrong_unit.assignments.as_mut().expect("assignments")[0].length_unit = "mm".to_string();
        let result =
            resolve_geometry_assignments(&wrong_unit, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(
            result.violations[0].code,
            "project-assignment-unit-mismatch"
        );
        assert!(result.artifacts.is_empty());

        let mut empty = base.clone();
        empty.assignments.as_mut().expect("assignments")[0].selector = MeshSelector::Box {
            min: [10.0, 10.0, 10.0],
            max: [11.0, 11.0, 11.0],
            tolerance: 0.0,
        };
        let result = resolve_geometry_assignments(&empty, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(result.violations[0].code, "mesh-assignment-empty-selection");
        assert!(!result.violations[0].fix.is_empty());
        assert!(result.artifacts.is_empty());

        let mut fragile = base.clone();
        fragile.assignments.as_mut().expect("assignments")[0].selector =
            MeshSelector::ExplicitFaceSet {
                faces: vec![0],
                fragility_acknowledged: false,
            };
        let result =
            resolve_geometry_assignments(&fragile, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(
            result.violations[0].code,
            "mesh-assignment-fragility-unacknowledged"
        );
        assert!(result.artifacts.is_empty());
    });
}

#[test]
fn g4_cancellation_is_atomic_at_the_l6_boundary() {
    let spec = one_top_spec();
    let artifact = spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut library = ImportedMeshLibrary::new();
    library.insert(&artifact, cube(), "m", Vec::new());
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        let result = resolve_geometry_assignments(&spec, &library, AssignmentLimits::DEFAULT, cx);
        assert_eq!(result.violations[0].code, "mesh-assignment-cancelled");
        assert!(result.artifacts.is_empty());
    });
}

#[test]
fn g0_g3_interface_audit_detects_seeded_omission_at_the_tolerance_boundary() {
    let spec = spec_without_interface();
    let artifact = spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut library = ImportedMeshLibrary::new();
    library.insert(&artifact, cube(), "m", named_groups());
    let gate = CancelGate::new_clock_free();

    with_cx(&gate, |cx| {
        let resolution =
            resolve_geometry_assignments(&spec, &library, AssignmentLimits::DEFAULT, cx);
        assert!(resolution.admissible(), "{:?}", resolution.violations);

        let near_miss = audit_interface_declarations(
            &spec,
            &library,
            AssignmentLimits::DEFAULT,
            InterfaceAuditLimits {
                proximity_tolerance: 0.99,
                max_triangle_pair_tests: 16,
            },
            cx,
        );
        assert!(near_miss.admissible(), "{near_miss:?}");
        assert!(near_miss.undeclared_contacts.is_empty());

        let detected = audit_interface_declarations(
            &spec,
            &library,
            AssignmentLimits::DEFAULT,
            InterfaceAuditLimits {
                proximity_tolerance: 1.0,
                max_triangle_pair_tests: 16,
            },
            cx,
        );
        assert_eq!(
            detected.violations[0].code,
            "project-interface-undeclared-contact"
        );
        assert_eq!(detected.undeclared_contacts.len(), 1);
        let contact = &detected.undeclared_contacts[0];
        assert_eq!(
            (
                contact.first_region.as_str(),
                contact.second_region.as_str()
            ),
            ("bottom", "top")
        );
        assert_eq!(contact.separation, 1.0);
        assert_eq!(contact.length_unit, "m");
        assert!(
            InterfaceDeclarationAudit::no_claim().contains("does not certify continuum contact")
        );
    });
}

#[test]
fn g0_g4_declared_pairs_are_exempt_and_resource_refusal_publishes_no_partial_list() {
    let declared_spec = named_spec();
    let artifact = declared_spec.geometry.as_ref().expect("geometry")[0].clone();
    let mut library = ImportedMeshLibrary::new();
    library.insert(&artifact, cube(), "m", named_groups());
    let gate = CancelGate::new_clock_free();

    with_cx(&gate, |cx| {
        let declared_resolution =
            resolve_geometry_assignments(&declared_spec, &library, AssignmentLimits::DEFAULT, cx);
        assert!(
            declared_resolution.admissible(),
            "{:?}",
            declared_resolution.violations
        );
        let declared = audit_interface_declarations(
            &declared_spec,
            &library,
            AssignmentLimits::DEFAULT,
            InterfaceAuditLimits {
                proximity_tolerance: 1.0,
                max_triangle_pair_tests: 0,
            },
            cx,
        );
        assert!(declared.admissible(), "{declared:?}");
        assert_eq!(declared.triangle_pair_tests, 0);

        let omitted_spec = spec_without_interface();
        let omitted_resolution =
            resolve_geometry_assignments(&omitted_spec, &library, AssignmentLimits::DEFAULT, cx);
        assert!(
            omitted_resolution.admissible(),
            "{:?}",
            omitted_resolution.violations
        );
        let bounded = audit_interface_declarations(
            &omitted_spec,
            &library,
            AssignmentLimits::DEFAULT,
            InterfaceAuditLimits {
                proximity_tolerance: 1.0,
                max_triangle_pair_tests: 1,
            },
            cx,
        );
        assert_eq!(
            bounded.violations[0].code,
            "project-interface-audit-resource-bound"
        );
        assert!(bounded.undeclared_contacts.is_empty());
        assert_eq!(bounded.triangle_pair_tests, 1);
    });
}
