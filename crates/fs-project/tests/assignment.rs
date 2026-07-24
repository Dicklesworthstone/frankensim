//! G0/G3/G4 coverage for the L6 project-to-fs-io assignment adapter.

use fs_blake3::hash_domain;
use fs_conduction::{ConductionMesh, ThermalInterfaces};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{AssignmentLimits, HalfSpaceSide, MeshSelector, NamedFaceGroup};
use fs_project::assignment::{
    ConductionInterfaceLimits, ConductionInterfaceResolution, InterfaceAuditLimits,
    InterfaceDeclarationAudit, audit_interface_declarations, resolve_conduction_interface_pairs,
};
use fs_project::{
    EntityDecl, GEOMETRY_ASSIGNMENT_REPORT_DOMAIN, GeometryArtifact, GeometryAssignment,
    GeometryResolution, ImportedMeshLibrary, ProjectSpec, geometry_source_identity,
    resolve_geometry_assignments,
};
use fs_rep_mesh::{Soup, TetComplex};

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

fn two_trace_conduction_fixture() -> (
    ProjectSpec,
    ImportedMeshLibrary,
    ConductionMesh,
    usize,
    usize,
) {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
    ];
    let complex = TetComplex::from_tets(positions.len(), vec![[0, 1, 2, 3], [4, 6, 5, 7]]);
    let mesh = ConductionMesh::new(complex, positions).expect("two valid duplicated-trace tets");
    let pair = ThermalInterfaces::coincident_face_pairs(&mesh)
        .expect("fixture has valid coincident traces")
        .pop()
        .expect("fixture has one coincident pair");
    let from_slot = [pair.side_a, pair.side_b]
        .into_iter()
        .find(|&slot| mesh.boundary()[slot].outward_normal[0] > 0.0)
        .expect("left trace points in +x");
    let to_slot = if from_slot == pair.side_a {
        pair.side_b
    } else {
        pair.side_a
    };

    let mut soup_positions = Vec::new();
    let mut triangles = Vec::new();
    for slot in [from_slot, to_slot] {
        let face = &mesh.boundary()[slot];
        let start = u32::try_from(soup_positions.len()).expect("fixture indices fit u32");
        let mut coordinates = face
            .vertices
            .map(|vertex| mesh.positions()[vertex as usize]);
        let first = [
            coordinates[1][0] - coordinates[0][0],
            coordinates[1][1] - coordinates[0][1],
            coordinates[1][2] - coordinates[0][2],
        ];
        let second = [
            coordinates[2][0] - coordinates[0][0],
            coordinates[2][1] - coordinates[0][1],
            coordinates[2][2] - coordinates[0][2],
        ];
        let normal = [
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        ];
        let alignment = normal[0] * face.outward_normal[0]
            + normal[1] * face.outward_normal[1]
            + normal[2] * face.outward_normal[2];
        if alignment < 0.0 {
            coordinates.swap(1, 2);
        }
        soup_positions.extend(
            coordinates
                .into_iter()
                .map(|point| Point3::new(point[0], point[1], point[2])),
        );
        triangles.push([start, start + 1, start + 2]);
    }
    let soup = Soup {
        positions: soup_positions,
        triangles,
    };
    let artifact = GeometryArtifact {
        role: "two-trace".to_string(),
        format: "fixture".to_string(),
        source_hash: 0x17_03,
        parser_version: "fixture/v1".to_string(),
    };
    let spec = ProjectSpec {
        geometry: Some(vec![artifact.clone()]),
        assignments: Some(vec![
            GeometryAssignment {
                artifact: artifact.role.clone(),
                target: "solid-a".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "SIDE_A".to_string(),
                },
                allow_overlap: true,
            },
            GeometryAssignment {
                artifact: artifact.role.clone(),
                target: "solid-b".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "SIDE_B".to_string(),
                },
                allow_overlap: true,
            },
            GeometryAssignment {
                artifact: artifact.role.clone(),
                target: "bondline".to_string(),
                length_unit: "m".to_string(),
                selector: MeshSelector::NamedGroup {
                    name: "BONDLINE".to_string(),
                },
                allow_overlap: true,
            },
        ]),
        assembly: Some(vec![
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
                name: "solid-a".to_string(),
                display: "Solid A".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "stack".to_string(),
                name: "solid-b".to_string(),
                display: "Solid B".to_string(),
                expect_id: None,
            },
            EntityDecl::Interface {
                parent: "rig".to_string(),
                name: "bondline".to_string(),
                display: "Bondline".to_string(),
                from: "solid-a".to_string(),
                to: "solid-b".to_string(),
                expect_id: None,
            },
        ]),
        ..ProjectSpec::default()
    };
    let mut library = ImportedMeshLibrary::new();
    library.insert(
        &artifact,
        soup,
        "m",
        vec![
            NamedFaceGroup {
                name: "SIDE_A".to_string(),
                faces: vec![0],
            },
            NamedFaceGroup {
                name: "SIDE_B".to_string(),
                faces: vec![1],
            },
            NamedFaceGroup {
                name: "BONDLINE".to_string(),
                faces: vec![0, 1],
            },
        ],
    );
    (spec, library, mesh, from_slot, to_slot)
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

#[test]
fn g0_conduction_interface_lowering_binds_exact_oriented_boundary_slots() {
    let (spec, library, mesh, from_slot, to_slot) = two_trace_conduction_fixture();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let resolution = resolve_conduction_interface_pairs(
            &spec,
            &library,
            AssignmentLimits::DEFAULT,
            ConductionInterfaceLimits::DEFAULT,
            &mesh,
            cx,
        );
        assert!(resolution.admissible(), "{:?}", resolution.violations);
        assert_eq!(resolution.source_faces_indexed, 4);
        assert_eq!(resolution.pairs.len(), 1);
        let pair = &resolution.pairs[0];
        assert_eq!(pair.interface, "bondline");
        assert_eq!(pair.from_region, "solid-a");
        assert_eq!(pair.to_region, "solid-b");
        assert_eq!(pair.from_boundary_slot, from_slot);
        assert_eq!(pair.to_boundary_slot, to_slot);
        assert_eq!(pair.face_pair().side_a, from_slot);
        assert_eq!(pair.face_pair().side_b, to_slot);
        assert_eq!(pair.from_source.face, 0);
        assert_eq!(pair.to_source.face, 1);
        assert_eq!(pair.interface_sources.len(), 2);
        assert!(
            resolution
                .render_table()
                .contains("bondline | solid-a slot")
        );
        assert!(
            ConductionInterfaceResolution::no_claim()
                .contains("does not authenticate the importer")
        );
    });
}

#[test]
fn g0_g4_conduction_interface_lowering_refuses_orientation_and_budget_atomically() {
    let (mut wrong_orientation, library, mesh, _, _) = two_trace_conduction_fixture();
    let assignments = wrong_orientation
        .assignments
        .as_mut()
        .expect("fixture assignments");
    assignments[1].selector = MeshSelector::NamedGroup {
        name: "SIDE_A".to_string(),
    };
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let refused = resolve_conduction_interface_pairs(
            &wrong_orientation,
            &library,
            AssignmentLimits::DEFAULT,
            ConductionInterfaceLimits::DEFAULT,
            &mesh,
            cx,
        );
        assert_eq!(
            refused.violations[0].code,
            "project-conduction-interface-orientation"
        );
        assert!(refused.pairs.is_empty());

        let (spec, library, mesh, _, _) = two_trace_conduction_fixture();
        let bounded = resolve_conduction_interface_pairs(
            &spec,
            &library,
            AssignmentLimits::DEFAULT,
            ConductionInterfaceLimits {
                max_source_faces: 2,
            },
            &mesh,
            cx,
        );
        assert_eq!(
            bounded.violations[0].code,
            "project-conduction-interface-resource-bound"
        );
        assert_eq!(bounded.source_faces_indexed, 2);
        assert!(bounded.pairs.is_empty());

        gate.request();
        let cancelled = resolve_conduction_interface_pairs(
            &spec,
            &library,
            AssignmentLimits::DEFAULT,
            ConductionInterfaceLimits::DEFAULT,
            &mesh,
            cx,
        );
        assert_eq!(
            cancelled.violations[0].code,
            "project-conduction-interface-cancelled"
        );
        assert_eq!(cancelled.source_faces_indexed, 0);
        assert!(cancelled.pairs.is_empty());
    });
}
