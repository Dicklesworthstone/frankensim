//! G0/G3/G4 coverage for the caller-tessellated STEP-to-estimated-SDF handoff.

use fs_evidence::{NumericalCertificate, NumericalKind, vv::UnitId};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_io::{
    ParsedStep, StepImportRefusal, StepMeshDefectKind, StepTessellatorIdentity,
    import_step_tessellation, parse_step,
};
use fs_rep_mesh::{MeshSdfError, Soup, shapes};

fn parsed_fixture() -> ParsedStep {
    let source = b"ISO-10303-21;\n\
        HEADER;\n\
        FILE_DESCRIPTION(('caller tessellation fixture'),'2;1');\n\
        FILE_NAME('fixture.step','2026-07-16T00:00:00',('fs-io'),\
        ('FrankenSim'),'fs-io','FrankenSim','');\n\
        FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\n\
        ENDSEC;\n\
        DATA;\n\
        #1=FIXTURE();\n\
        ENDSEC;\n\
        END-ISO-10303-21;\n";
    parse_step(source).expect("bounded fixture parses")
}

fn tessellator() -> StepTessellatorIdentity {
    StepTessellatorIdentity {
        name: "fixture-tessellator".to_string(),
        version: "1.0".to_string(),
        configuration_fingerprint: 0x51_7e_11_a7,
    }
}

fn deviation() -> NumericalCertificate {
    NumericalCertificate::estimate(0.0, 0.02)
}

fn length_unit() -> UnitId {
    UnitId::try_new("m").expect("fixture unit is admitted")
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_mode_cx(gate, ExecMode::Deterministic, f)
}

fn with_mode_cx<R>(gate: &CancelGate, mode: ExecMode, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x57_E9,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            mode,
        );
        f(&cx)
    })
}

#[test]
fn step_import_001_closed_mesh_yields_source_bound_estimate_receipt() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let first = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect("closed outward cube is admitted");
        let second = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect("repeated deterministic import is admitted");

        let receipt = first.receipt();
        assert_eq!(
            receipt.source_fingerprint(),
            parsed.receipt().source_fingerprint()
        );
        assert_eq!(
            receipt.canonical_layout_fingerprint(),
            parsed.receipt().canonical_layout_fingerprint()
        );
        assert_eq!(
            receipt.schema_identifiers(),
            &["CONFIG_CONTROL_DESIGN".to_string()]
        );
        assert_eq!(receipt.tessellator(), &tessellator());
        assert_eq!(receipt.length_unit(), &length_unit());
        assert_eq!(receipt.target_h(), 1.0);
        assert_eq!(receipt.tessellation_deviation(), deviation());
        assert_eq!(receipt.combined_numerical().kind, NumericalKind::Estimate);
        assert_eq!(receipt.combined_numerical().lo, 0.0);
        assert!(
            receipt.combined_numerical().hi > receipt.mesh_sdf_numerical().hi + deviation().hi,
            "combined upper bound is outward-rounded"
        );
        assert_eq!(first.evidence().numerical, receipt.combined_numerical());
        assert_eq!(first.evidence().qoi, receipt.combined_numerical().hi);
        assert_eq!(first.evidence().provenance, receipt.output_provenance());
        assert!(receipt.quality().passes_basic_orientation_checks());

        let json = receipt.to_json();
        assert_eq!(json, second.receipt().to_json());
        assert!(json.contains("\"authority\":\"estimate\""));
        assert!(json.contains("\"sign_confidence\":\"uncertified\""));
        assert!(json.contains("\"source_fingerprint_fnv1a64\""));
        assert!(json.contains("\"step_import_semantics\":\"step-tessellation-to-sdf-v1\""));
        assert!(json.contains("\"execution_mode\":\"deterministic\""));
        assert!(json.contains("\"tessellation_fingerprint_domain\""));
        assert!(json.contains("caller-supplied tessellation"));

        let different_deviation = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            NumericalCertificate::estimate(0.0, 0.03),
            length_unit(),
            1.0,
            cx,
        )
        .expect("alternate declared deviation remains admissible");
        assert_ne!(
            receipt.output_provenance(),
            different_deviation.receipt().output_provenance(),
            "changed numerical claim must move provenance"
        );

        let mut deformed = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        deformed.positions[6].x = 1.125;
        let deformed = import_step_tessellation(
            &parsed,
            deformed,
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect("closed deformed cube remains admissible");
        assert_ne!(
            receipt.source_tessellation_fingerprint(),
            deformed.receipt().source_tessellation_fingerprint(),
            "changed soup bits must move the tessellation identity"
        );
        assert_ne!(
            receipt.output_provenance(),
            deformed.receipt().output_provenance(),
            "changed soup must move output provenance"
        );
    });
}

#[test]
fn step_import_002_open_mesh_refuses_with_localized_boundary_edges() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let mut soup = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        soup.triangles.pop();
        let error = import_step_tessellation(
            &parsed,
            soup,
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("open mesh must not publish an SDF");

        let StepImportRefusal::MeshIntegrity {
            source_fingerprint,
            quality,
            defects,
            localized_truncated,
            repairs,
        } = error
        else {
            panic!("unexpected refusal variant");
        };
        assert_eq!(source_fingerprint, parsed.receipt().source_fingerprint());
        assert!(quality.boundary_edges > 0);
        assert!(!localized_truncated);
        assert!(defects.iter().any(|defect| {
            defect.kind == StepMeshDefectKind::BoundaryEdge
                && defect.edge.is_some()
                && defect.incident_faces[0].is_some()
                && defect.total_incidents == 1
        }));
        assert!(
            repairs
                .iter()
                .all(|receipt| receipt.defect != "boundary-hole"),
            "zero hole-fill budget must not conceal the leak"
        );
    });
}

#[test]
fn step_import_003_non_manifold_edge_is_deterministically_localized() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let soup = Soup {
            positions: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
            ],
            triangles: vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]],
        };
        let error = import_step_tessellation(
            &parsed,
            soup,
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("three faces on one edge must refuse");

        let StepImportRefusal::MeshIntegrity {
            quality, defects, ..
        } = error
        else {
            panic!("unexpected refusal variant");
        };
        assert!(quality.nonmanifold_edges > 0);
        let edge = defects
            .iter()
            .find(|defect| defect.kind == StepMeshDefectKind::NonManifoldEdge)
            .expect("non-manifold edge localized");
        assert_eq!(edge.edge, Some([0, 1]));
        assert_eq!(edge.total_incidents, 3);
        assert_eq!(edge.incident_faces, [Some(0), Some(1), Some(2)]);
        assert_eq!(edge.kind.label(), "non-manifold-edge");
    });
}

#[test]
fn step_import_004_disconnected_closed_vertex_link_refuses() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let soup = Soup {
            positions: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.0, 0.0, 1.0),
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(0.0, -1.0, 0.0),
                Point3::new(0.0, 0.0, -1.0),
                Point3::new(10.0, 0.0, 0.0),
                Point3::new(11.0, 0.0, 0.0),
                Point3::new(10.0, 1.0, 0.0),
            ],
            triangles: vec![
                [0, 2, 1],
                [0, 1, 3],
                [0, 3, 2],
                [1, 2, 3],
                [0, 4, 5],
                [0, 6, 4],
                [0, 5, 6],
                [4, 6, 5],
                [7, 8, 9],
            ],
        };
        let error = import_step_tessellation(
            &parsed,
            soup,
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("two closed face fans sharing only a vertex must refuse");
        let StepImportRefusal::MeshIntegrity {
            quality, defects, ..
        } = error
        else {
            panic!("unexpected refusal variant");
        };
        assert!(
            quality.boundary_edges > 0,
            "the isolated face keeps the basic edge gate red"
        );
        let defect = defects
            .iter()
            .find(|defect| defect.kind == StepMeshDefectKind::VertexLinkNonManifold)
            .expect("disconnected vertex link localized");
        assert_eq!(defect.vertex, Some(0));
        assert!(defect.total_incidents > 0);
    });
}

#[test]
fn step_import_005_hostile_handoff_inputs_fail_at_admission() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let admission = |soup, adapter, deviation, target_h| {
            import_step_tessellation(
                &parsed,
                soup,
                adapter,
                deviation,
                length_unit(),
                target_h,
                cx,
            )
            .expect_err("invalid handoff must refuse")
        };

        let mut non_finite = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        non_finite.positions[0].x = f64::NAN;
        assert!(matches!(
            admission(non_finite, tessellator(), deviation(), 1.0),
            StepImportRefusal::Admission { .. }
        ));

        let mut missing_vertex = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        missing_vertex.triangles[0][2] = u32::MAX;
        assert!(matches!(
            admission(missing_vertex, tessellator(), deviation(), 1.0),
            StepImportRefusal::Admission { .. }
        ));

        let mut blank_adapter = tessellator();
        blank_adapter.name.clear();
        assert!(matches!(
            admission(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                blank_adapter,
                deviation(),
                1.0,
            ),
            StepImportRefusal::Admission { .. }
        ));
        assert!(matches!(
            admission(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                tessellator(),
                NumericalCertificate {
                    kind: NumericalKind::Exact,
                    lo: 0.0,
                    hi: 0.1,
                },
                1.0,
            ),
            StepImportRefusal::Admission { .. }
        ));

        assert!(matches!(
            admission(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                tessellator(),
                NumericalCertificate::no_claim(),
                1.0,
            ),
            StepImportRefusal::Admission { .. }
        ));
        assert!(matches!(
            admission(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                tessellator(),
                deviation(),
                0.0,
            ),
            StepImportRefusal::Admission { .. }
        ));
    });
}

#[test]
fn step_import_006_safe_repairs_are_retained_without_hole_filling() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let mut soup = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        soup.positions.push(Point3::new(9.0, 9.0, 9.0));
        soup.triangles.push(soup.triangles[0]);
        soup.triangles.push([0, 0, 1]);
        let outcome = import_step_tessellation(
            &parsed,
            soup,
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect("duplicate and degenerate faces are removable");

        assert!(
            outcome
                .receipt()
                .repairs()
                .iter()
                .any(|receipt| receipt.defect == "duplicate-face")
        );
        assert!(
            outcome
                .receipt()
                .repairs()
                .iter()
                .any(|receipt| receipt.defect == "degenerate-face")
        );
        assert!(
            outcome
                .receipt()
                .repairs()
                .iter()
                .any(|receipt| receipt.defect == "unreferenced-vertex")
        );
        assert!(
            outcome
                .receipt()
                .repairs()
                .iter()
                .all(|receipt| receipt.defect != "boundary-hole")
        );

        let exhausted = import_step_tessellation(
            &parsed,
            Soup {
                positions: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                triangles: vec![[0, 0, 1]],
            },
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("repair-exhausted soup must retain a structured refusal");
        let StepImportRefusal::MeshIntegrity { repairs, .. } = exhausted else {
            panic!("unexpected refusal variant");
        };
        assert!(
            repairs
                .iter()
                .any(|receipt| receipt.defect == "degenerate-face"),
            "post-repair refusal retains the repair audit trail"
        );
    });
}

#[test]
fn step_import_007_pre_requested_cancellation_refuses_at_entry() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let error = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("pre-requested cancellation must refuse publication");
        assert!(matches!(
            error,
            StepImportRefusal::SdfBuild {
                error: MeshSdfError::Cancelled,
                ..
            }
        ));
    });
}

#[test]
fn step_import_008_outward_rounding_overflow_refuses_evidence() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed_fixture();
        let error = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            NumericalCertificate::estimate(0.0, f64::MAX),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("infinite outward-rounded upper bound must refuse");
        assert!(matches!(error, StepImportRefusal::Evidence { .. }));
    });
}

#[test]
fn step_import_009_fast_mode_refuses_d0_publication() {
    let gate = CancelGate::new_clock_free();
    with_mode_cx(&gate, ExecMode::Fast, |cx| {
        let parsed = parsed_fixture();
        let error = import_step_tessellation(
            &parsed,
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            tessellator(),
            deviation(),
            length_unit(),
            1.0,
            cx,
        )
        .expect_err("D0 publication must reject fast mode");
        assert!(matches!(
            error,
            StepImportRefusal::Admission { what, .. } if what.contains("deterministic")
        ));
    });
}
