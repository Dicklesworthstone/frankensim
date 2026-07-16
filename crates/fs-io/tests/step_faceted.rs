//! G0/G3/G4 coverage for strict native triangular FACETED_BREP decoding.

use fs_evidence::{NumericalKind, vv::UnitId};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::{
    STEP_FACETED_DECODER_VERSION, StepFacetedImportRefusal, StepFacetedLimits, StepFacetedProfile,
    StepFacetedRefusal, StepImportRefusal, decode_faceted_brep_with_limits, import_faceted_brep,
    parse_step,
};

fn tetra_source() -> String {
    "ISO-10303-21;\n\
     HEADER;\n\
     FILE_DESCRIPTION(('strict faceted fixture'),'2;1');\n\
     FILE_NAME('faceted.step','2026-07-16T00:00:00',('fs-io'),('FrankenSim'),'fs-io','FrankenSim','');\n\
     FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\n\
     ENDSEC;\n\
     DATA;\n\
     #60=FACETED_BREP('',#50);\n\
     #43=FACE('',(#42));\n\
     #21=POLY_LOOP('',(#1,#2,#4));\n\
     #2=CARTESIAN_POINT('',(1.0,0.0,0.0));\n\
     #42=FACE_OUTER_BOUND('',#41,.T.);\n\
     #12=FACE_OUTER_BOUND('',#11,.T.);\n\
     #31=POLY_LOOP('',(#1,#4,#3));\n\
     #4=CARTESIAN_POINT('',(0.0,0.0,1.0));\n\
     #50=CLOSED_SHELL('',(#43,#13,#33,#23));\n\
     #11=POLY_LOOP('',(#1,#3,#2));\n\
     #23=FACE('',(#22));\n\
     #1=CARTESIAN_POINT('',(0.0,0.0,0.0));\n\
     #32=FACE_OUTER_BOUND('',#31,.T.);\n\
     #41=POLY_LOOP('',(#2,#3,#4));\n\
     #13=FACE('',(#12));\n\
     #3=CARTESIAN_POINT('',(0.0,1.0,0.0));\n\
     #22=FACE_OUTER_BOUND('',#21,.T.);\n\
     #33=FACE('',(#32));\n\
     #70=UNRELATED_ENTITY($);\n\
     ENDSEC;\n\
     END-ISO-10303-21;\n"
        .to_string()
}

fn plane_tetra_source() -> String {
    tetra_source()
        .replace(
            "#13=FACE('',(#12));",
            "#13=FACE_SURFACE('',(#12),#113,.T.);",
        )
        .replace(
            "#23=FACE('',(#22));",
            "#23=FACE_SURFACE('',(#22),#123,.T.);",
        )
        .replace(
            "#33=FACE('',(#32));",
            "#33=FACE_SURFACE('',(#32),#133,.T.);",
        )
        .replace(
            "#43=FACE('',(#42));",
            "#43=FACE_SURFACE('',(#42),#143,.T.);",
        )
        .replace(
            "#70=UNRELATED_ENTITY($);",
            "#113=PLANE('',#114);\n\
             #114=AXIS2_PLACEMENT_3D('',#1,#115,#116);\n\
             #115=DIRECTION('',(0.0,0.0,-1.0));\n\
             #116=DIRECTION('',(1.0,0.0,0.0));\n\
             #123=PLANE('',#124);\n\
             #124=AXIS2_PLACEMENT_3D('',#1,#125,#126);\n\
             #125=DIRECTION('',(0.0,-1.0,0.0));\n\
             #126=DIRECTION('',(1.0,0.0,0.0));\n\
             #133=PLANE('',#134);\n\
             #134=AXIS2_PLACEMENT_3D('',#1,#135,#136);\n\
             #135=DIRECTION('',(-1.0,0.0,0.0));\n\
             #136=DIRECTION('',(0.0,1.0,0.0));\n\
             #143=PLANE('',#144);\n\
             #144=AXIS2_PLACEMENT_3D('',#2,#145,#146);\n\
             #145=DIRECTION('',(1.0,1.0,1.0));\n\
             #146=DIRECTION('',(0.0,1.0,0.0));\n\
             #70=UNRELATED_ENTITY($);",
        )
}

fn parsed(source: &str) -> fs_io::ParsedStep {
    parse_step(source.as_bytes()).expect("bounded FACETED_BREP fixture parses")
}

fn length_unit() -> UnitId {
    UnitId::try_new("m").expect("fixture unit is admitted")
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xfa_ce_7e_d0,
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

#[test]
fn step_faceted_001_unsorted_tetra_decodes_to_canonical_triangle_soup() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parsed(&tetra_source());
        let first = decode_faceted_brep_with_limits(&parsed, 60, StepFacetedLimits::default(), cx)
            .expect("strict tetra closure decodes");
        let second = decode_faceted_brep_with_limits(&parsed, 60, StepFacetedLimits::default(), cx)
            .expect("repeated strict decoding succeeds");

        assert_eq!(first.soup().positions.len(), 4);
        assert_eq!(first.soup().triangles.len(), 4);
        assert_eq!(
            first.soup().triangles,
            vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
        );
        assert_eq!(first.soup().positions, second.soup().positions);
        assert_eq!(first.soup().triangles, second.soup().triangles);

        let receipt = first.receipt();
        assert_eq!(receipt.profile(), StepFacetedProfile::ConfigControlDesign);
        assert_eq!(receipt.root_id(), 60);
        assert_eq!(receipt.shell_id(), 50);
        assert_eq!(receipt.vertex_count(), 4);
        assert_eq!(receipt.triangle_count(), 4);
        assert_eq!(receipt.bare_face_count(), 4);
        assert_eq!(receipt.plane_face_count(), 0);
        assert_eq!(receipt.reversed_bounds(), 0);
        assert_eq!(
            receipt.semantic_fingerprint(),
            second.receipt().semantic_fingerprint()
        );
        assert_eq!(
            receipt.coordinate_conversion().kind,
            NumericalKind::Estimate
        );
        assert_eq!(receipt.coordinate_conversion().lo, 0.0);
        assert!(receipt.coordinate_conversion().hi.is_finite());
        assert_eq!(
            receipt.materialization_deviation(),
            receipt.coordinate_conversion()
        );

        let json = receipt.to_json();
        assert!(json.contains(STEP_FACETED_DECODER_VERSION));
        assert!(json.contains("schema-declaration-gated-resource-decoding"));
        assert!(json.contains("\"triangles\":4"));
        assert!(json.contains("no full EXPRESS or AP conformance"));
    });
}

#[test]
fn step_faceted_002_false_bound_orientation_reverses_only_that_loop() {
    let source = tetra_source().replace(
        "#12=FACE_OUTER_BOUND('',#11,.T.);",
        "#12=FACE_OUTER_BOUND('',#11,.F.);",
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let decoded =
            decode_faceted_brep_with_limits(&parsed(&source), 60, StepFacetedLimits::default(), cx)
                .expect("false orientation remains in the admitted subset");
        assert_eq!(decoded.soup().triangles[0], [0, 1, 2]);
        assert_eq!(decoded.soup().triangles[1], [0, 1, 3]);
        assert_eq!(decoded.receipt().reversed_bounds(), 1);
    });
}

#[test]
fn step_faceted_003_shell_set_permutation_preserves_semantics() {
    let permuted = tetra_source().replace(
        "#50=CLOSED_SHELL('',(#43,#13,#33,#23));",
        "#50=CLOSED_SHELL('',(#23,#33,#13,#43));",
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let original = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("original shell decodes");
        let permuted = decode_faceted_brep_with_limits(
            &parsed(&permuted),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("permuted SET decodes");

        assert_eq!(original.soup().positions, permuted.soup().positions);
        assert_eq!(original.soup().triangles, permuted.soup().triangles);
        assert_eq!(
            original.receipt().semantic_fingerprint(),
            permuted.receipt().semantic_fingerprint()
        );
        assert_ne!(
            original.receipt().source_fingerprint(),
            permuted.receipt().source_fingerprint(),
            "source spelling remains separately identity-bearing"
        );
    });
}

#[test]
fn step_faceted_004_schema_gate_is_exact_and_non_authoritative() {
    let automotive = tetra_source().replace("CONFIG_CONTROL_DESIGN", "AUTOMOTIVE_DESIGN");
    let unsupported = tetra_source().replace(
        "CONFIG_CONTROL_DESIGN",
        "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF",
    );
    let ambiguous = tetra_source().replace(
        "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));",
        "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN','AUTOMOTIVE_DESIGN'));",
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let decoded = decode_faceted_brep_with_limits(
            &parsed(&automotive),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("exact AUTOMOTIVE_DESIGN declaration is admitted");
        assert_eq!(
            decoded.receipt().profile(),
            StepFacetedProfile::AutomotiveDesign
        );

        for source in [&unsupported, &ambiguous] {
            let error = decode_faceted_brep_with_limits(
                &parsed(source),
                60,
                StepFacetedLimits::default(),
                cx,
            )
            .expect_err("unsupported or ambiguous schema declarations refuse");
            assert!(matches!(error, StepFacetedRefusal::Schema { .. }));
        }
    });
}

#[test]
fn step_faceted_005_polygon_coordinate_and_resource_drift_refuse() {
    let quad = tetra_source().replace("(#1,#3,#2)", "(#1,#3,#2,#4)");
    let duplicate = tetra_source().replace("(#1,#3,#2)", "(#1,#3,#1)");
    let non_finite = tetra_source().replace("(1.0,0.0,0.0)", "(1.0E+999,0.0,0.0)");
    let malformed_face_surface =
        tetra_source().replace("#13=FACE('',(#12));", "#13=FACE_SURFACE('',(#12),#1,.T.);");
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for source in [&quad, &duplicate, &non_finite, &malformed_face_surface] {
            let error = decode_faceted_brep_with_limits(
                &parsed(source),
                60,
                StepFacetedLimits::default(),
                cx,
            )
            .expect_err("semantic drift must refuse");
            assert!(matches!(error, StepFacetedRefusal::Entity { .. }));
        }

        let error = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits {
                max_vertices: 3,
                ..StepFacetedLimits::default()
            },
            cx,
        )
        .expect_err("injected vertex limit must refuse before publication");
        assert!(matches!(
            error,
            StepFacetedRefusal::Resource {
                stage: "point-plan",
                ..
            }
        ));

        let error = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits {
                max_triangles: 3,
                ..StepFacetedLimits::default()
            },
            cx,
        )
        .expect_err("injected triangle limit must refuse the shell face SET");
        assert!(matches!(
            error,
            StepFacetedRefusal::Resource {
                stage: "shell-faces",
                ..
            }
        ));

        let error = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits {
                max_auxiliary_bytes: 768,
                ..StepFacetedLimits::default()
            },
            cx,
        )
        .expect_err("injected auxiliary-memory limit must refuse");
        assert!(matches!(
            error,
            StepFacetedRefusal::Resource {
                stage: "semantic-plan",
                ..
            }
        ));
    });
}

#[test]
fn step_faceted_006_pre_requested_cancellation_refuses_before_materialization() {
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        let error = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect_err("pre-requested cancellation must refuse");
        assert!(matches!(
            error,
            StepFacetedRefusal::Cancelled { stage: "entry", .. }
        ));
    });
}

#[test]
fn step_faceted_007_native_bridge_reuses_topology_quarantine() {
    let two_roots = tetra_source().replace(
        "#60=FACETED_BREP('',#50);\n",
        "#60=FACETED_BREP('',#50);\n#61=FACETED_BREP('',#50);\n",
    );
    let open = tetra_source().replace(
        "#50=CLOSED_SHELL('',(#43,#13,#33,#23));",
        "#50=CLOSED_SHELL('',(#13,#23,#33));",
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed_closed = parsed(&two_roots);
        let outcome = import_faceted_brep(&parsed_closed, 60, length_unit(), 1.0, cx)
            .expect("closed outward tetra reaches estimated-SDF publication");
        let alternate_root = import_faceted_brep(&parsed_closed, 61, length_unit(), 1.0, cx)
            .expect("second explicit root reaches the same materialized soup");
        assert_eq!(outcome.decoder_receipt().triangle_count(), 4);
        assert_eq!(
            outcome.import().receipt().tessellator().name,
            "fs-io-native-faceted-brep"
        );
        assert_eq!(
            outcome.import().receipt().tessellation_deviation(),
            outcome.decoder_receipt().materialization_deviation()
        );
        assert_eq!(
            outcome.import().receipt().source_tessellation_fingerprint(),
            alternate_root
                .import()
                .receipt()
                .source_tessellation_fingerprint(),
            "both roots deliberately materialize the same soup"
        );
        assert_ne!(
            outcome.decoder_receipt().semantic_fingerprint(),
            alternate_root.decoder_receipt().semantic_fingerprint(),
            "the caller-selected root remains semantic input"
        );
        assert_ne!(
            outcome
                .import()
                .receipt()
                .tessellator()
                .configuration_fingerprint,
            alternate_root
                .import()
                .receipt()
                .tessellator()
                .configuration_fingerprint,
            "native materializer identity binds the selected semantic closure"
        );
        assert_ne!(
            outcome.import().receipt().output_provenance(),
            alternate_root.import().receipt().output_provenance(),
            "downstream evidence provenance must bind the selected native root"
        );

        let error = import_faceted_brep(&parsed(&open), 60, length_unit(), 1.0, cx)
            .expect_err("open native closure must not bypass topology quarantine");
        assert!(matches!(
            error,
            StepFacetedImportRefusal::Import {
                decoder_receipt,
                error: StepImportRefusal::MeshIntegrity { quality, .. },
            } if decoder_receipt.root_id() == 60 && quality.boundary_edges > 0
        ));
    });
}

#[test]
fn step_faceted_008_plane_backed_faces_preserve_soup_and_bind_plane_semantics() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let bare = decode_faceted_brep_with_limits(
            &parsed(&tetra_source()),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("bare tetra decodes");
        let plane = decode_faceted_brep_with_limits(
            &parsed(&plane_tetra_source()),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("plane-backed tetra passes coplanarity and orientation checks");

        assert_eq!(plane.soup().positions, bare.soup().positions);
        assert_eq!(plane.soup().triangles, bare.soup().triangles);
        assert_eq!(plane.receipt().bare_face_count(), 0);
        assert_eq!(plane.receipt().plane_face_count(), 4);
        assert_eq!(
            plane.receipt().plane_consistency().kind,
            NumericalKind::Estimate
        );
        assert_eq!(plane.receipt().plane_consistency().lo, 0.0);
        assert!(plane.receipt().plane_consistency().hi.is_finite());
        assert!(
            plane.receipt().materialization_deviation().hi
                >= plane.receipt().coordinate_conversion().hi
        );
        assert_ne!(
            plane.receipt().semantic_fingerprint(),
            bare.receipt().semantic_fingerprint(),
            "the plane closure is provenance-bearing even when the soup agrees"
        );
        assert!(plane.receipt().to_json().contains("\"plane_faces\":4"));

        let reversed_surface_normal = plane_tetra_source()
            .replace(
                "#13=FACE_SURFACE('',(#12),#113,.T.);",
                "#13=FACE_SURFACE('',(#12),#113,.F.);",
            )
            .replace(
                "#115=DIRECTION('',(0.0,0.0,-1.0));",
                "#115=DIRECTION('',(0.0,0.0,1.0));",
            );
        let reversed = decode_faceted_brep_with_limits(
            &parsed(&reversed_surface_normal),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("same_sense false paired with the reversed plane normal remains valid");
        assert_eq!(reversed.soup().triangles, plane.soup().triangles);
        assert_ne!(
            reversed.receipt().semantic_fingerprint(),
            plane.receipt().semantic_fingerprint()
        );

        let default_axis = plane_tetra_source()
            .replace(
                "#13=FACE_SURFACE('',(#12),#113,.T.);",
                "#13=FACE_SURFACE('',(#12),#113,.F.);",
            )
            .replace(
                "#114=AXIS2_PLACEMENT_3D('',#1,#115,#116);",
                "#114=AXIS2_PLACEMENT_3D('',#1,$,$);",
            );
        decode_faceted_brep_with_limits(
            &parsed(&default_axis),
            60,
            StepFacetedLimits::default(),
            cx,
        )
        .expect("omitted placement axis uses the EXPRESS default positive Z normal");

        let outcome =
            import_faceted_brep(&parsed(&plane_tetra_source()), 60, length_unit(), 1.0, cx)
                .expect("validated plane-backed soup reaches the existing quarantine");
        assert_eq!(
            outcome.import().receipt().tessellation_deviation(),
            outcome.decoder_receipt().materialization_deviation()
        );
    });
}

#[test]
fn step_faceted_009_inconsistent_plane_placement_and_direction_refuse() {
    let off_plane = plane_tetra_source().replace(
        "#114=AXIS2_PLACEMENT_3D('',#1,#115,#116);",
        "#114=AXIS2_PLACEMENT_3D('',#4,#115,#116);",
    );
    let wrong_normal = plane_tetra_source().replace(
        "#115=DIRECTION('',(0.0,0.0,-1.0));",
        "#115=DIRECTION('',(0.0,1.0,0.0));",
    );
    let parallel_axes = plane_tetra_source().replace(
        "#116=DIRECTION('',(1.0,0.0,0.0));",
        "#116=DIRECTION('',(0.0,0.0,2.0));",
    );
    let short_direction = plane_tetra_source().replace(
        "#115=DIRECTION('',(0.0,0.0,-1.0));",
        "#115=DIRECTION('',(0.0,-1.0));",
    );
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for source in [&off_plane, &wrong_normal, &parallel_axes, &short_direction] {
            let error = decode_faceted_brep_with_limits(
                &parsed(source),
                60,
                StepFacetedLimits::default(),
                cx,
            )
            .expect_err("inconsistent plane geometry must refuse before soup publication");
            assert!(matches!(error, StepFacetedRefusal::Entity { .. }));
        }
    });
}
