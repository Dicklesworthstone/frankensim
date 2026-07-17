//! G0/G3/G4/G5 coverage for sealed strict FACETED_BREP re-emission.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::{
    STEP_FACETED_EXPORT_VERSION, StepFacetedExportLimits, StepFacetedExportMetadata,
    StepFacetedExportRefusal, StepFacetedLimits, StepFacetedProfile,
    decode_faceted_brep_with_limits, export_decoded_faceted_brep,
    export_decoded_faceted_brep_with_limits, parse_step,
};

fn tetra_source() -> String {
    "ISO-10303-21;\n\
     HEADER;\n\
     FILE_DESCRIPTION(('strict faceted export fixture'),'2;1');\n\
     FILE_NAME('source.step','2026-07-17T00:00:00',('fs-io'),('FrankenSim'),'fs-io','FrankenSim','');\n\
     FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\n\
     ENDSEC;\n\
     DATA;\n\
     #60=FACETED_BREP('',#50);\n\
     #43=FACE('',(#42));\n\
     #21=POLY_LOOP('',(#1,#2,#4));\n\
     #2=CARTESIAN_POINT('',(1.2345678901234567,-0.0,0.0));\n\
     #42=FACE_OUTER_BOUND('',#41,.T.);\n\
     #12=FACE_OUTER_BOUND('',#11,.T.);\n\
     #31=POLY_LOOP('',(#1,#4,#3));\n\
     #4=CARTESIAN_POINT('',(0.0,0.0,1.2345678901234567E+100));\n\
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

fn metadata() -> StepFacetedExportMetadata {
    StepFacetedExportMetadata::new(
        "supplier-part.step",
        "2026-07-17T12:34:56",
        "O'Connor",
        "FrankenSim",
        "controlled-export",
    )
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0xfa_ce_7e_e0,
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

fn assert_same_soup_bits(left: &fs_rep_mesh::Soup, right: &fs_rep_mesh::Soup) {
    assert_eq!(left.triangles, right.triangles);
    assert_eq!(left.positions.len(), right.positions.len());
    for (left, right) in left.positions.iter().zip(&right.positions) {
        assert_eq!(left.x.to_bits(), right.x.to_bits());
        assert_eq!(left.y.to_bits(), right.y.to_bits());
        assert_eq!(left.z.to_bits(), right.z.to_bits());
    }
}

#[test]
fn step_faceted_export_001_reemits_exact_geometry_with_nested_receipts() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parse_step(tetra_source().as_bytes()).expect("source syntax");
        let source = decode_faceted_brep_with_limits(&parsed, 60, StepFacetedLimits::default(), cx)
            .expect("sealed source resource");
        let exported = export_decoded_faceted_brep(&source, metadata(), cx)
            .expect("strict faceted resource export");

        assert_eq!(
            exported.receipt().writer_version(),
            STEP_FACETED_EXPORT_VERSION
        );
        assert_eq!(exported.receipt().metadata(), &metadata());
        assert_eq!(
            exported.receipt().source_decoder_receipt(),
            source.receipt()
        );
        assert_eq!(
            exported
                .receipt()
                .output_structure_receipt()
                .schema_identifiers(),
            &["CONFIG_CONTROL_DESIGN".to_string()]
        );
        assert_eq!(
            exported.receipt().output_decoder_receipt().profile(),
            StepFacetedProfile::ConfigControlDesign
        );
        assert_eq!(
            exported.receipt().output_decoder_receipt().vertex_count(),
            source.soup().positions.len()
        );
        assert_eq!(
            exported.receipt().output_decoder_receipt().triangle_count(),
            source.soup().triangles.len()
        );
        assert_eq!(
            exported
                .receipt()
                .output_decoder_receipt()
                .bare_face_count(),
            source.soup().triangles.len()
        );

        let replay = parse_step(exported.bytes()).expect("exported syntax reparses");
        let decoded = decode_faceted_brep_with_limits(
            &replay,
            exported.receipt().output_decoder_receipt().root_id(),
            StepFacetedLimits::default(),
            cx,
        )
        .expect("exported semantic resource replays");
        assert_same_soup_bits(source.soup(), decoded.soup());
        assert_eq!(decoded.soup().positions[1].y.to_bits(), (-0.0f64).to_bits());
    });
}

#[test]
fn step_faceted_export_002_is_byte_deterministic_and_preserves_profile() {
    let automotive = tetra_source().replace("CONFIG_CONTROL_DESIGN", "AUTOMOTIVE_DESIGN");
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parse_step(automotive.as_bytes()).expect("automotive source syntax");
        let source = decode_faceted_brep_with_limits(&parsed, 60, StepFacetedLimits::default(), cx)
            .expect("sealed automotive resource");
        let first = export_decoded_faceted_brep(&source, metadata(), cx).expect("first export");
        let replay = export_decoded_faceted_brep(&source, metadata(), cx).expect("replayed export");

        assert_eq!(first, replay);
        assert_eq!(first.bytes(), replay.bytes());
        assert_eq!(
            first.receipt().output_decoder_receipt().profile(),
            StepFacetedProfile::AutomotiveDesign
        );
        assert!(
            core::str::from_utf8(first.bytes())
                .expect("Part-21 output is ASCII")
                .contains("O''Connor"),
            "STEP string apostrophes remain canonically doubled"
        );
    });
}

#[test]
fn step_faceted_export_003_bounds_metadata_and_cancellation_refuse() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let parsed = parse_step(tetra_source().as_bytes()).expect("source syntax");
        let source = decode_faceted_brep_with_limits(&parsed, 60, StepFacetedLimits::default(), cx)
            .expect("sealed source resource");

        let mut instance_limits = StepFacetedExportLimits::default();
        instance_limits.syntax.max_instances = 17;
        assert!(matches!(
            export_decoded_faceted_brep_with_limits(&source, metadata(), instance_limits, cx,),
            Err(StepFacetedExportRefusal::Resource {
                stage: "instance-plan",
                ..
            })
        ));

        let mut vertex_limits = StepFacetedExportLimits::default();
        vertex_limits.faceted.max_vertices = 3;
        assert!(matches!(
            export_decoded_faceted_brep_with_limits(&source, metadata(), vertex_limits, cx,),
            Err(StepFacetedExportRefusal::Resource {
                stage: "source-vertices",
                ..
            })
        ));

        let plane_source = tetra_source()
            .replace(
                "#13=FACE('',(#12));",
                "#13=FACE_SURFACE('',(#12),#113,.T.);",
            )
            .replace(
                "#70=UNRELATED_ENTITY($);",
                "#113=PLANE('',#114);\n\
                 #114=AXIS2_PLACEMENT_3D('',#1,#115,#116);\n\
                 #115=DIRECTION('',(0.0,0.0,-1.0));\n\
                 #116=DIRECTION('',(1.0,0.0,0.0));\n\
                 #70=UNRELATED_ENTITY($);",
            );
        let plane_parsed = parse_step(plane_source.as_bytes()).expect("plane source syntax");
        let plane_source =
            decode_faceted_brep_with_limits(&plane_parsed, 60, StepFacetedLimits::default(), cx)
                .expect("sealed plane-backed source resource");
        assert!(matches!(
            export_decoded_faceted_brep(&plane_source, metadata(), cx),
            Err(StepFacetedExportRefusal::Admission { .. })
        ));

        let invalid_metadata = StepFacetedExportMetadata::new(
            "supplier.step",
            "2026-07-17",
            "non-ascii-autor-\u{00e9}",
            "FrankenSim",
            "",
        );
        assert!(matches!(
            export_decoded_faceted_brep(&source, invalid_metadata, cx),
            Err(StepFacetedExportRefusal::Syntax { stage: "write", .. })
        ));
    });

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        let active = CancelGate::new_clock_free();
        with_cx(&active, |decode_cx| {
            let parsed = parse_step(tetra_source().as_bytes()).expect("source syntax");
            let source = decode_faceted_brep_with_limits(
                &parsed,
                60,
                StepFacetedLimits::default(),
                decode_cx,
            )
            .expect("sealed source resource");
            assert!(matches!(
                export_decoded_faceted_brep(&source, metadata(), cx),
                Err(StepFacetedExportRefusal::Cancelled { stage: "entry" })
            ));
        });
    });
}
