//! G0/G3/G4 evidence for the CLI-layer quarantine -> assignment -> ledger path.

use fs_cli::{GeometryImportLimits, RawGeometryLibrary, import_project_geometry};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_io::{NamedFaceGroup, STEP_FACETED_DECODER_VERSION, parse_step, quarantine::import_mesh};
use fs_ledger::{EdgeRole, ExtensionTable, Ledger};
use fs_project::{
    Budgets, ConsequenceClass, Cooling, DecisionGate, EntityDecl, Envelope, GeometryArtifact,
    GeometryAssignment, HalfSpaceSide, MeshSelector, Metadata, OutputRequest, PowerDissipation,
    ProjectSpec, RequirementDirection, RequirementSeverity, RequirementSource,
    RequirementSourceKind, SafetyFactorPolicy, Seeds, SolverSettings, ThermalLimit, UnitsDoctrine,
    Versions,
};
use fs_qty::QtyAny;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x6a_03_01,
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

fn facet(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> String {
    format!(
        "facet normal 0 0 0\nouter loop\nvertex {} {} {}\nvertex {} {} {}\nvertex {} {} {}\nendloop\nendfacet\n",
        a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2],
    )
}

fn tetra_stl(dirty: bool, open: bool) -> Vec<u8> {
    let p0 = [0.0, 0.0, 0.0];
    let p1 = [1.0, 0.0, 0.0];
    let p2 = [0.0, 1.0, 0.0];
    let p3 = [0.0, 0.0, 1.0];
    let mut stl = String::from("solid enclosure\n");
    stl.push_str(&facet(p0, p2, p1));
    stl.push_str(&facet(p0, p1, p3));
    stl.push_str(&facet(p0, p3, p2));
    if !open {
        stl.push_str(&facet(p1, p2, p3));
    }
    if dirty {
        stl.push_str(&facet(p0, p2, p1));
        stl.push_str(&facet(p0, p0, p1));
    }
    stl.push_str("endsolid enclosure\n");
    stl.into_bytes()
}

fn tetra_step() -> Vec<u8> {
    b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('strict faceted fixture'),'2;1');
FILE_NAME('faceted.step','2026-07-23T00:00:00',('fs-cli'),('FrankenSim'),'fs-cli','FrankenSim','');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#60=FACETED_BREP('',#50);
#43=FACE('',(#42));
#21=POLY_LOOP('',(#1,#2,#4));
#2=CARTESIAN_POINT('',(1.0,0.0,0.0));
#42=FACE_OUTER_BOUND('',#41,.T.);
#12=FACE_OUTER_BOUND('',#11,.T.);
#31=POLY_LOOP('',(#1,#4,#3));
#4=CARTESIAN_POINT('',(0.0,0.0,1.0));
#50=CLOSED_SHELL('',(#43,#13,#33,#23));
#11=POLY_LOOP('',(#1,#3,#2));
#23=FACE('',(#22));
#1=CARTESIAN_POINT('',(0.0,0.0,0.0));
#32=FACE_OUTER_BOUND('',#31,.T.);
#41=POLY_LOOP('',(#2,#3,#4));
#13=FACE('',(#12));
#3=CARTESIAN_POINT('',(0.0,1.0,0.0));
#22=FACE_OUTER_BOUND('',#21,.T.);
#33=FACE('',(#32));
ENDSEC;
END-ISO-10303-21;
"
    .to_vec()
}

fn project_for_receipt(format: &str, source_hash: u64, parser_version: &str) -> ProjectSpec {
    let kelvin = |value| QtyAny::new(value, fs_project::spec::dims::TEMPERATURE);
    let watts = |value| QtyAny::new(value, fs_project::spec::dims::POWER);
    ProjectSpec {
        metadata: Some(Metadata {
            name: "import-reference".to_string(),
            created: "2026-07-23".to_string(),
            context_of_use: "geometry import conformance".to_string(),
            intended_decision: "exercise retained enclosure assignments".to_string(),
            decision_gate: DecisionGate::ScopingEstimate,
            consequence: ConsequenceClass::Advisory,
        }),
        versions: Some(Versions {
            schema: fs_project::FSIM_VERSION,
            constellation: "00".repeat(32),
            workspace: "11".repeat(20),
        }),
        seeds: Some(Seeds { master: 7 }),
        budgets: Some(Budgets {
            solve_time: QtyAny::new(60.0, fs_project::spec::dims::TIME),
            memory_bytes: 64 * 1024 * 1024,
            accuracy_rel: 0.01,
        }),
        capabilities: Some(vec!["geometry.import".to_string()]),
        units: Some(UnitsDoctrine {
            storage: "si-base".to_string(),
            display: "engineering".to_string(),
        }),
        geometry: Some(vec![GeometryArtifact {
            role: "enclosure".to_string(),
            format: format.to_string(),
            source_hash,
            parser_version: parser_version.to_string(),
        }]),
        assignments: Some(vec![GeometryAssignment {
            artifact: "enclosure".to_string(),
            target: "air".to_string(),
            length_unit: "m".to_string(),
            selector: MeshSelector::HalfSpace {
                normal: [1.0, 0.0, 0.0],
                offset: 1.0,
                side: HalfSpaceSide::AtMost,
                tolerance: 0.0,
            },
            allow_overlap: false,
        }]),
        assembly: Some(vec![
            EntityDecl::Assembly {
                name: "assembly".to_string(),
                display: "Assembly".to_string(),
                expect_id: None,
            },
            EntityDecl::Part {
                parent: "assembly".to_string(),
                name: "enclosure".to_string(),
                display: "Enclosure".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "enclosure".to_string(),
                name: "air".to_string(),
                display: "Internal air".to_string(),
                expect_id: None,
            },
        ]),
        materials: Some(Vec::new()),
        interface_cards: Some(Vec::new()),
        power: Some(vec![PowerDissipation {
            region: "air".to_string(),
            watts: watts(5.0),
            duty: 1.0,
        }]),
        cooling: Some(Cooling {
            fans: Vec::new(),
            vents: Vec::new(),
            leakage: watts(0.0),
        }),
        envelope: Some(Envelope {
            ambient_lo: kelvin(293.15),
            ambient_hi: kelvin(313.15),
            pressure: QtyAny::new(101_325.0, fs_project::spec::dims::PRESSURE),
        }),
        requirements: Some(vec![ThermalLimit {
            qoi: "temperature-max".to_string(),
            class: "surface".to_string(),
            region: "air".to_string(),
            direction: RequirementDirection::AtMost,
            limit: kelvin(353.15),
            margin: kelvin(5.0),
            source: RequirementSource {
                kind: RequirementSourceKind::UserDeclaration,
                document: "geometry-import-fixture".to_string(),
                version: "1".to_string(),
                locator: "temperature-max".to_string(),
            },
            safety_factor: SafetyFactorPolicy {
                factor: 1.0,
                source: RequirementSource {
                    kind: RequirementSourceKind::UserDeclaration,
                    document: "geometry-import-margin-policy".to_string(),
                    version: "1".to_string(),
                    locator: "factor".to_string(),
                },
            },
            severity: RequirementSeverity::ReliabilityDerating,
        }]),
        solver: Some(SolverSettings {
            fidelity: "auto".to_string(),
            tolerance_rel: 1e-6,
        }),
        outputs: Some(vec![OutputRequest {
            name: "temperature-max".to_string(),
            kind: "scalar".to_string(),
        }]),
    }
}

fn project_for(bytes: &[u8]) -> ProjectSpec {
    let receipt = import_mesh(bytes, "stl")
        .expect("fixture parses")
        .source_receipt;
    project_for_receipt("stl", receipt.source_hash, receipt.parser_version)
}

fn source_library(project: &ProjectSpec, bytes: Vec<u8>, unit: &str) -> RawGeometryLibrary {
    let artifact = &project.geometry.as_ref().expect("geometry")[0];
    let mut raw = RawGeometryLibrary::new();
    assert!(!raw.insert_mesh(
        artifact,
        "fixtures/enclosure.stl",
        bytes,
        unit,
        0,
        Vec::new(),
    ));
    raw
}

fn step_source_library(project: &ProjectSpec, bytes: Vec<u8>, target_h: f64) -> RawGeometryLibrary {
    let artifact = &project.geometry.as_ref().expect("geometry")[0];
    let mut raw = RawGeometryLibrary::new();
    assert!(!raw.insert_faceted_step(
        artifact,
        "fixtures/enclosure.step",
        bytes,
        60,
        "m",
        target_h,
        Vec::new(),
    ));
    raw
}

#[test]
fn g0_dirty_stl_promotes_assigns_and_retains_complete_lineage() {
    let bytes = tetra_stl(true, false);
    let project = project_for(&bytes);
    assert!(project.validate().is_empty());
    let raw = source_library(&project, bytes, "m");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let run = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect("dirty but repairable mesh imports")
    });

    assert_eq!(run.artifacts.len(), 1);
    assert!(run.assignment_table.contains("air | entity region:"));
    assert!(run.assignment_table.contains("| faces 4 |"));
    let retained = &run.artifacts[0];
    let receipt = ledger
        .get_artifact(&retained.promotion_receipt)
        .expect("receipt read")
        .expect("receipt exists");
    let receipt = String::from_utf8(receipt).expect("receipt utf8");
    assert!(receipt.contains("\"trust\":\"promoted\""));
    assert!(receipt.contains("\"class\":\"duplicate-face\""));
    assert_eq!(
        ledger
            .get_extension(ExtensionTable::Imports, &retained.import_record)
            .expect("imports read")
            .expect("imports row"),
        receipt
    );
    assert_eq!(
        ledger.op(run.op_id).expect("op read").expect("op").outcome,
        Some("ok".to_string())
    );
    assert!(
        ledger
            .edge_exists(run.op_id, &retained.raw_source, EdgeRole::In)
            .expect("raw lineage")
    );
    for hash in [
        retained.promotion_receipt,
        retained.promoted_mesh,
        retained.assignment_report,
        run.summary_artifact,
    ] {
        assert!(
            ledger
                .edge_exists(run.op_id, &hash, EdgeRole::Out)
                .expect("output lineage")
        );
    }
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g0_faceted_step_import_retains_both_receipts_and_repaired_mesh() {
    let bytes = tetra_step();
    let parsed = parse_step(&bytes).expect("STEP fixture parses");
    let project = project_for_receipt(
        "step",
        parsed.receipt().source_fingerprint(),
        STEP_FACETED_DECODER_VERSION,
    );
    assert!(project.validate().is_empty());
    let raw = step_source_library(&project, bytes, 1.0);
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let run = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect("strict faceted STEP imports")
    });

    assert!(run.assignment_table.contains("| faces 4 |"));
    let retained = &run.artifacts[0];
    let receipt = ledger
        .get_artifact(&retained.promotion_receipt)
        .expect("receipt read")
        .expect("receipt exists");
    let receipt = String::from_utf8(receipt).expect("receipt utf8");
    assert!(receipt.contains("frankensim.cli.faceted-step-import-receipt.v1"));
    assert!(receipt.contains("\"kind\":\"step-triangular-faceted-brep-receipt\""));
    assert!(receipt.contains("\"kind\":\"step-tessellation-to-sdf-receipt\""));
    assert!(receipt.contains("\"repaired_mesh\":{\"vertices\":4,\"triangles\":4}"));
    let promoted_mesh = ledger
        .get_artifact(&retained.promoted_mesh)
        .expect("mesh read")
        .expect("mesh exists");
    let promoted_mesh = String::from_utf8(promoted_mesh).expect("PLY utf8");
    assert!(promoted_mesh.contains("element vertex 4"));
    assert!(promoted_mesh.contains("element face 4"));
    let op = ledger.op(run.op_id).expect("op read").expect("op");
    assert!(op.ir.contains("\"kind\":\"faceted-step\""));
    assert!(op.ir.contains("\"root_id\":60"));
    assert!(op.versions.contains(STEP_FACETED_DECODER_VERSION));
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g3_repaired_retessellation_preserves_selector_stats() {
    let clean = tetra_stl(false, false);
    let dirty = tetra_stl(true, false);
    let clean_project = project_for(&clean);
    let dirty_project = project_for(&dirty);
    let clean_raw = source_library(&clean_project, clean, "m");
    let dirty_raw = source_library(&dirty_project, dirty, "m");
    let clean_ledger = Ledger::open(":memory:").expect("clean ledger");
    let dirty_ledger = Ledger::open(":memory:").expect("dirty ledger");
    let gate = CancelGate::new_clock_free();

    let clean_run = with_cx(&gate, |cx| {
        import_project_geometry(
            &clean_project,
            &clean_raw,
            &clean_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("clean source imports")
    });
    let dirty_run = with_cx(&gate, |cx| {
        import_project_geometry(
            &dirty_project,
            &dirty_raw,
            &dirty_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("repaired re-tessellation imports")
    });

    for (ledger, run) in [(&clean_ledger, &clean_run), (&dirty_ledger, &dirty_run)] {
        let report = ledger
            .get_artifact(&run.artifacts[0].assignment_report)
            .expect("assignment report read")
            .expect("assignment report exists");
        let report = String::from_utf8(report).expect("assignment report utf8");
        assert!(report.contains("\"face_count\":4"));
        assert!(report.contains("\"surface_area\":2.3660254037844384"));
        assert!(report.contains("\"enclosed_volume\":0.16666666666666666"));
        assert!(report.contains("\"bounds_min\":[0,0,0]"));
        assert!(report.contains("\"bounds_max\":[1,1,1]"));
    }
    assert!(clean_ledger.lint().expect("clean lint").is_clean());
    assert!(dirty_ledger.lint().expect("dirty lint").is_clean());
}

#[test]
fn g3_changed_source_bytes_refuse_and_record_the_attempt() {
    let expected = tetra_stl(false, false);
    let project = project_for(&expected);
    let changed = tetra_stl(false, true);
    let raw = source_library(&project, changed, "m");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("changed source identity must refuse")
    });

    assert_eq!(error.code, "cli-import-source-hash-mismatch");
    let recorded = error.recorded.expect("refusal retained");
    assert_eq!(recorded.raw_sources.len(), 1);
    assert!(recorded.receipt_artifacts.is_empty());
    assert!(recorded.promoted_meshes.is_empty());
    assert!(recorded.import_records.is_empty());
    assert_eq!(
        ledger
            .op(recorded.op_id)
            .expect("op read")
            .expect("op")
            .outcome,
        Some("error".to_string())
    );
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g3_open_mesh_refusal_retains_fs_io_receipt_and_fixes() {
    let bytes = tetra_stl(false, true);
    let project = project_for(&bytes);
    let raw = source_library(&project, bytes, "m");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("open mesh must not promote")
    });

    assert_eq!(error.code, "cli-import-promotion-refused");
    assert!(error.what.contains("non-manifold-or-open"));
    assert!(error.fix.contains("increase max_hole_edges"));
    let recorded = error.recorded.expect("refusal retained");
    let receipt_hash = recorded.receipt_artifacts[0];
    let receipt = ledger
        .get_artifact(&receipt_hash)
        .expect("receipt read")
        .expect("receipt exists");
    let receipt = String::from_utf8(receipt).expect("receipt utf8");
    assert!(receipt.contains("\"trust\":\"refused\""));
    assert_eq!(
        ledger
            .get_extension(ExtensionTable::Imports, &recorded.import_records[0],)
            .expect("imports read")
            .expect("imports row"),
        receipt
    );
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g3_named_groups_refuse_when_repair_moves_face_ordinals() {
    let bytes = tetra_stl(true, false);
    let project = project_for(&bytes);
    let artifact = &project.geometry.as_ref().expect("geometry")[0];
    let mut raw = RawGeometryLibrary::new();
    raw.insert_mesh(
        artifact,
        "fixtures/enclosure.stl",
        bytes,
        "m",
        0,
        vec![NamedFaceGroup {
            name: "DUPLICATE".to_string(),
            faces: vec![4],
        }],
    );
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("pre-repair ordinals must not be laundered onto the promoted soup")
    });

    assert_eq!(error.code, "cli-import-group-remap-unavailable");
    assert!(error.what.contains("ordinals no longer identify"));
    let recorded = error.recorded.expect("refusal retained");
    assert_eq!(recorded.receipt_artifacts.len(), 1);
    assert!(recorded.promoted_meshes.is_empty());
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g3_unit_mismatch_retains_completed_import_and_terminal_refusal() {
    let bytes = tetra_stl(false, false);
    let project = project_for(&bytes);
    let raw = source_library(&project, bytes, "mm");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("unit mismatch must refuse assignment")
    });

    assert_eq!(error.code, "project-assignment-unit-mismatch");
    assert!(error.fix.contains("coordinate unit"));
    let recorded = error.recorded.expect("post-import refusal retained");
    assert_eq!(recorded.raw_sources.len(), 1);
    assert_eq!(recorded.receipt_artifacts.len(), 1);
    assert_eq!(recorded.promoted_meshes.len(), 1);
    assert_eq!(recorded.import_records.len(), 1);
    assert_eq!(
        ledger
            .op(recorded.op_id)
            .expect("op read")
            .expect("op")
            .outcome,
        Some("error".to_string())
    );
    assert!(ledger.lint().expect("lint").is_clean());
}

#[test]
fn g3_dangling_assignment_refuses_before_import_side_effects() {
    let bytes = tetra_stl(false, false);
    let mut project = project_for(&bytes);
    project.assignments.as_mut().expect("assignments")[0].target = "missing-region".to_string();
    let raw = source_library(&project, bytes, "m");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("dangling assignment must fail project admission")
    });

    assert_eq!(error.code, "project-ref-unknown");
    assert!(error.fix.contains("declared names"));
    assert!(error.recorded.is_none());
    assert_eq!(ledger.table_count("artifacts").expect("artifacts"), 0);
    assert_eq!(ledger.table_count("ops").expect("ops"), 0);
}

#[test]
fn g5_independent_replays_produce_identical_content_identities() {
    let bytes = tetra_stl(false, false);
    let project = project_for(&bytes);
    let first_raw = source_library(&project, bytes.clone(), "m");
    let second_raw = source_library(&project, bytes, "m");
    let first_ledger = Ledger::open(":memory:").expect("first ledger");
    let second_ledger = Ledger::open(":memory:").expect("second ledger");
    let gate = CancelGate::new_clock_free();

    let first = with_cx(&gate, |cx| {
        import_project_geometry(
            &project,
            &first_raw,
            &first_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("first replay")
    });
    let second = with_cx(&gate, |cx| {
        import_project_geometry(
            &project,
            &second_raw,
            &second_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("second replay")
    });

    assert_eq!(first.project_hash, second.project_hash);
    assert_eq!(first.summary_artifact, second.summary_artifact);
    assert_eq!(first.assignment_table, second.assignment_table);
    assert_eq!(first.artifacts, second.artifacts);
}

#[test]
fn g5_step_sampling_policy_changes_frozen_operation_identity() {
    let bytes = tetra_step();
    let parsed = parse_step(&bytes).expect("STEP fixture parses");
    let project = project_for_receipt(
        "step",
        parsed.receipt().source_fingerprint(),
        STEP_FACETED_DECODER_VERSION,
    );
    let coarse_raw = step_source_library(&project, bytes.clone(), 1.0);
    let fine_raw = step_source_library(&project, bytes, 0.5);
    let coarse_ledger = Ledger::open(":memory:").expect("coarse ledger");
    let fine_ledger = Ledger::open(":memory:").expect("fine ledger");
    let gate = CancelGate::new_clock_free();

    let coarse = with_cx(&gate, |cx| {
        import_project_geometry(
            &project,
            &coarse_raw,
            &coarse_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("coarse replay")
    });
    let fine = with_cx(&gate, |cx| {
        import_project_geometry(
            &project,
            &fine_raw,
            &fine_ledger,
            GeometryImportLimits::DEFAULT,
            cx,
        )
        .expect("fine replay")
    });

    let coarse_ir = coarse_ledger
        .op(coarse.op_id)
        .expect("coarse op read")
        .expect("coarse op")
        .ir;
    let fine_ir = fine_ledger
        .op(fine.op_id)
        .expect("fine op read")
        .expect("fine op")
        .ir;
    assert_ne!(coarse_ir, fine_ir);
    assert!(coarse_ir.contains(&format!("{:016x}", 1.0f64.to_bits())));
    assert!(fine_ir.contains(&format!("{:016x}", 0.5f64.to_bits())));
    assert_ne!(coarse.summary_artifact, fine.summary_artifact);
}

#[test]
fn g4_precancelled_import_publishes_nothing() {
    let bytes = tetra_stl(false, false);
    let project = project_for(&bytes);
    let raw = source_library(&project, bytes, "m");
    let ledger = Ledger::open(":memory:").expect("ledger");
    let gate = CancelGate::new_clock_free();
    gate.request();

    let error = with_cx(&gate, |cx| {
        import_project_geometry(&project, &raw, &ledger, GeometryImportLimits::DEFAULT, cx)
            .expect_err("pre-cancelled import refuses")
    });

    assert_eq!(error.code, "cli-import-cancelled");
    assert!(error.recorded.is_none());
    assert_eq!(ledger.table_count("artifacts").expect("artifacts"), 0);
    assert_eq!(ledger.table_count("ops").expect("ops"), 0);
    assert_eq!(ledger.table_count("imports").expect("imports"), 0);
}
